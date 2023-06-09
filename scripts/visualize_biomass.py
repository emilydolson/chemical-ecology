import sys
import ast
import re
import matplotlib.pyplot as plt
import pandas as pd


def format_and_plot(world_states_str, plot_name):
    #convert string representation into list
    world_states = []
    for world_state in world_states_str:
        string_with_commas = world_state.replace('] [', '],[')
        string_with_commas = re.sub(r'(\d+)\s+', r'\1, ', string_with_commas)
        world_state = ast.literal_eval(string_with_commas)
        world_states.append(world_state)

    #reformat to plot
    biomasses = [[[0 for _ in range(len(world_states))] for _ in range(len(world_states[0][0]))] for _ in range(len(world_states[0]))]
    for world_state_i in range(len(world_states)):
        for cell_i in range(len(world_states[world_state_i])):
            for species_i in range(len(world_states[world_state_i][cell_i])):
                biomasses[cell_i][species_i][world_state_i] = world_states[world_state_i][cell_i][species_i]

    #plot
    figure, axis = plt.subplots(10, 10, figsize=(50,50))
    row = 0
    col = 0
    for cell in biomasses:
        for species_i in range(len(cell)):
            if row == 0 and col == 0:
                axis[row][col].plot(cell[species_i], label=species_i, alpha=0.75)
            else:
                axis[row][col].plot(cell[species_i], alpha=0.75)
        row += 1
        if row % 10 == 0:
            col += 1
            row = 0
    figure.legend(fontsize=32)
    plt.savefig(plot_name)
    plt.close()


def main(file_path):
    df = pd.read_csv(f'{file_path}/a-eco_data.csv')
    world_states_str = df['worldState'].tolist()
    format_and_plot(world_states_str, 'biomass_plot.png')

    df = pd.read_csv(f'{file_path}/a-eco_model_data.csv')
    df_repro = df.loc[df['worldType'] == 'Repro']
    df_soup = df.loc[df['worldType'] == 'Soup']

    soup_world_states_str = df_soup['stochasticWorldState'].tolist()
    format_and_plot(soup_world_states_str, 'biomass_plot_soup.png')
    repro_world_states_str = df_repro['stochasticWorldState'].tolist()
    format_and_plot(repro_world_states_str, 'biomass_plot_repro.png')


if __name__ == '__main__':
    main(sys.argv[1])