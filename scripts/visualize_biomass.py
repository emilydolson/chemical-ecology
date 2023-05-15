import sys
import ast
import re
import matplotlib.pyplot as plt
import pandas as pd

def main(file_path):
    df = pd.read_csv(f'{file_path}/a-eco_data.csv')
    world_states_str = df['worldState'].tolist()

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
                axis[row][col].plot(cell[species_i], label=species_i)
            else:
                axis[row][col].plot(cell[species_i])
        row += 1
        if row % 10 == 0:
            col += 1
            row = 0
    figure.legend(fontsize=32)
    plt.savefig('biomass_plot.png')
    plt.close()


if __name__ == '__main__':
    main(sys.argv[1])