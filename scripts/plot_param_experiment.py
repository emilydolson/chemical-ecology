import pandas as pd
import matplotlib.pyplot as plt
import sys


def main():
    measures = ['Growth_Rate', 'Biomass', 'Heredity', 'Invasion_Ability']
    # measures = ['Growth_Rate', 'Biomass', 'Heredity']
    fig, axs = plt.subplots(len(measures), sharex=True)
    fig.suptitle('Growth of Fitness Proxies')
    name = '_'.join(sys.argv[1:])
    for condition in sys.argv[1:]: 
        file_location = "experiments/" + condition
        df = pd.read_csv(f'{file_location}/all_data.csv')
        df = df.groupby('Time').mean()
        diffusion = round(df.loc[0]['DIFFUSION'], 4)
        seeding_prob = round(df.loc[0]['SEEDING_PROB'], 4)
        clear_prob = round(df.loc[0]['PROB_CLEAR'], 4)
        params = "Diffusion=" + str(diffusion) +"  Seeding=" + str(seeding_prob) + "   Clear=" + str(clear_prob)
        for ax, measure in zip(axs, measures):
            ax.plot(df.index, df['mean_' + measure], label=params)
            ax.set_ylabel(measure)
        axs[0].legend(prop={'size': 6})
        axs[-1].set(xlabel='Updates')
        # plt.plot(df.index, df['mean_Fitness'], label=params)

    plt.show()
    for condition in sys.argv[1:]:
        plt.savefig(f'experiments/{condition}/{name}Plot.png')


if __name__ == '__main__':
    main()