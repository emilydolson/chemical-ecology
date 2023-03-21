import plotly.express as px
import csv 
import pandas as pd
import plotly
import sys

def plot(dir):
    #measures = ["Biomass", "Growth_Rate", "Heredity", "Invasion_Ability", "Resiliance"]
    measures = ["Biomass"]
    for measure in measures:
        df = pd.read_csv(dir + "/" + measure + ".csv")
        fig = px.parallel_coordinates(df, color="Score",
                                    dimensions=['Seeding', 'Diffusion', 'Clear'],
                                    color_continuous_scale=px.colors.diverging.RdBu)
        fig.update_layout(title_text=dir + " " + measure, title_x=0.5, title_y=.025)
        plotly.offline.plot(fig, filename=f'/mnt/home/foreba10/chemical-ecology/searchresults/{dir}/{measure}.html')

def main():
    plot(sys.argv[1])

main()