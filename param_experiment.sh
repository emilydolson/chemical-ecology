#!/bin/bash

## arguments: name of sbatch script, name of condition
## This will automatically run the sbatch script passed in as an argument to this file, munge the data, and visualize the best linage
## Create all sbatch scripts following the templates in /config
## Run this script using "./param_experiment.sh config/param_test.sb param"

## Modified version of run_experiment. This way we dont have to worry about what individual SCRATCH directories look like, 
## and parameters can be individually plotted without clogging interaction graphs

mkdir experiments/$2
sbatch --wait $1 /mnt/home/$USER/chemical-ecology $2 
rm *.out
python3 scripts/munge_data.py "/mnt/scratch/$USER/chemical-ecology/data/$2/*" all_data
mv all_data.csv experiments/$2
python3 scripts/munge_scores.py "/mnt/scratch/$USER/chemical-ecology/data/$2/*" all_scores
mv all_scores.csv experiments/$2
python3 scripts/plot_param_experiment.py $2

## plot_param_experiment can be run with more than one argument to make combined graphs