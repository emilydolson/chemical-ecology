#!/bin/sh

#SBATCH -A ecode

## Email settings
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=foreba10@msu.edu

## FIRST ARG IS NAME OF EXPERIMENT
## SECOND ARG IS WHICH CLASS (class1 class2 etc)

## Job name settings
#SBATCH --job-name=help-me

## Time requirement in format "days-hours:minutes"
#SBATCH --time=6-0:00

## Memory requirement in megabytes. You might need to make this bigger.
#SBATCH --mem-per-cpu=5120                                                 

module load GNU/8.2.0-2.31.1

cd /mnt/home/$USER/chemical-ecology/searchresults
mkdir $1

cd /mnt/scratch/$USER/chemical-ecology/

## Make a directory for this condition
## (this line will throw an error once the directory exists, but
## that's harmless)
mkdir $1
cd $1

# Get executable
cp /mnt/home/$USER/chemical-ecology/chemical-ecology .

# Get config
cp -r /mnt/home/$USER/chemical-ecology/config .

mkdir search
cd search

# Get script
cp /mnt/home/$USER/chemical-ecology/scripts/search_variance.py .

# Run script
python3 search_variance.py $2 > outfile.txt

cp *.csv /mnt/home/$USER/chemical-ecology/searchresults/$1