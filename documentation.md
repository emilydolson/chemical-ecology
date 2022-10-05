
## Files / folders

analysis
* contains R script to read in generated data and plots it
* see example of plot on page 18 of the grant
* TODO: convert script to python

ci
* stuff for cookiecutter / versioning???

config
* proof of concept sbatch files and interaction matrices
* TODO: investigate what the sb files + interaction matrices generate and how they are proof of concepts

docs
* more just cookiecutter stuff??
* TODO: update all the .md files eventually

include/chemical-ecology
* where the bulk of the work is
* split into config and a-eco (where the actual work is)
* TODO: document document document

scripts
* contains a python script that extracts the path taken by the fittest lineage from the phylogeny_5000.csv file from each directory
* TODO: is this even being used and if so where/how

source
* contains files to run the native or web stuff
* TODO: document

tests
* ??

third-party
* TODO: are Catch and force-cover being used and if so how

web
* display files for UI
* TODO why does "make web" change the files in here that were uploaded to github and why does that break the interaction matrix


# a-eco 

**Particle**
* Currently unused. Likely to be deleted. 

**CellData**
* A struct with members: (doubles) fitness, heredity, growth_rate, synergy, equilib_growth_rate, (int) count, (vector<int>) species. 
* Used in CalcAllFitness(), return type of GetFitness()

**AEcoWorld**
* Main class in which all relevant functions are created

Important globals of AEcoWorld:
**world_t**
* Vector of vectors of ints. Defines a 2D array to create world objects

**interactions**
* Vector of vectors of doubles. Contains interaction matrix for the world

**N_TYPES**
* Number of types of species we will have in this world. Set by config in Setup()

**MAX_POP**
* Maximum population size for one type in one cell. Set by config in Setup()

Misc
* Our world object holds 100 vectors of size 10. This is becasue we have a 10x10 grid of cells in our world, each with 10 types. 
* Types can move to adjacent cells with some probability. 