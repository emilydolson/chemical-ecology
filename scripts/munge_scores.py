import sys
import csv
import glob
import os

# This script expects two command line arguments:
#   - a glob-style pattern indicating which directories to analyze
#         (remember to put quotation marks around it)
#   - the prefix filename to store data in 

def main():
    glob_pattern = sys.argv[1]
    outfilename = sys.argv[2]
    scores = [[], [], [], [], [], []]


    for dirname in glob.glob(glob_pattern):
        scores_file = dirname + "/scores.csv"
        if not (os.path.exists(scores_file)):
            print("Scores not found")
            continue

        with open(scores_file) as scores_file:
            reader = csv.reader(scores_file)
            data = list(reader)
            for i, score in enumerate(data[1]):
                scores[i].append(float(score))
        scores_file.close()

    avg_scores = []
    for fit_score in scores:
        avg_scores.append(sum(fit_score)/len(fit_score))

        
    with open(outfilename + ".csv", 'w') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(data[0])
        writer.writerow(avg_scores)
    outfile.close()


if __name__ == "__main__":
    main()