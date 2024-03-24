import pandas as pd
import subprocess 
import os
import argparse

df = pd.DataFrame(columns = ['Sample', 'min_depth', 'max_depth', 'median_depth', 'mean_depth'])

parser = argparse.ArgumentParser()

parser.add_argument(
    'suffix'
)

args = parser.parse_args()


depth_cmd = 'samtools depth _input_file_ > output_file.tsv'

counter = 0
for file in os.listdir():
    if file.endswith(args.suffix):
        cmd_to_run = depth_cmd.replace('_input_file_', file)
        cmd_to_run = cmd_to_run.replace('output_file', file.split('_')[0])

        subprocess.call(cmd_to_run, shell = True)

        tsv = pd.read_csv(file.split('_')[0] + '.tsv', sep = '\t')

        df.loc[counter, 'Sample'] = file.split('_')[0]
        df.loc[counter, 'min_depth'] = tsv.iloc[:, 2].min()
        df.loc[counter, 'max_depth'] = tsv.iloc[:, 2].max()
        df.loc[counter, 'median_depth'] = tsv.iloc[:, 2].median()
        df.loc[counter, 'mean_depth'] = tsv.iloc[:, 2].mean()
        counter += 1


df.to_csv('depth_report.csv')




