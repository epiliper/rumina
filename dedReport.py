import pandas as pd
import os
import subprocess

df = pd.DataFrame()


cov_cmd = 'samtools coverage _input_file_ > out_file_dedreport.tsv'

for file in os.listdir():
    if file.endswith('_dedUMI.bam'):
        cmd_to_run = cov_cmd.replace('_input_file_', file)
        cmd_to_run = cmd_to_run.replace('out_file', file.split('dedUMI.bam')[0])
        subprocess.call(cmd_to_run, shell = True)



i = 0
for file in os.listdir():
    if 'dedreport.tsv' in file:
        read_file = pd.read_csv(file, sep = '\t', header = 0) 
        if i == 0:
            print(read_file)
            df = pd.DataFrame(columns = read_file.columns.tolist().append('sample'))
        df = pd.concat([df, read_file])
        df = df.reset_index(drop=True)
        df.loc[i, 'sample'] = file.split('.')[0]
        i += 1


df.to_csv(
    'final_ded_report.csv'
)
