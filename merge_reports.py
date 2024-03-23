import pandas as pd
import os

df = pd.DataFrame()


i = 0
for file in os.listdir():
    if '.tsv' in file:
        read_file = pd.read_csv(file, sep = '\t', header = 0) 
        if i == 0:
            print(read_file)
            df = pd.DataFrame(columns = read_file.columns.tolist().append('sample'))
        df = pd.concat([df, read_file])
        df = df.reset_index(drop=True)
        df.loc[i, 'sample'] = file.split('.')[0]
        i += 1


df.to_csv(
    'final_report.csv'
)
