import pandas as pd 
import subprocess
import os
from shutil import which
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

columns = ['min umis per group', 'max umis per group', 'min_depth', 'max_depth', 'median_depth', 'mean_depth', 'coverage_percent','query_name']

def report_coverage(input, min_groupsize, max_groupsize):

    # infile = os.path.basename(input)
    # save_dir = os.path.dirname(os.path.basename(input))
    save_dir = os.path.dirname(input)
    outfile = os.path.basename(input).split('.bam')[0] + '_depth.tsv'
    query_name = os.path.basename(input).split('.bam')[0]

    with open(outfile, "w") as cov_file:
        subprocess.run([which("bedtools"), 'genomecov', '-d', '-ibam', os.path.abspath(input)], stdout = cov_file)

    df = pd.read_csv(outfile, sep='\t', names = ['reference', 'position', 'num_reads'])

    ## depth statistics
    min_depth = df['num_reads'].min()
    max_depth = df['num_reads'].max()
    median_depth = df['num_reads'].median()
    mean_depth = df['num_reads'].mean()

    ## coverage statistics
    num_positions = df['position'].max()
    coverage = 100 * (num_positions - len(df.loc[df['num_reads'] == 0])) / num_positions 

    data = [columns, 
            [min_groupsize, max_groupsize, min_depth, max_depth, median_depth, mean_depth, coverage, query_name]]

    report = pd.DataFrame(data)

    csv_name = os.path.join(
        save_dir,
        outfile.split('_depth.tsv')[0] + '_coverage.tsv'
    )

    # report.to_csv(os.path.join(save_dir, csv_name), sep = '\t', header = 0, index = None)
    report.to_csv(csv_name, sep = '\t', header = 0, index = None)
    
    os.remove(outfile)

def summarize_coverage(work_dir):

    if not os.path.isabs(work_dir):
        work_dir = os.path.abspath(work_dir)

    total_df = pd.DataFrame(columns = columns)
    
    final_file = os.path.join(work_dir, "COVERAGE_REPORT.csv")
    
    print(f"Combining coverage reports...\nCoverage csv name: {final_file}")

    for file in os.listdir(work_dir):
        if file.endswith('_coverage.tsv'):
            read_file = pd.read_csv(os.path.join(work_dir, file), sep = '\t', header = 0, index_col = False) 
            total_df = pd.concat([total_df, read_file])

    total_df.to_csv(final_file)

    for file in os.listdir(work_dir):
        if file.endswith('_coverage.tsv'):
            os.remove(os.path.join(work_dir, file))
    

def report_merged_coverage(input):
    work_path = os.path.dirname(input)

    minmax_file = os.path.join(
        work_path, 
        'minmax.txt'
    )

    df = pd.read_csv(minmax_file, sep = '\t', names = ['mins', 'maxes'])
    true_min = df['mins'].min()
    true_max = df['maxes'].max()

    report_coverage(input, true_min, true_max)

    os.remove(minmax_file)

    # os.remove(minmax_file)






    







    

