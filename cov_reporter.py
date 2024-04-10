import pandas as pd 
import argparse
import subprocess
import os

columns = ['min_depth', 'max_depth', 'median_depth', 'mean_depth', 'coverage_percent','query_name']

def report_coverage(work_dir, input, report_type):

    cov_cmd = 'bedtools genomecov -d -ibam input_file > output_file'

    infile = os.path.basename(input)
    save_dir = os.path.dirname(input)
    outfile = input.split('.bam')[0] + '_depth.tsv'
    query_name = infile.split('.bam')[0]

    cov_cmd = cov_cmd.replace('input_file', input)
    cov_cmd = cov_cmd.replace('output_file', outfile)

    subprocess.run(cov_cmd, shell = True)

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
            [min_depth, max_depth, median_depth, mean_depth, coverage, query_name]]

    report = pd.DataFrame(data)

    if report_type == 'original':
        csv_name = infile.split('.bam')[0] + '_coverage.tsv'
    elif report_type == 'dedup':
        csv_name = infile.split('.bam')[0] + '_coverage_dedup.tsv'

    report.to_csv(save_dir + '/' + csv_name, sep = '\t', header = 0, index = None)
    
    os.remove(outfile)

def summarize_coverage(report_type, work_dir):

    if not os.path.isabs(work_dir):
        work_dir = os.path.abspath(work_dir)

    total_df = pd.DataFrame(columns = columns)
    
    final_file = ''
    suffix = ''
    search_path = ''

    if report_type == 'original':
        suffix = '_coverage.tsv'
        search_path = work_dir + '/'
        final_file = search_path + 'COMPLETE_COVERAGE.csv'

    elif report_type == 'dedup':
        suffix = '_coverage_dedup.tsv'
        search_path = work_dir + '/cleaned/'
        final_file = search_path + 'DEDUP_COMPLETE_COVERAGE.csv'

    print(f"Combining coverage reports...\nCoverage csv name: {final_file}")
    print(search_path)

    for file in os.listdir(search_path):
        if file.endswith(suffix):
            read_file = pd.read_csv(search_path + file, sep = '\t', header = 0) 
            total_df = pd.concat([total_df, read_file])
            # df = df.append(read_file)

    total_df.to_csv(final_file)

    for file in os.listdir(search_path):
        if file.endswith(suffix):
            os.remove(search_path + file)
    



    







    

