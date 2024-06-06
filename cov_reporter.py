import pandas as pd 
import os
import warnings
import pysam

import pybedtools

warnings.simplefilter(action="ignore", category=FutureWarning)

columns = [
    'num_reads_input_file',
    'num_reads_output_file',
    'least_reads_group', 
    'min umis per group', 
    'most_reads_group', 
    'max umis per group', 
    'num_groups', 
    'min_depth', 
    'max_depth', 
    'median_depth', 
    'mean_depth', 
    'coverage_percent',
    'query_name'
]

def get_counts(infile, outfile):
    infile_reads = pysam.AlignmentFile(infile).count(until_eof = True)
    outfile_reads = pysam.AlignmentFile(outfile).count(until_eof = True)

    return infile_reads, outfile_reads

def report_coverage(input_reads, output_reads, input, min_group, min_groupsize, max_group, max_groupsize, num_groups):

    # infile = os.path.basename(input)
    save_dir = os.path.dirname(input)
    outfile = os.path.basename(input).split('.bam')[0] + '_depth.tsv'
    query_name = os.path.basename(input).split('.bam')[0]

    ## run bedtools genomecov via pybedtools API
    pybedtools.example_bedtool(os.path.abspath(input)).genome_coverage(d=True).saveas(outfile)

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
            [input_reads, output_reads, min_group, min_groupsize, max_group, max_groupsize, num_groups, min_depth, max_depth, median_depth, mean_depth, coverage, query_name]]

    report = pd.DataFrame(data)

    csv_name = os.path.join(
        save_dir,
        outfile.split('_depth.tsv')[0] + '_coverage.tsv'
    )

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
    
def report_merged_coverage(original_file, final_file):
    work_path = os.path.dirname(final_file)

    minmax_file = os.path.join(
        work_path, 
        'minmax.txt'
    )

    inreads, outreads = get_counts(original_file, final_file)

    df = pd.read_csv(minmax_file, sep = '\t', names = ['min_groups', 'mins', 'max_groups', 'maxes', 'num_groups'])
    true_min = int(df['mins'].min())
    true_min_group = df.iloc[df['mins'].idxmin()].min_groups
    true_max = int(df['maxes'].max())
    true_max_group = df.iloc[df['maxes'].idxmax()].max_groups
    num_groups = df['num_groups'].sum()

    report_coverage(inreads, outreads, final_file, true_min_group, true_min, true_max_group, true_max, num_groups)
    os.remove(minmax_file)
