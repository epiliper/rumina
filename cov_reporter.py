import pandas as pd 
import os
import warnings
import pysam

import pybedtools

warnings.simplefilter(action="ignore", category=FutureWarning)

COLUMNS = [
    'num_reads_input_file',
    'num_reads_output_file',
    'num_raw_umis',
    'num_total_groups',
    'num_passing_groups', 
    'least_reads_group', 
    'min reads per group', 
    'most_reads_group', 
    'max reads per group', 
    'min_depth', 
    'max_depth', 
    'median_depth', 
    'mean_depth', 
    'coverage_percent',
    'query_name'
]

def generate_report(original_file, final_file):
    work_path = os.path.dirname(final_file)
    minmax_file = os.path.join(
        work_path, 
        'minmax.txt'
    )

    inreads, outreads = get_counts(original_file, final_file)

    df = pd.read_csv(minmax_file, sep = '\t', names = ['min_groups', 'mins', 'max_groups', 'maxes', 'num_passing_groups', 'num_total_groups', 'num_umis'])
    true_min = int(df['mins'].min())
    true_min_group = df.iloc[df['mins'].idxmin()].min_groups
    true_max = int(df['maxes'].max())
    true_max_group = df.iloc[df['maxes'].idxmax()].max_groups
    num_passing_groups = df['num_passing_groups'].sum()
    num_total_groups = df['num_total_groups'].sum()
    num_umis = df['num_umis'].sum()

    report_coverage(inreads, outreads, final_file, num_umis, num_total_groups, num_passing_groups, true_min_group, true_min, true_max_group, true_max)
    os.remove(minmax_file)


def get_counts(infile, outfile):
    infile_reads = pysam.AlignmentFile(infile).count(until_eof = True)
    outfile_reads = pysam.AlignmentFile(outfile).count(until_eof = True)

    return infile_reads, outfile_reads

def report_coverage(input_reads, output_reads, input, num_umis, num_total_groups, num_passing_groups, min_group, min_groupsize, max_group, max_groupsize):

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

    data = [COLUMNS, 
            [input_reads, output_reads, num_umis, num_total_groups, num_passing_groups, min_group, min_groupsize, max_group, max_groupsize, min_depth, max_depth, median_depth, mean_depth, coverage, query_name]]

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

    total_df = pd.DataFrame(columns = COLUMNS)
    
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
    
