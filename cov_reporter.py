#! /usr/bin/env python3

import pandas as pd
import os
import warnings
import pysam
from concurrent.futures import ThreadPoolExecutor
import argparse

parser = argparse.ArgumentParser(description="A simple python wrapper over samtools coverage, with multithreading for parallel processing of references.")
parser.add_argument("infile", help = "bam file to get coverage reporting for")
parser.add_argument("--threads", default = os.cpu_count(), help = "number of threads to use", required = False)

warnings.simplefilter(action="ignore", category=FutureWarning)

def get_coverage(input_file, region):
    region_df = pd.DataFrame()

    cov_stats = pysam.coverage("-r", region, input_file)
    cov_stats = cov_stats.split("\n")
    cov_stats = zip(cov_stats[0].split("\t"), cov_stats[1].split("\t"))

    for col_val in cov_stats:
        region_df[col_val[0]] = [col_val[1]]

    return region_df


def generate_cov_depth_report(final_file):
    sort_and_index(final_file)
    cov_report = pd.DataFrame()
    with pysam.AlignmentFile(final_file, "rb") as outbam:
        references = outbam.references

    # use a thread for each region
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(get_coverage, final_file, region) for region in references
        ]
        for future in futures:
            cov_report = pd.concat([cov_report, future.result()], ignore_index=True)

    save_dir = os.path.dirname(final_file)
    cov_tsv = os.path.join(save_dir, final_file.split(".bam")[0] + "_coverage.tsv")
    cov_report.to_csv(cov_tsv, sep="\t", index=None)

def sort_and_index(output_file):
    print("sorting and indexing output BAM...\r")

    temp_file = output_file.split(".bam")[0] + "_s.bam"
    os.rename(output_file, temp_file)

    pysam.sort(f"-@ {args.threads}", temp_file, "-o", output_file)
    os.remove(temp_file)

    pysam.index(f"-@ {args.threads}", output_file)

if __name__ == "__main__":
    args = parser.parse_args()
    generate_cov_depth_report(args.infile)

