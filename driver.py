### This is an example file showing how the pipeline
## can be applied to multiple input bamfiles 
## within a single directory.
import os 
import subprocess
import argparse
import pandas 
from cov_reporter import report_coverage, summarize_coverage
from uclean import *
from args import init_args
import time

args = init_args()

## last characters of input file, used to identify relevant files
suffix = '.bam'
start_time = time.time()

## if input is a directory, process every bamfile within 
if os.path.isdir(args.input):
    for file in os.listdir(args.input):
        if file.endswith(suffix):
            ## make sure inputs aren't already-processed files
            if 'tagged' not in file and 'cleaned' not in file:
                file_to_clean = os.path.abspath(os.path.join(args.input, file))
                print(f"WORKING ON FILE: {file_to_clean}")
                tagged_bam = tag_bam(file_to_clean)
                bam_to_clean, blacklist = build_onesies(tagged_bam)
                clean_file = remove_onesies(bam_to_clean, blacklist)
                file_to_dedup, file_qc = check_cleaned(clean_file)
                file_to_report = dedup(file_to_dedup) 
                if args.report_coverage: 
                    report_coverage(args.input, file_to_clean, 'original')
                    report_coverage(args.input, file_to_report, 'dedup')

    ## assuming more than one input bamfile in the input directory, 
    ## compile coverage/depth report .tsvs into a single .csv
    if args.report_coverage:
        summarize_coverage('original', args.input)
        summarize_coverage('dedup', args.input)

## if input is single file, process it
elif os.path.isfile(args.input):
    tagged_bam = tag_bam(args.input)
    bam_to_clean, blacklist = build_onesies(tagged_bam)
    clean_file = remove_onesies(bam_to_clean, blacklist)
    # file_to_dedup, file_qc = check_cleaned(clean_file)
    # file_to_report = dedup(file_to_dedup) 

end_time = time.time()

execution_time = end_time - start_time
print("Execution time: ", execution_time)

