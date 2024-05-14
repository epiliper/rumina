import os 
from process import group_bam, build_onesies, remove_onesies, calculate_split, split_bam
from args import init_args
import time
from cov_reporter import report_coverage, summarize_coverage

args = init_args()

suffix = '.bam'
start_time = time.time()
split = False
window_size = 0
split_dirs = []

def process_dir(dir):
    for file in os.listdir(dir):
        if 'tagged' not in file and 'cleaned' not in file:
            file_to_clean = os.path.abspath(os.path.join(dir, file))
            print(f"WORKING ON FILE: {file_to_clean}")
            tagged_bam = group_bam(file_to_clean)
            bam_to_clean, blacklist = build_onesies(tagged_bam)
            clean_file = remove_onesies(bam_to_clean, blacklist)

def process_file(file):
    tagged_bam = group_bam(file)
    bam_to_clean, blacklist = build_onesies(tagged_bam)
    clean_file = remove_onesies(bam_to_clean, blacklist)

if os.path.isdir(args.input): 
    ## if a bam is split, a new directory of sub-bams is created
    ## and all sub-bams within are processed
    match args.split_window:

        ## calculate recommended split window size
        ## if zero, then just process all files in the dir
        case "auto":
            for file in os.listdir(args.input):
                window_size = calculate_split(file)
                if window_size == 0:
                    print("Processing file without splitting...")
                    process_dir(args.input)
                else:
                    split_dirs.append(split_bam(file, window_size))

            for dir in split_dirs:
                process_dir(dir)

        ## no splitting; process files normally
        case None:
            print("Processing directory without splitting...")
            process_dir(args.input)

        ## process all bams with a supplied split window size
        case args.split_window.isdigit():
            if int(args.split_window) == 0:
                print("Processing directory without splitting...")
                process_dir(args.input)

            else: 
                for file in os.listdir(args.input):
                    split_dirs.append(split_bam(file, args.split_window))

                for dir in split_dirs:
                    process_dir(dir)

        ## --split_window arg must be either "auto", a positive integer, or None (left blank)
        case _:
            print("Invalid value for split window size! Enter none, an integer, or 'auto'")
            exit(5)

if os.path.isfile(args.input):

    match args.split_window:

        ## calculate split
        case "auto":
            window_size = calculate_split(args.input)
            if window_size == 0:
                print("Processing file without splitting...")
                process_file(args.input)
            else:
                split_dir = split_bam(args.input, window_size)
                process_dir(split_dir)

        ## no split
        case None:
            print("Processing file without splitting...")
            process_file(args.input)

        ## split according to supplied --split_window value
        case args.split_window.isdigit():
            if int(args.split_window) == 0:
                print("Processing file without splitting...")
                process_file(args.input)
            else: 
                split_dir = split_bam(args.input, args.split_window)
                process_dir(split_dir)

        case _:
            print("Invalid value for split window size! Enter none, an integer, or 'auto'")
            exit(5)


end_time = time.time()

execution_time = end_time - start_time
print("Execution time: ", execution_time)
