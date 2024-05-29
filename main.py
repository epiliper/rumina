import os 
from process import group_bam, calculate_split, split_bam, merge_processed_splits

from args import init_args
import time
from shutil import rmtree

args = init_args()

suffix = '.bam'
start_time = time.time()
window_size = 0
split_dirs = []

def process_dir(dir, split):
    for file in os.listdir(dir):
        if 'tagged' not in file and 'cleaned' not in file:
            file_to_clean = os.path.abspath(os.path.join(dir, file))
            print(f"WORKING ON FILE: {file_to_clean}")
            tagged_bam = group_bam(file_to_clean, split)

def process_file(file):
    tagged_bam = group_bam(file, False)

## if input is a directory, process all bams within
if os.path.isdir(args.input): 

    match args.split_window:

        ## calculate recommended split window size
        ## if zero, then just process all files in the dir
        case "auto":
            for file in os.listdir(args.input):
                if file.endswith('.bam'):
                    window_size = calculate_split(file)
                    if window_size == 0:
                        print("Processing file without splitting...")
                        process_dir(args.input, split = False)
                    else:
                        split_dirs.append(split_bam(os.path.join(args.input, file), window_size))

            for dir in split_dirs:
                process_dir(dir, split = True)
                
            merge_processed_splits()

            [rmtree(dir) for dir in split_dirs]

        ## no splitting; process files normally
        case None:
            print("Processing directory without splitting...")
            process_dir(args.input, split = False)

        ## process all bams with a supplied split window size
        case x if x.isdigit():
            if int(x) == 0:
                print("Processing directory without splitting...")
                process_dir(args.input, split = False)

            else: 
                for file in os.listdir(args.input):
                    if file.endswith('.bam'):
                        split_dirs.append(split_bam(os.path.join(args.input, file), x))

                for dir in split_dirs:
                    clean_dir = process_dir(dir, split = True)

                    merge_processed_splits()

                [rmtree(dir) for dir in split_dirs]

        ## --split_window arg must be either "auto", a positive integer, or None (left blank)
        case _:
            print("Invalid value for split window size! Enter none, an integer, or 'auto'")
            exit(5)

## if input is just a file, process it
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
                process_dir(split_dir, split = True)
                merge_processed_splits()
                rmtree(split_dir)

        ## no split
        case None:
            print("Processing file without splitting...")
            process_file(args.input)

        ## split according to supplied --split_window value
        case x if x.isdigit():
            if int(x) == 0:
                print("Processing file without splitting...")
                process_file(args.input)
            else: 
                split_dir = split_bam(args.input, x)
                process_dir(split_dir, split = True)
                merge_processed_splits()
                rmtree(split_dir)

        case _:
            print("Invalid value for split window size! Enter none, an integer, or 'auto'")
            exit(5)

end_time = time.time()

execution_time = end_time - start_time
print("Execution time: ", execution_time)
