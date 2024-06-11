import os 
from process import group_bam, calculate_split, split_bam, merge_processed_splits, prepare_files

from args import init_args
import time
from shutil import rmtree

args = init_args()

suffix = '.bam'
start_time = time.time()
window_size = 0
split_dirs = []

def process_dir(dir, split):
    temp_bams = []

    for file in os.listdir(dir):
        if file.endswith('.bam') and 'tagged' not in file and 'cleaned' not in file:
            file_to_clean = os.path.abspath(os.path.join(dir, file))
            print(f"WORKING ON FILE: {file_to_clean}")
            tagged_bam = group_bam(file_to_clean, split)

def process_file(file):
    tagged_bam = group_bam(file, False)

## if input is a directory, process all bams within
if os.path.isdir(args.input): 

    temp_bams = prepare_files(args.input)

    match args.split_window:

        ## calculate recommended split window size for each file
        case "auto":
            for file in os.listdir(args.input):
                if file.endswith('.bam'):
                    window_size = calculate_split(os.path.join(args.input, file))

                    if window_size == 0:
                        print("Processing file without splitting...")
                        process_dir(args.input, split = False)
                        break
                    else:
                        file_split, split_dir = split_bam(os.path.join(args.input, file), window_size)
                        split_dirs.append(split_dir)

                        process_dir(dir, split = True)
                        merge_processed_splits(file)

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
                        file_split, split_dir = split_bam(os.path.join(args.input, file), args.split_window)
                        split_dirs.append(split_dir)

                        process_dir(split_dir, split = True)
                        
                        merge_processed_splits(file_split)
                [rmtree(dir) for dir in split_dirs]

        case _:
            print("Invalid value for split window size! Enter none, an integer, or 'auto'")
            exit(5)

    [os.remove(temp) for temp in temp_bams]

## if input is just a file, process it
if os.path.isfile(args.input):

    temp_bams = prepare_files(args.input)

    if not len(temp_bams) == 0:
        input = temp_bams[0]
    else: 
        input = args.input

    match args.split_window:

        ## calculate split
        case "auto":
            window_size = calculate_split(input)
            if window_size == 0:
                print("Processing file without splitting...")
                process_file(input)
            else:
                input, split_dir = split_bam(input, window_size)
                process_dir(split_dir, split = True)
                merge_processed_splits(input)
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
                input, split_dir = split_bam(args.input, x)
                print(input, split_dir)
                process_dir(split_dir, split = True)
                merge_processed_splits(input)
                rmtree(split_dir)

        case _:
            print("Invalid value for split window size! Enter none, an integer, or 'auto'")
            exit(5)

    [os.remove(temp) for temp in temp_bams]

end_time = time.time()

execution_time = end_time - start_time
print("Execution time: ", execution_time)
