import os 
from process import group_bam, build_onesies, remove_onesies, calculate_split, split_bam
from args import init_args
import time
from cov_reporter import report_coverage, summarize_coverage

args = init_args()

## last characters of input file, used to identify relevant files
suffix = '.bam'
start_time = time.time()

def process_dir(input_dir):
    for file in os.listdir(input_dir):
        if file.endswith(suffix):
            ## make sure inputs aren't already-processed files
            if 'tagged' not in file and 'cleaned' not in file:
                file_to_clean = os.path.abspath(os.path.join(input_dir, file))
                print(f"WORKING ON FILE: {file_to_clean}")
                tagged_bam = group_bam(file_to_clean)
                bam_to_clean, blacklist = build_onesies(tagged_bam)
                clean_file = remove_onesies(bam_to_clean, blacklist)

                if args.report_coverage:
                    report_coverage(clean_file)


## if input is a directory, process every bamfile within 
if os.path.isdir(args.input):
    split_dirs = []
    split = False
    for file in os.listdir(args.input):
        if file.endswith(suffix):

            match args.split_window:

                case "auto":
                    window_size = calculate_split(file)
                    if window_size == 0:
                        print("Processing file without splitting...")
                    else:
                        split_bam(file, window_size)

                case None:
                    print("Processing file without splitting...")

                case args.split_window.isdigit():
                    split_dirs.append(split_bam(file, args.split_window))
                    split = True

                case _:
                    print("Invalid value for split window size! Enter none, an integer, or 'auto'")
                    exit(5)


    for dir in split_dirs: 
        for file in dir:
            ## make sure inputs aren't already-processed files
            if 'tagged' not in file and 'cleaned' not in file:
                file_to_clean = os.path.abspath(os.path.join(args.input, file))
                print(f"WORKING ON FILE: {file_to_clean}")
                tagged_bam = group_bam(file_to_clean)
                bam_to_clean, blacklist = build_onesies(tagged_bam)
                clean_file = remove_onesies(bam_to_clean, blacklist)

                if args.report_coverage:
                    report_coverage(clean_file)

    ## assuming more than one input bamfile in the input directory, 
    ## compile coverage/depth report .tsvs into a single .csv
    if args.report_coverage:
        summarize_coverage('original', args.input)
        summarize_coverage('cleaned', args.input)

## if input is single file, process it
elif os.path.isfile(args.input):
    match args.split_window:

        case "auto":
            window_size = calculate_split(args.input)
            if window_size == 0:
                print("Processing file without splitting...")
            else:
                split_bam(args.input, window_size)

        case None:
            print("Processing file without splitting...")

        case _:
            if not args.split_window.isdigit():
                print("Invalid value for split window size! Enter none, an integer, or 'auto'")
            split_bam(args.input, args.split_window)
            exit(5)

    tagged_bam = group_bam(args.input)
    bam_to_clean, blacklist = build_onesies(tagged_bam)
    clean_file = remove_onesies(bam_to_clean, blacklist)

    if args.report_coverage:
        report_coverage(os.path.abspath(args.input))
        report_coverage(os.path.abspath(clean_file))

end_time = time.time()

execution_time = end_time - start_time
print("Execution time: ", execution_time)

