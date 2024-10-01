from process import (
    calculate_split,
    split_bam,
    merge_processed_splits,
    prepare_files,
    get_all_files,
    process_dir,
    process_file,
    print_file_info
)

from logo import r, M, C, print_logo
from args import init_args
import time
from shutil import rmtree

print_logo()
args = init_args()

suffix = ".bam"
start_time = time.time()
window_size = 0
split_dirs = []


print(f"{M}Initializing...{r}")
print(f"{M}============================={r}")
print("parameters:")
for arg in vars(args):
    print(f"{M}{arg}{r}: {getattr(args, arg)}")
temp_bams = prepare_files(get_all_files(args.input))

for i, file in enumerate(temp_bams, 1):
    print_file_info(len(temp_bams), i, file)
    match args.split_window:
        # calculate recommended split window size for each file
        case "auto":
            window_size = calculate_split(file)

            if window_size == 0:
                print("Processing file without splitting...")
                outbam = process_file(file, split=False)
            else:
                file_split, split_dir = split_bam(file, window_size)
                split_dirs.append(split_dir)

                process_dir(split_dir, split=True)
                outbam = merge_processed_splits(file)

        # no splitting; process files normally
        case None:
            print("Processing file without splitting...")
            outbam = process_file(file, split=False)

        # process all bams with a supplied split window size
        case x if x.isdigit():
            if int(x) == 0:
                print("Processing file without splitting...")
                outbam = process_file(file, split=False)

            else:
                file_split, split_dir = split_bam(file, args.split_window)
                split_dirs.append(split_dir)

                process_dir(split_dir, split=True)

                outbam = merge_processed_splits(file)

        case _:
            print(
                "Invalid value for split window size! Enter none, an integer, or 'auto'"
            )
            exit(5)

[rmtree(dir) for dir in split_dirs]

end_time = time.time()

execution_time = end_time - start_time
print("Execution time: ", execution_time)
