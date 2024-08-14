from process import (
    calculate_split,
    split_bam,
    merge_processed_splits,
    prepare_files,
    get_all_files,
    process_dir,
    process_file,
)

from logo import LOGO
from args import init_args
import time
from shutil import rmtree

print(LOGO)
args = init_args()

suffix = ".bam"
start_time = time.time()
window_size = 0
split_dirs = []


temp_bams = prepare_files(get_all_files(args.input))

for file in temp_bams:
    match args.split_window:
        # calculate recommended split window size for each file
        case "auto":
            # for file in os.listdir(args.input):
            window_size = calculate_split(file)

            if window_size == 0:
                print("Processing file without splitting...")
                process_file(file, split=False)
            else:
                file_split, split_dir = split_bam(file, window_size)
                split_dirs.append(split_dir)

                process_dir(split_dir, split=True)
                merge_processed_splits(file)

        # no splitting; process files normally
        case None:
            print("Processing file without splitting...")
            process_file(file, split=False)

        # process all bams with a supplied split window size
        case x if x.isdigit():
            if int(x) == 0:
                print("Processing file without splitting...")
                process_file(file, split=False)

            else:
                file_split, split_dir = split_bam(file, args.split_window)
                split_dirs.append(split_dir)

                process_dir(split_dir, split=True)

                merge_processed_splits(file)

        case _:
            print(
                "Invalid value for split window size! Enter none, an integer, or 'auto'"
            )
            exit(5)

[rmtree(dir) for dir in split_dirs]

end_time = time.time()

execution_time = end_time - start_time
print("Execution time: ", execution_time)
