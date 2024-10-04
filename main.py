from process import (
    calculate_split,
    prepare_files,
    get_all_files,
    process_file,
)

from logo import r, M, print_logo, print_file_info
from args import init_args
import time

print_logo()
args = init_args()

suffix = ".bam"
start_time = time.time()
window_size = 0

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
                outbam = process_file(file, None)
            else:
                process_file(file, window_size)

        # no splitting; process files normally
        case None:
            print("Processing file without splitting...")
            outbam = process_file(file, None)

        # process all bams with a supplied split window size
        case window_size if window_size.isdigit():
            if int(window_size) == 0:
                print("Processing file without splitting...")
                outbam = process_file(file, None)

            else:
                process_file(file, window_size)

        case _:
            print(
                "Invalid value for split window size! Enter none, an integer, or 'auto'"
            )
            exit(5)

end_time = time.time()

execution_time = end_time - start_time
print("Execution time: ", execution_time)
