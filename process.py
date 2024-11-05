import pysam
import subprocess
import os
from cov_reporter import generate_report
from logo import print_file_end

from args import init_args

args = init_args()

# Find parent directory of input .bam
work_path = args.input

if os.path.isabs(work_path):
    work_path = os.path.dirname(work_path)
else:
    work_path = os.path.abspath(os.path.dirname(args.input))

exec_path = os.path.dirname(os.path.abspath(__file__))


def process_file(file, split_window):
    group_bam(file, split_window)
    print_file_end()


def get_all_files(input):
    print(f"Working path: {work_path}")
    files_to_clean = []

    print("gathering files...\r")

    if os.path.isdir(input):
        for file in os.listdir(input):
            if file.endswith(".bam") or file.endswith(".sam"):
                if "rumina" not in file and "tagged" not in file:
                    files_to_clean.append(os.path.join(input, file))

    else:
        files_to_clean.append(input)

    return files_to_clean


def prepare_files(files):
    print("converting SAMs to temporary BAMs and sorting inputs...\r")

    temp_files = []

    for file in files:
        if file.endswith(".sam"):
            bam_name = file.split(".sam")[0] + "_temp.bam"
            pysam.sort(f"-@ {args.threads}", "-o", bam_name, file)
            temp_files.append(bam_name)

        elif file.endswith(".bam"):
            temp_files.append(file)

    return temp_files


# calculate split window from file size
def calculate_split(input):
    # get file size in megabytes
    size = int(os.stat(input).st_size / 1024**2)

    if size in range(0, 500):
        return 0

    elif size in range(500, 1_000):
        return 250

    elif size in range(1000, 10_000):
        return 100

    else:
        return 100


def merge_fr(tagged_file_name):
    print("Merging overlapping forward/reverse amplicons...")
    dupes = os.path.join(work_path, args.outdir, "barcodes.tsv")
    outfile = tagged_file_name.split(".bam")[0] + "_merged.bam"

    subprocess.run(
        [
            os.path.join(exec_path, "pair_merger/target/release/pair_merger"),
            "-l",
            dupes,
            "-i",
            tagged_file_name,
            "-o",
            outfile,
            "--threads",
            args.threads,
        ],
    )

    os.remove(dupes)
    return outfile


# assign UG tag for each group of clustered UMIs
def group_bam(input_file, split_window):
    suffix = "_rumina.bam"

    output_dir = os.path.join(work_path, args.outdir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    tagged_file_name = os.path.join(
        os.path.abspath(output_dir),
        os.path.basename(input_file).split(".bam")[0] + suffix,
    )

    if args.only_group:
        tagged_file_name = tagged_file_name.split(".bam")[0] + "_GROUP_ONLY.bam"

    tag_cmd = os.path.join(exec_path, "bam_processor/target/release/bam_processor")
    tag_cmd = [
        tag_cmd,
        input_file,
        tagged_file_name,
        args.separator,
        args.grouping_method,
        args.threads,
    ]

    if split_window:
        tag_cmd.append(str(split_window))

    if args.length:
        tag_cmd.append("--length")

    if args.only_group:
        tag_cmd.append("--only-group")

    if args.singletons:
        tag_cmd.append("--singletons")

    if args.merge_fr:
        tag_cmd.append("--track-umis")

    subprocess.run(tag_cmd)

    if args.merge_fr:
        merge_file = merge_fr(tagged_file_name)
        os.remove(tagged_file_name)
        tagged_file_name = merge_file

    if args.sort_outbam:
        sort_and_index(tagged_file_name)

    if not args.no_report:
        generate_report(input_file, tagged_file_name)

    return tagged_file_name


def sort_and_index(output_file):
    print("sorting and indexing output BAM...\r")

    temp_file = output_file.split(".bam")[0] + "_s.bam"
    os.rename(output_file, temp_file)

    # if args.only_group:
    #     output_file = output_file.split(".bam")[0] + "_GROUP_ONLY.bam"

    pysam.sort(f"-@ {args.threads}", temp_file, "-o", output_file)
    os.remove(temp_file)

    pysam.index(f"-@ {args.threads}", output_file)
    print("getting coverage/depth report...\r")

    # generate_report(input_file, output_file)
