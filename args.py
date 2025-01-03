import argparse
import os
import sys


def init_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="RUMINA",
        description="A pipeline to perform consensus-based error correction via UMI barcodes",
        usage="rumina [input] [grouping_method] [--separator [SEPARATOR]] [--split_window [SPLIT_WINDOW]] [--no_report]",
    )

    required_named = parser.add_argument_group("required named arguments")
    flags = parser.add_argument_group("flags")
    optional = parser.add_argument_group("optional arguments")

    required_named.add_argument(
        "input",
        help="""A .bam file or a directory of .bam files to be processed. Bam files must have UMIs present in the read QNAME.
""",
    )

    required_named.add_argument(
        "--grouping_method",
        required=True,
        type=str,
        choices=["raw", "directional", "acyclic"],
        help="""Specifies how/if to merge UMIs based on edit distance, to account for PCR mutations and NGS errors in UMI sequence.
Options are:
- directional: Merge UMIs via directinal clustering.
- raw: Treat each raw UMI as genuine; UMIs are not merged.

""",
    )

    required_named.add_argument(
        "--separator",
        action="store",
        type=str,
        required=True,
        help="""Specifies the character in the read QNAME delimiting the UMI barcode from the rest of the string. This is usually '_' or ':'. 

""",
    )

    flags.add_argument(
        "--cov_depth_report",
        action="store_true",
        help="""Calculate coverage and depth reporting on output files using 'bedtools genomecov'. This can add several minutes of runtime per file when working with large files\n""",
    )

    flags.add_argument(
        "--no_report",
        action="store_true",
        help="""Disables grouping statistics reporting. Not computationaly intensive, but reduces output file count.""",
    )

    flags.add_argument(
        "--length",
        action="store_true",
        help="if used, groups reads by length in addition to reference coordinate.",
    )

    flags.add_argument(
        "--sort_outbam",
        action="store_true",
        help="sort and index the output BAM file. Required for generating coverage/depth reports (see --no_report).",
    )

    flags.add_argument(
        "--only-group",
        action="store_true",
        help="if used, reads are grouped and tagged but not deduplicated or error corrected. This is useful if you want to manually review what grouping looks like.",
    )

    flags.add_argument(
        "--merge_pairs",
        help="merge forward and reverse reads that overlap, with the same (corrected) UMI. Specify a reference genome (in FASTA format) for realignment of merged reads",
    )

    flags.add_argument(
        "--min_overlap_bp",
        default="3",
        help="minimum number of bases for reads to overlap to be merged via --merge_pairs. Reads below this threshold will be discarded.",
    )

    flags.add_argument(
        "--halve_pairs",
        action="store_true",
        help="only use R1 for deduplication, and discard R2. Similar to UMI-tools.",
    )

    flags.add_argument("--singletons", action="store_true")

    optional.add_argument(
        "--split_window",
        nargs="?",
        default="auto",
        help="""dictates how to split input bam files into subfiles (for avoiding memory overflow). Options are:
- auto: calcluate the recommended subfile size (in basepairs along genome) based on input file size. If 'input' is a directory, this will be calculated for each file within the directory.
- integer from 1-500: split input files by fixed window size. If 'input' is a directory, this will be applied to each file within the directory.
- none (default): process the input bam files as one file. If your system has ~16GB of ram, this is suitable for bams containing up to ~15 million reads. Using this option with larger bams without additional RAM may result in memory overuse.

""",
    )

    optional.add_argument("--threads", action="store", type=str)
    optional.add_argument("--outdir", action="store", type=str, default="rumina_output")

    args = parser.parse_args()

    if args.threads:
        if int(args.threads) > os.cpu_count():
            sys.exit("Number of threads specified exceeds system threads. Exiting...")
    else:
        args.threads = str(os.cpu_count())

    if args.split_window is None:
        pass

    if args.cov_depth_report and not args.sort_outbam:
        sys.exit(
            "Cannot generate coverage and depth reports without sorting output files. Either use --sort_outbam or disable reporting with --no_report."
        )

    elif args.split_window.isdigit():
        if int(args.split_window) == 0:
            args.split_window = None

        elif int(args.split_window) < 0:
            sys.exit(
                "Split window size cannot be negative. Please use a positive value."
            )

    elif args.split_window == "auto":
        pass

    else:
        sys.exit(f"""
        Invalid value provided for --split_window: {args.split_window}. Options:\n
        - 'auto': Calculate a recommended window size, based on the size of the input bam file\n
        - any positive integer: a window size of your choosing\n

        If you don't want to split your input, don't use this argument.
        """)

    if args.merge_pairs:
        if args.halve_pairs:
            sys.exit(
                "Cannot use --merge_pairs and --halve_pairs simultaneously. Please pick one option."
            )

            if args.min_overlap_bp < 0:
                sys.exit(
                    "Cannot use negative value for --min_overlap_bp. Please pick a positive value"
                )

        if not any(ext in str(args.merge_pairs) for ext in [".fa", ".fasta"]):
            sys.exit(
                "Error: Reference file provided via --merge_pairs doesn't appear to be in FASTA format"
            )
        elif not os.path.exists(args.merge_pairs):
            sys.exit("Error: cannot find --merge_pairs reference FASTA")

    return args
