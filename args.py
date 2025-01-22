import argparse
import os
import sys

class CustomHelpFormatter(argparse.HelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings:
            return self._metavar_formatter(action, action.dest)(1)[0]
        else:
            parts = action.option_strings
            return ', '.join(parts)

def init_args():
    parser = argparse.ArgumentParser(
        # formatter_class=argparse.RawTextHelpFormatter,
        formatter_class=CustomHelpFormatter,
        prog="RUMINA",
        description="A pipeline to perform consensus-based error correction via UMI barcodes",
        usage="rumina [input] [grouping_method] [--separator [SEPARATOR]] [--split_window [SPLIT_WINDOW]] [--no_report]",
    )

    required_named = parser.add_argument_group("required named arguments")
    flags = parser.add_argument_group("flags")
    optional = parser.add_argument_group("optional arguments")

    required_named.add_argument(
        "input",
        help="A .bam file or a directory of .bam files to be processed.\n"
             "Bam files must have UMIs present in the read QNAME, AND be sorted and indexed with samtools.",
    )

    required_named.add_argument(
        "--grouping_method",
        required=True,
        type=str,
        choices=["raw", "directional", "acyclic"],
        help="Specifies UMI merging method for UMI error correction:\n"
             "- directional: Merge UMIs via directional clustering\n"
             "- acyclic: Same as directional, but networks are limited to a depth of one. May reduce overcorrection in barcode-diverse samples (e.g. metagenomics specimens)\n"
             "- raw: Treat each raw UMI as genuine; no merging",
    )

    required_named.add_argument(
        "--separator",
        action="store",
        type=str,
        required=True,
        help="Character in read QNAME delimiting UMI barcode. Usually '_' or ':'",
    )

    flags.add_argument(
        "--cov_depth_report",
        action="store_true",
        help="Calculate coverage/depth using 'bedtools genomecov (may add runtime for large files)",
    )

    flags.add_argument(
        "--length",
        action="store_true",
        help="Group reads by length in addition to reference coordinate",
    )

    flags.add_argument(
        "--sort_outbam",
        action="store_true",
        help="Sort and index output BAM file",
    )

    flags.add_argument(
        "--only-group",
        action="store_true",
        help="Group and tag reads without deduplication/error correction",
    )

    flags.add_argument(
        "--merge_pairs",
        help="Merge overlapping forward/reverse reads with same UMI. Requires reference genome FASTA for realignment"
    )

    flags.add_argument(
        "--min_overlap_bp",
        default="3",
        help="Minimum bases for read overlap (--merge_pairs)",
    )

    flags.add_argument(
        "--halve_pairs",
        action="store_true",
        help="Use only R1 for deduplication, discard R2, similar to UMI-tools",
    )

    flags.add_argument(
        "--singletons", 
        action="store_true",
        help="Process singleton reads"
    )

    # Optional arguments with improved formatting
    optional.add_argument(
        "--split_window",
        nargs="?",
        default="auto",
        help="BAM file splitting strategy: \n"
             "- 'auto (Default)': Recommended subfile size based on input\n"
             "- Integer (1-500): Fixed window size\n"
             "- 'none': Process entire file. May incur heavy memory usage.\n"
             "If you find an option to use too much memory, try increasing the split_window size.",
    )

    optional.add_argument(
        "--threads", 
        action="store", 
        type=str, 
        help="Number of threads to use"
    )
    
    optional.add_argument(
        "--outdir", 
        action="store", 
        type=str, 
        default="rumina_output",
        help="Output directory"
    )

    args = parser.parse_args()

    if args.threads:
        if int(args.threads) > os.cpu_count():
            sys.exit("Number of threads specified exceeds system threads. Exiting...")
    else:
        args.threads = str(os.cpu_count())

    if args.split_window is None:
        pass

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

        If you don't want to split your input (which may incur more memory usage), use this argument with 0.
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
