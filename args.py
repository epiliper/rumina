import argparse

def init_args():
    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawTextHelpFormatter,
        prog="RUMINA",
        description="A pipeline to perform consensus-based error correction via UMI barcodes",
        usage="rumina [input] [grouping_method] [--separator [SEPARATOR]] [--split_window [SPLIT_WINDOW]] [--delete_temps], [--report_coverage]"
    )

    required_named = parser.add_argument_group('required named arguments')
    flags = parser.add_argument_group('flags')
    optional = parser.add_argument_group('optional arguments')

    required_named.add_argument(
        'input',
        help = 
"""A .bam file or a directory of .bam files to be processed. Bam files must have UMIs present in the read QNAME.
"""
    )

    required_named.add_argument(
        "--grouping_method",
        required=True,
        type = str,
        choices=['raw', 'directional'],
        help = 
"""Specifies how/if to merge UMIs based on edit distance, to account for PCR mutations and NGS errors in UMI sequence.
Options are:
- directional: Merge UMIs via directinal clustering.
- raw: Treat each raw UMI as genuine; UMIs are not merged.

"""
    )

    required_named.add_argument(
        '--separator', action='store',
        type = str,
        required=True,
        help = 
"""Specifies the character in the read QNAME delimiting the UMI barcode from the rest of the string. This is usually '_' or ':'. 

"""
    )

    flags.add_argument(
        '--delete_temps',
        action = 'store_true',
        help = 
"""If specified, deletes pipeline-generated files needed temporarily for processing. Can save gigabytes of space when working with large files.

"""
    )

    flags.add_argument(
        '--no_report', action = 'store_true', 
    help=
"""If used, disables coverage and depth reporting on output files using 'bedtools genomecov'. This can save several minutes per file when working with large files\n"""
    )

    
    optional.add_argument('--split_window', nargs='?', default="auto",
                          help = 
"""dictates how to split input bam files into subfiles (for avoiding memory overflow). Options are:
- auto: calcluate the recommended subfile size (in basepairs along genome) based on input file size. If 'input' is a directory, this will be calculated for each file within the directory.
- integer from 1-500: split input files by fixed window size. If 'input' is a directory, this will be applied to each file within the directory.
- none (default): process the input bam files as one file. If your system has ~16GB of ram, this is suitable for bams containing up to ~15 million reads. Using this option with larger bams without additional RAM may result in memory overuse.

"""
                          )

    args = parser.parse_args()
    if args.split_window is None:
        pass

    elif args.split_window.isdigit():
        if int(args.split_window) == 0:
            args.split_window = None

        elif int(args.split_window) < 0:
            print("Split window size cannot be negative. Please use a positive value.")
            exit(7)

    elif args.split_window == "auto":
        pass

    else:
        print(
            f"""
            Invalid value provided for --split_window: {args.split_window}. Options:\n
            - 'auto': Calculate a recommended window size, based on the size of the input bam file\n
            - any positive integer: a window size of your choosing\n
            
            If you don't want to split your input, don't use this argument. 
            """
            )
        exit(6)
    return args
