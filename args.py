import argparse

def init_args():
    parser = argparse.ArgumentParser(
        prog="RUMINA",
        description="A pipeline to perform consensus-based error correction via UMI barcodes",
        usage="rumina [input] [grouping_method] [separator [SEPARATOR]] [--split_window [SPLIT_WINDOW]] [--delete_temps], [--report_coverage]"
    )

    parser.add_argument(
        'input')

    parser.add_argument(
        dest ='grouping_method',
        type = str,
        choices=['raw', 'directional'],
    )

    parser.add_argument(
        '--separator', action='store',
        type = str,
        required=True,
    )

    parser.add_argument(
        '--delete_temps',
        action = 'store_true')

    parser.add_argument(
        '--report_coverage', action = 'store_true')

    
    parser.add_argument('--split_window', nargs='?', default="auto")

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
