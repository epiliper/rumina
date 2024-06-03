import argparse

def init_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'input')

    parser.add_argument(
        'grouping_method',
        choices=['raw', 'directional']
    )

    parser.add_argument(
        '--delete_temps',
        action = 'store_true')

    parser.add_argument(
        '--report_coverage', action = 'store_true')

    parser.add_argument(
        '--separator', action='store'
    )
    
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
