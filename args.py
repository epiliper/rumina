import argparse


def init_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'input')

    parser.add_argument(
        '--delete_temps',
        action = 'store_true')

    parser.add_argument(
        '--report_coverage', action = 'store_true')

    parser.add_argument(
        '--separator', action='store'
    )

    args = parser.parse_args()

    return parser.parse_args()
