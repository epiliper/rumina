import argparse 
import binsearch
import pysam
from collections import Counter


parser = argparse.ArgumentParser()

parser.add_argument(
    'input'
)

def above_one(number):
    return number > 1

args = parser.parse_args()

def read_sam(input_file):

    samfile = pysam.AlignmentFile(input_file, 'rb')

    iter = samfile.fetch(until_eof = True)

    ug_list = []

    n = 0
    for read in iter:
        # tag_report = {"UG":read.tags["UG"]}
        ug = str(dict(read.tags)["UG"])

        ug_list.append(ug)

        print(ug)
        print(n)
        n += 1
        
    count_list = Counter(ug_list)

    # filtered_list = dict(filter(above_one, ug_list))
    filtered_list = dict(
        filter(lambda x: x[1] > 1, count_list.items())
    )

    print(f"Length of unfiltered bam: {len(count_list)}")
    print(f"Length of bam with once-observed UMI groups removed: {len(filtered_list)}")


    print(filtered_list)




    with open("test.txt", "w") as filter_file:
        for key in filtered_list.keys():
            filter_file.write(key + "\n")








read_sam(args.input)

