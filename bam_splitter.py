import pysam
import argparse
import numpy as np
import os
import subprocess
from itertools import pairwise

parser = argparse.ArgumentParser()

parser.add_argument('input_file')
parser.add_argument('output_dir')
parser.add_argument('chunk_size')
args = parser.parse_args() 

def get_max_pos(input_file):
    positions = set()
    for read in pysam.AlignmentFile(input_file, "rb"):
        positions.add(read.reference_start)

    max_pos = max(positions)
    return max_pos

def splitter(input_file, output_dir, chunk_size):
    chunk_size = int(chunk_size)
    output_dir = os.path.join(
        os.path.abspath(os.path.dirname(input_file)),
        output_dir
    )

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    print(f"outputting to {output_dir}")

    subfile_name = output_dir + '/' + input_file.split('bam')[0] + str({0}) + '.bam'
    max_pos = get_max_pos(input_file)
    ref_name = pysam.AlignmentFile(input_file, "rb").get_reference_name(0)
    split_cmd = 'samtools view -b -h {0} {1}:{2}-{3} > {4}'
    pysam.index(input_file)

    chunk_indices = np.arange(0, max_pos + 1, chunk_size).tolist()

    if len(chunk_indices) % 2 != 0: 
        leftover = max_pos % chunk_size 
        leftover = chunk_indices[-1] + leftover + 1
        chunk_indices.append(leftover)

    print(f"Length of reference{max_pos}")

    for i, index_pair in enumerate(pairwise(chunk_indices)):
        subfile = subfile_name.format(i)
        if i == 0:
            cmd_to_run = split_cmd.format(input_file, ref_name, index_pair[0], index_pair[1] - 1, subfile)
            subprocess.run(cmd_to_run, shell = True)
            print(cmd_to_run)
        elif i == 1:
            cmd_to_run = split_cmd.format(input_file, ref_name, index_pair[0], index_pair[1] - 1, subfile)
            subprocess.run(cmd_to_run, shell = True)
            print(cmd_to_run)
        elif i == -1: 
            cmd_to_run = split_cmd.format(input_file, ref_name, index_pair[0], index_pair[1], subfile)
            subprocess.run(cmd_to_run, shell = True)
            print(cmd_to_run)
        else:
            cmd_to_run = split_cmd.format(input_file, ref_name, index_pair[0], index_pair[1] - 1 , subfile)
            subprocess.run(cmd_to_run, shell = True)
            print(cmd_to_run)

splitter(args.input_file, args.output_dir, args.chunk_size)








