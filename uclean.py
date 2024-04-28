import argparse 
import pysam 
from collections import Counter 
import subprocess 
import os
import time
from cov_reporter import report_coverage, summarize_coverage

### import args
from args import init_args

args = init_args()

### Find parent directory of input .bam
work_path = args.input

if os.path.isabs(work_path):
    work_path = os.path.dirname(work_path)
else:
    work_path = os.path.abspath(
        os.path.dirname(args.input))

print(f"Working path: {work_path}")

### assign UG tag for each group of clustered UMIs
def tag_bam(input_file):

    # use umi_tools group to assign unique UG tag per UMI cluster
    tagged_file_name = input_file.split('.bam')[0] + '_tagged.bam'
    # tag_cmd = 'umi_tools group -I ' + input_file + " --output-bam --umi-separator=':' --paired -S " + tagged_file_name.split('.bam')[0] + '_temp.bam' 
    tag_cmd = 'bam_processor/target/release/bam_processor ' + input_file + ' ' + tagged_file_name.split('.bam')[0] + '_temp.bam ' + args.separator
    subprocess.run(tag_cmd, shell = True)

    # filter tagged bam to get only reads with UMI tag
    filter_cmd = 'samtools view -@ 6 -h -b -d UG ' + tagged_file_name.split('.bam')[0] + '_temp.bam > ' + tagged_file_name
    subprocess.run(filter_cmd, shell = True)

    # CLEAN 
    if args.delete_temps:
        print('DELETING FILE!')
        os.remove(tagged_file_name.split('.bam')[0] + '_temp.bam')

    return(tagged_file_name)

### read tagged bam and generate a list of UG tags that only appeared once. 
def build_onesies(input_file):
    name_of_txt = input_file.split('.bam')[0] + '_onesies.txt'

    print(f"TAGGED BAM: {input_file}")
    samfile = pysam.AlignmentFile(input_file, 'rb')
    iter = samfile.fetch(until_eof = True)
    ug_list = []
    n = 0
    for read in iter:
        # tag_report = {"UG":read.tags["UG"]}
        ug = str(dict(read.tags)["UG"])
        ug_list.append(ug)
        n += 1
        
    count_list = Counter(ug_list)

    filtered_list = dict(
        filter(lambda x: x[1] > 1, count_list.items())
    )
    print(f"Length of unfiltered bam: {len(count_list)}")
    print(f"Length of bam with once-observed UMI groups removed: {len(filtered_list)}")

    with open(name_of_txt, "w") as filter_file:
        for key in filtered_list.keys():
            filter_file.write(key + "\n")

    return input_file, name_of_txt

### Using file generated from build_onesies, 
### remove singlets
def remove_onesies(input_file, blacklist):
    file_to_clean = input_file

    output_file = os.path.abspath(input_file).split('/')[-1].split('.bam')[0]

    # save cleaned bamfiles to /cleaned folder within original bamfile directory
    output_dir = work_path + '/cleaned'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    clean_file_name = output_dir  + '/' + output_file + '_cleaned.bam'


    clean_cmd = 'samtools view -b -@ 6 -h -D UG:' + blacklist + ' ' + file_to_clean + ' > ' + clean_file_name
    subprocess.run(clean_cmd, shell = True)

    # CLEAN
    if args.delete_temps: 
        os.remove(file_to_clean)
        os.remove(blacklist)

    return clean_file_name

### check cleaned files to ensure no onesies remain
def check_cleaned(input_file):
    print("Checking...")

    samfile = pysam.AlignmentFile(input_file, 'rb')
    iter = samfile.fetch(until_eof = True)

    ug_list = []
    for read in iter:
        ug = str(dict(read.tags)['UG'])

        ug_list.append(ug)

    count_list = Counter(ug_list)
    print(count_list)

    qc = ''
    if 1 in count_list.keys():
        print('ERROR ERROR ERROR ERROR: ONESIE DETECTED, FILTERING FAILED!')

        qc = 'FAIL'
    else:
        print("Filtering successful; no onesies detected.")
        qc = 'PASS'

    return input_file, qc

def dedup(input_file):

    output_file = input_file.split('.bam')[0] + '_dedup.bam'
    
    subprocess.run('samtools index ' + input_file, shell = True)

    dedup_cmd = 'umi_tools dedup --paired --buffer-whole-contig -I input_file --umi-separator=":" -S output_file'
    dedup_cmd = dedup_cmd.replace('input_file', input_file)
    dedup_cmd = dedup_cmd.replace('output_file', output_file)
    subprocess.run(dedup_cmd, shell = True)

    subprocess.run('samtools index ' + output_file, shell = True)

    return output_file
