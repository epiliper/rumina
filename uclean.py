import argparse 
import pysam 
from collections import Counter 
import subprocess 
import os
import time

from cov_reporter import report_coverage, summarize_coverage

parser = argparse.ArgumentParser()

parser.add_argument(
    'input'
)

parser.add_argument(
    '--delete_temps',
    action = 'store_true'
)

parser.add_argument(
    '--report_coverage', action = 'store_true'
)

args = parser.parse_args()

def above_one(number):
    return number > 1


### Find parent directory of input .bam
work_path = args.input

if os.path.isabs(work_path):
    work_path = os.dirname(work_path)
else:
    work_path = os.path.abspath(
        os.path.dirname(args.input))

print(f"Working path: {work_path}")

### assign UG tag for each group of clustered UMIs
def tag_bam(input_file):

    # use umi_tools group to assign unique UG tag per UMI cluster
    tagged_file_name = input_file.split('.bam')[0] + '_tagged.bam'
    tag_cmd = 'umi_tools group -I ' + input_file + " --output-bam --umi-separator=':' --paired -S " + tagged_file_name.split('.bam')[0] + '_temp.bam' 
    subprocess.call(tag_cmd, shell = True)

    # filter tagged bam to get only reads with UMI tag
    filter_cmd = 'samtools view -@ 6 -h -b -d UG ' + tagged_file_name.split('.bam')[0] + '_temp.bam > ' + tagged_file_name
    subprocess.call(filter_cmd, shell = True)

    # CLEAN 
    if args.delete_temps:
        os.remove(tagged_file_name.split('.bam')[0] + '_temp.bam')

    return(tagged_file_name)

### read tagged bam and generate a list of UG tags that only appeared once. 
def build_onesies(input_file):

    print(f"TAGGED BAM: {input_file}")

    onesie_cmd = 'onesie_remover '
    onesie_cmd = onesie_cmd + input_file
    subprocess.call(onesie_cmd, shell = True)

    name_of_txt = input_file.split('.bam')[0] + '_onesies.txt'

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
    subprocess.call(clean_cmd, shell = True)

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
        print(ug)

    count_list = Counter(ug_list)
    print(count_list)

    qc = ''
    if 1 in count_list.keys():
        print('ERROR ERROR ERROR ERROR: ONESIE DETECTED, FILTERING FAILED! File name prefixed with "FAILED"')

        qc = 'FAIL'
    else:
        print("Filtering successful; no onesies detected.")
        qc = 'PASS'

    return input_file, qc

# driver code
start_time = time.time()
tagged_bam = tag_bam(args.input)
bam_to_clean, blacklist = build_onesies(tagged_bam)
clean_file = remove_onesies(bam_to_clean, blacklist)
# file_to_report, file_qc = check_cleaned(clean_file)

if args.report_coverage:
    # make_report(file_to_report, file_qc)
    report_coverage(args.input, 'original')
    # report_coverage(file_to_report, 'dedup')

    # summarize_coverage(file_to_report, 'original')
    # summarize_coverage(file_to_report, 'dedup')
    

end_time = time.time()

execution_time = end_time - start_time
print("Execution time: ", execution_time)

