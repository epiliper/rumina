import argparse 
import binsearch
import pysam
from collections import Counter
import subprocess


parser = argparse.ArgumentParser()

parser.add_argument(
    'input'
)

def above_one(number):
    return number > 1

args = parser.parse_args()

### assign UG tag for each group of clustered UMIs
def tag_bam(input_file):
    
    tagged_file_name = input_file.split('.bam')[0] + '_tagged.bam'
    tag_cmd = 'umi_tools group -I ' + input_file + " --output-bam --umi-separator=':' --paired -S " + tagged_file_name.split('.bam')[0] + '_temp.bam' 
    subprocess.call(tag_cmd, shell = True)
    filter_cmd = 'samtools view -h -b -d UG ' + tagged_file_name.split('.bam')[0] + '_temp.bam > ' + tagged_file_name
    subprocess.call(filter_cmd, shell = True)

    return(tagged_file_name)

### read tagged bam and generate a list of UG tags that only appeared once. 
### this list represents reads that only appeared once per amplicon, once per coordinate. 
### this list can be used to filter for reads that were observed more than once per UMI group
def build_onesies(input_file):

    print(f"TAGGED BAM: {input_file}")

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

    name_of_txt = input_file.split('.bam')[0] + '_onesies.txt'


    with open(name_of_txt, "w") as filter_file:
        for key in filtered_list.keys():
            filter_file.write(key + "\n")


    return input_file, name_of_txt

def remove_onesies(input_file, blacklist):
    file_to_clean = input_file
    clean_file_name = input_file.split('.bam')[0] + '_cleaned.bam'
    clean_cmd = 'samtools view -h -D UG:' + blacklist + ' ' + file_to_clean + ' > ' + clean_file_name
    subprocess.call(clean_cmd, shell = True)

# driver code
tagged_bam = tag_bam(args.input)
bam_to_clean, blacklist = build_onesies(tagged_bam)
remove_onesies(bam_to_clean, blacklist)

