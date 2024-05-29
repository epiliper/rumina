import pysam 
from collections import Counter 
import subprocess 
import os
from cov_reporter import report_coverage, summarize_coverage, report_merged_coverage

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


def calculate_split(input):
    size = int(os.stat(input).st_size / 1024**2)

    if size in range(0, 500):
        return 0

    elif size in range(500, 1_000):
        return 250

    elif size in range(1000, 10_000):
        return 100

    else:
        return 100

def split_bam(input, window_size):

    file = os.path.basename(input)
    split_dir = file.split('.bam')[0] + f"_{window_size}"

    print(f"Splitting {input} into {window_size}bp windows...")

    subprocess.run(["multibam/target/release/multibam", input, split_dir, str(window_size)])

    # return os.path.abspath(os.path.join(os.path.dirname(input), split_dir))
    return os.path.join(os.path.dirname(input), split_dir)

def merge_processed_splits():

    clean_dir = os.path.join(work_path, "cleaned")

    prefixes_for_merging = set()

    # prefixes_for_merging = [filename.split('.')[0] for filename in os.listdir(clean_dir)]
    for filename in os.listdir(clean_dir):
        if filename.endswith('.bam'):
            prefixes_for_merging.add(filename.split('.')[0])

    for prefix in prefixes_for_merging:
        splits = [os.path.join(clean_dir, file) for file in os.listdir(clean_dir) if file.startswith(prefix) and file.endswith('.bam')]
        final_file = os.path.join(clean_dir, prefix + "_final.bam")
        pysam.merge("-@ 6", final_file, *splits)
        for split in splits:
            os.remove(split)

        if args.report_coverage:
            report_merged_coverage(final_file)

### assign UG tag for each group of clustered UMIs
def group_bam(input_file, split):

    output_dir = os.path.join(work_path, 'cleaned')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    tagged_file_name = os.path.join(
        os.path.abspath(output_dir), 
        os.path.basename(input_file).split('.bam')[0] + '_cleaned.bam'
    )
    
    print(tagged_file_name)

    tag_cmd = 'snv/target/release/snv'

    subprocess.run([tag_cmd, input_file, tagged_file_name, args.separator])

    # pysam.view("-@ 6", "-h", "-b", "-d", "UG", os.path.abspath(temp_file_name), "-o", os.path.abspath(tagged_file_name), catch_stdout = False)

    if not split:

        min_groupsize = 0 
        max_groupsize = 0

        minmax_group_file = os.path.join(
            os.path.dirname(input_file),
            "minmax.txt"
        )

        with open(minmax_group_file) as minmax:
            min_and_max = minmax.readline().split('\t')
            min_groupsize = min_and_max[0]
            max_groupsize = min_and_max[1]

        os.remove(minmax_group_file)

        report_coverage(tagged_file_name, min_groupsize, max_groupsize)

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
        ug = (str(dict(read.tags)["UG"]), str(read.reference_start), str(read.qname.split(str(args.separator))[-1]))
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
            filter_file.write(key[0] +"\t" + key[1] + "\t" + key[2] + "\n")

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
        
    pysam.view("-b", "-@ 6", "-h", "-D", f"UG:{blacklist}", file_to_clean, "-o", clean_file_name, catch_stdout = False) 
    

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

def report():
    clean_path = os.path.join(work_path, "cleaned")
    for file in os.listdir(clean_path):
        if file.endswith('.bam'):
            report_coverage(os.path.join(clean_path, file))

def summarize():
    clean_path = os.path.join(work_path, "cleaned")
    summarize_coverage(clean_path)
    

    

