import pysam 
from collections import Counter 
import subprocess 
import os

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
    size = os.stat(input).st_size / 1024**2
    print(size)
    window_size = 0

    if size in range(0, 500):
        window_size = 0 

    elif size in range(500, 1_000):
        window_size = 250

    elif size in range(1000, 10_000):
        print("oof")
        window_size = 100

    else:
        window_size = 100

    return window_size 

def split_bam(input, window_size):

    file = os.path.basename(input)
    split_dir = file.split('.bam')[0] + f"_{window_size}"

    print(f"Splitting {input} into {window_size} base pair windows...")

    subprocess.run(["multibam/target/release/multibam", input, split_dir, str(window_size)])

    return os.path.abspath(os.path.join(os.path.dirname(input), split_dir))

### assign UG tag for each group of clustered UMIs
def group_bam(input_file):

    tagged_file_name = input_file.split('.bam')[0] + '_tagged.bam'
    temp_file_name = tagged_file_name.split('.bam')[0] + '_temp.bam'
    tag_cmd = 'bam_processor/target/release/bam_processor'
    subprocess.run([tag_cmd, input_file, temp_file_name, args.separator])

    pysam.view("-@ 6", "-h", "-b", "-d", "UG", os.path.abspath(temp_file_name), "-o", os.path.abspath(tagged_file_name), catch_stdout = False)

    # CLEAN 
    if args.delete_temps:
        os.remove(temp_file_name)

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
