import os 
import subprocess

cmd_clean = 'python3.11 uclean.py'

# directory = '/Users/eli/stuff/bioinformatics/CHONKYFILES/20240307_UMI3_SG_STW_manualupload/BAM'
directory = '/Volumes/GreningLab6/NEB_UMI_SG_3-NovaSeq/FASTQ/bam'

for file in os.listdir(directory):
    if file.endswith('_dedUMI.bam'):
        file_to_clean = os.path.abspath(os.path.join(directory, file))
        print(f"WORKING ON FILE: {file_to_clean}")
        subprocess.call(cmd_clean + ' ' + file_to_clean, shell = True)
        
