<p align="center">
    <img src="imgs/RUMINA_LOGO_light.png#gh-light-mode-only" width=50% height=50%>
    <img src="imgs/RUMINA_LOGO_dark.png#gh-dark-mode-only" width=50% height=50%>
</p>

## RUMINA - Rust-based Unique Molecular Identifier Network Analysis

### Dependencies: 
- :snake: python ≥ 3.12
- :crab: cargo (rust) ≥ 1.77.0


### Installation:

1. Clone this repository
2. `cd` into cloned repo
3. run `sh install.sh`

This will compile the rust components of the pipeline, set up a python virtual environment with the necessary packages, and create a script named `rumina` to enable running RUMINA from any directory. This script will be located in `$HOME/.cargo/bin/`


### Usage: 
`rumina <input (file or directory)> --grouping_method <grouping_method> --separator <separator> <optional args>`

an example command:<br>
`rumina example.bam --grouping_method directional --separator : --report_coverage --delete_temps`

---
The `input` to `rumina` can be a file or a directory; if a directory, all BAM files within (exlcuding pipeline products) are processed.


###  Arguments 
:small_blue_diamond: = mandatory, no default

##### `input` :small_blue_diamond:
The input file or directory. If a file, it must be: 

1. in BAM format
2. sorted (e.g. via `samtools sort`)

BAM indexes or any files associated with reference genomes are not required.

If the input is a directory, all BAM files within (excluding pipeline products) will be processed per the other arguments specified. 

##### `--grouping_method` :small_blue_diamond:

Specifies how/if to merge UMIs based on edit distance, to account for PCR mutations and NGS errors in UMI sequence. Options are: 
* **directional**: Merge UMIs via directional clustering. See *Amplicon* section for more details.
* **raw**: Treat each UMI as genuine; UMIs are not merged.

##### `--separator` :small_blue_diamond:
Specifies the character in the read QNAME delimiting the UMI barcode from the rest of the string. This is usually `_` or `:`.<br>

##### `--split_window` (default = None)
dictates how to split input bam files into subfiles (for avoiding memory overflow). <br><br> This is usually necessary for bam files containing above ~15 million reads, assuming 16GB of total system memory, and has been used to process BAM files containing up to 110 million reads. <br> 

Splitting happens along coordinates of the reference genome in the specified BAM file; If `--split_window 100` was used, reads for every 100bp stretch of the reference would be stored in a separate subfile. These subfiles would be processed individually and merged together into the final output file. Once the final file has been created, the subfiles are deleted.

Because reads are grouped per reference coordinate regardless of splitting, this option does not change underlying analysis.

Options are: 
* **auto**: calculate the recommended subfile size (in basepairs along genome) based on input file size. If `input` is a directory, this will be applied independently to each file within the directory
* **positive integer from 1 - 500**: split input files by a fixed window size. If `input` is a directory, this will be applied to each file within the directory. 
* **none** (default): process the input as one file without splitting into subfiles. If your system has ~16GB of RAM, this is suitable for BAMs containing up to ~15 million reads. Using this option with larger bams may result in memory overuse.

##### `--delete_temps` (optional, off by default) 
deletes pipeline-generated files needed temporarily for processing.<br>Can save gigabytes of space when working with large files 
##### `--report_coverage` (optional, off by default)
generates coverage and depth .csv reports on output files using `bedtools genomecov` via `pybedtools`. Also records the UMI groups with the least and most reads, respectively, as well as the number of UMI groups surviving initial filtering. Doing this with large bam files can increase the runtime by several minutes per file

#### Todo

- establish UMI clustering methods for shotgun sequencing methods
- make illustrations for:
    - general pipeline process from input to output
    - bam file splitting
- actually have an explanation in the README for what the pipeline does

