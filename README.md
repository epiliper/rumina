<p align="center">
    <img src="imgs/RUMINA_LOGO_light.png#gh-light-mode-only" width=50% height=50%>
    <img src="imgs/RUMINA_LOGO_dark.png#gh-dark-mode-only" width=50% height=50%>
</p>

## RUMINA - Rust-based Unique Molecular Identifier Network Analysis

- removing singleton UMIs (potential artifacts)
- deduplicating UMIs via UMI_tools 

### Installation:

1. Clone this repository
2. cd into `rumina`
3. run sh `install.sh`

This will compile the rust crates comprising the pipeline, set up a python virtual environment with the necessary packages, and create a script to run RUMINA from any directory: This script will named: `$HOME/.local/bin/rumina`.

After this, you should be able to run `rumina` in your shell of choice with the necessary arguments (see **Usage**).


### Usage: 
`rumina <input (file or directory)> --grouping_method <grouping_method> --separator <separator> <optional args>`

an example command:<br>
`rumina example.bam --grouping_method directional --separator : --report_coverage --delete_temps`

---
The `input` to `main.py` can be a file or a directory; in the latter case, all .bam files within a directory (exlcuding pipeline products) are processed. This is not yet parallelized. 

### Requirements: 
- python3.12+

###  Arguments 
:small_blue_diamond: = mandatory, no default

##### `input` :small_blue_diamond:
The input file or directory. If a file, it must be in .bam format. BAM indexes, or any files associated with the reference, are not required.

If the input is a directory, all .bam files within (excluding pipeline products) will be processed per the other arguments specified. 

##### `--grouping_method` :small_blue_diamond:

Specifies how/if to merge UMIs based on edit distance, to account for PCR mutations and NGS errors in UMI sequence. Options are: 
* **directional**: Merge UMIs via directional clustering. See *Amplicon* section for more details.
* **raw**: Treat each UMI as genuine; UMIs are not merged.

##### `--separator` :small_blue_diamond:
Specifies the character in the read QNAME delimiting the UMI barcode from the rest of the string. This is usually `_` or `:`.<br>


##### `--split_window` (default = None)
dictates how to split input bam files into subfiles (for avoiding memory overflow). Options are: 
* **auto**: calculate the recommended subfile size (in basepairs along genome) based on input file size. If `input` is a directory, this will be applied independently to each file within the directory
* **positive integer from 1 - 500**: split input files by a fixed window size. If `input` is a directory, this will be applied to each file within the directory. 
* **none** (default): process the input as one file without splitting into subfiles. If your system has ~16GB of RAM, this is suitable for bams containing up to 1e7 reads. Using this option with larger bams may result in memory overuse.

##### `--delete_temps` (optional, off by default) 
deletes pipeline-generated files needed temporarily for processing.<br>Can save gigabytes of space when working with large files 
##### `--report_coverage` (optional, off by default)
generates coverage and depth .csv reports on output files using `bedtools genomecov` via `pybedtools`. Doing this with large bam files can increase the runtime by several minutes per file

#### Todo

- establish UMI clustering methods for shotgun sequencing methods
- make illustrations for:
    - general pipeline process from input to output
    - bam file splitting
- actually have an explanation in the README for what the pipeline does

