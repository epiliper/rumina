<p align="center">
    <img src="imgs/RUMINA_LOGO_light.png#gh-light-mode-only" width=50% height=50%>
    <img src="imgs/RUMINA_LOGO_dark.png#gh-dark-mode-only" width=50% height=50%>
</p>

## RUMINA - Rust-based Unique Molecular Identifier Network Analysis

RUMINA is a performant pipeline for error correction of Next-Generation Sequencing (NGS) data using Unique Molecular Identifiers (UMIs), fit for shotgun-based and amplicon-based methodologies. 

RUMINA deduplicates reads associated with a single template molecule, using the coordinates of read alignment and UMI barcode sequence to perform correction of PCR and sequencing errors. The above strategy also allows for correction of errors in UMI sequences (directional clustering). 

This pipeline is tested for processing ~600 million reads in ~5 hours, at a rate of 120 million reads processed per hour (tested on 10-core M1 Max Mac Studio).


### Workflow:

<p align="center">
    <img src="imgs/workflow_light.png#gh-light-mode-only" width=100% height=100%>
    <img src="imgs/workflow_dark.png#gh-dark-mode-only" width=100% height=100%>
</p>


### Installation: 

dependencies: 
- [cargo](https://www.rust-lang.org/tools/install) v1.86.0+
#### option 1: install with cargo

```bash
cargo install rumina
```

#### option 2: compile from source

```bash
export RUSTFLAGS="-C target-cpu=native" 

git clone https://github.com/epiliper/rumina.git
cd rumina
cargo build --release 
mv target/release/rumina .
```

The binary will be located at `./rumina`. It's recommended that you move it somewhere to your `$PATH`, so you can run it from anywhere.  
NOTE: Using this option may yield performance gains, as the `target-cpu=native` flag is not used when making the release binaries.

#### option 3: release binaries

Navigate to releases and download the zip for your system's CPU architecture. Unzip and `cd` into the directory, and run `./rumina -h` to ensure it's working.   
It's recommended to move the binary to someplace in your `$PATH` for convenience.


### Usage: 
```rumina <input (file or directory)> -g <grouping_method> -s <separator> <optional args>```

an example command:<br>
```rumina example.bam -g directional -s : -x 100 -t 8``` 

---
The `input` to `rumina` can be a file or a directory; if a directory, all BAM files within (excluding pipeline products) are processed.

RUMINA will write output BAM files and reports to an output directory (`rumina_output` by default), which can be specified with `--outdir`. Output BAM files will be sorted and indexed.

###  Arguments 
:small_blue_diamond: = mandatory, no default

##### `input` :small_blue_diamond:
The input file or directory. If a file, it must be: 

- in BAM format
    - The UMI should be in the read QNAME field (see image under `--separator`). Illumina data base-called by [BCL Convert](https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/bcl-convert.html) should be formatted this way by default.
    - sorted and indexed with `samtools sort` and `samtools index` beforehand.


If the input is a directory, all BAM files within (excluding pipeline products) will be processed per the other arguments specified. 

##### `-g, --grouping_method` :small_blue_diamond:

Specifies how/if to merge UMIs based on edit distance, to account for PCR mutations and NGS errors in UMI sequence. Options are: 
* **directional**: Merge UMIs via directional clustering. See *Amplicon* section for more details. This is the best option for amplicon sequencing data.
* **acyclic**: Same as directional clustering, except UMI networks are limited to a depth of 1, i.e. UMIs that are predicted children cannot have UMIs as child nodes.
* **raw**: Treat each UMI as genuine; UMIs are not merged. This is the best option for metagenomics/shotgun sequencing data.

##### `-s, --separator` :small_blue_diamond:
Specifies the character in the read QNAME delimiting the UMI barcode from the rest of the string. This is usually `_` or `:`.<br>

<p align="center">
    <img src="imgs/barcode_light.png#gh-light-mode-only" width=75% height=75%>
    <img src="imgs/barcode_dark.png#gh-dark-mode-only" width=75% height=75%>
</p>


##### `-x, --split_window` (default = None)
dictates how to split input BAM files into subfiles (for avoiding memory overflow). <br><br> This is usually necessary for BAM files with high sequencing depth that would otherwise cause the program to overuse available memory.

Splitting happens along coordinates of the reference genome in the specified BAM file; If `--split_window 100` was used, reads for every 100bp stretch of the reference would be processed in separate batches, prior to being written to output. This applies to every reference genome present in the input alignment.

Options are: 
* **positive integer**: split input files by a fixed window size. If `input` is a directory, this will be applied to each file within the directory. This has been tested with values ranging from 50 - 500.
* **none** (default): process the input as one file without splitting by coordinate window. Using this option with larger BAMs may result in memory overuse.

##### `-t, --threads` (default = all) 

Specifies the number of threads RUMINA is allowed to use. Threads are used to parallelize processing of individual reference coordinates, and for I/O operations. 

By default, RUMINA will attempt to use all available threads.

##### `-l, --length` (optional)
if used, groups reads by length as well as coordinate. This is recommended for metagenomics data with high read depth, as this will group reads more stringently and likely produce more singleton groups. 

##### `--only-group` (optional)
if used, reads will be grouped (assigned a group-specific "UG" tag), but not deduplicated or error-corrected. This is useful if you want to manually check how grouping works with a given file.

##### `--outdir` (default = rumina_output)
The output directory, relative to the parent directory of the input files/directory, in which RUMINA's output will be stored.

#### Arguments for paired-end input

using either of these arguments will automatically treat the input as paired-end.

##### `--halve_pairs` (optional)

Use only R1 for deduplication, pairing deduplicated R1s with their associated R2s. This is similar to UMI-tools, in that R2 reads are not part of UMI clusters.

##### `--merge_pairs` (REF_FASTA, :small_blue_diamond:)
Use both R1 and R2 for deduplication, and merge overlapping forward/reverse reads with the same barcode after initial deduplication. Merged reads are then realigned to the reference genome, which should be supplied in FASTA format. This is untested with segmented genomes or eukaryotic genomes, and is under active development.

Forward/reverse pairs are merged only if they contain a minimum number of overlapping bases, which is controlled by the `--min_overlap_bp` argument. Forward/reverse pairs identified to have discordant sequences are discarded, and reads unable to be merged for other reasons are still written to output.

##### `--min_overlap_bp` (default = 3)
The minimum number of bases shared by two reads at the same reference coordinates for merging to occur in `--merge_pairs`. Reads not discordant in sequence but not meeting this threshold will not be merged, and instead both be written to the output file.
