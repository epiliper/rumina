<p align="center">
    <img src="imgs/RUMINA_LOGO_light.png#gh-light-mode-only" width=50% height=50%>
    <img src="imgs/RUMINA_LOGO_dark.png#gh-dark-mode-only" width=50% height=50%>
</p>

## RUMINA - Rust-based Unique Molecular Identifier Network Analysis

#### to run: 
`python3 main.py <input (file or directory)> <optional args>`

The `input` to `main.py` can be a file or a directory; in the latter case, all .bam files within a directory (excluding pipeline products) are processed (including those that require splitting).

If `--report_coverage` is used (see optional arguments) with a directory input, the reports for each bamfile will be summarized into a final .csv file (`COVERAGE_COMPLETE.csv`). 

#### requirements: 
- python 3
    - pandas
    - pysam v0.22.0+
- bedtools v2.31.1+

#### Input types 

- as of 2024-Mar-27, this pipeline only works on .BAM files. 

#### Optional arguments
| arg | explanation | 
| --- | --- | 
| `--delete_temps` | deletes pipeline-generated files needed temporarily for processing.<br>Can save gigabytes of space when working with large files |
| `--report_coverage` | generates per-sample and combined .csv reports on both *input* and *output* files using `bedtools genomecov`.<br>Doing this with large bam files can increase the runtime by several minutes per file | 
| `--split_window` | specifies how to split input bam files. Options: <br> - `auto` (default): calculate recommended split size based on size of input file <br> - any positive integer: arbitrary window size |

