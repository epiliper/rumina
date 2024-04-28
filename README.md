## RUMINA - Rust-based Unique Molecular Identifier Network Analysis

- removing singleton UMIs (potential artifacts)
- deduplicating UMIs via UMI_tools 

#### to run: 
<!-- `python3.11 <parent dir of uclean.py>/uclean.py <path_of_input.bam>` -->
`python3 driver.py <input (file or directory)> <optional args>`

The `input` to `driver.py` can be a file or a directory; in the latter case, all .bam files within a directory (exlcuding pipeline products) are processed. This is not yet parallelized. 

If `--report_coverage` is used (see optional arguments) with a directory input, the reports for each bamfile will be summarized into a final .csv file (`(DEDUP)_COVERAGE_REPORT.csv`). 

#### requirements: 
- python3+
- pysam v0.22.0+ 
- umi_tools v1.1.4+

#### Input types 

- as of 2024-Mar-27, this pipeline only works on .BAM files. .SAM files are not currently supported. 

#### Optional arguments
| arg | explanation | 
| --- | --- | 
| `--delete_temps` | deletes pipeline-generated files needed temporarily for processing.<br>Can save gigabytes of space when working with large files |
| `--report_coverage` | generates per-sample and combined .csv reports on both *input* and *output* files using `bedtools genomecov`.<br>Doing this with large bam files can increase the runtime by several minutes per file | 

#### Todo
- Attempt to include mapping quality stats in `--report_coverage`
- investigate faster ways of grouping/deduping bams
- parallelize processing multiple bamfiles (contingent on above task)
