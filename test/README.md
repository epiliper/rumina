# Testing RUMINA with sub-sampled iCLIP dataset

**NOTE:** The first part of this test (going from raw FASTQ to BAM) uses the [mm9 reference genome](https://genome.ucsc.edu/cgi-bin/hgGateway?db=mm9) and [minimap2 2.29-r1283](https://github.com/lh3/minimap2). If you don't wish to perform these steps, the output BAM (`SRR2057564_sub_ext.bam`) is already provided; skip to step 4 to avoid dependencies.

---
Test file: SRR2057564 sub-sampled to the first 12000 reads. The sub-sampled file (SRR2057564_sub.fastq.gz) is already provided.

0. If you don't wish to redownload the original SRA FASTQ, then clone this repository and use the included FASTQ and BAM files, skipping to step 2.

1. Sub-sampling was performed below:
```bash
gzcat SRR2057564.fastq.gz | head -n 48000 | pigz > SRR2057564_sub.fastq.gz
```

2. We next want to extract barcodes from the reads according to the layout specified in the [UMI-tools publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC5340976/)
```bash
rumina extract -i SRR2057564_sub.fastq.gz \
  -p NNNXXXXNN -o SRR2057564_sub_ext.fastq.gz
```

3. We next want to map to the mm9 reference genome to generate a BAM:
```bash
minimap2 -t 10 -L --secondary no -ax splice:sr \
  mm9.fa.gz SRR2057564_sub_ext.fastq.gz | samtools view -bS -F 2052 > SRR2057564_sub_ext

# remember to sort and index; RUMINA needs that
samtools sort SRR2057564_sub_ext > SRR2057564_sub_ext.bam
rm SRR2057564_sub_ext
samtools index SRR2057564_sub_ext.bam
```

4. We now have everything we need to run RUMINA on the file: barcodes are in read headers, and reads are mapped to the appropriate reference. Now let's run RUMINA: 
```bash
rumina dedup -i SRR2057564_sub_ext.bam -g directional -s _ --threads 1 --outdir ext_dir
```

Here, we are telling RUMINA to deduplicate the BAM file with the directional method, detect the UMI as the character sequence immediately after the last underscore ('_') in the header, use 1 thread, and output results to the `./ext_dir` directory.

You should find an output, deduplicated BAM in `./ext_dir/SRR2057564_sub_ext_RUMINA.bam`. An associated deduplication report can be found in the same directory. RUMINA also outputs `./rumina_group.log`; this reports detailed deduplication metrics and shouldn't be considered unless troubleshooting is required.

5. Alternatively, if you just want to see UMI clustering without deduplication, you can run the following:
```bash
rumina dedup -i SRR2057564_sub_ext.bam -g directional -s _ --threads 1 --outdir ext_dir_group --only-group
```

Notice the `--only-group` at the end. This will output to `./ext_dir_group/SRR2057564_sub_ext_RUMINA.bam`), but will only mark reads with the UMI cluster ID they belong to, not deduplicate them. Cluster ID can be found in the UG tag. 
