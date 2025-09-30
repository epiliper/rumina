echo "Testing: deduplication..."
rumina dedup -i SRR2057564_sub_ext.bam -g directional -s _ --threads 1 --outdir ext_dir

echo "Testing: only grouping..."
rumina dedup -i SRR2057564_sub_ext.bam -g directional -s _ --threads 1 --outdir ext_dir_group --only-group
