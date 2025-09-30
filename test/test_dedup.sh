outdir=$1
input_dir="${outdir}/RUMINA_TEST/input"
dedup_out="${outdir}/RUMINA_TEST/ext_dir"
group_out="${outdir}/RUMINA_TEST/ext_dir_group"

mkdir input_dir
mkdir -p dedup_out
mkdir -p group_out

# make input test file visible
cp SRR2057564_sub_ext.bam $input_dir/

echo "Testing: deduplication. Outputting to '${dedup_out}'..."
rumina dedup -i SRR2057564_sub_ext.bam -g directional -s _ --threads 1 --outdir $dedup_out

echo "Testing: only grouping. Outputting to '${group_out}'..."
rumina dedup -i SRR2057564_sub_ext.bam -g directional -s _ --threads 1 --outdir $group_out --only-group
