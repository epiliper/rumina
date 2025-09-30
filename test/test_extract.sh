outdir=$1
outdir="${outdir}/RUMINA_TEST"
input_dir="${outdir}/input"

mkdir -p $outdir
mkdir -p $input_dir

# make input test file visible
cp SRR2057564_sub.fastq.gz $input_dir/

echo "Testing: extract. Outputting to ${outdir}..."
rumina extract -i SRR2057564_sub.fastq.gz \
  -p NNNXXXXNN -o $outdir/SRR2057564_sub_ext.fastq.gz
