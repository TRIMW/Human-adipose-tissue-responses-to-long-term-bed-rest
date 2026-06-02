for fq in *.fastq; do
    echo submitting job for $fq
    sbatch -p short -t 0-0:40:0 --mem 1G --wrap "gzip -f $fq"
done
