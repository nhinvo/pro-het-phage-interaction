#!/bin/bash
#SBATCH --job-name bcftools-view
#SBATCH -p sched_mit_tami
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH -o logs/bcftools.%j.out
#SBATCH -e logs/bcftools.%j.err

eval "$(conda shell.bash hook)"  
conda activate bcftools-1.22

VCF_DIR=../../modified-widevariant-workflow/1-run-WV/1-Mapping/vcf

mkdir -p data/1.bcftools_out


# PROCESS VCF # 
for fpath in ${VCF_DIR}/*strain.vcf.gz; do
    # obtain file name 
    fname=$(basename "$fpath" | cut -d'.' -f1)
    echo "$(date): Processing file: ${fname}"

    # run bcftools to obtain all coords (after quality filtering by bcftools view)
    bcftools view \
        --include 'QUAL >= 30 && ALT != "." && INFO/FQ < -30' \
        ${fpath} | \
    bcftools query \
        --format '%CHROM\t%POS\t%REF\t%ALT\t%INFO/INDEL\n' \
        --output "data/1.bcftools_out/${fname}_strain.tsv"

    echo "$(date): Processed file: ${fname}"
    echo 
done


echo "$(date): Completed all bcftools runs"