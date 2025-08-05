#!/bin/bash
#SBATCH --job-name WV
#SBATCH -p sched_mit_chisholm
#SBATCH -n 1
#SBATCH --time=5-0
#SBATCH --mem=1G
#SBATCH -o logs/snakemake.%j.out
#SBATCH -e logs/snakemake.%j.err

eval "$(conda shell.bash hook)"  
conda activate widevariant-snakemake

echo "$(date): Begnning WV Snakemake pipeline..."
snakemake -v 
snakemake --profile .

echo "$(date): Done!"
