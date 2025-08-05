rule SRA_download:
    """
    Download raw read sample from SRA. 
    """
    output:
        r1 = temp(scratch_dict["SRA_dl"] / "{sample}_1.fastq"),
        r2 = temp(scratch_dict["SRA_dl"] / "{sample}_2.fastq"),
        sra = temp(scratch_dict["SRA_dl"] / "{sample}" / "{sample}.sra"),
    conda:
        "../envs/sra-tools.yaml"
    shell:
        """
        ACCESSION={wildcards.sample}
        OUTPUT_DIR=$(dirname {output.r1})

        echo prefetch ... 
        prefetch \
            $ACCESSION \
            --max-size 10000G \
            --force ALL \
            --verbose \
            --output-directory $OUTPUT_DIR

        echo fasterq-dump ... 
        fasterq-dump \
            $OUTPUT_DIR/$ACCESSION/ \
            --split-files \
            --force --temp $OUTPUT_DIR \
            --verbose \
            --threads {resources.cpus_per_task} \
            --outdir $OUTPUT_DIR
        """

