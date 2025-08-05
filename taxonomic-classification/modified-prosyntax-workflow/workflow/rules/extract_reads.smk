rule kaiju_name_extract:
    """
    Extract names of reads classified as Pro/Syn and obtain their full taxonomic classification. 
    """
    input:
        kaiju_name = scratch_dict["kaiju_names"] / "{sample}_names.out",
    output:
        read_name_file = temp(scratch_dict["prosyn_reads"]["read_name"] / "{sample}.txt"), 
        read_name_taxa_file = temp(scratch_dict["prosyn_reads"]["read_name_classification"] / "{sample}_classification.txt"), 
    shell:
        """
        # obtain name of reads whose classification contains Pro/Syn 
        grep -E "Prochlorococcus|Synechococcus" {input.kaiju_name} | cut -f2 > {output.read_name_file}

        # obtain name and full taxonomic classification of reads whose classification contains Pro/Syn 
        grep -E "Prochlorococcus|Synechococcus" {input.kaiju_name} | cut -f2,4 > {output.read_name_taxa_file}
        """


rule extract_fastq_reads:
    """
    Extract fastq of reads classified as Pro/Syn.
    """
    input: 
        # r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",  # use the thermus-removed file instead 
        r1 = scratch_dict["thermus_removed_reads"] / "{sample}_1_trimmed_no_thermus.fastq", 
        prosyn_read_name = scratch_dict["prosyn_reads"]["read_name"] / "{sample}.txt", 
    output:
        fwd_prosyn_reads = temp(scratch_dict["prosyn_reads"]["extracted_reads"] / "{sample}_fwd.fastq"), 
    conda:
        "../envs/seqtk.yaml"
    shell:
        """
        seqtk subseq {input.r1} {input.prosyn_read_name} > {output.fwd_prosyn_reads}
        """


rule extracted_fastq_to_fasta: 
    """
    Convert Pro/Syn fastq to fasta. 
    """
    input: scratch_dict["prosyn_reads"]["extracted_reads"] / "{sample}_fwd.fastq", 
    output: temp(scratch_dict["prosyn_reads"]["extracted_reads"] / "{sample}_fwd.fasta"), 
    conda: "../envs/seqtk.yaml"
    shell: "seqtk seq -A {input} > {output}"

