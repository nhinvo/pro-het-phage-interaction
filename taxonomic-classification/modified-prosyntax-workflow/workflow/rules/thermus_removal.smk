rule index_genome:
    """
    Index reference Thermus genome. 
    """
    input:
        thermus_genome = config["input"]["thermus_genome"], 
    output:
        touch(scratch_dict["genome_index_done"]),
    conda:
        "../envs/thermus_removal.yaml"
    shell:
        "bowtie2-build --threads {resources.tasks} {input.thermus_genome} {input.thermus_genome}"

        
rule obtain_thermus_reads:
    """
    Map reads to Thermus genome and obtain names of mapped reads. 
    """
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq.gz",
        thermus_genome = config["input"]["thermus_genome"], 
        indexing = scratch_dict["genome_index_done"]
    output:
        mapped_sam = temp(scratch_dict["mapped_thermus"] / "{sample}_mapped.sam"),
        mapped_bam = temp(scratch_dict["mapped_thermus"] / "{sample}_mapped.bam"),
        thermus_read_name = temp(scratch_dict["mapped_thermus"] / "{sample}_thermus_reads.txt"), 
    conda:
        "../envs/thermus_removal.yaml"  
    shell:
        """
        # map reads to Thermus genome 
        bowtie2 -x {input.thermus_genome} -1 {input.r1} -2 {input.r2} -S {output.mapped_sam} -p {resources.tasks}

        # convert to bam 
        samtools view -S -b -F 4 {output.mapped_sam} > {output.mapped_bam}

        # extract read name and remove duplicate names 
        samtools view {output.mapped_bam} | cut -f1 | sort | uniq > {output.thermus_read_name}
        """


rule remove_thermus_reads:
    """
    Remove reads mapped to Thermus genome. 
    
    -v: invert match (get reads not in thermus list)
    -i: ignore case 
    -f: read name file 
    """
    input: 
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq.gz",
        thermus_read_name = scratch_dict["mapped_thermus"] / "{sample}_thermus_reads.txt", 
    output:
        o1 = temp(scratch_dict["thermus_removed_reads"] / "{sample}_1_trimmed_no_thermus.fastq"), 
        o2 = temp(scratch_dict["thermus_removed_reads"] / "{sample}_2_trimmed_no_thermus.fastq"), 
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit grep -v -i -f {input.thermus_read_name} {input.r1} > {output.o1}
        seqkit grep -v -i -f {input.thermus_read_name} {input.r2}> {output.o2}
        """
