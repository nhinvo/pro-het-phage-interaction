# specifies with run trim rule to run first 
ruleorder: run_trim_PE_local > run_trim_PE_sra 


rule run_trim_PE_local:
    """
    Trims reads for local files that have forward and reverse read files. 
    """
    input:
        r1 = lambda wildcards: str(SAMPLE_TABLE.loc[wildcards.sample, 'forward read']),
        r2 = lambda wildcards: str(SAMPLE_TABLE.loc[wildcards.sample, 'reverse read']),
        ref = Path(config["input"]["adapter_file"]),
    output:
        o1 = temp(scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz"),
        o2 = temp(scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq.gz"),
    conda:
        "../envs/bbtools.yaml"
    shell:
        "bbduk.sh threads={resources.cpus_per_task} "
        "in1={input.r1} in2={input.r2} "
        "out1={output.o1} out2={output.o2} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1"


rule run_trim_PE_sra:
    """
    Trims reads for files that need SRA download first. 
    """
    input:
        r1 = scratch_dict["SRA_dl"] / "{sample}_1.fastq",
        r2 = scratch_dict["SRA_dl"] / "{sample}_2.fastq",
        ref = Path(config["input"]["adapter_file"]),
    output:
        o1 = temp(scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz"),
        o2 = temp(scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq.gz"),
    conda:
        "../envs/bbtools.yaml"
    shell:
        "bbduk.sh threads={resources.cpus_per_task} "
        "in1={input.r1} in2={input.r2} "
        "out1={output.o1} out2={output.o2} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1"
