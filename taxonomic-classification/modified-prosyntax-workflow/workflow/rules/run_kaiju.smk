rule kaiju_run:
    """
    Run kaiju on Thermus-removed read files.
    """
    input:
        # thermus-removed read files 
        r1 = scratch_dict["thermus_removed_reads"] / "{sample}_1_trimmed_no_thermus.fastq", 
        r2 = scratch_dict["thermus_removed_reads"] / "{sample}_2_trimmed_no_thermus.fastq", 
        # kaiju files 
        nodes = Path(config["input"]["nodes_file"]),
        fmi = Path(config["input"]["fmi_file"]),
    output:
        temp(scratch_dict["base_kaiju"] / "{sample}_kaiju.txt"), 
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju \
            -z {resources.cpus_per_task} \
            -m 11 -s 65 -E 0.05 -x \
            -e 5 -t {input.nodes} -f {input.fmi} \
            -i {input.r1} {input.r2} \
            -o {output}
        """


rule kaiju_name:
    """
    Returns [sample]_name.out file that contains cols: 
        - [classification_status, read_name, taxonid, full_taxa]
    
    Notes: 
        -p: print the full taxon path instead of just the taxon name.
        -u: omit unclassified reads (saves space).
    """
    input:
        kaiju = scratch_dict["base_kaiju"] / "{sample}_kaiju.txt",
        # kaiju files
        nodes = Path(config["input"]["nodes_file"]),
        names = Path(config["input"]["names_file"]),
    output:
        scratch_dict["kaiju_names"] / "{sample}_names.out",      
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju-addTaxonNames \
            -t {input.nodes} \
            -n {input.names} \
            -i {input.kaiju} \
            -p -u -o {output}
        """

rule kaiju_summary_genus:
    input:
        kaiju = scratch_dict["base_kaiju"] / "{sample}_kaiju.txt",
        # kaiju files
        nodes = Path(config["input"]["nodes_file"]),
        names = Path(config["input"]["names_file"]),
    output:
        scratch_dict["kaiju_summary"] / "{sample}_kaiju_summary_genus.tsv", 
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju2table \
            -t {input.nodes} \
            -n {input.names} \
            -r genus \
            {input.kaiju} \
            -o {output}
        """

rule kaiju_summary_family:
    input:
        kaiju = scratch_dict["base_kaiju"] / "{sample}_kaiju.txt",
        # kaiju files
        nodes = Path(config["input"]["nodes_file"]),
        names = Path(config["input"]["names_file"]),
    output:
        scratch_dict["kaiju_summary"] / "{sample}_kaiju_summary_family.tsv", 
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju2table \
            -t {input.nodes} \
            -n {input.names} \
            -r family \
            {input.kaiju} \
            -o {output}
        """

rule kaiju_summary_order:
    input:
        kaiju = scratch_dict["base_kaiju"] / "{sample}_kaiju.txt",
        # kaiju files
        nodes = Path(config["input"]["nodes_file"]),
        names = Path(config["input"]["names_file"]),
    output:
        scratch_dict["kaiju_summary"] / "{sample}_kaiju_summary_order.tsv", 
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju2table \
            -t {input.nodes} \
            -n {input.names} \
            -r order \
            {input.kaiju} \
            -o {output}
        """

rule kaiju_summary_class:
    input:
        kaiju = scratch_dict["base_kaiju"] / "{sample}_kaiju.txt",
        # kaiju files
        nodes = Path(config["input"]["nodes_file"]),
        names = Path(config["input"]["names_file"]),
    output:
        scratch_dict["kaiju_summary"] / "{sample}_kaiju_summary_class.tsv", 
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju2table \
            -t {input.nodes} \
            -n {input.names} \
            -r class \
            {input.kaiju} \
            -o {output}
        """

rule kaiju_summary_phylum:
    input:
        kaiju = scratch_dict["base_kaiju"] / "{sample}_kaiju.txt",
        # kaiju files
        nodes = Path(config["input"]["nodes_file"]),
        names = Path(config["input"]["names_file"]),
    output:
        scratch_dict["kaiju_summary"] / "{sample}_kaiju_summary_phylum.tsv", 
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju2table \
            -t {input.nodes} \
            -n {input.names} \
            -r phylum \
            {input.kaiju} \
            -o {output}
        """