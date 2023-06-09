import pandas as pd

samples = pd.read_table("config/samples.tsv", dtype = str).set_index("sample_name", drop=False)

def get_run_name(wildcards):
    return samples.loc[wildcards.sample]["run_name"]

def get_flowcell_name(wildcards):
    return samples.loc[wildcards.sample]["flowcell_name"]

def get_barcode_name(wildcards):
    return samples.loc[wildcards.sample]["barcode_name"]

configfile: "config/config.yml"

rule all:
    input:
        mapping_index=expand("mapped_reads/{sample.sample_name}_aln.bam.bai", sample=samples.itertuples()),
        mapping_depth=expand("mapped_reads/{sample.sample_name}_aln_depth.tsv", sample=samples.itertuples()),
        apping_coverage=expand("mapped_reads/{sample.sample_name}_aln_coverage.tsv", sample=samples.itertuples()),
        consensus_sequence=expand("consensus_sequences/{sample.sample_name}_consensus.fasta", sample=samples.itertuples()),

rule combine_reads:
    params:
        run_name=get_run_name,
        flowcell_name=get_flowcell_name,
        barcode_name=get_barcode_name
    output:
        "reads_raw/{sample}_fastq_pass.fastq"
    shell:
        "zcat /data/GridION/GridIONOutput/{params.run_name}/{params.flowcell_name}/*/fastq_pass/{params.barcode_name}/*.fastq.gz > {output}"

rule extract_amplicon:
    input:
        "reads_raw/{sample}_fastq_pass.fastq"
    params:
        primer_forward=config["primer_forward"],
        primer_reverse=config["primer_reverse"],
        amplicon_start=len(config["primer_forward"])+1,
        amplicon_end=-(len(config["primer_reverse"])+1)
    output:
        "reads_filtered/{sample}_fastq_pass_amplicon.fastq"
    shell:
        # extract amplicon by perfect primer match and remove primer parts
        "cat {input} | seqkit amplicon -F {params.primer_forward} -R {params.primer_reverse} -r {params.amplicon_start}:{params.amplicon_end} > {output}"

rule mapping_reads:
    input:
        reference=config["reference"],
        reads="reads_filtered/{sample}_fastq_pass_amplicon.fastq"
    output:
        "mapped_reads/{sample}_aln.bam"
    shell:
        "minimap2 -ax map-ont {input.reference} {input.reads} | samtools sort > {output}"

rule index_alignment:
    input:
        "mapped_reads/{sample}_aln.bam"
    output:
        "mapped_reads/{sample}_aln.bam.bai"
    shell:
        "samtools index {input}"

rule compute_coverage:
    input:
        "mapped_reads/{sample}_aln.bam"
    output:
        depth="mapped_reads/{sample}_aln_depth.tsv",
        coverage="mapped_reads/{sample}_aln_coverage.tsv"
    shell:
        """
        samtools depth -a {input} > {output.depth}

        # get coverage statistics and only keep rows where numreads > 0
        samtools coverage {input} | awk 'NR == 1 || $4 > 0' > {output.coverage}
        """

rule create_consensus:
    input:
        "mapped_reads/{sample}_aln.bam"
    params:
        min_depth_factor=config["min_depth_factor"],
        min_depth_reads=config["min_depth_reads"],
        het_fract=config["het_fract"]
    output:
        "consensus_sequences/{sample}_consensus_raw.fasta"
    shell:
        """
        mapped_reads=$(samtools view -c -F 4 {input})
        min_depth=$(( mapped_reads / {params.min_depth_factor} ))
        # choose whatever number is larger for the actual min_depth
        min_depth=$(( min_depth > {params.min_depth_reads} ? min_depth : {params.min_depth_reads} ))

        samtools consensus -f fasta -a -l 0 -m simple --ambig --use-qual --het-fract {params.het_fract} --min-depth $min_depth {input} -o {output}

        # replace reference accession number with sample id in the fasta header
        sed -i -r 's/^>.*(\\|.*\\|.*$)/>{wildcards.sample}\\1/g' {output}
        """

rule filter_consensus:
    input:
        "consensus_sequences/{sample}_consensus_raw.fasta"
    params:
        ns_max_p=config["ns_max_p"]
    output:
        "consensus_sequences/{sample}_consensus.fasta"
    shell:
        # Filter sequence with more than ns_max_p percentage of Ns and remove exact duplicates
        # use the stdout option to also write a (empty) file even if no sequence remains 
        "prinseq -fasta {input} -ns_max_p {params.ns_max_p} -derep 1 -line_width 0 -out_bad null -out_good stdout > {output}"
