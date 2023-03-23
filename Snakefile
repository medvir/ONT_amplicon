SAMPLES = ["1000463191", "1000463386", "1000465755", "1000548083", "1000463192", "1000465563", "1000536513", "100045744501", "100045751801", "100045756501", "100045771201"]

rule all:
    input:
        mapping_index=expand("mapped_reads/{sample}_aln.bam.bai", sample=SAMPLES),
        mapping_depth=expand("mapped_reads/{sample}_aln_depth.tsv", sample=SAMPLES),
        consensus_sequence=expand("consensus_sequences/{sample}_consensus.fasta", sample=SAMPLES),

rule extract_amplicon_from_reads:
    input:
        "data/samples/{sample}_fastq_pass.fastq"
    output:
        "filtered_reads/{sample}_fastq_pass_amplicon.fastq"
    shell:
        "cat {input} | seqkit amplicon -F 'CCAGCACTGACAGCAGYNGARAYNGG' -R 'TACTGGACCACCTGGNGGNAYRWACAT' -r 27:-28 > {output}"

rule mapping_reads:
    input:
        reference="data/references/Enterovirus_all_without_primer.fasta",
        reads="filtered_reads/{sample}_fastq_pass_amplicon.fastq"
    output:
        "mapped_reads/{sample}_aln.bam"
    shell:
        "minimap2 -ax map-ont {input.reference} {input.reads} | samtools sort > {output}"

rule samtools_index:
    input:
        "mapped_reads/{sample}_aln.bam"
    output:
        "mapped_reads/{sample}_aln.bam.bai"
    shell:
        "samtools index {input}"

rule mapping_depth:
    input:
        "mapped_reads/{sample}_aln.bam"
    output:
        depth="mapped_reads/{sample}_aln_depth.tsv",
        coverage="mapped_reads/{sample}_aln_coverage.tsv"
    run:
        shell("samtools depth -a {input} > {output.depth}")
        # get coverage statistics and only keep rows where numreads > 0
        shell("samtools coverage {input} | awk 'NR == 1 || $4 > 0' > {output.coverage}")

rule create_consensus:
    input:
        "mapped_reads/{sample}_aln.bam"
    output:
        "consensus_sequences/{sample}_consensus_raw.fasta"
    run:
        shell("samtools consensus -f fasta -a -l 0 -m simple --ambig --use-qual --het-fract 0.25 --min-depth 30 {input} -o {output}")
        # replace reference accession number with sample id in the fasta header
        shell("sed -i -r 's/^>.*(\\|.*\\|.*$)/>{wildcards.sample}\\1/g' {output}")

rule filter_consensus:
    input:
        "consensus_sequences/{sample}_consensus_raw.fasta"
    output:
        "consensus_sequences/{sample}_consensus.fasta"
    shell:
        # Filter sequence with more than ns_max_p percentage of Ns
        "prinseq -fasta {input} -ns_max_p 5 -line_width 0 -out_bad null -out_good stdout > {output}"
