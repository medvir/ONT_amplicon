SAMPLES = ["1000463191", "1000463386", "1000465755", "1000548083", "1000463192", "1000465563", "1000536513"]

rule all:
    input:
        mapping_index=expand("mapped_reads/{sample}_aln.bam.bai", sample=SAMPLES),
        mapping_depth=expand("mapped_reads/{sample}_aln_depth.tsv", sample=SAMPLES),
        consensus_sequence=expand("consensus_sequences/{sample}_consensus.fasta", sample=SAMPLES),
        target_consensus_sequence=expand("consensus_sequences_target/{sample}_consensus.fasta", sample=SAMPLES)

rule human_genome_indexing:
    input:
        "data/GRCh38/GRCh38_latest_genomic.fna"
    output:
        "data/GRCh38/GRCh38_latest_genomic.mmi"
    shell:
        "minimap2 -d {output} {input}"

rule human_genome_decontamination:
    input:
        human_genome_index="data/GRCh38/GRCh38_latest_genomic.mmi",
        sequencing_reads="data/samples/{sample}_fastq_pass.fastq"
    output:
        "filtered_reads/{sample}_fastq_pass_decont.fastq"
    shell:
        "minimap2 -ax map-ont {input.human_genome_index} {input.sequencing_reads} | samtools view -f 4 | samtools fastq > {output}"

rule quality_filtering:
    input:
        "filtered_reads/{sample}_fastq_pass_decont.fastq"
    output:
        "filtered_reads/{sample}_fastq_pass_decont_filtered.fastq"
    shell:
        "prinseq -fastq {input} -min_len 250 -max_len 450 -min_qual_mean 15 -lc_method entropy -lc_threshold 80 -out_bad null -out_good filtered_reads/{wildcards.sample}_fastq_pass_decont_filtered"

rule extract_amplicon_from_reads:
    input:
        "filtered_reads/{sample}_fastq_pass_decont_filtered.fastq"
    output:
        "filtered_reads/{sample}_fastq_pass_decont_filtered_amplicon.fastq"
    shell:
        "cat {input} | seqkit amplicon -F 'CCAGCACTGACAGCAGYNGARAYNGG' -R 'TACTGGACCACCTGGNGGNAYRWACAT' -r 27:-28 > {output}"

rule denovo_assembly:
    input:
        "filtered_reads/{sample}_fastq_pass_decont_filtered_amplicon.fastq"
    output:
        assembly_folder=directory("assembled_reads/{sample}_megahit"),
        assembly="assembled_reads/{sample}_contigs.fasta"
    run:
        shell("megahit -r {input} --min-count 3 --out-dir {output.assembly_folder}")
        shell("cp assembled_reads/{wildcards.sample}_megahit/final.contigs.fa {output.assembly}")

rule mapping_reads:
    input:
        contigs="assembled_reads/{sample}_contigs.fasta",
        reads="filtered_reads/{sample}_fastq_pass_decont_filtered.fastq"
    output:
        "mapped_reads/{sample}_aln.bam"
    shell:
        "minimap2 -ax map-ont {input.contigs} {input.reads} | samtools sort > {output}"

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
        "mapped_reads/{sample}_aln_depth.tsv"
    shell:
        "samtools depth -aa {input} > {output}"

rule create_consensus:
    input:
        "mapped_reads/{sample}_aln.bam"
    output:
        "consensus_sequences/{sample}_consensus.fasta"
    run:
        shell("samtools consensus -f fasta -m simple --ambig --use-qual --het-fract 0.25 --min-depth 30 {input} -o {output}")
        # replace k-mer number with sample id in the fasta header
        shell("sed -i 's/>.*_/>{wildcards.sample}_/g' {output}")

rule create_target:
    input:
        reference="data/target/{sample}_without_primer.fasta",
        reads="filtered_reads/{sample}_fastq_pass_decont_filtered.fastq"
    output:
        bam="mapped_reads_target/{sample}_aln.bam",
        bai="mapped_reads_target/{sample}_aln.bam.bai",
        depth="mapped_reads_target/{sample}_aln_depth.tsv",
        consensus="consensus_sequences_target/{sample}_consensus.fasta"
    run:
        shell("minimap2 -ax map-ont {input.reference} {input.reads} | samtools sort > {output.bam}")
        shell("samtools index {output.bam}")
        shell("samtools depth -aa {output.bam} > {output.depth}")
        shell("samtools consensus -f fasta -m simple --ambig --use-qual --het-fract 0.25 --min-depth 30 {output.bam} -o {output.consensus}")