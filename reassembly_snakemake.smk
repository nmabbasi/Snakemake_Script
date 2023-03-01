# Index reference genome in fasta format
rule index:
    input:
        "21.fasta"
    output:
        "21.fasta.bwt",
        "21.fasta.pac",
        "21.fasta.amb",
        "21.fasta.ann",
        "21.fasta.sa"
    shell:
        "/home/rohit/josemanuel/programs/bwa-0.6.2/bwa index {input}"

# Align reads to indexed reference genome using BWA
rule aln:
    input:
        "21.fasta",
        "input/P1.renamed.fasta",
        "input/P2.renamed.fasta"
    output:
        "21.fasta.{sample}.P1.sai",
        "21.fasta.{sample}.P2.sai"
    params:
        threads = 30
    shell:
        "/home/rohit/josemanuel/programs/bwa-0.6.2/bwa aln -t {params.threads} {input[0]} {input[1]} > {output[0]}",
        "/home/rohit/josemanuel/programs/bwa-0.6.2/bwa aln -t {params.threads} {input[0]} {input[2]} > {output[1]}"

# Generate SAM file from aligned reads
rule sampe:
    input:
        "21.fasta",
        "21.fasta.{sample}.P1.sai",
        "21.fasta.{sample}.P2.sai",
        "input/P1.renamed.fasta",
        "input/P2.renamed.fasta"
    output:
        "21_m_{sample}.sam"
    shell:
        "/home/rohit/josemanuel/programs/bwa-0.6.2/bwa sampe {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} > {output}"

# Extract sequences from SAM file based on list of read IDs
rule extract_fasta:
    input:
        "21_m_{sample}.sam",
        "input/{sample}_list.txt"
    output:
        "{sample}_final.fasta"
    shell:
        "perl extract.fasta.from.sam.using.list.pl -s {input[0]} -l {input[1]} -o {output}"

# Combine extracted sequences into one fasta file
rule combine_fasta:
    input:
        "2_fw_Tous12m_final.fasta",
        "2_fw_Tous25m_final.fasta"
    output:
        "all_reads.fasta"
    shell:
        "cat {input} > {output}"

# Assemble contigs from combined fasta file using SPAdes
rule spades:
    input:
        "all_reads.fasta",
        "2_freshwater.fasta"
    output:
        "assembly_try1/contigs.fasta"
    params:
        threads = 20,
        memory = 500
    shell:
        "/home/rohit/josemanuel/programs/SPAdes-3.6.1-Linux/bin/spades.py -o assembly_try1 --12 {input[0]} --trusted-contigs {input[1]} --only-assembler --careful -t {params.threads} -m {params.memory}"



