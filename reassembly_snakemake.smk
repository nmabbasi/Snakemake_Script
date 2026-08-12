REFERENCE = "21.fasta"
TRUSTED_CONTIGS = "2_freshwater.fasta"
SAMPLES = ["2_fw_Tous12m", "2_fw_Tous25m"]

READS = {
    "2_fw_Tous12m": ("input/P1.renamed.fasta", "input/P2.renamed.fasta"),
    "2_fw_Tous25m": ("input/P1.renamed.fasta", "input/P2.renamed.fasta"),
}

BWA = "bwa"
SPADES = "spades.py"
EXTRACT_SCRIPT = "extract.fasta.from.sam.using.list.pl"


rule all:
    input:
        "assembly_try1/contigs.fasta"


rule index:
    input:
        REFERENCE
    output:
        REFERENCE + ".bwt",
        REFERENCE + ".pac",
        REFERENCE + ".amb",
        REFERENCE + ".ann",
        REFERENCE + ".sa"
    shell:
        "{BWA} index {input}"


rule aln:
    input:
        ref=REFERENCE,
        p1=lambda wc: READS[wc.sample][0],
        p2=lambda wc: READS[wc.sample][1]
    output:
        p1=REFERENCE + ".{sample}.P1.sai",
        p2=REFERENCE + ".{sample}.P2.sai"
    threads:
        30
    shell:
        """
        {BWA} aln -t {threads} {input.ref} {input.p1} > {output.p1}
        {BWA} aln -t {threads} {input.ref} {input.p2} > {output.p2}
        """


rule sampe:
    input:
        ref=REFERENCE,
        p1_sai=REFERENCE + ".{sample}.P1.sai",
        p2_sai=REFERENCE + ".{sample}.P2.sai",
        p1=lambda wc: READS[wc.sample][0],
        p2=lambda wc: READS[wc.sample][1]
    output:
        "21_m_{sample}.sam"
    shell:
        "{BWA} sampe {input.ref} {input.p1_sai} {input.p2_sai} {input.p1} {input.p2} > {output}"


rule extract_fasta:
    input:
        sam="21_m_{sample}.sam",
        ids="input/{sample}_list.txt"
    output:
        "{sample}_final.fasta"
    shell:
        "perl {EXTRACT_SCRIPT} -s {input.sam} -l {input.ids} -o {output}"


rule combine_fasta:
    input:
        expand("{sample}_final.fasta", sample=SAMPLES)
    output:
        "all_reads.fasta"
    shell:
        "cat {input} > {output}"


rule spades:
    input:
        reads="all_reads.fasta",
        trusted=TRUSTED_CONTIGS
    output:
        "assembly_try1/contigs.fasta"
    threads:
        20
    params:
        memory=500
    shell:
        "{SPADES} -o assembly_try1 --12 {input.reads} --trusted-contigs {input.trusted} --only-assembler --careful -t {threads} -m {params.memory}"
