INPUT_FASTQ = "C.C.fastq.gz"
OUTPUT_PREFIX = "asm_CC.hifiasm.asm"
PRIMARY_GFA = OUTPUT_PREFIX + ".p_ctg.gfa"
PRIMARY_CONTIGS = OUTPUT_PREFIX + ".p_ctgl2.fasta"
BUSCO_OUTDIR = "output"
BLASTN_OUT = "blastn.tsv"
BLOBDIR = "tempo"
BLOBDB = BLOBDIR + "/dataset_0.blobDB.json"

TAXDUMP = "../../../../3-DATABASES/taxdump"
BUSCO_LINEAGE = "busco_downloads/lineages/insecta_odb10"
BUSCO_TABLE = BUSCO_OUTDIR + "/run_insecta_odb10/full_table.tsv"
READ_TO_CONTIG_BAM = "read_to_contig.bam"


rule all:
    input:
        BLOBDB


rule fastqc:
    input:
        INPUT_FASTQ
    output:
        html="fastqc_report.html"
    shell:
        "fastqc {input} --outdir ."


rule hifiasm:
    input:
        INPUT_FASTQ
    output:
        OUTPUT_PREFIX + ".p_ctg.gfa"
    threads:
        32
    shell:
        "hifiasm -l2 --primary -o {OUTPUT_PREFIX} -t {threads} {input}"


rule contigs:
    input:
        PRIMARY_GFA
    output:
        PRIMARY_CONTIGS
    shell:
        r"""awk '/^S/ {{print ">"$2"\n"$3}}' {input} | fold > {output}"""


rule busco:
    input:
        PRIMARY_CONTIGS
    output:
        directory(BUSCO_OUTDIR)
    threads:
        16
    shell:
        "busco -i {input} -o {output} -l {BUSCO_LINEAGE} -m genome --cpu {threads}"


rule blastn:
    input:
        PRIMARY_CONTIGS
    output:
        BLASTN_OUT
    threads:
        16
    log:
        "log_blastn.txt"
    shell:
        "blastn -query {input} -db nt -outfmt '6 qseqid staxids bitscore std' "
        "-max_target_seqs 10 -max_hsps 1 -num_threads {threads} "
        "-evalue 1e-25 -out {output} > {log} 2>&1"


rule blobtools:
    input:
        fasta=PRIMARY_CONTIGS,
        hits=BLASTN_OUT,
        busco=BUSCO_TABLE,
        cov=READ_TO_CONTIG_BAM
    output:
        BLOBDB
    threads:
        16
    log:
        "log_blob.txt"
    shell:
        "mkdir -p {BLOBDIR} && "
        "blobtools create --fasta {input.fasta} --cov {input.cov} "
        "--hits {input.hits} --taxdump {TAXDUMP} --busco {input.busco} "
        "--threads {threads} {BLOBDIR} > {log} 2>&1"
