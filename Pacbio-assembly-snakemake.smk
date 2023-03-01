import os

# Define input and output files
input_fastq = "C.C.fastq.gz"
output_asm = "asm_CC.hifiasm.asm"
output_contigs = "asm_CC.hifiasm.asm.p_ctgl2.fasta"
output_busco = "output"
output_blastn = "blastn.tsv"
output_blobtools = "tempo/dataset_0.blobDB.json"

# Define the number of threads to use
threads = 92

# Define the location of the taxdump folder
taxdump = "../../../../3-DATABASES/taxdump"

# Define the lineage file for BUSCO analysis
busco_lineage = "/home/biopatic/abbasi/P1--Cotesia-Pacbio/2-Canu-Assembly/CC/CC_hifi/3-Busco_analysis/busco_downloads/lineages/insecta_odb10"

# Define the rule to run FastQC on the input fastq file
rule fastqc:
    input:
        input_fastq
    output:
        "fastqc_report.html"
    shell:
        "fastqc {input} --outdir ."

# Define the rule to run hifiasm on the input fastq file
rule hifiasm:
    input:
        input_fastq
    output:
        output_asm
    shell:
        "hifiasm -l2 --primary -o {output} -t {threads} {input}"

# Define the rule to generate contigs from the hifiasm assembly
rule contigs:
    input:
        output_asm
    output:
        output_contigs
    shell:
        "awk '/^S/{{print \">{2}\\n\"{3}}}' {input}.p_ctg.gfa | fold > {output}"

# Define the rule to run BUSCO analysis on the contigs
rule busco:
    input:
        output_contigs
    output:
        directory(output_busco)
    shell:
        "busco -i {input} -o {output} -l {busco_lineage} -m genome"

# Define the rule to run BLASTn on the contigs
rule blastn:
    input:
        output_contigs
    output:
        output_blastn
    threads:
        180
    shell:
        "blastn -query {input} -db nt -outfmt '6 qseqid staxids bitscore std' "
        "-max_target_seqs 10 -max_hsps 1 -num_threads {threads} "
        "-evalue 1e-25 -out {output} &> log_blastn.txt"

# Define the rule to run BlobTools on the contigs
rule blobtools:
    input:
        output_contigs,
        "blastn.tsv"
    output:
        output_blobtools
    threads:
        192
    shell:
        "blobtools create --fasta {input[0]} --cov /home/biopatic/abbasi/P1--Cotesia-Pacbio/2-Assembly/2-Hafiasm/Hifiasm_on_CC/Inspector/Pacbio_CC_Hifiasm_out/read_to_contig.bam "
        "--hits {input[1]} --taxdump {taxdump} "
        "--busco ../../Hifiasm_on_CC/BUSCO/BUSCO_Hifiasm_Result/CC.hifiasm.p_ctgl0.busco.out/run_insecta_odb10/full_table.tsv "
        "--threads {threads} {output} &> log_blob"

