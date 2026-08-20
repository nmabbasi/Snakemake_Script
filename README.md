# Snakemake PacBio Assembly & Reassembly Pipelines 🐍

[![Snakemake](https://img.shields.io/badge/Snakemake-≥6.0.0-brightgreen.svg?style=flat)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)]()

This repository contains two robust Snakemake pipelines tailored for the genome assembly, quality assessment, and targeted reassembly of PacBio sequencing data (specifically HiFi reads).

## 🚀 Pipelines Overview

### 1. Primary Assembly Pipeline (`Pacbio-assembly-snakemake.smk`)
An automated end-to-end workflow for *de novo* genome assembly of PacBio HiFi reads and subsequent taxonomic/quality evaluation.

**Key Steps:**
1. **Quality Control:** Validates input FastQ reads using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. **De Novo Assembly:** Assembles long reads using [Hifiasm](https://github.com/chhylp123/hifiasm), optimized for diploid/haploid phasing.
3. **Contig Extraction:** Parses `.gfa` graphs to extract primary contigs (`.fasta`).
4. **Completeness Assessment:** Evaluates assembly completeness using [BUSCO](https://busco.ezlab.org/) against the `insecta_odb10` lineage (customizable).
5. **Taxonomic Profiling:** Queries contigs against the NCBI `nt` database using `blastn`.
6. **Contamination Screening:** Generates Taxon-Annotated GC-Coverage (TAGC) plots using [BlobTools](https://github.com/DRL/blobtools) by integrating the assembly, BAM coverage, BLAST hits, and BUSCO results.

### 2. Targeted Reassembly Pipeline (`reassembly_snakemake.smk`)
A highly specialized pipeline designed for extracting subset reads aligned to a specific region/reference and reassembling them using SPAdes.

**Key Steps:**
1. **Reference Indexing:** Indexes the target/reference genome using [BWA](https://github.com/lh3/bwa).
2. **Alignment:** Maps paired-end reads (`P1`, `P2`) to the reference genome via BWA-ALN / BWA-SAMPE.
3. **Read Extraction:** Extracts reads of interest from the resulting SAM files utilizing a custom Perl script (`extract.fasta.from.sam.using.list.pl`).
4. **Targeted Assembly:** Concatenates extracted reads and reassembles them *ab initio* utilizing [SPAdes](https://github.com/ablab/spades) (`--careful`, `--only-assembler`) incorporating trusted contigs.

---

## 💻 Installation & Requirements

Ensure you have [Conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your HPC or local machine.

### Dependencies
The pipelines rely on the following bioinformatics tools being available in your environment or `$PATH`:
- `fastqc`
- `hifiasm`
- `busco`
- `blastn`
- `blobtools`
- `bwa` (v0.6.2+)
- `spades.py` (v3.6.1+)
- `perl`

*Tip: It is highly recommended to manage these dependencies using a Conda environment.*

---

## ⚙️ Usage

### Running the Primary Assembly Pipeline

1. Edit the input variables at the top of `Pacbio-assembly-snakemake.smk`:
   ```python
   input_fastq = "your_reads.fastq.gz"
   threads = 92
   taxdump = "/path/to/taxdump"
   busco_lineage = "/path/to/lineage_odb10"
   ```
2. Execute the workflow:
   ```bash
   snakemake -s Pacbio-assembly-snakemake.smk --cores 92
   ```

### Running the Targeted Reassembly Pipeline

1. Ensure your inputs (FASTA references and subset lists) are accurately defined within the rule parameters of `reassembly_snakemake.smk`.
2. Execute the workflow:
   ```bash
   snakemake -s reassembly_snakemake.smk --cores 30
   ```

---

## 📬 Contact & Support

**Nasir Mahmood Abbasi, PhD**  
*Computational Biologist & Single-Cell Transcriptomics Expert*  
🌐 Portfolio: [nmabbasi.github.io](https://nmabbasi.github.io)  
📚 Educational Portal: [The Omics Hub](https://theomicshub.com)  

If you use these pipelines in your research, please consider citing or acknowledging this repository. For bug reports or feature requests, feel free to open an issue!
