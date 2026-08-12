# Snakemake PacBio Assembly

This repository contains two Snakemake workflows:

- `/home/runner/work/Snakemake-PacBio-Assembly/Snakemake-PacBio-Assembly/Pacbio-assembly-snakemake.smk`
- `/home/runner/work/Snakemake-PacBio-Assembly/Snakemake-PacBio-Assembly/reassembly_snakemake.smk`

## 1) `Pacbio-assembly-snakemake.smk`

Pipeline steps:
1. **FastQC** on `C.C.fastq.gz`
2. **hifiasm** assembly (`asm_CC.hifiasm.asm`)
3. Convert GFA to contigs fasta (`asm_CC.hifiasm.asm.p_ctgl2.fasta`)
4. **BUSCO** genome completeness analysis (`output/`)
5. **BLASTn** taxonomic hits (`blastn.tsv`)
6. **BlobTools** dataset generation (`tempo/dataset_0.blobDB.json`)

### Main hardcoded values used in the script
- Input reads: `C.C.fastq.gz`
- BUSCO lineage path: `/home/biopatic/abbasi/P1--Cotesia-Pacbio/2-Canu-Assembly/CC/CC_hifi/3-Busco_analysis/busco_downloads/lineages/insecta_odb10`
- Taxdump path: `../../../../3-DATABASES/taxdump`
- Coverage BAM used by BlobTools: `/home/biopatic/abbasi/P1--Cotesia-Pacbio/2-Assembly/2-Hafiasm/Hifiasm_on_CC/Inspector/Pacbio_CC_Hifiasm_out/read_to_contig.bam`

## 2) `reassembly_snakemake.smk`

Pipeline steps:
1. **bwa index** on `21.fasta`
2. **bwa aln** for paired inputs `input/P1.renamed.fasta` and `input/P2.renamed.fasta`
3. **bwa sampe** SAM generation
4. Read extraction with `extract.fasta.from.sam.using.list.pl`
5. Merge extracted fasta files into `all_reads.fasta`
6. **SPAdes** assembly into `assembly_try1/contigs.fasta`

### Main files expected by the script
- Reference: `21.fasta`
- Read lists: `input/{sample}_list.txt`
- Paired read files: `input/P1.renamed.fasta`, `input/P2.renamed.fasta`
- Trusted contigs for SPAdes: `2_freshwater.fasta`

## Run

Run either workflow with Snakemake:

```bash
snakemake -s Pacbio-assembly-snakemake.smk --cores <N>
```

or

```bash
snakemake -s reassembly_snakemake.smk --cores <N>
```

## Notes

- Both workflows currently include environment-specific absolute paths.
- Update tool paths and input/output file names before running in a new environment.

## Contact

- GitHub: [@nmabbasi](https://github.com/nmabbasi)
