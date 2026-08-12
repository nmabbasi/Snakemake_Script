# Snakemake PacBio Assembly

This repository contains two Snakemake workflows for PacBio-based assembly tasks.

## Workflows

- `/home/runner/work/Snakemake-PacBio-Assembly/Snakemake-PacBio-Assembly/Pacbio-assembly-snakemake.smk`  
  HiFi assembly and QC workflow using FastQC, Hifiasm, BUSCO, BLAST, and BlobTools.

- `/home/runner/work/Snakemake-PacBio-Assembly/Snakemake-PacBio-Assembly/reassembly_snakemake.smk`  
  Read extraction and re-assembly workflow using BWA, custom extraction script, and SPAdes.

## Requirements

Install and make available on `PATH`:

- `snakemake` (7+)
- `fastqc`
- `hifiasm`
- `busco`
- `blastn`
- `blobtools`
- `bwa`
- `spades.py`
- `perl` (for `extract.fasta.from.sam.using.list.pl`)

## Expected Inputs

### PacBio assembly workflow

- `C.C.fastq.gz`
- `read_to_contig.bam`
- `busco_downloads/lineages/insecta_odb10`
- `../../../../3-DATABASES/taxdump`

### Reassembly workflow

- `21.fasta`
- `2_freshwater.fasta`
- `input/P1.renamed.fasta`
- `input/P2.renamed.fasta`
- `input/2_fw_Tous12m_list.txt`
- `input/2_fw_Tous25m_list.txt`
- `extract.fasta.from.sam.using.list.pl`

## Usage

Run the PacBio assembly workflow:

```bash
snakemake -s Pacbio-assembly-snakemake.smk -j 16
```

Run the reassembly workflow:

```bash
snakemake -s reassembly_snakemake.smk -j 16
```

Preview planned jobs without running:

```bash
snakemake -s Pacbio-assembly-snakemake.smk -n
snakemake -s reassembly_snakemake.smk -n
```

## Outputs

### PacBio assembly workflow

- `tempo/dataset_0.blobDB.json` (final target)
- `asm_CC.hifiasm.asm.p_ctgl2.fasta`
- `blastn.tsv`
- `output/` (BUSCO output directory)

### Reassembly workflow

- `assembly_try1/contigs.fasta` (final target)
- `all_reads.fasta`
- `21_m_<sample>.sam`
- `<sample>_final.fasta`

## Notes

- Tool and database paths are intentionally defined at the top of each Snakefile for easier editing.
- If your CPU/memory differ, adjust per-rule thread values in the Snakefiles.
