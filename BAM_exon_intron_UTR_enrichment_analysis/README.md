# Splice Event and Exon Enrichment Analysis from Bulk RNA-seq
## Overview

This repository contains a reproducible workflow to detect splice junctions, infer unannotated (cryptic) exon usage, and quantify exon-level expression from bulk RNA-seq data using aligned BAM files and splice junction output. The pipeline is designed for high-resolution analysis of alternative splicing and exon inclusion within a defined genomic locus and is particularly suited for disease-associated splicing abnormalities.

The current implementation focuses on UBA1 splicing, motivated by studies of VEXAS syndrome and myelodysplastic syndromes (MDS), where aberrant splicing and exon usage play a critical pathogenic role.

The workflow integrates:

- STAR splice junction outputs (`SJ.out.tab`)

- BAM-level read counting

- Gene model annotation (GTF)

- Exon-level FPKM normalization

- Prediction of putative exons based on splice junction gaps

## Motivation

Alternative splicing is a major contributor to transcriptomic diversity and disease pathogenesis. While many tools detect differential splice junction usage, cryptic exon inclusion and exon-level expression changes are often under-characterized - particularly when they do not fully overlap annotated gene models.

This pipeline addresses that gap by:

1. Mapping observed splice junctions to annotated exons

2. Identifying gaps between junctions consistent with exon inclusion

3. Quantifying read support for these inferred exon regions

4. Normalizing exon-level signal relative to gene-level expression

5. Linking predicted exons to known isoforms when possible

This approach enables mechanistic interpretation of splicing abnormalities beyond junction-level statistics alone.

## Input Files and Formats
### 1. Aligned BAM files
### Example of STAR aligner required output:

`VEXAS_CD14Aligned.sortedByCoord.out.bam`

`VEXAS_CD14Aligned.sortedByCoord.out.bam.bai`

### 2. STAR splice junction files
`SJ.out.tab`

### 3. Gene annotation file (.gtf)

## Methods
### 1. Exon Coordinate Extraction

Annotated exon coordinates are extracted directly from the GTF file using `awk`. Exons are indexed by genomic coordinates and associated transcript (isoform) identifiers.

### 2. Splice Junction Filtering

Splice junctions are filtered to a user-defined genomic window (currently UBA1 on chrX:47190861–47215128). Junctions are processed sequentially to preserve positional order.

### 3. Predicted Exon Identification

A putative exon is inferred when:

- There is a gap between consecutive splice junctions

- The gap exceeds a 3 bp tolerance

- The region lies within the gene locus

This captures exon-like regions consistent with cryptic exon inclusion or altered splice site usage.

### 4. Read Quantification

For each predicted exon:

- Reads overlapping the exon region are counted using `samtools view -c`

- Total gene-level reads are calculated for normalization

### 5. Expression Normalization

Exon expression is reported as FPKM, normalized to:

- Exon length

- Total reads mapped to the gene locus

<img width="384" height="64" alt="Image" src="https://github.com/user-attachments/assets/3e3c2839-ca76-4950-80c4-1cc569bd9923" />

### 6. Isoform Mapping

Predicted exon coordinates are compared against annotated exon boundaries with a ±3 bp tolerance to determine isoform overlap, allowing assignment of exon usage to known transcripts when possible.

<img width="843" height="294" alt="Image" src="https://github.com/user-attachments/assets/ccc1224b-c132-4fd8-9169-19dc7fa5ca6b" />
