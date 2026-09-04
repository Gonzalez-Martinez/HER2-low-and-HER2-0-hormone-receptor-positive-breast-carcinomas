# HER2-low and HER2-0 hormone receptor-positive breast carcinomas

Single-nucleus RNA sequencing analysis of HER2-low and HER2-0 HR+ breast carcinomas.

## Project Overview
This repository contains the code and analysis workflows used to:
- Perform cell-type compositional analysis
- Compare HER2 groups (0, Low)
- Calculate inter-tumoral variability (CV)
- Generate publication-ready figures

## Analysis
All analyses were performed in R using:
- Seurat, SCTransform, DoubletFinder
- inferCNV
- speckle (propeller)
- DESeq2, edgeR, fgsea, genefu
- dplyr, tidyr, ggplot2, ggpubr

## Data availability

Raw sequencing data (FASTQ files) for the 20 tumours are deposited in the Gene Expression Omnibus under accession GSE318213.

### Obtaining the input data

1. Download the FASTQ files for all 20 samples from GEO (GSE318213).
2. Generate count matrices using Cell Ranger v9.0.0 with the 10x Genomics GRCh38 human reference, following the 10x Genomics protocols CG000632 (FFPE sample preparation) and CG000477 (single-nucleus library preparation).
3. Place the resulting count matrices in a `data/` directory at the repository root.
4. The Seurat objects used in the analyses are generated from these matrices by the scripts in this repository; they are not deposited separately.

## Author
Silvia González-Martínez
