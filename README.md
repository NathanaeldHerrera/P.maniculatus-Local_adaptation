# P. maniculatus RNA-Seq Pipeline. 

This repository is for data processing and analyses related to Bautista et al. (in review), "Contributions of genotypic specialization and adaptive plasticity in the evolved resistance to hypoxic cold stress in high-altitude deer mice".

Our analysis pipeline is similar to Schweizer et al., 2023. Gene regulatory changes underlie developmental plasticity in respiration and aerobic performance in highland deer mice. Molecular Ecology. 32:3483–3496 and have used some scripts with modification here (see: https://github.com/renaschweizer/pman_diaphragm_rnaseq). 

Here, our workflow for processing and analyzing raw RNA-Seq fastq data for the right ventricle of the heart and lung tissue to generate count data. R scripts are provided for each tissue where users will find step-by-step workflows to perform gene expression analyses in R, including quality control of sequence data, identifying regulatory modules in WGCNA, testing for correlations of WGCNA modules with phenotypes, and testing for effects of treatment/population on module expression.

In each section below, I describe the pipeline used (follow links below). Scripts referenced in each section can be found inside the respective folders.

Contact:ndh04c (at) gmail.com

Sections

1. [Processing raw RNAseq fastq read data](https://github.com/NathanaeldHerrera/Pman_rnaseq/blob/main/raw_read_processing_mapping_featureCounts)


