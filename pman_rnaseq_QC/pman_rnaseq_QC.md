## Processed read data QC
Here, we run PCA and generate a 95% data ellipse to identify 
technical outliers for each tissue and then summarize sequencing metrics.

For each tissue we have a table containing:

|          |           |
|----------|-----------|
|sample_ID | batch| perc_dups| mean_quality_score|	perc_bases_above_q30|	perc_gc| num_reads|	alignment_rate|	assignment_rate|

alignment rate is generated during mapping by HiSat2 and assignment rate comes from featureCounts


See R script:[pman_rnaseq_QC.R]([http://daehwankimlab.github.io/hisat2/](https://github.com/NathanaeldHerrera/Pman_rnaseq/blob/main/pman_rnaseq_QC.R))
