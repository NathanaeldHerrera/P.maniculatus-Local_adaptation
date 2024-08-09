### Load Packages ###########################
message = FALSE # hide package loading messages
rm(list=ls())
library(stats)
library(ggbiplot)
library(qwraps2)
library(report)
#==============================================================================
# We will then run principal components analysis and generate a 95% data ellipse to identify 
# technical outliers. 
# read_in_metrics ########################################################
run_date <- format(Sys.Date(), "%d%b%Y")

#metrics <- read.csv("./right_ventricle/rnaseq_Pman_RV_QCmetrics.csv",header=TRUE,stringsAsFactors=FALSE)
metrics <- read.csv("./lung/rnaseq_Pman_LUNG_QCmetrics.csv",header=TRUE,stringsAsFactors=FALSE)
head(metrics)
metrics$group_for_pca <- "samples"
metrics.pca <- prcomp(metrics[,c(4:6,8:9)], center = TRUE,scale. = TRUE) #don't include perc dups or # reads 
summary(metrics.pca)
# output to Console by sample ID
ggbiplot(metrics.pca, ellipse=TRUE, labels=metrics$sample_name,groups=metrics$group_for_pca,circle=FALSE,colour="black",ellipse.prob = 0.95) + scale_colour_manual(values=c("black")) + ggtitle("PCA of P. maniculatus lung RNASeq samples with 95% ellipse")
ggbiplot(metrics.pca, ellipse=TRUE, groups=metrics$group_for_pca,circle=FALSE,colour="black",ellipse.prob = 0.95) + scale_colour_manual(values=c("black")) + ggtitle("PCA of P. maniculatus lung RNASeq samples with 95% ellipse")

# remove 2 samples for RV and calculate averages for metrics
metrics_mod <- metrics[metrics$sample_name!=c("CC3MEBF"),]
metrics_mod <- metrics_mod[metrics_mod$sample_name!=c("KNG2304G"),]

# remove 5 samples for LUNG and calculate averages for metrics
metrics_mod <- metrics[metrics$sample_name!=c("AC1KNO"),]
metrics_mod <- metrics_mod[metrics_mod$sample_name!=c("CC2KNF"),]
metrics_mod <- metrics_mod[metrics_mod$sample_name!=c("AA4KNV"),]
metrics_mod <- metrics_mod[metrics_mod$sample_name!=c("CC1KNN"),]
metrics_mod <- metrics_mod[metrics_mod$sample_name!=c("CC4MEX"),]

set.seed(42)
options(qwraps2_markup = "markdown")

#LUNG_seq_summary <-
RV_seq_summary <-
  list("Mean Quality Score" =
         list("min"       = ~ min(mean_quality_score),
              "max"       = ~ max(mean_quality_score),
              "mean (sd)" = ~ qwraps2::mean_sd(mean_quality_score)),
       "Percent Bases Above Q30" =
         list("min"       = ~ min(perc_bases_above_q30),
              "median"    = ~ median(perc_bases_above_q30),
              "max"       = ~ max(perc_bases_above_q30),
              "mean (sd)" = ~ qwraps2::mean_sd(perc_bases_above_q30)),
       "GC (%)" =
         list("min"       = ~ min(perc_gc),
              "median"    = ~ median(perc_gc),
              "max"       = ~ max(perc_gc),
              "mean (sd)" = ~ qwraps2::mean_sd(perc_gc)),
       "Reads" =
         list("min"       = ~ min(reads),
              "median"    = ~ median(reads),
              "max"       = ~ max(reads),
              "mean (sd)" = ~ qwraps2::mean_sd(reads)),
       "Genome Alignment Rate (%)" =
         list("min"       = ~ min(align_rate),
              "max"       = ~ max(align_rate),
              "mean (sd)" = ~ qwraps2::mean_sd(align_rate)),
       "Assignment Rate (%)" =
         list("min"       = ~ min(assign_rate),
              "max"       = ~ max(assign_rate),
              "mean (sd)" = ~ qwraps2::mean_sd(assign_rate))
  )

### Overall
#table <- summary_table(metrics_mod, LUNG_seq_summary)
table <- summary_table(metrics_mod, RV_seq_summary)
table

# Run anova to check for batch effects
lm_qs <- lm(mean_quality_score ~ group, data=metrics)
lm_q30 <- lm(perc_bases_above_q30 ~ group, data=metrics)
lm_gc <- lm(perc_gc ~ group, data=metrics)
lm_align <- lm(align_rate ~ group, data=metrics)
lm_assign <- lm(assign_rate ~ group, data=metrics)

report(anova(lm_qs))
report(anova(lm_q30))
report(anova(lm_gc))
report(anova(lm_align))
report(anova(lm_assign))
