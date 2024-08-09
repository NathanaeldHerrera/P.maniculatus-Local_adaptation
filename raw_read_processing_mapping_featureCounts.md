## Processing raw RNAseq fastq read data
First, we will clean raw read data using [FastP](https://github.com/OpenGene/fastp)
```
for i in *R1_001.fastq.gz;
do
  name1=${i%-*};
  fastp -w 6 -i "$name1"-RV_R1_001.fastq.gz -I "$name1"-RV_R2_001.fastq.gz --out1 "$name1"_fastp_R1.fastq.gz --out2 "$name1"_fastp_R2.fastq.gz --unpaired1 "$name1"_fastp_1U.fastq.gz --unpaired2 "$name1"_fastp_2U.fastq.gz --detect_adapter_for_pe --cut_front --cut_front_window_size 5 --cut_front_mean_quality 20 -l 36 -j "$name1"_fastp.json -h "$name1"_fastp.html 2> "$name1".log
done
```
Collect FastP results for QC statistics using [Multiqc](https://multiqc.info/)
```
multiqc .
```
Now we can map and sort. We are using [Hisat2](http://daehwankimlab.github.io/hisat2/) for mapping.
We are using the P. maniculatus [assembly 2.1.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003704035.1/)
```
# Define reference:
PMAN_REF=/reference/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna

for i in *R1.fastq.gz;
do
        name1=$(echo $i | cut -d '_' -f 1);
if [ ! -f ./bams/"$name1"_Lung_Halign.bam ]
then
        hisat2 -p 14 --mp 2,0 -q -x $PMAN_REF -1 "$name1"_fastp_Lungs_R1.fastq.gz -2 "$name1"_fastp_Lungs_R2.fastq.gz -U "$name1"_fastp_Lungs_1U.fastq.gz --summary-file "$name1"_Lung_align.summ.txt --un "$name1"_Lung.unmapped.fq | samtools view -Sbo ./bams/"$name1"_Lung_Halign.bam
fi
# sort final bams
if [ ! -f ./bams/"$name1"_Lung_Halign_sort.bam ]
then
        samtools sort ./bams/"$name1"_Lung_Halign.bam -o ./bams/"$name1"_Lung_Halign_sort.bam
fi
done
```
Collect mapping QC metrics using [Qualimap](http://qualimap.conesalab.org/) bamQC.
```
for i in */*.bam;
do
        name1=${i%.*};
        qualimap bamqc -sd -bam "$name1".bam -nt 46
done
```
Finally, we generate count data using [featureCounts](https://subread.sourceforge.net/featureCounts.html). Each tissue (right ventricle and lung) is treated seperately.
```
# Set some environmental variables:
PMAN_GTF=/reference/GCF_003704035.1_HU_Pman_2.1.3_genomic.gtf
RV_OUT=/right_ventricle/counts/Pman_readcounts_RV_CountMM.txt
LUNG_OUT=/lung/counts/Pman_readcounts_LUNG_CountMM.txt

# Right Ventricle
featureCounts -p -O -F GTF -a $PMAN_GTF -o $RV_OUT *_Halign_sort.bam 2> featureCounts_RV_log.txt
# Lung
featureCounts -p -O -F GTF -a $PMAN_GTF -o $LUNG_OUT *_Halign_sort.bam 2> featureCounts_LUNG_log.txt
```
