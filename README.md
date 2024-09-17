# ChIP-Seq_pipeline
A ChIP-Seq pipeline

### Software requirement
- cutadapt
- bowtie
- samtools
- bedtools
- bedToBigWig
- macs
- picard
- HOMER

### Usage
```
sh ChIPseq_pipeline_cutadapt.sh input_fq1,input_fq2 out_dir out_label bowtie_index_dir n_cpu n_trim5 n_trim3 n_mismatch genome_name gtf_path
```
Here is an example of a figure using the results from this pipeline:  
![image](https://github.com/maxuying1218/ChIP-Seq_pipeline/blob/main/figures/Peak_Annotation.jpg)  
PCC and PCA plot can be drawn using peaks like this:  
![image](https://github.com/maxuying1218/ChIP-Seq_pipeline/blob/main/figures/PCC_PCA_example.jpg)  
PCA plot can be drawn using [2.PCA_cal_plot.R](https://github.com/maxuying1218/RNA-Seq_pipeline/blob/main/1.mapping_counting/2.PCA_cal_plot.R).


