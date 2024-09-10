# ChIP-Seq_pipeline
A ChIP-Seq pipeline

### Software requirement
- cutadapt
- bowtie
- samtools
- bedtools
- bedToBigWig
- macs
- removeDupBed(from Weilab)
- peakAnnotate(from Weilab)

### Usage
sh ChIPseq_pipeline.sh input_fq1,input_fq2 out_dir out_label bowtie_index_dir n_cpu n_trim5 n_trim3 n_mismatch genome_name gtf_path

### Note
removeDupBed and peakAnnotate used in this pipeline are customized scripts by Weilab members, which are stored at [1.perl_scripts_from_Weilab](./1.perl_scripts_from_Weilab).
