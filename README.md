# RNA-Seq Protocol: HISAT, Stringtie, and Ballgown

### Following protocol by Pertea, M. *et al*:
[Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown](https://www.ncbi.nlm.nih.gov/pubmed/27560171)

---

### Tools Used (insatlled using *conda*):
- **HISAT2**: version 2.1.0
- **SAMtools**: version 1.9
- **StringTie**: version
- **IGV**: version 2.5.2
- **R**: version 3.4.1
- **Ballgown**: version

```bash
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda stringtie
conda install -c bioconda bioconductor-ballgown
conda install -c bioconda igv
```

### Data used:
Obtain and extract files from FTP server: <ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol>

```bash
wget ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
tar xvzf chrX_data.tar.gz
```

### Align reads using HISAT2:
Align all reads in the FASTQ format located in the *chrX_data/samples* directory to the indexed reference genome located in the *chrX_data/indexes* directory.  This will generate SAM files (one for each pair of FASTQ samples).  
- --dta option is used to report alignments tailored for transcript assemblers.
- -x is to indicate the reference genome used

```bash
hisat2 --dta -x chrX_data/indexes/chrX_tran \
-1 chrX_data/samples/ERR188044_chrX_1.fastq.gz \
-2 chrX_data/samples/ERR188044_chrX_2.fastq.gz \
-S analysis/hisat2-alignment/ERR188044_chrX.sam
```

To align all of the reads automatically, run the hisat2-alignment.sh Bash script
```bash
bash hisat2-alignment.sh
```

### Sort and convert generated SAM files to BAM using Samtools
Using Samtools sort command. Optionally, delete original SAM files
- -@ indicates number of threads to use
- -o indicates to write final output to FILE rather than standard output

```bash
samtools sort -@ 1 ./analysis/hisat2-alignment/ERR188044_chrX.sam -o ./analysis/hisat2-alignment/ERR188044_chrX.sorted.bam
rm ./analysis/hisat2-alignment/ERR188044_chrX.sam
```

To sort and convert all of the SAM files automatically, run the samtools-sort.sh Bash script
```bash
bash samtools-sort.sh
```

### Index sorted BAM files for downstream analysis and visualization with IGV
Using Samtools index command.  Will generate a .bai index file for each BAM file

```bash
samtools index ./analysis/hisat2-alignment/ERR188044_chrX.sorted.bam -o ./analysis/hisat2-alignment/ERR188044_chrX.sorted.bam.bai
```

To index all of the BAM files automatically, run the samtools-index.sh Bash script
```bash
bash samtools-index.sh
```
