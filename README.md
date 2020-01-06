# RNA-Seq Protocol: HISAT, Stringtie, and Ballgown

### Following protocol by Pertea, M. *et al*:
[Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown](https://www.ncbi.nlm.nih.gov/pubmed/27560171)

---

### Tools Used (insatlled using *conda*):
- **HISAT2**: version 2.1.0
- **SAMtools**: version 1.9
- **StringTie**: version 2.0
- **Gffcompare**: version 0.11.2
- **IGV**: version 2.5.2
- **R**: version 3.6.1
- **Ballgown**: version 2.18.0

```bash
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda stringtie
conda install -c bioconda gffcompare
conda install -c bioconda bioconductor-ballgown
conda install -c bioconda igv
```

### Data used:
Obtain and extract files from FTP server: <ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol>

```bash
wget ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
tar xvzf chrX_data.tar.gz
```

---

### Align reads using **HISAT2**:
Align all reads in the FASTQ format located in the *chrX_data/samples* directory to the indexed reference genome located in the *chrX_data/indexes* directory.  This will generate SAM files (one for each pair of FASTQ samples).  
- -p indicates number of threads to use
- --dta option is used to report alignments tailored for transcript assemblers.
- -x is to indicate the reference genome used

```bash
hisat2 -p 1 --dta -x chrX_data/indexes/chrX_tran \
-1 chrX_data/samples/ERR188044_chrX_1.fastq.gz \
-2 chrX_data/samples/ERR188044_chrX_2.fastq.gz \
-S analysis/hisat2-alignment/ERR188044_chrX.sam
```

To align all of the reads automatically, run the *hisat2-alignment.sh* Bash script
```bash
bash hisat2-alignment.sh
```

---

### Sort and convert generated SAM files to BAM using **Samtools**:
Using *samtools sort* command. Optionally, delete original SAM files
- -@ indicates number of threads to use
- -o indicates to write final output to FILE rather than standard output

```bash
samtools sort -@ 1 ./analysis/hisat2-alignment/ERR188044_chrX.sam -o ./analysis/hisat2-alignment/ERR188044_chrX.sorted.bam
rm ./analysis/hisat2-alignment/ERR188044_chrX.sam
```

To sort and convert all of the SAM files automatically, run the *samtools-sort.sh* Bash script
```bash
bash samtools-sort.sh
```

To view a BAM file including the header, use *samtools view* command and pipe the output to *less*
```bash
samtools view -h analysis/hisat2-alignment/ERR188044_chrX.sorted.bam | less -S
```

---

### Index sorted BAM files for downstream analysis and visualization with **IGV**:
Using *samtools index* command.  Will generate a .bai index file for each BAM file

```bash
samtools index ./analysis/hisat2-alignment/ERR188044_chrX.sorted.bam -o ./analysis/hisat2-alignment/ERR188044_chrX.sorted.bam.bai
```

To index all of the BAM files automatically, run the *samtools-index.sh* Bash script
```bash
bash samtools-index.sh
```

---

### Assemble transcripts for each sample using **stringtie**:
- -p indicates number of threads to use (default: 1)
- -G indicates reference annotation to use (GTF/GFF)
- -o indicates output file name for assembled transcripts GTF (default: stdout)
- -l indicates name prefix for output transcripts (default: STRG)
- last, provide input BAM file

```bash
stringtie -p 1 -G chrX_data/genes/chrX.gtf -o ERR188044_chrX.gtf ./analysis/hisat2-alignment/ERR188044_chrX.bam
```

To assemble all transcripts, run *stringtie-assemble.sh* Bash script
```bash
bash stringtie-assemble.sh
```

To see number of records in each file (minus the header), run:
```bash
for f in analysis/stringtie-assemble/*.gtf; do grep -v '^#' $f | wc -l; done
```

---

### Merge transcripts from all samples
Generate a *mergelist.txt* file containing the names of the generated GTF files
```bash
for f in analysis/stringtie-assemble/*.gtf; do echo $f >> "analysis/stringtie-assemble/mergelist.txt"; done
```

Use *stringtie --merge* command to merge all of the GTF files
- -p indicates number of threads to use (default: 1)
- -G indicates reference annotation to include in the merging (GTF/GFF)
- -o indicates output file name for the merged transcripts GTF (default: stdout)
- -l indicates name prefix for output transcripts (default: MSTRG)
- last, provide generated mergelist text filw with names of all GTF files to merge

```bash
stringtie --merge -p 1 -G ./chrX_data/genes/chrX.gtf \
-o ./analysis/stringtie-assemble/stringtie_merged.gtf \
./analysis/stringtie-assemble/mergelist.txt
```

To generate the *mergelist.txt* and run the *stringtie --merge* command automatically, run the *stringtie-merge-assemble.sh* Bash script:
```bash
bash stringtie-merge-assemble.sh
```

---

### Compare transcripts with reference annotation using **gffcompare**:
- -r indicates reference annotation file (GTF/GFF)
- -G option tells gffcompare to compare all transcripts in the input transcripts.gtf file, even those that might be redundant
- -o indicates output prefix
- last, provide the merged GTF file
- redirect stderr to a log file

```bash
gffcompare -r chrX_data/genes/chrX.gtf -G \
-o analysis/stringtie-assemble/merged analysis/stringtie-assemble/stringtie-merged.gtf \
2> analysis/log/gffcompare/gffcompare.log
```

---

### Estimate transcript abundances using **stringtie**:
Use *stringtie* to estimate transcript abundances and create table counts for *Ballgown*.  It will generate several table files for each sample in a separate sample directory.
- -p number of threads (CPUs) to use (default: 1)
- -e only estimate the abundance of given reference transcripts (requires -G)
- -G reference annotation to use for guiding the assembly process (GTF/GFF3)
- -B enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)

```bash
stringtie -p 1 -e -B -G ./analysis/stringtie-assemble/stringtie-merged.gtf \
-o ./analysis/ballgown/ERR188044/ERR188044_chrX.gtf \
./analysis/hisat2-alignment/ERR188044_chrX.sorted.bam
```
