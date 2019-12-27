# RNA-Seq Protocol: HISAT, Stringtie, and Ballgown

### Following protocol by Pertea, M. *et al*:
[Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown](https://www.ncbi.nlm.nih.gov/pubmed/27560171)

---

### Tools Used (insatlled using *conda*):
- **HISAT2**: version 2.1.0
- **SAMtools**: version 1.9
- **StringTie**: version
- **IGV**: version 
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
- Obtain and extract files from FTP server: <ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol>

```bash
wget ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
tar xvzf chrX_data.tar.gz
```
