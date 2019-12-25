# RNA-Seq Protocol: HISAT, Stringtie, and Ballgown

## Following protocol by Pertea, M. *et al*:
[Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown](https://www.ncbi.nlm.nih.gov/pubmed/27560171)

---

### Tools Used:
- **HISAT2**: version 2.1.0, installed using *conda*
```bash
conda install -c bioconda hisat2
```

- **SAMtools**: version 1.9, installed using *conda*
```bash
conda install -c bioconda samtools
```

- **StringTie**: version, installed using *conda*
```bash
conda install -c bioconda stringtie
```

- **R**: version 3.4.1

- **Ballgown**: version, installed using 

### Data used:
Obtained from ftp site using **wget**: <ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz>
```bash
wget ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
```

Extract files:
```bash
tar xvzf chrX_data.tar.gz
```



