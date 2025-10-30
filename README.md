# Q40 DNA and RNA
The scripts related to comparing "Higher Q score (Q40) Sequencing Technology Enhances Accuracy and Cost-Efficiency in Genomic Events Detection"

The scripts used to analysis datasets are saved in branch "Analysis-scripts"
The scripts used to visualization are saved in branch "Visualization scripts"

## Software Runtime Environment/Software

### 1. Operating System

Linux Operating System
### 2. Programming Languages & Core Runtimes

Python (v3.8+)
R (v4.2+)
Bash/Shell
### 3. Bioinformatics Core Tools & Software Packages

**3.1. Data Quality Control & Preprocessing**
```
#script: Analysis scripts/fastQC.smk, Analysis scripts/fastq-downsample.sh
fastp (v0.19.6)
fastqc (v0.11.5)
seqtk (v1.4-r122)
```


**3.2. Sequence Alignment**

BWA (v0.7.17-r1188)
HISAT2 (v2.7.10)
SAMtools

**3.3. Variant Calling (DNA Level)**

GATK (v4.4.0, including Mutect2, HaplotypeCaller, FilterMutectCalls)
ASCAT (v3.1.2)
hap.py
som.py (for somatic variant evaluation)
VBT (for Mendelian consistency analysis)
BCFtools

**3.4. Transcriptome Analysis (RNA Level)**

StringTie (v2.1.7, for transcript assembly & quantification)
limma (v3.54.0)
edgeR (v3.40.2, for differential expression analysis)

**3.5. Quality Assessment & Visualization**

Qualimap (v2.0.0)
MultiQC (v1.18)

**3.6. R Bioinformatics Packages**

ASCAT (for copy number analysis)
limma, edgeR (listed above in 3.4)
ggplot2, ggunchained, etc. (for plotting & visualization)
### 4. Dependent Key Databases/Reference Files

Reference Genome: GRCh38.d1.vd1.fa (hg38)
Gene Annotation File: Homo_sapiens.GRCh38.109.gtf
High-Confidence Variant Reference Sets: (e.g., standard variant sets provided by the NIST, Quartet, and SEQC2 projects)

