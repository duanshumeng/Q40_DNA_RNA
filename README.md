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
#fastp (v0.19.6),fastqc (v0.11.5),seqtk (v1.4-r122)
Analysis scripts/fastQC.smk
Analysis scripts/fastq-downsample.sh
Analysis scripts/bam-downsample.sh
```

**3.2. Sequence Alignment**
```
#BWA from sentieon (v202112.05),HISAT2 (v2.7.10),SAMtools (v1.18)
#DNA alignment
Analysis scripts/wes-fastq2bam.wdl

#RNA alignment
Analysis scripts/rnaseq-fastq2bam.sh
```

**3.3. Variant Calling (DNA Level)** 

```
#GATK (v4.4.0, including Mutect2, HaplotypeCaller), ASCAT (v3.1.2), hap.py (for germline variant evaluation), som.py (for somatic variant evaluation), VBT (for Mendelian consistency analysis), BCFtools(v1.17)
#Germline SNV/Indel
Analysis scripts/Haplotyper_calling-germline.sh

#Somatic SNV/Indel
Analysis scripts/Mutect2_calling-somatic.sh

#Somatic CNV
Analysis scripts/ascat-WES.r
Analysis scripts/CNV-calling.sh
```

**3.4. Transcriptome Analysis (RNA Level)**

```
#StringTie (v2.1.7)
Analysis scripts/gene-quantification.sh
```

**3.5. Quality Assessment & Visualization**
```
#Qualimap (v2.0.0),MultiQC (v1.18)
#Data quality Assessment
Visualization scripts/Fig2-data_quality.R

#Data performance Assessment
Visualization scripts/Fig3-4-WES_combine.R
Visualization scripts/Fig5-RNAseq_combine.R
Visualization scripts/Fig5.1.ERCC_RNAseq.V2.R
Visualization scripts/Fig5.2.RNA_Quartet.R
Visualization scripts/Fig5.3.MAQC2-RNAseq.R
```

**3.6. R Bioinformatics Packages**

R (v4.2+),ASCAT (v3.1.2),limma (v3.54.0), edgeR (v3.40.2), ggplot2, ggunchained, etc. (for plotting & visualization)
### 4. Dependent Key Databases/Reference Files

Reference Genome: GRCh38.d1.vd1.fa (hg38)

Gene Annotation File: Homo_sapiens.GRCh38.109.gtf

High-Confidence Variant Reference Sets: (e.g., standard variant sets provided by the NIST, Quartet, and SEQC2 projects)

