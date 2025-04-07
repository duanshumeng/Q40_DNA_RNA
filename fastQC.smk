import os
import pandas as pd
# Define the paths
# Read sample names and VCF file paths from input.csv
samples = pd.read_csv("/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/LC_RM/RNA_results/NQI_RNAseq_info_20250226.txt",sep='\t')
print(samples)

sample_names = samples['Sample_name'].tolist()
print(sample_names)

def get_fastq(sample,read):
    print('============',sample,'============')
    fq = samples.loc[samples['Sample_name'] == sample, read].values[0]
    return fq

output_path = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/LC_RM/RNA_results/NQI_RNAseq" # You will pass this as a parameter or specify it here

TRIMMED_OUTPUT_PATH=f"{output_path}.clean_data"
print(TRIMMED_OUTPUT_PATH)
FASTQC_OUTPUT_PATH=f"{output_path}.fastqc"
print(FASTQC_OUTPUT_PATH)

 
if not os.path.exists(TRIMMED_OUTPUT_PATH):
    print('====mkdir file====',TRIMMED_OUTPUT_PATH)
    os.makedirs(TRIMMED_OUTPUT_PATH)


if not os.path.exists(FASTQC_OUTPUT_PATH):
    print('====mkdir file====',FASTQC_OUTPUT_PATH)
    os.makedirs(FASTQC_OUTPUT_PATH)

# Rule to read the input file and preprocess
rule all:
    input:
        expand(TRIMMED_OUTPUT_PATH+"/{Sample_name}_1.trimmed.fq.gz", Sample_name=sample_names), # Add more samples here
        expand(TRIMMED_OUTPUT_PATH+"/{Sample_name}_2.trimmed.fq.gz", Sample_name=sample_names),
        expand(FASTQC_OUTPUT_PATH+"/{Sample_name}.tar.gz", Sample_name=sample_names)


# Rule for trimming fastq files using fastp
rule trim_fastq:
    input:
        fastq1 = lambda wc: get_fastq(wc.Sample_name,'fq1'),
        fastq2 = lambda wc: get_fastq(wc.Sample_name,'fq2')
    output:
        trimmed_fastq1 = TRIMMED_OUTPUT_PATH+"/{Sample_name}_1.trimmed.fq.gz",
        trimmed_fastq2 = TRIMMED_OUTPUT_PATH+"/{Sample_name}_2.trimmed.fq.gz",
        html = TRIMMED_OUTPUT_PATH+"/{Sample_name}.html",
        json = TRIMMED_OUTPUT_PATH+"/{Sample_name}.json"
    conda:
        "mapping_V2"
    shell:
        """
        fastp --thread 16 -i {input.fastq1} -I {input.fastq2} \
        -o {output.trimmed_fastq1} -O {output.trimmed_fastq2} \
        -h {output.html} \
        -j {output.json}
        """

rule fastqc:
    input:
        fastq1 = TRIMMED_OUTPUT_PATH+"/{Sample_name}_1.trimmed.fq.gz",
        fastq2 = TRIMMED_OUTPUT_PATH+"/{Sample_name}_2.trimmed.fq.gz"
    output:
        out_dir = FASTQC_OUTPUT_PATH+"/{Sample_name}",
        out_file = FASTQC_OUTPUT_PATH+"/{Sample_name}.tar.gz"
    conda:
        "mapping_V2"
    shell:
        "mkdir -p {output.out_dir};fastqc {input.fastq1} {input.fastq2} -o {output.out_dir};tar -czvf {output.out_file} {output.out_dir}"