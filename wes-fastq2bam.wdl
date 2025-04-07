version 1.0

# 在此定义您的分析流程
workflow WES_Sentieon {
    
    # 流程输入文件和参数
    input {
        String sample_id
        File fastq1
        File fastq2
        String Seq_platform

        String ref_fasta = "GRCh38.d1.vd1.fa"
        File ref_dir = "oss://hcc1395/reference_genome/GRCh38.d1.vd1/"
        File regions = "oss://hcc1395/reference_genome/Bed_files/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_noy_cut.bed"

        File dbsnp_dir = "oss://hcc1395/reference_genome/GRCh38.d1.vd1/"
        File dbmills_dir = "oss://hcc1395/reference_genome/GRCh38.d1.vd1/"
        String dbsnp = "dbsnp_146.hg38.vcf"
        String db_mills = "Mills_and_1000G_gold_standard.indels.hg38.vcf"
    }
    
    # 调用子任务
    call SentieonFastqToBam {
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            sample_id=sample_id,
            ref_fasta=ref_fasta,
            ref_dir=ref_dir,
            Seq_platform=Seq_platform
    }


    call Sentieon_BQSR{
        input:
            ref_dir = ref_dir,
            dbsnp_dir = dbsnp_dir,
            dbmills_dir = dbmills_dir,
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            regions = regions,
            dbsnp = dbsnp,
            db_mills = db_mills,
            deduped_bam = SentieonFastqToBam.deduped_bam,
            deduped_bam_index = SentieonFastqToBam.deduped_bam_bai
    }

    output {
            File deduped_bam = SentieonFastqToBam.deduped_bam
            File deduped_bam_bai = SentieonFastqToBam.deduped_bam_bai
            Array[File] qc_metrics = SentieonFastqToBam.qc_metrics
            File recal_table = Sentieon_BQSR.recal_table
            File recal_post = Sentieon_BQSR.recal_post
            File recaled_bam = Sentieon_BQSR.recaled_bam
            File recaled_bam_index = Sentieon_BQSR.recaled_bam_index
            File recal_csv = Sentieon_BQSR.recal_csv
            File bqsrreport_pdf = Sentieon_BQSR.bqsrreport_pdf
            File vcf = Sentieon_BQSR.vcf
            File vcf_idx = Sentieon_BQSR.vcf_idx
    }

}

# 在此定义流程中使用的工具
task SentieonFastqToBam {
    # 工具输入文件和参数
    input {
        File fastq1
        File fastq2
        String sample_id
        String ref_fasta
        File ref_dir
        String Seq_platform
        

        ## Extra driver parameters
        #String qc_driver_args = ""
        String lc_driver_args = "--traverse_param=200000/10000"
        String dedup_driver_args = "--traverse_param=200000/10000"
        ## Extra algo parameters
        String bwa_args = "-Y -M"
        String bwa_chunk_size = "100000000"
        # String lc_args = ""
        String bam_option = "--bam_compression 1"

        Int cpu = 32
        String memory = "64G"
        String disks = "local-disk 250 cloud_ssd"
    }
    
    String out_bam = sample_id + ".dedup.bam"
    String out_bai = sample_id + ".dedup.bam.bai"

    # 工具运行命令
    command <<<
        set -exo pipefail
        
        mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
        export bwt_max_mem="$((mem_kb / 1024 / 1024 - 2))g"
        export MALLOC_CONF=lg_dirty_mult:-1
        export LD_PRELOAD=/usr/local/sentieon-genomics/lib/libjemalloc.so
        
        sentieon bwa mem -R "@RG\tID:~{sample_id}\tSM:~{sample_id}\tPL:~{Seq_platform}" ~{bwa_args} -K ~{bwa_chunk_size} -t ~{cpu} ~{ref_dir}/~{ref_fasta} ~{fastq1} ~{fastq2} \
        | sentieon util sort ~{bam_option} -i - -r ~{ref_dir}/~{ref_fasta} -t ~{cpu} -o ~{sample_id}.sorted.bam --sam2bam

        sentieon driver -r ~{ref_dir}/~{ref_fasta} -t ~{cpu} -i ~{sample_id}.sorted.bam ~{lc_driver_args} \
         --algo LocusCollector \
         ~{sample_id}.score.txt.gz

        sentieon driver -r ~{ref_dir}/~{ref_fasta} -t ~{cpu} -i ~{sample_id}.sorted.bam ~{dedup_driver_args} \
         --algo Dedup \
         --score_info ~{sample_id}.score.txt.gz \
         --metrics ~{sample_id}.dedup_metrics.txt \
         ~{bam_option} ~{out_bam}
    >>>

    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        software: "sentieon:202112.05"
    }

    # 工具运行输出结果
    output {
        File deduped_bam = out_bam
        File deduped_bam_bai = out_bai
        Array[File] qc_metrics = glob("*_metrics.txt")
    }

}



task Sentieon_BQSR{
    input{
        File ref_dir
        File dbsnp_dir
        File dbmills_dir

        String sample_id
        String ref_fasta
        String dbsnp
        String db_mills

        File deduped_bam
        File deduped_bam_index

        File regions

        Int cpu = 32
        String memory = "64G"
        String disks = "local-disk 250 cloud_ssd"
        #String SENTIEON_LICENSE = "172.25.164.226:8990"
        #String docker = "registry.cn-shanghai.aliyuncs.com/hcc1395_aliyun/sentieon-genomics:v202112.05"
    }

    command<<<
        set -o pipefail
        set -exo
        mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
        export bwt_max_mem="$((mem_kb / 1024 / 1024 - 2))g"
        export MALLOC_CONF=lg_dirty_mult:-1
        export LD_PRELOAD=/usr/local/sentieon-genomics/lib/libjemalloc.so
        
        nt=$(nproc)

        sentieon driver -t $nt \
        -r ~{ref_dir}/~{ref_fasta} -i ~{deduped_bam} \
        --interval ~{regions} \
        --algo QualCal \
        -k ~{dbsnp_dir}/~{dbsnp} -k ~{dbmills_dir}/~{db_mills} \
        ~{sample_id}_recal_data.table

        sentieon driver -t $nt \
        -r ~{ref_dir}/~{ref_fasta} -i ~{deduped_bam} \
        -q ~{sample_id}_recal_data.table \
        --algo QualCal \
        -k ~{dbsnp_dir}/~{dbsnp} -k ~{dbmills_dir}/~{db_mills} \
        ~{sample_id}_recal_data.table.post \
        --algo ReadWriter ~{sample_id}.sorted.deduped.recaled.bam

        sentieon driver -t $nt --algo QualCal \
        --plot --before ~{sample_id}_recal_data.table --after ~{sample_id}_recal_data.table.post ~{sample_id}_recal_data.csv

        sentieon plot bqsr -o ~{sample_id}_bqsrreport.pdf ~{sample_id}_recal_data.csv
        
        sentieon driver -t $nt \
        --interval ~{regions} -r ~{ref_dir}/~{ref_fasta} \
        -i ~{sample_id}.sorted.deduped.recaled.bam \
        --algo Haplotyper -d ~{dbsnp_dir}/~{dbsnp} \
        ~{sample_id}.Haplotyper.vcf
    >>>

    runtime{
        cpu:cpu
        memory:memory
        disks:disks
        #docker:docker
        software: "sentieon:202112.05"
    }

    output{
        File recal_table = "~{sample_id}_recal_data.table"
        File recal_post = "~{sample_id}_recal_data.table.post"
        File recaled_bam = "~{sample_id}.sorted.deduped.recaled.bam"
        File recaled_bam_index = "~{sample_id}.sorted.deduped.recaled.bam.bai"
        File recal_csv = "~{sample_id}_recal_data.csv"
        File bqsrreport_pdf = "~{sample_id}_bqsrreport.pdf"
        File vcf = "~{sample_id}.Haplotyper.vcf"
        File vcf_idx = "~{sample_id}.Haplotyper.vcf.idx"
    }
}