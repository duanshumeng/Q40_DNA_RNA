#source activate mapping

input_dir="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Fastq/RNAseq/Fastp"

nt=24
idx="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/ERCC_RNA/ercc_reference/SRM2374_Sequence_v1"
pen_noncansplice=0
pen_cansplice=0
pen_intronlen="G,-8,1"
min_intronlen=20
max_intronlen=500000
maxins=500
minins=0

insert_size=8000

for i in `ls $input_dir | grep 'D6\|M8' | grep ".fastq.gz" | grep R1`
do
    echo $i
    r1=$(basename $i)
    echo $r1
    r2=${r1/R1/R2}
    echo $r2
    sample_id=$(echo $r1 | awk -F '.' '{print $1}')
    echo $sample_id
    set -o pipefail
    set -e
    echo "hisat2 -t -p $nt -x ${idx} -1 ${input_dir}/${r1} -2 ${input_dir}/${r2}  -S $sample_id.sam --pen-cansplice ${pen_cansplice} --pen-noncansplice ${pen_noncansplice} --pen-canintronlen ${pen_intronlen} --min-intronlen ${min_intronlen} --max-intronlen ${max_intronlen} --maxins ${maxins} --minins ${minins} --rna-strandness RF --un-conc-gz ${sample_id}_un.fq.gz"
    hisat2 -t -p $nt -x ${idx} -1 ${input_dir}/${r1} -2 ${input_dir}/${r2}  -S $sample_id.sam --pen-cansplice ${pen_cansplice} --pen-noncansplice ${pen_noncansplice} --pen-canintronlen ${pen_intronlen} --min-intronlen ${min_intronlen} --max-intronlen ${max_intronlen} --maxins ${maxins} --minins ${minins} --rna-strandness RF --un-conc-gz ${sample_id}_un.fq.gz
       
    
    samtools view -bS ${sample_id}.sam > ${sample_id}.bam
    samtools sort -m 1000000000 ${sample_id}.bam -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    samtools view -bs 42.1 ${sample_id}.sorted.bam > ${sample_id}.sorted.percent.bam
    samtools stats -i ${insert_size} ${sample_id}.sorted.bam |grep ^IS|cut -f 2- > ${sample_id}.ins_size
    
done

echo "fastq2bam finished!!!!"