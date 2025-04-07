export PATH="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/miniconda3/bin:$PATH"
source activate mapping_V2
cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/MGI_Q40/MGI_WES
input_path="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/MGI_Q40/MGI_WES/Bam_downsample"
smp=$(ls -d $input_path/* | grep 'HCC1395BL' | grep -v 'bai')
for i in ${smp[@]}
do
    n=$(basename $i)
    #echo $n
    t=${n/HCC1395BL/HCC1395}
    #echo $t
    normal=$input_path/$n
    tumor=$input_path/$t
    #ls -l $normal
    #ls  -l $tumor
    bash shell_scripts/Call_SNV_INDEL.sh -i $normal -I $tumor -o Call.variants -f /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Reference/bwa_index/GRCh38.d1.vd1.fa -d /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Software/annovar &
done
wait
echo "================END================"

# 运行GATK Mutect2
gatk Mutect2 \
    -R "${ref_fasta}" \
    -I "${tumor_bam}" \
    -I "${normal_bam}" \
    --native-pair-hmm-threads ${nt} \
    -normal ${Normal_sample} \
    --tumor-sample ${Tumor_sample} \
    -O ${out_dir}/${Tumor_name}/${Tumor_name}.vcf

## 过滤vcf文件
gatk FilterMutectCalls \
    -V ${out_dir}/${Tumor_name}/${Tumor_name}.vcf \
    -R ${ref_fasta} \
    -O ${out_dir}/${Tumor_name}/${Tumor_name}.filter.vcf