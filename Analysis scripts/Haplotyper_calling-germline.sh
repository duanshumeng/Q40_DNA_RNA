export PATH="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/miniconda3/bin:$PATH"
source activate mapping_V2
cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/MGI_Q40/MGI_WES
input_path="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/MGI_Q40/MGI_WES/Bam_downsample"
echo $input_path
smp=$(ls -d $input_path/* | grep -v 'Q30_HCC1395\|Q40_HCC1395' | grep 'bam$')
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
    bash shell_scripts/Call_SNV_INDEL.germline.sh -i $normal  -o Call.germline -f /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Reference/bwa_index/GRCh38.d1.vd1.fa -d /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Software/annovar &
done
wait
echo "================END================"

# 运行GATK HaplotypeCaller
gatk HaplotypeCaller \
    -R "${ref_fasta}" \
    -I "${normal_bam}" \
    --native-pair-hmm-threads ${nt} \
    --dbsnp ${dbsnp} \
    -O ${out_dir}/${Normal_name}.Haplotyper.vcf
