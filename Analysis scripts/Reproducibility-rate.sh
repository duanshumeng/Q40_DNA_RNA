#NIST
source activate hcc1395
lst=('30X' '60X' '90X' '120X' '25X' '20X' '15X' '10X')
for i in ${lst[@]}; do echo $i; mkdir Haplotyper_WES_${i}; cp Haplotyper_WES/*${i}.Haplotyper.vcf* Haplotyper_WES_${i}; done

nohup bash JI_batch.sh &
dirs=('Haplotyper_WES_120X/' 'Haplotyper_WES_90X/' 'Haplotyper_WES_60X/' 'Haplotyper_WES_30X/' 'Haplotyper_WES_25X/' 'Haplotyper_WES_20X/' 'Haplotyper_WES_15X/' 'Haplotyper_WES_10X/')
for i in ${dirs[@]}
do
echo $i
#python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/Jaccard_bedtools.py -i ${i} -o JI -s no -c no
python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/NIST_F1_and_Reproducibility.py -i $i -o ${i}Rep -s no -R /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/NIST_high_confidence/HG001_GRCh38_1_22_v4.2.1_benchmark.spikein.sorted.bed
done


#HCC1395
cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Data_Performance/HCC1395_WES
for i in ${lst[@]}; do echo $i; mkdir TNseq_${i}; cp TNseq/*${i}.TNseq.vcf* TNseq_${i}; done

nohup bash Rep_batch.sh &
dirs=('TNseq_120X/' 'TNseq_90X/' 'TNseq_60X/' 'TNseq_30X/')
for i in ${dirs[@]}
do
echo $i
python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/HCC1395_F1_and_Reproducibility.py -i $i -o ${i}Rep  -s no -R /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Validation_by_WES_element/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_noy_cut.HCR.bed
done

#Quartet
#D5/D6/F7/M8 

lst=('30X' '60X' '90X' '120X' '25X' '20X' '15X' '10X')
for i in ${lst[@]}; do echo $i; mkdir D5_${i}; cp *${i}.Haplotyper.vcf* D5_${i}; done
for i in ${lst[@]}; do echo $i; mkdir D6_${i}; cp *${i}.Haplotyper.vcf* D6_${i}; done
for i in ${lst[@]}; do echo $i; mkdir F7_${i}; cp *${i}.Haplotyper.vcf* F7_${i}; done
for i in ${lst[@]}; do echo $i; mkdir M8_${i}; cp *${i}.Haplotyper.vcf* M8_${i}; done
lst=('30X' '60X' '90X' '120X' '25X' '20X' '15X' '10X')
smp=('D5' 'D6' 'F7' 'M8')
for i in ${lst[@]}
do
           for j in ${smp[@]}
           do 
           echo $i
           python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/Quartet_F1_and_Reproducibility.py -i ${j}_${i} -o ${j}_${i}Rep  -s no -R /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Validation_by_WES_element/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_noy_cut.HCR.bed
done
