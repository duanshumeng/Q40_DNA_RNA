#构建索引：
source activate happy_env
for i in `ls *.vcf`; do echo $i; bcftools view -Oz $i -o ${i}.gz; bcftools index -t ${i}.gz;done

#合并vcf
lst=('30X' '60X' '90X' '120X' '25X' '20X' '15X' '10X')
quality=('Q30' 'Q40')
for i in ${lst[@]}
do
for j in ${quality[@]}
do
echo ${j}_${i}
rtg vcfmerge --force-merge-all -o MCR/${j}_1.${i}.family.vcf.gz ${j}_D5_1.${i}.Haplotyper.vcf.gz ${j}_D6_1.${i}.Haplotyper.vcf.gz ${j}_F7_1.${i}.Haplotyper.vcf.gz ${j}_M8_1.${i}.Haplotyper.vcf.gz
rtg vcfmerge --force-merge-all -o MCR/${j}_2.${i}.family.vcf.gz ${j}_D5_2.${i}.Haplotyper.vcf.gz ${j}_D6_2.${i}.Haplotyper.vcf.gz ${j}_F7_2.${i}.Haplotyper.vcf.gz ${j}_M8_2.${i}.Haplotyper.vcf.gz
rtg vcfmerge --force-merge-all -o MCR/${j}_3.${i}.family.vcf.gz ${j}_D5_3.${i}.Haplotyper.vcf.gz ${j}_D6_3.${i}.Haplotyper.vcf.gz ${j}_F7_3.${i}.Haplotyper.vcf.gz ${j}_M8_3.${i}.Haplotyper.vcf.gz
done
done

#运行vbt
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/MGI_Q40/MGI_WES/Call.germlinevcf/Quartet_VCF/MCR

for i in `ls *.family.vcf.gz`
do
echo $i
family_name=${i/.family.vcf.gz/}
family_vcf=${i}
rm -r VBT_D5
rm -r VBT_D6
## 指定测试数据批次
nt=16
d5=$(zcat ${i} | grep '#CHROM' | awk -F ' ' '{print $10}')
d6=$(zcat ${i} | grep '#CHROM' | awk -F ' ' '{print $11}')
f7=$(zcat ${i} | grep '#CHROM' | awk -F ' ' '{print $12}')
m8=$(zcat ${i} | grep '#CHROM' | awk -F ' ' '{print $13}')
## 构建D5、F7、M8 ped文件
echo -e "${family_name}\t${m8}\t0\t0\t2\t-9\n${family_name}\t${f7}\t0\t0\t1\t-9\n${family_name}\t${d5}\t${f7}\t${m8}\t2\t-9" > ${family_name}.D5.ped

mkdir VBT_D5
## 计算D5家系孟德尔遗传率
/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/public/software/VBT-TrioAnalysis/vbt mendelian -ref /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Reference/GRCh38.d1.vd1.fa -mother ${family_vcf} -father ${family_vcf} -child ${family_vcf} -pedigree ${family_name}.D5.ped -outDir VBT_D5 -out-prefix ${family_name}.D5 --output-violation-regions -thread-count $nt

## 文件重命名
cat VBT_D5/${family_name}.D5_trio.vcf > ${family_name}.D5.vcf

## 构建D6、F7、M8 ped文件

echo -e "${family_name}\t${m8}\t0\t0\t2\t-9\n${family_name}\t${f7}\t0\t0\t1\t-9\n${family_name}\t${d6}\t${f7}\t${m8}\t2\t-9" > ${family_name}.D6.ped

mkdir VBT_D6
## 计算D6家系的孟德尔遗传率
/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/public/software/VBT-TrioAnalysis/vbt mendelian -ref /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Reference/GRCh38.d1.vd1.fa -mother ${family_vcf} -father ${family_vcf} -child ${family_vcf} -pedigree ${family_name}.D6.ped -outDir VBT_D6 -out-prefix ${family_name}.D6 --output-violation-regions -thread-count $nt

cat VBT_D6/${family_name}.D6_trio.vcf > ${family_name}.D6.vcf

## 文件重命名
cat ${family_name}.D5.vcf | grep -v '##' > ${family_name}.D5.txt
cat ${family_name}.D6.vcf | grep -v '##' > ${family_name}.D6.txt
zcat ${family_vcf} | grep -v '##' | awk '
                BEGIN { OFS = "\t" }
                NF > 2 && FNR > 1 { 
                    for ( i=9; i<=NF; i++ ) { 
                        split($i,a,":") ;$i = a[1];
                    } 
                } 
                { print }
' > ${family_name}.consensus.txt

python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Software/Quartet-DNA-QC-report/merge_two_family_with_genotype.py -LCL5 ${family_name}.D5.txt  -LCL6 ${family_name}.D6.txt -genotype ${family_name}.consensus.txt -family ${family_name}
done

