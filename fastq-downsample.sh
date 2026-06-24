cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/MGI_Q40/MGI_WES/Fastq_10G
for i in `cat trimed_fastq.txt`
do echo $i
smp=$(basename $i)
echo ${smp/trimmed/10G}
seqtk sample -s 12345 $i 33300000 > ${smp/trimmed/10G}
done
wait
echo "seqtk has beed finished..."


#RNAseq fastq to 10G
cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Fastq/RNAseq/Fastp
for i in `ls *.clean.fastq.gz`
do 
echo  $i
smp=$i
seqtk sample -s 12345 $i 33300000 > ${smp/clean/10G}
done
