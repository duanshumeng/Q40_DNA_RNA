#source activate mapping

input_dir="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/ERCC_RNA/Bam_files"
minimum_isoform_abundance=0.01
minimum_length_allowed_for_the_predicted_transcripts=200
maximum_fraction_of_muliplelocationmapped_reads=0.95
Junctions_no_spliced_reads=10
gtf="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/ERCC_RNA/ercc_reference/SRM2374_ambiguities_resolved_v1.GTF"
mkdir ballgown
nt=24
for i in `ls $input_dir/*.bam`
do
    echo $i
    bam=$(basename $i)
    echo $bam
    sample_id=$(echo $bam | awk -F '.' '{print $1}')
    echo $sample_id
    stringtie -e \
    -B \
    -p $nt \
    -f ${minimum_isoform_abundance} \
    -m ${minimum_length_allowed_for_the_predicted_transcripts} \
    -a ${Junctions_no_spliced_reads} \
    -M ${maximum_fraction_of_muliplelocationmapped_reads} \
    -G ${gtf} \
    --rf \
    -o ballgown/${sample_id}/${sample_id}.gtf \
    -C ${sample_id}.cov.ref.gtf \
    -A ${sample_id}.gene.abundance.txt \
    $input_dir/${bam}

done

cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/ERCC_RNA/Stringtie_quant
for i in `ls -d */`; do echo "${i/\//} ${i}${i/\//}.gtf" >> sample_list.txt; done

awk -F ' ' -v OFS='\t' '{print $1,$2}' sample_list.txt > sample_list_1.txt
python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/miniconda3/envs/rnaseq_mapping/bin/prepDE.py \
-i sample_list.txt  \
-g gene_count_matrix.csv  \
-t transcript_count_matrix.csv