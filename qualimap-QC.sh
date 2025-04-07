#qualimap 质控
for i in `ls *.bam`
do
        echo $i
        qualimap bamqc -bam $i -gff /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Fastq/Illumina/Fastq_raw/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_noy_cut.bed -outformat PDF:HTML -nt 24 -nr 500 -nw 1500 -outdir $i.bamqc --java-mem-size=64G
done