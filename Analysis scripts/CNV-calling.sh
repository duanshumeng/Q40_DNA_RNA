source activate ascat_env
cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Data_Performance/HCC1395_WES/ASCAT

# The allele distribution within the targeted region of WES was obtained first
Rscript ascat-WES-preTarget.r hcc1395-WES-preTarget.txt /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Validation_by_WES_element/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_noy_cut.bed /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Data_Performance/HCC1395_WES/ASCAT hg38

# Copy number variation analysis
Rscript ascat-WES.r hcc1395-WES.txt /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Validation_by_WES_element/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_noy_cut.bed /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Data_Performance/HCC1395_WES/ASCAT hg38

#Extract copy number
for i in `ls *.Rdata`; do Rscript ascat_grep.r $i ASCAT_csv ${i/.ASCAT_objects.Rdata/.ASCAT}; done

