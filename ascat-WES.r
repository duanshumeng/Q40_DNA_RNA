library(ASCAT)
library(getopt)

Sys.getenv("R_HOME")

rhome = Sys.getenv("R_HOME")
print('R_HOMR')
print(rhome)
args <- commandArgs()
print(args)
rhome <- args[1]
worksheet <- args[6]
bed_file = args[7]
workdir = args[8]
hg_version = args[9]
#setwd("/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Data_Performance/HCC1395_WES/ASCAT")
#worksheet="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Data_Performance/HCC1395_WES/ASCAT/myWorksheet.txt"
#tumor_bam = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Data_Performance/HCC1395_WES/WES_bams/WES_HCC1395_3_ILM_10G.30X.bam"
#normal_bam= "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Data_Performance/HCC1395_WES/WES_bams/WES_HCC1395BL_3_ILM_10G.30X.bam"
#tumor_name="WES_HCC1395_3_ILM_10G.30X"
#normal_name="WES_HCC1395BL_3_ILM_10G.30X"
#sample_name="WES_HCC1395_3_ILM_10G.30X"
#data_type="WES"
#output_dir="./"
#gender = "XX"
#hg_version = "hg38"
#bed = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Validation_by_WES_element/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_noy_cut.bed"
loci_prefix <- "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Software/ASCAT/WES/G1000_lociAll_hg38_chr/G1000_loci_hg38_chr"
alleles_prefix <- "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Software/ASCAT/WES/G1000_allelesAll_hg38/G1000_alleles_hg38_chr"
gc_file <- '/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Software/ASCAT/WES/GC_G1000_hg38.txt'
rt_file <- '/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Software/ASCAT/WES/RT_G1000_hg38.txt'
allelecounter ="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/miniconda3/envs/ascat_env/bin/alleleCounter" 


#ascat.prepareTargetedSeq(
#  Worksheet = worksheet, # A tab-separated file with specific information. Check format using ?ascat.prepareTargetedSeq
#  Workdir = workdir,
#  alleles.prefix = alleles_prefix,
#  BED_file = bed_file,
#  allelecounter_exe = allelecounter,
#  genomeVersion = hg_version,
#  is_chr_based=T,
#  nthreads = 24)
smp_info <- read.table(worksheet,sep='\t',header=T)
head(smp_info)

for (i in rownames(smp_info)){
    
    normal_name = smp_info[i,'Normal_ID']
    normal_bam = smp_info[i,'Normal_file']
    tumor_name = smp_info[i,'Tumor_ID']
    tumor_bam = smp_info[i,'Tumor_file']
    gender = smp_info[1,'Gender']
    hg_version = smp_info[1,'hg_version']
    print(paste0('======Processing==',normal_name,"====",tumor_name,"====="))
      ascat.prepareHTS(
      tumourseqfile = tumor_bam,
      normalseqfile = normal_bam,
      tumourname = tumor_name,
      normalname = normal_name,
      allelecounter_exe = allelecounter,
      alleles.prefix = paste0(workdir,"/alleleData/Cleaned/alleleData_chr"),
      loci.prefix =paste0(workdir, "/alleleData/Cleaned/loci_chr"),
      gender = gender,
      genomeVersion = hg_version,
      #chrom_names = paste0('chr',c(1:22, "X")),
      nthreads = 24,
      tumourLogR_file = paste0(tumor_name, "_LogR.txt"),
      tumourBAF_file = paste0(tumor_name, "_BAF.txt"),
      normalLogR_file = paste0(normal_name, "_LogR.txt"),
      normalBAF_file = paste0(normal_name, "_BAF.txt"))

    print('======ascat.prepareHTS==========Finished.......')

    ascat.bc = ascat.loadData(Tumor_LogR_file = paste0(tumor_name, "_LogR.txt"),
                              Tumor_BAF_file = paste0(tumor_name, "_BAF.txt"),
                              Germline_LogR_file = paste0(normal_name, "_LogR.txt"),
                              Germline_BAF_file = paste0(normal_name, "_BAF.txt"), 
                              gender = gender, genomeVersion = hg_version, isTargetedSeq=T)
    print('=======ascat.bc======Finished.......')

    ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
    ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_file, replictimingfile = rt_file)
    ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
    ascat.bc = ascat.aspcf(ascat.bc, penalty=25)
    ascat.plotSegmentedData(ascat.bc)
    ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
    QC = ascat.metrics(ascat.bc,ascat.output)
    save(ascat.bc, ascat.output, QC, file =paste0(tumor_name,'.ASCAT_objects.Rdata'))
    print('=======All======Finished.......')
    
    }


pattern <- "alleleFrequencies"
files <- list.files()
files_to_remove <- files[grep(pattern, files)]
file.remove(files_to_remove)
