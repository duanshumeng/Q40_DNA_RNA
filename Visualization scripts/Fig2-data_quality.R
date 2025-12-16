library(ggunchained)
source('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/Quartet_Multiomics_Ratio_Code/utils/theme_nature.r')
colors = c('#4472CA','#E69F00')
names(colors) <- c('Q40','Q30')
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/Data_quality')
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(gghalves)

colors = c('#3783BB','#91569F','#C94741')
names(colors) <- c('Q30','Q35','Q40')

my_theme <- theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 0),
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.x = element_text(color = "black",face = "bold",size = 16),
        axis.title.y = element_text(color = "black",face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16),
        legend.text = element_text(face = "bold",size = 14,colour = "black"),
        legend.position = "bottom",legend.background = element_blank(),
        panel.grid.minor  = element_line(colour = NA),
        panel.grid.major.x   = element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA)
  )+
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(colour = "black",face = "bold",size = 16)) + 
  theme(strip.background.y = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.y = element_text(colour = "black",face = "bold",size = 16)) + 
  theme(strip.placement = "inside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))

shape_values = c(21,22,24,23)
names(shape_values) = c('HCC1395/BL','NA12878','Quartet','MAQC')

jitter_plot <- function(df,d_type,y_lab,seq_type,ylim,y_position){
  print(ylim[2]-1)
  if (d_type != 'Duplicates'){
    df['Duplicates']= df[d_type]
  }
  p1 <- ggplot(df,aes(x=Quality,y=Duplicates,color = Quality,fill=Quality))+
    geom_half_violin(
      position = position_nudge(x = 0.1, y = 0),
      side = 'R', adjust = 1.2, trim = T,alpha=0.5
    )+
    geom_half_point(
      aes(fill=Quality,color = Quality),
      position = position_nudge(x = -0.35, y = 0), 
      size = 3, range_scale = 0.9,alpha=0.5,shape = 21
    )+
    geom_boxplot(
      outlier.shape = NA,  # hide outliers
      width = 0.1,alpha=0.5
    )+
    geom_signif(comparisons = list(c("Q30","Q40")),# groups to compare
                map_signif_level = T, # show significance with *
                test = t.test, ## method
                textsize = 10,
                y_position = y_position,
                size=0.8,color="black")+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+scale_shape_manual(values = shape_values)+
    my_theme+labs(y=y_lab,x=seq_type)+
    ylim(ylim)
  return(p1)
  
}

## DNA -----
# Duplicates
fastqc <- read.csv("WES_fastQC/multiqc_fastqc.txt",sep = '\t')
fastqc.stats <- read.csv("WES_fastQC/multiqc_general_stats.txt",sep='\t')
dup.df = fastqc.stats[,c('Sample','FastQC_mqc.generalstats.fastqc.percent_duplicates')]
dup.df$Quality <- lapply(dup.df$Sample,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
dup.df$source <- dup.df$Sample
dup.df$Sample <- str_split(dup.df$source,'\_',simplify = T)[,1]
dup.df$Sample[dup.df$Sample %in% c('D5','D6','F7','M8')] = 'Quartet'
dup.df$Sample <- gsub('HCC1395BL','HCC1395',dup.df$Sample)
colnames(dup.df) <- c('Sample','Duplicates','Quality','source')
dup.df$Sample <- gsub('HCC1395','HCC1395/BL',dup.df$Sample) %>% gsub('NIST','NA12878',.)
dup.df <- subset(dup.df,Quality !='Q35')
dup.dna <- jitter_plot(dup.df,'Duplicates','Duplication rates (%)','WES',c(10,40),c(37,38));dup.dna

dup.df %>%
  group_by(Quality) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(Duplicates, na.rm = TRUE),   # calculate mean
    SD_Duplicates = sd(Duplicates, na.rm = TRUE)       # calculate sd
  )

# Base Quality (ILM & ELE)
bq.df <- read.csv('WES_fastQC/fastqc_per_base_sequence_quality_plot.tsv',sep = '\t')
bq.df.sub <- colMeans(bq.df[,-1]) %>% as.data.frame()
bq.df.sub$batch <- lapply(rownames(bq.df.sub),function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()

bq.df.sub$Sample <- str_split(rownames(bq.df.sub),'\_',simplify = T)[,1]
bq.df.sub$Sample[bq.df.sub$Sample %in% c('D5','D6','F7','M8')] = 'Quartet'
bq.df.sub$Sample <- gsub('HCC1395BL','HCC1395',bq.df.sub$Sample)
bq.df.sub$Sample <- gsub('HCC1395','HCC1395/BL',bq.df.sub$Sample) %>% gsub('NIST','NA12878',.)
colnames(bq.df.sub) <- c('Base_quality','Quality','Sample')

bq.df.sub$Read <- lapply(rownames(bq.df.sub),function(x){
  p = ifelse(grepl('R1',x),'R1','R2')
  return(p)
}) %>% unlist()

size_values <- c(5,3)
names(size_values) <- c('R1','R2')
jitter_plot.bq <- function(df,d_type,y_lab,seq_type){
  if (d_type != 'Duplicates'){
    df['Duplicates']= df[d_type]
  }
  p1 <- ggplot(df,aes(x=Quality,y=Duplicates,color = Quality))+
    geom_half_violin(aes(fill=Quality),
                     position = position_nudge(x = 0.1, y = 0),
                     side = 'R', adjust = 1.2, trim = FALSE,alpha=0.5
    )+
    geom_half_point(
      aes(size = Read,fill=Quality,alpha=0.1),
      position = position_nudge(x = -0.35, y = 0), 
      range_scale = 0.9,shape=21)+
    geom_boxplot(aes(fill=Quality),
                 outlier.shape = NA,  # hide outliers
                 width = 0.1,alpha=0.5
    )+
    geom_signif(comparisons = list(c("Q30","Q40")),# groups to compare
                map_signif_level = T, # show significance with *
                test = t.test, ## method
                textsize = 10,
                y_position = c(43),# y-position of the bar
                size=0.8,color="black")+
    scale_color_manual(values = colors)+scale_size_manual(values = size_values)+
    scale_fill_manual(values = colors)+scale_shape_manual(values = shape_values)+
    my_theme+labs(y=y_lab,x=seq_type)+ylim(c(35,44))
  return(p1)
}

bq.df.sub %>%
  group_by(Quality) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(Base_quality, na.rm = TRUE),   # calculate mean
    SD_Duplicates = sd(Base_quality, na.rm = TRUE)       # calculate sd
  )
# proportion of data reaching Q40 (97.61905 %)
sum(round(subset(bq.df.sub,Quality=='Q40')$value) >= 40)/(length(subset(bq.df.sub,Quality=='Q40')$value)) * 100
bq.df.sub <- subset(bq.df.sub,Quality !='Q35')
bq.dna <- jitter_plot.bq(bq.df.sub,"Base_quality",'Phred quality score','WES');bq.dna

## RNA --------
# Duplicates ----------
# Duplicates (ELE & ILM)
fastqc.stats.rna <- read.csv("RNAseq_fastQC/multiqc_general_stats.txt",sep='\t')
dup.df.rna = fastqc.stats.rna[,c('Sample','FastQC_mqc.generalstats.fastqc.percent_duplicates')]
dup.df.rna$Quality <- lapply(dup.df.rna$Sample,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
  return(p)
}) %>% unlist()
dup.df.rna$source <- dup.df.rna$Sample
dup.df.rna$Sample <- str_split(dup.df.rna$source,'\_',simplify = T)[,1]
dup.df.rna$Sample[dup.df.rna$Sample %in% c('D5','D6','F7','M8','M81D63','M83D61')] = 'Quartet'
dup.df.rna$Sample[dup.df.rna$Sample %in% c('A','B')] = 'MAQC'
colnames(dup.df.rna) <- c('Sample','Duplicates','Quality','source')
dup.df.rna <- subset(dup.df.rna,Quality !='Q35')
dup.rna <- jitter_plot(dup.df.rna,'Duplicates','Duplication rates (%)','RNA-seq',c(30,48),c(46,46));dup.rna

dup.df.rna %>%
  group_by(Quality) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(Duplicates, na.rm = TRUE),   # calculate mean
    SD_Duplicates = sd(Duplicates, na.rm = TRUE)       # calculate sd
  )

# Base quality -----------
# Base quality (ELE & ILM)
bq.df.rna <- read.csv('RNAseq_fastQC/fastqc_per_base_sequence_quality_plot.tsv',sep = '\t')
bq.df.sub.rna <- colMeans(bq.df.rna[,-1]) %>% as.data.frame()
bq.df.sub.rna$batch <- lapply(rownames(bq.df.sub.rna),function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()

bq.df.sub.rna$Sample <- str_split(rownames(bq.df.sub.rna),'\_',simplify = T)[,1]
bq.df.sub.rna$Sample[bq.df.sub.rna$Sample %in% c('D5','D6','F7','M8','M81D63','M83D61')] = 'Quartet'
bq.df.sub.rna$Sample[bq.df.sub.rna$Sample %in% c('A','B')] = 'MAQC'

colnames(bq.df.sub.rna) <- c('Base_quality','Quality','Sample')

bq.df.sub.rna$Read <- lapply(rownames(bq.df.sub.rna),function(x){
  p = ifelse(grepl('R1',x),'R1','R2')
  return(p)
}) %>% unlist()

bq.df.sub.rna <- subset(bq.df.sub.rna,Quality !='Q35')
bq.rna <- jitter_plot.bq(bq.df.sub.rna,"Base_quality",'Phred quality score','RNA-seq');bq.rna

# Bam --------
# DNA --------
# ILM & ELE
mq.df <- read.csv('WES_bamQC/multiqc_qualimap_bamqc_genome_results.txt',sep = '\t')
mq.df$percentage_aligned <- mq.df$mapped_reads/mq.df$total_reads*100
mq.df.sub <- mq.df[,c('Sample','mean_coverage','mean_mapping_quality','percentage_aligned')]
mq.df.sub$Quality <- lapply(mq.df.sub$Sample,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
mq.df.sub$source <- mq.df.sub$Sample

mq.df.sub$Sample <- str_split(mq.df.sub$source,'\_',simplify = T)[,1]
mq.df.sub$Sample[mq.df.sub$Sample %in% c('D5','D6','F7','M8','M81D63','M83D61')] = 'Quartet'
mq.df.sub$Sample <- gsub('HCC1395BL','HCC1395',mq.df.sub$Sample)
mq.df.sub$Sample <- gsub('NIST','NA12878',mq.df.sub$Sample) %>% gsub('HCC1395','HCC1395/BL',.)

mq.df.sub <- subset(mq.df.sub ,Quality !='Q35')
mq.dna <- jitter_plot(mq.df.sub,'mean_mapping_quality','Mapping quality (%)','WES',c(58.35,58.55),c(58.49,58.45));mq.dna

aln.dna <- jitter_plot(mq.df.sub,'percentage_aligned','Aligned reads (%)','WES',c(99.85,100),c(99.97,99.95));aln.dna

# mean_mapping_quality
mq.df.sub %>%
  group_by(Quality) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(mean_mapping_quality, na.rm = TRUE),   # calculate mean
    SD_Duplicates = sd(mean_mapping_quality, na.rm = TRUE)       # calculate sd
  )

# percentage_aligned
mq.df.sub %>%
  group_by(Quality) %>%
  dplyr::summarise(
    Mean_percentage_aligned = mean(percentage_aligned, na.rm = TRUE) %>% round(2),   # calculate mean
    SD_percentage_aligned = sd(percentage_aligned, na.rm = TRUE)       # calculate sd
  )

#Aligned_reads proportion
# ELE & ILM--
mq.df.rna <- read.csv('RNAseq_RNAseqQC/qualimap_rnaseq_genome_results.txt',sep = '\t')
mq.df.sub.rna <- mq.df.rna[,c('Sample','not_aligned','total_alignments','X5_3_bias')]
mq.df.sub.rna$Quality <- lapply(mq.df.sub.rna$Sample,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
mq.df.sub.rna$source <- mq.df.sub.rna$Sample

mq.df.sub.rna$Sample <- str_split(mq.df.sub.rna$source,'\_',simplify = T)[,1]
mq.df.sub.rna$Sample[mq.df.sub.rna$Sample %in% c('D5','D6','F7','M8','M81D63','M83D61')] = 'Quartet'
mq.df.sub.rna$Sample[mq.df.sub.rna$Sample %in% c('A','B')] = 'MAQC'
mq.df.sub.rna$percentage_aligned <- ((mq.df.sub.rna$total_alignments - mq.df.sub.rna$not_aligned)/mq.df.sub.rna$total_alignments)*100 %>% round(2)
mq.df.sub.rna <- subset(mq.df.sub.rna,Quality !='Q35')
aln.rna <- jitter_plot(mq.df.sub.rna,'percentage_aligned','Aligned reads (%)','RNA-seq',c(96,98.5),c(98.1,97.8));aln.rna

mq.df.sub.rna[,c('Sample','X5_3_bias','Quality')]
mq.df.sub.rna$bias_5_3 <- mq.df.sub.rna$X5_3_bias - 1
mq.df.sub.rna$source <- paste0(mq.df.sub.rna$Quality,'',mq.df.sub.rna$Sample,'',mq.df.sub.rna$source) %>% 
  gsub('_ELE.10G','',.) %>% gsub('_ILM.10G','',.)
mq.df.sub.rna <- mq.df.sub.rna[order(mq.df.sub.rna$source), ]
mq.df.sub.rna$bias_type <- ifelse(mq.df.sub.rna$bias_5_3 > 0,'bias_5','bias_3')

size_values.bias = c(5,3)
names(size_values.bias) <- c('bias_5','bias_3')

mq.df.sub.rna$bias_abs <- abs(mq.df.sub.rna$bias_5_3)

bias_rna <- jitter_plot(mq.df.sub.rna,'bias_abs','5/3 bias','RNA-seq',c(0,0.12),c(0.105,0.089));bias_rna

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}

rna_legend <- g_legend(bq.rna)  %>% plot()
dna_legend <- g_legend(bq.dna) %>% plot()

Fig.QC.v1.png = ((bq.dna+dup.dna)|(mq.dna+aln.dna))/((bq.rna+dup.rna)|(aln.rna+bias_rna))+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & 
  theme(legend.position='',                       # remove legend
        plot.tag = element_text(color = "black",size = 18,face = "bold"))

# save main figure
ggsave('Fig.QC_main.pdf',  Fig.QC.v1.png, width = 8.15*1.5, height = 5.7*1.5, dpi = 300)
ggsave('Fig.QC_main.tiff', Fig.QC.v1.png, width = 8.15*1.5, height = 5.7*1.5, dpi = 300)

# supplementary figure
Fig.QC_sup.png = ((aln.rna + bias_rna) / (mq.dna + aln.dna)) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect') & 
  theme(legend.position = '',                     # remove legend
        plot.tag = element_text(color = "black", size = 18, face = "bold"))

ggsave('Fig.QC_sup.pdf',  Fig.QC_sup.png, width = 8.15*1.5, height = 5.7*1.5, dpi = 300)
ggsave('Fig.QC_sup.tiff', Fig.QC_sup.png, width = 8.15*1.5, height = 5.7*1.5, dpi = 300)