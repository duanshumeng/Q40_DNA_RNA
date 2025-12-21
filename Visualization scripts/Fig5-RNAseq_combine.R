library(edgeR)
library(cowplot)
library(scales)
library(RColorBrewer)
library(grDevices)
library(grid)
library(ggh4x)
library(patchwork)
library(dplyr)
library(stringr)
library(ggpattern)
library(gg.gap)

source('~/PGx_lab/HCC1395/ElementVSIllumina/文章/Scripts/Visualization scripts/Statistics_analysis.R')
#colors = c('#3783BB','#91569F','#C94741')
#names(colors) <- c('Q30','Q35','Q40')

#colors = c('#3783BB','#C94741')

colors=c("#709AE1FF","#F05C3BFF")
names(colors) <- c('Q30','Q40')

color_gradient <- colorRampPalette(c("#709AE1FF", "#FFFFFF")) 
colors <- color_gradient(6)
color_gradient <- colorRampPalette(c("#F05C3BFF", "#FFFFFF")) 
colors <- color_gradient(6)

rna_color <- c("#709AE1","#8CAEE7","#A9C2ED","#C5D6F3","#E2EAF9","#F05C3B","#F37C62","#F69D89","#F9BDB0", "#FCDED7")
names(rna_color) <- c('Q30-10','Q30-8','Q30-6','Q30-4','Q30-2','Q40-10','Q40-8','Q40-6','Q40-4','Q40-2')

smp_colors = c('#F2E3EB' ,'#D6DEFF' ,'#FFEDC2')
names(smp_colors) <- c('Quartet','MAQC','ERCC')

my_theme_1 <- theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 0),
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.x = element_text(color = "black",face = "bold",size = 16),
        axis.title.y = element_text(color = "black",face = "bold",size = 16),
        #axis.ticks.x = element_text(color = "black",face = "bold",size = 14),
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

my_theme <- theme_bw()+
  theme(
    axis.text.x = element_text(face = "bold",color = "black"),
    axis.text.y = element_text(face = "bold",color = "black"),
    axis.title  = element_text(face = "bold",color = "black"),
    panel.grid.minor  =  element_blank(),
    legend.background =  element_blank(),
    legend.title = element_text(face = "bold",color = "black"),
    legend.text = element_text(color = "black")
  )
setwd('/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/Data_performance_V2/RNAseq/')

load('MAQC/MAQC.rnaseq.RData')
load('Quartet/Quartet.rnaseq.RData')
load('ERCC/ERCC.rnaseq.RData')


#SNR---------
snr_df.maqc$Sample <- 'MAQC'
snr_df.quartet$Sample <- 'Quartet'
snr_df <- rbind(snr_df.maqc,snr_df.quartet)
snr_df$SNR <- snr_df$. %>% as.numeric() %>% round(.,2)
snr_df$Platform <- str_split(rownames(snr_df),'\\ ',simplify = T)[,1]
snr_df$Quality <- sub('ILM','Q30',snr_df$Platform) %>% sub('ELE','Q40',.)  %>% sub('Mix','Q35',.)
snr_df$Size <- str_split(rownames(snr_df),'\\ ',simplify = T)[,2]
snr_df$Size <- sub('G1|G','',snr_df$Size) %>% as.numeric()
snr_df$group <- paste0(snr_df$Quality,' ',snr_df$Sample)
snr_df <- subset(snr_df,Platform!='MGI')
snr_df <- subset(snr_df,Quality !='Q35')
snr_df$color <- paste0(snr_df$Quality,'-',snr_df$Size)
library(ggplot2)
library(dplyr)
library(colorspace)


get_bar_pattern <- function(df,d_type,y_title,y_lim=c(0,1),x_lim=seq(0,10,2),smp_colors = c('#F2E3EB' ,'#D6DEFF' ,'#FFEDC2')){
  df$SNR <- df[,d_type]
  p=ggplot(df, aes(x = Size, y = SNR,fill=color)) +
    #geom_bar_pattern(alpha=0.5,stat = "identity", width=1.5,
                     #position=position_dodge(width=1.5), colour = 'black', pattern_fill = 'white') +
    geom_bar(alpha=0.8,stat = "identity", width=1.5,position=position_dodge(width=1.5), colour = 'black')+
    facet_grid2(.~Sample,scales = "free",space = "free_x",
                strip  = strip_nested(
                  background_x = elem_list_rect(fill = smp_colors,color = NA)))+
    theme(legend.key.size = unit(1.5, 'cm'))+
    scale_fill_manual(values = rna_color)+
    labs(x="Size (Gb)",y=y_title)+
    scale_x_continuous(breaks=seq(0,10,2))+my_theme
  return(p)
  } 


snr.p <- get_bar_pattern(df=snr_df,d_type='SNR',y_title='SNR');snr.p

subset(snr_df,Sample=='Quartet') %>%
  group_by(Quality) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(SNR, na.rm = TRUE),  
    SD_Duplicates = sd(SNR, na.rm = TRUE)      
  )

jitter_plot <- function(df,d_type,y_lab,seq_type,ylim,y_position){
  print(ylim[2]-1)
  if (d_type != 'Duplicates'){
    df['Duplicates']= df[d_type]
  }
  p1 <- ggplot(df,aes(x=Quality,y=Duplicates,fill=Quality))+
    #geom_violin(cex=1)+
    geom_half_point(
      aes(fill=Quality),
      position = position_nudge(x = -0.35, y = 0), 
      size = 3,color='white', range_scale = 0.9,shape = 21
    )+
    geom_boxplot(
      outlier.shape = NA,  
      width = 0.5,alpha=0.5
    )+
    #geom_boxplot(cex=1)+
    #geom_jitter(
    #  aes(shape = Sample,fill=Quality), 
    #  position = position_jitter(0.1),
    #  size = 5,alpha=0.5,
    #) +
    geom_signif(comparisons = list(c("Q30","Q40")),
                map_signif_level = T, 
                test = t.test, 
                textsize = 5,
                y_position = y_position,
                #y_position = c(ylim[2]-1),
                #tip_length = c(c(0.05,0.05)),
                size=0.8,color="black")+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+
    my_theme+labs(y=y_lab,x=seq_type)+
    ylim(ylim)
  return(p1)
  
}
snr.p.quartet <- jitter_plot(subset(snr_df,Sample=='Quartet'),'SNR','SNR','Quartet',c(25,38),c(35,37));snr.p.quartet
snr.p.maqc <- jitter_plot(subset(snr_df,Sample=='MAQC'),'SNR','SNR','MAQC',c(18,21),c(20,22));snr.p.maqc

Fig.5d.snr.p.quartet.test <- get_test_type_with_effectsize_RNAseq(subset(snr_df,Sample=='Quartet'))
Fig.S4d.snr.p.maqc.test <- get_test_type_with_effectsize_RNAseq(subset(snr_df,Sample=='MAQC'))


g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}
legend <- g_legend(snr.p)  %>% plot()


##PCC-----------
pcc_df.ercc <- pcc_df.ercc %>% unlist() %>% as.data.frame()
pcc_df.ercc$Sample <- 'ERCC'
pcc_df.maqc$Sample <- 'MAQC'
pcc_df.quartet$Sample <- 'Quartet'


pcc_df <- rbind(pcc_df.maqc,pcc_df.quartet) %>% rbind(.,pcc_df.ercc)
pcc_df$PCC <- pcc_df$. %>% as.numeric() 
pcc_df$Platform <- str_split(rownames(pcc_df),'\\ ',simplify = T)[,1]
pcc_df$Quality <- sub('ILM','Q30',pcc_df$Platform) %>% sub('ELE','Q40',.)  %>% sub('Mix','Q35',.)
pcc_df$Size <- str_split(rownames(pcc_df),'\\ ',simplify = T)[,2]
pcc_df$Size <- sub('G.cor1|G.cor|G','',pcc_df$Size) %>% as.numeric()
pcc_df$group <- paste0(pcc_df$Quality,' ',pcc_df$Sample)
pcc_df <- subset(pcc_df,Platform!='MGI')
#pcc_df$PCC <- round(pcc_df$PCC,3)
pcc_df <- subset(pcc_df,Quality !='Q35')
pcc_df$color <- paste0(pcc_df$Quality,'-',pcc_df$Size)
subset(pcc_df,Sample=='Quartet') %>%
  group_by(Quality) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(PCC, na.rm = TRUE),   
    SD_Duplicates = sd(PCC, na.rm = TRUE)       
  )

pcc.p <- get_bar_pattern(df=pcc_df,d_type='PCC',y_title='PCC',smp_colors=c('#FFEDC2','#F2E3EB' ,'#D6DEFF'));pcc.p

pcc.p1 = gg.gap(plot=pcc.p,
                  segments=c(0.25,0.8),
                  tick_width = c(0.25,0.1),
                  #rel_heights=c(1),
                  ylim=c(0,1))+theme(legend.position="bottom");pcc.p1

pcc.p.quartet <- jitter_plot(subset(pcc_df,Sample=='Quartet'),'PCC','PCC','Quartet',c(0.85,0.95),c(0.91,0.92));pcc.p.quartet
pcc.p.maqc <- jitter_plot(subset(pcc_df,Sample=='MAQC'),'PCC','PCC','MAQC',c(0.9,0.95),c(0.94,0.94));pcc.p.maqc
pcc.p.ercc <- jitter_plot(subset(pcc_df,Sample=='ERCC'),'PCC','PCC','ERCC',c(0.95,0.98),c(0.97,0.97));pcc.p.ercc

pcc_df$SNR <- pcc_df$PCC
Fig.5b.pcc.p.quartet.test <- get_test_type_with_effectsize_RNAseq(subset(pcc_df,Sample=='Quartet'))
Fig.S4a.pcc.p.maqc.test <- get_test_type_with_effectsize_RNAseq(subset(pcc_df,Sample=='MAQC'))
Fig.S4b.pcc.p.ercc.test <- get_test_type_with_effectsize_RNAseq(subset(pcc_df,Sample=='ERCC'))


get_lineplot.MCC <- function(df,d_type,y_title,y_lim=c(0,1),x_lim=seq(0,10,2),smp_colors=c('#F2E3EB' ,'#D6DEFF' ,'#D6DEFF','#D6DEFF')){
  df$SNR <- df[,d_type]
  p1 <- ggplot(df, aes(x = Size, y = SNR,fill=color)) +
    #geom_bar_pattern(alpha=0.5,stat = "identity", width=1.5,
    #position=position_dodge(width=1.5), colour = 'black', pattern_fill = 'white') +
    geom_bar(alpha=0.8,stat = "identity", width=1.5,position=position_dodge(width=1.5), colour = 'black')+
    facet_grid2(.~Sample,scales = "free",space = "free_x",
                strip  = strip_nested(
                  background_y = elem_list_rect(fill =  
                                                  # '#4472CA','#E69F00'
                                                  c('#D7EAE4', '#9FE5D2'
                                                    # ,'#FFDD33' ,'#527AFF','#FFDD33' ,'#527AFF'
                                                    #'#E69F00','#4472CA','#E69F00','#4472CA'
                                                  ),color = NA),
                  background_x = elem_list_rect(fill = smp_colors,color = NA)))+
    scale_fill_manual(values = rna_color)+
    scale_x_continuous(breaks =x_lim)+
    theme_bw()+
    my_theme+
    labs(x="Size (G)",y=y_title)+
    #expand_limits(y = y_lim)+
    guides(fill=guide_legend(title = "",nrow = 1,byrow = FALSE))
  
  p1;return(p1)
  
}
#MCC----------------
MCC.df.maqc$Source <- MCC.df.maqc$Sample
colnames(MCC.df.maqc) <- sub('F1.score','MCC',colnames(MCC.df.maqc))
MCC.df.maqc$Sample <- 'MAQC'
MCC.df.maqc$Compare <- 'A/B'

#MCC.df.quartet$Source <- MCC.df.quartet$Sample
colnames(MCC.df.quartet) <- sub('MCC.score','MCC',colnames(MCC.df.quartet))
MCC.df.quartet$Sample <- 'Quartet'
MCC.df <- rbind(MCC.df.maqc,MCC.df.quartet)
MCC.df$Platform <- str_split(MCC.df$Source,'\\ ',simplify = T)[,1]
MCC.df$Size <- str_split(MCC.df$Source,'\\ ',simplify = T)[,2]
MCC.df$Size <- sub('G','',MCC.df$Size) %>% as.numeric()
MCC.df$Quality <- sub('ILM','Q30',MCC.df$Platform) %>% sub('ELE','Q40',.)  %>% sub('Mix','Q35',.)
MCC.df$Sample <- paste0(MCC.df$Sample,' ',MCC.df$Compare)
MCC.df$group <- paste0(MCC.df$Quality,' ',MCC.df$Sample)
MCC.df <- subset(MCC.df,Platform!='MGI')
MCC.df <- subset(MCC.df,Quality !='Q35')
MCC.df$color <- paste0(MCC.df$Quality,'-',MCC.df$Size)

subset(MCC.df,grepl('Quartet',MCC.df$Sample)) %>%
  group_by(Quality) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(MCC, na.rm = TRUE),   
    SD_Duplicates = sd(MCC, na.rm = TRUE)       
  )

MCC.p <- get_lineplot.MCC(df=MCC.df,d_type='MCC',y_title='MCC');MCC.p

MCC.p1 = gg.gap(plot=MCC.p,
                segments=c(0.1,0.5),
                tick_width = c(0.1,0.1),
                #rel_heights=c(1),
                ylim=c(0,0.7))+theme(legend.position="bottom");MCC.p1

MCC.df$smp <- str_split(MCC.df$Sample,'\\ ',simplify = T)[,1]
mcc.p.quartet <- jitter_plot(subset(MCC.df,smp=='Quartet'),'MCC','MCC','Quartet',c(0.6,0.8),c(0.75,0.76));mcc.p.quartet
mcc.p.maqc <- jitter_plot(subset(MCC.df,smp=='MAQC'),'MCC','MCC','MAQC',c(0.6,0.8),c(0.75,0.76));mcc.p.maqc

MCC.df$SNR <- MCC.df$MCC
rownames(MCC.df) <- paste0(MCC.df$Source,MCC.df$Platform,MCC.df$Sample)

Fig.5f.mcc.p.quartet.test <- get_test_type_with_effectsize_RNAseq(subset(MCC.df,smp=='Quartet'))
Fig.S4c.mcc.p.maqc.test <- get_test_type_with_effectsize_RNAseq(subset(MCC.df,smp=='MAQC'))

#Save test results---------------
library(openxlsx)
test_file <- ls(pattern = '*Fig*') %>% sort()
#test_file <- test_file[-grep("get_test_type", test_file)]
wb <- createWorkbook()

for (var_name in test_file) {
  print(var_name)
  df <- get(var_name)
=
  sheet_name <- substr(var_name, 1, 31) 

  sheet_name <- gsub("[\\[\\]\\*\\?\\/\\\\]", "_", sheet_name)

  addWorksheet(wb, sheet_name)
  
  writeData(wb, sheet = sheet_name, df)
  
  cat("Add sheet:", sheet_name, "\n")
}

saveWorkbook(wb, "RNAseq_Statistic_test_type.xlsx", overwrite = TRUE)



rna_main.p1 <- (pcc.p1+theme(axis.title.x = element_blank()))/(MCC.p1+theme(axis.title.x = element_blank()))/snr.p+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect',heights = c(2,2,1)) & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 16,face = "bold"))

rna_main.p1

rna_main.p2 <- (pcc.p.quartet+theme(axis.title.x = element_blank()))+(mcc.p.quartet+theme(axis.title.x = element_blank()))+snr.p.quartet+
  rna_main.p2+plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='',
                                                           plot.tag = element_text(color = "black",size = 18,face = "bold"))
rna_main.p2
ggsave("Fig.rna.main-1.pdf",rna_main.p1,width=8.15*1.5, height=5.7*2.5,dpi = 300)
ggsave("Fig.rna.main-2.pdf",rna_main.p2,width=8.15*1.5, height=5.7*0.7,dpi = 300)


#Gene count----
snr_gene_df.maqc$Sample <- 'MAQC'
snr_gene_df.quartet$Sample <- 'Quartet'
snr_gene_df <- rbind(snr_gene_df.maqc,snr_gene_df.quartet)
snr_gene_df$Count <- snr_gene_df$. %>% as.numeric() 
snr_gene_df$Platform <- str_split(rownames(snr_gene_df),'\\ ',simplify = T)[,1]
snr_gene_df$Quality <- sub('ILM','Q30',snr_gene_df$Platform) %>% sub('ELE','Q40',.)  %>% sub('Mix','Q35',.)
snr_gene_df$Size <- str_split(rownames(snr_gene_df),'\\ ',simplify = T)[,2]
snr_gene_df$Size <- sub('G1|G','',snr_gene_df$Size) %>% as.numeric()
snr_gene_df$group <- paste0(snr_gene_df$Quality,' ',snr_gene_df$Sample)
snr_gene_df <- subset(snr_gene_df,Platform!='MGI')
snr_gene_df <- subset(snr_gene_df,Quality !='Q35')
snr_gene_df$color <- paste0(snr_gene_df$Quality,'-',snr_gene_df$Size)
count.p <- get_bar_pattern(df=snr_gene_df,d_type='Count',y_title='Gene Count');count.p

count.p1 = gg.gap(plot=count.p,
                  segments=c(10000,20000),
                  tick_width = c(20000,5000),
                  #rel_heights=c(1),
                  ylim=c(0,35000))+theme(legend.position="bottom");count.p1


#CV--------------

get_boxplot <- function(df,d_type,y_title,y_lim=c(0,1),smp_colors=c('#F2E3EB' ,'#D6DEFF' ,'#FFEDC2')){
  df$SNR <- df[,d_type]
  p1 <- ggplot(data=df,aes(x=Size,y=SNR)) +
    geom_boxplot(cex=0.5,alpha=0.7,aes(fill=color),outlier.size = 0,outlier.stroke = 0)+
    facet_grid2(.~Sample,scales = "free",space = "free_x",
                strip  = strip_nested(
                  background_y = elem_list_rect(fill =  
                                                  # '#4472CA','#E69F00'
                                                  c('#D7EAE4', '#9FE5D2'
                                                    # ,'#FFDD33' ,'#527AFF','#FFDD33' ,'#527AFF'
                                                    #'#E69F00','#4472CA','#E69F00','#4472CA'
                                                  ),color = NA),
                  background_x = elem_list_rect(fill = smp_colors,color = NA)))+
    scale_color_manual(values = rna_color)+
    scale_fill_manual(values = rna_color)+
    theme_bw()+
    my_theme+
    labs(x="Size (G)",y=y_title)+xlim(c('2','4','6','8','10'))+ylim(y_lim)+
    guides(fill=guide_legend(title = "",nrow = 1,byrow = FALSE))
  
  p1;return(p1)
  
}

df_cv.maqc$Sample <- 'MAQC'
df_cv.maqc <- df_cv.maqc[-1,]
df_cv.quartet$Sample <- 'Quartet'
df_cv.quartet <- df_cv.quartet[-1,]
cv_df <- rbind(df_cv.maqc,df_cv.quartet)
cv_df$CV <- cv_df$value %>% as.numeric() 
cv_df$Platform <- str_split(cv_df$Quality,'\\ ',simplify = T)[,1]
cv_df$Size <- str_split(cv_df$Quality,'\\ ',simplify = T)[,2]
cv_df$Size <- sub('G|G1','',cv_df$Size) %>% as.character()
cv_df$Quality <- sub('ILM','Q30',cv_df$Platform) %>% sub('ELE','Q40',.)  %>% sub('Mix','Q35',.)
cv_df <- subset(cv_df,Platform!='MGI')
cv_df <- subset(cv_df,Quality != 'Q35')
cv_df$color <- paste0(cv_df$Quality,'-',cv_df$Size)
cv.p <- get_boxplot(cv_df,d_type='CV',y_title='CV',y_lim=c(0,15));cv.p

#ERCC AUC---------
auc_df.ercc$Sample <- 'ERCC'
auc_df.ercc= auc_df.ercc[-1,]
auc_df.ercc$Platform <- str_split(auc_df.ercc$source,'\\ ',simplify = T)[,1]
auc_df.ercc$Size <- str_split(auc_df.ercc$source,'\\ ',simplify = T)[,2]
auc_df.ercc$Size <- sub('G','',auc_df.ercc$Size) %>% as.numeric()
auc_df.ercc$Quality <- sub('ILM','Q30',auc_df.ercc$Platform) %>% sub('ELE','Q40',.)  %>% sub('Mix','Q35',.)
auc_df.ercc$Sample <- paste0(auc_df.ercc$Sample,' ',auc_df.ercc$Ratio)
auc_df.ercc$group <- paste0(auc_df.ercc$Quality,' ',auc_df.ercc$Sample)
auc_df.ercc <- subset(auc_df.ercc,Platform!='MGI')
auc_df.ercc <- subset(auc_df.ercc,Quality != 'Q35')
auc_df.ercc$color <- paste0(auc_df.ercc$Quality,'-',auc_df.ercc$Size)
auc.p <- get_lineplot.MCC(df=auc_df.ercc,d_type = 'AUC',y_title='AUC',smp_colors=c('#FFEDC2' ,'#FFEDC2' ,'#FFEDC2'));auc.p


#ERCC LODR----------
lodr_df.ercc$Sample <- 'ERCC'
lodr_df.ercc = lodr_df.ercc[-1,]
lodr_df.ercc$Platform <- str_split(lodr_df.ercc$source,'\\ ',simplify = T)[,1]
lodr_df.ercc$Size <- str_split(lodr_df.ercc$source,'\\ ',simplify = T)[,2]
lodr_df.ercc$Size <- sub('G','',lodr_df.ercc$Size) %>% as.numeric()
lodr_df.ercc$Quality <- sub('ILM','Q30',lodr_df.ercc$Platform) %>% sub('ELE','Q40',.)  %>% sub('Mix','Q35',.)
lodr_df.ercc$Sample <- paste0(lodr_df.ercc$Sample,' ',lodr_df.ercc$Ratio)
lodr_df.ercc$group <- paste0(lodr_df.ercc$Quality,' ',lodr_df.ercc$Sample)
lodr_df.ercc$Estimate <- sub('<','',lodr_df.ercc$Estimate)
lodr_df.ercc$Estimate <- lodr_df.ercc$Estimate %>% as.numeric()
lodr_df.ercc$Estimate[is.infinite(lodr_df.ercc$Estimate)] <- 1000
lodr_df.ercc <- subset(lodr_df.ercc,Platform!='MGI')
lodr_df.ercc <- subset(lodr_df.ercc,Quality != 'Q35')
lodr_df.ercc$color <- paste0(lodr_df.ercc$Quality,'-',lodr_df.ercc$Size)
lodr.p <- get_lineplot.MCC(df=lodr_df.ercc,d_type = 'Estimate',y_title='LODR',smp_colors=c('#FFEDC2' ,'#FFEDC2' ,'#FFEDC2'));lodr.p

rna_extend.p = (pcc.p.maqc|pcc.p.ercc|mcc.p.maqc|snr.p.maqc)/((auc.p/lodr.p)|(count.p/cv.p))+plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect',heights = c(1,2)) & theme(legend.position = '',
                                        plot.tag = element_text(color = "black",size = 16,face = "bold"))


rna_extend.p
ggsave("Fig.rna.sup.pdf",rna_extend.p,width=8.15*1.5, height=5.7*1.5,dpi = 300)
ggsave("Fig.rna.sup.tiff",rna_extend.p,width=8.15, height=5.7*1.5,dpi = 300)


## WES---
wes_qc_index <- read.csv('/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/Data_quality/WES.QC.index.summary.csv')
rnaseq_qc_index <- read.csv('/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/Data_quality/RNAseq.QC.index.summary.csv')

wes_index <- read.csv('/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/Data_performance/WES.index.summary.csv')
cnv_index <- read.csv('/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/Data_performance/CNV.index.summary.csv')
rnaseq_index <- read.csv('/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/Data_performance_V2/RNAseq/RNAseq_index_summary.csv')
cnv_index$Sample <- 'HCC1395/BL'
colnames(cnv_index)[4]='cnPCC'
wes_index <- merge(wes_index,cnv_index,by=c('Quality','Sample'),all=TRUE)
wes_qc_index$Sample <- gsub('HCC1395','HCC1395/BL',wes_qc_index$Sample) %>% gsub('NIST','NA12878',.)

wes_index.all <- merge(wes_index,wes_qc_index,by=c('Quality','Sample'),all=TRUE)



library(tidyverse)
library(patchwork)



library(ggradar)
normalize_custom_range <- function(x, new_min, new_max) {
  return (((x - min(x)) / (max(x) - min(x))) * (new_max - new_min) + new_min)
}



index_radar <- function(idx_df,idx_type,split_smp=FALSE){
  radar_p.rnaseq <- list()
  if (split_smp){
    for (i in unique(idx_df$Sample)){
      print(i)
      df <- subset(idx_df,Sample==i)
      df_clean <- df[, colSums(is.na(df)) < nrow(df)]
      df_clean <- df_clean[,-c(2)]
      df_clean[,-1] <- data.frame(lapply(df_clean[,-1], function(x) normalize_custom_range(x, 0,1)))
      print(df_clean)
      df.p <- ggradar(df_clean, 
                      background.circle.transparency = 0,
                      group.colours =colors,
                      grid.min=0,
                      grid.max=1,
                      plot.title =i,
                      fill =T,
                      fill.alpha = 0.3,
                      group.point.size = 5)
      radar_p.rnaseq[[i]] <- df.p
    } 
  } else {
    df <- idx_df
    df_clean <- df[, colSums(is.na(df)) < nrow(df)]
    df_clean <- df_clean[,-c(1)]
    df_clean[,-1] <- data.frame(lapply(df_clean[,-1], function(x) normalize_custom_range(x, 0,1)))
    df.p <- ggradar(df_clean, 
                    background.circle.transparency = 0,
                    group.colours =colors[-2],
                    grid.min=0,
                    grid.max=1,
                    plot.title =idx_type,
                    fill =T,
                    fill.alpha = 0.3,
                    group.point.size = 5)
    radar_p.rnaseq[[idx_type]] <- df.p
  }
  return(radar_p.rnaseq)
}

#WES 
wes_index.all$Dup <- 1/wes_index.all$Dup
wes_index.all <- subset(wes_index.all,Quality != 'Q30_M' & Quality != 'Q40_M')
wes_index.all$Quality <- str_split(wes_index.all$Quality,'\\_',simplify = T)[,1]
wes_index.all <- subset(wes_index.all,Quality != 'Q30_M' & Quality != 'Q40_M')
wes_index.all$Quality <- gsub('_I','',wes_index.all$Quality) %>% gsub('_E','',.) 
wes_radar <- index_radar(wes_index.all,'WES',split_smp=TRUE);wes_radar

wes_index.all <- subset(wes_index.all, select = -Reg)


#qc_index

T1_df_avg <- wes_index.all %>%
  group_by(Quality) %>%
  dplyr::summarise(across(F1_S:Ali, ~mean(. , na.rm = TRUE)))

T1_df_long <- T1_df_avg %>%
  pivot_longer(cols = F1_S:Ali, names_to = "Feature", values_to = "Value")

ggplot(T1_df_long, aes(x = Feature, y = Quality, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") + 
  theme_minimal() +
  labs(title = "Heatmap of Feature Values by Quality", 
       x = "Feature", 
       y = "Quality", 
       fill = "Mean Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#RNAseq QC
rnaseq_qc_index$Dup <- 1/rnaseq_qc_index$Dup
rnaseq_qc_index$Bias <- 1/abs(rnaseq_qc_index$Bias-1)
rnaseq_qc_index$Bias[is.infinite(rnaseq_qc_index$Bias)] <- 300
rnaseq_qc_index <- subset(rnaseq_qc_index,Quality != 'Q30_M')
rnaseq_qc_index$Quality <- str_split(rnaseq_qc_index$Quality,'\\_',simplify = T)[,1]
ranseq_qc_radar <- index_radar(rnaseq_qc_index,'RNAseq-QC');ranseq_qc_radar

#RNAseq performance
rnaseq_index <- rnaseq_index[,c('Sample','Quality','SNR','PCC','MCC','AUC','LODR')]
rnaseq_index$LODR <- 1/rnaseq_index$LODR

rnaseq_index <- subset(rnaseq_index,Quality != 'Q40_M')
rnaseq_index$Quality <- str_split(rnaseq_index$Quality,'\\_',simplify = T)[,1]
ranseq_radar <- index_radar(rnaseq_index,'WES-PM',split_smp=TRUE);ranseq_radar

#WES performance
wes_index <- wes_index[,c('Sample','Quality','F1_G','F1_S','Rep','MCR','PCC','F1_L','Reg','JI','cnPCC')]
wes_index <- subset(wes_index,Quality != 'Q30_M' & Quality != 'Q40_M')
wes_index$Quality <- str_split(wes_index$Quality,'\\_',simplify = T)[,1]
wes_radar <- index_radar(wes_index,'WES-PM',split_smp=TRUE);wes_radar


CNV_Q40_all.main <- (CNV_Q40_JI.all|CNV_Q40_PCC.all)/(ranseq_qc_radar[[1]]|wes_qc_radar[[1]])/(ranseq_radar[['Quartet']]|wes_radar[['HCC1395/BL']])+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect',heights = c(1.5,0.8,0.8)) & theme(legend.position='bottom',
                                                                 plot.tag = element_text(color = "black",size = 18,face = "bold"))
ggsave('Fig.CNV.main.pdf',CNV_Q40_all.main,width=8.15*1.5, height=5.7*1.5,dpi = 300)

index_sup <- (ranseq_radar[['MAQC']]|ranseq_radar[['ERCC']])/(wes_radar[['Quartet']]|wes_radar[['NA12878']])+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect',heights = c(1,1)) & theme(legend.position='bottom',
                                                         plot.tag = element_text(color = "black",size = 18,face = "bold"))
index_sup

ggsave('Fig.index.sup.pdf',index_sup,width=8.15, height=5.7,dpi = 300)


