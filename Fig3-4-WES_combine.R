rm(list = ls())
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/Data_performance')
source('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/Quartet_Multiomics_Ratio_Code/utils/theme_nature.r')
library(ggplot2)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(stringr)
library(ggpubr)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggunchained) 
library(ggh4x)
library(ggpattern)
library(rstatix)
my_theme <- theme_bw()+
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

my_theme_2 <- theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x =element_text(color = "black",size = 10,face = "bold",angle = 0),
        axis.text.y = element_text(color = "black",size = 10),
        axis.title.x = element_text(color = "black",face = "bold",size = 10),
        axis.title.y = element_text(color = "black",face = "bold",size = 10),
        #axis.ticks.x = element_text(color = "black",face = "bold",size = 14),
        legend.title = element_text(face = "bold",size = 10),
        legend.text = element_text(face = "bold",size = 10,colour = "black"),
        legend.position = "bottom",legend.background = element_blank(),
        panel.grid.minor  = element_line(colour = NA),
        panel.grid.major.x   = element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA)
  )+
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(colour = "black",face = "bold",size = 10)) + 
  theme(strip.background.y = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.y = element_text(colour = "black",face = "bold",size = 10)) + 
  theme(strip.placement = "inside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))
colors = c('#3783BB','#91569F','#C94741')
names(colors) <- c('Q30','Q35','Q40')

smp_colors = c('#F2E3EB' ,'#D6DEFF' ,'#FFEDC2')
names(smp_colors) <- c('NA12878','Quartet','HCC1395/BL')
#变

get_boxcompare = function(tnseq.df.sub,m_type, y_lim,ct_1395.sum=data.frame(),d_type='F1.score'){
  tnseq.df.sub <- subset(tnseq.df.sub,type == m_type)
  tnseq.df.sub['F1.score'] <- tnseq.df.sub[d_type]
  tnseq.df.test <- tnseq.df.sub %>%
    group_by(type, Depth) %>%
    filter(n_distinct(F1.score) > 1, n_distinct(Quality) > 1) %>% 
    t_test(F1.score ~ Quality, p.adjust.method = 'none', comparisons = c("Q30", "Q40")) %>%
    mutate(p.adj.signif = p.adjust(p, method = "none")) %>%  # Adjust p-value and calculate significance
    add_significance() %>%  # Adds significance symbols based on adjusted p-values
    select(-df, -statistic, -p)
  tnseq.df.test
  tnseq.df.test <- tnseq.df.test %>% add_xy_position(x = "Depth")
  Max_F1_score <- tnseq.df.sub %>%
    group_by(type, Depth) %>%
    dplyr::summarise(y.position = max(F1.score, na.rm = TRUE)+0.01)
  print(Max_F1_score)
  print(tnseq.df.test)
  tnseq.df.test<- merge(tnseq.df.test,Max_F1_score,by=c('type','Depth'))
  tnseq.df.test$y.position <- tnseq.df.test$y.position.y
  
  #tnseq.df.test$y.position <- ifelse(tnseq.df.test$group2=='Q40',tnseq.df.test$y.position-0.01,tnseq.df.test$y.position-0.03)
  tnseq.df.sub$group <- paste0(tnseq.df.sub$Quality,'_',tnseq.df.sub$Depth)
  box_p <- ggplot(tnseq.df.sub, aes(x = Depth, y = F1.score)) +
    geom_rect(xmin = 0.5, xmax = 1.5,
              ymin = 0, ymax = 11, fill ="#E9E9E9") +
    geom_rect(xmin = 2.5, xmax = 3.5,
              ymin = 0, ymax = 11, fill ="#E9E9E9")+
    geom_rect(xmin = 4.5, xmax = 5.5,
              ymin = 0, ymax = 2000000, fill ="#E9E9E9") +
    geom_rect(xmin = 6.5, xmax = 7.5,
              ymin = 0, ymax = 2000000, fill ="#E9E9E9")+
    geom_boxplot(aes(fill=Quality,color=Quality),position = position_dodge(width = 0.8),
                 outlier.shape = NA,  # 隐藏离群点
                 width = 0.6,alpha=0.3,size = 0.5
    )+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                size = 3, alpha=0.5,aes(color = Quality))+
    #geom_line(aes(color = Quality,group=Quality),position = position_dodge(width = 0.5),alpha=0.5,size=1)+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+
    #facet_grid(type~.)+
    ylim(y_lim)+ylab(paste0(d_type,' ',m_type))+ggtitle(unique(tnseq.df.sub$Sample))+
    my_theme+
    stat_pvalue_manual(
      tnseq.df.test, label = "p.signif", 
      step.increase = 0.01,size = 5
    )+xlab('')
  if (length(ct_1395.sum) > 0){
    ct_1395_plot <- ggplot(subset(ct_1395.sum,type==m_type), aes(x = Depth, y = count_median,fill=Quality)) +
      geom_rect(xmin = 0.5, xmax = 1.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9") +
      geom_rect(xmin = 2.5, xmax = 3.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9")+
      geom_rect(xmin = 4.5, xmax = 5.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9") +
      geom_rect(xmin = 6.5, xmax = 7.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9")+
      geom_bar(stat = "identity", position = "dodge", width = 0.7)+
      geom_errorbar(aes(ymin = count_median - count_sd, ymax = count_median + count_sd), position = position_dodge(0.7), width = 0.2)+
      scale_fill_manual(values = colors)+
      my_theme_2+theme(legend.position = "none")+ylab('Count')+xlab('')
    if (unique(tnseq.df.sub$Sample)=='HCC1395/BL'){
      box_p = box_p +annotation_custom(grob=ggplotGrob(ct_1395_plot),ymin = y_lim[1], ymax=y_lim[1]+0.2, xmin=2, xmax=4.5)
    } else {
      box_p = box_p +annotation_custom(grob=ggplotGrob(ct_1395_plot),ymin = y_lim[1], ymax=y_lim[1]+0.3, xmin=3, xmax=8.5)
    }
  }

  return(box_p)
}
#变异位点数量统计-----
#ELE & ILM
df_count <- read.csv('ELE_ILM/WES_VCF_counts.vcf_stats.csv')
df_count$Quality <- lapply(str_split(df_count$source,'\\_',simplify = T)[,4],function(x){
  p = ifelse(grepl('ILM',x),'Q30_I','Q40_E')
  return(p)
}) %>% unlist()

df_count$Sample <- str_split(df_count$source,'\\_',simplify = T)[,2]
df_count$Sample[df_count$Sample %in% c('D5','D6','F7','M8')] = 'Quartet'
df_count$Source <- paste0(df_count$Sample,'.',df_count$Quality,'.',df_count$source)
df_count$Type <- gsub('SNVs','SNV',df_count$type) %>% gsub('Indels','INDEL',.)
df_count<- subset(df_count,source !='WES_D6_2_ELE_10G.120X.Haplotyper')
df_count$Quality <- as.factor(df_count$Quality)
df_count$Sample <- as.factor(df_count$Sample)
df_count$Rep <- paste0('Rep',str_split(df_count$source,'\\_',simplify = T)[,3])
df_count$Quality <- str_split(df_count$Quality,'\\_',simplify = T)[,1]

df_count_median <- df_count %>%
  dplyr::group_by(Type,Sample, Quality) %>%
  dplyr::summarize(
    across(
      count,
      list(
        mean = ~mean(.),
        sd = ~sd(.)
      )
    )
  )

y_title = 'Variant counts'
df_count.p <- ggplot(data=df_count_median,aes(x=Quality,y=count_mean,fill=Quality)) +
  geom_bar(stat="summary",fun=mean,position=position_dodge(0.9),alpha=0.8)+
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd),color='black',width=.2,
                position=position_dodge(0.9))+
  stat_summary(geom = "text", aes(label = paste(round(..y.., digits = 0))),
               position = position_dodge(width = 1), vjust = -0.9,hjust=0.5,size = 3,color='black')+
  facet_grid2(Type~Sample,scales = "free",space = "free_x",
              strip  = strip_nested(
                background_y = elem_list_rect(fill =  
                                                # '#4472CA','#E69F00'
                                                c('#D7EAE4', '#9FE5D2'
                                                  #,'#FFDD33' ,'#527AFF','#FFDD33' ,'#527AFF'
                                                ),color = NA),
                background_x = elem_list_rect(fill = 
                                                c('#F2E3EB' ,'#D6DEFF' ,'#FFEDC2',
                                                  '#F2E3EB','#F2E3EB',
                                                  '#D6DEFF','#D6DEFF',
                                                  '#FFEDC2','#FFEDC2'),color = NA)))+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)),
                     labels = scales::unit_format(unit = "K", scale = 1e-3))+
  my_theme+
  labs(x="",y=y_title)+
  #expand_limits(y = y_lim)+
  guides(fill=guide_legend(title = "",nrow = 1,byrow = FALSE));df_count.p


#不同深度下的变异位点准确性F1--------------
get_dotbar <- function(df,d_type,y_title,y_lim=c(0,1),xlab="Depth (X)"){
  df[,'mean'] = df[,paste0(d_type,'_mean')]
  df[,'sd'] = df[,paste0(d_type,'_sd')]
  #if (length(unique(df$Sample)) > 1){
  #  SNV_Q40 <-   ggplot(data=df,aes(x=Depth,y=mean,pattern=Sample,fill=Quality)) +geom_bar_pattern(stat="identity",position=position_dodge(0.8),width = 0.8,color='black',alpha=0.8)
  #} else {
    SNV_Q40 <- ggplot(data=df,aes(x=Depth,y=mean,fill=Quality))+geom_bar(stat="identity",position=position_dodge(0.8),width = 0.8,color='black',alpha=0.8)
  #}
  
  SNV_Q40 <- SNV_Q40+ geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),color='black',width=0.2,
                                    position=position_dodge(0.8))+

    
    facet_grid2(Sample~Type,scales = "free",space = "free_x",
                strip  = strip_nested(
                  background_x = elem_list_rect(fill =  
                                                  # '#4472CA','#E69F00'
                                                  c('#D7EAE4', '#9FE5D2'
                                                    # ,'#FFDD33' ,'#527AFF','#FFDD33' ,'#527AFF'
                                                    #'#E69F00','#4472CA','#E69F00','#4472CA'
                                                  ),color = NA),
                  background_y = elem_list_rect(fill = smp_colors,color = NA)))+
    scale_fill_manual(values = colors)+
    my_theme+
    labs(x=xlab,y=y_title)+
    expand_limits(y = y_lim)+
    guides(fill=guide_legend(title = "",nrow = 1,byrow = FALSE))
  if (xlab == "Depth (X)"){
    xbreaks <- c(0, 5, 10, 15, 20, 25, 30, 60, 90, 120)
    SNV_Q40 <- SNV_Q40+scale_x_discrete(breaks = xbreaks)
  } else {
    SNV_Q40 <- SNV_Q40+theme(axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 30,hjust = 1))
  }
  SNV_Q40;return(SNV_Q40)
}

#高置信区间-------
#NIST-NA12878
#ELE & ILM
nist.df <- read.csv("ELE_ILM/F1score_stats_total_NIST.csv")
nist.df.sub <- nist.df[,c('Type','source','METRIC.F1_Score')]
nist.df.sub$Type <- gsub('SNP','SNV',nist.df.sub$Type)
nist.df.sub$Depth <- str_split(nist.df.sub$source,'\\.',simplify = T)[,2]
nist.df.sub$F1.score <- nist.df.sub$METRIC.F1_Score
nist.df.sub$Quality <- lapply(nist.df.sub$source,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
nist.df.sub <- subset(nist.df.sub,Quality!='Q35')
nist.df.sub$Sample <- 'NA12878'
nist.df.sub$Depth <- factor(nist.df.sub$Depth,levels = c('10X','15X','20X','25X','30X','60X','90X','120X'))
#突变数量----
ct_nist <- read.csv('ELE_ILM/Haplotyper_WES.vcf_stats.NIST.csv')
ct_nist$source <- gsub('.Haplotyper','',ct_nist$source)

ct_nist$Quality <- lapply(ct_nist$source,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
ct_nist <- subset(ct_nist,Quality!='Q35')
ct_nist$Depth <- str_split(ct_nist$source,'\\.',simplify = T)[,2]
ct_nist.sum <- ct_nist %>%
  dplyr::group_by(Depth, Quality, type) %>%
  dplyr::summarize(
    across(
      count,
      list(
        median = ~median(.,na.rm = T),
        sd = ~sd(.,na.rm = T)
      )
    )
  )
ct_nist.sum$Depth <- factor(ct_nist.sum$Depth,levels = c('10X','15X','20X','25X','30X','60X','90X','120X'))
ct_nist.sum$type <- gsub('Indels','INDEL',ct_nist.sum$type) %>% gsub('SNVs','SNV',.)
nist.df.sub$type = nist.df.sub$Type
nist_snv <- get_boxcompare(tnseq.df.sub=nist.df.sub,m_type='SNV',y_lim=c(0.5,1.05),ct_1395.sum=ct_nist.sum,d_type='F1.score');nist_snv
#hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.3,0.9),ct_1395.sum);hcc_indel
nist_indel <- get_boxcompare(nist.df.sub,'INDEL',c(0.5,1.02),ct_nist.sum,'F1.score');nist_indel
tnseq.df.sub=nist.df.sub
m_type='SNV'
y_lim=c(0.5,1.05)
ct_1395.sum=ct_nist.sum
d_type='F1.score'

## Quartet----
#ELE & ILM
q.f1 <- read.csv("ELE_ILM/F1score_stats_total_Quartet.csv",sep = ',')
q.f1$source <- gsub('_10G_','.',q.f1$source)
q.f1.sub <- q.f1[,c('Type','source','METRIC.F1_Score')]
q.f1.sub$Type <- gsub('SNP','SNV',q.f1.sub$Type)
q.f1.sub$Depth <- str_split(q.f1.sub$source,'\\.',simplify = T)[,2]
q.f1.sub$F1.score <- q.f1.sub$METRIC.F1_Score
q.f1.sub$Platform <- str_split(q.f1.sub$source,'\\_',simplify = T)[,3]
q.f1.sub$Quality <- lapply(q.f1.sub$source,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()

q.f1.sub <- subset(q.f1.sub,Quality!='Q35')
q.f1.sub <- q.f1.sub[!apply(q.f1.sub, 1, function(x) any(grepl("WES_D6_2_ELE_10G.120X", x))), ]
q.f1.sub$Sample='Quartet'
q.f1.sub$Depth <- factor(q.f1.sub$Depth,levels = c('10X','15X','20X','25X','30X','60X','90X','120X'))

subset(q.f1.sub,Depth=='20X') %>%
  group_by(Quality,Type) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(F1.score, na.rm = TRUE),   # 计算均值
    SD_Duplicates = sd(F1.score, na.rm = TRUE)       # 计算标准差
  )
#突变数量----
ct_quartet <- data.frame()
quartet_lst = c('D5.vcf_stats.csv','D6.vcf_stats.csv','F7.vcf_stats.csv','M8.vcf_stats.csv')
for (i in quartet_lst){
  df <- read.csv(paste0('ELE_ILM/',i))
  ct_quartet <- rbind(ct_quartet,df)
}
ct_quartet$source <- gsub('.Haplotyper','',ct_quartet$source)

ct_quartet$Quality <- lapply(ct_quartet$source,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
ct_quartet <- subset(ct_quartet,Quality!='Q35')
ct_quartet$Depth <- str_split(ct_quartet$source,'\\.',simplify = T)[,2]
ct_quartet.sum <- ct_quartet %>%
  dplyr::group_by(Depth, Quality, type) %>%
  dplyr::summarize(
    across(
      count,
      list(
        median = ~median(.,na.rm = T),
        sd = ~sd(.,na.rm = T)
      )
    )
  )
ct_quartet.sum$Depth <- factor(ct_quartet.sum$Depth,levels = c('10X','15X','20X','25X','30X','60X','90X','120X'))
ct_quartet.sum$type <- gsub('Indels','INDEL',ct_quartet.sum$type) %>% gsub('SNVs','SNV',.)

q.f1.sub$type <- q.f1.sub$Type
quartet_snv <- get_boxcompare(q.f1.sub,'SNV',c(0.6,1.03),ct_nist.sum);quartet_snv
#hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.3,0.9),ct_1395.sum);hcc_indel
quartet_indel <- get_boxcompare(q.f1.sub,'INDEL',c(0.5,1.0),ct_nist.sum);quartet_indel


#全基因组范围内--
#ELE & ILM
q.mcr <- read.csv("ELE_ILM/MCR_stats_total_Quartet.csv")
q.mcr$source <- gsub('_10G_','.',q.mcr$source)
q.mcr$Depth <- str_split(q.mcr$source,'\\.',simplify = T)[,2]
q.mcr$Platform <- str_split(q.mcr$source,'\\_',simplify = T)[,3]
q.mcr$Quality <- lapply(q.mcr$source,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
q.mcr$Family <- gsub('_10G_','.',q.mcr$Family)
q.mcr$Type <- str_split(q.mcr$Family,'\\.',simplify = T)[,3]

q.mcr <- subset(q.mcr,source != 'WES_2_ELE.120X')
q.mcr <- subset(q.mcr,Quality !='Q35')
q.mcr$Sample <- 'Quartet'
q.mcr$MCR <- q.mcr$Mendelian_Concordance_Rate

q.mcr$type <- q.mcr$Type
q.mcr$Depth <- factor(q.mcr$Depth,levels = c('10X','15X','20X','25X','30X','60X','90X','120X'))
quartet_snv.mcr <- get_boxcompare(q.mcr,'SNV',c(0.4,1.03),d_type = 'MCR');quartet_snv.mcr
#hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.3,0.9),ct_1395.sum);hcc_indel
quartet_indel.mcr <- get_boxcompare(q.mcr,'INDEL',c(0.4,0.9),d_type = 'MCR');quartet_indel.mcr

subset(q.mcr,Depth=='30X') %>%
  group_by(Quality,Type) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(MCR, na.rm = TRUE),   # 计算均值
    SD_Duplicates = sd(MCR, na.rm = TRUE)       # 计算标准差
  )


get_boxcompare_rmsk = function(tnseq.df.sub,m_type, y_lim,ct_1395.sum=data.frame(),d_type='F1.score',xmax=5.5){
  tnseq.df.sub <- subset(tnseq.df.sub,type == m_type)
  tnseq.df.sub['F1.score'] <- tnseq.df.sub[d_type]
  
  #tnseq.df.test$y.position <- ifelse(tnseq.df.test$group2=='Q40',tnseq.df.test$y.position-0.01,tnseq.df.test$y.position-0.03)
  tnseq.df.sub$group <- paste0(tnseq.df.sub$Quality,'_',tnseq.df.sub$Depth)
  if (length(ct_1395.sum) > 0){
    ct_1395_plot <- ggplot(subset(ct_1395.sum,type==m_type), aes(x = Depth, y = count_median,fill=Quality)) +
      geom_rect(xmin = 0.5, xmax = 1.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9") +
      geom_rect(xmin = 2.5, xmax = 3.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9")+
      geom_rect(xmin = 4.5, xmax = 5.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9") +
      geom_rect(xmin = 6.5, xmax = 7.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9")+
      geom_bar(stat = "identity", position = "dodge", width = 0.7)+
      geom_errorbar(aes(ymin = count_median - count_sd, ymax = count_median + count_sd), position = position_dodge(0.7), width = 0.2)+
      scale_fill_manual(values = colors)+ylab('')+
      my_theme_2+theme(axis.text.x =element_text(color = "black",size = 10,face = "bold",angle = 15,hjust = 1))+theme(legend.position = "none")+ylab('Count')+xlab('')
    if (unique(tnseq.df.sub$Sample)=='HCC1395/BL' & m_type=='INDEL'){
        box_p <- ggplot(tnseq.df.sub, aes(x = Depth, y = F1.score)) +
          geom_rect(xmin = 0.5, xmax = 1.5,
                    ymin = 0, ymax = 11, fill ="#E9E9E9") +
          geom_rect(xmin = 2.5, xmax = 3.5,
                    ymin = 0, ymax = 11, fill ="#E9E9E9")+
          geom_rect(xmin = 4.5, xmax = 5.5,
                    ymin = 0, ymax = 2000000, fill ="#E9E9E9") +
          geom_rect(xmin = 6.5, xmax = 7.5,
                    ymin = 0, ymax = 2000000, fill ="#E9E9E9")+
          geom_boxplot(aes(fill=Quality,color=Quality),
                       outlier.shape = NA,  # 隐藏离群点
                       width = 0.6,alpha=0.3,size = 0.5
          )+
          geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                      size = 3, alpha=0.5,aes(color = Quality))+
          #geom_line(aes(color = Quality,group=Quality),position = position_dodge(width = 0.5),alpha=0.5,size=1)+
          scale_color_manual(values = colors)+
          scale_fill_manual(values = colors)+
          #facet_grid(type~.)+
          ylim(y_lim)+ylab(paste0(d_type,' ',m_type))+ggtitle(unique(tnseq.df.sub$Sample))+
          my_theme+theme(axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 15,hjust = 1))+xlab('')
        box_p = box_p +annotation_custom(grob=ggplotGrob(ct_1395_plot),ymin = y_lim[1], ymax=y_lim[1]+0.3, xmin=1, xmax=xmax)
      } else {
        tnseq.df.test <- tnseq.df.sub %>%
          group_by(type, Depth) %>%
          filter(n_distinct(F1.score) > 1, n_distinct(Quality) > 1) %>% 
          t_test(F1.score ~ Quality, p.adjust.method = 'none', comparisons = c("Q30", "Q40")) %>%
          mutate(p.adj.signif = p.adjust(p, method = "none")) %>%  # Adjust p-value and calculate significance
          add_significance() %>%  # Adds significance symbols based on adjusted p-values
          select(-df, -statistic, -p)
        print(tnseq.df.test)
        tnseq.df.test <- tnseq.df.test %>% add_xy_position(x = "Depth")
        Max_F1_score <- tnseq.df.sub %>%
          group_by(type, Depth) %>%
          dplyr::summarise(y.position = max(F1.score, na.rm = TRUE))
        tnseq.df.test<- merge(tnseq.df.test,Max_F1_score,by=c('type','Depth'))
        tnseq.df.test$y.position <- tnseq.df.test$y.position.y
        box_p <- ggplot(tnseq.df.sub, aes(x = Depth, y = F1.score)) +
          geom_rect(xmin = 0.5, xmax = 1.5,
                    ymin = 0, ymax = 11, fill ="#E9E9E9") +
          geom_rect(xmin = 2.5, xmax = 3.5,
                    ymin = 0, ymax = 11, fill ="#E9E9E9")+
          geom_rect(xmin = 4.5, xmax = 5.5,
                    ymin = 0, ymax = 2000000, fill ="#E9E9E9") +
          geom_rect(xmin = 6.5, xmax = 7.5,
                    ymin = 0, ymax = 2000000, fill ="#E9E9E9")+
          geom_boxplot(aes(fill=Quality,color=Quality),
                       outlier.shape = NA,  # 隐藏离群点
                       width = 0.6,alpha=0.3,size = 0.5
          )+
          geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                      size = 3, alpha=0.5,aes(color = Quality))+
          #geom_line(aes(color = Quality,group=Quality),position = position_dodge(width = 0.5),alpha=0.5,size=1)+
          scale_color_manual(values = colors)+
          scale_fill_manual(values = colors)+
          #facet_grid(type~.)+
          ylim(y_lim)+ylab(paste0(d_type,' ',m_type))+ggtitle(unique(tnseq.df.sub$Sample))+
          my_theme+theme(axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 15,hjust = 1))+
          stat_pvalue_manual(
            tnseq.df.test, label = "p.signif", 
            step.increase = 0.01,size = 5
          )+xlab('')
        box_p = box_p +annotation_custom(grob=ggplotGrob(ct_1395_plot),ymin = y_lim[1], ymax=y_lim[1]+0.4, xmin=2, xmax=xmax)
        
      }
  }
  return(box_p)
}
#复杂基因组区域-------
#NIST---
#ELE & ILM
df_nist = read.csv("ELE_ILM/F1score_stats_total_NIST_rmsk.csv")
df_nist.sub <- df_nist[,c('Type','source','METRIC.F1_Score','QUERY.TOTAL')]
df_nist.sub$Type <- gsub('SNP','SNV',df_nist.sub$Type)
df_nist.sub$Rmsk_class <- str_split(df_nist.sub$source,'\\.',simplify = T)[,2]
df_nist.sub$F1.score <- df_nist.sub$METRIC.F1_Score
#df_nist.sub$Platform <- str_split(df_nist.sub$source,'\\_',simplify = T)[,4]
df_nist.sub$Quality <- lapply(str_split(df_nist.sub$source,'\\_',simplify = T)[,4],function(x){
  p = ifelse(grepl('ILM',x),'Q30','Q40')
  return(p)
}) %>% unlist()

df_nist.sub$Sample <- 'NA12878'
df_nist.sub$count <- df_nist.sub$QUERY.TOTAL

df_nist.sub$type <- df_nist.sub$Type
#ct_quartet$Depth <- str_split(ct_quartet$source,'\\.',simplify = T)[,2]
df_nist.sub.rmsk <- df_nist.sub %>%
  dplyr::group_by(Rmsk_class, Quality, type) %>%
  dplyr::summarize(
    across(
      count,
      list(
        median = ~median(.,na.rm = T),
        sd = ~sd(.,na.rm = T)
      )
    )
  )
df_nist.sub.rmsk$Depth <- df_nist.sub.rmsk$Rmsk_class
df_nist.sub$Depth <- df_nist.sub$Rmsk_class
tnseq.df.sub = df_nist.sub
m_type = 'SNV'
y_lim=c(0.4,1.03)
ct_1395.sum=df_nist.sub.rmsk
d_type='F1.score'

nist_snv.rmsk <- get_boxcompare_rmsk(df_nist.sub,'SNV',c(0.5,1.08),df_nist.sub.rmsk,d_type = 'F1.score',xmax = 6.5);nist_snv.rmsk
#hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.3,0.9),ct_1395.sum);hcc_indel
nist_indel.rmsk <- get_boxcompare_rmsk(df_nist.sub,'INDEL',c(0.3,1.1),df_nist.sub.rmsk,d_type = 'F1.score');nist_indel.rmsk

#df_nist.sub$Lib <- 'Lib1'

#Quartet F1---
#ELE & ILM
df_quartet = read.csv("ELE_ILM/Quartet.rmsk.F1.csv")
df_quartet.sub <- df_quartet[,c('Type','source','METRIC.F1_Score','QUERY.TOTAL')]
df_quartet.sub$Type <- gsub('SNP','SNV',df_quartet.sub$Type)
df_quartet.sub$Rmsk_class <- str_split(df_quartet.sub$source,'\\.',simplify = T)[,2]
df_quartet.sub$F1.score <- df_quartet.sub$METRIC.F1_Score
df_quartet.sub$Quality <- lapply(str_split(df_quartet.sub$source,'\\_',simplify = T)[,4],function(x){
  p = ifelse(grepl('ILM',x),'Q30','Q40')
  return(p)
}) %>% unlist()

df_quartet.sub$Sample <- 'Quartet'
df_quartet.sub$count <- df_quartet.sub$QUERY.TOTAL
df_quartet.sub$type <- df_quartet.sub$Type
df_quartet.sub$Depth <- df_quartet.sub$Rmsk_class
#ct_quartet$Depth <- str_split(ct_quartet$source,'\\.',simplify = T)[,2]
ct_quartet.rmsk <- df_quartet.sub %>%
  dplyr::group_by(Rmsk_class, Quality, type) %>%
  dplyr::summarize(
    across(
      count,
      list(
        median = ~median(.,na.rm = T),
        sd = ~sd(.,na.rm = T)
      )
    )
  )
ct_quartet.rmsk$Depth <- ct_quartet.rmsk$Rmsk_class
quartet_snv.rmsk <- get_boxcompare_rmsk(df_quartet.sub,'SNV',c(0.5,1.08),ct_quartet.rmsk,d_type = 'F1.score',xmax = 6.5);quartet_snv.rmsk
#hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.3,0.9),ct_1395.sum);hcc_indel
quartet_indel.rmsk <- get_boxcompare_rmsk(df_quartet.sub,'INDEL',c(0.4,1.1),ct_quartet.rmsk,d_type = 'F1.score');quartet_indel.rmsk


#Quartet MCR---
#ELE & ILM
df_quartet.mcr = read.csv("ELE_ILM/Quartet.rmsk.MCR.csv")
df_quartet.mcr$type <- str_split(df_quartet.mcr$Family,'\\.',simplify = T)[,3]
df_quartet.mcr$Rmsk_class <- str_split(df_quartet.mcr$Family,'\\.',simplify = T)[,2]
df_quartet.mcr$Sample <- 'Quartet'
df_quartet.mcr$Quality <- lapply(str_split(df_quartet.mcr$Family,'\\_',simplify = T)[,3],function(x){
  p = ifelse(grepl('ILM',x),'Q30','Q40')
  return(p)
}) %>% unlist()
df_quartet.mcr$MCR = df_quartet.mcr$Mendelian_Concordance_Rate

df_quartet.mcr <- subset(df_quartet.mcr,!grepl('WES_2_ELE_10G',df_quartet.mcr$Family)) 


df_quartet.mcr$Depth <- df_quartet.mcr$Rmsk_class
quartet_snv.rmsk.mcr <- get_boxcompare_rmsk(df_quartet.mcr,'SNV',c(0.5,1.05),ct_quartet.rmsk,d_type = 'MCR');quartet_snv.rmsk.mcr
#hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.3,0.9),ct_1395.sum);hcc_indel
quartet_indel.rmsk.mcr <- get_boxcompare_rmsk(df_quartet.mcr,'INDEL',c(0.2,1),ct_quartet.rmsk,d_type = 'MCR');quartet_indel.rmsk.mcr


subset(df_quartet.mcr,Sample=='Quartet') %>%
  group_by(Quality,Depth,type) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(MCR, na.rm = TRUE),   # 计算均值
    SD_Duplicates = sd(MCR, na.rm = TRUE)       # 计算标准差
  )



#HCC1395.高置信区间------
#ELE_ILM
tnseq.df <- read.csv("ELE_ILM/HCC1395_F1score_stats_total.SEQC2.csv")
#tnseq.df <- read.csv("ELE_ILM/HCC1395_F1score_stats_total.csv")
#tnseq.df <- read.csv("ELE_ILM/HCC1395_F1score_stats_total.strelka2.csv")
tnseq.df.sub <- tnseq.df[,c('type','source','F1.score','recall','precision')]
tnseq.df.sub$type <- gsub('SNVs','SNV',tnseq.df.sub$type) %>% gsub('indels','INDEL',.)
tnseq.df.sub$Depth <- str_split(tnseq.df.sub$source,'\\.',simplify = T)[,2]
#tnseq.df.sub$Platform <- str_split(tnseq.df.sub$source,'\\_',simplify = T)[,4]
#突变数量----
ct_1395 <- read.csv('ELE_ILM/Mutect2.vcf_stats.csv')
#ct_1395 <- read.csv('ELE_ILM/strelka_vcf.vcf_stats.csv')
ct_1395$source <- gsub('.filter','',ct_1395$source)
#ct_1395_raw <- read.csv('ELE_ILM/Mutect2.vcf_stats.raw.csv')
#ct_1395 <- merge(ct_1395,ct_1395_raw,by=c('source','type'),suffixes = c(".filter",".raw"))
ct_1395$Quality <- lapply(ct_1395$source,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
ct_1395 <- subset(ct_1395,Quality != 'Q35')
ct_1395$Depth <- str_split(ct_1395$source,'\\.',simplify = T)[,2]
ct_1395.sum <- ct_1395 %>%
  dplyr::group_by(Depth, Quality, type) %>%
  dplyr::summarize(
    across(
      count,
      list(
        median = ~median(.,na.rm = T),
        sd = ~sd(.,na.rm = T)
      )
    )
  )
ct_1395.sum$Depth <- factor(ct_1395.sum$Depth,levels = c('30X','60X','90X','120X'))
ct_1395.sum$type <- gsub('Indels','INDEL',ct_1395.sum$type) %>% gsub('SNVs','SNV',.)
ct_1395.sum <- subset(ct_1395.sum,Quality != 'Q35')
tnseq.df.sub$Quality <- lapply(tnseq.df.sub$source,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
tnseq.df.sub <- subset(tnseq.df.sub,Quality !='Q35')
tnseq.df.sub$Type <- tnseq.df.sub$type
tnseq.df.sub$Sample <- 'HCC1395/BL'
tnseq.df.sub$Depth <- factor(tnseq.df.sub$Depth,levels = c('30X','60X','90X','120X'))

subset(tnseq.df.sub,Depth=='60X') %>% 
  group_by(Quality,type) %>% 
  dplyr::summarise(F1.score_mean = mean(F1.score,na.rm = T),
                   F1.score_sd = sd(F1.score,na.rm = T))

subset(tnseq.df.sub,Depth=='90X') %>% 
  group_by(Quality,type) %>% 
  dplyr::summarise(F1.score_mean = mean(F1.score,na.rm = T),
                   F1.score_sd = sd(F1.score,na.rm = T))
  

hcc_snv <- get_boxcompare(tnseq.df.sub,'SNV',c(0.45,0.9),ct_1395.sum);hcc_snv
#hcc_snv <- get_boxcompare(tnseq.df.sub,'SNV',c(0.7,1.0),ct_1395.sum);hcc_snv
hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.3,0.8),ct_1395.sum);hcc_indel
#hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.4,0.85),ct_1395.sum);hcc_indel

hcc_f1.p <- (hcc_snv/hcc_indel)+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))
hcc_f1.p.strelka <- (hcc_snv/hcc_indel)+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))

ggsave('Fig.hcc_f1.pdf',hcc_f1.p,width=8.15, height=5.7*1.5,dpi = 300)
ggsave('Fig.hcc_f1.strelka.pdf',hcc_f1.p.strelka,width=8.15, height=5.7*1.5,dpi = 300)


##HCC1395 复杂基因组区域-----
#HCC1395---
#ELE & ILM
df_hcc1395 = read.csv("ELE_ILM/F1score_stats_total_HCC1395_rmsk.SEQC2.csv")
df_hcc1395$QUERY.TOTAL <- df_hcc1395$total.query
df_hcc1395.sub <- df_hcc1395[,c('type','source','F1.score','QUERY.TOTAL')]
df_hcc1395.sub$Type <- df_hcc1395.sub$type <- gsub('SNVs','SNV',df_hcc1395.sub$type) %>% gsub('indels','INDEL',.)
df_hcc1395.sub$Rmsk_class <- str_split(df_hcc1395.sub$source,'\\.',simplify = T)[,2]
#df_hcc1395.sub$Platform <- str_split(df_hcc1395.sub$source,'\\_',simplify = T)[,4]
df_hcc1395.sub <- subset(df_hcc1395.sub,!grepl('Mix',df_hcc1395.sub$source))
df_hcc1395.sub$Quality <- lapply(str_split(df_hcc1395.sub$source,'\\_',simplify = T)[,4],function(x){
  p = ifelse(grepl('ILM',x),'Q30','Q40')
  return(p)
}) %>% unlist()

df_hcc1395.sub$Sample <- 'HCC1395/BL'
df_hcc1395.sub$count <- df_hcc1395.sub$QUERY.TOTAL
df_hcc1395.sub$type <- df_hcc1395.sub$Type
df_hcc1395.sub <- subset(df_hcc1395.sub,Type != 'records')
df_hcc1395.sub$Depth <- df_hcc1395.sub$Rmsk_class
rmsk_df_median.hcc <- df_hcc1395.sub %>%
  dplyr::group_by(type, Quality,Rmsk_class,Sample) %>%
  dplyr::summarize(
    across(
      count,
      list(
        median = ~median(.,na.rm=T),
        sd = ~sd(.,na.rm=T)
      )
    )
  )
rmsk_df_median.hcc$Depth <- rmsk_df_median.hcc$Rmsk_class
tnseq.df.sub = df_hcc1395.sub
m_type = 'INDEL'
y_lim=c(0.1,1.03)
ct_1395.sum=rmsk_df_median.hcc
d_type='F1.score'

hcc_snv.rmsk <- get_boxcompare_rmsk(df_hcc1395.sub,'SNV',c(0.2,1.1),rmsk_df_median.hcc,xmax = 4.5);hcc_snv.rmsk
#hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.3,0.9),ct_1395.sum);hcc_indel
hcc_indel.rmsk <- get_boxcompare_rmsk(df_hcc1395.sub,'INDEL',c(0.6,1),rmsk_df_median.hcc,xmax = 0.6);hcc_indel.rmsk


#WES靶向捕获区域中的重复区域
rep_df <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/文章/Q40_repeat-information.csv')
colnames(rep_df)[1] <- 'Repeat type'
rep_df.long <- melt(rep_df)
colnames(rep_df.long) <- c('Repeat type','Lib','Region')
rep_df.p1 <- ggplot(data=subset(rep_df.long,Lib=='Lib1'),aes(x=`Repeat type`,y=Region)) +
  geom_bar(stat="summary",fun=mean,position=position_dodge(0.9),alpha=1,color='gray',fill='#91569F')+
  my_theme+
  theme(axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 30,hjust = 1))+
  labs(x=" ",y='Region(bp)')+
  #expand_limits(y = y_lim)+
  guides(fill=guide_legend(title = "",nrow = 1,byrow = FALSE));rep_df.p1

rep_df.p1



#HCC1395 低频突变区域------
#低频突变的F1 score-----
#ELE & ILM
lowvaf.df <- read.csv("ELE_ILM/F1score_stats_total.LowVAF.csv")
lowvaf.df.sub <- lowvaf.df[,c('type','source','F1.score')]
lowvaf.df.sub$type <- gsub('SNVs','SNV',lowvaf.df.sub$type) %>% gsub('indels','INDEL',.)
lowvaf.df.sub$Depth <- str_split(lowvaf.df.sub$source,'\\.',simplify = T)[,2]
#lowvaf.df.sub$Platform <- str_split(lowvaf.df.sub$source,'\\_',simplify = T)[,4]
lowvaf.df.sub <- subset(lowvaf.df.sub,!grepl('Mix',lowvaf.df.sub$source))
lowvaf.df.sub$Quality <- lapply(str_split(lowvaf.df.sub$source,'\\_',simplify = T)[,4],function(x){
  p = ifelse(grepl('ILM',x),'Q30','Q40')
  return(p)
}) %>% unlist()
lowvaf.df.sub$Type <- lowvaf.df.sub$type

hcc_snv.lowVAF <- get_boxcompare(tnseq.df.sub=lowvaf.df.sub,m_type='SNV',y_lim=c(0.2,0.7),d_type = 'F1.score');hcc_snv.lowVAF
#hcc_indel <- get_boxcompare(tnseq.df.sub,'INDEL',c(0.3,0.9),ct_1395.sum);hcc_indel
hcc_indel.lowVAF <- get_boxcompare(lowvaf.df.sub,'INDEL',c(0.2,0.7),d_type = 'F1.score');hcc_indel.lowVAF

lowvaf.df.sub %>% 
  group_by(Quality,Type) %>% 
  dplyr::summarise(F1.score_mean = mean(F1.score,na.rm = T),
                   F1.score_sd = sd(F1.score,na.rm = T))

#补充图中需要放的图----
#HCC1395 VAF与标准集的相关性------
ele.vaf <- read.csv('ELE_ILM/HCC1395.VAF.ELE.csv')
breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
# Create labels for the intervals
labels <- c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]", "(0.6,0.8]", "(0.8,1]")
ele.vaf$Interval <- cut(ele.vaf$VAF.Mean, breaks = breaks, labels = labels, include.lowest = TRUE)
ele.vaf$Quality <- 'Q40'
ilm.vaf <- read.csv('ELE_ILM/HCC1395.VAF.ILM.csv')
ilm.vaf$Interval <- cut(ilm.vaf$VAF.Mean, breaks = breaks, labels = labels, include.lowest = TRUE)
ilm.vaf$Quality <- 'Q30'

vaf_df <- rbind(ele.vaf[c('VAF.STD','Interval','Quality','tag','VAF.Mean','CV')],ilm.vaf[c('VAF.STD','Interval','Quality','tag','VAF.Mean','CV')])
vaf_df$Sample='HCC1395/BL'

ref_vcf <- read.table('ELE_ILM/high-confidence_sSNV_sIndel_v1.HCR.spikein.vcf')

ref_vcf$tag <- paste0(ref_vcf$V1,'_',ref_vcf$V2,'_',ref_vcf$V4,'_',ref_vcf$V5)
ref_vcf$VAF.ref <- str_split(ref_vcf$V8,'\\;',simplify = T)[,5] %>% gsub('VAF=','',.) %>% as.numeric()
breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
# Create labels for the intervals
labels <- c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]", "(0.6,0.8]", "(0.8,1]")
ref_vcf$Interval <- cut(ref_vcf$VAF.ref, breaks = breaks, labels = labels, include.lowest = TRUE)

ref_vcf.all <- merge(vaf_df[c('VAF.Mean','Quality','tag')],ref_vcf[c('tag','VAF.ref','Interval')],by="tag")
#ref_vcf.ele <- merge(ele.vaf[c('VAF.Mean','Quality','tag')],ref_vcf[c('tag','VAF.ref','Interval')],by="tag")

#ref_vcf.ilm <- merge(ilm.vaf[c('VAF.Mean','Interval','Quality','tag')],ref_vcf[c('tag','VAF.ref')],by="tag")

library(plyr)
cor.df <- ddply(ref_vcf.all, .(Interval,Quality), summarise, correlation = cor.test(VAF.Mean, VAF.ref)$estimate)
#cor.ele <- ddply(ref_vcf.ele, .(Interval), summarise, correlation = cor.test(VAF.Mean, VAF.ref)$estimate)
#cor.ilm <- ddply(ref_vcf.ilm,.(Interval), summarise, correlation = cor.test(VAF.Mean, VAF.ref)$estimate)

#cor.ele$Quality <- 'Q40'
#cor.ilm$Quality <- 'Q30'

cor.df$Sample <- 'HCC1395/BL'
cor.df$Rmsk_class <- cor.df$Interval
cor.df$Sample='HCC1395/BL'

y_title='VAF.PCC'
y_lim = c(0,1.05)
cor.p.somatic <- ggplot(data=cor.df,aes(x=Interval,y=correlation,fill=Quality)) +
  geom_bar(stat="summary",fun=mean,position=position_dodge(0.9),alpha=0.8)+
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),color='black',width=.2,
  #              position=position_dodge(0.9))+
  stat_summary(geom = "text", aes(label = paste(round(..y.., digits = 2))),
               position = position_dodge(width = 1), vjust = -0.9,hjust=0.5,size = 3,color='black')+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  my_theme+
  theme(axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 30,hjust = 1))+
  labs(x="VAF Interval",y=y_title)+
  expand_limits(y = y_lim)+
  guides(fill=guide_legend(title = "",nrow = 1,byrow = FALSE));cor.p.somatic


#不同深度下的变异位点可重复性Reproducibility--------------
#NIST
#ELE & ILM
df_nist = read.csv('ELE_ILM/NIST_Reproducibility_stats_total.csv') 
df_nist$Reproducibility <- (df_nist$TRUTH.TP/df_nist$TRUTH.TOTAL + df_nist$TRUTH.TP/df_nist$QUERY.TOTAL)/2
df_nist$Source <- gsub('_10G','',df_nist$source)
df_nist$Depth <- str_split(df_nist$Source,'\\.',simplify = T)[,2]
df_nist$Sample <- 'NA12878'
add_batch.e <- function(df){
  df[(str_count(df$source,'ELE')==2),'Quality'] = 'Q40'
  df[(str_count(df$source,'ILM')==2),'Quality'] = 'Q30'
  df[(str_count(df$source,'ILM')==1),'Quality'] = 'Q40_vs_Q30'
  return(df)
}
df_nist <- add_batch.e(df_nist)
df_nist$type <- gsub('SNP','SNV',df_nist$Type) 

df_nist$Depth <- factor(df_nist$Depth,levels = c('10X','15X','20X','25X','30X','60X','90X','120X'))

#Quartet
#ELE & ILM
df_quartet = read.csv('ELE_ILM/Quartet_Reproducibility_stats_total.csv')
df_quartet$Reproducibility <- (df_quartet$TRUTH.TP/df_quartet$TRUTH.TOTAL + df_quartet$TRUTH.TP/df_quartet$QUERY.TOTAL)/2
df_quartet$type <- gsub('SNP','SNV',df_quartet$Type) 
df_quartet$Source <- gsub('_10G','',df_quartet$source)

df_quartet$Depth <- str_split(df_quartet$Source,'\\.',simplify = T)[,2]
df_quartet$Sample <- 'Quartet'
df_quartet <- add_batch.e(df_quartet)
df_quartet <- df_quartet[!apply(df_quartet, 1, function(x) any(grepl("D6_2_ELE_10G.120X.Haplotyper", x))), ]

#HCC1395
#ELE & ILM
hcc1395_dirs <- dir(path = 'ELE_ILM',pattern = 'HCC1395_F1score_stats_total.')
df_lst = list()
for (d in hcc1395_dirs){
  print(d)
  dep = str_split(d,'\\.',simplify = T)[,2]
  df = read.csv(paste0('ELE_ILM/',d))
  df$type <- gsub('indels','INDEL',df$type) %>% gsub('SNVs','SNV',.)
  add_batch <- function(df){
    df[(str_count(df$source,'ELE')==2),'Quality'] = 'Q40'
    df[(str_count(df$source,'ILM')==2),'Quality'] = 'Q30'
    df[(str_count(df$source,'ILM')==1),'Quality'] = 'Q40_vs_Q30'
    return(df)
  }
  
  df <- add_batch(df)
  df$Depth <- dep
  head(df) %>% print()
  df_lst[[d]] = df
}
df_hcc1395 = do.call(rbind,df_lst)
df_hcc1395$Depth <- factor(df_hcc1395$Depth,levels = c('30X','60X','90X','120X'))

df_hcc1395$Sample <- 'HCC1395/BL'
#rownames(df_nist) <- NULL
#rownames(df_hcc1395) <- NULL
rep_all <- rbind(df_nist[,c('type','Reproducibility','Quality','Sample','Depth')],
                 df_quartet[,c('type','Reproducibility','Quality','Sample','Depth')]) %>% 
  rbind(.,df_hcc1395[,c('type','Reproducibility','Quality','Sample','Depth')])
rep_all <- subset(rep_all,Quality!='Q40_vs_Q30' & Quality!='Q40_vs_Q30_M' & type %in% c('INDEL','SNV'))
rep_all$Depth <- gsub('X','',rep_all$Depth) %>% as.numeric()
rep_all$Type <- rep_all$type
rep_all_median <- rep_all %>%
  dplyr::group_by(Type, Quality, Depth,Sample) %>%
  dplyr::summarize(
    across(
      Reproducibility,
      list(
        mean = ~mean(.),
        sd = ~sd(.)
      )
    )
  )

rep_all_median <- subset(rep_all_median,Type != 'records')
rep_all_median$Depth <- factor(rep_all_median$Depth,levels=c(0,5,10,15,20,25,30,60,90,120))
#rep_all_median$group <- paste0(rep_all_median$Quality,'-',rep_all_median$Sample,'-',rep_all_median$Type,'-',rep_all_median$Lib)

rep_all_median$Sample <- gsub('NIST','NA12878',rep_all_median$Sample)
rep.p.somatic <- get_dotbar(subset(rep_all_median,Sample=='HCC1395/BL'),'Reproducibility','Reproducibility',y_lim=c(0.5,1))
rep.p.somatic
rep.p.germline <- get_dotbar(subset(rep_all_median,Sample!='HCC1395/BL'),'Reproducibility','Reproducibility',y_lim=c(0.58,1))
rep.p.germline

#

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

germline_legend <- g_legend(f1.p.germline)  %>% plot()


#nist_snv,nist_indel,quartet_snv,quartet_indel,quartet_snv.mcr,quartet_indel.mcr,nist_snv.rmsk,nist_indel.rmsk,quartet_snv.rmsk,quartet_indel.rmsk

germline_main <- ((quartet_snv.mcr|quartet_indel.mcr)/((quartet_snv|quartet_indel))/((quartet_snv.rmsk.mcr|quartet_indel.rmsk.mcr)))+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect',heights = c(1,1,1)) & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))
germline_main
ggsave('Fig.germline-1.pdf',germline_main,width=8.15*1.5, height=5.7*2.5,dpi = 300)

#hcc_snv,hcc_indel,hcc_snv.rmsk,hcc_indel.rmsk,hcc_snv.lowVAF,hcc_indel.lowVAF
somatic_main <- (hcc_snv|hcc_indel)/((hcc_snv.rmsk+hcc_indel.rmsk+hcc_snv.lowVAF+
                                        hcc_indel.lowVAF)+plot_layout(guides='collect',widths = c(3,1,1,1)))/
  (CNV_Q40_JI.all|CNV_Q40_PCC.all)+
                 plot_annotation(tag_levels = 'a') & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))
somatic_main
ggsave('Fig.somatic_1.pdf',somatic_main,width=8.15*1.5, height=5.7*2.5,dpi = 300)


sup_wes <-(nist_snv|nist_indel)/(nist_snv.rmsk|nist_indel.rmsk)/(rep.p.germline)/(rep.p.somatic|cor.p.somatic|rep_df.p1)+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))

ggsave('Fig.wes_sup.pdf',sup_wes,width=8.15*1.5, height=5.7*3.5,dpi = 300)



#
nist_mut <- list()
for (i in list.files(path = 'ELE_ILM/',pattern = '*10X.Haplotyper.bed')){
  type='ELE'
  if (grepl('Mix',i)){
    type='Mix'
  } else if (grepl('ILM',i)){
    type='ILM'
  }
  df <- read.csv(paste0('ELE_ILM/',i),sep = '\t',header = F)
  #df[type] = paste(df$V1,df$V2,df$V4,df$V5,sep = ':')
  nist_mut[[type]] = paste(df$V1,df$V2,df$V4,df$V5,sep = ':')
}

combine_df <- data.frame(ID = c(nist_mut[['Mix']],nist_mut[['ELE']],nist_mut[['ILM']]),
  Lab = c(rep('Mix',length(nist_mut[['Mix']])),rep('ELE',length(nist_mut[['ELE']])),rep('ILM',length(nist_mut[['ILM']]))))
intersect(nist_mut[['Mix']],nist_mut[['ELE']]) %>% length()

library(ggVennDiagram)
library(ggplot2)
library("animation")


combine_df <- na.omit(combine_df)
venn_plot <- function(df,p_name){
  list_df = list()
  for (i in df[,2]){
    list_df[[i]] = df[df[,2]==i,1]
  }
  
  vp = ggVennDiagram(list_df,
                     label = "count",
                     edge_lty = "dashed",
                     edge_size = 0.5,
                     set_size = 5,
                     label_txtWidth = 1,
                     label_size = 5,
                     scale_size = 0.5)+scale_fill_distiller(palette = "RdBu")+
    theme_void()+
    #labs(title = paste0("Reference datasets ","(",p_name,")"))+
    theme(legend.position = 'right',
          plot.margin=unit(c(1, 1, 1, 1),'cm'),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15),)+scale_x_continuous(expand = expansion(mult = .2))
  ggsave(paste0(p_name,'_','venn_plot.png'),plot = vp)
  ggsave(paste0(p_name,'_','venn_plot.pdf'),plot = vp)
  
}


venn_plot(combine_df,'NIST_10X')
