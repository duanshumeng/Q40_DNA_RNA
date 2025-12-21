rm(list = ls())
setwd('/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/Data_performance')
source('/Users/duanshumeng/PGx_lab/HCC1395/Quartet_Multiomics_Ratio_Code/utils/theme_nature.r')

# Load function: get_test_type_with_effectsize
source("~/Gx_lab/HCC1395/ElementVSIllumina/Scripts/Visualization scripts/Statistics_analysis.R")

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
  theme(
    axis.text.x = element_text(face = "bold",color = "black"),
    axis.text.y = element_text(face = "bold",color = "black"),
    axis.title  = element_text(face = "bold",color = "black"),
    panel.grid.minor  =  element_blank(),
    legend.background =  element_blank(),
    legend.title = element_text(face = "bold",color = "black"),
    legend.text = element_text(color = "black")
  )

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

colors=c("#709AE1FF","#F05C3BFF")
#colors = c('#3783BB','#91569F','#C94741')
names(colors) <- c('Q30','Q40')

smp_colors = c('#F2E3EB' ,'#D6DEFF' ,'#FFEDC2')
names(smp_colors) <- c('RM8398','Quartet','HCC1395/BL')


# Function arguments example:
# lowvaf.df.sub,'InDel',c(0,0.5),d_type = 'F1.score'
# tnseq.df.sub=lowvaf.df.sub
# m_type='InDel'
# y_lim=c(0,0.5)
# d_type = 'F1.score'
get_boxcompare = function(tnseq.df.sub,m_type, y_lim,ct_1395.sum=data.frame(),d_type='F1.score'){
  tnseq.df.sub <- subset(tnseq.df.sub,type == m_type)
  tnseq.df.sub['F1.score'] <- tnseq.df.sub[d_type]
  tnseq.df.test <- tnseq.df.sub %>%
    group_by(type, Depth) %>%
    filter(n_distinct(F1.score) > 1, n_distinct(Quality) > 1) %>% 
    t_test(F1.score ~ Quality, p.adjust.method = 'none', comparisons = c("Q30", "Q40"),paired=FALSE) %>%  #t_test(F1.score ~ Quality, p.adjust.method = 'none', comparisons = c("Q30", "Q40"),paired=TRUE) %>%
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
    geom_boxplot(aes(fill=Quality),position = position_dodge(width = 0.8),
                 outlier.shape = NA,  
                 width = 0.6,alpha=0.35,size = 0.5
    )+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                size = 3, shape=21,colour='white',aes(fill = Quality))+
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
      my_theme+theme(legend.position = "none")+ylab('Count')+xlab('')
    #if (unique(tnseq.df.sub$Sample)=='HCC1395/BL'){
    #  box_p = box_p +annotation_custom(grob=ggplotGrob(ct_1395_plot),ymin = y_lim[1], ymax=y_lim[1]+0.2, xmin=2, xmax=4.5)
    #} else {
    #  box_p = box_p +annotation_custom(grob=ggplotGrob(ct_1395_plot),ymin = y_lim[1], ymax=y_lim[1]+0.3, xmin=3, xmax=8.5)
    #}
  }
  
  return(box_p)
}
# Variant count statistics -----
# ELE & ILM
df_count <- read.csv('ELE_ILM/WES_VCF_counts.vcf_stats.csv')
df_count$Quality <- lapply(str_split(df_count$source,'\\_',simplify = T)[,4],function(x){
  p = ifelse(grepl('ILM',x),'Q30_I','Q40_E')
  return(p)
}) %>% unlist()

df_count$Sample <- str_split(df_count$source,'\\_',simplify = T)[,2]
df_count$Sample[df_count$Sample %in% c('D5','D6','F7','M8')] = 'Quartet'
df_count$Source <- paste0(df_count$Sample,'.',df_count$Quality,'.',df_count$source)
df_count$Type <- gsub('SNVs','SNV',df_count$type) %>% gsub('Indels','InDel',.)
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


# Variant accuracy (F1 score) at different depths --------------
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

# High-confidence region -------
# NIST-RM8398
# ELE & ILM
nist.df <- read.csv("ELE_ILM/F1score_stats_total_NIST.csv")
nist.df.sub <- nist.df[,c('Type','source','METRIC.F1_Score')]
nist.df.sub$Type <- gsub('SNP','SNV',nist.df.sub$Type) %>% gsub('INDEL','InDel',.)
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
nist.df.sub$Sample <- 'RM8398'
nist.df.sub$Depth <- factor(nist.df.sub$Depth,levels = c('10X','15X','20X','25X','30X','60X','90X','120X'))
# Mutation count ----
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
ct_nist.sum$type <- gsub('Indels','InDel',ct_nist.sum$type) %>% gsub('SNVs','SNV',.)
nist.df.sub$type = nist.df.sub$Type
nist_snv <- get_boxcompare(tnseq.df.sub=nist.df.sub,m_type='SNV',y_lim=c(0.5,1.05),d_type='F1.score');nist_snv

Fig.S1a.nist.snv.test <- get_test_type_with_effectsize(subset(nist.df.sub,Type=='SNV'))

#hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.3,0.9),ct_1395.sum);hcc_indel
nist_indel <- get_boxcompare(tnseq.df.sub=nist.df.sub,m_type='InDel',y_lim=c(0.5,1.02),d_type='F1.score');nist_indel
Fig.S1b.nist.indel.test <- get_test_type_with_effectsize(subset(nist.df.sub,Type=='InDel'))

#p_qq <- plot_all_qqplots(subset(nist.df.sub,Type=='SNV'))

#ggsave(file='/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/文章/GB投稿/Reviewer_comments_reply/Figures/NIST.QQ_plot.pdf',p_qq,width = 8.15,height = 5.7)

# Example of function call:
# tnseq.df.sub=nist.df.sub
# m_type='SNV'
# y_lim=c(0.5,1.05)
# ct_1395.sum=ct_nist.sum
# d_type='F1.score'

## Quartet ----
# ELE & ILM
q.f1 <- read.csv("ELE_ILM/F1score_stats_total_Quartet.csv",sep = ',')
q.f1$source <- gsub('_10G_','.',q.f1$source)
q.f1.sub <- q.f1[,c('Type','source','METRIC.F1_Score')]
q.f1.sub$Type <- gsub('SNP','SNV',q.f1.sub$Type) %>% gsub('INDEL','InDel',.)
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
    Mean_Duplicates = mean(F1.score, na.rm = TRUE),   # Calculate mean
    SD_Duplicates = sd(F1.score, na.rm = TRUE)       # Calculate standard deviation
  )
# Mutation count ----
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
ct_quartet.sum$type <- gsub('Indels','InDel',ct_quartet.sum$type) %>% gsub('SNVs','SNV',.)

plot_count_F1 <- function(ct_1395_data, m_type, bar_width = 0.7, error_width = 0.2) {
  plot_data <- subset(ct_1395_data, type == m_type)
  
  # If data is empty, return an empty plot with a warning
  if (nrow(plot_data) == 0) {
    warning("No data found for type: ", m_type)
    return(ggplot() + 
             labs(title = paste("No data for type:", m_type)) + 
             my_theme)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Depth, y = count_median, fill = Quality)) +
    # Add gray background rectangles (based on original code logic)
    geom_rect(xmin = 0.5, xmax = 1.5, ymin = 0, ymax = max(plot_data$count_median) * 1.1, fill = "#E9E9E9", alpha = 0.3) +
    geom_rect(xmin = 2.5, xmax = 3.5, ymin = 0, ymax = max(plot_data$count_median) * 1.1, fill = "#E9E9E9", alpha = 0.3) +
    geom_rect(xmin = 4.5, xmax = 5.5, ymin = 0, ymax = max(plot_data$count_median) * 1.1, fill = "#E9E9E9", alpha = 0.3) +
    geom_rect(xmin = 6.5, xmax = 7.5, ymin = 0, ymax = max(plot_data$count_median) * 1.1, fill = "#E9E9E9", alpha = 0.3) +
    
    # Add bar plot
    geom_bar(stat = "identity", position = position_dodge(width = bar_width), width = bar_width) +
    
    # Add error bars
    geom_errorbar(
      aes(ymin = count_median - count_sd, ymax = count_median + count_sd),
      position = position_dodge(width = bar_width),
      width = error_width
    ) +
    
    # Set colors
    scale_fill_manual(values = colors) +
    
    # Apply custom theme
    my_theme +
    
    # Remove legend
    theme(legend.position = "none") +
    
    # Axis labels
    labs(
      x = "Depth",
      y = "Count",
      title = paste("Count Median by Depth for", m_type)
    )
  
  return(p)
}

quartet_snv.count <- plot_count_F1(ct_1395_data=ct_quartet.sum,m_type='SNV',bar_width = 0.7, error_width = 0.2)
quartet_indel.count <- plot_count_F1(ct_1395_data=ct_quartet.sum,m_type='InDel',bar_width = 0.7, error_width = 0.2)


q.f1.sub$type <- q.f1.sub$Type
#q.f1.sub$source <- str_split(q.f1.sub$source,'\\.WES',simplify = T)[,1]
q.f1.sub$source <- gsub('_ELE|_ILM','',q.f1.sub$source)
# Example function call:
# tnseq.df.sub=nist.df.sub,m_type='InDel',y_lim=c(0.5,1.02),d_type='F1.score'
quartet_snv <- get_boxcompare(tnseq.df.sub=q.f1.sub,m_type='SNV',y_lim=c(0.6,1.03),d_type = 'F1.score');quartet_snv
#hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.3,0.9),ct_1395.sum);hcc_indel
quartet_indel <- get_boxcompare(q.f1.sub,'InDel',c(0.5,1.0));quartet_indel

Fig.3c.quartet_snv.test <- get_test_type_with_effectsize(subset(q.f1.sub,Type=='SNV'))
Fig.3d.quartet_indel.test <- get_test_type_with_effectsize(subset(q.f1.sub,Type=='InDel'))
quartet_snv.QQ <- plot_all_qqplots(subset(q.f1.sub,Type=='SNV'))

subset(q.f1.sub,Depth=='25X') %>%
  group_by(Quality,Type) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(F1.score, na.rm = TRUE),   # Calculate mean
    SD_Duplicates = sd(F1.score, na.rm = TRUE)       # Calculate standard deviation
  )


#ggsave(file='/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/文章/GB投稿/Reviewer_comments_reply/Figures/Quartet_QQ_plot.pdf',quartet_snv.QQ,width = 8.15,height = 5.7)

# Genome-wide analysis --
# ELE & ILM
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

q.mcr$type <- q.mcr$Type %>% gsub('INDEL','InDel',.)
q.mcr$Depth <- factor(q.mcr$Depth,levels = c('10X','15X','20X','25X','30X','60X','90X','120X'))
quartet_snv.mcr <- get_boxcompare(q.mcr,'SNV',c(0.4,1.03),d_type = 'MCR');quartet_snv.mcr
#hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.3,0.9),ct_1395.sum);hcc_indel
quartet_indel.mcr <- get_boxcompare(q.mcr,'InDel',c(0.4,0.9),d_type = 'MCR');quartet_indel.mcr

q.mcr$F1.score <- q.mcr$MCR
Fig.3a.q_mcr.snv.test <- get_test_type_with_effectsize(subset(q.mcr,type=='SNV'))
Fig.3b.q_mcr.indel.test <- get_test_type_with_effectsize(subset(q.mcr,type=='InDel'))

subset(q.mcr,Depth=='30X') %>%
  group_by(Quality,Type) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(MCR, na.rm = TRUE),   # Calculate mean
    SD_Duplicates = sd(MCR, na.rm = TRUE)       # Calculate standard deviation
  )


# Wilcoxon rank-sum test for 120x:
# SNV
wilcox.test(subset(q.mcr,type=='SNV'&Depth=='120X'&Quality=='Q40')$MCR, subset(q.mcr,type=='SNV'&Depth=='120X'&Quality=='Q30')$MCR,
            paired = FALSE, alternative = "two.sided", 
            exact = FALSE, conf.int = FALSE)

calculate_unpaired_cohens_d(subset(q.mcr,type=='SNV'&Depth=='120X'&Quality=='Q30')$MCR,subset(q.mcr,type=='SNV'&Depth=='120X'&Quality=='Q40')$MCR)
# InDel
wilcox.test(subset(q.mcr,type=='InDel'&Depth=='120X'&Quality=='Q40')$MCR, subset(q.mcr,type=='InDel'&Depth=='120X'&Quality=='Q30')$MCR,
            paired = FALSE, alternative = "two.sided", 
            exact = FALSE, conf.int = FALSE)

calculate_unpaired_cohens_d(subset(q.mcr,type=='InDel'&Depth=='120X'&Quality=='Q30')$MCR,subset(q.mcr,type=='InDel'&Depth=='120X'&Quality=='Q40')$MCR)

get_MCR_count_plot <- function(df, v_type) {
  error_type = "sd"
  # error_type: specifies error type, 'sd' for standard deviation, 'se' for standard error
  
  df_long <- df %>%
    pivot_longer(
      cols = c("Total_Variants", "Mendelian_Concordant_Variants"),
      names_to = "Variant_Type",
      values_to = "Count"
    )
  
  # --- Key modification: calculate mean and error for each group (assuming input is raw replicate data) ---
  # Note: If your df is already summarized (with error columns), skip this step and use your data directly
  # Here we assume df_long has multiple observations per group, and we need to calculate mean and error
  plot_summary <- df_long %>%
    group_by(Quality, Depth, Variant_Type) %>%
    summarise(
      mean_count = mean(Count, na.rm = TRUE),
      sd_count = sd(Count, na.rm = TRUE),
      se_count = sd_count / sqrt(n()), # Calculate standard error
      n = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      # Determine error bar limits based on selected error type
      error_lower = mean_count - ifelse(error_type == "sd", sd_count, se_count),
      error_upper = mean_count + ifelse(error_type == "sd", sd_count, se_count)
    )
  
  variant_colors <- c("steelblue", "lightcoral")
  
  # Create clustered bar chart (using calculated mean)
  p <- ggplot(plot_summary, aes(x = Quality, y = mean_count, fill = Variant_Type)) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      size = 0.3
    ) +
    # --- New: error bars (geom_errorbar) ---
    geom_errorbar(
      aes(ymin = error_lower, ymax = error_upper),
      position = position_dodge(width = 0.8),
      width = 0.25,        # Width of error bar end lines
      size = 0.4,          # Error bar line thickness
      color = "black"      # Error bar color
    ) +
    facet_wrap(~ Depth, ncol = 4, scales = "free_y") +
    scale_fill_manual(
      values = variant_colors,
      labels = c("Mendelian Concordant Variants", "Total Variants")
    ) +
    labs(
      x = "Sequencing Quality",
      y = "Mean Variant Count", # Y-axis label changed to "Mean"
      fill = "Variant Type",
      title = v_type
    ) +
    my_theme +
    # Adjust text labels: display mean
    geom_text(
      aes(label = format(round(mean_count), big.mark = ",")),
      position = position_dodge(width = 0.8),
      vjust = -1.5,         # Move text further up to make room for error bars
      size = 2.8,
      color = "black"
    ) +
    # Adjust y-axis range to leave more space for error bars and value labels
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  return(p)
}
q.mcr.count <- q.mcr[,c('Total_Variants','Mendelian_Concordant_Variants','Quality','Type','Depth')]
quartet_snv.mcr.count <- get_MCR_count_plot(subset(q.mcr.count,Type=='SNV'),'SNV')
quartet_indel.mcr.count <- get_MCR_count_plot(subset(q.mcr.count,Type=='INDEL'),'InDel')

quartet_var_count <- quartet_snv.mcr.count/quartet_indel.mcr.count/(quartet_snv.count+quartet_indel.count)+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect',widths = c(2,2,1)) & theme(legend.position='bottom',
                                                          plot.tag = element_text(color = "black",size = 16,face = "bold"))
ggsave('Fig.wes_sup.germline-count.pdf',quartet_var_count,width=8.15, height=5.7*2,dpi = 300)


get_boxcompare_rmsk = function(tnseq.df.sub,m_type, y_lim,ct_1395.sum=data.frame(),d_type='F1.score',xmax=5.5){
  tnseq.df.sub <- subset(tnseq.df.sub,type == m_type)
  tnseq.df.sub['F1.score'] <- tnseq.df.sub[d_type]
  
  #tnseq.df.test$y.position <- ifelse(tnseq.df.test$group2=='Q40',tnseq.df.test$y.position-0.01,tnseq.df.test$y.position-0.03)
  tnseq.df.sub$group <- paste0(tnseq.df.sub$Quality,'_',tnseq.df.sub$Depth)
  if (unique(tnseq.df.sub$Sample)=='HCC1395/BL' & m_type=='InDel'){
    box_p <- ggplot(tnseq.df.sub, aes(x = Depth, y = F1.score)) +
      geom_rect(xmin = 0.5, xmax = 1.5,
                ymin = 0, ymax = 11, fill ="#E9E9E9") +
      geom_rect(xmin = 2.5, xmax = 3.5,
                ymin = 0, ymax = 11, fill ="#E9E9E9")+
      geom_rect(xmin = 4.5, xmax = 5.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9") +
      geom_rect(xmin = 6.5, xmax = 7.5,
                ymin = 0, ymax = 2000000, fill ="#E9E9E9")+
      geom_boxplot(aes(fill=Quality),
                   outlier.shape = NA,  # Hide outliers
                   width = 0.6,alpha=0.35,size = 0.5
      )+
      geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 3,shape=21,colour='white' ,aes(fill = Quality))+
      #geom_line(aes(color = Quality,group=Quality),position = position_dodge(width = 0.5),alpha=0.5,size=1)+
      scale_color_manual(values = colors)+
      scale_fill_manual(values = colors)+
      #facet_grid(type~.)+
      ylim(y_lim)+ylab(paste0(d_type,' ',m_type))+ggtitle(unique(tnseq.df.sub$Sample))+
      my_theme+theme(axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 15,hjust = 1))+xlab('')
    #box_p = box_p +annotation_custom(grob=ggplotGrob(ct_1395_plot),ymin = y_lim[1], ymax=y_lim[1]+0.3, xmin=1, xmax=xmax)
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
      geom_boxplot(aes(fill=Quality),
                   outlier.shape = NA,  # Hide outliers
                   width = 0.6,alpha=0.35,size = 0.5
      )+
      geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 3,shape=21,colour='white',aes(fill = Quality))+
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
    #box_p = box_p +annotation_custom(grob=ggplotGrob(ct_1395_plot),ymin = y_lim[1], ymax=y_lim[1]+0.4, xmin=2, xmax=xmax)
    
  }
  return(box_p)
}
# Complex genomic regions -------
get_rmsk_count_compare <- function(rmsk_count,Sample,v_type){
  trend_data <- rmsk_count %>%
    dplyr::group_by(Rmsk_class, Quality) %>%
    dplyr::summarise(
      mean_QUERY_TOTAL = mean(QUERY.TOTAL),
      mean_TRUTH_TP = mean(TRUTH.TP),
      .groups = 'drop'
    )
  
  p_rmsk = ggplot(rmsk_count, aes(x = TRUTH.TP, y = QUERY.TOTAL)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
               size = 3, shape=21,colour='white',aes(fill = Quality))+
    facet_wrap(~ Rmsk_class, scales = "free") + # Facet by repeat type
    labs(
      title = Sample,
      subtitle = v_type,
      x = "TRUTH.TP",
      y = "QUERY.TOTAL",
      color = "Quality"
    ) +
    scale_fill_manual(values = colors)+
    my_theme+
    theme(legend.position = "bottom")
  return(p_rmsk)
}

# NIST ---
# ELE & ILM
df_nist = read.csv("ELE_ILM/F1score_stats_total_NIST_rmsk.csv")
df_nist.sub <- df_nist[,c('Type','source','METRIC.F1_Score','QUERY.TOTAL','TRUTH.TP','TRUTH.FN','TRUTH.TOTAL','QUERY.FP','QUERY.UNK')]
df_nist.sub$Type <- gsub('SNP','SNV',df_nist.sub$Type)
df_nist.sub$Rmsk_class <- str_split(df_nist.sub$source,'\\.',simplify = T)[,2]
df_nist.sub$F1.score <- df_nist.sub$METRIC.F1_Score
#df_nist.sub$Platform <- str_split(df_nist.sub$source,'\\_',simplify = T)[,4]
df_nist.sub$Quality <- lapply(str_split(df_nist.sub$source,'\\_',simplify = T)[,4],function(x){
  p = ifelse(grepl('ILM',x),'Q30','Q40')
  return(p)
}) %>% unlist()

df_nist.sub$Sample <- 'RM8398'
df_nist.sub$count <- df_nist.sub$QUERY.TOTAL

df_nist.sub$type <- df_nist.sub$Type %>% gsub('INDEL','InDel',.)
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
# Example of function call:
# tnseq.df.sub = df_nist.sub
# m_type = 'SNV'
# y_lim=c(0.4,1.03)
# ct_1395.sum=df_nist.sub.rmsk
# d_type='F1.score'

nist_snv.rmsk <- get_boxcompare_rmsk(df_nist.sub,'SNV',c(0.5,1.08),df_nist.sub.rmsk,d_type = 'F1.score',xmax = 6.5);nist_snv.rmsk
#hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.3,0.9),ct_1395.sum);hcc_indel
nist_indel.rmsk <- get_boxcompare_rmsk(df_nist.sub,'InDel',c(0.3,1.1),df_nist.sub.rmsk,d_type = 'F1.score');nist_indel.rmsk

Fig.S1c.nist.rmsk.snv.test <- get_test_type_with_effectsize(subset(df_nist.sub,type=='SNV'))
Fig.S1d.nist.rmsk.indel.test <- get_test_type_with_effectsize(subset(df_nist.sub,type=='InDel'))


df_nist.sub.count <- df_nist.sub[,c('Rmsk_class','Quality','QUERY.TOTAL','TRUTH.TP','type')]
df_nist.sub.count.snv <- get_rmsk_count_compare(subset(df_nist.sub.count,type=='SNV'),'RM8398','SNV')
df_nist.sub.count.indel <- get_rmsk_count_compare(subset(df_nist.sub.count,type=='InDel'),'RM8398','SNV')

#df_nist.sub$Lib <- 'Lib1'

# Quartet F1 ---
# ELE & ILM
df_quartet = read.csv("ELE_ILM/Quartet.rmsk.F1.csv")
df_quartet.sub <- df_quartet[,c('Type','source','METRIC.F1_Score','QUERY.TOTAL','TRUTH.TP')]
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
df_quartet.sub <- subset(df_quartet.sub,!grepl('WES_D6_2',df_quartet.sub$source)) 
df_quartet.sub$type <- gsub('INDEL','InDel',df_quartet.sub$type)

quartet_snv.rmsk <- get_boxcompare_rmsk(df_quartet.sub,'SNV',c(0.5,1.08),ct_quartet.rmsk,d_type = 'F1.score',xmax = 6.5);quartet_snv.rmsk
#hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.3,0.9),ct_1395.sum);hcc_indel
quartet_indel.rmsk <- get_boxcompare_rmsk(df_quartet.sub,'InDel',c(0.5,1.1),ct_quartet.rmsk,d_type = 'F1.score');quartet_indel.rmsk

Fig.3e.q.quartet_snv.rmsk.test <- get_test_type_with_effectsize(subset(df_quartet.sub,type=='SNV'))
Fig.3f.q.quartet_indel.rmsk.test <- get_test_type_with_effectsize(subset(df_quartet.sub,type=='InDel'))

quartet_snv.rmsk.count <- df_quartet.sub[,c('QUERY.TOTAL','TRUTH.TP','type','Rmsk_class','Quality')]

quartet_snv.rmsk.count.snv <- get_rmsk_count_compare(subset(quartet_snv.rmsk.count,type=='SNV'),'Quartet','SNV (F1.score)')
quartet_snv.rmsk.count.indel <- get_rmsk_count_compare(subset(quartet_snv.rmsk.count,type=='InDel'),'Quartet','InDel (F1.score)')

# Quartet MCR ---
# ELE & ILM
df_quartet.mcr = read.csv("ELE_ILM/Quartet.rmsk.MCR.csv")
df_quartet.mcr$type <- str_split(df_quartet.mcr$Family,'\\.',simplify = T)[,3] %>% gsub('INDEL','InDel',.)
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
#hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.3,0.9),ct_1395.sum);hcc_indel
quartet_indel.rmsk.mcr <- get_boxcompare_rmsk(df_quartet.mcr,'InDel',c(0.2,1),ct_quartet.rmsk,d_type = 'MCR');quartet_indel.rmsk.mcr

df_quartet.mcr$F1.score <- df_quartet.mcr$MCR
df_quartet.mcr$source <- df_quartet.mcr$Family
Fig.3e.q.rmsk.mcr.snv.test <- get_test_type_with_effectsize(subset(df_quartet.mcr,type=='SNV'))
Fig.3f.q.rmsk.mcr.indel.test <- get_test_type_with_effectsize(subset(df_quartet.mcr,type=='InDel'))

df_quartet.mcr.rsmk_count <- df_quartet.mcr[,c('Rmsk_class','Quality','Total_Variants','Mendelian_Concordant_Variants','type')]
df_quartet.mcr.rsmk_count$QUERY.TOTAL <- df_quartet.mcr.rsmk_count$Total_Variants
df_quartet.mcr.rsmk_count$TRUTH.TP <- df_quartet.mcr.rsmk_count$Mendelian_Concordant_Variants
df_quartet.mcr.rsmk_count.snv <- get_rmsk_count_compare(subset(df_quartet.mcr.rsmk_count,type=='SNV'),'Quartet','SNV (MCR)')
df_quartet.mcr.rsmk_count.indel <- get_rmsk_count_compare(subset(df_quartet.mcr.rsmk_count,type=='InDel'),'Quartet','InDel (MCR)')


for (v_type in unique(df_quartet.mcr$type)){
  for (rmsk in unique(df_quartet.mcr$Rmsk_class)){
    print(paste0(v_type,' ',rmsk))
    df_rmsk <- subset(df_quartet.mcr,type==v_type & Rmsk_class==rmsk)
    wilcox.test(MCR ~ Quality, data = df_rmsk) %>% print()
    calculate_unpaired_cohens_d(subset(df_rmsk,Quality=='Q30')$MCR,subset(df_rmsk,Quality=='Q40')$MCR) %>% print()
    
  }
}

results_list <- list()

for (v_type in unique(df_quartet.mcr$type)){
  for (rmsk in unique(df_quartet.mcr$Rmsk_class)){
    
    cat("Processing:", v_type, "-", rmsk, "\n")
    
    df_rmsk <- subset(df_quartet.mcr, type == v_type & Rmsk_class == rmsk)
    
    if(nrow(df_rmsk) < 2) {
      cat("  Skipped: insufficient data\n")
      next
    }
    
    wilcox_result <- wilcox.test(MCR ~ Quality, data = df_rmsk)
    
    cohens_d_result <- calculate_unpaired_cohens_d(
      subset(df_rmsk, Quality == 'Q30')$MCR,
      subset(df_rmsk, Quality == 'Q40')$MCR
    )
    
    result_row <- data.frame(
      Variant_Type = v_type,
      Rmsk_Class = rmsk,
      W_statistic = wilcox_result$statistic,
      Wilcoxon_p_value = wilcox_result$p.value,
      Cohens_d = cohens_d_result$cohens_d,
      CI_lower = cohens_d_result$ci_lower,
      CI_upper = cohens_d_result$ci_upper,
      Q30_n = sum(df_rmsk$Quality == 'Q30'),
      Q40_n = sum(df_rmsk$Quality == 'Q40'),
      Total_n = nrow(df_rmsk)
    )
    results_list[[paste(v_type, rmsk, sep = "_")]] <- result_row
  }
}

results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL
print(results_df)


subset(df_quartet.mcr,Sample=='Quartet') %>%
  group_by(Quality,Depth,type) %>%
  dplyr::summarise(
    Mean_Duplicates = mean(MCR, na.rm = TRUE),   # Calculate mean
    SD_Duplicates = sd(MCR, na.rm = TRUE)       # Calculate standard deviation
  )



# HCC1395. High-confidence region ------
# ELE_ILM
tnseq.df <- read.csv("ELE_ILM/HCC1395_F1score_stats_total.SEQC2.csv")
#tnseq.df <- read.csv("ELE_ILM/HCC1395_F1score_stats_total.csv")
#tnseq.df <- read.csv("ELE_ILM/HCC1395_F1score_stats_total.strelka2.csv")
tnseq.df.sub <- tnseq.df[,c('type','source','F1.score','recall','precision')]
tnseq.df.sub$type <- gsub('SNVs','SNV',tnseq.df.sub$type) %>% gsub('indels','InDel',.)
tnseq.df.sub$Depth <- str_split(tnseq.df.sub$source,'\\.',simplify = T)[,2]
#tnseq.df.sub$Platform <- str_split(tnseq.df.sub$source,'\\_',simplify = T)[,4]
# Mutation count ----
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
ct_1395.sum$type <- gsub('Indels','InDel',ct_1395.sum$type) %>% gsub('SNVs','SNV',.)
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


hcc_snv <- get_boxcompare(tnseq.df.sub,'SNV',c(0.45,0.9));hcc_snv
#hcc_snv <- get_boxcompare(tnseq.df.sub,'SNV',c(0.7,1.0),ct_1395.sum);hcc_snv
hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.3,0.8));hcc_indel
#hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.4,0.85),ct_1395.sum);hcc_indel
Fig.4a.hcc_snv.test <- get_test_type_with_effectsize(subset(tnseq.df.sub,Type=='SNV'))
Fig.4b.hcc_indel.test <- get_test_type_with_effectsize(subset(tnseq.df.sub,Type=='InDel'))

#write.csv(Fig.4a.hcc_snv.test,file = 'Fig.4a.hcc_snv.test.csv',row.names = F)
#write.csv(Fig.4b.hcc_indel.test,file = 'Fig.4b.hcc_indel.test.csv',row.names = F)

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


## HCC1395 complex genomic regions ------
# HCC1395 ---
# ELE & ILM
df_hcc1395 = read.csv("ELE_ILM/F1score_stats_total_HCC1395_rmsk.SEQC2.csv")
df_hcc1395$QUERY.TOTAL <- df_hcc1395$total.query
df_hcc1395$TRUTH.TP <- df_hcc1395$tp
df_hcc1395.sub <- df_hcc1395[,c('type','source','F1.score','QUERY.TOTAL','TRUTH.TP')]
df_hcc1395.sub$Type <- df_hcc1395.sub$type <- gsub('SNVs','SNV',df_hcc1395.sub$type) %>% gsub('indels','InDel',.)
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
# Example of function call:
# tnseq.df.sub = df_hcc1395.sub
# m_type = 'InDel'
# y_lim=c(0.1,1.03)
# ct_1395.sum=rmsk_df_median.hcc
# d_type='F1.score'

hcc_snv.rmsk <- get_boxcompare_rmsk(tnseq.df.sub=df_hcc1395.sub,m_type='SNV',y_lim=c(0.2,1.1),xmax = 4.5);hcc_snv.rmsk
#hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.3,0.9),ct_1395.sum);hcc_indel
hcc_indel.rmsk <- get_boxcompare_rmsk(df_hcc1395.sub,'InDel',c(0.6,1),xmax = 0.6);hcc_indel.rmsk

Fig.4c.hcc_snv.rmsk.test <- get_test_type_with_effectsize(subset(df_hcc1395.sub,Type=='SNV'))
Fig.4d.hcc_indel.rmsk.test <- get_test_type_with_effectsize(subset(df_hcc1395.sub,Type=='InDel'))

df_hcc1395.sub.rmsk_count <- df_hcc1395.sub[,c('type','Quality','Rmsk_class','QUERY.TOTAL','TRUTH.TP')]
df_hcc1395.sub.rmsk_count.snv <- get_rmsk_count_compare(subset(df_hcc1395.sub.rmsk_count,type=='SNV'),'','SNV')
df_hcc1395.sub.rmsk_count.indel <- get_rmsk_count_compare(subset(df_hcc1395.sub.rmsk_count,type=='InDel'),'','InDel')


# WES targeted capture region repetitive elements
rep_df <- read.csv('/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/文章/Q40_repeat-information.csv')
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



# HCC1395 low-frequency variant region ------
# Low-frequency variant F1 score ------
# ELE & ILM
#lowvaf.df <- read.csv("ELE_ILM/F1score_stats_total.LowVAF.csv")
lowvaf.df <- read.csv("ELE_ILM/F1score_stats_total-lowVAF.SEQC2.csv")
lowvaf.df <- subset(lowvaf.df,grepl('120X',lowvaf.df$source))
#F1score_stats_total-lowVAF.hcc1395.csv
#lowvaf.df <- read.csv("/Users/duanshumeng/fsdownload/F1score_stats_total-lowVAF.hcc1395.csv")
lowvaf.df.sub <- lowvaf.df[,c('type','source','F1.score','recall','precision')]
lowvaf.df.sub$type <- gsub('SNVs','SNV',lowvaf.df.sub$type) %>% gsub('indels','InDel',.)
lowvaf.df.sub$Depth <- str_split(lowvaf.df.sub$source,'\\.',simplify = T)[,2]
#lowvaf.df.sub$Platform <- str_split(lowvaf.df.sub$source,'\\_',simplify = T)[,4]
lowvaf.df.sub <- subset(lowvaf.df.sub,!grepl('Mix',lowvaf.df.sub$source))
lowvaf.df.sub$Quality <- lapply(str_split(lowvaf.df.sub$source,'\\_',simplify = T)[,4],function(x){
  p = ifelse(grepl('ILM',x),'Q30','Q40')
  return(p)
}) %>% unlist()
lowvaf.df.sub$Type <- lowvaf.df.sub$type
lowvaf.df.sub$Depth <- factor(lowvaf.df.sub$Depth,levels=c('30X','60X','90X','120X'))

hcc_snv.lowVAF <- get_boxcompare_lowVAF_indel(tnseq.df.sub=lowvaf.df.sub,m_type='SNV',y_lim=c(0,0.7),d_type = 'F1.score');hcc_snv.lowVAF
#hcc_indel <- get_boxcompare(tnseq.df.sub,'InDel',c(0.3,0.9),ct_1395.sum);hcc_indel
hcc_indel.lowVAF <- get_boxcompare_lowVAF_indel(lowvaf.df.sub,'InDel',c(0,0.5),d_type = 'F1.score');hcc_indel.lowVAF

Fig.4e.hcc_snv.lowVAF.test <- get_test_type_with_effectsize(subset(lowvaf.df.sub,Type=='SNV'))
Fig.4f.hcc_indel.lowVAF.test <- get_test_type_with_effectsize(subset(lowvaf.df.sub,Type=='InDel'))


get_boxcompare_lowVAF_indel <- function(tnseq.df.sub,m_type,y_lim,d_type = 'F1.score'){
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
  tnseq.df.test$xmin  <- 0.8
  tnseq.df.test$xmax  <- 1.2
  tnseq.df.test$x <- 1.0
  #tnseq.df.test$y.position <- ifelse(tnseq.df.test$group2=='Q40',tnseq.df.test$y.position-0.01,tnseq.df.test$y.position-0.03)
  tnseq.df.sub$group <- paste0(tnseq.df.sub$Quality,'_',tnseq.df.sub$Depth)
  box_p <- ggplot(tnseq.df.sub, aes(x = Depth, y = F1.score)) +
    geom_boxplot(aes(fill=Quality),position = position_dodge(width = 0.8),
                 outlier.shape = NA,  # Hide outliers
                 width = 0.6,alpha=0.35,size = 0.5
    )+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                size = 3, shape=21,colour='white',aes(fill = Quality))+
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
  return(box_p)
}

hcc_snv.recall.lowVAF <- get_boxcompare(tnseq.df.sub=lowvaf.df.sub,m_type='SNV',y_lim=c(0,0.5),d_type = 'recall');hcc_snv.recall.lowVAF
hcc_indel.recall.lowVAF <- get_boxcompare(lowvaf.df.sub,'InDel',c(0,0.5),d_type = 'recall');hcc_indel.recall.lowVAF


lowvaf.df.sub %>% 
  group_by(type,Quality) %>%
  dplyr::summarise(f1_mean=mean(F1.score))

hcc_snv.precision.lowVAF <- get_boxcompare_lowVAF_indel(tnseq.df.sub=lowvaf.df.sub,m_type='SNV',y_lim=c(0,0.5),d_type = 'precision');hcc_snv.precision.lowVAF
hcc_indel.precision.lowVAF <- get_boxcompare_lowVAF_indel(lowvaf.df.sub,'InDel',c(0,0.8),d_type = 'precision');hcc_indel.precision.lowVAF

hcc_snv.recall.lowVAF+hcc_indel.recall.lowVAF+hcc_snv.precision.lowVAF+hcc_indel.precision.lowVAF

df=subset(lowvaf.df.sub,Type!='records')[,c('F1.score','recall','precision','Depth','Type','Quality')]
library(tidyverse)
library(ggplot2)
library(ggrepel)
# Data formatting: convert to long format
df_long <- df %>%
  pivot_longer(
    cols = c(F1.score, recall, precision),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Depth = factor(Depth, levels = c("30X", "60X", "90X", "120X")),
    Type = factor(Type,levels=c('SNV','InDel')),
    Quality = factor(Quality,levels=c('Q30','Q40')),
    Metric = factor(Metric, levels = c("F1.score", "recall", "precision"),
                    labels = c("F1 Score", "Recall", "Precision"))
  )

# Calculate mean for each depth and metric (for connecting lines)
df_mean <- df_long %>%
  group_by(Depth, Metric,Type,Quality) %>%
  dplyr::summarise(Mean_Value = mean(Value), .groups = "drop")

# Create three-line plot
low_vaf.p <- ggplot() +
  # Add original data points
  geom_point(data = df_long, 
             aes(x = Depth, y = Value, color = Metric),
             position = position_jitter(width = 0.1), 
             alpha = 0.6, size = 2.5) +
  
  # Add mean connecting lines
  geom_line(data = df_mean, 
            aes(x = Depth, y = Mean_Value, color = Metric, group = Metric),
            linewidth = 1.2) +
  
  # Add mean points
  geom_point(data = df_mean, 
             aes(x = Depth, y = Mean_Value, color = Metric), 
             size = 4, shape = 18) +
  
  # Add mean labels
  geom_text_repel(data = df_mean,
                  aes(x = Depth, y = Mean_Value, label = sprintf("%.2f", Mean_Value), color = Metric),
                  size = 4, fontface = "bold", box.padding = 0.5) +
  
  # Set colors and theme
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) + # Red-Blue-Green
  labs(
    x = "Depth",
    y = "Value") +
  my_theme+
  # Add background grid for better readability
  geom_hline(yintercept = seq(0.1, 0.4, by = 0.05), color = "grey92", linetype = "dashed")+facet_grid(Quality~Type)

ggsave('Fig.somatic_sup2.pdf',low_vaf.p,width=8.15*1.5, height=5.7*1.5,dpi = 300)
ggsave('Fig.somatic_sup2.tiff',low_vaf.p,width=8.15*1.5, height=5.7*1.5,dpi = 300)

lowvaf.df.sub %>% 
  group_by(Quality,Type) %>% 
  dplyr::summarise(F1.score_mean = mean(F1.score,na.rm = T),
                   F1.score_sd = sd(F1.score,na.rm = T))

# Supplementary figures needed ----
# HCC1395 VAF correlation with reference dataset ------
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
  geom_bar(stat="summary",fun=mean,position=position_dodge(0.9),alpha=0.9,color='black')+
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


# Variant reproducibility at different depths --------------
# NIST
# ELE & ILM
df_nist = read.csv('ELE_ILM/NIST_Reproducibility_stats_total.csv') 
df_nist$Reproducibility <- (df_nist$TRUTH.TP/df_nist$TRUTH.TOTAL + df_nist$TRUTH.TP/df_nist$QUERY.TOTAL)/2
df_nist$Source <- gsub('_10G','',df_nist$source)
df_nist$Depth <- str_split(df_nist$Source,'\\.',simplify = T)[,2]
df_nist$Sample <- 'RM8398'
add_batch.e <- function(df){
  df[(str_count(df$source,'ELE')==2),'Quality'] = 'Q40'
  df[(str_count(df$source,'ILM')==2),'Quality'] = 'Q30'
  df[(str_count(df$source,'ILM')==1),'Quality'] = 'Q40_vs_Q30'
  return(df)
}
df_nist <- add_batch.e(df_nist)
df_nist$type <- gsub('SNP','SNV',df_nist$Type)  %>% gsub('INDEL','InDel',.)

df_nist$Depth <- factor(df_nist$Depth,levels = c('10X','15X','20X','25X','30X','60X','90X','120X'))

# Quartet
# ELE & ILM
df_quartet = read.csv('ELE_ILM/Quartet_Reproducibility_stats_total.csv')
df_quartet$Reproducibility <- (df_quartet$TRUTH.TP/df_quartet$TRUTH.TOTAL + df_quartet$TRUTH.TP/df_quartet$QUERY.TOTAL)/2
df_quartet$type <- gsub('SNP','SNV',df_quartet$Type) %>% gsub('INDEL','InDel',.)
df_quartet$Source <- gsub('_10G','',df_quartet$source)

df_quartet$Depth <- str_split(df_quartet$Source,'\\.',simplify = T)[,2]
df_quartet$Sample <- 'Quartet'
df_quartet <- add_batch.e(df_quartet)
df_quartet <- df_quartet[!apply(df_quartet, 1, function(x) any(grepl("D6_2_ELE_10G.120X.Haplotyper", x))), ]

# HCC1395
# ELE & ILM
hcc1395_dirs <- dir(path = 'ELE_ILM',pattern = 'HCC1395_F1score_stats_total.')
df_lst = list()
for (d in hcc1395_dirs){
  print(d)
  dep = str_split(d,'\\.',simplify = T)[,2]
  df = read.csv(paste0('ELE_ILM/',d))
  df$type <- gsub('indels','InDel',df$type) %>% gsub('SNVs','SNV',.)
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
rep_all <- subset(rep_all,Quality!='Q40_vs_Q30' & Quality!='Q40_vs_Q30_M' & type %in% c('InDel','SNV'))
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

rep_all_median$Sample <- gsub('NIST','RM8398',rep_all_median$Sample)
rep.p.somatic <- get_dotbar(subset(rep_all_median,Sample=='HCC1395/BL'),'Reproducibility','Reproducibility',y_lim=c(0.5,1))
rep.p.somatic
rep.p.germline <- get_dotbar(subset(rep_all_median,Sample!='HCC1395/BL'),'Reproducibility','Reproducibility',y_lim=c(0.58,1))
rep.p.germline


g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

germline_legend <- g_legend(f1.p.germline)  %>% plot()


# CNV ------
get_dotbar <- function(df,d_type,m_type,y_title,y_lim,scale_y_continuous=FALSE){
  df[,'mean'] = df[,paste0(d_type,'_mean')]
  df[,'sd'] = df[,paste0(d_type,'_sd')]
  df <- subset(df,cnv_type == m_type)
  CNV_Q40 =ggplot(data=df,aes(x=Depth,y=mean,fill=Quality))+
    geom_rect(xmin = 0.5, xmax = 1.5,
              ymin = 0, ymax = 30, fill ="#E9E9E9") +
    geom_rect(xmin = 2.5, xmax = 3.5,
              ymin = 0, ymax = 30, fill ="#E9E9E9")+
    geom_bar(stat="identity",position=position_dodge(0.8),width = 0.8,color='black',alpha=0.9)+
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),color='black',width=0.2,
                  position=position_dodge(0.8))+
    scale_x_discrete(breaks = c("30X", "60X", "90X", "120X"))+
    #facet_grid2(.~cnv_type,scales = "free",space = "free_x",
    #            strip  = strip_nested(
    #              background_x = elem_list_rect(fill =  
    #                                              c('#E2F3EB' ,'#D6DEFF' ,'#FFEDC2'),color = NA),
    #              background_y = elem_list_rect(fill =
    #                                              c('#F2E3EB'),color = NA)))+
    scale_fill_manual(values = colors)+
    theme_bw()+
    my_theme+
    labs(x="",y=y_title)+
    expand_limits(y = y_lim)+
    guides(fill=guide_legend(title = "",nrow = 1,byrow = FALSE))
  
  if (scale_y_continuous){
    CNV_Q40 = CNV_Q40+ scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)),
                                          labels = scales::unit_format(unit = "K", scale = 1e-3))
  }
  
  CNV_Q40;return(CNV_Q40)
  
}
# CNV Region --------------------
# ELE & ILM
all_region <- data.frame(bed=1,region=1,cnv_type=1,Quality=1)

for (i in dir('ELE_ILM/Somatic_CNV',pattern = '*bed.jaccard.txt')){
  print(i)
  cnv_type = str_split(i,'\\.',simplify = T)[,3]
  reg_df <- read.table(paste0('ELE_ILM/Somatic_CNV/',i),header = T)
  bed_reg <- data.frame(bed = c(reg_df$bed1,reg_df$bed2),
                        region = c(reg_df$bed1_region,reg_df$bed2_region)) %>% distinct(.keep_all = FALSE)
  bed_reg$bed <- gsub('WES_','',bed_reg$bed) %>% gsub('_10G','',.)
  bed_reg$Quality 
  bed_reg$Quality <- lapply(bed_reg$bed,function(x){
    p='Q40'
    if (grepl('ILM',x)){
      p='Q30'
    } else if (grepl('Mix',x)){
      p = 'Q35'
    }
    return(p)
  }) %>% unlist()
  bed_reg$cnv_type = cnv_type
  all_region = rbind(all_region,bed_reg)
}

all_region <- all_region[-1,]
all_region$Depth <- str_split(all_region$bed,'\\_',simplify = T)[,4] %>% factor(., levels=c("30X", "60X", "90X", "120X"))

all_region$region_raw <- all_region$region
all_region$region <- all_region$region_raw/1000000
all_region_median <- all_region %>%
  dplyr::group_by(cnv_type, Quality, Depth) %>%
  dplyr::summarize(
    across(
      region,
      list(
        mean = ~mean(.),
        sd = ~sd(.)
      )
    )
  )

all_region_median$Depth <- factor(all_region_median$Depth,levels=c("30X", "60X", "90X", "120X"))
all_region_median$Sample <- 'HCC1395/BL'

all_region_median <- subset(all_region_median,Quality != 'Q35')
CNV_Q40_RG.gain = get_dotbar(all_region_median,"region",'gain',"Gain Region (Mbp)",c(0, 5));CNV_Q40_RG.gain
CNV_Q40_RG.loss = get_dotbar(all_region_median,"region",'loss',"Loss Region (Mbp)",c(0, 5));CNV_Q40_RG.loss
CNV_Q40_RG.loh = get_dotbar(all_region_median,"region",'loh',"Loss Region (Mbp)",c(0, 5));CNV_Q40_RG.loh


get_boxcompare_cnv = function(tnseq.df.sub,d_type,m_type, y_lim,ct_1395.sum=data.frame()){
  if (m_type=='all'){
    tnseq.df.sub = tnseq.df.sub
    tnseq.df.sub$cnv_type = 'all'
  } else {
    tnseq.df.sub <- subset(tnseq.df.sub,cnv_type == m_type)
  }
  tnseq.df.sub[,'jaccard'] = tnseq.df.sub[,d_type]
  #tnseq.df.test <- tnseq.df.sub %>%
  #  group_by(cnv_type, Depth) %>%
  #  rstatix::t_test(
  #    jaccard ~ Quality,p.adjust.method = 'none',comparisons =list(c("Q30","Q40"))) %>%
  #  select(-df, -statistic, -p) # Remove details
  tnseq.df.test <- tnseq.df.sub %>%
    group_by(cnv_type, Depth) %>%
    filter(n_distinct(jaccard) > 1, n_distinct(Quality) > 1) %>% 
    t_test(jaccard ~ Quality, p.adjust.method = 'none', comparisons = c("Q30", "Q40")) %>%
    mutate(p.adj.signif = p.adjust(p, method = "none")) %>%  # Adjust p-value and calculate significance
    add_significance() %>%  # Adds significance symbols based on adjusted p-values
    select(-df, -statistic, -p)
  tnseq.df.test
  #tnseq.df.test <- tnseq.df.test %>% add_xy_position(x = "Depth")
  #if (d_type=='jaccard'){
  #tnseq.df.test$y.position <- ifelse(tnseq.df.test$group2=='Q40',tnseq.df.test$y.position-0.1,tnseq.df.test$y.position-0.2)
  #} 
  #tnseq.df.sub$group <- paste0(tnseq.df.sub$Quality,'_',tnseq.df.sub$Depth)
  tnseq.df.test <- tnseq.df.test %>% add_xy_position(x = "Depth")
  Max_F1_score <- tnseq.df.sub %>%
    group_by(cnv_type, Depth) %>%
    dplyr::summarise(y.position = max(jaccard, na.rm = TRUE))
  tnseq.df.test<- merge(tnseq.df.test,Max_F1_score,by=c('cnv_type','Depth'))
  tnseq.df.test$y.position <- tnseq.df.test$y.position.y
  if (length(ct_1395.sum) == 0){
    box_p <- ggplot(tnseq.df.sub, aes(x = Depth, y = jaccard)) +
      geom_rect(xmin = 0.5, xmax = 1.5,
                ymin = 0, ymax = 11, fill ="#E9E9E9") +
      geom_rect(xmin = 2.5, xmax = 3.5,
                ymin = 0, ymax = 11, fill ="#E9E9E9")+
      # Overlay narrow violin plot
      geom_violin(
        aes(fill = Quality),
        width = 0.6,              # Narrow violin width
        scale = "width",          # Uniform width scaling
        position = position_dodge(width = 0.8),  # Align with boxplot
        alpha = 0.35,
        size = 0.5,
        trim = TRUE
      )+
      geom_boxplot(
        aes(fill = Quality),
        #   fill = 'white',
        width = 0.1,              # Wider boxplot
        position = position_dodge(width = 0.8), 
        outlier.shape = NA,       # Hide outliers
        alpha = 0.9,
        size = 0.5
      ) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 3, shape=21,colour='white',aes(fill = Quality))+
      stat_summary(aes(color=Quality),fun = median, geom = "point",size = 1.5,position = position_dodge(width = 1)) +  # Add median
      #stat_summary(fun = median, geom = "text", color = "black", 
      #             aes(label = round(..y.., 2)), vjust = -1)+
      #geom_line(aes(color = Quality,group=Quality),position = position_dodge(width = 0.5),alpha=0.5,size=1)+
      scale_color_manual(values = colors)+
      scale_fill_manual(values = colors)+
      #facet_grid(type~.)+
      ylim(y_lim)+ylab(paste0(d_type,' ',m_type))+xlab('')+
      my_theme+ stat_pvalue_manual(
        tnseq.df.test, label = "p.signif", 
        step.increase = 0.01,size = 5
      )
  } else {
    ct_1395_plot <- ggplot(subset(ct_1395.sum,cnv_type==m_type), aes(x = Depth, y = region_mean,fill=Quality)) +
      geom_rect(xmin = 0.5, xmax = 1.5,
                ymin = 0, ymax = 200000000, fill ="#E9E9E9") +
      geom_rect(xmin = 2.5, xmax = 3.5,
                ymin = 0, ymax = 200000000, fill ="#E9E9E9")+
      geom_bar(stat = "identity", position = "dodge", width = 0.7)+
      geom_errorbar(aes(ymin = region_mean - region_sd, ymax = region_mean + region_sd), position = position_dodge(0.7), width = 0.2)+
      scale_fill_manual(values = colors)+xlab('')+
      my_theme+theme(legend.position = "none")+ylab('Region (Mbp)')
    box_p <- ggplot(tnseq.df.sub, aes(x = Depth, y = jaccard)) +
      geom_rect(xmin = 0.5, xmax = 1.5,
                ymin = 0, ymax = 11, fill ="#E9E9E9") +
      geom_rect(xmin = 2.5, xmax = 3.5,
                ymin = 0, ymax = 11, fill ="#E9E9E9")+
      # Overlay narrow violin plot
      geom_violin(
        aes(fill = Quality),
        width = 0.6,              # Narrow violin width
        scale = "width",          # Uniform width scaling
        position = position_dodge(width = 0.8),  # Align with boxplot
        alpha = 0.35,
        size = 0.5,
        trim = TRUE
      )+
      geom_boxplot(
        aes(fill = Quality),
        #fill = 'white',
        width = 0.1,              # Wider boxplot
        position = position_dodge(width = 0.8), 
        outlier.shape = NA,       # Hide outliers
        alpha = 0.9,
        size = 0.5
      ) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 3, shape=21,colour='white',aes(fill = Quality))+
      stat_summary(aes(color=Quality),fun = median, geom = "point",size = 1.5,position = position_dodge(width = 1)) +  # Add median
      #stat_summary(fun = median, geom = "text", color = "black", 
      #             aes(label = round(..y.., 2)), vjust = -1)+
      #geom_line(aes(color = Quality,group=Quality),position = position_dodge(width = 0.5),alpha=0.5,size=1)+
      scale_color_manual(values = colors)+
      scale_fill_manual(values = colors)+
      #facet_grid(type~.)+
      ylim(y_lim)+ylab(paste0('Jaccard ',m_type))+xlab('')+
      my_theme+ stat_pvalue_manual(
        tnseq.df.test, label = "p.signif", 
        step.increase = 0.01,size = 5
      )
    #+annotation_custom(grob=ggplotGrob(ct_1395_plot),ymin = y_lim[1], ymax=y_lim[1]+0.7, xmin=2.5, xmax=4.5)
  }
  return(box_p)
}

#CNV_Q40_Fig1 = get_boxcompare_cnv(all_region_median,"region","CNV Region (bp)",c(0.5, 1),scale_y_continuous=TRUE);CNV_Q40_Fig1

# CNV Jaccard (Reproducibility) ------------------

# ELE & ILM
all_jaccard <- data.frame(Depth=1,jaccard=1,cnv_type=1,Quality=1,source=1)

for (i in dir('ELE_ILM/Somatic_CNV',pattern = '*bed.jaccard.txt')){
  cnv_type = str_split(i,'\\.',simplify = T)[,3]
  q_type = str_split(i,'\\.',simplify = T)[,1]
  reg_df <- read.table(paste0('ELE_ILM/Somatic_CNV/',i),header = T)
  reg_df$bed1 <- gsub('WES_','',reg_df$bed1) %>% gsub('_10G','',.)
  reg_df$bed2 <- gsub('WES_','',reg_df$bed2) %>% gsub('_10G','',.)
  reg_df$Depth1 = str_split(reg_df$bed1,'\\_',simplify = T)[,4]
  reg_df$Depth2 = str_split(reg_df$bed2,'\\_',simplify = T)[,4]
  reg_df$source <- paste0(reg_df$bed1,'.vs.',reg_df$bed2)
  reg_df$Quality <- 'Q40'
  if (grepl('ILM',q_type)){
    reg_df$Quality='Q30'
  } else if (grepl('Mix',q_type)){
    reg_df$Quality = 'Q35'
  }
  
  reg_df$Depth <- ifelse(reg_df$Depth1==reg_df$Depth2,reg_df$Depth1,'')
  reg_df$cnv_type = cnv_type
  reg_df = subset(reg_df,Depth!="",c('jaccard','Quality','cnv_type','Depth','source'))
  all_jaccard = rbind(all_jaccard,reg_df)
}

all_jaccard = all_jaccard[-1,]
all_jaccard$Depth <- factor(all_jaccard$Depth,levels=c("30X", "60X", "90X", "120X"))
all_jaccard <- subset(all_jaccard,jaccard!=1)
all_jaccard <- subset(all_jaccard,Quality!='Q35')
CNV_Q40_JI.gain = get_boxcompare_cnv(all_jaccard,'jaccard',"gain",c(0.2, 1.1),all_region_median);CNV_Q40_JI.gain
CNV_Q40_JI.loss = get_boxcompare_cnv(all_jaccard,'jaccard',"loss",c(0, 1.1),all_region_median);CNV_Q40_JI.loss
CNV_Q40_JI.loh = get_boxcompare_cnv(all_jaccard,'jaccard',"loh",c(0, 1.2),all_region_median);CNV_Q40_JI.loh
all_jaccard <- subset(all_jaccard,Quality != 'Q35')
CNV_Q40_JI.all = get_boxcompare_cnv(all_jaccard,'jaccard',"all",c(0, 1.0));CNV_Q40_JI.all

all_jaccard$F1.score <- all_jaccard$jaccard
Fig.S3a.CNV_JI.gain.test <- get_test_type_with_effectsize(subset(all_jaccard,cnv_type=='gain'))
Fig.S3c.CNV_JI.loss.test <- get_test_type_with_effectsize(subset(all_jaccard,cnv_type=='loss'))
Fig.S3e.CNV_JI.loh.test <- get_test_type_with_effectsize(subset(all_jaccard,cnv_type=='loh'))
Fig.4g.CNV_JI.all.test <- get_test_type_with_effectsize(all_jaccard)
# At 30× sequencing depth, the paired data structure was incomplete due to failed CNV detection in some Q30 samples. 
# Therefore, comparisons between Q30 and Q40 platforms were performed using the Wilcoxon rank-sum test for independent samples.

all_jaccard.30X <- subset(all_jaccard,Depth=='30X')
wilcox.test(jaccard ~ Quality, data = all_jaccard.30X)
calculate_unpaired_cohens_d(subset(all_jaccard.30X,Quality=='Q30')$jaccard,subset(all_jaccard.30X,Quality=='Q40')$jaccard)

wilcox.test(jaccard ~ Quality, data = subset(all_jaccard.30X,cnv_type=='gain'))
calculate_unpaired_cohens_d(subset(all_jaccard.30X,Quality=='Q30' & cnv_type=='gain')$jaccard,subset(all_jaccard.30X,Quality=='Q40' & cnv_type=='gain')$jaccard)

wilcox.test(jaccard ~ Quality, data = subset(all_jaccard.30X,cnv_type=='loss'))
calculate_unpaired_cohens_d(subset(all_jaccard.30X,Quality=='Q30' & cnv_type=='loss')$jaccard,subset(all_jaccard.30X,Quality=='Q40'& cnv_type=='loss')$jaccard)

wilcox.test(jaccard ~ Quality, data = subset(all_jaccard.30X,cnv_type=='loh'))
calculate_unpaired_cohens_d(subset(all_jaccard.30X,Quality=='Q30' & cnv_type=='loh')$jaccard,subset(all_jaccard.30X,Quality=='Q40' & cnv_type=='loh')$jaccard)


all_jaccard_median <- all_jaccard %>%
  dplyr::group_by(cnv_type, Quality, Depth) %>%
  dplyr::summarize(
    across(
      jaccard,
      list(
        mean = ~mean(.),
        sd = ~sd(.)
      )
    )
  )

all_jaccard_median$Depth <- factor(gsub('X','',all_jaccard_median$Depth),levels=c(30, 60, 90, 120))
all_jaccard_median$Sample <- 'HCC1395/BL'


all_jaccard %>% 
  group_by(Quality,Depth)  %>%
  dplyr::summarise(jaccard_mean = median(jaccard,na.rm=T))

#CNV_Q40_Fig2 = get_dotbar(all_jaccard_median,"jaccard","Jaccard Index",c(0.5, 1));CNV_Q40_Fig2

# CNV copy number correlation (Reproducibility) ------------------
# ELE & ILM
all_pcc <- data.frame(Sample=1,cnPCC=1,cnType=1,source=1)

for (i in dir('ELE_ILM/Somatic_CNV',pattern = '*.cnPCC.*')){
  cn_pcc <- read.csv(paste0('ELE_ILM/Somatic_CNV/',i),header =TRUE)
  print(head(cn_pcc))
  cn_pcc$source <- cn_pcc$Sample
  all_pcc <- rbind(all_pcc,cn_pcc)
}

all_pcc <- all_pcc[-1,]
all_pcc$Sample <- gsub('_10G','',all_pcc$Sample)
all_pcc$Quality <- lapply(all_pcc$Sample,function(x){
  p='Q40'
  if (grepl('ILM',x)){
    p='Q30'
  } else if (grepl('Mix',x)){
    p = 'Q35'
  }
  return(p)
}) %>% unlist()
all_pcc$Depth <- str_split(all_pcc$Sample,"[_.]",simplify = T)[,2]


#all_pcc = all_pcc[-1,]
all_pcc$cnv_type <- all_pcc$cnType


all_pcc$Sample <- 'HCC1395/BL'

all_pcc_median <- all_pcc %>%
  dplyr::group_by(cnv_type, Quality, Depth,Sample) %>%
  dplyr::summarize(
    across(
      cnPCC,
      list(
        mean = ~mean(.),
        sd = ~sd(.)
      )
    )
  )

all_pcc_median$Depth <- factor(gsub('X','',all_pcc_median$Depth),levels=c(30, 60, 90, 120))
all_pcc <- subset(all_pcc,Quality !='Q35')
#CNV_Q40_Fig3 = get_dotbar(all_pcc_median,"cnPCC","PCC of Copy Number",c(0.5, 1));CNV_Q40_Fig3
all_pcc$cnv_type <- all_pcc$cnType
all_pcc$Depth <- factor(all_pcc$Depth,levels=c("30X", "60X", "90X", "120X"))
CNV_Q40_PCC.gain_loss = get_boxcompare_cnv(all_pcc,'cnPCC',"gain/loss",c(0.6, 1.05));CNV_Q40_PCC.gain_loss
CNV_Q40_PCC.loh = get_boxcompare_cnv(all_pcc,'cnPCC',"loh",c(0.5, 1.05));CNV_Q40_PCC.loh
CNV_Q40_PCC.all =  get_boxcompare_cnv(all_pcc,'cnPCC',"all",c(0.5, 1.05));CNV_Q40_PCC.all
all_pcc %>% 
  group_by(Quality,Depth)  %>%
  dplyr::summarise(jaccard_mean = median(cnPCC,na.rm=T))


all_pcc$F1.score <- all_pcc$cnPCC
Fig.S3g.CNV_PCC.gain_loss.test <- get_test_type_with_effectsize(subset(all_pcc,cnv_type=='gain/loss'))
Fig.S3h.CNV_PCC.loh.test <- get_test_type_with_effectsize(subset(all_pcc,cnv_type=='loh'))
Fig.4h.CNV_PCC.all.test <- get_test_type_with_effectsize(all_pcc)

all_pcc.30X <- subset(all_pcc,Depth=='30X')
wilcox.test(cnPCC ~ Quality, data = all_pcc.30X)
calculate_unpaired_cohens_d(subset(all_pcc.30X,Quality=='Q30')$cnPCC,subset(all_pcc.30X,Quality=='Q40')$cnPCC)

wilcox.test(cnPCC ~ Quality, data = subset(all_pcc.30X,cnv_type=='gain'))
calculate_unpaired_cohens_d(subset(all_pcc.30X,Quality=='Q30' & cnv_type=='gain')$cnPCC,subset(all_pcc.30X,Quality=='Q40' & cnv_type=='gain')$cnPCC)

wilcox.test(cnPCC ~ Quality, data = subset(all_pcc.30X,cnv_type=='loss'))
calculate_unpaired_cohens_d(subset(all_pcc.30X,Quality=='Q30' & cnv_type=='loss')$cnPCC,subset(all_pcc.30X,Quality=='Q40'& cnv_type=='loss')$cnPCC)

#Save statistical results --------------
library(openxlsx)
test_file <- ls(pattern = '*Fig*') %>% sort()
#test_file <- test_file[-grep("get_test_type", test_file)]
wb <- createWorkbook()

# Loop through each variable
for (var_name in test_file) {
  print(var_name)
  # Get the variable object (dataframe)
  df <- get(var_name)
  
  # Use variable name as sheet name (note: Excel sheet names have length limits and cannot contain special characters)
  # Can truncate or clean variable name to meet Excel requirements
  sheet_name <- substr(var_name, 1, 31)  # Excel sheet name max 31 characters
  
  # Replace possible invalid characters
  sheet_name <- gsub("[\\[\\]\\*\\?\\/\\\\]", "_", sheet_name)
  
  # Add worksheet
  addWorksheet(wb, sheet_name)
  
  # Write data to worksheet
  writeData(wb, sheet = sheet_name, df)
  
  cat("Added worksheet:", sheet_name, "\n")
}

# Save Excel file
saveWorkbook(wb, "WES_Statistic_test_type.xlsx", overwrite = TRUE)

#Save Figures----------
CNV_Q40_all <- (((CNV_Q40_JI.gain+CNV_Q40_RG.gain+plot_layout(widths = c(1,0.5)))/
                   (CNV_Q40_JI.loss+CNV_Q40_RG.loss+plot_layout(widths = c(1,0.5)))/
                   (CNV_Q40_JI.loh+CNV_Q40_RG.loh+plot_layout(widths = c(1,0.5))))+(CNV_Q40_PCC.gain_loss|CNV_Q40_PCC.loh))+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect',heights = c(1,1,1,0.8)) & theme(legend.position='bottom',
                                                               plot.tag = element_text(color = "black",size = 18,face = "bold"))
CNV_Q40_all

ggsave('Fig.CNV.pdf',CNV_Q40_all,width=8.15, height=5.7*2,dpi = 300)

# Germline main figure: nist_snv,nist_indel,quartet_snv,quartet_indel,quartet_snv.mcr,quartet_indel.mcr,nist_snv.rmsk,nist_indel.rmsk,quartet_snv.rmsk,quartet_indel.rmsk

germline_main <- ((quartet_snv.mcr|quartet_indel.mcr)/((quartet_snv|quartet_indel))/((quartet_snv.rmsk|quartet_indel.rmsk)))+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect',heights = c(1,1,1)) & theme(legend.position='',
                                                           plot.tag = element_text(color = "black",size = 16,face = "bold"))
germline_main
ggsave('Fig.germline-final-revised.pdf',germline_main,width=8.15, height=5.7*1.5,dpi = 300)




# Somatic main figure: hcc_snv,hcc_indel,hcc_snv.rmsk,hcc_indel.rmsk,hcc_snv.lowVAF,hcc_indel.lowVAF
somatic_main <- (hcc_snv|hcc_indel)/((hcc_snv.rmsk+hcc_indel.rmsk+hcc_snv.lowVAF+
                                        hcc_indel.lowVAF)+plot_layout(guides='collect',widths = c(2,1,2,2)))/
  (CNV_Q40_JI.all|CNV_Q40_PCC.all)+
  plot_annotation(tag_levels = 'a') & theme(legend.position='',
                                            plot.tag = element_text(color = "black",size = 16,face = "bold"))
somatic_main
ggsave('Fig.somatic-final.pdf',somatic_main,width=8.15, height=5.7*1.5,dpi = 300)




sup_wes.germline <-(nist_snv|nist_indel)/((quartet_snv.rmsk.mcr|quartet_indel.rmsk.mcr)/(nist_snv.rmsk|nist_indel.rmsk))/(rep.p.germline)+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 16,face = "bold"))
ggsave('Fig.wes_sup.germline-final.pdf',sup_wes.germline,width=8.15, height=5.7*1.5,dpi = 300)

sup_wes.germline_v2 <- (df_nist.sub.count.snv|df_nist.sub.count.indel)/(quartet_snv.rmsk.count.snv|quartet_snv.rmsk.count.indel)/
  (df_quartet.mcr.rsmk_count.snv|df_quartet.mcr.rsmk_count.indel)+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 16,face = "bold"))

ggsave('Fig.wes_sup.germline-rmsk.pdf',sup_wes.germline_v2,width=8.15*1.5, height=5.7*2.5,dpi = 300)


sup_wes.somatic <- rep.p.somatic/(cor.p.somatic|rep_df.p1)/(df_hcc1395.sub.rmsk_count.snv+df_hcc1395.sub.rmsk_count.indel+plot_layout(widths = c(2,1)))+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='',
                                        plot.tag = element_text(color = "black",size = 16,face = "bold"))

ggsave('Fig.wes_sup.somatic-final-revised.pdf',sup_wes.somatic,width=8.15, height=5.7*1.5,dpi = 300)




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
