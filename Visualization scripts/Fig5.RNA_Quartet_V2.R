#Quartet RNAseq 分析------
setwd("/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/Data_performance_V2/RNAseq/Quartet")
source("/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/quartet-rseqc-report/exp2qcdt/R/exp2qcdt.R")
source("/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/quartet-rseqc-report/exp2qcdt/R/config.R")
source("/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/quartet-rseqc-report/exp2qcdt/R/multiple_group_output.R")
source("/Users/duanshumeng/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/quartet-rseqc-report/exp2qcdt/R/one_group_output.R")
source('~/PGx_lab/HCC1395/ElementVSIllumina/Scripts/Visualization scripts/Statistics_analysis.R')
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggthemes)
library(edgeR)
library(cowplot)
library(scales)
library(RColorBrewer)
library(grDevices)
library(grid)
library(dtplyr)
library(plyr)
library(ggpubr)
library(gghalves)
load("Quartet.rnaseq.RData")
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
#colors = c('#4472CA','#E69F00','#CC79A7')
#names(colors) <- c('Q40_E','Q30','Q40_M')
#colors = c('#3783BB','#91569F','#C94741')
#names(colors) <- c('Q30','Q35','Q40')
colors=c("#709AE1FF","#F05C3BFF")
names(colors) <- c('Q30','Q40')
#SNR-----
output_snr_res <- function(exp_table_file, count_file, meta_file){
  dt_fpkm <- fread(exp_table_file)
  dt_fpkm_log <- data.table(apply(dt_fpkm[, !'gene_id'], 2, function(x)(log2(x + 0.01))))
  dt_fpkm_log[, gene_id := dt_fpkm$gene_id]
  
  dt_counts <- fread(count_file)
  change_cols <- colnames(dt_counts[, !'gene_id'])
  dt_counts[, (change_cols):= lapply(.SD, as.numeric), .SDcols = change_cols]
  dt_meta <- fread(meta_file)
  dt_detect_gene <- do.call(cbind, lapply(unique(dt_meta$sample), function(x){
    detect_res <- apply(dt_counts[, dt_meta[sample == x][['library']], with = F], 1, function(x){length(which(x >= 3)) >= 2})
    return(
      detect_res = detect_res
    ) 
  }))
  gene_list_snr <- dt_counts[['gene_id']][apply(dt_detect_gene, 1, function(x){any(x)})]
  exp_design = (dt_meta[, .(library, group = sample)] %>% setkey(., library))
  dt_fpkm_f <- dt_fpkm_log[gene_list_snr, on = .(gene_id)]
  dt_fpkm_zscore <- data.table(t(apply(dt_fpkm_f[, dt_meta$library, with = F], 1, function(x){(x - mean(x))/sd(x)})))
  #print(dt_fpkm_zscore)
  dt_fpkm_zscore <- na.omit(dt_fpkm_zscore)
  pca_list <- get_pca_list(dt_fpkm_zscore, exp_design, dt_meta)
  return(pca_list)
}


get_SNR_plot <- function(exp_table_file, count_file, meta_file,type){
  dt_snr <- output_snr_res(exp_table_file, count_file, meta_file)
  snr_gene_num <- dt_snr$gene_num[1]
  dt_snr$SNR %>% unique() %>% print()
  ## figure of pca with snr
  pt_snr <- ggplot(dt_snr, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 3, show.legend = FALSE) +
    theme_few() +
    guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 1, title.position = "top")) +
    scale_fill_manual(values = c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745","#f5f3a5","#e1c5ee")) +
    scale_color_manual(values = c("#2f5c85", "#7BC8A4", "#FFC65D", "#F16745","#f5f3a5","#e1c5ee")) +
    theme(plot.title = element_text(hjust = 0.5,size = 14,face = "bold")) + 
    #coord_fixed()+
    labs(
      title = paste(type,"\n"," SNR: ", dt_snr$SNR[1], ' (N = ', dt_snr$gene_num[1], ')', sep = ""),
      x = paste("PC1 (", dt_snr$PC1_ratio, "%)", sep = ""),
      y = paste("PC1 (", dt_snr$PC2_ratio, "%)", sep = ""))+theme(plot.title = element_text(size = 10))
  return(list(unique(dt_snr$SNR),pt_snr,snr_gene_num))
}
sizes=c('2G','4G','6G','8G','10G')
snr_lst <- list()
snr_plots <- list()
snr_gene_lst <- list()
for (i in sizes){
  print(i)
  for (a in c('ILM','ELE')){
    print(paste0('quartet_fpkm_',a,'.',i,'.csv'))
    snr <- get_SNR_plot(paste0('quartet_fpkm_',a,'.',i,'.csv'),paste0('quartet_count_',a,'.',i,'.csv'),paste0('metadata_',a,'.',i,'.csv'),paste0(a,' ',i))
    snr_lst[[paste0(a,' ',i)]] = snr[[1]]
    snr_plots[[paste0(a,' ',i)]] = snr[[2]]
    snr_gene_lst[[paste0(a,' ',i)]] = snr[[3]]
  }
}
snr_df.quartet <- snr_lst %>% unlist() %>% as.data.frame()
library(gridExtra)
Quartet_snr_detail <- grid.arrange(grobs = snr_plots, ncol = 4)
#ggsave('Quartet.SNR_addlegend.pdf',Quartet_snr_detail,width=5.7*2, height=5.7,dpi = 300)
ggsave('Quartet.SNR.pdf',Quartet_snr_detail,width=5.7*1.5, height=5.7*1.2,dpi = 300)
#Gene count-------
snr_gene_df.quartet <- snr_gene_lst %>% unlist() %>% as.data.frame()

#PCC & RMSE----
get_RMSE <- function(predictions,observed){
  differences <- predictions - observed
  # Step 2
  squared_differences <- differences^2
  # Step 3
  mean_squared_difference <- mean(squared_differences)
  # Step 4
  RMSE <- sqrt(mean_squared_difference)
  return(RMSE)
}
pcc_lst = list()
rmse_df = data.frame(gene=1,
                     RMSE=1,
                     source=1)
for (i in sizes){
  print(i)
  for (a in c('ILM','ELE','Mix')){
    print(paste0(a,'_',i,'/performance_assessment/logfc_cor_ref_test.txt'))
    ele_expr <- read.table(paste0(a,'_',i,'/performance_assessment/logfc_cor_ref_test.txt'),header = TRUE)
    correlation = cor.test(ele_expr$meanlogFC_test, ele_expr$meanlogFC_ref,method = "pearson")$estimate
    ele_rmse <- ddply(ele_expr, .(gene), summarise, RMSE = get_RMSE(meanlogFC_test,meanlogFC_ref))
    ele_rmse$source <- paste0(a,' ',i)
    #rmse = get_RMSE(ele_expr$meanlogFC_test, ele_expr$meanlogFC_ref)
    pcc_lst[[paste0(a,' ',i)]] = correlation
    #rmse_lst[[paste0(a,' ',i)]] = rmse
    rmse_df <- rbind(rmse_df,ele_rmse)
    #print(correlation)
  }
}
pcc_df.quartet <- pcc_lst %>% unlist() %>% as.data.frame()
rmse_df.quartet <- rmse_df

#CV----
get_cv <- function(input_file,quality){
  df = fread(input_file)
  calculate_cv <- function(row) {
    row[row == 0] <- NA
    #print(row)
    mean_row <- mean(row,na.rm = T)
    sd_row <- sd(row,na.rm = T)
    cv <- (sd_row / mean_row) * 100
    return(cv)
  }
  #print(head(df))
  #df = ilm.count_low[1:10,]
  for (i in c('D5','D6','F7','M8')){
    col_names <- colnames(df)[!grepl('CV',colnames(df))]
    aim_cols <- col_names[grepl(i,col_names)]
    cv = paste0(i)
    #print(i)
    #print(aim_cols)
    df[,cv]<- apply(df[,..aim_cols,with = FALSE],1,calculate_cv)
  }
  df_sub<-df[,c('D5','D6','F7','M8')] %>% melt()
  df_sub$Quality = quality
  return(df_sub)
}
df_cv = data.frame(variable=1,
                   value=1,
                   Quality=1)
for (i in sizes){
  print(i)
  for (a in c('ILM','ELE','Mix')){
      ele.cv <- get_cv(paste0('quartet_fpkm_',a,'.',i,'.csv'),paste0(a,' ',i))
      df_cv=rbind(df_cv,ele.cv)
  }
}
df_cv.quartet <- df_cv

#CV of quantified genes-------
df_cv_count = data.frame(variable=1,
                   value=1,
                   Quality=1,
                   gene_count=1)
gene_lst = list()
get_cv_count <- function(df,quality){
  calculate_cv <- function(row) {
    row[row == 0] <- NA
    #print(row)
    mean_row <- mean(row,na.rm = T)
    sd_row <- sd(row,na.rm = T)
    cv <- (sd_row / mean_row) * 100
    return(cv)
  }
  #print(head(df))
  #df = ilm.count_low[1:10,]
  for (i in c('D5','D6','F7','M8')){
    col_names <- colnames(df)[!grepl('CV',colnames(df))]
    aim_cols <- col_names[grepl(i,col_names)]
    cv = paste0(i)
    #print(i)
    #print(aim_cols)
    df[,cv]<- apply(df[,..aim_cols,with = FALSE],1,calculate_cv)
  }
  df_sub<-df[,c('D5','D6','F7','M8')] %>% melt()
  df_sub$Quality = quality
  return(df_sub)
}
for (i in sizes){
  print(i)
  for (a in c('ILM','ELE','Mix')){
    print(paste0('quartet_fpkm_',a,'.',i,'.csv'))
    dt_fpkm <- fread(paste0('quartet_fpkm_',a,'.',i,'.csv'))
    dt_counts <- fread(paste0('quartet_count_',a,'.',i,'.csv'))
    change_cols <- colnames(dt_counts[, !'gene_id'])
    dt_counts[, (change_cols):= lapply(.SD, as.numeric), .SDcols = change_cols]
    dt_meta <- fread(paste0('metadata_',a,'.',i,'.csv'))
    dt_detect_gene <- do.call(cbind, lapply(unique(dt_meta$sample), function(x){
      detect_res <- apply(dt_counts[, dt_meta[sample == x][['library']], with = F], 1, function(x){length(which(x >= 3)) >= 2})
      return(
        detect_res = detect_res
      ) 
    }))
    gene_list_snr <- dt_counts[['gene_id']][apply(dt_detect_gene, 1, function(x){any(x)})]
    gene_lst[[paste0(a,' ',i)]] = gene_list_snr
    dt_fpkm_10G <- fread(paste0('quartet_fpkm_',a,'.','10G','.csv'))
    dt_fpkm_f <- dt_fpkm_10G[gene_list_snr, on = .(gene_id)]
    ele.cv <- get_cv_count(dt_fpkm_f,paste0(a,' ',i))
    ele.cv$gene_count <- length(gene_list_snr)
    df_cv_count=rbind(df_cv_count,ele.cv)
  }
}


samples=c('D5','D6','F7','M8','T1','T2')
fpkm_cv_df <- data.frame(group=1,
                         cv_value=1,
                         Quality=1,
                         Size=1)
row_cv <- function(x) {
  apply(x, 1, function(y) sd(y) / mean(y) * 100)
}
for (a in c('ILM','ELE')){
  for (b in sizes){
    list1 <- gene_lst[[paste0(a," ",b)]] %>% unlist()
    gene_lists.grouped <- gene_lists %>% dplyr::group_by(Gene) %>% dplyr::summarise(List = paste(List, collapse = "-"), .groups = 'drop')
    dt_fpkm <- fread(paste0('quartet_fpkm_',a,'.',b,'.csv'))
    dt_fpkm <- dt_fpkm[gene_id %in% list1]
    for (i in samples){
      dt_fpkm[, (i) := rowMeans(.SD), .SDcols = patterns(i)]
    }
    breaks <- c(0, 2, 4, 6, 8, 10, Inf)
    labels <- c("[0,2)", "[2,4)", "[4,6)", "[6,8)", "[8,10)", "≥10")
    
    aim_cols=colnames(dt_fpkm)[grepl('F7',colnames(dt_fpkm))]

    dt_fpkm.F7 <- dt_fpkm[,.SD,.SDcols = aim_cols]
    dt_fpkm.F7[, group := cut(F7, breaks = breaks, labels = labels, right = FALSE, include.lowest = TRUE)
    ][, cv_value := row_cv(.SD), .SDcols = 1:3]
    dt_fpkm.F7[, Quality := ifelse(a == "ILM", "Q30", "Q40")]
    dt_fpkm.F7[,Size := b]
    fpkm_cv_df = rbind(fpkm_cv_df,dt_fpkm.F7[,c('group','cv_value','Quality','Size')])
  }
}
fpkm_cv_df[,-1]

fpkm_cv_df[fpkm_cv_df$Size=='10G',] %>% 
  dplyr::group_by(group,Quality) %>%
  dplyr::summarise(mean_cv = mean(cv_value,na.rm = T))

fpkm_cv_compare <- ggplot(data = fpkm_cv_df[fpkm_cv_df$Size=='10G',] ,       
       aes(x=Quality, y=cv_value, fill=Quality)) +
  #geom_half_violin(side="r")+
  #geom_half_boxplot(side = "r") +    
  # geom_half_point_panel(side = "l")+
  geom_half_violin(side = "r", color=NA, alpha=0.35) +    
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, 
                    outlier.size = 1,
                    outlier.stroke = 0.2,
                    outlier.shape = NA, 
                    width=0.2, linewidth=0.5) +    
  geom_half_point_panel(aes(fill = Quality),side = "l", 
                        range_scale = .85,
                        shape=21, size=1.5, color="white")+
  scale_fill_manual(values = c(pal_simpsons()(16)[c(2,8)])) +
  #scale_fill_simpsons()+
  theme_bw()+
  theme(
    axis.text.x = element_text(face = "bold",color = "black"),
    axis.text.y = element_text(face = "bold",color = "black"),
    axis.title  = element_text(face = "bold",color = "black"),
    panel.grid.minor  =  element_blank(),
    legend.background =  element_blank(),
    legend.title = element_text(face = "bold",color = "black"),
    legend.text = element_text(color = "black")
  )+
  scale_y_continuous(expand = expansion(mult=c(0.1,0.1)))+
  stat_compare_means(aes(group = Quality),
                     label="p.signif",
                     method="wilcox.test",
                     show.legend = F)+
  facet_wrap(.~group,scales = "free",ncol = 6)

# remove NA
df_clean <- fpkm_cv_df[fpkm_cv_df$Size=='10G',] %>% filter(!is.na(cv_value))
Fig.5g.fpkm_cv_compare.test <- get_test_type_with_effectsize_fpkm(df_clean)
#Save test results---------------
library(openxlsx)
test_file <- ls(pattern = '*Fig*') %>% sort()
#test_file <- test_file[-grep("get_test_type", test_file)]
wb <- createWorkbook()

for (var_name in test_file) {
  print(var_name)

  df <- get(var_name)

  sheet_name <- substr(var_name, 1, 31)  

  sheet_name <- gsub("[\\[\\]\\*\\?\\/\\\\]", "_", sheet_name)

  addWorksheet(wb, sheet_name)

  writeData(wb, sheet = sheet_name, df)
  
  cat("Add sheet:", sheet_name, "\n")
}


saveWorkbook(wb, "RNAseq_Statistic_test_type.xlsx", overwrite = TRUE)

library(data.table)

get_fpkm_plot <- function(dt_fpkm.F7,gene_lists.grouped,pname){
  gene_lists.grouped$Group <- str_split(gene_lists.grouped$List,'\\-',simplify = T)[,1]
  dt_fpkm <- unique(dt_fpkm)
  colnames(gene_lists.grouped) <- c('gene_id','Group_raw','Group')

  dt_fpkm.log2 <- as.data.frame(log2(dt_fpkm[,-1]+0.001))
  rownames(dt_fpkm.log2) <- dt_fpkm$gene_id
  dt_fpkm.log2$gene_id <- rownames(dt_fpkm.log2)
  dt_fpkm.log2 <- merge(dt_fpkm.log2,gene_lists.grouped[,c('gene_id','Group')],by='gene_id',all.x = T)
  dt_fpkm.log2$Mean_log2FPKM <- rowMeans(dt_fpkm.log2[,-c(1,20)])
  dt_fpkm.log2<- na.omit(dt_fpkm.log2)
  fpkm.p <- ggplot(dt_fpkm.log2, aes(x=Group, y=Mean_log2FPKM,fill=Group)) +
    geom_boxplot(alpha=0.5,width=0.45,
                 position=position_dodge(width=0.8),
                 size=0.75,outlier.colour = NA)+
    my_theme+
    xlim(c("2G", "4G",'6G','8G','10G'))+
    theme(legend.position = 'none',)+xlab(paste0(pname,' Data volumes'))
  return(fpkm.p)
}

df_cv_count.quartet <- df_cv_count
library(tidyr)
library(UpSetR)
upset_lst = list()
fpkmP_lst = list()
fpkm_lst = list()
for (a in c('ILM','ELE','Mix')){
  list1 <- gene_lst[[paste0(a," 2G")]] %>% unlist()
  list2 <- gene_lst[[paste0(a," 4G")]] %>% unlist()
  list3 <- gene_lst[[paste0(a," 6G")]] %>% unlist()
  list4 <- gene_lst[[paste0(a," 8G")]] %>% unlist()
  list5 <- gene_lst[[paste0(a," 10G")]] %>% unlist()
  gene_lists <- data.frame(
    Gene = c(list1,list2,list3,list4,list5),
    List = factor(rep(sizes, times = c(length(list1), length(list2), length(list3), length(list4), length(list5))))
  ) %>% as.data.frame()
  gene_lists.grouped <- gene_lists %>% dplyr::group_by(Gene) %>% dplyr::summarise(List = paste(List, collapse = "-"), .groups = 'drop')
  dt_fpkm <- fread(paste0('quartet_fpkm_',a,'.','10G','.csv'))
  fpkm_lst[[a]] <- 
  fpkm_p <- get_fpkm_plot(dt_fpkm,gene_lists.grouped,a)
  fpkmP_lst[[a]] <- fpkm_p
  str(gene_lists)
  gene_matrix <- gene_lists %>%
    dplyr::count(Gene, List) %>%
    pivot_wider(names_from = List, values_from = n, values_fill = list(n = 0))
  
  row.names(gene_matrix) <- gene_matrix$Gene
  gene_matrix <- gene_matrix %>% dplyr::select(-Gene)
  colnames(gene_matrix) <- factor(colnames(gene_matrix),levels = c('2G','4G','6G','8G','10G'))
  up_p <- upset(as.data.frame(gene_matrix), order.by = c("freq"), decreasing = c(TRUE),
                mb.ratio = c(0.55,0.45),
                point.size = 3,
                line.size = 1.2,
                number.angles = 15,
                sets.x.label = paste0("Counts ",a), 
                sets.bar.color=brewer.pal(5,"Set1"),
                shade.color="red",
                text.scale=c(1.2,1.2,1.2,1.2,1.2,1.2))
  upset_lst[[a]] <- up_p
}
library(gridExtra)
Quartet_fpkm_detail <- grid.arrange(grobs = fpkmP_lst, ncol = 1)
ggsave('Quartet.FPKM.pdf',Quartet_fpkm_detail,width=5.7, height=5.7*1.2,dpi = 300)

pdf("Quartet.upset.ILM.pdf",width = 5.7, height = 5.7*0.6,pointsize = 16)
#png("Quartet.upset.ILM.png",width = 570*3, height = 570*1.5,res=300)
upset_lst[[1]]
dev.off()
pdf("Quartet.upset.ELE.pdf",width = 5.7, height = 5.7*0.6,pointsize = 16)
#png("Quartet.upset.ELE.png",width = 570*3, height = 570*1.5,res=300)
upset_lst[[2]]
dev.off()
pdf("Quartet.upset.Mix.pdf",width = 5.7, height = 5.7*0.6,pointsize = 16)
#png("Quartet.upset.MGI.png",width = 570*3.5, height = 570*1.5,res=300)
upset_lst[[3]]
dev.off()

#MCC----
quartet_ref <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/quartet-rseqc-report/exp2qcdt/data/ref_data_fc_value.csv',
                        header = TRUE)
quartet_ref$Feature <- quartet_ref$gene
quartet_ref$logFC <- quartet_ref$meanlogFC


get_DEGtype <- function(df){
  DEGtype <- apply(df,1,function(x){
    #print(x['adj.P.Val']%>% as.numeric())
    #print(x['logFC'] %>% as.numeric())
    if (x['P.Value'] %>% as.numeric() >= 0.05){
      return('non-DEG')
    }else{
      if (x['logFC'] %>% as.numeric() > 1){
        return('up-regulate')
      } else if (x['logFC'] %>% as.numeric()  < -1){
        return('down-regulate')
      } else {
        return('non-DEG')
      }
    }
  }) %>% unlist()
}

get_mcc <- function(tn,tp,fn,fp){
  print(paste(tn,tp,fn,fp,sep = ' '))
  
  sqrt_term = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc = (tp*tn - fp*fn)/sqrt_term
  return(mcc)
}

if (TURE){get_f1 <- function(tn,tp,fn,fp){
  print(paste(tn,tp,fn,fp,sep = ' '))
  
  precision = tp/(tp+fp)
  recall = tp/(tp+fn)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  return(f1_score)
}
}
get_mcc_all <- function(test_expr,ref_expr,smp){
  TN_non_non = intersect(intersect(subset(test_expr,DEGtype=='non-DEG')$gene_id,ref_expr$gene),
                         subset(ref_expr,DEGtype=='non-DEG')$gene) %>% length() 
  FN_non_up = intersect(intersect(subset(test_expr,DEGtype=='non-DEG')$gene_id,ref_expr$gene),
                        subset(ref_expr,DEGtype=='up-regulate')$gene) %>% length() 
  
  FN_non_down = intersect(intersect(subset(test_expr,DEGtype=='non-DEG')$gene_id,ref_expr$gene),
                          subset(ref_expr,DEGtype=='down-regulate')$gene) %>% length()
  
  FP_up_non = intersect(intersect(subset(test_expr,DEGtype=='up-regulate')$gene_id,ref_expr$gene),
                        subset(ref_expr,DEGtype=='non-DEG')$gene) %>% length() 
  TP_up_up = intersect(intersect(subset(test_expr,DEGtype=='up-regulate')$gene_id,ref_expr$gene),
                       subset(ref_expr,DEGtype=='up-regulate')$gene) %>% length()
  FP_up_down = intersect(intersect(subset(test_expr,DEGtype=='up-regulate')$gene_id,ref_expr$gene),
                         subset(ref_expr,DEGtype=='down-regulate')$gene) %>% length()
  
  FP_down_non = intersect(intersect(subset(test_expr,DEGtype=='down-regulate')$gene_id,ref_expr$gene),
                          subset(ref_expr,DEGtype=='non-DEG')$gene) %>% length() 
  
  FP_down_up = intersect(intersect(subset(test_expr,DEGtype=='down-regulate')$gene_id,ref_expr$gene),
                         subset(ref_expr,DEGtype=='up-regulate')$gene) %>% length() 
  
  TP_down_down = intersect(intersect(subset(test_expr,DEGtype=='down-regulate')$gene_id,ref_expr$gene),
                           subset(ref_expr,DEGtype=='down-regulate')$gene) %>% length() 
  
  #FN_not_up = intersect(setdiff(test_expr$gene_id,ref_expr$gene),
  #                      subset(ref_expr,DEGtype=='up-regulate')$gene) %>% length() 
  
  #FN_not_down = intersect(setdiff(test_expr$gene_id,ref_expr$gene),
  #                       subset(ref_expr,DEGtype=='down-regulate')$gene) %>% length() 
  #grep("FN", ls(), value = TRUE)
  FN = FN_non_down + FN_non_up 
  #+ FN_not_down+FN_not_up
  #grep("TN", ls(), value = TRUE)
  TN = TN_non_non
  #grep('FP',ls(),value = TRUE)
  FP = FP_down_non + FP_down_up + FP_up_down+FP_up_non
  #grep('TP',ls(),value = TRUE)
  TP = TP_down_down+TP_up_up
  df = data.frame(TN=TN,
                  TP=TP,
                  FN=FN,
                  FP=FP,
                  Sample=smp)
  #mcc <- get_mcc(TN*0.01,TP*0.01,FN*0.01,FP*0.01)
  return(df)
  
}
mcc_df = data.frame(TN=1,
                   TP=1,
                   FN=1,
                   FP=1,
                   Sample=1,
                   Compare=1)


for (i in sizes){
  print(i)
  for (a in c('ILM','ELE','Mix')){
    ele_expr <- read.table(paste0(a,'_',i,'/performance_assessment/logfc_test.txt'),header = TRUE)
    ele_expr$DEGtype <- get_DEGtype(ele_expr)
    for (b in quartet_ref$compare %>% unique()){
      quartet_ref.sub <- subset(quartet_ref,compare==b)
      ele_expr.sub <- subset(ele_expr,compare==b)
      f1 = get_mcc_all(ele_expr.sub,quartet_ref.sub,paste0(a,' ',i))
      f1$Compare <- b
      mcc_df <- rbind(mcc_df,f1)
    }

  }
}
MCC.df <- mcc_df[-1,]
MCC.df$MCC.score <- apply(MCC.df,1,function(x){
  #print(x)
  get_mcc(x['TN'] %>% as.numeric(),
         x['TP']%>% as.numeric(),
         x['FN']%>% as.numeric(),
         x['FP']%>% as.numeric())
}) %>% unlist()
MCC.df.quartet <- MCC.df
#load("Quartet.rnaseq.RData")

#The ratio distribution of Quartet and MAQC reference datasets
quartet_ref <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/quartet-rseqc-report/exp2qcdt/data/ref_data_fc_value.csv',
                        header = TRUE)
quartet_ref$Feature <- quartet_ref$gene
quartet_ref$logFC <- quartet_ref$meanlogFC

quartet_ref$logFC %>% abs() %>% max()

maqc_ref <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/quartet-rseqc-report/exp2qcdt/data/ref_data_fc_value.csv',
                        header = TRUE)
quartet_ref$Feature <- quartet_ref$gene
quartet_ref$logFC <- quartet_ref$meanlogFC

maqc_ref <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/MAQC2_sampleA_sampleB/MAQC_taq_reference.csv')

maqc_ref$DEGtype <- gsub("Non-DEG","non-DEG",maqc_ref$DEG)
maqc_ref$gene <- maqc_ref$Gene.Name
maqc_ref$logFC <- maqc_ref$log2FC
maqc_ref$Feature <- maqc_ref$Gene.Name

maqc_ref$logFC %>% abs() %>% max(na.rm = T)

all_ref <- data.frame(logFC_abs=c(quartet_ref$logFC,maqc_ref$logFC) %>% abs(),
                      DEGtype=c(quartet_ref$DEGtype,maqc_ref$DEGtype),
                      Sample=c(rep('Quartet',nrow(quartet_ref)),rep('MAQC',nrow(maqc_ref))))


logFC_dist <- ggplot(all_ref, aes(x = logFC_abs,fill=Sample)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.5, 
                 position = "identity", alpha = 0.5,color='gray')+
  my_theme+xlim(c(0,10))+scale_fill_manual(values = c('#FF4500' ,'#8A2BE2'))

## |logFC|< 1 and other
data_summary <- all_ref %>%
  dplyr::mutate(Status = ifelse(logFC_abs < 1, "|logFC| < 1", "|logFC| >= 1")) %>%
  dplyr::group_by(Status,Sample) %>%
  dplyr::summarise(Count = n(), .groups = 'drop')


data_summary <- na.omit(data_summary)
logFC_dist.2 <- ggplot(data_summary, aes(x = Sample, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill",alpha=0.7,width = 0.5,color='gray') +
  labs( y = "Proportion",
       fill = "") +
  scale_y_continuous(labels = scales::percent) +
  my_theme

logFC_all <- logFC_dist + logFC_dist.2

ggsave('logFC_abs_dist.png',logFC_all,width=8.15, height=5.7,dpi = 300)

## 0<|logFC|<1
quartet_ref_lowdiff <- subset(quartet_ref,abs(quartet_ref$logFC)<1)

subset(quartet_ref,abs(quartet_ref$logFC)<1) %>% dim()

quartet_ref %>% dim()

# Identify genes with higher expression in M8 than in D6 and genes with lower expression in M8 than in D6
quartet_ref.M8_D6 <- subset(quartet_ref, compare == 'M8/D6' & DEGtype != 'non-DEG')

# Supplement the 1:3 and 3:1 fitting situations for M8 and D3 ----
T1_df = data.frame(value = 1, Quality = 1, Size = 1)
T2_df = data.frame(value = 1, Quality = 1, Size = 1)
for (i in sizes) {
  print(i)
  for (a in c('ILM', 'ELE', 'Mix')) {
    quality = 'Q40'
    if (a == 'ILM') {
      quality = 'Q30'
    } else if (a == 'Mix') {
      quality = 'Q35'
    }
    print(paste0('quartet_fpkm_', a, '.', i, '.csv'))
    df_fpkm <- read.csv(paste0('quartet_fpkm_', a, '.', i, '.csv'))
    df_fpkm$gene <- df_fpkm$gene_id
    df_fpkm_diff <- merge(df_fpkm, quartet_ref.M8_D6['gene'], by = 'gene', all = F)
    df_fpkm.sub <- df_fpkm_diff[, grepl('M8|D6|T1|T2', colnames(df_fpkm_diff))]
    for (b in c('M8', 'D6', 'T1', 'T2')) {
      df_fpkm.sub[b] <- rowMeans(df_fpkm.sub[, grepl(b, colnames(df_fpkm.sub))])
    }
    
    # M8 > D6: M8 > T1 > T2 > D6
    # D6 > M8: D6 > T2 > T1 > M8
    # Theoretical value of T1-D6/M8-D6 is 0.75
    df_fpkm.sub <- subset(df_fpkm.sub, (df_fpkm.sub$M8 > df_fpkm.sub$D6 & df_fpkm.sub$T1 > df_fpkm.sub$T2 & df_fpkm.sub$T2 > df_fpkm.sub$D6) |
                            (df_fpkm.sub$D6 > df_fpkm.sub$T2 & df_fpkm.sub$T2 > df_fpkm.sub$T1 & df_fpkm.sub$T1 > df_fpkm.sub$M8))[c('M8', 'D6', 'T1', 'T2')]
    df_fpkm.sub.T1 <- df_fpkm.sub
    df_fpkm.sub.T1['value'] <- (log2(df_fpkm.sub.T1['T1']) - log2(df_fpkm.sub.T1['D6'])) / (log2(df_fpkm.sub.T1['M8']) - log2(df_fpkm.sub.T1['D6']))
    replace(df_fpkm.sub.T1[, 'value'], is.infinite(df_fpkm.sub.T1[, 'value']), 1) %>% median(., na.rm = T)
    df_fpkm.sub.T1$Quality = quality
    df_fpkm.sub.T1$Size = i
    T1_df <- rbind(T1_df, df_fpkm.sub.T1[, c('value', 'Quality', 'Size')])
    
    # Theoretical value of T2-D6/M8-D6 is 0.25
    # D6 > M8, T2 < D6, M8 < T2
    df_fpkm.sub.T2 <- df_fpkm.sub
    df_fpkm.sub.T2['value'] <- (log2(df_fpkm.sub.T2['T2']) - log2(df_fpkm.sub.T2['D6'])) / (log2(df_fpkm.sub.T2['M8']) - log2(df_fpkm.sub.T2['D6']))
    replace(df_fpkm.sub.T2[, 'value'], is.infinite(df_fpkm.sub.T2[, 'value']), 1) %>% median(., na.rm = T)
    df_fpkm.sub.T2$Quality = quality
    df_fpkm.sub.T2$Size = i
    T2_df <- rbind(T2_df, df_fpkm.sub.T2[, c('value', 'Quality', 'Size')])
  }
}
T1_df <- T1_df[-1, ]
T1_df$RMSE <- abs(T1_df$value - 0.75)
T1_df$Size <- gsub('G', '', T1_df$Size)
T1_df$Size <- factor(T1_df$Size, levels = c('2', '4', '6', '8', '10'))
T1_df <- subset(T1_df, Quality != 'Q35')

T1_df.test <- T1_df %>%
  dplyr::group_by(Size) %>%
  dplyr::filter(n_distinct(value) > 1, n_distinct(Quality) > 1) %>% 
  rstatix::wilcox_test(value ~ Quality, p.adjust.method = 'none', comparisons = c("Q30", "Q40")) %>%
  rstatix::mutate(p.adj.signif = ifelse(is.na(p), NA, p.adjust(p, method = "none"))) %>% 
  rstatix::add_significance() %>%  
  select(-statistic, -p)
T1_df.test
T1_df.test <- T1_df.test %>% rstatix::add_xy_position(x = "Size")
T1_df <- subset(T1_df, Quality != 'Q35')

get_rmse_plot <- function(T1_df, trd_type, hline = 0.75) {
  T1_df.sum <- T1_df %>%
    dplyr::group_by(Size, Quality) %>%
    dplyr::summarize(
      count = n(),
      across(
        value,
        list(
          median = ~median(., na.rm = T),
          sd = ~sd(., na.rm = T)
        )
      )
    )
  plot1 <- ggplot(T1_df.sum, aes(x = Size, y = count, fill = Quality)) +
    geom_rect(xmin = 0.5, xmax = 1.5,
              ymin = 0, ymax = 11, fill = "#E9E9E9") +
    geom_rect(xmin = 2.5, xmax = 3.5,
              ymin = 0, ymax = 11, fill = "#E9E9E9") +
    geom_rect(xmin = 4.5, xmax = 5.5,
              ymin = 0, ymax = 11, fill = "#E9E9E9") +
    geom_bar(stat = "identity", position = "dodge", width = 0.7, color = 'black') +
    scale_fill_manual(values = colors) +
    my_theme + theme(legend.position = "none")
  plot1
  
  T1_plot <- ggplot(T1_df, aes(x = Size, y = value)) +
    geom_rect(xmin = 0.5, xmax = 1.5,
              ymin = 0, ymax = 11, fill = "#E9E9E9") +
    geom_rect(xmin = 2.5, xmax = 3.5,
              ymin = 0, ymax = 11, fill = "#E9E9E9") +
    geom_rect(xmin = 4.5, xmax = 5.5,
              ymin = 0, ymax = 11, fill = "#E9E9E9") +
    geom_violin(
      aes(fill = Quality),
      width = 0.6,
      alpha = 0.5,
      scale = "width",
      position = position_dodge(width = 0.8),
      size = 0.5,
      trim = TRUE
    ) +
    geom_boxplot(
      aes(fill = Quality),
      width = 0.1,
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      size = 0.5
    ) +
    geom_hline(yintercept = hline, linetype = "solid", color = "gray") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    ylab(paste0(trd_type)) +
    ylim(c(0, 1)) +
    my_theme + annotation_custom(grob = ggplotGrob(plot1),
                                 ymin = 0.1, ymax = 0.5, xmin = 2.5, xmax = 5.5)
  return(T1_plot)
}

T1_p <- get_rmse_plot(T1_df, '(T1-D6)/(M8-D6)'); T1_p

T2_df <- T2_df[-1, ]
T2_df$RMSE <- abs(T2_df$value - 0.25)
T2_df$Size <- gsub('G', '', T2_df$Size)
T2_df$Size <- factor(T2_df$Size, levels = c('2', '4', '6', '8', '10'))
T2_df <- subset(T2_df, Quality != 'Q35')
T2_p <- get_rmse_plot(T2_df, '(T2-D6)/(M8-D6)', hline = 0.25); T2_p

T2_df.test <- T2_df %>%
  dplyr::group_by(Size) %>%
  dplyr::filter(n_distinct(value) > 1, n_distinct(Quality) > 1) %>% 
  rstatix::wilcox_test(value ~ Quality, p.adjust.method = 'none', comparisons = c("Q30", "Q40")) %>%
  rstatix::mutate(p.adj.signif = ifelse(is.na(p), NA, p.adjust(p, method = "none"))) %>% 
  rstatix::add_significance() %>%  
  select(-statistic, -p)
T2_df.test

T1_df$Type = '(T1-D6)/(M8-D6)'
T2_df$Type = '(T2-D6)/(M8-D6)'
T_df <- rbind(T1_df, T2_df)

get_rmse_plot_v2 <- function(T_df, hline = 0) {
  p <- ggplot(data = T_df,
              aes(x = Quality, y = value, fill = Quality)) +
    geom_half_violin(side = "r", color = NA, alpha = 0.35) +
    geom_half_boxplot(side = "r", errorbar.draw = FALSE,
                      outlier.size = 1,
                      outlier.stroke = 0.2,
                      outlier.shape = NA,
                      width = 0.2, linewidth = 0.5) +
    geom_half_point_panel(aes(fill = Quality), side = "l",
                          range_scale = .85,
                          shape = 21, size = 1.5, color = "white") +
    scale_fill_manual(values = c(pal_simpsons()(16)[c(2, 8)])) +
    theme_bw() +
    theme(
      axis.text.x = element_text(face = "bold", color = "black"),
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      panel.grid.minor = element_blank(),
      legend.background = element_blank(),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    stat_compare_means(aes(group = Quality),
                       label = "p.signif",
                       method = "t.test",
                       show.legend = F) +
    facet_wrap(Type ~ Size, scales = "free", ncol = 5)
  if (hline != 0) {
    p = p + geom_hline(yintercept = hline, linetype = "solid", color = "gray") + ylim(c(0, 1))
  }
  return(p)
}

T1_tr_value <- get_rmse_plot_v2(subset(T_df, Type == "(T1-D6)/(M8-D6)"), hline = 0.75); T1_tr_value
T2_tr_value <- get_rmse_plot_v2(subset(T_df, Type == "(T2-D6)/(M8-D6)"), hline = 0.25); T2_tr_value

ls_all <- ls()
save(df_cv_count.quartet, df_cv.quartet, MCC.df.quartet, pcc_df.quartet, rmse_df.quartet, snr_df.quartet, snr_gene_df.quartet,
     T1_df, T2_df,
     file = "Quartet.rnaseq.RData")

load('Quartet.rnaseq.RData')
T_all <- (T1_tr_value / T2_tr_value) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect') & theme(legend.position = '',
                                          plot.tag = element_text(color = "black", size = 16, face = "bold"))
ggsave('Fig.rna-tr.pdf', T_all, width = 8.15, height = 5.7, dpi = 300)
ggsave('Fig.rna-fpkm.pdf', fpkm_cv_compare + theme(legend.position = ''), width = 8.15, height = 5.7 * 0.4, dpi = 300)
