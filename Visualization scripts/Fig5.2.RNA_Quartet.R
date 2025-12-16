setwd("/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq")
source("./quartet-rseqc-report/exp2qcdt/R/exp2qcdt.R")
source("./quartet-rseqc-report/exp2qcdt/R/config.R")
source("./quartet-rseqc-report/exp2qcdt/R/multiple_group_output.R")
source("./quartet-rseqc-report/exp2qcdt/R/one_group_output.R")
source('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/Quartet_Multiomics_Ratio_Code/utils/theme_nature.r')
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

colors = c('#4472CA','#E69F00')
names(colors) <- c('Q40','Q30')
#Differentially Expression-------------
ele_expr <- read.table('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/ELE/performance_assessment/logfc_cor_ref_test.txt',
                       header = TRUE)

ilm_expr <- read.table('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/ILM/performance_assessment/logfc_cor_ref_test.txt',
                       header = TRUE)

breaks <- c(0, 0.99, 1.99,2.99,3.99, 4.99, 5.99,10)
# Create labels for the intervals
labels <- c("[0,1)", "[1,2)", "[2,3)", "[3,4)","[4,5)","[5,6)", ">=6")
ele_expr$Interval <- cut(abs(ele_expr$meanlogFC_ref), breaks = breaks, labels = labels, include.lowest = TRUE)
ilm_expr$Interval <- cut(abs(ilm_expr$meanlogFC_ref), breaks = breaks, labels = labels, include.lowest = TRUE)

#Gene Counts----
ele_expr$group <- paste0(ele_expr$Interval,'.',ele_expr$compare)
ele_count <- table(ele_expr$group) %>% as.data.frame()
ele_count$Quality <- 'Q40'
ele_count$Interval <- str_split(ele_count$Var1,'\\.',simplify = T)[,1]
ele_count$Compare <- str_split(ele_count$Var1,'\\.',simplify = T)[,2]

ilm_expr$group <- paste0(ilm_expr$Interval,'.',ilm_expr$compare)
ilm_count <- table(ilm_expr$group) %>% as.data.frame()
ilm_count$Quality <- 'Q30'
ilm_count$Interval <- str_split(ilm_count$Var1,'\\.',simplify = T)[,1]
ilm_count$Compare <- str_split(ilm_count$Var1,'\\.',simplify = T)[,2]

quartet_count.all <- rbind(ele_count,ilm_count)
quartet_count.all<- quartet_count.all[,c('Interval','Freq','Quality','Compare')]
quartet_count.all$Sample <- 'Quartet'

gene_counts = ggplot(quartet_count.all,aes(x=Interval,y=Freq,fill=Quality,color=Quality))+geom_bar(stat = "identity", 
                                                                             position = "dodge",color="black",alpha = 0.5)+facet_grid(Compare~.)+
  scale_color_manual(values = colors)+
  stat_summary(geom = "text", aes(label = paste(round(..y.., digits = 0))),
               position = position_dodge(width = 0.7), vjust = -0.7,hjust=0,size = 3)+ylim(0,max(quartet_count.all$Freq)+2000)+
  scale_fill_manual(values = colors)+theme_nature_border()+labs(x='Ratio interval ',y='Gene Counts')+
  theme(legend.position = 'none')


#PCC of fold change----
library(plyr)
fc_q40 <- ggplot(ele_expr, aes(x = meanlogFC_test, y = meanlogFC_ref)) + theme_classic()+
  geom_point(aes(color=Interval),size=3,alpha=0.75)+ scale_color_brewer(palette = "Dark2")+scale_fill_brewer(palette = "Dark2")+
  geom_smooth(aes(color = Interval, fill = Interval), 
              method = "lm", fullrange = TRUE) +stat_cor(method = "pearson",size = 3)+
  theme_nature_border()+facet_grid(rows=vars(Interval),cols = vars(compare))+theme(legend.position="none")+
  labs(y = "logFC (Reference)", x = "LogFC (Q40)")


fc_q30 <-ggplot(ilm_expr, aes(x = meanlogFC_test, y = meanlogFC_ref)) + theme_classic()+
  geom_point(aes(color=Interval),size=3,alpha=0.75)+ scale_color_brewer(palette = "Dark2")+scale_fill_brewer(palette = "Dark2")+
  geom_smooth(aes(color = Interval, fill = Interval), 
              method = "lm", fullrange = TRUE) +stat_cor(method = "pearson",size = 3)+
  theme_nature_border()+facet_grid(rows=vars(Interval),cols = vars(compare))+theme(legend.position="none")+
  labs(y = "logFC (Reference)", x = "LogFC (Q30)")
p.fc <- ggarrange(fc_q40, fc_q30,ncol = 2)
quartet_FC.cor <- p.fc

cor.ele <- ddply(ele_expr, .(group), summarise, correlation = cor.test(meanlogFC_test, meanlogFC_ref,method = "pearson")$estimate)

check_group_size <- function(group_data) {
  if (length(group_data) < 4) {  # 假设每个组中至少需要30个观测值
    return(FALSE)
  } else {
    return(TRUE)
  }
}

cor.ilm <- ddply(ilm_expr, .(group), function(sub_data) {
  if (check_group_size(sub_data$group)) {
    correlation <- cor.test(sub_data$meanlogFC_test, sub_data$meanlogFC_ref, method = "pearson")$estimate
  } else {
    correlation <- NA  
  }
  return(data.frame(correlation = correlation))
})

cor.ele$Quality <- 'Q40'
cor.ilm$Quality <- 'Q30'

quartet_cor.df <- rbind(cor.ele,cor.ilm)
quartet_cor.df$Interval <- str_split(quartet_cor.df$group,'\\.',simplify = T)[,1]
quartet_cor.df$Compare <- str_split(quartet_cor.df$group,'\\.',simplify = T)[,2]

quartet_cor.df<- subset(quartet_cor.df,quartet_cor.df$correlation > 0)
quartet_cor.df$Sample <- 'Quartet'

p_line.ratio <- function(test,d_type,comp=F){
  if (d_type != 'Reproducibility'){
    test$Reproducibility = test[,d_type]
  } 
  p1 = ggplot(test,aes(x=Interval,y=Reproducibility))+theme_classic()+
    geom_point(size = 4,alpha=0.7,aes(fill = Quality,color =Quality))+scale_fill_manual(values = colors)+
    geom_line(aes(color =Quality,group=Quality))+scale_color_manual(values = colors)+
    theme_nature_border()+ggtitle(paste0(''))+
    labs(x="Ratio Interval", y=d_type) 
  if (comp){
    p1=p1+facet_grid(Compare~.)
  }
  return(p1)
}

cor.ratio <- p_line.ratio(quartet_cor.df,'correlation',TRUE)

#RMSE---
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
rmse.ele <- ddply(ele_expr, .(Interval), summarise, RMSE = get_RMSE(meanlogFC_ref,meanlogFC_test))
rmse.ele$Quality <- 'Q40'
rmse.ilm <- ddply(ilm_expr, .(Interval), summarise, RMSE = get_RMSE(meanlogFC_ref,meanlogFC_test))
rmse.ilm$Quality <- 'Q30'

quartet_rmse <- rbind(rmse.ele,rmse.ilm)
quartet_rmse$Sample <- 'Quartet'
rmse.ratio <- p_line.ratio(quartet_rmse,'RMSE')

#SNR----
output_snr_res <- function(exp_table_file, count_file, meta_file,gene_list_snr){
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
  
  exp_design = (dt_meta[, .(library, group = sample)] %>% setkey(., library))
  dt_fpkm_f <- dt_fpkm_log[gene_list_snr, on = .(gene_id)]
  dt_fpkm_zscore <- data.table(t(apply(dt_fpkm_f[, dt_meta$library, with = F], 1, function(x){(x - mean(x))/sd(x)})))
  #print(dt_fpkm_zscore)
  dt_fpkm_zscore <- na.omit(dt_fpkm_zscore)
  pca_list <- get_pca_list(dt_fpkm_zscore, exp_design, dt_meta)
  return(pca_list)
}


get_SNR_plot <- function(exp_table_file, count_file, meta_file,gene_list_snr,type){
  dt_snr <- output_snr_res(exp_table_file, count_file, meta_file,gene_list_snr)
  snr_gene_num <- dt_snr$gene_num[1]
  dt_snr$SNR %>% unique() %>% print()
  ## figure of pca with snr
  pt_snr <- ggplot(dt_snr, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 3, show.legend = FALSE) +
    theme_few() +
    guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 1, title.position = "top")) +
    scale_fill_manual(values = c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745")) +
    scale_color_manual(values = c("#2f5c85", "#7BC8A4", "#FFC65D", "#F16745")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(
      title = paste(type,"\n"," SNR: ", dt_snr$SNR[1], ' (N = ', dt_snr$gene_num[1], ')', sep = ""),
      x = paste("PC1 (", dt_snr$PC1_ratio, "%)", sep = ""),
      y = paste("PC1 (", dt_snr$PC2_ratio, "%)", sep = ""))+theme(plot.title = element_text(size = 10))
  return(list(unique(dt_snr$SNR),pt_snr))
}

snr_lst <- list()
snr_plots <- list()
for (i in labels){
  print(i)
  aim_genes.ele <- subset(ele_expr,Interval==i)$gene
  aim_genes.ilm <- subset(ilm_expr,Interval==i)$gene
  snr.ele = get_SNR_plot("quartet_fpkm_ele_D5M8.csv", "quartet_count_ele_D5M8.csv", "metadata_ele.csv",aim_genes.ele,paste0('Q40 ',i))
  snr.ilm = get_SNR_plot("quartet_fpkm_ilm_D5M8.csv", "quartet_count_ilm_D5M8.csv", "metadata_ilm.csv",aim_genes.ilm,paste0('Q30 ',i))
  snr_lst[[paste0('Q40 ',i)]] = snr.ele[[1]]
  snr_lst[[paste0('Q30 ',i)]] = snr.ilm[[1]]
  snr_plots[[paste0('Q30 ',i)]] = snr.ilm[[2]]
  snr_plots[[paste0('Q40 ',i)]] = snr.ele[[2]]
}

library(gridExtra)
snr_detail <- grid.arrange(grobs = snr_plots, ncol = 4)
quartet_snr_detail <- snr_detail

quartet_snr_df <- snr_lst %>% unlist() %>% as.data.frame()
quartet_snr_df$SNR <- quartet_snr_df$.
quartet_snr_df$Quality <- str_split(rownames(quartet_snr_df),' ',simplify = T)[,1]
quartet_snr_df$Interval <- str_split(rownames(quartet_snr_df),' ',simplify = T)[,2]
snr.ratio <- p_line.ratio(quartet_snr_df,'SNR')
snr.ratio$Sample <- 'Quartet'
cor.ratio + snr.ratio

save(quartet_cor.df,quartet_count.all,quartet_rmse,quartet_snr_df,quartet_snr_detail, quartet_FC.cor,
     file = "Quartet.rnaseq.RData")

#CV----
get_cv <- function(df_1,quality,exp_type,ratio_df){
  gene_lst = subset(ratio_df,Interval==i)$gene
  df = subset(df_1,gene_id %in% gene_lst)
  calculate_cv <- function(row) {
    row[row == 0] <- NA
    print(row)
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
  df_sub$Exp_type = exp_type
  return(df_sub)
}
ele.count = fread("quartet_fpkm_ele_D5M8.csv")
ilm.count = fread("quartet_fpkm_ilm_D5M8.csv")

df_cv = data.frame(variable=1,
                   value=1,
                   Quality=1,
                   Exp_type=1)

for (b in labels){
    ele.cv <- get_cv(ele.count,'Q40',b,ele_expr)
    ilm.cv <- get_cv(ilm.count,'Q30',b,ilm_expr)
    df_cv=rbind(df_cv,ele.cv)
    df_cv=rbind(df_cv,ilm.cv)
}

df_cv = df_cv[-1,]
df_cv$Compare <- df_cv$variable
df_cv$FPKM.CV <- df_cv$value
df_cv$Interval <- df_cv$Exp_type
fpkm.cv <- ggplot(df_cv,aes(x=Interval,y=FPKM.CV,fill=Quality,color=Quality))+geom_split_violin(trim=F,alpha = 0.5,width = 1,color="black")+
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175))+scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+theme_nature_border()+labs(x='Ratio Interval',y='FPKM.CV')

fig6 <- ((gene_counts+ggtitle('a')+theme(plot.title = element_text(hjust =0, vjust=0,size = 15, face = "bold")))|
    (cor.ratio+ggtitle('b')+theme(plot.title = element_text(hjust =0, vjust=0,size = 15, face = "bold"))))+
  ((fpkm.cv+ggtitle('c')+theme(plot.title = element_text(hjust =0, vjust=0,size = 15, face = "bold")))/
     (snr.ratio+ggtitle('d')+theme(plot.title = element_text(hjust =0, vjust=0,size = 15, face = "bold"))))+plot_layout(guides = 'collect')

ggsave('Fig6.Ratio.png',fig6,width=8.15*2, height=5.7*2,dpi = 300)




#M8/D6 Titration difference-------
quartet_color <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745","#f5f3a5","#e1c5ee")
names(quartet_color) <- c('D5','D6','F7','M8','T1','T2')
#The coincidence of Titration differenc (TCR)

get_FC <- function(dt_counts,dt_meta){
  dt_fc_test <- do.call(rbind, lapply(list(c('M8', 'D6'), c('T1', 'D6'), c('T2', 'D6')), function(x){
    compare_name <- paste(x[1], '/', x[2], sep = '')
    
    ### at least tow replicate counts >= 3 
    dt_detect_gene <-  data.table(apply(dt_counts[, dt_meta[x[1], on = .(sample)][['library']], with = F], 1, function(x){length(which(x >= 3)) >= 2}),
                                  apply(dt_counts[, dt_meta[x[2], on = .(sample)][['library']], with = F], 1, function(x){length(which(x >= 3)) >= 2}))
    gene_list_com <- dt_counts[['GeneName']][apply(dt_detect_gene, 1, function(x){all(x)})]
    
    dt_count_compare <- dt_counts[, c('GeneName',
                                      dt_meta[c(x[1], x[2]), on = .(sample)][['library']]), with = F]
    dt_count_compare_com <- dt_count_compare[gene_list_com, on = .(GeneName)]
    group_compare <- dt_meta[c(x[1], x[2]), on = .(sample)][['sample']]
    group_compare[group_compare == x[1]] <- 'ZZZ'
    group_compare[group_compare == x[2]] <- 'AAA'
    
    # DEG analysis was perform with edger pakcages
    deg_output <- DEGanalysis(dt_count_compare_com[, !'GeneName'], group_compare)
    deg_output[, GeneName := dt_count_compare_com[as.numeric(deg_output$gene), 'GeneName']]
    deg_output[, compare := compare_name]
    deg_output_d <- deg_output[,c('GeneName', 'compare', 'logFC'), with = FALSE]
    colnames(deg_output_d) <- c('gene', 'compare', 'meanlogFC')
    return(deg_output_d)
  }))
  return(dt_fc_test)
}

ele_diff <- get_FC(fpkm_df.ele %>% data.table::data.table(),
       dt_meta.ele)

ilm_diff <- get_FC(fpkm_df.ilm %>% data.table::data.table(),
                   dt_meta.ilm)
model_func <- function(x, k1, k2) {
  log2(k1 + k2 * 2^x)
}

get_k <- function(z,t_type){
  k1 = z/(z+3) 
  k2 = 3*z/(3*z+1)
  if (t_type=='T1'){
    return(c(k1=k1,k2=k2))
  } else {
    return(c(k1=k2,k2=k1))
  }
  
}
initial.t1 <-get_k(0.965,'T1')
initial.t2 <-get_k(0.941,'T2')
x_range <- seq(-6, max(6), length.out = 10000)
theoretical.t1 <- data.frame(x = x_range, y = model_func(x_range, initial.t1[[1]], initial.t1[[2]]))
theoretical.t2 <- data.frame(x = x_range, y = model_func(x_range, initial.t2[[1]], initial.t2[[2]]))
#print(theoretical_data)

#return k1 and k2
get_k <- function(z,t_type){
  k1 = z/(z+3) 
  k2 = 3*z/(3*z+1)
  if (t_type=='T1'){
    return(c(k1=k1,k2=k2))
  } else {
    return(c(k1=k2,k2=k1))
  }
  
}


get_fit <- function(data,initial){
  fit <- nls(y ~ model_func(x, k1, k2), data = data, start = initial)
  k1=coef(fit)[1]
  k2=coef(fit)[2]
  return(list(fit,k1,k2))
}

get_fit.plot <- function(data){

  initial.t1 <-get_k(0.965,'T1')
  initial.t2 <-get_k(0.941,'T2')
  x_range <- seq(-6, max(6), length.out = 100)
  theoretical.t1 <- data.frame(x = x_range, y = model_func(x_range, initial.t1[[1]], initial.t1[[2]]))
  theoretical.t2 <- data.frame(x = x_range, y = model_func(x_range, initial.t2[[1]], initial.t2[[2]]))
  
  print(initial.t1)
  #Q40
  df.q40 <- data[,colnames(data)[grepl('Q40',colnames(data))]]
  colnames(df.q40) <- colnames(data)[grepl('Q40',colnames(data))] %>% gsub('.Q40','',.) %>% gsub('.t1','',.)
  t1.q40 <- get_fit(df.q40,initial.t1)
  fit.t1.q40 <- t1.q40[[1]]
  k1.t1.q40=coef(fit.t1.q40)[1]
  k2.t1.q40=coef(fit.t1.q40)[2]
  #Q30
  df.q30 <- data[,colnames(data)[grepl('Q30',colnames(data))]]
  colnames(df.q30) <- colnames(data)[grepl('Q30',colnames(data))] %>% gsub('.Q30','',.) %>% gsub('.t1','',.)
  t1.q30 <- get_fit(df.q30,initial.t1)
  fit.t1.q30 <- t1.q30[[1]]
  k1.t1.q30=coef(fit.t1.q30)[1]
  k2.t1.q30=coef(fit.t1.q30)[2]
  

  p.t1 <- ggplot() +
    geom_point(data=data, aes(x = x.Q40, y = y.t1.Q40),color='#4472CA',alpha=0.5) +  
    geom_line(data= theoretical.t1,aes(x = x, y = y), color = '#9FE61F', linetype = "dashed")+
    geom_smooth(data=data, aes(x = x.Q40, y = y.t1.Q40),method = "nls",method.args=list(formula =y ~ log2(k1 + k2 * 2^x),
                                                start=list(k1=k1.t1.q40,k2=k2.t1.q40)),se = FALSE, color = '#4472CA', linetype = "dashed")+
    geom_point(data=data, aes(x = x.Q30, y = y.t1.Q30),color='#E69F00',alpha=0.5) + 
    geom_smooth(data=data, aes(x = x.Q30, y = y.t1.Q30),method = "nls",method.args=list(formula =y ~ log2(k1 + k2 * 2^x),
                                                start=list(k1=k1.t1.q30,k2=k2.t1.q30)),
                se = FALSE, color = '#E69F00', linetype = "dashed")+
    theme_nature_border()+coord_fixed()+
    labs(x = "M8.vs.D6 (Log2FC)", y = "T1.vs.D6 (Log2FC)", title = paste0(''))
  
  #Q40
  df.q40 <- data[,colnames(data)[grepl('Q40',colnames(data))]]
  colnames(df.q40) <- colnames(data)[grepl('Q40',colnames(data))] %>% gsub('.Q40','',.) %>% gsub('.t2','',.)
  t2.q40 <- get_fit(df.q40,initial.t2)
  fit.t2.q40 <- t2.q40[[1]]
  k1.t2.q40=coef(fit.t2.q40)[1]
  k2.t2.q40=coef(fit.t2.q40)[2]
  #Q30
  df.q30 <- data[,colnames(data)[grepl('Q30',colnames(data))]]
  colnames(df.q30) <- colnames(data)[grepl('Q30',colnames(data))] %>% gsub('.Q30','',.) %>% gsub('.t2','',.)
  t2.q30 <- get_fit(df.q30,initial.t2)
  fit.t2.q30 <- t2.q30[[1]]
  k1.t2.q30=coef(fit.t2.q30)[1]
  k2.t2.q30=coef(fit.t2.q30)[2]
  
  p.t2 <-  ggplot() +
    geom_point(data=data, aes(x = x.Q40, y = y.t2.Q40),color='#4472CA',alpha=0.5) +  
    geom_line(data= theoretical.t2,aes(x = x, y = y), color = '#9FE61F', linetype = "dashed")+
    geom_smooth(data=data, aes(x = x.Q40, y = y.t2.Q40),method = "nls",method.args=list(formula =y ~ log2(k1 + k2 * 2^x),
                                                                                        start=list(k1=k1.t2.q40,k2=k2.t2.q40)),
                se = FALSE, color = '#4472CA', linetype = "dashed")+
    geom_point(data=data, aes(x = x.Q30, y = y.t2.Q30),color='#E69F00',alpha=0.5) + 
    geom_smooth(data=data, aes(x = x.Q30, y = y.t2.Q30),method = "nls",method.args=list(formula =y ~ log2(k1 + k2 * 2^x),
                                                                                        start=list(k1=k1.t2.q30,k2=k2.t2.q30)),
                se = FALSE, color = '#E69F00', linetype = "dashed")+
    theme_nature_border()+coord_fixed()+
    labs(x = "M8.vs.D6 (Log2FC)", y = "T2.vs.D6 (Log2FC)", title = paste0(''))
  
  
  # 使用 predict() 函数和拟合模型生成预测的 y 值
  predicted_y.t1.q30 <- predict(fit.t1.q30, newdata = data.frame(x = data$x.Q30))
  predicted_y.t1.q40 <- predict(fit.t1.q40, newdata = data.frame(x = data$x.Q40))
  
  predicted_y.t2.q30 <- predict(fit.t2.q30, newdata = data.frame(x = data$x.Q30))
  predicted_y.t2.q40 <- predict(fit.t2.q40, newdata = data.frame(x = data$x.Q40))

  data$y.t1.Q40.pred <- predicted_y.t1.q40
  data$y.t2.Q40.pred <- predicted_y.t2.q40
  data$y.t1.Q30.pred <- predicted_y.t1.q30
  data$y.t2.Q30.pred <- predicted_y.t2.q30
  lst = list(p.t1,p.t2,data)
  return(lst)
}

library(tidyr)
ele_diff.wide <- pivot_wider(ele_diff, names_from = compare, values_from = meanlogFC) %>% as.data.frame()
ele_diff.wide <- na.omit(ele_diff.wide)

colnames(ele_diff.wide) <- c('gene','x.Q40','y.t1.Q40','y.t2.Q40')

#ele.lst <- get_fit(ele_diff.wide,'Q40')
#ele.lst[[1]]

ilm_diff.wide <- pivot_wider(ilm_diff, names_from = compare, values_from = meanlogFC) %>% as.data.frame()
ilm_diff.wide <- na.omit(ilm_diff.wide)

colnames(ilm_diff.wide) <- c('gene','x.Q30','y.t1.Q30','y.t2.Q30')
#ilm.lst <- get_fit(ilm_diff.wide,'Q30')
#ele.lst

diff.wide <- merge(ele_diff.wide,ilm_diff.wide,by='gene')

fit_lst <- get_fit.plot(diff.wide)

library(patchwork)
p_all.fit = fit_lst[[1]] + fit_lst[[2]]

diff.pred <- fit_lst[[3]]


diff.pred.RMSE <- ddply(diff.pred, .(), summarise, RMSE.Q40.t1 = get_RMSE(y.t1.Q40 ,y.t1.Q40.pred),
      RMSE.Q30.t1 = get_RMSE(y.t1.Q30 ,y.t1.Q30.pred),
      RMSE.Q40.t2 = get_RMSE(y.t2.Q40 ,y.t2.Q40.pred),
      RMSE.Q30.t2 = get_RMSE(y.t2.Q30 ,y.t2.Q30.pred))



#SNR--
fpkm_df <- read.csv('RNA_quartet_fpkm.csv')
#D53D61 -> M8:D6 = 3:1 -> T1
#D51D63 -> M8:D6 = 1:3 -> T2

#fpkm_df_sub <- fpkm_df[,colnames(fpkm_df)[grepl('M8|D6',colnames(fpkm_df))]]

fpkm_df_sub <- lapply(fpkm_df, function(x) ifelse(x == 0, NA, x)) %>% as.data.frame()
fpkm_df_sub$gene <- fpkm_df$gene
colnames(fpkm_df_sub) <- colnames(fpkm_df_sub) %>% gsub('FPKM.','',.) %>% gsub('.10G','',.) %>% gsub('D53D61','T1',.) %>% gsub('D51D63','T2',.)
get_SNR_plot <- function(dt_fpkm,dt_counts,dt_meta,platform,d_type){
  change_cols <- colnames(dt_fpkm[, !'GeneName'])
  dt_fpkm[, (change_cols):= lapply(.SD, as.numeric), .SDcols = change_cols]
  dt_counts[, (change_cols):= lapply(.SD, as.numeric), .SDcols = change_cols]
  
  dt_fpkm_log <- data.table(apply(dt_fpkm[, !'GeneName'], 2, function(x)(log2(x + 0.01))))
  dt_fpkm_log[, GeneName := dt_fpkm$GeneName]
  
  calc_signoise_ratio <- function(pca_prcomp, exp_design) {
    
    pcs <- as.data.frame(predict(pca_prcomp))
    dt_perc_pcs <- data.table(PCX = 1:nrow(pcs),
                              Percent = summary(pca_prcomp)$importance[2,],
                              AccumPercent = summary(pca_prcomp)$importance[3,])
    
    dt_dist <- data.table(ID.A = rep(rownames(pcs), each = nrow(pcs)),
                          ID.B = rep(rownames(pcs), time = nrow(pcs)))
    
    dt_dist$Group.A <- exp_design[dt_dist$ID.A]$group
    dt_dist$Group.B <- exp_design[dt_dist$ID.B]$group
    
    dt_dist[, Type := ifelse(ID.A == ID.B, "Same",
                             ifelse(Group.A == Group.B, "Intra", "Inter"))]
    dt_dist[, Dist := dt_perc_pcs[1]$Percent * (pcs[ID.A, 1] - pcs[ID.B, 1]) ^ 2 + dt_perc_pcs[2]$Percent * (pcs[ID.A, 2] - pcs[ID.B, 2]) ^ 2]
    
    dt_dist_stats <- dt_dist[, .(Avg.Dist = mean(Dist)), by = .(Type)]
    setkey(dt_dist_stats, Type)
    signoise <- dt_dist_stats["Inter"]$Avg.Dist / dt_dist_stats["Intra"]$Avg.Dist
    signoise_db <- 10*log10(signoise)
    return(signoise_db)
  }
  
  get_pca_list <- function(expr_mat_forsignoise, exp_design, dt_meta) {
    pca_prcomp = prcomp(t(expr_mat_forsignoise), scale = F)
    pcs = predict(pca_prcomp) %>% data.frame()
    pcs$library = row.names(pcs)
    pcs_add_meta = merge(pcs, dt_meta, by = "library")
    PC1_ratio = round(summary(pca_prcomp)$importance[2, 1] * 100, digits = 2)
    PC2_ratio = round(summary(pca_prcomp)$importance[2, 2] * 100, digits = 2)
    #PC3_ratio = round(summary(pca_prcomp)$importance[2, 3] * 100, digits = 2)
    PC3_ratio = NA
    SNR = format(round(calc_signoise_ratio(pca_prcomp, exp_design = exp_design), digits = 3), nsmall = 3) 
    gene_num = dim(expr_mat_forsignoise)[1]
    pca_list = cbind(pcs_add_meta, PC1_ratio, PC2_ratio, PC3_ratio, SNR, gene_num)
    return(pca_list)
  }
  
  output_snr_res <- function(dt_fpkm_log, dt_counts, dt_meta){
    dt_detect_gene <- do.call(cbind, lapply(unique(dt_meta$sample), function(x){
      detect_res <- apply(dt_counts[, dt_meta[sample == x][['library']], with = F], 1, function(x){length(which(x >= 3)) >= 2})
      return(
        detect_res = detect_res
      ) 
    }))
    
    gene_list_snr <- dt_counts[['GeneName']][apply(dt_detect_gene, 1, function(x){any(x)})]
    exp_design = (dt_meta[, .(library, group = sample)] %>% setkey(., library))
    dt_fpkm_log <- aggregate(.~GeneName,dt_fpkm_log,mean) %>% data.table()
    dt_fpkm_f <- dt_fpkm_log[gene_list_snr, on = .(GeneName)]
    dt_fpkm_zscore <- data.table(t(apply(dt_fpkm_f[, dt_meta$library, with = F], 1, function(x){(x - mean(x))/sd(x)})))
    dt_fpkm_zscore <- na.omit(dt_fpkm_zscore)
    pca_list <- get_pca_list(dt_fpkm_zscore, exp_design, dt_meta)
    return(pca_list)
  }
  
  dt_snr <- output_snr_res(dt_fpkm_log, dt_counts, dt_meta)
  snr_gene_num <- dt_snr$gene_num[1]
  
  ## figure of pca with snr
  print(dt_snr)
  pt_snr <- ggplot(dt_snr, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 3, show.legend = FALSE) +
    theme_few() +
    guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 1, title.position = "top")) +
    scale_fill_manual(values = quartet_color) +
    scale_color_manual(values = quartet_color) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(
      title = paste0(platform,'',d_type,"\nSNR: ", dt_snr$SNR[1], ' (N = ', dt_snr$gene_num[1], ')'),
      x = paste("PC1 (", dt_snr$PC1_ratio, "%)", sep = ""),
      y = paste("PC1 (", dt_snr$PC2_ratio, "%)", sep = ""))+theme_nature_border()
  
  ## output snr table
  dt_snr$PC1 <- round(dt_snr$PC1 , digits = 3)
  dt_snr$PC2 <- round(dt_snr$PC2 , digits = 3)
  fwrite(dt_snr, file = paste(platform, ".pca_with_snr.txt", sep = ""), sep = "\t")
  return(list(unique(dt_snr$SNR),pt_snr))
}



fpkm_df.ele <- fpkm_df_sub[,grepl('ELE',colnames(fpkm_df_sub))]
fpkm_df.ele$GeneName <- fpkm_df_sub$gene

fpkm_df.ilm <- fpkm_df_sub[,grepl('ILM',colnames(fpkm_df_sub))]
fpkm_df.ilm$GeneName <- fpkm_df_sub$gene

dt_meta.ilm = data.frame(group=c(colnames(fpkm_df.ilm)[-19]),
                         library=c(colnames(fpkm_df.ilm)[-19]),
                         sample=str_split(colnames(fpkm_df.ilm)[-19],'\\_',simplify = T)[,2]) %>% data.table()

dt_meta.ele = data.frame(group=c(colnames(fpkm_df.ele)[-19]),
                         library=c(colnames(fpkm_df.ele)[-19]),
                         sample=str_split(colnames(fpkm_df.ele)[-19],'\\_',simplify = T)[,2]) %>% data.table()


ele.pca <- get_SNR_plot(fpkm_df.ele %>% data.table::data.table(),
                        fpkm_df.ele %>% data.table::data.table(),
                        dt_meta.ele,
                        'Q40','')

ilm.pca <- get_SNR_plot(fpkm_df.ilm %>% data.table::data.table(),
                        fpkm_df.ilm %>% data.table::data.table(),
                        dt_meta.ilm,
                        'Q30','')

ele.pca[[2]] + ilm.pca[[2]]

snr_lst <- list()
snr_plots <- list()
for (i in labels){
  print(i)
  aim_genes.ele <- subset(ele_expr,Interval==i)$gene
  aim_genes.ilm <- subset(ilm_expr,Interval==i)$gene
  snr.ele = get_SNR_plot(subset(fpkm_df.ele,GeneName %in% aim_genes.ele) %>% data.table::data.table(),
                         subset(fpkm_df.ele,GeneName %in% aim_genes.ele) %>% data.table::data.table(),
                         dt_meta.ele,
                         'Q40',i)
  snr.ilm = get_SNR_plot(subset(fpkm_df.ilm,GeneName %in% aim_genes.ilm) %>% data.table::data.table(),
                         subset(fpkm_df.ilm,GeneName %in% aim_genes.ilm) %>% data.table::data.table(),
                         dt_meta.ilm,
                         'Q30',i)
  snr_lst[[paste0('Q40 ',i)]] = snr.ele[[1]]
  snr_lst[[paste0('Q30 ',i)]] = snr.ilm[[1]]
  snr_plots[[paste0('Q30 ',i)]] = snr.ilm[[2]]
  snr_plots[[paste0('Q40 ',i)]] = snr.ele[[2]]
}


library(gridExtra)
Quartet_snr_detail <- grid.arrange(grobs = snr_plots, ncol = 4)

quartet_snr_df <- snr_lst %>% unlist() %>% as.data.frame()
quartet_snr_df$SNR <- quartet_snr_df$.
quartet_snr_df$Quality <- str_split(rownames(quartet_snr_df),' ',simplify = T)[,1]
quartet_snr_df$Interval <- str_split(rownames(quartet_snr_df),' ',simplify = T)[,2]
quartet_snr_df$Sample <- 'Quartet'

save(Quartet_snr_detail,ele.pca,ilm.pca,quartet_snr_df,diff.wide,diff.pred,diff.pred.RMSE,p_all.fit,file = 'Quartet_titration.RData')

ggsave('Quartet.SNR.png',Quartet_snr_detail,width=8.15*1.5, height=5.7*2.5,dpi = 300)

#Filter differentially expressed genes：
#t-test P < 0.05 & |log2FC| > 1
df = expr.ele[1:10,]
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

ele_expr$DEGtype <- get_DEGtype(ele_expr)
ilm_expr$DEGtype <- get_DEGtype(ilm_expr)
get_mcc <- function(tn,tp,fn,fp){
  print(paste(tn,tp,fn,fp,sep = ' '))
  
  sqrt_term = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc = (tp*tn - fp*fn)/sqrt_term
  return(mcc)
}

get_f1 <- function(tn,tp,fn,fp){
  print(paste(tn,tp,fn,fp,sep = ' '))
  
  precision = tp/(tp+fp)
  recall = tp/(tp+fn)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  return(f1_score)
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
f1_df = data.frame(TN=1,
                   TP=1,
                   FN=1,
                   FP=1,
                   Sample=1)
quartet_lodr_lst = list()
for (comp in unique(ele_expr$compare)){
  a = str_split(comp,'\\/',simplify = T)[,1]
  b = str_split(comp,'\\/',simplify = T)[,2]
  #ELE--
  quartet_ref.sub <- subset(quartet_ref,compare==comp)
  ele_expr.sub <- subset(ele_expr,compare==comp)
  ele_count <- quartet_count[,colnames(quartet_count)[(grepl(a,colnames(quartet_count)) | grepl(b,colnames(quartet_count))) & grepl('ELE',colnames(quartet_count))]]
  ele_count$GeneName <- quartet_count$gene
  ele_count$MnSignal <- rowMeans(ele_count[,-7])
  rownames(ele_expr.sub) <- ele_expr.sub$gene_id
  ele_expr.sub$Feature = ele_expr.sub$gene_id
  ele_expr.sub <- merge(ele_expr.sub,quartet_ref.sub[,c('meanlogFC','Feature')],by='Feature')
  ele_expr.sub$Var = abs(ele_expr.sub$logFC-ele_expr.sub$meanlogFC)
  print(paste0('ELE ',comp))
  f1.ele = get_mcc_all(ele_expr.sub,quartet_ref.sub,paste0('ELE ',comp))
  f1_df = rbind(f1_df,f1.ele)
  #ele_lodr <- get_lodr(ele_expr.sub,ele_count,quartet_ref.sub)
  #lodr_lst[[paste0('ELE ',comp)]] = ele_lodr
  
  #ILM---
  ilm_expr.sub <- subset(ilm_expr,compare==comp)
  ilm_count <- quartet_count[,colnames(quartet_count)[(grepl(a,colnames(quartet_count)) | grepl(b,colnames(quartet_count))) & grepl('ILM',colnames(quartet_count))]]
  ilm_count$GeneName <- quartet_count$gene
  ilm_count$MnSignal <- rowMeans(ilm_count[,-7])
  rownames(ilm_expr.sub) <- ilm_expr.sub$gene_id
  #ilm_lodr <- get_lodr(ilm_expr.sub,ilm_count,quartet_ref.sub)
  #lodr_lst[[paste0('ILM ',comp)]] = ilm_lodr
  print(paste0('ILM ',comp))
  #mcc.ilm = get_mcc_all(ilm_expr.sub,quartet_ref.sub)
  #mcc_lst[[paste0('ILM ',comp)]] = mcc.ilm
  f1.ilm = get_mcc_all(ilm_expr.sub,quartet_ref.sub,paste0('ILM ',comp))
  f1_df = rbind(f1_df,f1.ilm)
  
  if (F){
    for (deg in unique(ele_expr.sub$DEGtype)){
      print(deg)
      ele_expr.sub.sub <- subset(ele_expr.sub,DEGtype == deg)
      ilm_expr.sub.sub <- subset(ilm_expr.sub,DEGtype == deg)
      quartet_ref.sub.sub <- subset(quartet_ref.sub,DEGtype == deg)
      ele.lodr <- get_lodr(ele_expr.sub.sub,ele_count,quartet_ref.sub.sub)
      ilm.lodr <- get_lodr(ilm_expr.sub.sub,ilm_count,quartet_ref.sub.sub)
      quartet_lodr_lst[[paste0('ELE ',comp,' ',deg)]] = ele.lodr[[1]]
      quartet_lodr_lst[[paste0('ILM ',comp,' ',deg)]] = ilm.lodr[[1]]
    }
  }
  
}


quartet_lodr.df <- quartet_lodr_lst[[1]]
quartet_lodr.df$source <- names(quartet_lodr_lst)[1]
for (i in names(quartet_lodr_lst)[2:length(quartet_lodr_lst)]){
  print(i)
  df = quartet_lodr_lst[[i]]
  df$source = i
  quartet_lodr.df = rbind(quartet_lodr.df,df)
}

quartet_lodr.df$Quality = ifelse(grepl('ELE',quartet_lodr.df$source),'Q40','Q30')
quartet_lodr.df$DEGtype <- str_split(quartet_lodr.df$source,'\\ ',simplify = T)[,3]
quartet_lodr.df$LODR <- gsub('<','',quartet_lodr.df$Estimate) %>% as.numeric()
quartet_lodr.df$Ratio <- ifelse(quartet_lodr.df$Ratio==0,'<0.01',quartet_lodr.df$Ratio)
quartet_lodr.df$Sample <- 'Quartet'
quartet_lodr.df$group <- quartet_lodr.df$Quality
quartet_lodr.df$Compare <- str_split(quartet_lodr.df$source,'\\ ',simplify = T)[,2]


quartet_lodr.df$LODR <- ifelse(is.infinite(quartet_lodr.df$LODR),10,quartet_lodr.df$LODR)


quartet.non <- p_line.ratio(subset(quartet_lodr.df,DEGtype=='non-DEG'),'LODR',comp=T)

quartet.up <- p_line.ratio(subset(quartet_lodr.df,DEGtype=='up-regulate'),'LODR',comp=T)

quartet.down <- p_line.ratio(subset(quartet_lodr.df,DEGtype=='down-regulate'),'LODR',comp=T)

maqc.up <- p_line.ratio(subset(maqc_lodr.df,DEGtype=='up-regulate'),'LODR')
maqc.down <- p_line.ratio(subset(maqc_lodr.df,DEGtype=='down-regulate'),'LODR')

p.lodr <- (maqc.up|maqc.down)/quartet.up/quartet.down

F1.df <- f1_df[-1,]
F1.df <- rbind(F1.df,ilm_F1)
F1.df <- rbind(F1.df,ele_F1)
F1.df$F1.score <- apply(F1.df,1,function(x){
  #print(x)
  get_f1(x['TN'] %>% as.numeric(),
         x['TP']%>% as.numeric(),
         x['FN']%>% as.numeric(),
         x['FP']%>% as.numeric())
}) %>% unlist()
F1.df$Sample = ifelse(!grepl('MAQC',F1.df$Sample),paste0('Quartet ',F1.df$Sample),F1.df$Sample)
F1.df$source = F1.df$Sample

F1.df$Quality <- ifelse(grepl('ELE',F1.df$source),'Q40','Q30')
F1.df$Sample <-gsub('ELE |ILM ','',F1.df$source)

F1.df$F1.score <- as.numeric(F1.df$F1.score)
F1_score <- ggplot(F1.df,aes(Sample,F1.score,color=Quality,fill=Quality))+
  geom_bar(stat="summary",fun=mean,position=position_dodge(width = 1))+ 
  stat_summary(geom = "text", aes(label = paste(round(..y.., digits =2))),
               position = position_dodge(width = 1), vjust = -0.7,hjust=0.5,size = 3)+
  theme_bw()+
  scale_fill_manual(values = colors)+scale_color_manual(values = colors)+
  labs(x = NULL, y = 'F1.score', col = "Quality",size=12,face="bold") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(face="bold",color = "black"),
        axis.title = element_text(face="bold",size=12,color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face="bold",color = "black"))

save(maqc_lodr_lst,mcc_lst,quartet_lodr_lst,quartet_lodr.df,maqc_lodr.df,F1.df,file='LODR.F1score.Rdata')

ggsave('F1_score.quartet.maqc.png',F1_score,width=8.15*0.8, height=5.7*0.8,dpi = 300)


