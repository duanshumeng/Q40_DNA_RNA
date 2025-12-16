library(stringr)
library(HTqPCR)
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/MAQC2_sampleA_sampleB')
rna_count <- read.csv('../gene_count_matrix.ALL.csv',header = T)
taq_count <- read.table('MAQC_TAQ_16r1044_ANA_gene_20150804.txt',header = T)
taq_count$GeneName %>% length()

gene_id <- data.frame(GeneName=str_split(rna_count$gene_id,'\\|',simplify = T)[,1],
                      GeneId=str_split(rna_count$gene_id,'\\|',simplify = T)[,2])

rna_count$GeneName <- str_split(rna_count$gene_id,'\\|',simplify = T)[,2]

rna_count.sub <- rna_count[,grepl('RNA_A|RNA_B|GeneName',colnames(rna_count))]

rna_count.ele <- rna_count.sub[,grepl('_ELE|GeneName',colnames(rna_count.sub))]

rna_count.ilm <- rna_count.sub[,grepl('_ILM|GeneName',colnames(rna_count.sub))]

subset(taq_count,GeneName=='NGFR')
subset(rna_count.ele,GeneName=='NGFR')


#PCC with Taqman quantification(RC)------------------
library(reshape2)
library(dplyr)
#Taq
polr2a.count <- subset(taq_count,GeneName=='POLR2A') %>% melt()
polr2a.count$Sample <- str_split(polr2a.count$variable,'\\_',simplify = T)[,4] %>% gsub('1|2|3|4','',.)

polr2a.mean <- polr2a.count %>% group_by(Sample) %>% summarize(count.mean=mean(value))

taq_count_AB = taq_count[,grepl('GeneName|_A|_B',colnames(taq_count))]


for (x in colnames(taq_count_AB)){
  if (grepl('_A',x)){
    print(x)
    m = polr2a.mean[polr2a.mean$Sample=='A',2] %>% as.numeric()
    taq_count_AB[paste0(x,' norm')] = taq_count_AB[x] - m
  } else if (grepl('_B',x)){
    m = polr2a.mean[polr2a.mean$Sample=='B',2] %>% as.numeric()
    taq_count_AB[paste0(x,' norm')] = taq_count_AB[x] - m
  }
  
}

taq_count_AB.norm <- taq_count_AB[,grepl('norm',colnames(taq_count_AB))]
rownames(taq_count_AB.norm) <- taq_count_AB$GeneName

colnames(taq_count_AB.norm) <- gsub('.norm','',colnames(taq_count_AB.norm))

taq_count_AB.norm = taq_count[,grepl('_A|_B',colnames(taq_count))]
rownames(taq_count_AB.norm) <- taq_count$GeneName
smp_list <- str_split(colnames(taq_count_AB.norm),'\\_',simplify = T)[,4] %>% gsub('1|2|3|4','',.)
names(smp_list) <- colnames(taq_count_AB.norm)

taq_count_AB.norm$A.Mean <- apply(taq_count_AB.norm,1,function(x){
  A.mean = mean(x[grepl('_A',names(x))])
  return(A.mean)
})

taq_count_AB.norm$B.Mean <- apply(taq_count_AB.norm,1,function(x){
  A.mean = mean(x[grepl('_B',names(x))])
  return(A.mean)
})

taq_count_AB.norm$logFC <- apply(taq_count_AB.norm,1,function(x){
  return(log2(x['A.Mean']/x['B.Mean']))
})

taq_expr <- taq_count_AB.norm[,c('A.Mean','B.Mean','logFC')]

if (F){
  list <- model.matrix(~factor(smp_list)+0)
  colnames(list) <- c('A','B')
  
  taq.fit <- lmFit(taq_count_AB.norm[,1:8],list)
  
  df.matrix <- makeContrasts(A - B , levels = list)
  fit <- contrasts.fit(taq.fit, df.matrix)
  fit <- eBayes(fit)
  taq_expr <- topTable(fit,n = Inf, adjust = "fdr")
}


#rna-seq
library(plyr)

get_expr <- function(df_count){
  if (F){
    polr2a.count <- subset(df_count,GeneName=='POLR2A') %>% melt()
    polr2a.count$Sample <- str_split(polr2a.count$variable,'\\_',simplify = T)[,2] %>% gsub('1|2|3|4','',.) %>% as.factor()
    polr2a.mean <- polr2a.count %>% group_by(Sample) %>% dplyr::summarize(count.mean=mean(value))
    
    
    for (x in colnames(df_count)){
      if (grepl('_A',x)){
        print(x)
        m = polr2a.mean[polr2a.mean$Sample=='A',2] %>% as.numeric()
        df_count[paste0(x,' norm')] = df_count[x] /m
      } else if (grepl('_B',x)){
        m = polr2a.mean[polr2a.mean$Sample=='B',2] %>% as.numeric()
        df_count[paste0(x,' norm')] = df_count[x] /m
      }
      
    }
    df_count.norm <- df_count[,grepl('norm|GeneName',colnames(df_count))]
    colnames(df_count.norm) <- colnames(df_count.norm)  %>% gsub(' norm','',.)
    df_count.mean <- aggregate(.~GeneName,df_count.norm,mean)
  }

  df_count.mean <- aggregate(.~GeneName,df_count,mean)
  rownames(df_count.mean) <- df_count.mean$GeneName 
  df_count.mean = df_count.mean[,-1] %>% log()
  df_count.mean = df_count.mean[,-1]
  smp_list <- str_split(colnames(df_count.mean),'\\_',simplify = T)[,2]
  list <- model.matrix(~factor(smp_list)+0)
  colnames(list) <- c('A','B')
  df.fit <- lmFit(df_count.mean, list)
  
  df.matrix <- makeContrasts(A - B , levels = list)
  fit <- contrasts.fit(df.fit, df.matrix)
  fit <- eBayes(fit)
  tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
  return(tempOutput)
}

expr.ele <- get_expr(rna_count.ele)
expr.ele$Genename <- rownames(expr.ele)
expr.ilm <- get_expr(rna_count.ilm)
expr.ilm$Genename <- rownames(expr.ilm)
taq_expr$Genename <- rownames(taq_expr)

save(expr.ele, expr.ilm,taq_expr, file = "MAQC.rnaseq.exp.RData")
load("MAQC.rnaseq.exp.RData")

breaks <- c(0, 0.99, 1.99,2.99,3.99, 4.99, 5.99,10)
# Create labels for the intervals
labels <- c("[0,1)", "[1,2)", "[2,3)", "[3,4)","[4,5)","[5,6)", ">=6")


#PCC of TaqMam and Element sequencing
expr.ele.taq <- merge(expr.ele[,c("logFC","Genename")],taq_expr[,c("logFC","Genename")],by = c("Genename"),suffixes=c(".ELE",".TaqMan"))
expr.ele.taq$Interval <- cut(abs(expr.ele.taq$logFC.TaqMan), breaks = breaks, labels = labels, include.lowest = TRUE)

check_group_size <- function(group_data) {
  if (length(group_data) < 4) {  
    return(FALSE)
  } else {
    return(TRUE)
  }
}
expr.ele.taq <- na.omit(expr.ele.taq)
cor.ele.taq <- ddply(expr.ele.taq, .(Interval), function(sub_data) {
  if (check_group_size(sub_data$Interval)) {
    correlation <- cor.test(sub_data$logFC.ELE, sub_data$logFC.TaqMan, method = "pearson",na.rm=TRUE)$estimate
  } else {
    correlation <- NA  
  }
  return(data.frame(correlation = correlation))
})
cor.ele.taq$Quality <- 'Q40'


colors = c('#4472CA','#E69F00')
names(colors) <- c('ELE','ILM')
get_VAF.ele <- function(df,platform){
  p_vaf_genome <- ggplot(df, aes(x = `logFC.ELE`, y = `logFC.TaqMan`)) + theme_classic()+
    geom_point(size=3,alpha=0.5,color = '#4472CA')+ 
    geom_smooth(method = "lm", se = FALSE, color = 'black') +stat_cor(method = "pearson",p.accuracy = 0.001,size = 5)+
    theme_nature_border()+coord_fixed()+xlim(c(-15,15))+ylim(c(-15,15))+
    labs(x = "logFC.ELE", y = "logFC.TaqMan", title = paste0(platform))
  return(p_vaf_genome)
}

expr.ele.taq.p <- get_VAF.ele(expr.ele.taq,'ELE.vs.TaqMan')

#PCC of TaqMan and Illumina sequencing
expr.ilm.taq <- merge(expr.ilm[,c("logFC","Genename")],taq_expr[,c("logFC","Genename")],by = c("Genename"),suffixes=c(".ILM",".TaqMan"))
expr.ilm.taq$Interval <- cut(abs(expr.ilm.taq$logFC.TaqMan), breaks = breaks, labels = labels, include.lowest = TRUE)

expr.ilm.taq <- na.omit(expr.ilm.taq)
cor.ilm.taq <- ddply(expr.ilm.taq, .(Interval), function(sub_data) {
  if (check_group_size(sub_data$Interval)) {
    correlation <- cor.test(sub_data$logFC.ILM, sub_data$logFC.TaqMan, method = "pearson",na.rm=TRUE)$estimate
  } else {
    correlation <- NA 
  }
  return(data.frame(correlation = correlation))
})
cor.ilm.taq$Quality <- 'Q30'

maqc_cor.df <- rbind(cor.ilm.taq,cor.ele.taq)
cor.ratio.tag <- p_line.ratio(maqc_cor.df,'correlation')

#Total count----
maqc_count.ilm <- table(expr.ilm.taq$Interval) %>% as.data.frame()
maqc_count.ilm$Quality = 'Q30'
maqc_count.ele <- table(expr.ele.taq$Interval) %>% as.data.frame()
maqc_count.ele$Quality = 'Q40'
maqc_count.all<- rbind(maqc_count.ilm,maqc_count.ele)
colnames(maqc_count.all) <- c('Interval','Freq','Quality')
maqc_count.all$Sample <- 'MAQC'

#RMSE----
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

rmse.ele <- ddply(expr.ele.taq, .(Interval), summarise, RMSE = get_RMSE(logFC.TaqMan,logFC.ELE))
rmse.ele$Quality <- 'Q40'
rmse.ilm <- ddply(expr.ilm.taq, .(Interval), summarise, RMSE = get_RMSE(logFC.TaqMan,logFC.ILM))
rmse.ilm$Quality <- 'Q30'

maqc_rmse <- rbind(rmse.ele,rmse.ilm)

rmse.ratio <- p_line.ratio(maqc_rmse,'RMSE')


#PCA分析(SNR)-----------------------
## obtain SNR results
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
    PC3_ratio = round(summary(pca_prcomp)$importance[2, 3] * 100, digits = 2)
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
  pt_snr <- ggplot(dt_snr, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 3, show.legend = FALSE) +
    theme_few() +
    guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 1, title.position = "top")) +
    scale_fill_manual(values = c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745")) +
    scale_color_manual(values = c("#2f5c85", "#7BC8A4", "#FFC65D", "#F16745")) +
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
library(data.table)
rna_fpkm <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/RNA_seqc_fpkm.csv',header = T)
rna_fpkm$GeneName <- str_split(rna_fpkm$gene,'\\|',simplify = T)[,1]

rna_fpkm.ele <- rna_fpkm[,grepl('_ELE|GeneName',colnames(rna_fpkm))]

rna_fpkm.ilm <- rna_fpkm[,grepl('_ILM|GeneName',colnames(rna_fpkm))]

colnames(rna_fpkm.ele)[-7] <- gsub('FPKM.','',colnames(rna_fpkm.ele)[-7]) %>% gsub('.10G','',.)
dt_meta.ele = data.frame(group=c(colnames(rna_fpkm.ele)[-7]),
                         library=c(colnames(rna_fpkm.ele)[-7]),
                         sample=str_split(colnames(rna_fpkm.ele)[-7],'\\_',simplify = T)[,2]) %>% data.table()

rna_fpkm.ele <- merge(rna_fpkm.ele,gene_id)
rna_fpkm.ele$GeneName  <- rna_fpkm.ele$GeneId
rna_fpkm.ele <- rna_fpkm.ele[,-8]

colnames(rna_fpkm.ilm)[-7] <- gsub('FPKM.','',colnames(rna_fpkm.ilm)[-7]) %>% gsub('.10G','',.)
dt_meta.ilm = data.frame(group=c(colnames(rna_fpkm.ilm)[-7]),
                         library=c(colnames(rna_fpkm.ilm)[-7]),
                         sample=str_split(colnames(rna_fpkm.ilm)[-7],'\\_',simplify = T)[,2]) %>% data.table()

rna_fpkm.ilm <- merge(rna_fpkm.ilm,gene_id)
rna_fpkm.ilm$GeneName  <- rna_fpkm.ilm$GeneId
rna_fpkm.ilm <- rna_fpkm.ilm[,-8]

snr_lst <- list()
snr_splot <- list()
for (i in labels){
  print(i)
  gene.ele = subset(expr.ele.taq,Interval==i)$Genename
  gene.ilm = subset(expr.ilm.taq,Interval==i)$Genename
  ele.pca <- get_SNR_plot(subset(rna_fpkm.ele,rna_fpkm.ele$GeneName %in% gene.ele) %>% data.table::data.table(),
                          subset(rna_count.ele,rna_count.ele$GeneName %in% gene.ele) %>% data.table::data.table(),
                          dt_meta.ele,
                          'Q40',i)
  snr_lst[[paste0('Q40 ',i)]] = ele.pca[[1]]
  snr_splot[[paste0('Q40 ',i)]] = ele.pca[[2]]
  ilm.pca <- get_SNR_plot(subset(rna_fpkm.ilm,rna_fpkm.ilm$GeneName %in% gene.ele) %>% data.table::data.table(),
                          subset(rna_count.ilm,rna_count.ilm$GeneName %in% gene.ele) %>% data.table::data.table(),
                          dt_meta.ilm,
                          'Q30',i)
  snr_lst[[paste0('Q30 ',i)]] = ilm.pca[[1]]
  snr_splot[[paste0('Q30 ',i)]] = ilm.pca[[2]]
  
}
maqc_snr_df <- snr_lst %>% unlist() %>% as.data.frame()
maqc_snr_df$SNR <- maqc_snr_df$.
maqc_snr_df$Quality <- str_split(rownames(maqc_snr_df),' ',simplify = T)[,1]
maqc_snr_df$Interval <- str_split(rownames(maqc_snr_df),' ',simplify = T)[,2]

snr.ratio <- p_line.ratio(maqc_snr_df,'SNR')
snr_detail.taq <- grid.arrange(grobs = snr_splot,ncol = 4)

maqc_cor.df$Sample <- 'MAQC'
maqc_rmse$Sample <- 'MAQC'
maqc_snr_df$Sample <- 'MAQC'
save(maqc_cor.df, maqc_count.all,maqc_rmse,maqc_snr_df, file = "MAQC.rnaseq.RData")
load("MAQC.rnaseq.RData")
