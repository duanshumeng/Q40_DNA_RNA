library(plyr)

setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/ERCC_RNAseq')
load("ELE.M8D6.count.M8.D6.RData")
load("ILM.M8D6.count.M8.D6.RData")
#Count----
ele_sd <- ELE.M8D6.count.M8.D6.exDat$Results$maDatAll
#ilm_sd <- ILM.M8D6.count.M8.D6.exDat$Results$dynRangeDat
ilm_sd <- ILM.M8D6.count.M8.D6.exDat$Results$maDatAll
ele_sd$Quality <- 'Q40'
ilm_sd$Quality <- 'Q30'
ercc_sd <- rbind(ele_sd,ilm_sd)
#ercc_sd<- subset(ercc_sd,Ratio !='1:1')
#x_lim <- c('4:1','1:2','1:1.5')

ercc_sd$group <- paste0(ercc_sd$Ratio,' ',ercc_sd$Quality)
ercc_count <- table(ercc_sd$group) %>% as.data.frame()
ercc_count$Quality <- str_split(ercc_count$Var1,' ',simplify = T)[,2]
ercc_count$Ratio <- str_split(ercc_count$Var1,' ',simplify = T)[,1]

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

get_log2FC <- function(ratio){
  ratio <- as.character(ratio)
  ratio_1 <- as.numeric(strsplit(ratio,':')[[1]][1]) / as.numeric(strsplit(ratio,':')[[1]][2])
  log2_ratio <- log2(ratio_1)
  return(log2_ratio)
}

ercc_sd <- ddply(ercc_sd, .(Feature), mutate, log2FC = get_log2FC(Ratio))

ercc_rmse <- ddply(ercc_sd, .(group), summarise, RMSE = get_RMSE(log2FC,M.Ave))
ercc_rmse$Ratio <- str_split(ercc_rmse$group,' ',simplify = T)[,1]
ercc_rmse$Quality <- str_split(ercc_rmse$group,' ',simplify = T)[,2]

#AUC----
ele_auc <- ELE.M8D6.count.M8.D6.exDat$Results$AUCdat
ilm_auc <- ILM.M8D6.count.M8.D6.exDat$Results$AUCdat

ele_auc$Quality <- 'Q40'
ilm_auc$Quality <- 'Q30'

ercc_auc <- rbind(ele_auc,ilm_auc)

#LODR----
ele_lodr <- ELE.M8D6.count.M8.D6.exDat$Results$lodr.res.ERCC
ilm_lodr <- ILM.M8D6.count.M8.D6.exDat$Results$lodr.res.ERCC
ele_lodr$Quality = 'Q40'
ilm_lodr$Quality = 'Q30'
ercc_lodr <- rbind(ele_lodr,ilm_lodr)
ercc_lodr$Estimate <- as.numeric(ercc_lodr$Estimate)
ercc_lodr$Estimate[is.infinite(ercc_lodr$Estimate)] <- 450

save(ercc_sd, ercc_auc,ercc_count,ercc_lodr,ercc_rmse, file = "/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/ERCC.rnaseq.RData")
