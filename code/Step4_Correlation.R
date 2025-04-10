library(Hmisc)
library(corrplot)
library(dummies)

# Load data
metab_info = read.csv('../output/knn_vsn_aftercombat_meta.csv', row.names=1)
clinical = read.csv('../data/complete_CI.csv',row.names=1)
name_map<-read.csv("../data/rename_df_reor.csv")
cobre_pd = read.csv('../data/cobre_pd.csv')
rownames(cobre_pd) = cobre_pd$Sample_Name
top_df <- read.table("../output/ranked_features_with_shap.txt",header = FALSE,sep = "\t",stringsAsFactors = FALSE)
colnames(top_df) <- c("Feature", "Value")
top_formatted <- top_df$Feature




clinical$Mat_Ethnicity = as.factor(clinical$Mat_Ethnicity)
levels(clinical$Mat_Ethnicity) = c('White','Asian','NHPI')
clinical$Pat_Ethnicity = as.factor(clinical$Pat_Ethnicity)
levels(clinical$Pat_Ethnicity) = c('White','Asian','NHPI')
clinical$OBS_NORMAL_con = clinical$OBS_NORMAL_con

comb_pd = merge(cobre_pd,clinical,by=0)
rownames(comb_pd) = comb_pd$Row.names
sample_eth_anno = comb_pd[,c('Row.names','Mat_Ethnicity.y','Pat_Ethnicity.y','OBS_NORMAL_con')]
sample_eth_anno = sample_eth_anno[,-c(1)]
colnames(sample_eth_anno) = c('Mat_Ethnicity','Pat_Ethnicity','BMI')


# Correct Metabolites name and filter metabolites data for top SHAP.
map_vector <- setNames(name_map$lili, name_map$clean)
current_names <- names(metab_info)
new_names <- ifelse(current_names %in% names(map_vector), 
                    map_vector[current_names], 
                    current_names)
names(metab_info) <- new_names

common_cols <- intersect(top_formatted, colnames(metab_info))



metab_eth_all <- merge(sample_eth_anno, metab_info[, common_cols], by = 0)
rownames(metab_eth_all) = metab_eth_all$Row.names
metab_eth_all_relevel = metab_eth_all[,-c(1)]
metab_eth_all_relevel$Mat_Ethnicity = as.numeric(factor(metab_eth_all_relevel$Mat_Ethnicity, levels = c("Asian" ,"White" ,"NHPI")))
metab_eth_all_relevel$Pat_Ethnicity = as.numeric(factor(metab_eth_all_relevel$Pat_Ethnicity, levels = c("Asian" ,"White" ,"NHPI")))


df_encoded <- dummy.data.frame(metab_eth_all, names = c("Mat_Ethnicity", "Pat_Ethnicity"), sep = "_")



# Run Spearman Correlation
spearman_corr = rcorr(as.matrix(metab_eth_all_relevel),type="spearman")
spearman_corr$P[is.na(spearman_corr$P)] = 0
png('../Figures/Fig5_correlation.png', 
    width = 8, height = 8, units = 'in', res = 600)

corrplot(
  spearman_corr$r, 
  method = "circle", 
  type = "upper", 
  order = "original", 
  tl.col = "black", 
  tl.srt = 45, 
  tl.cex = 0.5, 
  p.mat = spearman_corr$P, 
  sig.level = 0.05, 
  insig = "blank", 
  number.cex = 0.3, 
  addrect = 4
)

dev.off()