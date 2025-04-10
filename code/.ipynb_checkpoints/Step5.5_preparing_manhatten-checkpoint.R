load("../output/multinomresult_Momage_BMI.RData")
Asian<-sapply(p_values, function(x) x[[1]])
Mixed<-sapply(p_values, function(x) x[[2]])
White<-sapply(p_values, function(x) x[[3]])
df<-data.frame(Asian=Asian,Mixed=Mixed,White=White)
df_adjusted<-as.data.frame(lapply(df,function(x) p.adjust(x,method='fdr')))
rownames(df_adjusted)<-row.names(df)
df_adjusted_reor<-df_adjusted[order(df_adjusted$White),]
selected_WH <-df_adjusted_reor[which(df_adjusted_reor$White<0.05),]
selected_AS <-df_adjusted_reor[which(df_adjusted_reor$Asian<0.05),]   
selected_Mix<-df_adjusted_reor[which(df_adjusted_reor$Mixed<0.05),] 
snplist <- rownames(selected_WH)
snplist <- gsub("\\.", "-", snplist)
snplist <- substr(snplist, 1, nchar(snplist) - 2)
write.csv(selected_WH,"../output/adjp_val_multinomial_momage.csv")
                                  
                
                                  
                                  
                                  
bim <- read.table("../data/PE_maternal_09_10_control_region_filtered.bim", header = FALSE, tringsAsFactors = FALSE,col.names = c("CHR", "SNP", "cm", "BP","allele1", "allele2"))

df_adjusted_reor<-read.csv("../data/df_adjusted_reor.csv",row.names = 1)
new_names <- substr(rownames(df_adjusted_reor), 1, nchar(rownames(df_adjusted_reor)) - 2)
new_names <- gsub("\\.", "-", new_names)
rownames(df_adjusted_reor) <- new_names

all_snps<-  bim[, c("CHR", "SNP", "BP","allele2","allele1")]
all_snps$White <-df_adjusted_reor$White[match(all_snps$snp, rownames(df_adjusted_reor))]
all_snps$Mixed <-df_adjusted_reor$Mixed[match(all_snps$snp, rownames(df_adjusted_reor))]
all_snps$Asian <-df_adjusted_reor$Asian[match(all_snps$snp, rownames(df_adjusted_reor))]                                  
all_snps<-na.omit(all_snps)

all_snps_filtered <- all_snps[all_snps$White < 0.05, ]
                                  
write.csv(all_snps_filtered,"../output/sigsnps_info.csv")

                                  
                                  
                                  
all_snps$p <-df_adjusted_reor$White[match(all_snps$snp, rownames(df_adjusted_reor))]
names(all_snps)[names(all_snps) == "p"]   <- "P"
write.csv(all_snps,"../output/all_snps_adjp_val_multinomial_momage.csv")
                                  
                                  

target_snps <- snplist
selected_snps <- bim[bim$snp %in% target_snps, c("snp", "pos")]
write.csv(selected_snps,"slected_snps.csv")                                  