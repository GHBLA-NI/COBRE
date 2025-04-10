setwd("/nfs/dcmb-lgarmire/boweil/ED/GITHUB_CODE/data")

library(dplyr)
library(stringr)
library(stats)
library(corrplot)
library(pathifier)
library(princurve)
library(caret)
library(gdata)

source('../lilikoipack/stdev.R')
source('../lilikoipack/lilikoi_MetaTOpathway.R')
source('../lilikoipack/lilikoi_PDSfun.R')
source('../lilikoipack/pathifier.R')
source('../lilikoipack/getpathway.R') 
source('../lilikoipack/score.R')
source('../lilikoipack/scorepath.R')

load('../lilikoipack/sysdata.rda',verbose=TRUE)
Delete_items<-read.csv("kegg_supergroup_anno.csv")
Delete_items<-as.data.frame(Delete_items)
Exclude <- Delete_items[
  Delete_items[,"group"] == "Global and overview maps" |
  Delete_items[,"supergroup"] == "Drug Development" |
  Delete_items[,"supergroup"] == "Human Diseases" |
  Delete_items[,"name"] %in% c("Porphyrin and chlorophyll metabolism", 
                                  "Protein digestion and absorption", 
                                  "Mineral absorption"), 
]


keep_items<-Delete_items[-which(Delete_items$id %in% Exclude$id),]




index_seq_shap <- read.table("../output/index_sequence_update.txt",header = FALSE,sep = "\n",stringsAsFactors = FALSE)
rename_df_reor<-read.csv("../data/rename_df_reor.csv")
Metadata_orig <- read.csv('../output/knn_vsn_aftercombat_meta.csv', row.names=1)
Metadata_orig<-t(Metadata_orig)

new_names <- rownames(Metadata_orig)
matches <- match(new_names, rename_df_reor$clean)
new_names[!is.na(matches)] <- rename_df_reor$lili[matches[!is.na(matches)]]
rownames(Metadata_orig) <- new_names

selected_rows <- as.character(index_seq_shap$V1)
Metadata <- Metadata_orig[selected_rows, , drop = FALSE]
stopifnot(identical(rownames(Metadata), selected_rows))
index_rows <- as.character(index_seq_shap$V1)
rename_df_reor$clean <- as.character(rename_df_reor$clean)
matched_reor <- rename_df_reor[match(index_rows, rename_df_reor$clean), ]
new_names <- matched_reor$lili
new_names[is.na(new_names)] <- rownames(Metadata)[is.na(new_names)]
rownames(Metadata) <- new_names
dataSet <- list()
dataSet$cmpd<-(rownames(Metadata))
Metadata<-t(Metadata)
ClinicalInformation<-read.csv('complete_CI.csv',row.names=1)
#Make sure that the row matches
ClinicalInformation<-ClinicalInformation[match(rownames(Metadata),rownames(ClinicalInformation)),]
Metadata<-cbind(Metadata,ClinicalInformation$OBS_NORMAL_con)
Metadata<-cbind(Metadata,ClinicalInformation$label)
Metadata <- cbind(Metadata,
                  OBS_NORMAL_con = ClinicalInformation$OBS_NORMAL_con,
                  label = ClinicalInformation$label)
newList <- list("Metadata"=Metadata, "dataSet"=dataSet)





metabolites.list<-metabolites.list[names(metabolites.list)%in%str_to_title(keep_items$name)]
pathway.list<-pathway.list[(pathway.list%in%str_to_title(keep_items$name))]
metabolites.list[["One Carbon Metabolism"]]<-c('HMDB0001409','HMDB0001227','HMDB01846','HMDB01396','HMDB00108','HMDB00043','HMDB00123','HMDB0304213','HMDB15532','HMDB00244','HMDB00972','HMDB01354','HMDB00097','HMDB00121','HMDB00121','HMDB01056','HMDB01354','HMDB00696','HMDB00939','HMDB01185','HMDB02174','HMDB00187','HMDB02274','HMDB00676','HMDB01491','HMDB00742')
pathway.list<-c(pathway.list,"One Carbon Metabolism")


dataSet<-newList$dataSet
convertResults<-lilikoi.MetaTOpathway('name')



table<-convertResults$table
convertResults$table<-convertResults$table[!(convertResults$table$Match=="NA"),]
tab<-convertResults$table[!(is.na(convertResults$table$pathway)),]
Metabolite_pathway_table<-convertResults$table
Metadata<-newList$Metadata
Metadata<-as.data.frame(Metadata)
convert<-convertResults$table



PDSmatrix=lilikoi.PDSfun(convert)


df <- as.data.frame(do.call(rbind, PDSmatrix$score))
rownames(df) <- names(PDSmatrix$score)
pds_matrix<-df

ClinicalInformation <- ClinicalInformation[rownames(ClinicalInformation) != "MT65_rep", ]


# Create lists to store both p-values and correlation coefficients
p_values <- list()
cor_coefficients <- list()

for (i in 1:nrow(pds_matrix)){
  pathway_score <- pds_matrix[i,]
  cor_test <- cor.test(t(pathway_score), ClinicalInformation$OBS_NORMAL_con, method='kendall')
  p_values[i] <- cor_test$p.value
  cor_coefficients[i] <- cor_test$estimate  
}

# Convert to named vectors for easier interpretation
names(p_values) <- rownames(pds_matrix)
names(cor_coefficients) <- rownames(pds_matrix)
p_val_adj<-p.adjust(p_values,method='fdr')




sig_pathways <- rownames(pds_matrix)[p_values < 0.05]
sig_pathways_adj <- rownames(pds_matrix)[p_val_adj < 0.05]


sig_pathway_direction <- data.frame(
  Pathway = sig_pathways,
  Correlation = unlist(cor_coefficients[which(p_values < 0.05)]),
  Direction = ifelse(unlist(cor_coefficients[which(p_values < 0.05)]) > 0, "Positive", "Negative"),
  P_value = unlist(p_values[which(p_values < 0.05)])
)

sig_pathway_direction <- sig_pathway_direction[order(abs(sig_pathway_direction$Correlation), decreasing=TRUE),]



gene_df <- do.call(rbind, lapply(names(PDSmatrix$genesinpathway), function(pathway) {
  genes <- PDSmatrix$genesinpathway[[pathway]]

  genes <- as.vector(genes)
  data.frame(Pathway = pathway, Metabolites = genes, stringsAsFactors = FALSE)
}))

sig_gene_df <- gene_df[gene_df$Pathway %in% sig_pathways, ]



sig_pathway_direction
sig_gene_df