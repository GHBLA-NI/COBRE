library(limma)
library(sva)
library(mice)
library(impute)
library(tidyverse)
library(car)

# Data loading
ClinicalInformation <- read.csv("../data/ClinicalInformation.csv", row.names = 1)
merged123WO <- read.csv("../data/merged123WO.csv", row.names = 1)


# KNN imputation
knn_imputed <- impute.knn(as.matrix(merged123WO))$data
completed_data_knn <- as.data.frame(knn_imputed)


# normalization (vsn)
mat_vsn<- normalizeVSN(as.matrix(t(completed_data_knn)))
data_v <- as.data.frame(mat_vsn)
data_v <- t(data_v)




# Combat
common_samples <- intersect(rownames(data_v), rownames(ClinicalInformation))
data_v_matched <- data_v[common_samples, , drop = FALSE]
clinical_matched <- ClinicalInformation[common_samples, , drop = FALSE]
data_v_new <- cbind(clinical_matched, data_v_matched)


mypheno_v <- data_v_new[, c("chip", "Mat_Ethnicity", "Pat_Ethnicity",
                            "Maternal_Age", "Paternal_Age", "Gravidity", "Parity", "OBS_NORMAL_con")]
mypheno_v$Mat_Ethnicity <- factor(mypheno_v$Mat_Ethnicity)
mypheno_v$Pat_Ethnicity <- factor(mypheno_v$Pat_Ethnicity)
mypheno_v$Gravidity <- factor(mypheno_v$Gravidity)
mypheno_v$Parity <- factor(mypheno_v$Parity)

sva_data_format_v <- data.matrix(t(data_v_new[-c(1:10)]))
modcombat_v <- model.matrix(~1, data = mypheno_v)
batch_v <- mypheno_v$chip

combat_edata_v <- ComBat(dat = sva_data_format_v, batch = batch_v,
                         mod = modcombat_v, par.prior = TRUE, prior.plots = FALSE)


data_v_new_combat <- data.frame(
  Batch = data_v_new$chip,
  ID = data_v_new$ID,
  Mat_Ethnicity = data_v_new$Mat_Ethnicity,
  Pat_Ethnicity = data_v_new$Pat_Ethnicity,
  Maternal_Age = data_v_new$Maternal_Age,
  Paternal_Age = data_v_new$Paternal_Age,
  Gravidity = data_v_new$Gravidity,
  Parity = data_v_new$Parity,
  OBS_NORMAL_con = data_v_new$OBS_NORMAL_con,
  t(combat_edata_v),
  check.names = FALSE
)

metab_common_v <- data_v_new_combat[, -(1:9)]

Ftab_v <- matrix(NA, ncol = ncol(metab_common_v), nrow = 8)
rownames(Ftab_v) <- c('chip', 'Mat_Ethnicity', 'Pat_Ethnicity', 'Maternal_Age', 
                      'Paternal_Age', 'Gravidity', 'Parity', 'OBS_NORMAL_con')



# ANOVA
for(i in 1:ncol(metab_common_v)) {
  fit <- lm(metab_common_v[, i] ~ factor(Batch) + factor(Mat_Ethnicity) + factor(Pat_Ethnicity) +
              Maternal_Age + Paternal_Age + factor(Gravidity) + factor(Parity) + OBS_NORMAL_con,
            data = data_v_new_combat, na.action = na.omit)
  
  aovfit <- car::Anova(fit, type = 3, singular.ok = TRUE)
  Fvals <- aovfit$`F value`
  Fvals <- Fvals[2:(length(Fvals)-1)]
  Ftab_v[, i] <- Fvals
}

Ftab_v
F_stats <- rowMeans(Ftab_v, na.rm = TRUE)
result_v <- data.frame(
  Confounder = rownames(Ftab_v),
  FStatistic = F_stats
)

result_v$Confounder[result_v$Confounder == "chip"] <- "Batch"

write.csv(data_v_new_combat, "../output/knn_vsn_aftercombat_meta.csv", row.names = TRUE)
write.csv(result_v, "../output/SOV.csv", row.names = FALSE)




# SOV

png("../Figures/Fig1_SOV.png", width = 8, height = 6, units = "in", res = 600)

ggplot(result_v, aes(x = reorder(Confounder, -FStatistic), y = FStatistic, fill = Confounder)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ggtitle("SOV for VSN normalization after combat") +
  xlab("Potential Confounder") +
  ylab("F Statistic") +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()