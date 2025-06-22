library(dplyr)
library(ggplot2)
library(ggpubr)

metab_info <- read.csv("../output/knn_vsn_aftercombat_meta.csv", row.names = 1)
clinical    <- readRDS("../data/complete_ci.rds")

clinical$Mat_Ethnicity <- factor(clinical$Mat_Ethnicity,
                                  levels = c("White", "Asian", "NHPI"))

# Preprocessing
metab_info <- metab_info[, -(1:10)]
metab_info <- t(metab_info)

metab_eth_all <- merge(
  t(metab_info),
  clinical[, c("ID", "label", "Mat_Ethnicity")],
  by.x = "row.names",
  by.y = "ID"
)
rownames(metab_eth_all) <- metab_eth_all$Row.names
metab_eth_all <- metab_eth_all[, -1]
colnames(metab_eth_all)[ncol(metab_eth_all)] <- "Mat_Ethnicity.y"

#  Wilcoxon p-values and Plots 
metabolite_list <- c("PC_ae_C44_6", "X15Anhydrosorbitol", "Glycine", "LThreonine")
results_list <- list()

for (metabolite in metabolite_list) {
  tmp <- metab_eth_all %>%
    select(all_of(metabolite), Mat_Ethnicity.y, label) %>%
    rename(Value = !!metabolite) %>%
    filter(Mat_Ethnicity.y %in% c("Asian", "White")) %>%
    filter(!is.na(Value))
  
  v1 <- tmp$Value[tmp$Mat_Ethnicity.y == "Asian"]
  v2 <- tmp$Value[tmp$Mat_Ethnicity.y == "White"]
  wtest <- wilcox.test(v1, v2)
  
  results_list[[metabolite]] <- tibble(
    Metabolite = metabolite,
    Comparison = "Asian_vs_White",
    wilcox_pvalue = signif(wtest$p.value, 4)
  )
  
  stats_tbl <- tmp %>%
    group_by(Mat_Ethnicity.y, label) %>%
    summarise(
      mean_val = mean(Value, na.rm = TRUE),
      sem_val  = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  all_neg <- all(stats_tbl$mean_val < 0)
  if (all_neg) {
    stats_tbl <- stats_tbl %>%
      mutate(plot_mean = abs(mean_val), plot_sem = sem_val)
    axis_transform <- scale_y_continuous(labels = function(x) -x)
  } else {
    stats_tbl <- stats_tbl %>%
      mutate(plot_mean = mean_val, plot_sem = sem_val)
    axis_transform <- NULL
  }
  
  p <- ggplot(stats_tbl, aes(
    x = Mat_Ethnicity.y, y = plot_mean, fill = label
  )) +
    geom_bar(stat = "identity",
             position = position_dodge(0.7),
             width = 0.6,
             color = "black") +
    geom_errorbar(
      aes(ymin = plot_mean - plot_sem, ymax = plot_mean + plot_sem),
      position = position_dodge(0.7),
      width = 0.2
    ) +
    scale_fill_manual(values = c("Normal" = "#D73027", "Obese" = "#4575B4")) +
    labs(
      title = metabolite,
      x = NULL,
      y = paste0("Mean ", metabolite, " (± SEM)"),
      fill = "Label"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 14)
    )
  
  if (!is.null(axis_transform)) {
    p <- p + axis_transform
  }
  
  p <- p + stat_compare_means(
    comparisons   = list(c("Asian", "White")),
    method        = "wilcox.test",
    label         = "p.signif",
    size          = 6,
    linewidth     = 1.2,
    tip.length    = 0.02,
    step.increase = 0.1,
    hide.ns       = TRUE
  )
  

  out_dir <- "../Figures"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ggsave(file.path(out_dir, paste0("Figure5_", metabolite, ".png")),
         plot = p, width = 6, height = 5, dpi = 600, bg = "white")
}


pvalue_table <- bind_rows(results_list)

pvalue_table