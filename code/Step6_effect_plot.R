library(Hmisc)
library(corrplot)
library(dummies)
library(ggplot2)
library(effects)
library(dplyr)
library(car)

vsn_rm_all = read.csv('../output/knn_vsn_aftercombat_meta.csv', row.names=1)
meta_info <- vsn_rm_all[, c(3, 4, 9:ncol(vsn_rm_all))]
meta_info$Mat_Ethnicity <- factor(meta_info$Mat_Ethnicity)
meta_info$Pat_Ethnicity <- factor(meta_info$Pat_Ethnicity)



generate_ethnicity_interaction_plot <- function(metabolite, data, output_file = NULL) {
  form <- as.formula(paste0("OBS_NORMAL_con ~ ", metabolite, "*Mat_Ethnicity"))
  model <- lm(form, data = data)
  effect_plot <- allEffects(model)
  effect_data <- as.data.frame(effect_plot[[paste0(metabolite, ":Mat_Ethnicity")]])
  
  pg_plot <- ggplot(effect_data, aes_string(x = metabolite, y = "fit", 
                                            color = "Mat_Ethnicity", fill = "Mat_Ethnicity", group = "Mat_Ethnicity")) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = NA, show.legend = FALSE) +
    theme(
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      panel.background = element_rect(fill = 'white'),
      plot.background = element_rect(fill = 'white', color = NA),
      legend.background = element_rect(fill = 'white'),
      legend.box.background = element_rect(fill = 'white'),
      panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"),
      panel.grid.minor = element_line(size = 0.10, linetype = 'solid', colour = "gray")
    ) +
    labs(
      title = paste("Effect of Maternal Ethnicity and", metabolite, "level on BMI"),
      x = paste(metabolite, "Level"),
      y = "BMI value",
      color = "Ethnicity"
    ) +
    scale_color_discrete(labels = c("Caucasian", "Asian", "NHPI")) +
    scale_fill_discrete(labels = c("Caucasian", "Asian", "NHPI"))
  
  
  if (is.null(output_file)) {
    output_file <- paste0(tolower(metabolite), ".png")
  }
  
  ggsave(filename = output_file, plot = pg_plot, width = 10, height = 7, dpi = 600)
  message("Plot saved to: ", output_file)
}






generate_ethnicity_interaction_plot(
  metabolite = "PC_ae_C44_6",
  data = meta_info,
  output_file = "../Figures/Fig6B_effect_PC.png"
)



generate_ethnicity_interaction_plot(
  metabolite = "LThreonine",
  data = meta_info,
  output_file = "../Figures/Fig6A_effect_LT.png"
)