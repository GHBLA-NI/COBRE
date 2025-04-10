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


generate_ethnicity_interaction_plot <- function(metabolite, ethnicity = c("mat", "pat", "both"), data, output_dir) {
  ethnicity <- match.arg(ethnicity)
  if (ethnicity == "mat") {
    form <- as.formula(paste0("OBS_NORMAL_con ~ ", metabolite, "*Mat_Ethnicity"))
  } else if (ethnicity == "pat") {
    form <- as.formula(paste0("OBS_NORMAL_con ~ ", metabolite, "*Pat_Ethnicity"))
  } else if (ethnicity == "both") {
    form <- as.formula(paste0("OBS_NORMAL_con ~ ", metabolite, "*Mat_Ethnicity + ", metabolite, "*Pat_Ethnicity"))
  }

  model <- lm(form, data = data)
  effect_plot <- allEffects(model)

  if (ethnicity == "mat") {
    effect_data <- as.data.frame(effect_plot[[paste0(metabolite, ":Mat_Ethnicity")]])
    aes_group <- "Mat_Ethnicity"
    ethnicity_label <- "Maternal Ethnicity"
    filename_suffix <- "mat"
  } else if (ethnicity == "pat") {
    effect_data <- as.data.frame(effect_plot[[paste0(metabolite, ":Pat_Ethnicity")]])
    aes_group <- "Pat_Ethnicity"
    ethnicity_label <- "Paternal Ethnicity"
    filename_suffix <- "pat"
  } else {
    generate_ethnicity_interaction_plot(metabolite, "mat", data, output_dir)
    generate_ethnicity_interaction_plot(metabolite, "pat", data, output_dir)
    return(invisible(NULL))
  }
  
  pg_plot <- ggplot(effect_data, aes_string(x = metabolite, y = "fit", color = aes_group, fill = aes_group, group = aes_group)) +
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
      title = paste("Effect of", ethnicity_label, "and", metabolite, "level on BMI"),
      x = paste(metabolite, "Level"),
      y = "BMI value",
      color = "Ethnicity"
    ) +
    scale_color_discrete(labels = c("Caucasian", "Asian", "NHPI")) +
    scale_fill_discrete(labels = c("Caucasian", "Asian", "NHPI"))
  
  print(pg_plot)
  
  output_file <- file.path(output_dir, paste0("interaction_", tolower(metabolite), "_", filename_suffix, ".png"))
  ggsave(filename = output_file, plot = pg_plot, width = 10, height = 7, dpi = 600)
  
  message("Plot saved to: ", output_file)
}


generate_ethnicity_interaction_plot(
  metabolite = "PC_ae_C44_6",
  ethnicity = "mat",
  data = meta_info,
  output_dir = "../Figures/"
)



generate_ethnicity_interaction_plot(
  metabolite = "LThreonine",
  ethnicity = "mat",
  data = meta_info,
  output_dir = "../Figures/"
)