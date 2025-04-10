library(ggplot2)
library(dplyr)
library(ggrepel)

snps<-read.csv("../output/slected_snps.csv",row.names = 1)
pvalue<-read.csv("../output/adjp_val_multinomial_momage.csv",row.names = 1)

snps$chr <- ifelse(snps$pos >= 186798081 & snps$pos <= 186958113, 1,
                   ifelse(snps$pos >= 71069743  & snps$pos <= 140999928, 9,
                          ifelse(snps$pos >= 23285775  & snps$pos <= 30386399, 15,
                                 ifelse(snps$pos >= 102839    & snps$pos <= 28327676, 16, NA))))

new_names <- substr(rownames(pvalue), 1, nchar(rownames(pvalue)) - 2)
new_names <- gsub("\\.", "-", new_names)
rownames(pvalue) <- new_names

snps$p <- pvalue$White[match(snps$snp, rownames(pvalue))]
names(snps)[names(snps) == "pos"] <- "bp"
names(snps)[names(snps) == "chr"] <- "CHR"
names(snps)[names(snps) == "bp"]  <- "BP"
names(snps)[names(snps) == "snp"] <- "SNP"
names(snps)[names(snps) == "p"]   <- "P"





allsnps<-read.csv("../output/all_snps_adjp_val_multinomial_momage.csv",row.names = 1)

# For PC: BP between 186798081 and 186958113 OR BP between 23285775 and 30386399
PC <- subset(allsnps, (BP >= 186798081 & BP <= 186958113) | 
               (BP >= 23285775  & BP <= 30386399))

# For LT: BP between 71069743 and 140999928 OR BP between 102839 and 28327676
LT <- subset(allsnps, (BP >= 71069743  & BP <= 140999928) | 
               (BP >= 102839    & BP <= 28327676))



# --- PREPARE THE DATA (using PC) ---
PC$CHR <- as.numeric(ifelse(PC$CHR == "X", "23", PC$CHR))
PC <- na.omit(PC)
don <- PC %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  arrange(CHR) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(CHR, tot) %>%
  left_join(PC, by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(
    BPcum = BP + tot,
    logp = -log10(P)
  ) %>%
  ungroup()


threshold <- -log10(0.05)  
non_sig_odd  <- don %>% filter(logp < threshold & CHR %% 2 == 1)
non_sig_even <- don %>% filter(logp < threshold & CHR %% 2 == 0)

sig_points <- don %>% filter(logp >= threshold)

chr_colors <- setNames(rainbow(23), as.character(1:23))
x_breaks <- seq(min(don$BPcum), max(don$BPcum), length.out = 23)
x_labels <- c(as.character(1:22), "X")

# --- CREATE THE PLOT ---
p <- ggplot() +
  geom_point(data = non_sig_odd, aes(x = BPcum, y = logp),
             color = "grey90", size = 1.3) +
  geom_point(data = non_sig_even, aes(x = BPcum, y = logp),
             color = "grey40", size = 1.3) +
  geom_point(data = sig_points, aes(x = BPcum, y = logp, color = factor(CHR)),
             size = 2) +
  geom_label_repel(data = sig_points,
                   aes(x = BPcum, y = logp, label = SNP, color = factor(CHR)),max.overlaps = 10,
                   size = 2, show.legend = FALSE) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_color_manual(values = chr_colors) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
  labs(x = "Chromosome", y = expression(-log[10](P)),
       title = "Manhattan Plot (PC(O-44:6))") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = "../Figures/Fig6C_ManhattanPC.png",plot = p, width = 12, height = 8, dpi = 300)






# --- PREPARE THE DATA (using LT) ---
LT$CHR <- as.numeric(ifelse(LT$CHR == "X", "23", LT$CHR))
LT <- na.omit(LT)

don <- LT %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  arrange(CHR) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(CHR, tot) %>%
  left_join(LT, by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(
    BPcum = BP + tot,
    logp = -log10(P)
  ) %>%
  ungroup()

threshold <- -log10(0.05)

non_sig_odd  <- don %>% filter(logp < threshold & CHR %% 2 == 1)
non_sig_even <- don %>% filter(logp < threshold & CHR %% 2 == 0)
sig_points   <- don %>% filter(logp >= threshold)

chr_colors <- setNames(rainbow(23), as.character(1:23))

x_breaks <- seq(min(don$BPcum), max(don$BPcum), length.out = 23)
x_labels <- c(as.character(1:22), "X")

q <- ggplot() +
  geom_point(data = non_sig_odd, aes(x = BPcum, y = logp),
             color = "grey90", size = 1.3) +
  geom_point(data = non_sig_even, aes(x = BPcum, y = logp),
             color = "grey40", size = 1.3) +
  geom_point(data = sig_points, aes(x = BPcum, y = logp, color = factor(CHR)),
             size = 2) +
  geom_label_repel(data = sig_points,
                   aes(x = BPcum, y = logp, label = SNP, color = factor(CHR)),max.overlaps = 50,
                   size = 2, show.legend = FALSE) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_color_manual(values = chr_colors) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
  labs(x = "Chromosome", y = expression(-log[10](P)),
       title = "Manhattan Plot (L-Threonine)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = "../Figures/Fig6D_ManhattanLT.png",
       plot = q, width = 12, height = 8, dpi = 300)
