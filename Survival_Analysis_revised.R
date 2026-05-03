#########################################################################
# Overall Survival Status, Overall Survival (months)
# K-M curve and Log Rank Test

# load the required packages
library(readr)
library(survival)
library(survminer)
library(dplyr)
library(survMisc)
library(patchwork)
library(ggplot2)

# read in the data
all_clinical <- read_tsv("tcga_clinical/tcga_all_clinical_new.tsv", quote = "", col_names = TRUE,
                         name_repair = "minimal")
TP53_mutation <- read_tsv("tcga_clinical/TP53_mutations_new.tsv", quote = "", col_names = TRUE,
                          name_repair = "minimal")
CTNNB1_mutation <- read_tsv("tcga_clinical/CTNNB1_mutations_new.tsv", quote = "", col_names = TRUE,
                            name_repair = "minimal")
JAK1_mutation <- read_tsv("tcga_clinical/JAK1_mutations_new.tsv", quote = "", col_names = TRUE,
                          name_repair = "minimal")
P4HA1_mutation <- read_tsv("tcga_clinical/P4HA1_mutations_new.tsv", quote = "", col_names = TRUE,
                           name_repair = "minimal")
TLE1_mutation <- read_tsv("tcga_clinical/TLE1_mutations_new.tsv", quote = "", col_names = TRUE,
                          name_repair = "minimal")

# extract all "Patient ID" lists from each of the gene mutation datasets
TP53_mutation_patients <- TP53_mutation$`Patient ID`
CTNNB1_mutation_patients <- CTNNB1_mutation$`Patient ID`
JAK1_mutation_patients <- JAK1_mutation$`Patient ID`
P4HA1_mutation_patients <- P4HA1_mutation$`Patient ID`
TLE1_mutation_patients <- TLE1_mutation$`Patient ID`

# create a new mutation indicator column in the complete clinical dataset 
# for each type of gene mutation
all_clinical <- all_clinical %>%
  mutate(TP53_mutation = ifelse(`Patient ID` %in% TP53_mutation_patients, "Yes", "No"),
         CTNNB1_mutation = ifelse(`Patient ID` %in% CTNNB1_mutation_patients, 
                                  "Yes", "No"),
         JAK1_mutation = ifelse(`Patient ID` %in% JAK1_mutation_patients, "Yes", "No"),
         P4HA1_mutation = ifelse(`Patient ID` %in% P4HA1_mutation_patients, "Yes", "No"),
         TLE1_mutation = ifelse(`Patient ID` %in% TLE1_mutation_patients, "Yes", "No"))

# create a new column for the overall survival status
all_clinical <- all_clinical %>%
  mutate(OS_status = ifelse(`Overall Survival Status` == "1:DECEASED", 1, 0))
# create a new column for the overall survival time in months
all_clinical <- all_clinical %>%
  mutate(OS_time = `Overall Survival (Months)`)

# plot the K-M curve and conduct the log rank test for each type of gene mutation
# TP53 mutation 
km_fit_TP53 <- survfit(Surv(OS_time, OS_status) ~ TP53_mutation, data = all_clinical)
ggsurvplot(km_fit_TP53, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for TP53 Mutation",
           xlab = "Overall Survival Time (Months)",
           ylab = "Survival Probability")
log_rank_TP53 <- survdiff(Surv(OS_time, OS_status) ~ TP53_mutation, data = all_clinical)
log_rank_TP53_pvalue <- 1 - pchisq(log_rank_TP53$chisq, 
                                   df = length(log_rank_TP53$n) - 1)
# CTNNB1 mutation
km_fit_CTN <- survfit(Surv(OS_time, OS_status) ~ CTNNB1_mutation, data = all_clinical)
ggsurvplot(km_fit_CTN, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for CTNNB1 Mutation",
           xlab = "Overall Survival Time (Months)",
           ylab = "Survival Probability")
log_rank_CTN <- survdiff(Surv(OS_time, OS_status) ~ CTNNB1_mutation, data = all_clinical)
log_rank_CTN_pvalue <- 1 - pchisq(log_rank_CTN$chisq, df = length(log_rank_CTN$n) - 1)
# JAK1 mutation
km_fit_JAK1 <- survfit(Surv(OS_time, OS_status) ~ JAK1_mutation, data = all_clinical)
ggsurvplot(km_fit_JAK1, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for JAK1 Mutation",
           xlab = "Overall Survival Time (Months)",
           ylab = "Survival Probability")
log_rank_JAK1 <- survdiff(Surv(OS_time, OS_status) ~ JAK1_mutation, data = all_clinical)
log_rank_JAK1_pvalue <- 1 - pchisq(log_rank_JAK1$chisq, df = length(log_rank_JAK1$n) - 1)
# P4HA1 mutation
km_fit_P4HA1 <- survfit(Surv(OS_time, OS_status) ~ P4HA1_mutation, data = all_clinical)
ggsurvplot(km_fit_P4HA1, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for P4HA1 Mutation",
           xlab = "Overall Survival Time (Months)",
           ylab = "Survival Probability")
log_rank_P4HA1 <- survdiff(Surv(OS_time, OS_status) ~ P4HA1_mutation, 
                           data = all_clinical)
log_rank_P4HA1_pvalue <- 1 - pchisq(log_rank_P4HA1$chisq, 
                                    df = length(log_rank_P4HA1$n) - 1)
# TLE1 mutation
km_fit_TLE1 <- survfit(Surv(OS_time, OS_status) ~ TLE1_mutation, data = all_clinical)
ggsurvplot(km_fit_TLE1, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for TLE1 Mutation",
           xlab = "Overall Survival Time (Months)",
           ylab = "Survival Probability")
log_rank_TLE1 <- survdiff(Surv(OS_time, OS_status) ~ TLE1_mutation, data = all_clinical)
log_rank_TLE1_pvalue <- 1 - pchisq(log_rank_TLE1$chisq, df = length(log_rank_TLE1$n) - 1)

# create a summary table for the log rank test results
log_rank_results <- data.frame(
  Gene_Mutation = c("TP53", "CTNNB1", "JAK1", "P4HA1", "TLE1"),
  p_value = c(log_rank_TP53_pvalue, log_rank_CTN_pvalue, log_rank_JAK1_pvalue,
              log_rank_P4HA1_pvalue, log_rank_TLE1_pvalue)
)
# print the summary table
print(log_rank_results)


#################################################################################
# Disease Free Status, Disease Free (Months)
# K-M curve and Log Rank Test

# create a new column for the disease free status
all_clinical <- all_clinical %>%
  mutate(DFS_status = ifelse(`Disease Free Status` == "1:Recurred/Progressed", 1, 0))
# create a new column for the disease free time in months
all_clinical <- all_clinical %>%
  mutate(DFS_time = `Disease Free (Months)`)

# plot the K-M curve and conduct the log rank test for each type of gene mutation
# TP53 mutation
km_fit_TP53 <- survfit(Surv(DFS_time, DFS_status) ~ TP53_mutation, data = all_clinical)
ggsurvplot(km_fit_TP53, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for TP53 Mutation",
           xlab = "Disease Free Time (Months)",
           ylab = "Survival Probability")
log_rank_TP53 <- survdiff(Surv(DFS_time, DFS_status) ~ TP53_mutation, 
                          data = all_clinical)
log_rank_TP53_pvalue <- 1 - pchisq(log_rank_TP53$chisq, 
                                   df = length(log_rank_TP53$n) - 1)
# CTNNB1 mutation
km_fit_CTN <- survfit(Surv(DFS_time, DFS_status) ~ CTNNB1_mutation, data = all_clinical)
ggsurvplot(km_fit_CTN, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for CTNNB1 Mutation",
           xlab = "Disease Free Time (Months)",
           ylab = "Survival Probability")
log_rank_CTN <- survdiff(Surv(DFS_time, DFS_status) ~ CTNNB1_mutation, 
                         data = all_clinical)
log_rank_CTN_pvalue <- 1 - pchisq(log_rank_CTN$chisq, 
                                  df = length(log_rank_CTN$n) - 1)
# JAK1 mutation
km_fit_JAK1 <- survfit(Surv(DFS_time, DFS_status) ~ JAK1_mutation, data = all_clinical)
ggsurvplot(km_fit_JAK1, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for JAK1 Mutation",
           xlab = "Disease Free Time (Months)",
           ylab = "Survival Probability")
log_rank_JAK1 <- survdiff(Surv(DFS_time, DFS_status) ~ JAK1_mutation, 
                          data = all_clinical)
log_rank_JAK1_pvalue <- 1 - pchisq(log_rank_JAK1$chisq, 
                                   df = length(log_rank_JAK1$n) - 1)
# P4HA1 mutation
km_fit_P4HA1 <- survfit(Surv(DFS_time, DFS_status) ~ P4HA1_mutation, 
                        data = all_clinical)
ggsurvplot(km_fit_P4HA1, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for P4HA1 Mutation",
           xlab = "Disease Free Time (Months)",
           ylab = "Survival Probability")
log_rank_P4HA1 <- survdiff(Surv(DFS_time, DFS_status) ~ P4HA1_mutation, 
                           data = all_clinical)
log_rank_P4HA1_pvalue <- 1 - pchisq(log_rank_P4HA1$chisq, 
                                    df = length(log_rank_P4HA1$n) - 1)
# TLE1 mutation
km_fit_TLE1 <- survfit(Surv(DFS_time, DFS_status) ~ TLE1_mutation, data = all_clinical)
ggsurvplot(km_fit_TLE1, data = all_clinical, pval = TRUE, conf.int = TRUE,
           title = "K-M Curve for TLE1 Mutation",
           xlab = "Disease Free Time (Months)",
           ylab = "Survival Probability")
log_rank_TLE1 <- survdiff(Surv(DFS_time, DFS_status) ~ TLE1_mutation, 
                          data = all_clinical)
log_rank_TLE1_pvalue <- 1 - pchisq(log_rank_TLE1$chisq, 
                                   df = length(log_rank_TLE1$n) - 1)
# create a summary table for the log rank test results
log_rank_results <- data.frame(
  Gene_Mutation = c("TP53", "CTNNB1", "JAK1", "P4HA1", "TLE1"),
  p_value = c(log_rank_TP53_pvalue, log_rank_CTN_pvalue, log_rank_JAK1_pvalue,
              log_rank_P4HA1_pvalue, log_rank_TLE1_pvalue))
# print the summary table
print(log_rank_results)



#########################################################################
# UPDATED Figure 5: K-M Plots of TP53 Mutation for Overall Survival and Disease-Free Survival
# K-M plot of TP53 Mutation with customized themes

# make consistent mutation labeling
all_clinical$TP53_mutation <- factor(
  all_clinical$TP53_mutation,
  levels = c("No", "Yes"),
  labels = c("Wild-type", "Mutated")
)

base_font <- 16   # ONE font size for all text elements

# ---------- Overall Survival ----------
fit_os <- survfit(Surv(OS_time, OS_status) ~ TP53_mutation, data = all_clinical)
lr_os  <- survdiff(Surv(OS_time, OS_status) ~ TP53_mutation, data = all_clinical)
p_os   <- 1 - pchisq(lr_os$chisq, df = length(lr_os$n) - 1)

g_os <- ggsurvplot(
  fit_os, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  # no panel title
  xlab = "Overall survival (months)",
  ylab = "Survival probability",
  legend.title = expression(italic("TP53") * " mutation"),
  legend.labs  = levels(all_clinical$TP53_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = base_font)
)

# put legend in upper-right whitespace & enforce Arial + single size
g_os$plot <- g_os$plot +
  theme(
    text = element_text(family = "Arial", size = base_font),
    axis.title = element_text(size = base_font),
    axis.text = element_text(size = base_font),
    legend.title = element_text(size = base_font),
    legend.text = element_text(size = base_font),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 20)
  )

# italic P in lower-left
p_os_label <- paste0("italic(P) == ", formatC(p_os, format = "f", digits = 3))
g_os$plot <- g_os$plot +
  annotate(
    "text",
    x = 5, y = 0.12,
    label = p_os_label,
    parse = TRUE,
    hjust = 0,
    family = "Arial",
    size = 5.6    # ~16 pt; matches other text
  )

# ---------- Disease-Free Survival ----------
fit_dfs <- survfit(Surv(DFS_time, DFS_status) ~ TP53_mutation, data = all_clinical)
lr_dfs  <- survdiff(Surv(DFS_time, DFS_status) ~ TP53_mutation, data = all_clinical)
p_dfs   <- 1 - pchisq(lr_dfs$chisq, df = length(lr_dfs$n) - 1)

g_dfs <- ggsurvplot(
  fit_dfs, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  xlab = "Disease-free survival (months)",
  ylab = "Survival probability",
  legend.title = expression(italic("TP53") * " mutation"),
  legend.labs  = levels(all_clinical$TP53_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = base_font)
)

g_dfs$plot <- g_dfs$plot +
  theme(
    text = element_text(family = "Arial", size = base_font),
    axis.title = element_text(size = base_font),
    axis.text = element_text(size = base_font),
    legend.title = element_text(size = base_font),
    legend.text = element_text(size = base_font),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 20)
  )

p_dfs_label <- paste0("italic(P) == ", formatC(p_dfs, format = "f", digits = 4))
g_dfs$plot <- g_dfs$plot +
  annotate(
    "text",
    x = 5, y = 0.12,
    label = p_dfs_label,
    parse = TRUE,
    hjust = 0,
    family = "Arial",
    size = 5.6    # ~16 pt; matches other text
  )

# ---------- Combine panels with labels A/B ----------
fig5 <- g_os$plot + g_dfs$plot +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(family = "Arial", size = base_font),
    plot.tag = element_text(family = "Arial", face = "bold", size = 18),
    # top pixel & left pixel alignment (approximate)
    plot.tag.position = c(0, 0.982)
  )

# show
fig5

# save
# ggsave("Fig5_TP53_OS_DFS_revised.png", fig5, width = 12, height = 6, dpi = 300)


#########################################################################
# UPDATED Figure Supplementary A2: K-M Plots of CTNNB1, JAK1, and P4HA1 Mutations 
# Make each mutation column use consistent labels
all_clinical$CTNNB1_mutation <- factor(all_clinical$CTNNB1_mutation,
                                       levels = c("No","Yes"),
                                       labels = c("Wild-type","Mutated"))
all_clinical$JAK1_mutation <- factor(all_clinical$JAK1_mutation,
                                     levels = c("No","Yes"),
                                     labels = c("Wild-type","Mutated"))
all_clinical$P4HA1_mutation <- factor(all_clinical$P4HA1_mutation,
                                      levels = c("No","Yes"),
                                      labels = c("Wild-type","Mutated"))

# single base font size for ALL text elements
base_font <- 16
p_text_size <- base_font / 2.845   # convert pt -> ggplot size
legend_pos <- c(0.95, 0.95)

# ============== CTNNB1 (OS) ==============
fit_os_CTNNB1 <- survfit(Surv(OS_time, OS_status) ~ CTNNB1_mutation, data = all_clinical)
lr_os_CTNNB1  <- survdiff(Surv(OS_time, OS_status) ~ CTNNB1_mutation, data = all_clinical)
p_os_CTNNB1   <- 1 - pchisq(lr_os_CTNNB1$chisq, df = length(lr_os_CTNNB1$n) - 1)

g1 <- ggsurvplot(
  fit_os_CTNNB1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = FALSE,   # will add p manually
  xlab = "Overall survival (months)",
  ylab = "Survival probability",
  legend.title = expression(italic("CTNNB1") * " mutation"),
  legend.labs = levels(all_clinical$CTNNB1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = base_font)
)

p1 <- g1$plot +
  theme(
    text = element_text(family = "Arial", size = base_font),
    axis.title = element_text(size = base_font),
    axis.text = element_text(size = base_font),
    legend.title = element_text(size = base_font),
    legend.text = element_text(size = base_font),
    legend.position = legend_pos,
    legend.justification = c("right","top"),
    legend.background = element_blank(),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 20)
  ) +
  annotate(
    "text",
    x = 0, y = 0.09,
    label = paste0("italic(P) == ", formatC(p_os_CTNNB1, format = "f", digits = 3)),
    parse = TRUE,
    hjust = 0,
    family = "Arial",
    size = p_text_size
  )

# ============== CTNNB1 (DFS) ==============
fit_dfs_CTNNB1 <- survfit(Surv(DFS_time, DFS_status) ~ CTNNB1_mutation, data = all_clinical)
lr_dfs_CTNNB1  <- survdiff(Surv(DFS_time, DFS_status) ~ CTNNB1_mutation, data = all_clinical)
p_dfs_CTNNB1   <- 1 - pchisq(lr_dfs_CTNNB1$chisq, df = length(lr_dfs_CTNNB1$n) - 1)

g2 <- ggsurvplot(
  fit_dfs_CTNNB1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = FALSE,
  xlab = "Disease-free survival (months)",
  ylab = "Survival probability",
  legend.title = expression(italic("CTNNB1") * " mutation"),
  legend.labs = levels(all_clinical$CTNNB1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = base_font)
)

p2 <- g2$plot +
  theme(
    text = element_text(family = "Arial", size = base_font),
    axis.title = element_text(size = base_font),
    axis.text = element_text(size = base_font),
    legend.title = element_text(size = base_font),
    legend.text = element_text(size = base_font),
    legend.position = legend_pos,
    legend.justification = c("right","top"),
    legend.background = element_blank(),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 20)
  ) +
  annotate(
    "text",
    x = 0, y = 0.09,
    label = paste0("italic(P) == ", formatC(p_dfs_CTNNB1, format = "f", digits = 4)),
    parse = TRUE,
    hjust = 0,
    family = "Arial",
    size = p_text_size
  )

# ============== JAK1 (OS) ==============
fit_os_JAK1 <- survfit(Surv(OS_time, OS_status) ~ JAK1_mutation, data = all_clinical)
lr_os_JAK1  <- survdiff(Surv(OS_time, OS_status) ~ JAK1_mutation, data = all_clinical)
p_os_JAK1   <- 1 - pchisq(lr_os_JAK1$chisq, df = length(lr_os_JAK1$n) - 1)

g3 <- ggsurvplot(
  fit_os_JAK1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = FALSE,
  xlab = "Overall survival (months)",
  ylab = "Survival probability",
  legend.title = expression(italic("JAK1") * " mutation"),
  legend.labs = levels(all_clinical$JAK1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = base_font)
)

p3 <- g3$plot +
  theme(
    text = element_text(family = "Arial", size = base_font),
    axis.title = element_text(size = base_font),
    axis.text = element_text(size = base_font),
    legend.title = element_text(size = base_font),
    legend.text = element_text(size = base_font),
    legend.position = legend_pos,
    legend.justification = c("right","top"),
    legend.background = element_blank(),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 20)
  ) +
  annotate(
    "text",
    x = 0, y = 0.09,
    label = paste0("italic(P) == ", formatC(p_os_JAK1, format = "f", digits = 3)),
    parse = TRUE,
    hjust = 0,
    family = "Arial",
    size = p_text_size
  )

# ============== JAK1 (DFS) ==============
fit_dfs_JAK1 <- survfit(Surv(DFS_time, DFS_status) ~ JAK1_mutation, data = all_clinical)
lr_dfs_JAK1  <- survdiff(Surv(DFS_time, DFS_status) ~ JAK1_mutation, data = all_clinical)
p_dfs_JAK1   <- 1 - pchisq(lr_dfs_JAK1$chisq, df = length(lr_dfs_JAK1$n) - 1)

g4 <- ggsurvplot(
  fit_dfs_JAK1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = FALSE,
  xlab = "Disease-free survival (months)",
  ylab = "Survival probability",
  legend.title = expression(italic("JAK1") * " mutation"),
  legend.labs = levels(all_clinical$JAK1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = base_font)
)

p4 <- g4$plot +
  theme(
    text = element_text(family = "Arial", size = base_font),
    axis.title = element_text(size = base_font),
    axis.text = element_text(size = base_font),
    legend.title = element_text(size = base_font),
    legend.text = element_text(size = base_font),
    legend.position = legend_pos,
    legend.justification = c("right","top"),
    legend.background = element_blank(),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 20)
  ) +
  annotate(
    "text",
    x = 0, y = 0.09,
    label = paste0("italic(P) == ", formatC(p_dfs_JAK1, format = "f", digits = 4)),
    parse = TRUE,
    hjust = 0,
    family = "Arial",
    size = p_text_size
  )

# ============== P4HA1 (OS) ==============
fit_os_P4HA1 <- survfit(Surv(OS_time, OS_status) ~ P4HA1_mutation, data = all_clinical)
lr_os_P4HA1  <- survdiff(Surv(OS_time, OS_status) ~ P4HA1_mutation, data = all_clinical)
p_os_P4HA1   <- 1 - pchisq(lr_os_P4HA1$chisq, df = length(lr_os_P4HA1$n) - 1)

g5 <- ggsurvplot(
  fit_os_P4HA1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = FALSE,
  xlab = "Overall survival (months)",
  ylab = "Survival probability",
  legend.title = expression(italic("P4HA1") * " mutation"),
  legend.labs = levels(all_clinical$P4HA1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = base_font)
)

p5 <- g5$plot +
  theme(
    text = element_text(family = "Arial", size = base_font),
    axis.title = element_text(size = base_font),
    axis.text = element_text(size = base_font),
    legend.title = element_text(size = base_font),
    legend.text = element_text(size = base_font),
    legend.position = legend_pos,
    legend.justification = c("right","top"),
    legend.background = element_blank(),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 20)
  ) +
  annotate(
    "text",
    x = 0, y = 0.09,
    label = paste0("italic(P) == ", formatC(p_os_P4HA1, format = "f", digits = 3)),
    parse = TRUE,
    hjust = 0,
    family = "Arial",
    size = p_text_size
  )

# ============== P4HA1 (DFS) ==============
fit_dfs_P4HA1 <- survfit(Surv(DFS_time, DFS_status) ~ P4HA1_mutation, data = all_clinical)
lr_dfs_P4HA1  <- survdiff(Surv(DFS_time, DFS_status) ~ P4HA1_mutation, data = all_clinical)
p_dfs_P4HA1   <- 1 - pchisq(lr_dfs_P4HA1$chisq, df = length(lr_dfs_P4HA1$n) - 1)

g6 <- ggsurvplot(
  fit_dfs_P4HA1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = FALSE,
  xlab = "Disease-free survival (months)",
  ylab = "Survival probability",
  legend.title = expression(italic("P4HA1") * " mutation"),
  legend.labs = levels(all_clinical$P4HA1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = base_font)
)

p6 <- g6$plot +
  theme(
    text = element_text(family = "Arial", size = base_font),
    axis.title = element_text(size = base_font),
    axis.text = element_text(size = base_font),
    legend.title = element_text(size = base_font),
    legend.text = element_text(size = base_font),
    legend.position = legend_pos,
    legend.justification = c("right","top"),
    legend.background = element_blank(),
    plot.margin = margin(t = 15, r = 15, b = 10, l = 20)
  ) +
  annotate(
    "text",
    x = 0, y = 0.09,
    label = paste0("italic(P) == ", formatC(p_dfs_P4HA1, format = "f", digits = 4)),
    parse = TRUE,
    hjust = 0,
    family = "Arial",
    size = p_text_size
  )

# ---------- Arrange as 3 rows × 2 columns (OS | DFS) with tags A–F ----------
fig_3x2 <- (p1 | p2) /
           (p3 | p4) /
           (p5 | p6)

fig_3x2 <- fig_3x2 +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(family = "Arial", size = base_font),
    plot.tag = element_text(family = "Arial", face = "bold", size = 18),
    # A–F at top-left, aligned with top of image, shifted slightly left
    plot.tag.position = c(0, 0.982)
  )

# show
fig_3x2

# save
# ggsave("Figure_CTNNB1_JAK1_P4HA1_3x2_OS_DFS_revised.png",
       # fig_3x2, width = 15, height = 20, dpi = 300)



###############################################################################
# Multivariate Cox Proportional Hazards Regression
#
# Three models are fitted for each of OS and DFS:
#   Model 1 – KM log-rank univariate (Hazard Ratios extracted from univariate coxph)
#   Model 2 – Joint model: gene mutations only (no clinical covariates)
#   Model 3 – Fully adjusted model: gene mutations + clinical covariates
###############################################################################

# ── 0. Ensure gene mutation columns are factors, Wild-type as reference ───────
# Reference: standard practice in survival analysis; Bradburn et al. (2003)
# Br J Cancer 89(2):232-238 — reference category should be the unexposed group.
mutation_vars <- c("TP53_mutation", "CTNNB1_mutation", "JAK1_mutation",
                   "P4HA1_mutation")

for (v in mutation_vars) {
  all_clinical[[v]] <- factor(all_clinical[[v]],
                              levels = c("Wild-type", "Mutated"))
}

# ── 1. Prepare clinical covariates ────────────────────────────────────────────

# --- 1a. Diagnosis Age: continuous numeric -----------------------------------
# Reference: treated as a continuous covariate; each 1-year increase in age
# is modelled as a linear change in log-hazard (Kleinbaum & Klein, 2012,
# Survival Analysis: A Self-Learning Text, 3rd ed., Springer).
all_clinical$Diagnosis_Age <- as.numeric(all_clinical$`Diagnosis Age`)

# --- 1b. Sex: binary factor, Female as reference ----------------------------
# Reference: binary covariate coded as factor. Female set as reference
# because it is typically the lower-risk group in HCC
# (McGlynn et al., 2021, Nat Rev Gastroenterol Hepatol 18:613-629).
all_clinical$Sex <- factor(all_clinical$Sex,
                           levels = c("Female", "Male"))

# --- 1c. Ethnicity Category: binary factor -----------------------------------
# Reference: binary factor; "NOT HISPANIC OR LATINO" set as reference
# (the majority group in TCGA-LIHC; encoding follows TCGA data dictionary).
all_clinical$Ethnicity_Category <- factor(
  all_clinical$`Ethnicity Category`,
  levels = c("NOT HISPANIC OR LATINO", "HISPANIC OR LATINO")
)

# --- 1d. Race Category: nominal factor, WHITE as reference ------------------
# Reference: nominal factor with four levels. WHITE set as reference
# (largest subgroup in TCGA-LIHC). Rare categories (AMERICAN INDIAN OR
# ALASKA NATIVE) are retained; sparse cells may widen CIs.
# See Sondak et al. approach in TCGA pan-cancer papers.
all_clinical$Race_Category <- factor(
  all_clinical$`Race Category`,
  levels = c("WHITE", "ASIAN", "BLACK OR AFRICAN AMERICAN",
             "AMERICAN INDIAN OR ALASKA NATIVE")
)

# --- 1e. AJCC Pathologic Stage: ordinal factor, Stage I as reference ---------
# Reference: collapsed to main stages (I/II/III/IV) following AJCC 8th ed.
# Sub-stages (IIIA/IIIB/IIIC, IVA/IVB) are merged to their parent stage
# to avoid sparse-cell problems (Groome et al., 2001, CA Cancer J Clin).
# Stage I is the lowest-risk reference category.
stage_raw <- all_clinical$`AJCC Pathologic Stage`
stage_collapsed <- dplyr::case_when(
  stage_raw %in% c("Stage I")                        ~ "Stage I",
  stage_raw %in% c("Stage II")                       ~ "Stage II",
  stage_raw %in% c("Stage III", "Stage IIIA",
                   "Stage IIIB", "Stage IIIC")        ~ "Stage III",
  stage_raw %in% c("Stage IV", "Stage IVA",
                   "Stage IVB")                       ~ "Stage IV",
  TRUE                                                ~ NA_character_
)
all_clinical$AJCC_Stage <- factor(stage_collapsed,
                                  levels = c("Stage I", "Stage II",
                                             "Stage III", "Stage IV"))

# --- 1f. AJCC Pathologic T-Stage: ordinal factor, T1 as reference ------------
# Reference: sub-stages collapsed (T2a/T2b → T2; T3a/T3b → T3).
# T1 is the lowest-risk reference. TX (unknown) recoded to NA.
# Rationale: same sparse-cell logic as for overall stage above.
t_raw <- all_clinical$`AJCC Pathologic T-Stage`
t_collapsed <- dplyr::case_when(
  t_raw %in% c("T1")           ~ "T1",
  t_raw %in% c("T2", "T2a", "T2b") ~ "T2",
  t_raw %in% c("T3", "T3a", "T3b") ~ "T3",
  t_raw %in% c("T4")           ~ "T4",
  TRUE                         ~ NA_character_   # TX → NA
)
all_clinical$AJCC_T <- factor(t_collapsed,
                              levels = c("T1", "T2", "T3", "T4"))

# --- 1g. AJCC Pathologic N-Stage: factor, N0 as reference -------------------
# Reference: N0 (no nodal involvement) is the reference.
# NX (cannot be assessed) recoded to NA to avoid biasing the estimate.
n_raw <- all_clinical$`AJCC Pathologic N-Stage`
all_clinical$AJCC_N <- factor(
  dplyr::case_when(n_raw == "N0" ~ "N0",
                   n_raw == "N1" ~ "N1",
                   TRUE          ~ NA_character_),   # NX → NA
  levels = c("N0", "N1")
)

# --- 1h. AJCC Pathologic M-Stage: factor, M0 as reference -------------------
# Reference: M0 (no distant metastasis) is the reference.
# MX (cannot be assessed) recoded to NA.
m_raw <- all_clinical$`AJCC Pathologic M-Stage`
all_clinical$AJCC_M <- factor(
  dplyr::case_when(m_raw == "M0" ~ "M0",
                   m_raw == "M1" ~ "M1",
                   TRUE          ~ NA_character_),   # MX → NA
  levels = c("M0", "M1")
)

# --- 1i. Prior Malignancy: binary factor, No as reference -------------------
# Reference: binary variable from TCGA clinical data (TRUE/FALSE).
# Converted to factor with No as reference.
all_clinical$Prior_Malignancy <- factor(
  ifelse(all_clinical$`Prior Malignancy` == TRUE, "Yes", "No"),
  levels = c("No", "Yes")
)

# --- 1j. Prior Treatment: binary factor, No as reference --------------------
# Reference: same coding as Prior Malignancy above.
all_clinical$Prior_Treatment <- factor(
  ifelse(all_clinical$`Prior Treatment` == TRUE, "Yes", "No"),
  levels = c("No", "Yes")
)

# --- 1k. Fraction Genome Altered: continuous, log-transformed ----------------
# Reference: continuous covariate [0,1] representing the proportion of the
# genome with copy-number alterations. Log-transformation applied here
# since the column distribution is right-skewed.
# (Taylor et al., 2018, Cancer Cell 33:676-689 — TCGA pan-cancer CNA study)
all_clinical$Fraction_Genome_Altered_log <- log1p(
  as.numeric(all_clinical$`Fraction Genome Altered`)
)

# --- 1l. Mutation Count: continuous, log-transformed ------------------------
# Reference: right-skewed count variable; log1p transformation applied to
# reduce the influence of extreme hypermutators and satisfy the log-linear
# hazard assumption (Alexandrov et al., 2013, Nature 500:415-421).
all_clinical$Mutation_Count_log <- log1p(
  as.numeric(all_clinical$`Mutation Count`)
)

# --- 1m. TMB (nonsynonymous): continuous, log-transformed -------------------
# Reference: TMB is highly correlated with Mutation Count but measured on a
# per-Mb basis. Log1p transformation applied for the same reason as above.
# Note: including both TMB and Mutation Count may introduce collinearity;
# consider retaining only one in the final model
# (Xue et al., 2007, Lifetime Data Anal 13:333-350 — collinearity in Cox).
all_clinical$TMB_log <- log1p(as.numeric(all_clinical$`TMB (nonsynonymous)`))

# Verify collinearity between 'Mutation_Count_log' and 'TMB_log':
cat("\n--- Collinearity check: Mutation_Count_log vs TMB_log ---\n")
cat(sprintf("Pearson r = %.4f \n",
            cor(all_clinical$Mutation_Count_log, all_clinical$TMB_log,
                use = "complete.obs")))

# In TCGA-LIHC, TMB (nonsynonymous) and Mutation Count have a Pearson
# correlation of r = 1.000 (log1p-transformed: r = 0.976), confirming they
# are effectively the same measurement on different scales. Including both
# produces a non-identifiable model (Xue et al., 2007, Lifetime Data Anal
# 13:333-350). 
# Potential Solution: Mutation_Count_log can be retained as the more biologically
# interpretable quantity; meanwhile, TMB_log should be excluded from Model 3.
# (or the other way around?)

# ── 2. Helper functions ───────────────────────────────────────────────────────

# Extract HR table for only the four gene mutation rows from a coxph object
extract_gene_hr <- function(cox_model, genes = c("TP53","CTNNB1","JAK1","P4HA1")) {
  cf  <- summary(cox_model)$coefficients
  ci  <- confint(cox_model)
  # row names contain e.g. "TP53_mutationMutated"
  idx <- grep("_mutationMutated$", rownames(cf))
  data.frame(
    Gene     = genes,
    HR       = exp(cf[idx, "coef"]),
    CI_lower = exp(ci[idx, 1]),
    CI_upper = exp(ci[idx, 2]),
    p_value  = cf[idx, "Pr(>|z|)"],
    row.names = NULL
  )
}

# ── 3. Model 1: Univariate Cox (one gene at a time) ──────────────────────────
# This mirrors the log-rank test already shown in K-M section but returns HR.
genes <- c("TP53_mutation","CTNNB1_mutation","JAK1_mutation",
           "P4HA1_mutation")
gene_labels <- c("TP53","CTNNB1","JAK1","P4HA1")

uni_os_rows <- lapply(seq_along(genes), function(i) {
  f   <- as.formula(paste0("Surv(OS_time, OS_status) ~ ", genes[i]))
  fit <- coxph(f, data = all_clinical)
  cf  <- summary(fit)$coefficients
  ci  <- confint(fit)
  data.frame(Gene     = gene_labels[i],
             HR       = exp(cf[1, "coef"]),
             CI_lower = exp(ci[1, 1]),
             CI_upper = exp(ci[1, 2]),
             p_value  = cf[1, "Pr(>|z|)"])
})
uni_os_tbl <- do.call(rbind, uni_os_rows)

uni_dfs_rows <- lapply(seq_along(genes), function(i) {
  f   <- as.formula(paste0("Surv(DFS_time, DFS_status) ~ ", genes[i]))
  fit <- coxph(f, data = all_clinical)
  cf  <- summary(fit)$coefficients
  ci  <- confint(fit)
  data.frame(Gene     = gene_labels[i],
             HR       = exp(cf[1, "coef"]),
             CI_lower = exp(ci[1, 1]),
             CI_upper = exp(ci[1, 2]),
             p_value  = cf[1, "Pr(>|z|)"])
})
uni_dfs_tbl <- do.call(rbind, uni_dfs_rows)

cat("\n===== Model 1 – Univariate Cox: Overall Survival =====\n")
print(uni_os_tbl)
cat("\n===== Model 1 – Univariate Cox: Disease-Free Survival =====\n")
print(uni_dfs_tbl)

# ── 4. Model 2: Joint model – gene mutations only ────────────────────────────
cox_os_m2 <- coxph(
  Surv(OS_time, OS_status) ~
    TP53_mutation + CTNNB1_mutation + JAK1_mutation +
    P4HA1_mutation,
  data = all_clinical
)

cox_dfs_m2 <- coxph(
  Surv(DFS_time, DFS_status) ~
    TP53_mutation + CTNNB1_mutation + JAK1_mutation +
    P4HA1_mutation,
  data = all_clinical
)

cat("\n===== Model 2 – Joint (mutations only): OS =====\n")
print(summary(cox_os_m2))
cat("\n===== Model 2 – Joint (mutations only): DFS =====\n")
print(summary(cox_dfs_m2))

m2_os_tbl  <- extract_gene_hr(cox_os_m2)
m2_dfs_tbl <- extract_gene_hr(cox_dfs_m2)

# ── 5. Model 3: Fully adjusted – gene mutations + clinical covariates ─────────
cox_os_m3 <- coxph(
  Surv(OS_time, OS_status) ~
    # Gene mutations
    TP53_mutation + CTNNB1_mutation + JAK1_mutation +
    P4HA1_mutation +
    # Patient demographics
    Diagnosis_Age + Sex + Ethnicity_Category + Race_Category +
    # Tumour stage (overall stage used as primary staging variable;
    # T/N/M sub-stages are components of overall stage so including all four
    # together would introduce collinearity — see Groome et al., 2001, CA Cancer J Clin)
    AJCC_Stage +
    # Genomic burden
    Fraction_Genome_Altered_log + Mutation_Count_log + TMB_log +
    # Treatment history
    Prior_Malignancy + Prior_Treatment,
  data = all_clinical
)

cox_dfs_m3 <- coxph(
  Surv(DFS_time, DFS_status) ~
    TP53_mutation + CTNNB1_mutation + JAK1_mutation +
    P4HA1_mutation +
    Diagnosis_Age + Sex + Ethnicity_Category + Race_Category +
    AJCC_Stage +
    Fraction_Genome_Altered_log + Mutation_Count_log + TMB_log +
    Prior_Malignancy + Prior_Treatment,
  data = all_clinical
)

cat("\n===== Model 3 – Fully Adjusted: OS =====\n")
print(summary(cox_os_m3))
cat("\n===== Model 3 – Fully Adjusted: DFS =====\n")
print(summary(cox_dfs_m3))

m3_os_tbl  <- extract_gene_hr(cox_os_m3)
m3_dfs_tbl <- extract_gene_hr(cox_dfs_m3)

# ── 6. PH assumption checks for fully adjusted models ─────────────────────────
cat("\n===== PH Assumption Check (cox.zph) – Model 3 OS =====\n")
print(cox.zph(cox_os_m3))
cat("\n===== PH Assumption Check (cox.zph) – Model 3 DFS =====\n")
print(cox.zph(cox_dfs_m3))


# ── 7. Three-model comparison table ───────────────────────────────────────────
# Side-by-side HR (95% CI) and p-value for each gene across all three models.

fmt_hr <- function(hr, lo, hi, pv) {
  sig_star <- ifelse(pv < 0.001, "***",
              ifelse(pv < 0.01,  "**",
              ifelse(pv < 0.05,  "*", "")))
  sprintf("%.2f (%.2f\u2013%.2f)%s", hr, lo, hi, sig_star)
}

build_comparison_tbl <- function(uni, m2, m3) {
  data.frame(
    Gene = uni$Gene,
    Univariate          = mapply(fmt_hr, uni$HR, uni$CI_lower,
                                 uni$CI_upper, uni$p_value),
    Mutations_Only      = mapply(fmt_hr, m2$HR,  m2$CI_lower,
                                 m2$CI_upper,  m2$p_value),
    Fully_Adjusted      = mapply(fmt_hr, m3$HR,  m3$CI_lower,
                                 m3$CI_upper,  m3$p_value),
    stringsAsFactors = FALSE
  )
}

cat("\n")
cat("=================================================================\n")
cat("  HR (95% CI) Comparison Table – Overall Survival\n")
cat("  * p<0.05  ** p<0.01  *** p<0.001\n")
cat("=================================================================\n")
print(build_comparison_tbl(uni_os_tbl, m2_os_tbl, m3_os_tbl),
      row.names = FALSE)

cat("\n")
cat("=================================================================\n")
cat("  HR (95% CI) Comparison Table – Disease-Free Survival\n")
cat("  * p<0.05  ** p<0.01  *** p<0.001\n")
cat("=================================================================\n")
print(build_comparison_tbl(uni_dfs_tbl, m2_dfs_tbl, m3_dfs_tbl),
      row.names = FALSE)




# Comments:

# ISSUE 1: NOTE ON P4HA1 — COMPLETE SEPARATION:
# P4HA1 has only 1 mutated patient in TCGA-LIHC, and that patient experienced
# zero events (OS and DFS). This constitutes complete separation: the model
# cannot estimate a finite hazard ratio and returns HR ≈ 0 with CI [0, Inf].
# This limitation is documented here and should be stated in the Methods.
# Potential Solution: 
# P4HA1 can be EXCLUDED from Models 2 and 3 (Cox regression) but can be
# retained in Model 1 (univariate) and the K-M section for descriptive purposes.
# Reference: Heinze & Schemper (2002), Stat Med 21:2409-2419.

# ISSUE 2: NOTE ON AMERICAN INDIAN OR ALASKA NATIVE — COMPLETE SEPARATION:
# Similar complete separation problem as above (n = 2). 
# Retaining a 2-patient category produces HR = 0 or Inf with unbounded CIs,
# rendering the estimate meaningless (Heinze & Schemper, 2002, Stat Med
# 21:2409-2419 — complete separation in Cox regression).
# Potential Solution: 
# Merge AMERICAN INDIAN OR ALASKA NATIVE (n = 2) into a
# combined "Other" category together with BLACK OR AFRICAN AMERICAN (n = 17)
# to resolve complete separation caused by near-zero cell size.

# ISSUE 3:
# Prior Treatment has only 2 TRUE cases out of 379 (0.5%), producing complete
# separation and meaningless estimates (HR ~ 0, CI [0, Inf]).
# Per Peduzzi et al. (1995, J Clin Epidemiol 48:1503-1510), a minimum of ~10
# events per variable is required for stable Cox estimation; this variable
# has far fewer.
# Potential Solution: 
# Exclude Prior treatment from our fully-adjusted model (Model 3)







