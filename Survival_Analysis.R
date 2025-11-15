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

# read in the data
all_clinical <- read_tsv("tcga_clinical/tcga_all_clinical.tsv", quote = "", col_names = TRUE,
                         name_repair = "minimal")
TP53_mutation <- read_tsv("tcga_clinical/TP53_mutations.tsv", quote = "", col_names = TRUE,
                          name_repair = "minimal")
CTNNB1_mutation <- read_tsv("tcga_clinical/CTNNB1_mutations.tsv", quote = "", col_names = TRUE,
                            name_repair = "minimal")
JAK1_mutation <- read_tsv("tcga_clinical/JAK1_mutations.tsv", quote = "", col_names = TRUE,
                          name_repair = "minimal")
P4HA1_mutation <- read_tsv("tcga_clinical/P4HA1_mutations.tsv", quote = "", col_names = TRUE,
                           name_repair = "minimal")
TLE1_mutation <- read_tsv("tcga_clinical/TLE1_mutations.tsv", quote = "", col_names = TRUE,
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
# Figure 5: K-M Plots of TP53 Mutation for Overall Survival and Disease-Free Survival
# Revise the K-M plot of TP53 Mutation with customized themes
# make consistent mutation labeling
all_clinical$TP53_mutation <- factor(
  all_clinical$TP53_mutation,
  levels = c("No", "Yes"),
  labels = c("Wild-type", "Mutated")
)

# ---------- Overall Survival ----------
fit_os <- survfit(Surv(OS_time, OS_status) ~ TP53_mutation, data = all_clinical)
lr_os  <- survdiff(Surv(OS_time, OS_status) ~ TP53_mutation, data = all_clinical)
p_os   <- 1 - pchisq(lr_os$chisq, df = length(lr_os$n) - 1)

g_os <- ggsurvplot(
  fit_os, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = sprintf("P = %.3f", p_os),
  title = expression("K-M Curve for " * italic("TP53") * " Mutation – Overall Survival"),
  xlab = "Overall survival (months)", 
  ylab = "Survival probability",
  legend.title = expression(italic("TP53") * " mutation"),
  legend.labs = levels(all_clinical$TP53_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = 14)
)
g_os$plot <- g_os$plot + theme(legend.position = "top")

# ---------- Disease-Free Survival ----------
fit_dfs <- survfit(Surv(DFS_time, DFS_status) ~ TP53_mutation, data = all_clinical)
lr_dfs  <- survdiff(Surv(DFS_time, DFS_status) ~ TP53_mutation, data = all_clinical)
p_dfs   <- 1 - pchisq(lr_dfs$chisq, df = length(lr_dfs$n) - 1)

g_dfs <- ggsurvplot(
  fit_dfs, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = sprintf("P = %.4f", p_dfs),
  title = expression("K-M Curve for " * italic("TP53") * " Mutation – Disease-Free Survival"),
  xlab = "Disease-free survival (months)", 
  ylab = "Survival probability",
  legend.title = expression(italic("TP53") * " mutation"),
  legend.labs = levels(all_clinical$TP53_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic(base_size = 14)
)
g_dfs$plot <- g_dfs$plot + theme(legend.position = "top")

# ---------- Combine panels ----------
fig5 <- g_os$plot + g_dfs$plot +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 18, face = "plain"))

# show
fig5

# save the figure
# ggsave("Fig5_TP53_OS_DFS.png", fig5, width = 15, height = 6, dpi = 300)


#########################################################################
# Figure Supplementary A2: K-M Plots of CTNNB1, JAK1, and P4HA1 Mutations 
# Make each mutation column use consistent labels
all_clinical$CTNNB1_mutation <- factor(all_clinical$CTNNB1_mutation, levels = c("No", "Yes"),
                                       labels = c("Wild-type","Mutated"))
all_clinical$JAK1_mutation   <- factor(all_clinical$JAK1_mutation, levels = c("No", "Yes"),
                                       labels = c("Wild-type","Mutated"))
all_clinical$P4HA1_mutation  <- factor(all_clinical$P4HA1_mutation,  levels = c("No", "Yes"),
                                       labels = c("Wild-type","Mutated"))

# ---------------- CTNNB1 (OS) ----------------
fit_os_CTNNB1 <- survfit(Surv(OS_time, OS_status) ~ CTNNB1_mutation, data = all_clinical)
lr_os_CTNNB1  <- survdiff(Surv(OS_time, OS_status) ~ CTNNB1_mutation, data = all_clinical)
p_os_CTNNB1   <- 1 - pchisq(lr_os_CTNNB1$chisq, df = length(lr_os_CTNNB1$n) - 1)

p1 <- ggsurvplot(
  fit_os_CTNNB1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = sprintf("P = %.3f", p_os_CTNNB1),
  title = expression("K-M curve for " * italic("CTNNB1") * " Mutation – Overall Survival"),
  xlab = "Overall survival (months)", ylab = "Survival probability",
  legend.title = expression(italic("CTNNB1") * " mutation"),
  legend.labs = levels(all_clinical$CTNNB1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE, ggtheme = theme_classic(base_size = 14)
)$plot + theme(legend.position = "top")

# ---------------- CTNNB1 (DFS) ----------------
fit_dfs_CTNNB1 <- survfit(Surv(DFS_time, DFS_status) ~ CTNNB1_mutation, data = all_clinical)
lr_dfs_CTNNB1  <- survdiff(Surv(DFS_time, DFS_status) ~ CTNNB1_mutation, data = all_clinical)
p_dfs_CTNNB1   <- 1 - pchisq(lr_dfs_CTNNB1$chisq, df = length(lr_dfs_CTNNB1$n) - 1)

p2 <- ggsurvplot(
  fit_dfs_CTNNB1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = sprintf("P = %.4f", p_dfs_CTNNB1),
  title = expression("K-M curve for " * italic("CTNNB1") * " Mutation – Disease-Free Survival"),
  xlab = "Disease-free survival (months)", ylab = "Survival probability",
  legend.title = expression(italic("CTNNB1") * " mutation"),
  legend.labs = levels(all_clinical$CTNNB1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE, ggtheme = theme_classic(base_size = 14)
)$plot + theme(legend.position = "top")

# ---------------- JAK1 (OS) ----------------
fit_os_JAK1 <- survfit(Surv(OS_time, OS_status) ~ JAK1_mutation, data = all_clinical)
lr_os_JAK1  <- survdiff(Surv(OS_time, OS_status) ~ JAK1_mutation, data = all_clinical)
p_os_JAK1   <- 1 - pchisq(lr_os_JAK1$chisq, df = length(lr_os_JAK1$n) - 1)

p3 <- ggsurvplot(
  fit_os_JAK1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = sprintf("P = %.3f", p_os_JAK1),
  title = expression("K-M curve for " * italic("JAK1") * " Mutation – Overall Survival"),
  xlab = "Overall survival (months)", ylab = "Survival probability",
  legend.title = expression(italic("JAK1") * " mutation"),
  legend.labs = levels(all_clinical$JAK1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE, ggtheme = theme_classic(base_size = 14)
)$plot + theme(legend.position = "top")

# ---------------- JAK1 (DFS) ----------------
fit_dfs_JAK1 <- survfit(Surv(DFS_time, DFS_status) ~ JAK1_mutation, data = all_clinical)
lr_dfs_JAK1  <- survdiff(Surv(DFS_time, DFS_status) ~ JAK1_mutation, data = all_clinical)
p_dfs_JAK1   <- 1 - pchisq(lr_dfs_JAK1$chisq, df = length(lr_dfs_JAK1$n) - 1)

p4 <- ggsurvplot(
  fit_dfs_JAK1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = sprintf("P = %.4f", p_dfs_JAK1),
  title = expression("K-M curve for " * italic("JAK1") * " Mutation – Disease-Free Survival"),
  xlab = "Disease-free survival (months)", ylab = "Survival probability",
  legend.title = expression(italic("JAK1") * " mutation"),
  legend.labs = levels(all_clinical$JAK1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE, ggtheme = theme_classic(base_size = 14)
)$plot + theme(legend.position = "top")

# ---------------- P4HA1 (OS) ----------------
fit_os_P4HA1 <- survfit(Surv(OS_time, OS_status) ~ P4HA1_mutation, data = all_clinical)
lr_os_P4HA1  <- survdiff(Surv(OS_time, OS_status) ~ P4HA1_mutation, data = all_clinical)
p_os_P4HA1   <- 1 - pchisq(lr_os_P4HA1$chisq, df = length(lr_os_P4HA1$n) - 1)

p5 <- ggsurvplot(
  fit_os_P4HA1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = sprintf("P = %.3f", p_os_P4HA1),
  title = expression("K-M curve for " * italic("P4HA1") * " Mutation – Overall Survival"),
  xlab = "Overall survival (months)", ylab = "Survival probability",
  legend.title = expression(italic("P4HA1") * " mutation"),
  legend.labs = levels(all_clinical$P4HA1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE, ggtheme = theme_classic(base_size = 14)
)$plot + theme(legend.position = "top")

# ---------------- P4HA1 (DFS) ----------------
fit_dfs_P4HA1 <- survfit(Surv(DFS_time, DFS_status) ~ P4HA1_mutation, data = all_clinical)
lr_dfs_P4HA1  <- survdiff(Surv(DFS_time, DFS_status) ~ P4HA1_mutation, data = all_clinical)
p_dfs_P4HA1   <- 1 - pchisq(lr_dfs_P4HA1$chisq, df = length(lr_dfs_P4HA1$n) - 1)

p6 <- ggsurvplot(
  fit_dfs_P4HA1, data = all_clinical,
  conf.int = TRUE, censor = TRUE,
  pval = sprintf("P = %.4f", p_dfs_P4HA1),
  title = expression("K-M curve for " * italic("P4HA1") * " Mutation – Disease-Free Survival"),
  xlab = "Disease-free survival (months)", ylab = "Survival probability",
  legend.title = expression(italic("P4HA1") * " mutation"),
  legend.labs = levels(all_clinical$P4HA1_mutation),
  xlim = c(0, 120), break.time.by = 20,
  risk.table = FALSE, ggtheme = theme_classic(base_size = 14)
)$plot + theme(legend.position = "top")

# -------- Arrange as 3 rows × 2 columns (OS | DFS) --------
fig_3x2 <- (p1 | p2) /
  (p3 | p4) /
  (p5 | p6)

fig_3x2 <- fig_3x2 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "plain"))

# show
fig_3x2

# save
# ggsave("Figure_CTNNB1_JAK1_P4HA1_3x2_OS_DFS.png", fig_3x2,
# width = 15, height = 20, dpi = 300)






#########################################################################
# UPDATED Figure 5: K-M Plots of TP53 Mutation for Overall Survival and Disease-Free Survival
# Revise the K-M plot of TP53 Mutation with customized themes

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
    plot.tag.position = c(0, 1)
  )

# show
fig5

# save
ggsave("Fig5_TP53_OS_DFS_revised.png", fig5, width = 12, height = 6, dpi = 300)


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
    plot.tag.position = c(0, 1)
  )

# show
fig_3x2

# save
ggsave("Figure_CTNNB1_JAK1_P4HA1_3x2_OS_DFS_revised.png",
       fig_3x2, width = 15, height = 20, dpi = 300)

