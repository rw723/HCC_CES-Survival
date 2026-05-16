# ============================================================
# TCGA-LIHC survival analysis
# ============================================================

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(survival)
library(ggplot2)
library(readr)
library(survminer)
library(flextable)
library(officer)
library(cancereffectsizeR)

# -----------------------------
# 1. Inputs
# -----------------------------

genes_to_test <- c(
  "TP53",
  "CTNNB1",
  "JAK1",
  "P4HA1",
  "STAT3",
  "IL6ST"
)

tcga_clin_file <- "hcc_tcga_gdc_clinical_data.tsv"

out_dir <- "TCGA_LIHC_survival_outputs"
km_dir <- file.path(out_dir, "KM_plots")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(km_dir, showWarnings = FALSE, recursive = TRUE)

harmonize_tcga_patient_id <- function(x) {
  substr(as.character(x), 1, 12)
}

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  )
}

fmt_num <- function(x, digits = 2) {
  ifelse(
    is.na(x),
    NA_character_,
    format(round(x, digits), nsmall = digits, scientific = FALSE)
  )
}

# -----------------------------
# 2. Load and prepare TCGA clinical survival data
# -----------------------------

tcga_clin <- fread(tcga_clin_file, sep = "\t", quote = "") %>%
  as.data.frame()

tcga_surv <- tcga_clin %>%
  mutate(
    patient_id_raw = `Patient ID`,
    sample_id = `Sample ID`,
    patient_id = harmonize_tcga_patient_id(`Patient ID`),
    
    OS_time = as.numeric(`Overall Survival (Months)`),
    OS_status = case_when(
      `Overall Survival Status` == "1:DECEASED" ~ 1L,
      `Overall Survival Status` == "0:LIVING" ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    DFS_time = as.numeric(`Disease Free (Months)`),
    DFS_status = case_when(
      `Disease Free Status` == "1:Recurred/Progressed" ~ 1L,
      `Disease Free Status` == "0:DiseaseFree" ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    age = as.numeric(`Diagnosis Age`),
    sex = factor(Sex),
    race = factor(`Race Category`),
    ethnicity = factor(`Ethnicity Category`),
    
    stage_raw = `AJCC Pathologic Stage`,
    stage_simple = case_when(
      str_detect(stage_raw, "Stage I($|[^I])|Stage IA|Stage IB") ~ "I",
      str_detect(stage_raw, "Stage II($|[^I])|Stage IIA|Stage IIB") ~ "II",
      str_detect(stage_raw, "Stage III") ~ "III",
      str_detect(stage_raw, "Stage IV") ~ "IV",
      TRUE ~ NA_character_
    ),
    stage_simple = factor(stage_simple, levels = c("I", "II", "III", "IV"))
  ) %>%
  arrange(patient_id, desc(`Sample Type` == "Primary Tumor")) %>%
  distinct(patient_id, .keep_all = TRUE)

cat("TCGA clinical rows:", nrow(tcga_clin), "\n")
cat("TCGA unique survival patients:", nrow(tcga_surv), "\n")
cat(
  "OS usable:",
  sum(!is.na(tcga_surv$OS_time) &
        !is.na(tcga_surv$OS_status) &
        tcga_surv$OS_time > 0),
  "\n"
)
cat(
  "DFS usable:",
  sum(!is.na(tcga_surv$DFS_time) &
        !is.na(tcga_surv$DFS_status) &
        tcga_surv$DFS_time > 0),
  "\n"
)

# -----------------------------
# 3. Load TCGA MAF
# -----------------------------

tcga_maf <- preload_maf(
  maf = "./lihc_tcga_hg38.maf",
  refset = "ces.refset.hg38",
  keep_extra_columns = TRUE
)

tcga_maf_anno <- as.data.frame(tcga_maf)

required_cols <- c(
  "Unique_Patient_Identifier",
  "Hugo_Symbol",
  "Variant_Classification"
)

missing_cols <- setdiff(required_cols, names(tcga_maf_anno))

if (length(missing_cols) > 0) {
  stop(
    "tcga_maf_anno is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("Filtered TCGA MAF records:", nrow(tcga_maf_anno), "\n")
cat("Variant_Classification values:\n")
print(sort(table(tcga_maf_anno$Variant_Classification, useNA = "ifany"), decreasing = TRUE))

# -----------------------------
# 4. Define nonsynonymous mutation classes
# -----------------------------

nonsyn_classes <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "Translation_Start_Site",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Splice_Site",
  "Splice_Region"
)

present_nonsyn_classes <- intersect(
  nonsyn_classes,
  unique(tcga_maf_anno$Variant_Classification)
)

cat("Nonsynonymous classes present in this MAF:\n")
print(present_nonsyn_classes)

# -----------------------------
# 5. Build gene-level mutation matrix
# -----------------------------

tcga_surv_mut_long <- tcga_maf_anno %>%
  mutate(
    patient_id = harmonize_tcga_patient_id(Unique_Patient_Identifier),
    gene = Hugo_Symbol,
    variant_class = Variant_Classification
  ) %>%
  filter(
    gene %in% genes_to_test,
    variant_class %in% nonsyn_classes
  ) %>%
  transmute(
    patient_id,
    gene,
    variant_class,
    sample_id = Unique_Patient_Identifier,
    mutated = 1L
  ) %>%
  distinct(patient_id, gene, .keep_all = TRUE)

tcga_surv_mut_counts <- tcga_surv_mut_long %>%
  count(gene, name = "n_mutated_post_filter_nonsyn") %>%
  arrange(desc(n_mutated_post_filter_nonsyn))

cat("Mutation-positive counts from filtered MAF annotation:\n")
print(tcga_surv_mut_counts)

tcga_surv_mut_wide <- tcga_surv_mut_long %>%
  select(patient_id, gene, mutated) %>%
  pivot_wider(
    names_from = gene,
    values_from = mutated,
    values_fill = 0,
    names_prefix = "mut_"
  )

for (g in genes_to_test) {
  col <- paste0("mut_", g)
  if (!col %in% names(tcga_surv_mut_wide)) {
    tcga_surv_mut_wide[[col]] <- 0L
  }
}

# -----------------------------
# 6. Merge mutation matrix with clinical survival table
# -----------------------------

surv_data <- tcga_surv %>%
  left_join(tcga_surv_mut_wide, by = "patient_id")

for (g in genes_to_test) {
  col <- paste0("mut_", g)
  if (!col %in% names(surv_data)) {
    surv_data[[col]] <- 0L
  }
  surv_data[[col]][is.na(surv_data[[col]])] <- 0L
}

# -----------------------------
# 7. Explicit N used for OS and DFS
# -----------------------------

n_used_summary <- tibble(
  Endpoint = c("OS", "DFS"),
  `N used` = c(
    sum(!is.na(surv_data$OS_time) &
          !is.na(surv_data$OS_status) &
          surv_data$OS_time > 0),
    sum(!is.na(surv_data$DFS_time) &
          !is.na(surv_data$DFS_status) &
          surv_data$DFS_time > 0)
  ),
  `Number of events` = c(
    sum(!is.na(surv_data$OS_time) &
          !is.na(surv_data$OS_status) &
          surv_data$OS_time > 0 &
          surv_data$OS_status == 1),
    sum(!is.na(surv_data$DFS_time) &
          !is.na(surv_data$DFS_status) &
          surv_data$DFS_time > 0 &
          surv_data$DFS_status == 1)
  )
)

cat("\nN used for survival analyses:\n")
print(n_used_summary)

# -----------------------------
# 8. Mutation counts after clinical merge
# -----------------------------

surv_gene_counts <- map_dfr(genes_to_test, function(g) {
  mut_col <- paste0("mut_", g)
  
  tibble(
    gene = g,
    n_total_clinical = nrow(surv_data),
    n_mutated = sum(surv_data[[mut_col]] == 1, na.rm = TRUE),
    n_wildtype = sum(surv_data[[mut_col]] == 0, na.rm = TRUE)
  )
})

cat("Mutation counts after clinical merge:\n")
print(surv_gene_counts)

# -----------------------------
# 9. Event-count diagnostic
# -----------------------------

surv_gene_event_counts <- map_dfr(genes_to_test, function(g) {
  mut_col <- paste0("mut_", g)
  
  tibble(
    gene = g,
    n_mutated_total = sum(surv_data[[mut_col]] == 1, na.rm = TRUE),
    
    OS_n_analyzed = sum(!is.na(surv_data$OS_time) &
                          !is.na(surv_data$OS_status) &
                          surv_data$OS_time > 0),
    OS_events_mutated = sum(
      surv_data[[mut_col]] == 1 &
        surv_data$OS_status == 1 &
        !is.na(surv_data$OS_time) &
        surv_data$OS_time > 0,
      na.rm = TRUE
    ),
    OS_events_wildtype = sum(
      surv_data[[mut_col]] == 0 &
        surv_data$OS_status == 1 &
        !is.na(surv_data$OS_time) &
        surv_data$OS_time > 0,
      na.rm = TRUE
    ),
    
    DFS_n_analyzed = sum(!is.na(surv_data$DFS_time) &
                           !is.na(surv_data$DFS_status) &
                           surv_data$DFS_time > 0),
    DFS_events_mutated = sum(
      surv_data[[mut_col]] == 1 &
        surv_data$DFS_status == 1 &
        !is.na(surv_data$DFS_time) &
        surv_data$DFS_time > 0,
      na.rm = TRUE
    ),
    DFS_events_wildtype = sum(
      surv_data[[mut_col]] == 0 &
        surv_data$DFS_status == 1 &
        !is.na(surv_data$DFS_time) &
        surv_data$DFS_time > 0,
      na.rm = TRUE
    ),
    
    recommendation = case_when(
      n_mutated_total < 5 ~ "Too few mutated cases; report descriptively only",
      OS_events_mutated < 3 & DFS_events_mutated < 3 ~ "Too few mutation-positive events; interpret Cox cautiously",
      TRUE ~ "Run KM/log-rank and Cox with caution"
    )
  )
})

cat("Event-count diagnostic:\n")
print(surv_gene_event_counts)

# -----------------------------
# 10. Survival-analysis function
# -----------------------------

run_survival_one <- function(data, gene, endpoint = c("OS", "DFS"),
                             min_mutated = 5,
                             min_total_events = 5,
                             min_mutated_events_for_reliable_cox = 3) {
  
  endpoint <- match.arg(endpoint)
  mut_col <- paste0("mut_", gene)
  
  if (!mut_col %in% names(data)) {
    stop("Mutation column not found: ", mut_col)
  }
  
  if (endpoint == "OS") {
    time_col <- "OS_time"
    status_col <- "OS_status"
  } else {
    time_col <- "DFS_time"
    status_col <- "DFS_status"
  }
  
  df <- data %>%
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[status_col]]),
      .data[[time_col]] > 0,
      !is.na(.data[[mut_col]])
    ) %>%
    mutate(
      time = as.numeric(.data[[time_col]]),
      status = as.integer(.data[[status_col]]),
      mut_status = as.integer(.data[[mut_col]])
    )
  
  n_total <- nrow(df)
  n_mut <- sum(df$mut_status == 1, na.rm = TRUE)
  n_wt <- sum(df$mut_status == 0, na.rm = TRUE)
  events_mut <- sum(df$status == 1 & df$mut_status == 1, na.rm = TRUE)
  events_wt <- sum(df$status == 1 & df$mut_status == 0, na.rm = TRUE)
  total_events <- events_mut + events_wt
  
  base_row <- tibble(
    gene = gene,
    endpoint = endpoint,
    n_total = n_total,
    n_cox = NA_integer_,
    n_mutated = n_mut,
    n_wildtype = n_wt,
    events_mutated = events_mut,
    events_wildtype = events_wt,
    logrank_p = NA_real_,
    cox_hr = NA_real_,
    cox_ci_low = NA_real_,
    cox_ci_high = NA_real_,
    cox_p = NA_real_,
    ph_p = NA_real_,
    model_type = NA_character_,
    covariates = NA_character_,
    note = NA_character_
  )
  
  if (n_mut < min_mutated || n_wt < min_mutated || total_events < min_total_events) {
    return(
      base_row %>%
        mutate(note = "Insufficient mutation-positive cases, wild-type cases, or total events")
    )
  }
  
  surv_obj <- Surv(df$time, df$status)
  
  lr <- tryCatch(
    survdiff(surv_obj ~ mut_status, data = df),
    error = function(e) NULL
  )
  
  logrank_p <- if (!is.null(lr)) {
    1 - pchisq(lr$chisq, length(lr$n) - 1)
  } else {
    NA_real_
  }
  
  candidate_covars <- c("age", "sex", "stage_simple", "race", "ethnicity")
  candidate_covars <- candidate_covars[candidate_covars %in% names(df)]
  
  usable_covars <- candidate_covars[
    sapply(candidate_covars, function(v) {
      x <- df[[v]]
      sum(!is.na(x)) >= 30 && length(unique(na.omit(x))) >= 2
    })
  ]
  
  if (length(usable_covars) > 0) {
    cox_vars <- c("time", "status", "mut_status", usable_covars)
    
    df_cox_adj <- df %>%
      select(all_of(cox_vars)) %>%
      drop_na()
    
    n_events_cox_adj <- sum(df_cox_adj$status == 1)
    n_params_adj <- 1 + length(usable_covars)
  } else {
    df_cox_adj <- NULL
    n_events_cox_adj <- 0
    n_params_adj <- Inf
  }
  
  if (length(usable_covars) == 0 ||
      n_events_cox_adj < 10 * n_params_adj ||
      nrow(df_cox_adj) < 50) {
    
    df_cox <- df %>%
      select(time, status, mut_status) %>%
      drop_na()
    
    formula_use <- Surv(time, status) ~ mut_status
    covar_text <- "mutation status only"
    model_type <- "unadjusted Cox"
    
  } else {
    
    df_cox <- df_cox_adj
    
    formula_use <- as.formula(
      paste("Surv(time, status) ~ mut_status +", paste(usable_covars, collapse = " + "))
    )
    
    covar_text <- paste(c("mutation status", usable_covars), collapse = ", ")
    model_type <- "covariate-adjusted Cox"
  }
  
  cox_fit <- tryCatch(
    coxph(formula_use, data = df_cox, x = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(cox_fit)) {
    return(
      base_row %>%
        mutate(
          logrank_p = logrank_p,
          n_cox = nrow(df_cox),
          model_type = model_type,
          covariates = covar_text,
          note = "Cox model failed"
        )
    )
  }
  
  cox_sum <- summary(cox_fit)
  
  if (!"mut_status" %in% rownames(cox_sum$coefficients)) {
    return(
      base_row %>%
        mutate(
          logrank_p = logrank_p,
          n_cox = nrow(df_cox),
          model_type = model_type,
          covariates = covar_text,
          note = "Mutation coefficient not estimable"
        )
    )
  }
  
  hr <- cox_sum$coefficients["mut_status", "exp(coef)"]
  p <- cox_sum$coefficients["mut_status", "Pr(>|z|)"]
  ci_low <- cox_sum$conf.int["mut_status", "lower .95"]
  ci_high <- cox_sum$conf.int["mut_status", "upper .95"]
  
  ph <- tryCatch(
    cox.zph(cox_fit),
    error = function(e) NULL
  )
  
  ph_p <- if (!is.null(ph) && "mut_status" %in% rownames(ph$table)) {
    ph$table["mut_status", "p"]
  } else {
    NA_real_
  }
  
  note_text <- if (events_mut < min_mutated_events_for_reliable_cox) {
    "Few mutation-positive events; interpret Cox cautiously"
  } else {
    "OK"
  }
  
  tibble(
    gene = gene,
    endpoint = endpoint,
    n_total = n_total,
    n_cox = nrow(df_cox),
    n_mutated = n_mut,
    n_wildtype = n_wt,
    events_mutated = events_mut,
    events_wildtype = events_wt,
    logrank_p = logrank_p,
    cox_hr = hr,
    cox_ci_low = ci_low,
    cox_ci_high = ci_high,
    cox_p = p,
    ph_p = ph_p,
    model_type = model_type,
    covariates = covar_text,
    note = note_text
  )
}

# -----------------------------
# 11. Run survival analyses
# -----------------------------

survival_summary <- expand_grid(
  gene = genes_to_test,
  endpoint = c("OS", "DFS")
) %>%
  mutate(
    result = map2(
      gene,
      endpoint,
      ~ run_survival_one(surv_data, .x, .y)
    )
  ) %>%
  select(result) %>%
  unnest(result) %>%
  group_by(endpoint) %>%
  mutate(
    cox_p_BH = {
      out <- rep(NA_real_, length(cox_p))
      ok <- !is.na(cox_p)
      out[ok] <- p.adjust(cox_p[ok], method = "BH")
      out
    }
  ) %>%
  ungroup() %>%
  mutate(
    HR_95CI = ifelse(
      is.na(cox_hr) | is.na(cox_ci_low) | is.na(cox_ci_high),
      NA_character_,
      paste0(
        fmt_num(cox_hr, 2), " (",
        fmt_num(cox_ci_low, 2), "–",
        fmt_num(cox_ci_high, 2), ")"
      )
    )
  )

cat("Full survival summary:\n")
print(survival_summary, n = Inf)

cat("\nN used by endpoint from survival_summary:\n")
print(
  survival_summary %>%
    distinct(endpoint, n_total) %>%
    arrange(endpoint)
)

# -----------------------------
# 12. Clean publication table
# -----------------------------

survival_table_clean <- survival_summary %>%
  transmute(
    Gene = gene,
    Endpoint = endpoint,
    `N used` = n_total,
    `N used in Cox` = n_cox,
    `N mutated` = n_mutated,
    `N wild-type` = n_wildtype,
    `Events mutated` = events_mutated,
    `Events wild-type` = events_wildtype,
    `Log-rank p-value` = fmt_p(logrank_p),
    `Cox HR (95% CI)` = HR_95CI,
    `Cox p-value` = fmt_p(cox_p),
    `Cox BH p-value` = fmt_p(cox_p_BH),
    `PH test p-value` = fmt_p(ph_p)
  ) %>%
  mutate(
    across(everything(), as.character),
    across(
      everything(),
      ~ ifelse(is.na(.x) | .x == "" | .x == "NA (NA–NA)", "NA", .x)
    )
  )

cat("Clean survival table:\n")
print(survival_table_clean, n = Inf)

# -----------------------------
# 13. Kaplan-Meier plots
# -----------------------------

if (!requireNamespace("ragg", quietly = TRUE)) {
  stop("Please install 'ragg' first: install.packages('ragg')")
}

format_logrank_p <- function(fit, data) {
  pdat <- survminer::surv_pvalue(fit, data = data)
  pval <- pdat$pval
  
  if (is.null(pval) || length(pval) == 0 || is.na(pval)) {
    return("Log-rank\np = NA")
  } else if (pval < 0.0001) {
    return("Log-rank\np < 0.0001")
  } else {
    return(paste0("Log-rank\np = ", formatC(pval, format = "f", digits = 4)))
  }
}

choose_xmax_months <- function(x) {
  xmax <- max(x, na.rm = TRUE)
  
  if (xmax <= 24) {
    24
  } else if (xmax <= 36) {
    36
  } else if (xmax <= 48) {
    48
  } else if (xmax <= 60) {
    60
  } else if (xmax <= 72) {
    72
  } else if (xmax <= 96) {
    96
  } else {
    120
  }
}

choose_break_months <- function(xmax) {
  if (xmax <= 24) {
    6
  } else if (xmax <= 60) {
    12
  } else {
    24
  }
}

plot_km_one <- function(data, gene, endpoint = c("OS", "DFS"),
                        min_mutated_for_plot = 5) {
  
  endpoint <- match.arg(endpoint)
  mut_col <- paste0("mut_", gene)
  
  if (!mut_col %in% colnames(data)) {
    message("Skipping ", gene, " ", endpoint, ": mutation column not found -> ", mut_col)
    return(NA_character_)
  }
  
  if (endpoint == "OS") {
    time_col <- "OS_time"
    status_col <- "OS_status"
    endpoint_short <- "OS"
  } else {
    time_col <- "DFS_time"
    status_col <- "DFS_status"
    endpoint_short <- "DFS"
  }
  
  df <- data %>%
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[status_col]]),
      .data[[time_col]] > 0,
      !is.na(.data[[mut_col]])
    ) %>%
    mutate(
      time = as.numeric(.data[[time_col]]),
      status = as.integer(.data[[status_col]]),
      Mutation = factor(
        ifelse(.data[[mut_col]] == 1, "Mutated", "Wild-type"),
        levels = c("Wild-type", "Mutated")
      )
    )
  
  if (nrow(df) == 0) {
    message("Skipping ", gene, " ", endpoint, ": no valid data")
    return(NA_character_)
  }
  
  if (length(unique(df$Mutation)) < 2) {
    message("Skipping ", gene, " ", endpoint, ": only one mutation group present")
    return(NA_character_)
  }
  
  n_mut <- sum(df$Mutation == "Mutated", na.rm = TRUE)
  
  if (n_mut < min_mutated_for_plot) {
    message("Skipping ", gene, " ", endpoint, ": n_mutated < ", min_mutated_for_plot)
    return(NA_character_)
  }
  
  fit <- survival::survfit(
    survival::Surv(time, status) ~ Mutation,
    data = df
  )
  
  p_label <- format_logrank_p(fit, df)
  
  x_max <- choose_xmax_months(df$time)
  break_by <- choose_break_months(x_max)
  x_breaks <- seq(0, x_max, by = break_by)
  
  col_wt <- "#3A3A3A"
  col_mut <- "#71CAD0"
  
  p <- survminer::ggsurvplot(
    fit = fit,
    data = df,
    
    conf.int = TRUE,
    conf.int.style = "ribbon",
    conf.int.alpha = 0.20,
    censor = TRUE,
    censor.shape = 124,
    censor.size = 2.4,
    size = 1.0,
    palette = c(col_wt, col_mut),
    
    risk.table = "nrisk_cumcensor",
    risk.table.title = "No. at risk (no. censored)",
    risk.table.height = 0.32,
    risk.table.col = "strata",
    risk.table.y.text = FALSE,
    risk.table.y.text.col = TRUE,
    fontsize = 4.8,
    
    xlim = c(0, x_max),
    break.time.by = break_by,
    axes.offset = FALSE,
    
    title = paste0(endpoint_short, " by ", gene, " mutation status"),
    xlab = "Time from diagnosis (months)",
    ylab = "Survival probability",
    
    legend.title = NULL,
    legend.labs = c("Wild-type", "Mutated"),
    legend = c(0.98, 0.98),
    
    pval = FALSE,
    
    ggtheme = ggplot2::theme_classic(base_size = 12),
    tables.theme = ggplot2::theme_classic(base_size = 10)
  )
  
  p$plot <- p$plot +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max),
      breaks = x_breaks,
      expand = ggplot2::expansion(mult = c(0.005, 0.035))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.02),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0.005, 0.015))
    ) +
    ggplot2::annotate(
      "text",
      x = 0.04 * x_max,
      y = 0.08,
      label = p_label,
      hjust = 0,
      vjust = 0,
      size = 4.1
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20, hjust = 0.5, face = "plain"),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 11, colour = "black"),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.4
      ),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key.size = grid::unit(0.45, "cm"),
      legend.margin = ggplot2::margin(3, 5, 3, 5),
      legend.text = ggplot2::element_text(size = 10.5),
      legend.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 8, r = 28, b = 4, l = 14)
    )
  
  p$table <- p$table +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max),
      breaks = x_breaks,
      expand = ggplot2::expansion(mult = c(0.005, 0.035))
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, hjust = 0, face = "bold"),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 11, face = "bold", colour = "black"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 4, r = 28, b = 10, l = 14)
    )
  
  for (i in seq_along(p$table$layers)) {
    if (inherits(p$table$layers[[i]]$geom, "GeomText")) {
      p$table$layers[[i]]$aes_params$fontface <- "bold"
      p$table$layers[[i]]$aes_params$size <- 4.2
    }
  }
  
  out_png <- file.path(km_dir, paste0("KM_", gene, "_", endpoint, ".png"))
  
  ragg::agg_png(
    filename = out_png,
    width = 10.2,
    height = 7.4,
    units = "in",
    res = 300
  )
  
  print(
    p,
    surv.plot.height = 0.68,
    risk.table.height = 0.32,
    newpage = FALSE
  )
  
  grDevices::dev.off()
  
  return(out_png)
}

km_outputs <- tidyr::expand_grid(
  gene = genes_to_test,
  endpoint = c("OS", "DFS")
) %>%
  mutate(
    file = map2_chr(
      gene,
      endpoint,
      ~ {
        x <- plot_km_one(
          data = surv_data,
          gene = .x,
          endpoint = .y,
          min_mutated_for_plot = 5
        )
        ifelse(is.na(x), NA_character_, x)
      }
    )
  )

print(km_outputs)

# -----------------------------
# 14. Word-ready table
# -----------------------------

ft_survival <- survival_table_clean %>%
  flextable() %>%
  set_caption(
    caption = "S1 Table. Survival analysis summary for selected driver mutations in TCGA-LIHC."
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  
  compose(
    j = "Log-rank p-value",
    part = "header",
    value = as_paragraph(
      "Log-rank ",
      as_i("p"),
      "-value"
    )
  ) %>%
  compose(
    j = "Cox p-value",
    part = "header",
    value = as_paragraph(
      "Cox ",
      as_i("p"),
      "-value"
    )
  ) %>%
  compose(
    j = "Cox BH p-value",
    part = "header",
    value = as_paragraph(
      "Cox BH ",
      as_i("p"),
      "-value"
    )
  ) %>%
  compose(
    j = "PH test p-value",
    part = "header",
    value = as_paragraph(
      "PH test ",
      as_i("p"),
      "-value"
    )
  ) %>%
  
  bg(
    i = ~ Endpoint == "OS",
    bg = "white",
    part = "body"
  ) %>%
  bg(
    i = ~ Endpoint == "DFS",
    bg = "gray80",
    part = "body"
  ) %>%
  
  align(align = "center", part = "all") %>%
  align(
    j = "Gene",
    align = "left",
    part = "all"
  ) %>%
  fontsize(size = 8, part = "all") %>%
  padding(padding = 3, part = "all") %>%
  autofit()

ft_survival

doc <- read_docx() %>%
  body_add_flextable(ft_survival)

print(
  doc,
  target = file.path(out_dir, "TCGA_LIHC_survival_table.docx")
)

# -----------------------------
# 15. Save analysis-ready object
# -----------------------------

saveRDS(
  list(
    surv_data = surv_data,
    n_used_summary = n_used_summary,
    mutation_counts = surv_gene_counts,
    event_counts = surv_gene_event_counts,
    survival_summary = survival_summary,
    survival_table_clean = survival_table_clean,
    km_outputs = km_outputs
  ),
  file = file.path(out_dir, "TCGA_survival_analysis_results_final.rds")
)

cat("\nDone. Outputs saved in:", out_dir, "\n")

