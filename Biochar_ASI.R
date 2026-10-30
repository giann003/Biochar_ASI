############################################################
## Robust BIOCHAR ASI workflow avoiding row-mismatch errors
## Georgios Giannopoulos (c) 2026
##
## Description: Complete analytical pipeline for biochar 
## multi-criteria evaluation, regulatory metal screening, 
## multivariate profiling (NMDS, Mantel), and Aggregated 
## Suitability Index (ASI) calculation using multi-criteria weights.
############################################################

# Set working directory (adjust as needed for local environment)
setwd("~/Fenix")

############################################################
## 1. Directory Structure & Output Folders Setup
############################################################
out_dir <- "output"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

############################################################
## 2. Package Loading
############################################################
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(vegan)
library(cluster)
library(Hmisc)
library(corrplot)
library(scales)

############################################################
## 3. Workflow Pipeline Flowchart Generation
############################################################
grViz("
digraph biochar_flow {
  graph [layout = dot, rankdir = TB]

  node [shape = rectangle,
        style = filled,
        fillcolor = \"#f7f7f7\",
        color = \"#555555\",
        fontname = \"Helvetica\",
        fontsize = 10]

  n1  [label = \"(1) DATA IMPORT & STRUCTURE\"]
  n2  [label = \"(2) TRACE METAL SCREENING\"]
  n3  [label = \"(3) FILTERED DATASET (BC_filtered)\"]
  n4  [label = \"(4) DESCRIPTIVE STATISTICS\"]
  n5  [label = \"(5) CORRELATION ANALYSIS\"]
  n6  [label = \"(6) NMDS (Gower)\"]
  n7  [label = \"(7) MANTEL TESTS\"]
  n8  [label = \"(8) CRITERIA DEFINITION\"]
  n9  [label = \"(9) ASI WEIGHTING & SCORING\"]
  n10 [label = \"(10) RANKING & OUTPUT\"]
  
  # Main vertical path
  n1  -> n2
  n2  -> n3
  n3  -> n4
  
  # Branch from descriptive statistics to correlations and NMDS
  n4  -> n5
  n4  -> n6
  
  # Correlation and NMDS feed into Mantel tests
  n5  -> n7
  n6  -> n7

  # Mantel tests feed into criteria/AHP blocks
  n7  -> n8
  n8  -> n9
  n9  -> n10

  { rank = same; n5; n6 }
}
") -> g

# Export flowchart as vector PDF to the figures folder
DiagrammeRsvg::export_svg(g) |>
  charToRaw() |>
  rsvg::rsvg_pdf(file.path(fig_dir, "biochar_flowchart_ASI_pipeline.pdf"))

############################################################
## 4. Data Import & Preprocessing
############################################################
input_file <- "BC_data.txt"

bc_raw <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(bc_raw) <- trimws(colnames(bc_raw))

# Standardize key column names and types
names(bc_raw)[1] <- "BC"
bc_raw$BC <- trimws(as.character(bc_raw$BC))
bc_raw$Sample <- trimws(as.character(bc_raw$Sample))
bc_raw$Material <- trimws(as.character(bc_raw$Material))

# Convert numeric measurement columns safely
char_cols <- c("BC", "Sample", "Material")
num_cols <- setdiff(names(bc_raw), char_cols)
bc_raw[num_cols] <- lapply(bc_raw[num_cols], function(x) suppressWarnings(as.numeric(x)))

# Basic integrity checks
stopifnot(!anyDuplicated(bc_raw$BC))
stopifnot(nrow(bc_raw) == 36)

############################################################
## 5. Helper Functions
############################################################

# Min-Max Normalization function
normalize_min_max <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(0.5, length(x)))
  if (higher_better) (x - rng[1]) / (rng[2] - rng[1]) else (rng[2] - x) / (rng[2] - rng[1])
}

# Standard Error calculation
standard_error <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x) / sqrt(length(x))
}

# Descending rank helper
rank_desc <- function(x) rank(-x, ties.method = "average", na.last = "keep")

# Formatted text table writer for reports
write_table_txt <- function(df, file, title = NULL, digits = 4) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  if (!is.null(title)) {
    writeLines(title, con)
    writeLines(paste(rep("=", nchar(title)), collapse = ""), con)
    writeLines("", con)
  }
  df2 <- as.data.frame(df) # Convert tibble to data frame to support row.names = FALSE
  num_cols <- sapply(df2, is.numeric)
  df2[num_cols] <- lapply(df2[num_cols], round, digits = digits)
  capture.output(print(df2, row.names = FALSE), file = con)
}

# Unique key verification
assert_unique_key <- function(df, key = "BC", obj_name = deparse(substitute(df))) {
  if (anyDuplicated(df[[key]]) > 0) {
    stop(paste("Duplicate key values found in", obj_name, "for key", key))
  }
}

############################################################
## 6. Stage 1: Regulatory Trace Element Screening (set critical values)
############################################################
metal_limits <- c(Pb = 120, Cd = 1.0, Cu = 100, Ni = 30, Hg = 1.0, Zn = 400, Cr = 90, As = 13)
metal_cols <- names(metal_limits)

exceed_matrix <- sapply(metal_cols, function(m) bc_raw[[m]] > metal_limits[m])
exceed_matrix <- as.data.frame(exceed_matrix)

excluded_flag <- apply(exceed_matrix, 1, any, na.rm = TRUE)

excluded_table <- bc_raw %>%
  mutate(
    Exceeded = apply(exceed_matrix, 1, function(z) {
      bad <- names(which(z))
      if (length(bad) == 0) "" else paste(bad, collapse = ", ")
    })
  ) %>%
  filter(excluded_flag) %>%
  select(BC, Sample, Material, Temp, Time, HTSI, all_of(metal_cols), Exceeded)

# Filtered working dataset excluding non-compliant samples
bc <- bc_raw %>%
  filter(!excluded_flag) %>%
  arrange(BC)

assert_unique_key(bc)

# Output screening results to tables folder
write.csv(excluded_table, file.path(tab_dir, "BC_excluded_samples.csv"), row.names = FALSE)
write_table_txt(excluded_table, file.path(tab_dir, "BC_excluded_samples.txt"),
                "Excluded biochar samples exceeding trace-element thresholds")

############################################################
## 7. Stage 2A: Spearman Correlograms by Material Group
############################################################
corr_vars <- c(
  "Time","Temp","HTSI","BSI","HMC","NUTC","pH","EC","WHC","CEC",
  "N_%","S_%","P","K","Mg","Ca","Fe","Pb","Cd","Cu","Ni","Hg","Zn","Cr","As",
  "C_%","H_%","O_%"
)
corr_vars <- intersect(corr_vars, names(bc))

# Custom p-value matrix function for correlation tests
cor_pmat <- function(mat, method = "spearman") {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA_real_, n, n)
  diag(p.mat) <- 0
  colnames(p.mat) <- colnames(mat)
  rownames(p.mat) <- colnames(mat)
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i < j) {
        tmp <- suppressWarnings(cor.test(mat[, i], mat[, j], method = method, exact = FALSE))
        p.mat[i, j] <- tmp$p.value
        p.mat[j, i] <- tmp$p.value
      }
    }
  }
  p.mat
}

for (mat in sort(unique(bc$Material))) {
  dfm <- bc %>%
    filter(Material == mat) %>%
    select(all_of(corr_vars)) %>%
    drop_na()
  
  if (nrow(dfm) >= 3 && ncol(dfm) >= 2) {
    cor_mat <- cor(dfm, method = "spearman")
    p_mat   <- cor_pmat(dfm, method = "spearman")
    
    ord <- corrplot::corrMatOrder(cor_mat, order = "hclust")
    cor_mat_ord <- cor_mat[ord, ord]
    p_mat_ord   <- p_mat[ord, ord]
    
    # Save correlogram plot to the figures folder
    pdf(file.path(fig_dir, paste0("BC_correlogram_", mat, ".pdf")), width = 12, height = 10)
    corrplot(
      cor_mat_ord,
      method = "color",
      type = "upper",
      tl.col = "black",
      tl.cex = 0.7,
      addCoef.col = "black",
      number.cex = 0.5,
      p.mat = p_mat_ord,
      sig.level = 0.05,
      insig = "blank"
    )
    title(main = paste("Spearman correlogram -", mat))
    dev.off()
    
    corr_long <- expand.grid(
      Var1 = rownames(cor_mat),
      Var2 = colnames(cor_mat),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        spearman_r = as.vector(cor_mat),
        p_value    = as.vector(p_mat)
      )
    
    # Save correlation text table to the tables folder
    write_table_txt(
      corr_long,
      file.path(tab_dir, paste0("BC_spearman_correlation_", mat, ".txt")),
      title = paste("Spearman correlation results -", mat),
      digits = 4
    )
  }
}

############################################################
## 8. Summary Statistics Calculation
############################################################
summary_vars <- c(
  "HTSI","BSI","HMC","NUTC","pH","EC","WHC","CEC",
  "N_%","S_%","P","K","Mg","Ca","Fe","Pb","Cd","Cu","Ni","Hg","Zn","Cr","As",
  "C_%","H_%","O_%"
)
summary_vars <- intersect(summary_vars, names(bc))

summary_table <- bc %>%
  group_by(Material) %>%
  summarise(
    across(all_of(summary_vars),
           list(mean = ~mean(.x, na.rm = TRUE), se = ~standard_error(.x)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

# Export summary statistics to the tables folder
write.csv(summary_table, file.path(tab_dir, "BC_summary_statistics_by_material.csv"), row.names = FALSE)
write_table_txt(summary_table, file.path(tab_dir, "BC_summary_statistics_by_material.txt"),
                "Basic statistics by Material (mean and standard error)", 4)

############################################################
## 9. Multivariate Analysis: NMDS (Gower Distance)
############################################################
nmds_vars <- c(
  "pH","EC","WHC","CEC","N_%","S_%","P","K","S","Mg","Ca",
  "Fe","Pb","Cd","Cu","Ni","Hg","Zn","Cr","As","C_%","H_%","O_%"
)
nmds_vars <- intersect(nmds_vars, names(bc))

nmds_input <- bc %>%
  select(BC, Material, Temp, Time, HTSI, BSI, all_of(nmds_vars)) %>%
  drop_na()

if (nrow(nmds_input) < 3) {
  stop("Not enough complete cases for NMDS.")
}

gower_dist <- cluster::daisy(
  nmds_input %>% select(all_of(nmds_vars)),
  metric = "gower"
)

set.seed(123)
nmds_res <- metaMDS(
  gower_dist,
  k = 2,
  try = 50,
  trymax = 100,
  autotransform = FALSE
)
nmds_stress <- nmds_res$stress

site_scores <- as.data.frame(scores(nmds_res, display = "sites")) %>%
  tibble::rownames_to_column("row_id") %>%
  mutate(row_id = as.integer(row_id)) %>%
  left_join(
    nmds_input %>%
      mutate(row_id = dplyr::row_number()) %>%
      select(row_id, BC, Material, Temp, Time, HTSI, BSI),
    by = "row_id"
  ) %>%
  select(-row_id)

envfit_res <- envfit(
  nmds_res,
  nmds_input %>% select(Temp, Time),
  permutations = 999
)

vec_scores <- as.data.frame(scores(envfit_res, display = "vectors")) %>%
  rownames_to_column("Var")

p_nmds <- ggplot(site_scores, aes(NMDS1, NMDS2)) +
  geom_point(
    aes(shape = Material, fill = BSI),
    size = 3.5,
    stroke = 0.9,
    color = "black"
  ) +
  scale_shape_manual(values = c(DV = 21, RC = 23)) +
  scale_fill_gradient(name = "HTSI", low = "yellow", high = "blue") +
  geom_segment(
    data = vec_scores,
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.22, "cm")),
    color = "black"
  ) +
  geom_text(
    data = vec_scores,
    aes(x = NMDS1, y = NMDS2, label = Var),
    inherit.aes = FALSE,
    color = "black",
    vjust = -0.5
  ) +
  coord_equal() +
  theme_bw() +
  labs(
    title = paste0("FENIX Biochar Samples (Gower NMDS, stress = ", round(nmds_stress, 4), ")"),
    x = "NMDS1",
    y = "NMDS2"
  )

# Save NMDS plot to the figures folder
ggsave(
  file.path(fig_dir, "BC_NMDS_material.pdf"),
  p_nmds,
  width = 8,
  height = 6,
  device = cairo_pdf
)

# Save NMDS results to the tables folder
write_table_txt(
  tibble(Stress = nmds_stress),
  file.path(tab_dir, "BC_NMDS_results.txt"),
  "NMDS stress (Gower distance)",
  6
)

############################################################
## 10. Mantel Tests by Material Group
############################################################
mantel_vars <- c(
  "pH","EC","WHC","CEC","N_%","S_%","P","K","S","Mg","Ca",
  "Fe","Pb","Cd","Cu","Ni","Hg","Zn","Cr","As","C_%","H_%","O_%"
)
mantel_vars <- intersect(mantel_vars, names(bc))

run_mantel_gower <- function(df_group, response_col, meas_vars, mat_name) {
  needed_cols <- c(response_col, meas_vars)
  needed_cols <- intersect(needed_cols, names(df_group))
  
  dsub <- df_group %>%
    dplyr::select(all_of(needed_cols)) %>%
    drop_na()
  
  if (!(response_col %in% names(dsub)) || nrow(dsub) < 3 || length(meas_vars) < 2) {
    return(tibble(
      Material = mat_name,
      Variable = response_col,
      n = nrow(dsub),
      mantel_r = NA_real_,
      p_value = NA_real_
    ))
  }
  
  meas_present <- intersect(meas_vars, names(dsub))
  if (length(meas_present) < 2) {
    return(tibble(
      Material = mat_name,
      Variable = response_col,
      n = nrow(dsub),
      mantel_r = NA_real_,
      p_value = NA_real_
    ))
  }
  
  dist_resp <- cluster::daisy(data.frame(resp = dsub[[response_col]]), metric = "gower")
  dist_meas <- cluster::daisy(dsub[, meas_present, drop = FALSE], metric = "gower")
  
  res <- vegan::mantel(
    as.dist(dist_resp),
    as.dist(dist_meas),
    method = "pearson",
    permutations = 999
  )
  
  tibble(
    Material = mat_name,
    Variable = response_col,
    n = nrow(dsub),
    mantel_r = unname(res$statistic),
    p_value = res$signif
  )
}

materials <- sort(unique(bc$Material))
mantel_list <- list()

for (mat in materials) {
  df_group <- bc %>% filter(Material == mat)
  
  mantel_list[[paste0(mat, "_BSI")]]  <- run_mantel_gower(df_group, "BSI",  mantel_vars, mat)
  mantel_list[[paste0(mat, "_HTSI")]] <- run_mantel_gower(df_group, "HTSI", mantel_vars, mat)
}

mantel_results <- dplyr::bind_rows(mantel_list)

if (!all(c("Material","Variable","n","mantel_r","p_value") %in% names(mantel_results))) {
  stop("Mantel results table was not created correctly. Check group inputs and required columns.")
}

mantel_results <- mantel_results %>%
  arrange(Material, Variable)

# Export Mantel test results to the tables folder
write.csv(
  mantel_results,
  file.path(tab_dir, "BC_Mantel_results_by_material.csv"),
  row.names = FALSE
)

write_table_txt(
  mantel_results,
  file.path(tab_dir, "BC_Mantel_results_by_material.txt"),
  title = "Mantel test results by Material for BSI and HTSI",
  digits = 4
)

############################################################
## 11. Multi-Criteria and Aggregate Suitability Index
############################################################

# Piecewise scoring functions and class assignment for criteria
score_WHC <- function(x) case_when(
  is.na(x) ~ NA_real_, x > 2 ~ 1.00, x >= 1.5 ~ 0.80, x >= 1.0 ~ 0.55, x >= 0 ~ 0.20, TRUE ~ 0.00
)
class_WHC <- function(x) case_when(
  is.na(x) ~ NA_character_, x > 2 ~ "excellent", x >= 1.5 ~ "very_good", x >= 1.0 ~ "good", x >= 0 ~ "low", TRUE ~ "outside_range"
)

score_pH <- function(x) case_when(
  is.na(x) ~ NA_real_, x >= 9.5 & x <= 12.5 ~ 1.00, x >= 7.5 & x < 9.5 ~ 0.70, x >= 4.0 & x < 7.5 ~ 0.20, TRUE ~ 0.00
)
class_pH <- function(x) case_when(
  is.na(x) ~ NA_character_, x >= 9.5 & x <= 12.5 ~ "excellent", x >= 7.5 & x < 9.5 ~ "good", x >= 4.0 & x < 7.5 ~ "low", TRUE ~ "outside_range"
)

score_EC <- function(x) case_when(
  is.na(x) ~ NA_real_, x >= 2 & x <= 4 ~ 1.00, x < 2 ~ 0.80, x > 4 & x <= 8 ~ 0.40, x > 8 ~ 0.00, TRUE ~ 0.00
)
class_EC <- function(x) case_when(
  is.na(x) ~ NA_character_, x >= 2 & x <= 4 ~ "very_good", x < 2 ~ "good", x > 4 & x <= 8 ~ "moderate", x > 8 ~ "avoid", TRUE ~ "outside_range"
)

score_CEC <- function(x) case_when(
  is.na(x) ~ NA_real_, x > 30 ~ 1.00, x >= 15 & x <= 29 ~ 0.70, x >= 10 & x <= 14 ~ 0.40, x >= 0 ~ 0.15, TRUE ~ 0.00
)
class_CEC <- function(x) case_when(
  is.na(x) ~ NA_character_, x > 30 ~ "excellent", x >= 15 & x <= 29 ~ "good", x >= 10 & x <= 14 ~ "moderate", x >= 0 ~ "low", TRUE ~ "outside_range"
)

score_HMC <- function(x) case_when(
  is.na(x) ~ NA_real_, x >= 0.3 & x <= 0.7 ~ 1.00, x < 0.3 ~ 0.55, x > 0.7 ~ 0.10, TRUE ~ 0.00
)
class_HMC <- function(x) case_when(
  is.na(x) ~ NA_character_, x < 0.3 ~ "single_risk", x >= 0.3 & x <= 0.7 ~ "balanced", x > 0.7 ~ "elevated_risk", TRUE ~ "outside_range"
)

score_NUTC <- function(x) case_when(
  is.na(x) ~ NA_real_, x >= 0.6 & x <= 0.8 ~ 1.00, x >= 0.3 & x < 0.6 ~ 0.70, x > 0.8 ~ 0.35, x < 0.3 ~ 0.15, TRUE ~ 0.00
)
class_NUTC <- function(x) case_when(
  is.na(x) ~ NA_character_, x < 0.3 ~ "deficiency", x >= 0.3 & x < 0.6 ~ "acceptable", x >= 0.6 & x <= 0.8 ~ "optimum", x > 0.8 ~ "limiting", TRUE ~ "outside_range"
)

score_BSI <- function(x) case_when(
  is.na(x) ~ NA_real_, x >= 70 ~ 1.00, x >= 50 ~ 0.70, x >= 30 ~ 0.40, x >= 0 ~ 0.15, TRUE ~ 0.00
)
class_BSI <- function(x) case_when(
  is.na(x) ~ NA_character_, x >= 70 ~ "excellent", x >= 50 ~ "good", x >= 30 ~ "moderate", x >= 0 ~ "low", TRUE ~ "outside_range"
)

bc_ASI <- bc %>%
  mutate(
    HMC_norm_raw  = normalize_min_max(HMC, TRUE),
    NUTC_norm_raw = normalize_min_max(NUTC, TRUE),
    WHC_class  = class_WHC(WHC),
    pH_class   = class_pH(pH),
    EC_class   = class_EC(EC),
    CEC_class  = class_CEC(CEC),
    HMC_class  = class_HMC(HMC_norm_raw),
    NUTC_class = class_NUTC(NUTC_norm_raw),
    BSI_class  = class_BSI(BSI)
  ) %>%
  select(BC, Sample, Material, Temp, Time, HTSI,
         WHC, pH, EC, CEC, HMC, NUTC, BSI,
         WHC_class, pH_class, EC_class, CEC_class, HMC_class, NUTC_class, BSI_class,
         HMC_norm_raw, NUTC_norm_raw)

assert_unique_key(bc_ASI)

criteria_scores <- bc_ASI %>%
  transmute(
    BC,
    WHC  = score_WHC(WHC),
    pH   = score_pH(pH),
    EC   = score_EC(EC),
    CEC  = score_CEC(CEC),
    HMC  = score_HMC(HMC_norm_raw),
    NUTC = score_NUTC(NUTC_norm_raw),
    BSI  = score_BSI(BSI)
  )

assert_unique_key(criteria_scores)

criteria_ranks <- criteria_scores %>%
  transmute(
    BC,
    WHC_rank  = rank_desc(WHC),
    pH_rank   = rank_desc(pH),
    EC_rank   = rank_desc(EC),
    CEC_rank  = rank_desc(CEC),
    HMC_rank  = rank_desc(HMC),
    NUTC_rank = rank_desc(NUTC),
    BSI_rank  = rank_desc(BSI)
  )

criteria_norm <- criteria_scores %>%
  mutate(across(-BC, ~ normalize_min_max(.x, higher_better = TRUE)))

# Define ASI criteria importance by order
# Adjust importance accordingly
criteria_weight_rank <- c(pH = 1, EC = 1, BSI = 2, WHC = 3, CEC = 4, NUTC = 5, HMC = 6)
w_raw <- 1 / criteria_weight_rank
criteria_weights <- w_raw / sum(w_raw)

criteria_rank_table <- tibble(
  Criterion = names(criteria_weights),
  Importance_rank = unname(criteria_weight_rank[names(criteria_weights)]),
  Weight = as.numeric(criteria_weights)
) %>% arrange(Importance_rank, Criterion)

# Export criteria weight table to the tables folder
write.csv(criteria_rank_table, file.path(tab_dir, "BC_ASI_criteria_rank_table.csv"), row.names = FALSE)
write_table_txt(criteria_rank_table, file.path(tab_dir, "BC_ASI_criteria_rank_table.txt"),
                "ASI criteria rank table", 6)

score_mat <- as.matrix(criteria_norm %>% select(all_of(names(criteria_weights))))
ahp_tbl <- tibble(
  BC = criteria_norm$BC,
  ASI_score = as.vector(score_mat %*% criteria_weights)
)

assert_unique_key(ahp_tbl)

############################################################
## 12. Final Ranking & Output Generation
############################################################

all_samples_table <- bc_ASI %>%
  left_join(ahp_tbl, by = "BC") %>%
  left_join(criteria_ranks, by = "BC") %>%
  arrange(desc(ASI_score)) %>%
  mutate(ASI_rank = row_number()) %>%
  select(
    BC, Sample, ASI_rank, ASI_score,
    Material, Temp, Time, HTSI,
    WHC, pH, EC, CEC, HMC, NUTC, BSI, BSI_class,
    WHC_class, pH_class, EC_class, CEC_class, HMC_class, NUTC_class,
    WHC_rank, pH_rank, EC_rank, CEC_rank, HMC_rank, NUTC_rank, BSI_rank
  )

stopifnot(nrow(all_samples_table) == nrow(bc_ASI))

top5_table <- all_samples_table %>%
  slice(1:5) %>%
  select(
    BC, Sample, ASI_rank, ASI_score,
    Material, Temp, Time, HTSI,
    WHC, pH, EC, CEC, HMC, NUTC, BSI, BSI_class
  )

# Export final ranked output tables and CSVs to the tables folder
write.csv(all_samples_table, file.path(tab_dir, "BC_all_ASI_samples_raw_values.csv"), row.names = FALSE)
write.csv(top5_table, file.path(tab_dir, "BC_top5_ASI_samples_raw_values.csv"), row.names = FALSE)
write_table_txt(all_samples_table, file.path(tab_dir, "BC_all_ASI_samples_raw_values.txt"),
                "All ASI-scored biochar samples", 4)
write_table_txt(top5_table, file.path(tab_dir, "BC_top5_ASI_samples_raw_values.txt"),
                "Top 5 ASI biochar samples", 4)

cat("\nAnalysis completed successfully.\n")
cat("Figures saved in:", fig_dir, "\n")
cat("Tables and CSVs saved in:", tab_dir, "\n")
cat("Remaining samples after exclusion:", nrow(bc), "\n")
cat("NMDS stress:", round(nmds_stress, 6), "\n")
