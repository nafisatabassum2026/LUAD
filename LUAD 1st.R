getwd()
dir.create("data/raw", recursive = TRUE)
dir.create("data/processed", recursive = TRUE)
dir.create("scripts", recursive = TRUE)
dir.create("results/figures", recursive = TRUE)
dir.create("results/tables", recursive = TRUE)
file.create("data/raw/.gitkeep")
file.create("data/processed/.gitkeep")
file.create("results/figures/.gitkeep")
file.create("results/tables/.gitkeep")
writeLines(c(
  "TCGA-LUAD (PanCancer Atlas) - Mutation counts from cBioPortal",
  "Checked on: 23 August,2026",
  "",
  "EGFR:",
  "  Total mutations: 86",
  "  Driver: 69",
  "  VUS: 17",
  "",
  "KRAS:",
  "  Total mutations: 173",
  "  Driver: 171",
  "  VUS: 2"
), "data/group_sizes_note_cbioportal.txt")
getwd()
list.files("data/raw")
install.packages("data.table")
library(data.table)
install.packages('R.utils')
# RNA-seq counts (still log2(count+1) — don't transform yet, just look at it)
expr <- fread("data/raw/TCGA-LUAD.star_counts.tsv.gz")

# Gene ID → gene symbol mapping
probemap <- fread("data/raw/gencode.v36.annotation.gtf.gene.probemap")

# Mutation data (ensemble WXS calls)
mut <- fread("data/raw/TCGA-LUAD.somaticmutation_wxs.tsv.gz")

# Clinical/phenotype data
clin <- fread("data/raw/TCGA-LUAD.clinical.tsv.gz")

# Survival data
surv <- fread("data/raw/TCGA-LUAD.survival.tsv.gz")
dim(expr)
dim(probemap)
dim(mut)
dim(clin)
dim(surv)
# Expression: first column should be gene ID, rest are sample columns
expr[1:5, 1:5]

# Probemap: should show gene ID → gene symbol mapping
head(probemap)

# Mutation: THIS is the one we need to understand carefully
colnames(mut)
head(mut)

# Clinical
colnames(clin)
head(clin, 3)


# Survival
colnames(surv)
head(surv, 3)
library(dplyr)

# Look at what values "effect" actually contains, so we filter correctly
table(mut$effect)

# Look at what "callers" looks like
head(mut$callers, 20)

library(dplyr)
library(stringr)

# Step 1: Define which "effect" categories count as functionally significant
functional_effects <- c(
  "missense_variant", "frameshift_variant", "stop_gained", "stop_lost",
  "start_lost", "inframe_deletion", "inframe_insertion",
  "splice_acceptor_variant", "splice_donor_variant"
)

# Step 2: Mark functionally significant mutations
mut <- mut %>%
  mutate(is_functional = str_detect(effect, paste(functional_effects, collapse = "|")))

# Step 3: Mark high-confidence calls (2+ callers agree)
mut <- mut %>%
  mutate(n_callers = str_count(callers, ";") + 1,
         is_confident = n_callers >= 2)

# Step 4: Filter for EGFR and KRAS
egfr_mut <- mut %>% filter(gene == "EGFR", is_functional, is_confident)
kras_mut <- mut %>% filter(gene == "KRAS", is_functional, is_confident)

# Step 5: Get unique patient sample IDs
egfr_samples <- unique(egfr_mut$sample)
kras_samples <- unique(kras_mut$sample)

length(egfr_samples)
length(kras_samples)

# All tumor samples from the expression data
# Barcode format: TCGA-XX-XXXX-01A → the "01" right after the 3rd hyphen = Primary Tumor
all_samples <- colnames(expr)[-1]   # remove the first column (Ensembl_ID)
sample_type <- substr(all_samples, 14, 15)
all_tumor_samples <- all_samples[sample_type == "01"]

length(all_tumor_samples)   # sanity check: should be close to ~535 or so

# Wild-type = tumor samples with NEITHER EGFR nor KRAS mutation
# Note: mutation file sample IDs may have a slightly different suffix than expression
# file sample IDs (e.g. "-01A" vs "-01") — we handle that by trimming to 15 characters
egfr_samples_trim <- substr(egfr_samples, 1, 15)
kras_samples_trim <- substr(kras_samples, 1, 15)
all_tumor_trim    <- substr(all_tumor_samples, 1, 15)

wt_samples <- all_tumor_trim[!(all_tumor_trim %in% c(egfr_samples_trim, kras_samples_trim))]

length(wt_samples)

mutation_groups <- list(
  egfr = egfr_samples_trim,
  kras = kras_samples_trim,
  wt   = wt_samples
)

saveRDS(mutation_groups, "data/processed/mutation_groups.rds")

# Also save a simple text summary for your Methods chapter
writeLines(c(
  paste("EGFR-mutant samples:", length(mutation_groups$egfr)),
  paste("KRAS-mutant samples:", length(mutation_groups$kras)),
  paste("Wild-type samples:", length(mutation_groups$wt))
), "data/group_sizes_note.txt")

mutation_groups <- list(
  egfr = egfr_samples_trim,
  kras = kras_samples_trim,
  wt   = wt_samples
)

saveRDS(mutation_groups, "data/processed/mutation_groups.rds")

writeLines(c(
  "TCGA-LUAD mutation-based patient stratification",
  "Method: functional mutation filter (missense/frameshift/stop/splice-site variants,",
  "detected by >=2 independent variant callers)",
  "",
  paste("EGFR-mutant samples:", length(mutation_groups$egfr)),
  paste("KRAS-mutant samples:", length(mutation_groups$kras)),
  paste("Wild-type samples:", length(mutation_groups$wt)),
  paste("Total tumor samples:", length(mutation_groups$egfr) + length(mutation_groups$kras) + length(mutation_groups$wt)),
  "",
  "Cross-check against cBioPortal (TCGA-LUAD PanCancer Atlas, driver mutations only):",
  "EGFR driver: 69 | KRAS driver: 171"
), "data/group_sizes_note.txt")

library(dplyr)

# Add a group label to the clinical data using your mutation_groups list
clin <- clin %>%
  mutate(
    sample_trim = substr(submitter_id.samples, 1, 15),
    mutation_group = case_when(
      sample_trim %in% mutation_groups$egfr ~ "EGFR",
      sample_trim %in% mutation_groups$kras ~ "KRAS",
      sample_trim %in% mutation_groups$wt   ~ "WT",
      TRUE ~ NA_character_
    )
  )

# Check sex distribution across groups (EGFR-mutant should skew more female)
table(clin$mutation_group, clin$gender.demographic)

# Check smoking history across groups (EGFR-mutant should skew more never/light-smoker)
clin %>%
  filter(!is.na(mutation_group)) %>%
  group_by(mutation_group) %>%
  summarise(mean_cigarettes_per_day = mean(cigarettes_per_day.exposures, na.rm = TRUE))

library(dplyr)

clin <- clin %>%
  mutate(
    sample_trim = substr(sample, 1, 15),
    mutation_group = case_when(
      sample_trim %in% mutation_groups$egfr ~ "EGFR",
      sample_trim %in% mutation_groups$kras ~ "KRAS",
      sample_trim %in% mutation_groups$wt   ~ "WT",
      TRUE ~ NA_character_
    )
  )

# Check sex distribution across groups
table(clin$mutation_group, clin$gender.demographic)

# Check smoking history across groups
clin %>%
  filter(!is.na(mutation_group)) %>%
  group_by(mutation_group) %>%
  summarise(mean_cigarettes_per_day = mean(cigarettes_per_day.exposures, na.rm = TRUE))
table(clin$mutation_group, clin$tobacco_smoking_history)
grep("smoking", colnames(clin), value = TRUE, ignore.case = TRUE)
grep("smok|tobacco|pack_years|cigarette", colnames(clin), value = TRUE, ignore.case = TRUE)
clin %>%
  filter(!is.na(mutation_group)) %>%
  group_by(mutation_group) %>%
  summarise(
    n = n(),
    n_never_smoker_proxy = sum(is.na(pack_years_smoked.exposures)),
    pct_never_smoker_proxy = round(100 * n_never_smoker_proxy / n, 1)
  )
writeLines(c(
  "Clinical characterization of mutation-based groups (sanity check)",
  "",
  "Sex distribution (% female):",
  "  EGFR: 68% | KRAS: 53% | WT: 51%",
  "",
  "Never-smoker proxy (NA in pack_years_smoked.exposures):",
  "  EGFR: 63% | KRAS: 25.5% | WT: 30.4%",
  "",
  "Both findings are consistent with established LUAD literature",
  "(EGFR-mutant tumors more common in female, never-smoker patients),",
  "supporting the validity of the mutation-based stratification."
), "data/group_sizes_note.txt")



library(dplyr)

# Convert to matrix, gene IDs as row names
expr_mat <- as.matrix(expr[, -1, with = FALSE])   # drop the Ensembl_ID column
rownames(expr_mat) <- expr$Ensembl_ID

dim(expr_mat)   # sanity check — should be 60660 genes x 589 samples

expr_mat <- as.matrix(expr[, -1])
rownames(expr_mat) <- expr$Ensembl_ID

dim(expr_mat)

expr_mat_raw <- round(2^expr_mat - 1)

# Sanity check — values should now look like real counts (whole numbers, often large)
expr_mat_raw[1:5, 1:5]
colnames(expr_mat_raw) <- substr(colnames(expr_mat_raw), 1, 15)

# Check for duplicate column names after trimming
sum(duplicated(colnames(expr_mat_raw)))
# Find which column positions are duplicates (excluding the first occurrence)
dup_cols <- duplicated(colnames(expr_mat_raw))

# Keep only non-duplicate columns
expr_mat_raw <- expr_mat_raw[, !dup_cols]

# Confirm
dim(expr_mat_raw)
sum(duplicated(colnames(expr_mat_raw)))   # should now be 0


library(dplyr)

sample_ids <- colnames(expr_mat_raw)

group_label <- case_when(
  sample_ids %in% mutation_groups$egfr ~ "EGFR",
  sample_ids %in% mutation_groups$kras ~ "KRAS",
  sample_ids %in% mutation_groups$wt   ~ "WT",
  TRUE ~ NA_character_
)

table(group_label, useNA = "always")
keep_samples <- !is.na(group_label)

expr_final <- expr_mat_raw[, keep_samples]
group_final <- group_label[keep_samples]

dim(expr_final)        # should be 60660 genes x 516 samples (65+139+312)
table(group_final)      # should match 65 / 139 / 312

library(DESeq2)

coldata <- data.frame(
  sample = colnames(expr_final),
  group = factor(group_final, levels = c("WT", "EGFR", "KRAS"))  # WT = reference/baseline
)
rownames(coldata) <- coldata$sample

dds <- DESeqDataSetFromMatrix(
  countData = expr_final,
  colData = coldata,
  design = ~ group
)

# Filter out very low-expression genes (standard practice, speeds things up too)
dds <- dds[rowSums(counts(dds)) > 10, ]

dim(dds)   # check how many genes remain after filtering


library(DESeq2)

# Step 3a: Build the sample metadata table
coldata <- data.frame(
  sample = colnames(expr_final),
  group = factor(group_final, levels = c("WT", "EGFR", "KRAS"))
)
rownames(coldata) <- coldata$sample

# Step 3b: Create the DESeq2 dataset object — THIS is what creates "dds"
dds <- DESeqDataSetFromMatrix(
  countData = expr_final,
  colData = coldata,
  design = ~ group
)

# Step 3c: NOW filtering works, because dds exists
dds <- dds[rowSums(counts(dds)) > 10, ]

dim(dds)

dds <- DESeq(dds)

resultsNames(dds)


res_egfr <- results(dds, contrast = c("group", "EGFR", "WT"))
summary(res_egfr)

res_egfr_ordered <- res_egfr[order(res_egfr$padj), ]
head(res_egfr_ordered, 10)


res_kras <- results(dds, contrast = c("group", "KRAS", "WT"))
summary(res_kras)

res_kras_ordered <- res_kras[order(res_kras$padj), ]
head(res_kras_ordered, 10)


library(dplyr)

# Clean up probemap: strip version suffix (.17, .6 etc.) to match expr's IDs if needed
# First check: does res_egfr's rownames format match probemap$id exactly?
head(rownames(res_egfr))
head(probemap$id)
library(dplyr)

res_egfr_df <- as.data.frame(res_egfr) %>%
  tibble::rownames_to_column("Ensembl_ID") %>%
  left_join(probemap %>% select(id, gene), by = c("Ensembl_ID" = "id")) %>%
  arrange(padj)

res_kras_df <- as.data.frame(res_kras) %>%
  tibble::rownames_to_column("Ensembl_ID") %>%
  left_join(probemap %>% select(id, gene), by = c("Ensembl_ID" = "id")) %>%
  arrange(padj)

head(res_egfr_df, 10)
head(res_kras_df, 10)

# Save the full DESeq2 object (so you never have to re-run DESeq2 from scratch)
saveRDS(dds, "data/processed/dds_object.rds")

# Save annotated results as CSV — these are actual thesis result tables
write.csv(res_egfr_df, "results/tables/DEG_EGFR_vs_WT.csv", row.names = FALSE)
write.csv(res_kras_df, "results/tables/DEG_KRAS_vs_WT.csv", row.names = FALSE)
