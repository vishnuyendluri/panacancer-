## =============================================================================
## REG4 PAN-CANCER STATE-DOMINANCE ANALYSIS  —  complete pipeline (Steps 1–10)
## Data: TCGA PanCanAtlas (UCSC Xena): EB++ batch-corrected RNA-seq, TCGA-CDR
##       clinical/survival, ABSOLUTE purity, PanCanAtlas molecular subtypes.
##
## HOW TO RUN (on IU Quartz, login node — internet needed for downloads):
##   module load r
##   mkdir -p ~/reg4_pancancer && cd ~/reg4_pancancer
##   Rscript reg4_full_pipeline.R > pipeline_log.txt 2>&1
##   (or inside R:  source("reg4_full_pipeline.R", echo = TRUE))
##
## OUTPUT: results/*.csv (tables), results/figures/*.pdf (figures),
##         data/ (cached downloads), pipeline_log.txt (all printed results)
## Runtime: a few minutes (downloads are cached after the first run).
## =============================================================================

## ---------------------------------------------------------------------------
## 0. SETUP: packages (installed to personal library), folders, gene sets
## ---------------------------------------------------------------------------
Sys.setenv(TZ = "America/New_York")
options(timeout = 600)
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
pkgs <- c("UCSCXenaTools", "data.table", "survival", "ppcor", "ggplot2", "patchwork")
miss <- setdiff(pkgs, rownames(installed.packages()))
if (length(miss)) install.packages(miss, lib = lib, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(UCSCXenaTools); library(data.table); library(survival)
  library(ppcor); library(ggplot2); library(patchwork)
})
dir.create("data", showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

sets <- list(
  secretory        = c("MUC2","TFF3","SPINK4","FCGBP","AGR2"),
  epithelial       = c("KRT8","KRT18","KRT19","EPCAM"),
  gi_lineage       = c("CDX2","VIL1","KLF4","HNF4A"),
  stromal_emt      = c("COL1A1","COL1A2","COL3A1","COL5A1","FN1","SPARC","POSTN","FAP",
                       "PDGFRB","THY1","DCN","LUM","VIM","ZEB1","SNAI2","TWIST1"),
  immune_cytotoxic = c("PTPRC","CD2","CD3D","CD3E","CD8A","LCK","GZMA","GZMB","PRF1",
                       "IFNG","CXCL9","CXCL10"),
  nk               = c("KLRD1","KLRF1","NCR1","NCR3","SH2D1B","KIR2DL4","KIR3DL1","NCAM1"))
genes <- unique(c("REG4", unlist(sets)))            # 50 genes; REG4 is in NO score
sets$epi_secretory <- c(sets$secretory, sets$epithelial)

gi     <- c("COAD","READ","PAAD","STAD","ESCA_EAC","CHOL")
nonepi <- c("LAML","DLBC","THYM","SARC","SKCM","UVM","GBM","LGG","PCPG","TGCT","MESO","UCS")
gi5    <- c("READ","COAD","STAD","PAAD","ESCA_EAC")
grp_col <- c("GI adeno" = "#b2182b", "Other epithelial" = "#4393c3", "Non-epithelial" = "grey55")
na_ct  <- list(estimate = NA_real_, p.value = NA_real_)
zscore <- function(x) { s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x))); (x - mean(x, na.rm = TRUE)) / s }

## ---------------------------------------------------------------------------
## STEP 1. Expression: fetch only the 50 needed genes (not the 1.9 GB matrix)
## ---------------------------------------------------------------------------
f_expr <- "data/expr_subset.rds"
if (!file.exists(f_expr)) {
  m <- fetch_dense_values("https://pancanatlas.xenahubs.net",
         "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena",
         identifiers = genes, time_limit = 600)
  saveRDS(m, f_expr)
}
m <- readRDS(f_expr)
cat("\n== STEP 1 ==\n"); print(dim(m)); print(setdiff(genes, rownames(m)))
print(summary(as.vector(m)))
if (max(m, na.rm = TRUE) > 50) m <- log2(m + 1)      # data is already log2; safety check
na <- rowSums(is.na(m)); print(na[na > 0])

## ---------------------------------------------------------------------------
## STEP 2. Clinical (TCGA-CDR), purity (ABSOLUTE), sample QC, cohorts/groups
## ---------------------------------------------------------------------------
expr <- as.data.table(t(m), keep.rownames = "sample")
expr <- expr[substr(sample, 14, 15) %in% c("01","03")]   # primary tumours (03 = LAML)
expr[, patient := substr(sample, 1, 12)]
expr <- expr[!duplicated(patient)]                       # one tumour per patient

cdr <- XenaGenerate(subset = XenaHostNames == "pancanAtlasHub") |>
  XenaFilter(filterDatasets = "Survival_SupplementalTable_S1_20171025") |>
  XenaQuery() |> XenaDownload(destdir = "data") |> XenaPrepare()
clin <- as.data.table(cdr)
setnames(clin, c("cancer type abbreviation","age_at_initial_pathologic_diagnosis",
                 "ajcc_pathologic_tumor_stage","histological_type"),
         c("cancer","age","stage","histology"), skip_absent = TRUE)
clin <- clin[, .(sample, cancer, age = as.numeric(age), stage, histology,
                 OS, OS.time, DSS, DSS.time, PFI, PFI.time)]
dat <- merge(expr, clin, by = "sample")

pur_file <- "data/absolute_purity.txt"
if (!file.exists(pur_file))
  download.file("https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
                pur_file, mode = "wb")
pur <- fread(pur_file)[, .(sample = substr(array, 1, 15), purity = as.numeric(purity))]
dat <- merge(dat, pur[!duplicated(sample)], by = "sample", all.x = TRUE)

dat[, cohort := cancer]
dat[cancer == "ESCA", cohort := ifelse(grepl("adeno", histology, ignore.case = TRUE),
                                       "ESCA_EAC", "ESCA_ESCC")]
dat[, group := fcase(cohort %in% gi, "GI adeno", cohort %in% nonepi, "Non-epithelial",
                     default = "Other epithelial")]
cat("\n== STEP 2 ==\n"); print(nrow(dat)); print(mean(!is.na(dat$purity)))
print(dat[, .N, by = .(group, cohort)][order(group, -N)], nrows = 40)
print(dat[cancer == "ESCA", .N, by = .(cohort, histology)])

## ---------------------------------------------------------------------------
## STEP 3. Scores (mean of within-cohort z-scores) + REG4 detection
## ---------------------------------------------------------------------------
for (s in names(sets)) {
  g <- intersect(sets[[s]], names(dat))
  dat[, (paste0("S_", s)) := rowMeans(sapply(.SD, zscore), na.rm = TRUE),
      by = cohort, .SDcols = g]
}
dat[, REG4z := as.numeric(scale(REG4)), by = cohort]
dat[, SECz  := as.numeric(scale(S_secretory)), by = cohort]
st <- toupper(dat$stage)
st <- ifelse(grepl("^STAGE [IV]+", st), sub("^STAGE ([IV]+).*$", "\\1", st), NA)
dat[, stage2 := fifelse(st %in% c("I","II"), 0, fifelse(st %in% c("III","IV"), 1, NA_real_))]

detect <- dat[, .(n = .N, REG4_median = round(median(REG4, na.rm = TRUE), 2),
                  REG4_detect = round(mean(REG4 > 1, na.rm = TRUE), 2)),
              by = .(group, cohort)][order(-REG4_median)]
fwrite(detect, "results/00_REG4_detection_by_cohort.csv")
cat("\n== STEP 3 ==\n"); print(detect, nrows = 40)
print(round(cor(dat[, .(S_epi_secretory, S_stromal_emt, S_immune_cytotoxic, S_nk)],
                method = "spearman", use = "pairwise"), 2))

## ---------------------------------------------------------------------------
## STEP 4 (Q1). REG4 vs each marker gene, per cohort (raw + purity-adjusted)
## ---------------------------------------------------------------------------
q1_genes <- c(sets$secretory, sets$epithelial, sets$gi_lineage)
rows <- list()
for (co in unique(dat$cohort)) {
  d <- dat[cohort == co]
  for (g in q1_genes) {
    x <- d$REG4; y <- d[[g]]; ok <- !is.na(x) & !is.na(y); ok2 <- ok & !is.na(d$purity)
    r  <- tryCatch(suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)),
                   error = function(e) na_ct)
    pr <- tryCatch(suppressWarnings(pcor.test(x[ok2], y[ok2], d$purity[ok2], method = "spearman")$estimate),
                   error = function(e) NA_real_)
    rows[[length(rows) + 1]] <- data.table(group = d$group[1], cohort = co, gene = g,
      rho = unname(r$estimate), p = r$p.value, rho_purity = unname(pr), n = sum(ok))
  }
}
q1 <- rbindlist(rows)
q1[, fdr := p.adjust(p, "BH")]
q1 <- merge(q1, detect[, .(cohort, REG4_detect)], by = "cohort")
fwrite(q1, "results/01_Q1_gene_correlations.csv")

summ <- q1[REG4_detect >= 0.2, .(
  secretory     = median(rho[gene %in% sets$secretory], na.rm = TRUE),
  epithelial    = median(rho[gene %in% sets$epithelial], na.rm = TRUE),
  gi_lineage    = median(rho[gene %in% sets$gi_lineage], na.rm = TRUE),
  secretory_pur = median(rho_purity[gene %in% sets$secretory], na.rm = TRUE)),
  by = .(group, cohort)][order(group, -secretory)]
num <- c("secretory","epithelial","gi_lineage","secretory_pur")
summ[, (num) := lapply(.SD, round, 2), .SDcols = num]
gi_wide <- dcast(q1[group == "GI adeno"], gene ~ cohort, value.var = "rho")
gi_wide[, (names(gi_wide)[-1]) := lapply(.SD, round, 2), .SDcols = -1]
cat("\n== STEP 4 (Q1) ==\n"); print(q1[n < 10, .(cohort, gene, n)])
print(summ, nrows = 40); print(gi_wide[match(q1_genes, gene)])

## ---------------------------------------------------------------------------
## STEP 5 (Q3 + key analysis). Axis correlations + dominant axis per cohort
## ---------------------------------------------------------------------------
axes <- c("S_epi_secretory","S_secretory","S_stromal_emt","S_immune_cytotoxic","S_nk")
rows <- list()
for (co in unique(dat$cohort)) {
  d <- dat[cohort == co]; x <- d$REG4
  for (a in axes) {
    y <- d[[a]]; ok <- !is.na(x) & !is.na(y); ok2 <- ok & !is.na(d$purity)
    r  <- tryCatch(suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)), error = function(e) na_ct)
    pc <- tryCatch(suppressWarnings(pcor.test(x[ok2], y[ok2], d$purity[ok2], method = "spearman")), error = function(e) na_ct)
    rows[[length(rows) + 1]] <- data.table(group = d$group[1], cohort = co, axis = sub("^S_", "", a),
      rho = unname(r$estimate), p = r$p.value, rho_pur = unname(pc$estimate), p_pur = pc$p.value)
  }
  ok <- complete.cases(x, d$S_nk, d$S_immune_cytotoxic)          # NK beyond immune
  pc <- tryCatch(suppressWarnings(pcor.test(x[ok], d$S_nk[ok], d$S_immune_cytotoxic[ok], method = "spearman")),
                 error = function(e) na_ct)
  rows[[length(rows) + 1]] <- data.table(group = d$group[1], cohort = co, axis = "nk_given_immune",
    rho = unname(pc$estimate), p = pc$p.value, rho_pur = NA_real_, p_pur = NA_real_)
}
ax <- rbindlist(rows)
ax[, fdr := p.adjust(p, "BH"), by = axis][, fdr_pur := p.adjust(p_pur, "BH"), by = axis]
ax <- merge(ax, detect[, .(cohort, REG4_detect)], by = "cohort")
fwrite(ax, "results/02_axis_correlations.csv")

## Rule: top |rho| >= 0.3, FDR < 0.05, and beats runner-up by >= 0.1 (3 axes)
call_dom <- function(tab, rcol, qcol, epi = "epi_secretory") {
  t3 <- tab[axis %in% c(epi, "stromal_emt", "immune_cytotoxic")]
  t3[, r := get(rcol)][, q := get(qcol)]
  t3[order(-abs(r)), .(top = axis[1], top_rho = round(r[1], 2), margin = round(abs(r[1]) - abs(r[2]), 2),
     dominant = if (!is.na(r[1]) && q[1] < 0.05 && abs(r[1]) >= 0.3 && abs(r[1]) - abs(r[2]) >= 0.1) axis[1] else "none"),
     by = .(group, cohort, REG4_detect)]
}
dom     <- call_dom(ax, "rho", "fdr")
dom_pur <- call_dom(ax, "rho_pur", "fdr_pur")
dom_sec <- call_dom(ax, "rho", "fdr", epi = "secretory")
dom[, dom_purity  := dom_pur$dominant[match(cohort, dom_pur$cohort)]]
dom[, dom_secOnly := dom_sec$dominant[match(cohort, dom_sec$cohort)]]
dom[REG4_detect < 0.2, c("dominant", "dom_purity", "dom_secOnly") := "not expressed"]
nk <- ax[axis == "nk_given_immune", .(cohort, nk_partial = round(rho, 2), nk_fdr = signif(fdr, 2))]
dom <- merge(dom, nk, by = "cohort")
fwrite(dom, "results/02_REG4_dominant_axis.csv")
cat("\n== STEP 5 (dominance) ==\n")
print(dom[order(group, -abs(top_rho)), .(group, cohort, top, top_rho, margin, dominant,
                                         dom_purity, dom_secOnly, nk_partial)], nrows = 40)

## ---------------------------------------------------------------------------
## STEP 6 (Q2). GI vs non-GI enrichment tests + Figure 1 (dominance heatmap)
## ---------------------------------------------------------------------------
dom[, epi_dom := dominant == "epi_secretory"][, GI := group == "GI adeno"]
expr_dom <- dom[REG4_detect >= 0.2]
e <- ax[axis == "epi_secretory" & REG4_detect >= 0.2]
f_all  <- fisher.test(table(dom$GI, dom$epi_dom))
f_expr <- fisher.test(table(expr_dom$GI, expr_dom$epi_dom))
f_pur  <- fisher.test(table(expr_dom$GI, expr_dom$dom_purity == "epi_secretory"))
w      <- wilcox.test(rho ~ (group == "GI adeno"), data = e)
q2 <- data.table(test = c("Fisher all","Fisher expressed","Fisher purity-adj","Wilcoxon rho"),
                 OR = round(c(f_all$estimate, f_expr$estimate, f_pur$estimate, NA), 1),
                 p  = signif(c(f_all$p.value, f_expr$p.value, f_pur$p.value, w$p.value), 2))
fwrite(q2, "results/02_Q2_enrichment_tests.csv")
cat("\n== STEP 6 (Q2) ==\n"); print(q2)
print(e[, .(median_rho = round(median(rho), 2), n = .N), by = group])
print(q1[cohort == "PRAD", .(gene, rho = round(rho, 2))])
print(ax[axis == "nk_given_immune" & fdr < 0.05, .(cohort, rho = round(rho, 2), fdr = signif(fdr, 2))])

## Figure 1: dominance heatmap (grey = REG4 not expressed; * = purity-robust call)
lv <- rev(dom[order(group, REG4_detect < 0.2, -top_rho)]$cohort)
hm <- ax[axis %in% c("epi_secretory","stromal_emt","immune_cytotoxic","nk")]
hm[, axis := factor(axis, levels = c("epi_secretory","stromal_emt","immune_cytotoxic","nk"),
     labels = c("Epithelial/\nsecretory","Stromal/\nEMT","Immune/\ncytotoxic","NK"))]
hm[, expressed := REG4_detect >= 0.2]
hm[, label := ifelse(fdr < 0.05 & expressed, sprintf("%.2f", rho), "")]
hm[, cohort := factor(cohort, levels = lv)]
calls <- dom[, .(cohort = factor(cohort, levels = lv), group,
  call = factor(dominant, levels = c("epi_secretory","immune_cytotoxic","stromal_emt","none","not expressed"),
                labels = c("Epithelial/secretory","Immune","Stromal","None","Not expressed")),
  robust = ifelse(dom_purity == dominant & !dominant %in% c("none","not expressed"), "*", ""))]
p1 <- ggplot(hm, aes(axis, cohort)) +
  geom_tile(aes(fill = ifelse(expressed, rho, NA)), colour = "white") +
  geom_text(aes(label = label), size = 2.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1),
                       na.value = "grey88", name = "Spearman rho") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 9) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"),
        panel.grid = element_blank(), legend.position = "bottom")
p2 <- ggplot(calls, aes(x = "Dominant\naxis", y = cohort, fill = call)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = robust), size = 4, colour = "white") +
  scale_fill_manual(values = c("Epithelial/secretory" = "#b2182b", "Immune" = "#7b3294",
    "Stromal" = "#1b7837", "None" = "grey70", "Not expressed" = "grey92"), name = NULL, drop = FALSE) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 9) +
  theme(axis.text.y = element_blank(), strip.text = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom") + guides(fill = guide_legend(ncol = 1))
ggsave("results/figures/Fig1B_REG4_dominance_heatmap.pdf",
       p1 + p2 + plot_layout(widths = c(4, 0.6)), width = 7, height = 9.5)

## ---------------------------------------------------------------------------
## STEP 7 (Q4). Mucinous histology
## ---------------------------------------------------------------------------
focus <- c("COAD","READ","STAD","ESCA_EAC","PAAD","CHOL","LUAD","OV","UCEC","BRCA","BLCA")
dat[, mucinous := grepl("mucin|signet|colloid", histology, ignore.case = TRUE)]
ht <- dat[cohort %in% focus, .(n = .N, REG4_med = round(median(REG4), 2),
          sec_med = round(median(S_secretory), 2)), by = .(cohort, histology)][order(cohort, -REG4_med)]
mt <- dat[cohort %in% focus, {
  nm <- sum(mucinous)
  if (nm >= 3 && nm < .N) {
    .(n_muc = nm, n_other = .N - nm,
      REG4_muc = round(median(REG4[mucinous]), 2), REG4_other = round(median(REG4[!mucinous]), 2),
      p_REG4 = signif(wilcox.test(REG4 ~ mucinous)$p.value, 2),
      sec_muc = round(median(S_secretory[mucinous]), 2), sec_other = round(median(S_secretory[!mucinous]), 2),
      p_sec = signif(wilcox.test(S_secretory ~ mucinous)$p.value, 2))
  } else .(n_muc = nm, n_other = .N - nm, REG4_muc = NA_real_, REG4_other = NA_real_, p_REG4 = NA_real_,
           sec_muc = NA_real_, sec_other = NA_real_, p_sec = NA_real_)
}, by = cohort]
nm_cor <- dat[cohort %in% focus & !mucinous, .(n = .N,
  rho_nonmucinous = round(cor(REG4, S_secretory, method = "spearman", use = "complete.obs"), 2)), by = cohort]
fwrite(ht, "results/03_Q4_histology_table.csv"); fwrite(mt, "results/03_Q4_mucinous_tests.csv")
fwrite(nm_cor, "results/03_Q4_nonmucinous_correlations.csv")
cat("\n== STEP 7 (Q4) ==\n"); print(ht, nrows = 80); print(mt); print(nm_cor)

## ---------------------------------------------------------------------------
## STEP 8 (Q5). Does REG4 add prognostic information beyond the secretory state?
## ---------------------------------------------------------------------------
run_cox <- function(d, endpoint) {
  d <- copy(d)
  d[, time := get(paste0(endpoint, ".time"))][, event := as.numeric(get(endpoint))]
  covs <- c("S_epi_secretory", "age")
  if (mean(!is.na(d$stage2)) >= 0.7) covs <- c(covs, "stage2")
  d <- d[complete.cases(d[, c("time","event","REG4z", covs), with = FALSE]) & time > 0]
  if (nrow(d) < 50 || sum(d$event) < 15) return(NULL)
  uni <- summary(coxph(Surv(time, event) ~ REG4z, data = d))$coefficients
  f0 <- as.formula(paste("Surv(time, event) ~", paste(covs, collapse = " + ")))
  m0 <- coxph(f0, data = d); m1 <- coxph(update(f0, . ~ . + REG4z), data = d)
  s1 <- summary(m1)$coefficients
  data.table(endpoint, n = nrow(d), events = sum(d$event), stage_adj = "stage2" %in% covs,
    HR_uni = round(exp(uni["REG4z","coef"]), 2), p_uni = signif(uni["REG4z","Pr(>|z|)"], 2),
    HR_adj = round(exp(s1["REG4z","coef"]), 2), p_adj = signif(s1["REG4z","Pr(>|z|)"], 2),
    HR_epi = round(exp(s1["S_epi_secretory","coef"]), 2), p_epi = signif(s1["S_epi_secretory","Pr(>|z|)"], 2),
    p_LRT = signif(anova(m0, m1)[["Pr(>|Chi|)"]][2], 2),
    C_base = round(unname(m0$concordance["concordance"]), 3),
    C_full = round(unname(m1$concordance["concordance"]), 3))
}
cox_cohorts <- c(gi5, "BLCA", "PRAD")
surv <- rbindlist(lapply(cox_cohorts, function(co) rbindlist(lapply(c("PFI","OS"), function(e) {
  r <- run_cox(dat[cohort == co], e); if (!is.null(r)) r[, cohort := co]; r }))), fill = TRUE)
setcolorder(surv, "cohort")
surv[, fdr_LRT := signif(p.adjust(p_LRT, "BH"), 2)]
fwrite(surv, "results/04_Q5_survival.csv")

gid <- dat[cohort %in% gi5 & !is.na(PFI.time) & PFI.time > 0 & !is.na(age)]
pm0 <- coxph(Surv(PFI.time, PFI) ~ S_epi_secretory + age + strata(cohort), data = gid)
pm1 <- update(pm0, . ~ . + REG4z)
cat("\n== STEP 8 (Q5) ==\n"); print(surv)
print(round(summary(pm1)$coefficients, 3)); print(anova(pm0, pm1))

## ---------------------------------------------------------------------------
## STEP 9 (Q6). REG4 x secretory four-state analysis
## ---------------------------------------------------------------------------
q6_cohorts <- c("COAD","READ","PAAD","STAD","ESCA_EAC","BLCA","PRAD")
dat[, reg4_hi := REG4 > median(REG4, na.rm = TRUE), by = cohort]
dat[, sec_hi  := S_secretory > median(S_secretory, na.rm = TRUE), by = cohort]
dat[, group4 := factor(paste0(ifelse(reg4_hi, "R+", "R-"), ifelse(sec_hi, "S+", "S-")),
                       levels = c("R-S-","R+S-","R-S+","R+S+"))]
sizes <- dcast(dat[cohort %in% q6_cohorts, .N, by = .(cohort, group4)], cohort ~ group4, value.var = "N")
feat <- dat[cohort %in% q6_cohorts, .(immune = round(median(S_immune_cytotoxic), 2),
  nk = round(median(S_nk), 2), stromal = round(median(S_stromal_emt), 2),
  purity = round(median(purity, na.rm = TRUE), 2)), by = .(cohort, group4)][order(cohort, group4)]
res <- rbindlist(lapply(q6_cohorts, function(co) {
  d <- dat[cohort == co & !is.na(purity)]
  rbindlist(lapply(c("S_immune_cytotoxic","S_nk","S_stromal_emt"), function(y) {
    cf <- summary(lm(as.formula(paste(y, "~ REG4z * SECz + purity")), data = d))$coefficients
    data.table(cohort = co, outcome = sub("^S_", "", y),
      beta_REG4 = round(cf["REG4z","Estimate"], 3), p_REG4 = signif(cf["REG4z","Pr(>|t|)"], 2),
      beta_int = round(cf["REG4z:SECz","Estimate"], 3), p_int = signif(cf["REG4z:SECz","Pr(>|t|)"], 2))
  }))
}))
res[, fdr_REG4 := signif(p.adjust(p_REG4, "BH"), 2)][, fdr_int := signif(p.adjust(p_int, "BH"), 2)]
s4 <- rbindlist(lapply(q6_cohorts, function(co) {
  d <- dat[cohort == co & !is.na(PFI.time) & PFI.time > 0 & !is.na(age)]
  if (sum(d$PFI) < 15 || min(table(d$group4)) < 5) return(NULL)
  cf <- summary(coxph(Surv(PFI.time, PFI) ~ group4 + age, data = d))$coefficients
  g <- grep("^group4", rownames(cf))
  data.table(cohort = co, contrast = sub("group4", "", rownames(cf)[g]),
             HR = round(exp(cf[g, "coef"]), 2), p = signif(cf[g, "Pr(>|z|)"], 2))
}))
fwrite(sizes, "results/05_Q6_group_sizes.csv"); fwrite(feat, "results/05_Q6_group_features.csv")
fwrite(res, "results/05_Q6_interaction.csv");   fwrite(s4, "results/05_Q6_survival_groups.csv")
cat("\n== STEP 9 (Q6) ==\n"); print(sizes); print(feat, nrows = 40); print(res, nrows = 30); print(s4, nrows = 30)

## ---------------------------------------------------------------------------
## STEP 9b. Molecular subtype check (MSI / EBV) for the colorectal immune link
## ---------------------------------------------------------------------------
sub <- XenaGenerate(subset = XenaHostNames == "pancanAtlasHub") |>
  XenaFilter(filterDatasets = "TCGASubtype") |> XenaQuery() |> XenaDownload(destdir = "data") |> XenaPrepare()
sub <- as.data.table(sub); setnames(sub, 1, "sampleID")
sub <- sub[, .(sample = substr(sampleID, 1, 15), gi_subtype = Subtype_Selected)]
dat <- merge(dat, sub[!duplicated(sample)], by = "sample", all.x = TRUE)

sub_tab <- dat[cohort %in% c("COAD","READ","STAD"), .(n = .N, REG4 = round(median(REG4), 2),
    secretory = round(median(S_secretory), 2), immune = round(median(S_immune_cytotoxic), 2)),
    by = .(cohort, gi_subtype)][order(cohort, gi_subtype)]
msi_chk <- rbindlist(lapply(c("COAD","READ","STAD"), function(co) {
  d <- dat[cohort == co & !is.na(purity) & !is.na(gi_subtype)]
  rbindlist(lapply(c("S_immune_cytotoxic","S_nk"), function(y) {
    b0 <- summary(lm(as.formula(paste(y, "~ REG4z + SECz + purity")), data = d))$coefficients
    b1 <- summary(lm(as.formula(paste(y, "~ REG4z + SECz + purity + gi_subtype")), data = d))$coefficients
    data.table(cohort = co, outcome = sub("^S_", "", y), n = nrow(d),
      beta_before = round(b0["REG4z","Estimate"], 3), p_before = signif(b0["REG4z","Pr(>|t|)"], 2),
      beta_after = round(b1["REG4z","Estimate"], 3), p_after = signif(b1["REG4z","Pr(>|t|)"], 2))
  }))
}))
fwrite(sub_tab, "results/06_subtype_medians.csv"); fwrite(msi_chk, "results/06_subtype_adjusted_immune.csv")
cat("\n== STEP 9b (subtype) ==\n"); print(sub_tab); print(msi_chk)
saveRDS(dat, "data/dat_final.rds")

## ---------------------------------------------------------------------------
## STEP 10. Supporting figures (Fig 1A, 2, 3, 4, S1–S2)
## ---------------------------------------------------------------------------
theme_set(theme_classic(base_size = 9) +
          theme(strip.background = element_blank(), strip.text = element_text(face = "bold")))

## Fig 1A: REG4 expression across 34 cancers
ord <- dat[, .(m = median(REG4, na.rm = TRUE)), by = cohort][order(-m)]$cohort
f1a <- ggplot(dat, aes(factor(cohort, levels = ord), REG4, fill = group)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = grp_col, name = NULL) +
  labs(x = NULL, y = "REG4 expression (log2)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top")
ggsave("results/figures/Fig1A_REG4_expression_pancancer.pdf", f1a, width = 8, height = 3.5)

## Fig 2: gene-level heatmap + REG4 vs secretory scatter
fig2_cohorts <- c(gi5, "CHOL", "BLCA", "PRAD", "LUAD", "BRCA", "ESCA_ESCC")
mod <- data.table(gene = q1_genes, module = rep(c("Secretory","Epithelial","GI lineage TF"), c(5, 4, 4)))
h2 <- merge(q1[cohort %in% fig2_cohorts], mod, by = "gene")
h2[, cohort := factor(cohort, levels = fig2_cohorts)]
h2[, gene := factor(gene, levels = rev(mod$gene))]
h2[, module := factor(module, levels = c("Secretory","Epithelial","GI lineage TF"))]
f2a <- ggplot(h2, aes(cohort, gene, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(fdr < 0.05, sprintf("%.2f", rho), "")), size = 2.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "Spearman rho") +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 6.5, linewidth = 0.5) +
  labs(x = NULL, y = NULL, subtitle = "GI adenocarcinomas  |  non-GI comparators") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank())
sc_cohorts <- c("READ","COAD","STAD","PAAD","PRAD","BRCA")
sc <- dat[cohort %in% sc_cohorts][, cohort := factor(cohort, levels = sc_cohorts)]
lab <- sc[, .(rho = cor(REG4, S_secretory, method = "spearman", use = "complete.obs")), by = cohort]
f2b <- ggplot(sc, aes(S_secretory, REG4)) +
  geom_point(aes(colour = group), size = 0.5, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.4) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = sprintf("rho = %.2f", rho)),
            hjust = -0.1, vjust = 1.5, size = 2.8) +
  scale_colour_manual(values = grp_col, guide = "none") +
  facet_wrap(~ cohort, nrow = 2, scales = "free") +
  labs(x = "Secretory score (MUC2, TFF3, SPINK4, FCGBP, AGR2)", y = "REG4 (log2)")
ggsave("results/figures/Fig2_REG4_secretory_program.pdf",
       f2a / f2b + plot_layout(heights = c(1.3, 1)) + plot_annotation(tag_levels = "A"),
       width = 7.5, height = 10)

## Fig 3: mucinous vs non-mucinous; breast negative control; non-mucinous correlations
m_cohorts <- c("COAD","READ","STAD","PAAD","BRCA")
md <- dat[cohort %in% m_cohorts][, cohort := factor(cohort, levels = m_cohorts)]
md[, histo := factor(ifelse(mucinous, "Mucinous", "Non-mucinous"), levels = c("Non-mucinous","Mucinous"))]
pv <- md[, .(pR = wilcox.test(REG4 ~ histo)$p.value, pS = wilcox.test(S_secretory ~ histo)$p.value), by = cohort]
pv[, labR := paste0("p = ", signif(pR, 2))][, labS := paste0("p = ", signif(pS, 2))]
hist_col <- c("Non-mucinous" = "grey75", "Mucinous" = "#d6604d")
box_mucin <- function(yvar, ylab, plab) {
  ggplot(md, aes(histo, .data[[yvar]], fill = histo)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.35) +
    geom_text(data = pv, aes(x = 1.5, y = Inf, label = .data[[plab]]), inherit.aes = FALSE, vjust = 1.3, size = 2.6) +
    facet_wrap(~ cohort, nrow = 1) +
    scale_fill_manual(values = hist_col, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = ylab) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
f3a <- box_mucin("REG4", "REG4 (log2)", "labR")
f3b <- box_mucin("S_secretory", "Secretory score", "labS")
nm <- dat[cohort %in% c(gi5, "CHOL","BLCA","OV","BRCA","LUAD","UCEC") & !mucinous,
          .(rho = cor(REG4, S_secretory, method = "spearman", use = "complete.obs")), by = .(cohort, group)]
f3c <- ggplot(nm, aes(reorder(cohort, rho), rho, fill = group)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = grp_col, name = NULL) +
  labs(x = NULL, y = "Spearman rho, REG4 vs secretory score\n(non-mucinous tumours only)") +
  theme(legend.position = "bottom")
ggsave("results/figures/Fig3_mucinous_and_breast_control.pdf",
       ((f3a / f3b) | f3c) + plot_layout(widths = c(3, 1)) + plot_annotation(tag_levels = "A"),
       width = 10, height = 6.5)

## Fig 4: molecular subtype (MSI / EBV) and subtype-adjusted immune effects
sd <- dat[cohort %in% c("COAD","READ","STAD") & !is.na(gi_subtype)]
sd[, tumour := ifelse(cohort == "STAD", "Gastric (STAD)", "Colorectal (COAD + READ)")]
sd[, subtype := factor(sub("^GI\\.", "", gi_subtype), levels = c("CIN","GS","HM-indel","HM-SNV","EBV"),
                       labels = c("CIN","GS","MSI (HM-indel)","POLE (HM-SNV)","EBV"))]
sub_col <- c("CIN" = "#4393c3", "GS" = "#92c5de", "MSI (HM-indel)" = "#d6604d",
             "POLE (HM-SNV)" = "#f4a582", "EBV" = "#5aae61")
box_sub <- function(yvar, ylab) {
  ggplot(sd, aes(subtype, .data[[yvar]], fill = subtype)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.35) +
    facet_grid(~ tumour, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = sub_col, guide = "none") +
    labs(x = NULL, y = ylab) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
f4a <- box_sub("REG4", "REG4 (log2)")
f4b <- box_sub("S_immune_cytotoxic", "Immune/cytotoxic score")
cf_rows <- list()
for (co in c("COAD","READ","STAD")) {
  d <- dat[cohort == co & !is.na(purity) & !is.na(gi_subtype)]
  for (y in c("S_immune_cytotoxic","S_nk")) {
    for (adj in c("Before", "After subtype adjustment")) {
      f <- paste(y, "~ REG4z + SECz + purity", if (adj != "Before") "+ gi_subtype" else "")
      mm <- lm(as.formula(f), data = d); ci <- confint(mm)["REG4z", ]
      cf_rows[[length(cf_rows) + 1]] <- data.table(cohort = co,
        outcome = ifelse(y == "S_nk", "NK", "Immune/cytotoxic"), model = adj,
        beta = unname(coef(mm)["REG4z"]), lo = unname(ci[1]), hi = unname(ci[2]))
    }
  }
}
cfd <- rbindlist(cf_rows)
cfd[, model := factor(model, levels = c("Before", "After subtype adjustment"))]
f4c <- ggplot(cfd, aes(cohort, beta, colour = model)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.3) +
  coord_flip() + facet_wrap(~ outcome) +
  scale_colour_manual(values = c("Before" = "black", "After subtype adjustment" = "#d6604d"), name = NULL) +
  labs(x = NULL, y = "REG4 effect (beta per SD)\nadjusted for secretory score + purity") +
  theme(legend.position = "bottom")
ggsave("results/figures/Fig4_molecular_subtype.pdf",
       (f4a / f4b / f4c) + plot_layout(heights = c(1, 1, 0.9)) + plot_annotation(tag_levels = "A"),
       width = 7, height = 10)

## Supp S1: survival forest plot (PFI); Supp S2: purity robustness
cox_ci <- function(d, endpoint = "PFI") {
  d <- copy(d)
  d[, time := get(paste0(endpoint, ".time"))][, event := as.numeric(get(endpoint))]
  covs <- c("S_epi_secretory", "age")
  if (mean(!is.na(d$stage2)) >= 0.7) covs <- c(covs, "stage2")
  d <- d[complete.cases(d[, c("time","event","REG4z", covs), with = FALSE]) & time > 0]
  if (nrow(d) < 50 || sum(d$event) < 15) return(NULL)
  mu <- coxph(Surv(time, event) ~ REG4z, data = d)
  ma <- coxph(as.formula(paste("Surv(time, event) ~ REG4z +", paste(covs, collapse = " + "))), data = d)
  rbind(
    data.table(model = "Unadjusted", HR = exp(coef(mu)[["REG4z"]]),
               lo = exp(confint(mu)["REG4z", 1]), hi = exp(confint(mu)["REG4z", 2])),
    data.table(model = "Adjusted for secretory state", HR = exp(coef(ma)[["REG4z"]]),
               lo = exp(confint(ma)["REG4z", 1]), hi = exp(confint(ma)["REG4z", 2])))
}
fs <- rbindlist(lapply(cox_cohorts, function(co) {
  r <- cox_ci(dat[cohort == co]); if (!is.null(r)) r[, cohort := co]; r }))
fs[, cohort := factor(cohort, levels = rev(cox_cohorts))]
fs[, model := factor(model, levels = c("Unadjusted", "Adjusted for secretory state"))]
s1 <- ggplot(fs, aes(cohort, HR, colour = model)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.3) +
  coord_flip() + scale_y_log10() +
  scale_colour_manual(values = c("Unadjusted" = "grey50", "Adjusted for secretory state" = "#b2182b"), name = NULL) +
  labs(x = NULL, y = "Hazard ratio for PFI per SD of REG4 (log scale)") +
  theme(legend.position = "bottom")
pr <- ax[axis == "epi_secretory" & REG4_detect >= 0.2]
s2 <- ggplot(pr, aes(rho, rho_pur, colour = group)) +
  geom_abline(linetype = 2, linewidth = 0.3) +
  geom_point(size = 1.5) +
  geom_text(aes(label = cohort), size = 2.3, vjust = -0.8, show.legend = FALSE) +
  scale_colour_manual(values = grp_col, name = NULL) +
  coord_equal(xlim = c(-0.35, 0.85), ylim = c(-0.35, 0.85)) +
  labs(x = "Spearman rho, REG4 vs epithelial/secretory score (raw)",
       y = "Partial rho (adjusted for ABSOLUTE purity)") +
  theme(legend.position = "bottom")
ggsave("results/figures/FigS1_S2_survival_and_purity.pdf",
       (s1 | s2) + plot_annotation(tag_levels = list(c("S1", "S2"))), width = 10, height = 4.8)

## ---------------------------------------------------------------------------
## DONE
## ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
cat("\n== DONE ==\nTables:\n");  print(list.files("results", pattern = "\\.csv$"))
cat("Figures:\n");               print(list.files("results/figures"))
