## === User parameters ===
study_name <- "example"
condition_1 <- "Controls"
condition_2 <- "Treated"
organism <- "Homo sapiens" # "Homo sapiens" or "Mus musculus"
input_pattern <- "\\.tabular$"
cond1_dir <- "input/cond1"
cond2_dir <- "input/cond2"
outdir <- "output"

# Optional: specify additional covariates file (CSV with samples as rows)
covariates_file <- NULL  # e.g., "metadata/covariates.csv" or NULL if none

# Create output directories
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(file.path(outdir, "plots"))) dir.create(file.path(outdir, "plots"), recursive = TRUE)
if (!dir.exists(file.path(outdir, "dea_results"))) dir.create(file.path(outdir, "dea_results"), recursive = TRUE)

set.seed(123)

## === Packages ===
pkgs <- c("edgeR", "limma", "sva", "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db", "ggplot2")
for (p in pkgs) {
  if (!suppressWarnings(requireNamespace(p, quietly = TRUE))) {
    message("Package '", p, "' not installed. Please install before running the pipeline.")
  }
}
library(edgeR)
library(limma)
library(sva)
library(AnnotationDbi)
library(ggplot2)

if (organism == "Mus musculus") {
  library(org.Mm.eg.db)
  orgdb <- org.Mm.eg.db
} else if (organism == "Homo sapiens") {
  library(org.Hs.eg.db)
  orgdb <- org.Hs.eg.db
} else {
  stop("organism must be either 'Mus musculus' or 'Homo sapiens'")
}

## === Helper functions ===
load_counts_from_files <- function(files, pattern = NULL) {
  if (length(files) == 0) stop("No files provided.")
  mats <- lapply(files, function(f) {
    df <- read.table(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    genes <- as.character(df[[1]])
    counts <- as.numeric(df[[2]])
    sample_name <- if (!is.null(pattern)) sub(pattern, "", basename(f)) else basename(f)
    list(genes = genes, counts = counts, sample = sample_name)
  })
  genes_ref <- mats[[1]]$genes
  for (m in mats[-1]) {
    if (!identical(genes_ref, m$genes)) stop("Gene IDs/order mismatch across files.")
  }
  mat <- do.call(cbind, lapply(mats, function(m) m$counts))
  colnames(mat) <- sapply(mats, function(m) m$sample)
  rownames(mat) <- genes_ref
  return(mat)
}

detect_id_type <- function(ids) {
  if (all(grepl("^[0-9]+$", ids))) return("ENTREZID")
  if (all(grepl("^ENS", ids))) return("ENSEMBL")
  return("OTHER")
}

annotate_results <- function(df, idtype = "ENTREZID", orgdb) {
  keys <- rownames(df)
  keytype <- ifelse(idtype == "OTHER", "SYMBOL", idtype)
  out <- tryCatch(
    AnnotationDbi::select(orgdb, keys = keys, columns = c("SYMBOL", "ENTREZID", "GENENAME"), keytype = keytype),
    error = function(e) NULL
  )
  if (is.null(out)) return(df)
  out <- out[!duplicated(out[[keytype]]), ]
  rownames(out) <- out[[keytype]]
  ann <- out[rownames(df), c("SYMBOL", "ENTREZID", "GENENAME")]
  df$SYMBOL <- ann$SYMBOL
  df$ENTREZID <- ann$ENTREZID
  df$GENENAME <- ann$GENENAME
  return(df)
}

## === Load sample files ===
cond1_files <- list.files(cond1_dir, pattern = input_pattern, full.names = TRUE)
cond2_files <- list.files(cond2_dir, pattern = input_pattern, full.names = TRUE)
all_files <- c(cond1_files, cond2_files)

message("Loading count tables...")
counts_mat <- load_counts_from_files(all_files, pattern = "\\.tabular$")

# Ensure group factor levels are stable and in desired order
group <- factor(c(rep("ctr", length(cond1_files)), rep("case", length(cond2_files))), levels = c("ctr", "case"))

dge <- DGEList(counts = counts_mat, group = group)

## === Filtering lowly expressed genes ===
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
#message(sprintf("Kept %d genes after filtering.", nrow(dge)))

dge <- calcNormFactors(dge)

## === Load covariates if provided ===
covariates <- NULL
if (!is.null(covariates_file) && file.exists(covariates_file)) {
  covariates <- read.csv(covariates_file, row.names = 1, check.names = FALSE)
  # keep only samples that are present in dge (columns = samples)
  covariates <- covariates[colnames(dge), , drop = FALSE]
  message("Loaded covariates: ", paste(colnames(covariates), collapse = ", "))
}

## === Design matrix ===
# base design (no intercept) with explicit levels
design_base <- model.matrix(~0 + group)
colnames(design_base) <- levels(group)

if (!is.null(covariates) && ncol(covariates) > 0) {
  # ensure covariates are in a data.frame with proper column names
  df_cov <- as.data.frame(covariates, stringsAsFactors = FALSE)
  design <- model.matrix(~0 + group + ., data = data.frame(group = group, df_cov))
} else {
  design <- design_base
}

# Remove empty / NA column names and make syntactically valid names
bad_cols <- which(is.na(colnames(design)) | colnames(design) == "" | grepl('^\\s+$', colnames(design)))
if (length(bad_cols) > 0) {
  design <- design[, -bad_cols, drop = FALSE]
}
colnames(design) <- make.names(colnames(design))

## === voom ===
v <- voom(dge, design = design, plot = FALSE)

## === Optional SVA ===
mod <- design
# null model should be intercept only (number of rows = samples)
mod0 <- model.matrix(~1, data = data.frame(group = group))

n.sv_est <- tryCatch(num.sv(v$E, mod, method = "leek"), error = function(e) 0)
if (n.sv_est > 0) {
  svobj <- sva(v$E, mod, mod0, n.sv = n.sv_est)
  # add surrogate variables with sensible column names
  sv_mat <- svobj$sv
  colnames(sv_mat) <- paste0("SV", seq_len(ncol(sv_mat)))
  design <- cbind(design, sv_mat)
  message("Added ", n.sv_est, " surrogate variables to design.")
}

# Clean again after adding SVs
colnames(design) <- make.names(colnames(design))

## === Fit model ===
fit <- lmFit(v, design)

# Ensure contrast uses the cleaned design column names
# Check that 'case' and 'ctr' are present
if (!all(c("case", "ctr") %in% colnames(design))) {
  stop("'case' and 'ctr' columns must be present in the design matrix. Current columns: ", paste(colnames(design), collapse = ", "))
}

contrast.matrix <- makeContrasts(case - ctr, levels = colnames(design))
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res <- topTable(fit2, number = Inf, sort.by = "P")
res$FDR <- p.adjust(res$P.Value, method = "BH")
idtype <- detect_id_type(rownames(res))
res_annot <- annotate_results(res, idtype = idtype, orgdb = orgdb)

## === Save only raw results ===
write.csv(res_annot, file = file.path(outdir, "dea_results", paste0(study_name, "_raw.csv")), row.names = TRUE)

## === QC plots ===
mds <- plotMDS(dge, plot = FALSE)
mds.df <- data.frame(MDS1 = mds$x, MDS2 = mds$y, group = group, sample = colnames(dge))
p <- ggplot(mds.df, aes(x = MDS1, y = MDS2, color = group, label = sample)) +
  geom_point(size = 3) + geom_text(vjust = -1, size = 3) + theme_bw() +
  ggtitle(paste0(study_name, " - MDS plot"))
ggsave(filename = file.path(outdir, "plots", paste0(study_name, "_MDSplot.pdf")), plot = p, width = 6, height = 6)

writeLines(capture.output(sessionInfo()), con = file.path(outdir, paste0(study_name, "_sessionInfo.txt")))
message("Pipeline finished. Raw results saved to ", file.path(outdir, "dea_results"))

