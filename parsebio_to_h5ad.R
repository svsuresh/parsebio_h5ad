library(fastMatMR)
library(Matrix)
library(anndata)

data_dir <- "/Users/suresh/Desktop/1M_PBMC_T1D_Parse_bio"
genes_file <- file.path(data_dir, "all_genes_1M_PBMC.csv")
metadata_file <- file.path(data_dir, "cell_metadata_1M_PBMC.csv")
mtx_file <- file.path(data_dir, "DGE_1M_PBMC.mtx")
output_file <- file.path(data_dir, "1M_PBMC_T1D_parsebio.h5ad")

genes <- read.csv(genes_file)
metadata <- read.csv(metadata_file)
mat <- as(fmm_to_sparse_Matrix(mtx_file), "CsparseMatrix")

cat(sprintf("Matrix: %d x %d\n", nrow(mat), ncol(mat)))
cat(sprintf("Genes: %d\n", nrow(genes)))
cat(sprintf("Cells: %d\n", nrow(metadata)))

# Orient to cells x genes
if (nrow(mat) == nrow(genes) && ncol(mat) == nrow(metadata)) {
  mat <- t(mat)
}

rownames(mat) <- metadata$bc_wells
colnames(mat) <- genes$gene_id

obs <- metadata[, colnames(metadata) != "bc_wells", drop = FALSE]
rownames(obs) <- metadata$bc_wells

var <- genes
rownames(var) <- genes$gene_id

adata <- AnnData(X = mat, obs = obs, var = var)
write_h5ad(adata, output_file)
cat(sprintf("Saved: %s\n", output_file))
