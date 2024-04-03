# run examples from localVariance() function documentation
library(SpotSweeper)
library(SpatialExperiment)

# load example data
spe <- STexampleData::Visium_humanDLPFC()

# change from gene id to gene names
rownames(spe) <- rowData(spe)$gene_name

# show column data before SpotSweepR
colnames(colData(spe))

# drop out-of-tissue spots
spe <- spe[, spe$in_tissue == 1]
spe <- spe[, !is.na(spe$ground_truth)]

# Identifying the mitochondrial transcripts in our SpatialExperiment.
is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]

# Calculating QC metrics for each spot using scuttle
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
colnames(colData(spe))

spe <- localVariance(spe,
    features = "subsets_Mito_percent",
    n_neighbors = 36,
    name = "local_mito_variance_k36"
)


# === Tests ===
test_that("example objects have correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("examples give correct number of colData", {
  expect_equal(length(colnames(colData(spe))), 15)
})

