# run examples from localOutliers() function documentation
library(SpatialExperiment)

# load example data
spe <- STexampleData::Visium_humanDLPFC()

# change from gene id to gene names
rownames(spe) <- rowData(spe)$gene_name

# drop out-of-tissue spots
spe <- spe[, spe$in_tissue == 1]
spe <- spe[, !is.na(spe$ground_truth)]

# Identifying the mitochondrial transcripts in our SpatialExperiment.
is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]

# Calculating QC metrics for each spot using scuttle
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
colnames(colData(spe))

# Identifying local outliers using SpotSweeper
spe <- localOutliers(spe,
                     metric = "sum",
                     direction = "lower",
                     log = TRUE
)

spe <- localOutliers(spe,
                     metric = "detected",
                     direction = "lower",
                     log = TRUE
)

spe <- localOutliers(spe,
                     metric = "subsets_Mito_percent",
                     direction = "higher",
                     log = FALSE
)

# combine all outliers into "local_outliers" column
spe$local_outliers <- as.logical(spe$sum_outliers) |
  as.logical(spe$detected_outliers) |
  as.logical(spe$subsets_Mito_percent_outliers)



# === Tests ===
test_that("example objects have correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("correct outliers were found", {
  expect_equal(sum(as.logical(spe$sum_outliers)), 7)
  expect_equal(sum(as.logical(spe$detected_outliers)), 9)
  expect_equal(sum(as.logical(spe$subsets_Mito_percent_outliers)), 2)
  expect_equal(sum(as.logical(spe$local_outliers)), 11)
})
