#' localOutliers Function
#'
#' This function detects local outliers in spatial transcriptomics data based on
#' standard quality control metrics, such as library size, unique genes, and
#' mitochondrial ratio. Local outliers are defined as spots with low/high
#' quality metrics compared to their surrounding neighbors, based on a modified
#' z-score statistic.
#'
#' @param spe SpatialExperiment object
#' @param metric colData QC metric to use for outlier detection
#' @param direction Direction of outlier detection (higher, lower, or both)
#' @param n_neighbors Number of nearest neighbors to use for outlier detection
#' @param samples Column name in colData to use for sample IDs
#' @param log Logical indicating whether to log1p transform the features
#' (default is TRUE)
#' @param cutoff Cutoff for outlier detection (default is 3)
#'
#' @return SpatialExperiment object with updated colData containing outputs
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom spatialEco outliers
#'
#' @export localOutliers
#'
#' @examples
#' library(SpotSweeper)
#' library(SpatialExperiment)
#'
#' # load example data
#' spe <- STexampleData::Visium_humanDLPFC()
#'
#' # change from gene id to gene names
#' rownames(spe) <- rowData(spe)$gene_name
#'
#' # drop out-of-tissue spots
#' spe <- spe[, spe$in_tissue == 1]
#' spe <- spe[, !is.na(spe$ground_truth)]
#'
#' # Identifying the mitochondrial transcripts in our SpatialExperiment.
#' is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]
#'
#' # Calculating QC metrics for each spot using scuttle
#' spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
#' colnames(colData(spe))
#'
#' # Identifying local outliers using SpotSweeper
#' spe <- localOutliers(spe,
#'                      metric = "sum",
#'                      direction = "lower",
#'                      log = TRUE
#' )
#'
localOutliers <- function(
        spe, metric = "detected",
        direction = "lower", n_neighbors = 36, samples = "sample_id",
        log = TRUE, cutoff = 3) {

  # ===== Validity checks =====
  # Check if 'spe' is a valid object with required components
  if (!("SpatialExperiment" %in% class(spe))) {
    stop("Input data must be a SpatialExperiment object.")
  }

  # Validate 'direction'
  if (!direction %in% c("lower", "higher", "both")) {
    stop("'direction' must be one of 'lower', 'higher', or 'both'.")
  }

  # Check 'n_neighbors' is a positive integer
  if (!is.numeric(n_neighbors) ||
      n_neighbors <= 0 ||
      n_neighbors != round(n_neighbors)) {
    stop("'n_neighbors' must be a positive integer.")
  }

  # Check 'cutoff' is a numeric value
  if (!is.numeric(cutoff)) {
    stop("'cutoff' must be a numeric value.")
  }

  # ===== Start function =====
  # log transform specified metric
  if (log) {
      metric_log <- paste0(metric, "_log")
      colData(spe)[metric_log] <- log1p(colData(spe)[[metric]])
      metric_to_use <- metric_log
  } else {
      metric_to_use <- metric
  }

  # Get a list of unique sample IDs
  unique_sample_ids <- unique(colData(spe)[[samples]])

  # Initialize list to store each columnData
  columnData_list <- sapply(unique_sample_ids, FUN = function(x) NULL)

  # Loop through each unique sample ID
  for (sample in unique_sample_ids) {
    # Subset the data for the current sample
    spe_subset <- spe[, colData(spe)[[samples]] ==
                        sample]

    # Create a list of spatial coordinates and qc features
    columnData <- colData(spe_subset)

    # Find nearest neighbors
    dnn <- BiocNeighbors::findKNN(spatialCoords(spe_subset),
                                  k = n_neighbors, warn.ties = FALSE
    )$index

    # get neighborhood metrics
    neighborhoods <- lapply(seq_len(nrow(dnn)), function(i) {
      indices <- dnn[i, ]
      indices <- indices[indices != 0]
      indices <- c(i, indices)

      columnData[indices, metric_to_use]
    })

    # Compute modified-z and return the middle spot
    mod_z_matrix <- vapply(neighborhoods, function(x) {
      spatialEco::outliers(x)[1]
    }, numeric(1))

    # find outliers based on cutoff, store in colData
    metric_outliers <- paste0(metric, "_outliers")
    columnData[metric_outliers] <- switch(direction,
                                     higher = sapply(
                                       mod_z_matrix,
                                       function(x) x > cutoff
                                     ),
                                     lower = sapply(
                                       mod_z_matrix,
                                       function(x) x < -cutoff
                                     ),
                                     both = sapply(
                                       mod_z_matrix,
                                       function(x) {
                                         x > cutoff | x <
                                           -cutoff
                                       }
                                     )
    )

    # add z-scores to colData
    metric_z <- paste0(metric, "_z")
    columnData[metric_z] <- mod_z_matrix[]

    # Store the modified columnData dataframe in the list
    columnData_list[[sample]] <- columnData
  }

  # rbind the list of dataframes
  columnData_aggregated <- do.call(rbind, columnData_list)

  # replace column data
  colData(spe) <- columnData_aggregated

  return(spe)
}
