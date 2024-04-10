#' localVariance Function
#'
#' This function does calculates the local variance based on kNN.
#'
#' @param spe SpatialExperiment object with the following columns in colData:
#'        sample_id, sum_umi, sum_gene
#' @param n_neighbors Number of nearest neighbors to use for variance
#'        calculation
#' @param features Features to use for variance calculation
#' @param samples Column in colData to use for sample ID
#' @param log Whether to log1p transform the features
#' @param name Name of the new column to add to colData
#'
#' @return SpatialExperiment object with feature variance added to colData
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocNeighbors findKNN
#' @importFrom MASS rlm
#'
#' @export localVariance
#'
#' @examples
#'
#' # for more details see extended example in vignettes
#' library(SpotSweeper)
#' library(SpatialExperiment)
#' library(escheR)
#'
#' # load example data
#' spe <- STexampleData::Visium_humanDLPFC()
#'
#' # change from gene id to gene names
#' rownames(spe) <- rowData(spe)$gene_name
#'
#' # show column data before SpotSweepR
#' colnames(colData(spe))
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
#' spe <- localVariance(spe,
#'     features = "subsets_Mito_percent",
#'     n_neighbors = 36,
#'     name = "local_mito_variance_k36"
#'     )
#'
#' plotQC(spe, metric="local_mito_variance_k36")
#'
localVariance <- function(spe, n_neighbors = 36,
                          features = c("expr_chrM_ratio"),
                          samples = "sample_id", log = FALSE, name = NULL) {

    # ===== Validity checks =====
    # Check if 'spe' is a valid object with required components
    if (!("SpatialExperiment" %in% class(spe))) {
      stop("Input data must be a SpatialExperiment object.")
    }

    # Validate 'features' is a character vector
    if (!is.character(features)) {
      stop("'features' must be a character vector.")
    }

    # validate 'features' are valid colData
    if (!all(features %in% colnames(colData(spe)))) {
      stop("All features must be present in colData.")
    }

    # validate samples is valid colData
    if (!samples %in% colnames(colData(spe))) {
      stop("Samples column must be present in colData.")
    }

    # Check 'n_neighbors' is a positive integer
    if (!is.numeric(n_neighbors) ||
        n_neighbors <= 0 ||
        n_neighbors != round(n_neighbors)) {
      stop("'n_neighbors' must be a positive integer.")
    }

    # ===== start function =====
    # log1p transform specified features
    if (log) {
      features_log <- lapply(features, function(feature) {
        feature_log <- paste0(feature, "_log")
        colData(spe)[[feature_log]] <- log1p(colData(spe)[[feature]])
        return(feature_log)  # Return the new feature name
      })

      features_to_use <- c(features, unlist(features_log))
    } else {
      features_to_use <- features
    }

    # Get a list of unique sample IDs
    unique_sample_ids <- unique(colData(spe)[[samples]])

    # Initialize list to store each columnData dataframe
    columnData_list <- sapply(unique_sample_ids, FUN = function(x) NULL)

    # Loop through each unique sample ID
    for (sample_id in seq_along(unique_sample_ids)) {
        # Subset the data for the current sample
        sample <- unique_sample_ids[sample_id]
        spe_subset <- subset(spe, , sample_id == sample)

        # Create a list of spatial coordinates and qc features
        columnData <- colData(spe_subset)
        columnData$coords <- spatialCoords(spe_subset)

        # Find nearest neighbors
        dnn <- BiocNeighbors::findKNN(spatialCoords(spe_subset),
            k = n_neighbors,
            warn.ties = FALSE
        )$index

        #  === Get local variance ===
        # Initialize a matrix to store variance for each feature
        var_matrix <- matrix(NA, nrow(columnData), length(features_to_use))
        colnames(var_matrix) <- features_to_use

        mean_matrix <- matrix(NA, nrow(columnData), length(features_to_use))
        colnames(mean_matrix) <- features_to_use

        # Loop through each row in the nearest neighbor index matrix
        for (i in seq_len(nrow(dnn))) {
            dnn.idx <- dnn[i, ]
            for (j in seq_along(features_to_use)) {
                neighborhood <- columnData[c(i, dnn.idx[dnn.idx != 0]), ][[features_to_use[j]]]

                var_matrix[i, j] <- var(neighborhood, na.rm = TRUE)[1]
                mean_matrix[i, j] <- mean(neighborhood, na.rm = TRUE)[1]
            }
        }

        # Handle non-finite values
        var_matrix[!is.finite(var_matrix)] <- 0
        mean_matrix[!is.finite(mean_matrix)] <- 0

        # == Regress out mean-variance bias ==
        for (feature_idx in seq_along(features_to_use)) {
            # Prepare data.frame for current feature
            mito_var_df <- data.frame(
                mito_var = log2(var_matrix[, feature_idx]),
                mito_mean = log2(mean_matrix[, feature_idx])
            )

            # Perform robust linear regression (IRLS) of variance vs mean
            fit.irls <- MASS::rlm(mito_var ~ mito_mean, data = mito_var_df)

            # Get residuals and update the variance matrix
            resid.irls <- resid(fit.irls)

            # Replace original variance values with residuals
            var_matrix[, feature_idx] <- resid.irls
        }

        # add local variance to columnData dataframe
        if (!is.null(name)) {
            columnData[name] <- var_matrix[, j]
        } else {
            feature_var <- paste0(features[j], "_var")
            columnData[feature_var] <- var_matrix[, j]
        }

        # Store the modified columnData dataframe in the list
        columnData_list[[sample_id]] <- columnData
    }

    # rbind the list of dataframes
    columnData_aggregated <- do.call(rbind, columnData_list)

    # replace SPE column data with aggregated data
    colData(spe) <- columnData_aggregated

    return(spe)
}
