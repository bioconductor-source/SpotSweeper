#' Identify and annotate artifacts in spatial transcriptomics data
#'
#' This function identifies and annotates potential artifacts in spatial
#' transcriptomics data. Artifacts are detected based on local mito variance,
#' and the results are added to the original
#' SpatialExperiment (sce) object.
#'
#' @param spe A SingleCellExperiment object.
#' @param mito_percent The column name representing the mitochondrial percent.
#' Default is 'expr_chrM_ratio'.
#' @param mito_sum The column name representing sum mitochondrial expression.
#' Default is 'expr_chrM'.
#' @param samples The column name representing sample IDs. Default is
#' 'sample_id'.
#' @param n_rings The number of rings for local mito variance calculation.
#' Default is 5.
#' @param log Logical, indicating whether to log1p transform mito_percent.
#' Default is TRUE.
#' @param name Prefix for the local variance column names. Default is
#' 'artifact'.
#' @param var_output Logical, indicating whether to include local variances in
#' the output. Default is TRUE.
#'
#' @return Returns the modified SingleCellExperiment object with artifact
#' annotations.
#'
#' @seealso
#' \code{\link{localVariance}}
#'
#' @import SingleCellExperiment
#' @import SpatialExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom stats prcomp kmeans resid var
#'
#'
#' @examples
#' library(SpotSweeper)
#' library(SpatialExperiment)
#' library(escheR)
#'
#' data(DLPFC_artifact)
#' spe <- DLPFC_artifact
#'
#' # find artifacts
#' spe <- findArtifacts(spe,
#'     mito_percent = "expr_chrM_ratio",
#'     mito_sum = "expr_chrM",
#'     n_rings = 5,
#'     name = "artifact"
#' )
#'
#' plotQC(spe, metric="expr_chrM_ratio", outliers="artifact")
#'
#' @export
findArtifacts <- function(
        spe, mito_percent = "expr_chrM_ratio",
        mito_sum = "expr_chrM", samples = "sample_id", n_rings = 5,
         log = TRUE, name = "artifact", var_output = TRUE) {

    # ===== Validity checks =====
    if (!("SpatialExperiment" %in% class(spe))) {
      stop("Input data must be a SpatialExperiment object.")
    }

    if (!all(mito_percent %in% colnames(colData(spe)))) {
      stop("mito_percent must be present in colData.")
    }

    if (!mito_sum %in% colnames(colData(spe))) {
      stop("mito_sum must be present in colData.")
    }

    if (!samples %in% colnames(colData(spe))) {
      stop("Samples column must be present in colData.")
    }

    if (!is.numeric(n_rings) ||
        n_rings <= 0 ||
        n_rings != round(n_rings)) {
      stop("'n_rings' must be a positive integer.")
    }

    # get unique sample IDs
    unique_sample_ids <- unique(colData(spe)[[samples]])

    # Initialize a list to store spe for each sample
    columnData_list <- sapply(unique_sample_ids, FUN = function(x) NULL)

    for (sample in unique_sample_ids) {
        # subset by sample
        spe.temp <- spe[, colData(spe)[[samples]] == sample]

        # ======= Calculate local mito variance ========
        # Use vapply to iterate over rings_seq
        var_matrix <- vapply(seq_len(n_rings), function(i) {
          # Calculate n_neighbors for the current ring
          n_neighbors <- 3 * i * (i + 1)
          tmp.name <- paste0("k", n_neighbors)

          # Apply local variance calculation
          spe.temp <<- localVariance(spe.temp,
                                    metric = mito_percent,
                                    n_neighbors = n_neighbors,
                                    name = tmp.name,
                                    log=log)

          # Extract and return the column corresponding to tmp.name from colData
          colData(spe.temp)[[tmp.name]]
        }, numeric(length(spe.temp[[1]])))


        # ========== PCA and clustering ==========
        var_matrix <- cbind(
            var_matrix,
            colData(spe.temp)[[mito_percent]],
            colData(spe.temp)[[mito_sum]]
        )

        var_df <- data.frame(var_matrix)

        # Run PCA and add to reduced dims
        pc <- prcomp(var_df, center = TRUE, scale. = TRUE)
        rownames(pc$x) <- colnames(spe.temp) # assign rownames to avoid error
        reducedDim(spe.temp, "PCA_artifacts") <- pc$x

        # Cluster using kmeans and add to temp sce
        clus <- kmeans(pc$x, centers = 2, nstart = 25)
        spe.temp$Kmeans <- clus$cluster

        # =========== Artifact annotation ===========

        # calculate average local variance of the two clusters
        clus1_mean <- mean(colData(spe.temp)[[paste0("k",18
        )]][spe.temp$Kmeans == 1])
        clus2_mean <- mean(colData(spe.temp)[[paste0("k",18
        )]][spe.temp$Kmeans == 2])

        artifact_clus <- which.min(c(clus1_mean, clus2_mean))

        # create a new $artifact column; if clus1 mean < clus2 - rename Kmeans
        # 1 to 'TRUE'
        spe.temp$artifact <- FALSE
        spe.temp$artifact[spe.temp$Kmeans == artifact_clus] <- TRUE
        spe.temp$Kmeans <- NULL

        # add to sample list
        columnData_list[[sample]] <- colData(spe.temp)
    }

    # rbind the list of dataframes
    columnData_aggregated <- do.call(rbind, columnData_list)

    # replace SPE column data with aggregated data
    colData(spe) <- columnData_aggregated
    reducedDims(spe) <- reducedDims(spe.temp)

    # return sce
    return(spe)
}
