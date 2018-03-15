#' Raw zebrafish data
#' 
#' @name zebra_r
#' @title Raw zebrafish data
#' @description This dataset contains raw RNA-sequencing read counts from embryonic zebrafish olfactory sensory neurons. Treatment samples were treated with gallein and control samples were not treated with gallein.
#' @docType data
#' @format a \code{RData} instance, 1 row per gene
#' @details \itemize{
  #' \item ID gene name
  #' \item C.1 control replicate 1 raw read counts
  #' \item C.3 control replicate 2 raw read counts
  #' \item C.5 control replicate 3 raw read counts
  #' \item T.9 gallein-treated replicate 1 raw read counts
  #' \item T.11 gallein-treated replicate 2 raw read counts
  #' \item T.13 gallein-treated replicate 3 raw read counts
  #' }
  #'
#' @docType data
#' @keywords datasets
#' @name zebra_r
#' @usage data(zebra_r)
#' @format A data frame with 20,865 rows and 7 variables
#' @references
#' Risso D, Ngai J, Speed TP, Dudoit S (2014) Normalization of RNA-seq data using factor analysis of control genes or samples. Nature Biotechnol 32:896–902
NULL
