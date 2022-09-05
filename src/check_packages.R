check_packages <- function(cran_packages = c(""),
                           bioc_packages = c("")) {
  # Description:
  #   If package is not installed, the function installs the package
  # Arguments:
  #   cran_packages: packages available at CRAN.
  #   bioc_packages: packages available at bioconductor
  if (!"" %in% cran_packages) {
    lapply(
      cran_packages,
      FUN = function(x) {
        if (!x %in% rownames(installed.packages())) {
          stop(paste0("ERROR: package ", x," is not found."))
        }
      }
    )
  }
  if (!"" %in% bioc_packages) {
    lapply(
      bioc_packages,
      FUN = function(x) {
        if (!x %in% rownames(installed.packages())) {
          stop(paste0("ERROR: package ", x," is not found."))
        }
      }
    )
  }
  paste0("Passed dependency check.")
}