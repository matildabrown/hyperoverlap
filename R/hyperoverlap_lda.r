#' Hyperoverlap visualisation using linear discriminant analysis (LDA)
#'
#' @param x An \code{\link{hyperoverlap-class}} object.
#' @param return.plot Logical. If TRUE, data are plotted using \code{plot()}.
#' @param visualise3d Logical. If FALSE, data are projected onto two axes (LDA1, residualPCA1). If TRUE, data are projected onto three axes (LDA1, residualPCA1, residualPCA2)
#'
#' @usage hyperoverlap_lda(x, return.plot=TRUE, visualise3d=FALSE)
#'
#' @return Returns a dataframe with columns "Entity", "LDA1", "residualPCA1",
#' "residualPCA2" (if \code{visualise3d = TRUE})
#' @seealso \code{\link{hyperoverlap_detect}}
#'
#'
#' @examples
#' \dontrun{
#' #using iris dataset reduced to two species
#' data = iris[which(iris$Species!=("versicolor")),]
#' x = hyperoverlap_detect(data[1:4], data$Species)
#' hyperoverlap_lda(x)
#'}
#' @details This function provides a way to visualise overlap (or non-overlap)
#' between classes of high dimensional data. For inspection, it is useful to
#' use the base graphics package (implemented by return.plot=TRUE). The
#' transformed coordinates of each point are also returned as a dataframe,
#' which can be plotted with user-defined parameters.
#'
#'
#' @export


hyperoverlap_lda <- function (x, return.plot = TRUE, visualise3d = FALSE)
{
  occ <- x@occurrences
  if (is.data.frame(occ)==FALSE) occ <- as.data.frame(occ)
  n <- ncol(occ) - 1
  colnames(occ)[1] = "Entity"
  occ$Entity=factor(occ$Entity)
  for (i in 2:ncol(occ)){
    occ[,i] <- as.numeric(occ[,i])
  }
  lda1 <- MASS::lda(formula = Entity ~ ., data = occ)
  ord <- matrix(nrow = n, ncol = n)
  ord[, 1] <- lda1$scaling
  for (i in 2:n){
    ord[,i] <- stats::runif(n, -1, 1)
  }

  orth <- matlib::GramSchmidt(ord, normalize = TRUE)

  while (length(which(orth==0))>0){
    for (i in 2:n){
      ord[,i] <- stats::runif(n, -1, 1)
    }
    orth <- matlib::GramSchmidt(ord, normalize = TRUE)
  }

  tran <- occ
  cols <- c(colnames(occ)[1], "LDA1")
  for (i in 2:n) {
    cols <- c(cols, paste0("O", i - 1))
  }
  tran[,2:ncol(tran)] <- NA
  colnames(tran) <- cols
  for (i in 2:(n + 1)) {
    m <- matrix(nrow = nrow(occ), ncol = n)
    for (j in 1:n) {
      m[, j] <- (occ[, i] * orth[j, (i - 1)])
    }
    tran[, (i)] <- rowSums(m)
  }
  if (visualise3d == FALSE) {
    pca <- stats::princomp(tran[, 3:(n + 1)])
    tran2 <- tran[1:3]
    colnames(tran2)[3] <- "residualPCA"
    m <- matrix(nrow = nrow(occ), ncol = (n - 1))
    for (j in 1:(n - 1)) {
      m[, j] <- (occ[, (j + 2)] * pca$loadings[j, 1])
    }
    tran2[, 3] <- rowSums(m)
    if (return.plot == TRUE) {
      graphics::plot(tran2[, 2:3], col = c("red", "blue")[as.factor(tran2$Entity)])
    }
    else {
      pca <- stats::princomp(tran[, 3:(n + 1)])
      tran2 <- tran[1:4]
      colnames(tran2)[3:4] <- c("residualPCA1", "residualPCA2")
      m <- matrix(nrow = nrow(occ), ncol = (n - 1))
      for (j in 1:(n - 1)) {
        m[, j] <- (occ[, (j + 2)] * pca$loadings[j,
                                                 1])
      }
      tran2[, 3] <- rowSums(m)
      for (j in 1:(n - 1)) {
        m[, j] <- (occ[, (j + 2)] * pca$loadings[j,
                                                 2])
      }
      tran2[, 4] <- rowSums(m)
      if (return.plot == TRUE)
        rgl::plot3d(tran2[, 2:4], col = c("red", "blue")[as.factor(tran2$Entity)])
    }
  }
  return(tran2)
}
