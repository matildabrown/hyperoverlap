#' Pairwise overlap detection in n-dimensional space of multiple entities using support vector machines (SVMs)
#'
#' This function is a wrapper for \code{\link{hyperoverlap_detect}} for pairwise overlap detection between multiple entities.
#'
#' @usage hyperoverlap_set(x, y, kernel = "polynomial",kernel.degree = 3, cost = 1000,
#' stoppage.threshold = 0.2, write.to.file = FALSE,
#' path = NULL,
#' sample.dimensionality.omit = "FALSE")
#'
#' @param x A matrix or data.frame containing the variables of interest for both entities.
#' @param y A vector of labels.
#' @param kernel Character. Either "linear" or "polynomial" (default = "polynomial").
#' @param kernel.degree Parameter needed for \code{kernel = polynomial} (default = 3).
#' @param cost Specifies the SVM margin 'hardness'. Default value is 1000, but can be increased for improved accuracy (although this increases runtimes and memory usage).
#' @param stoppage.threshold Numeric. If the number of points misclassified using a linear hyperplane exceeds this proportion of the number of observations, non-linear separation is not attempted. Must be between 0 and 1 (default = 0.2).
#' @param write.to.file Logical. If TRUE, each \code{\link{hyperoverlap-class}} object is saved as a .rds file.
#' @param path Character. Path to write .rds files to. Ignored if \code{write.to.file=FALSE}
#' @param sample.dimensionality.omit Logical. If TRUE, omits any entity pairs with a combined sample size less than n+1, where n is the number of dimensions (see details).
#'
#' @return A long-form matrix with the following columns:
#' entity1,
#' entity2,
#' shape,
#' polynomial.order (if \code{kernel="polynomial"}),
#' result,
#' number.of.points.misclassified.
#'
#'If specified, individual \code{Hyperoverlap-class} objects are written to file.
#'
#'@details In n dimensions, any set of points up to n+1 points can be separated using a linear hyperplane. This may produce an artefactual non-overlap result.
#'The \code{sample.dimensionality.omit} parameter gives two options for dealing with these pairs when they form part of a larger analysis.
#'If \code{sample.dimensionality.omit = "TRUE"}, this pair is removed from the analysis (result = NA).
#'If \code{sample.dimensionality.omit = "FALSE"}, the pair is included, but a warning is printed.
#'
#' @examples
#'
#' data(iris)
#' hyperoverlap.iris.set = hyperoverlap_set(iris[1:3],iris$Species, kernel="linear")
#'
#' @export


hyperoverlap_set <-  function(x,y,kernel="polynomial",kernel.degree=3, cost=1000,stoppage.threshold=0.2, write.to.file=FALSE, path=NULL, sample.dimensionality.omit="FALSE"){

  #if data are incomplete, break
  if (length(which(is.na(x)==TRUE))>0){
    stop("Missing data. NAs detected.")
  }
  if (length(which(x==""))>0){
    stop("Missing data. Empty elements detected.")
  }

  if (length(which(is.na(y)==TRUE))>0){
    stop("Missing labels. NAs detected.")
  }

  if (length(which(y==""))>0){
    stop("Missing labels. Empty elements detected.")
  }


  if (length(y)!=nrow(x)){
    stop("Number of labels does not match number of points. ")
  }

  dat=cbind(y,x)

  entities <- unique(y)


  if (kernel=="linear"){
    results <-  list(NULL,NULL,NULL,NULL,NULL)
    names(results) <-  c("entity1", "entity2", "shape", "result", "number.of.points.misclassified")
  } else {
    results = list(NULL,NULL,NULL,NULL,NULL,NULL)
    names(results) <-  c("entity1", "entity2", "shape", "polynomial.order", "result", "number.of.points.misclassified")
  }

  if(write.to.file==TRUE){
    if(!dir.exists(path)) dir.create(path)
  }

  for (i in 1:(length(entities))){
    dat.i <-  dat[which(dat[,1]==entities[i]),]


      for (j in 1:length(entities)){

        if(j==i){
          results[[1]] <- c(results[[1]], levels(entities)[i])
          results[[2]] <- c(results[[2]], levels(entities)[i])
          results[[3]] <- c(results[[3]], " ")

          if (kernel=="linear") {
            results[[4]] <- c(results[[4]], NA)
            results[[5]] <- c(results[[5]], NA)
          } else {
            results[[4]] <- c(results[[4]], NA)
            results[[5]] <- c(results[[5]], NA)
            results[[6]] <- c(results[[6]], NA)
          }
        }

        if (j>i){

      dat.j <-  dat[which(dat[,1]==entities[j]),]
      dat.ij  <-  rbind(dat.i,dat.j)
      if (nrow(dat.ij)<ncol(dat.ij)) {
        if(sample.dimensionality.omit == "FALSE") {
        rlang::warn(paste0("Dimensionality and low number of points mean that a boundary can always be found between ",entities[i], " and ",entities[j],"."))
        }
        if(sample.dimensionality.omit == "TRUE"){break(paste0("Dimensionality and low number of points mean that a boundary can always be found between ",entities[i], " and ",entities[j],". This pair has been excluded."))}
      }
      hyperoverlap.ij  <-  hyperoverlap_detect(dat.ij[,-1],dat.ij[,1],kernel=kernel,kernel.degree=kernel.degree, cost=cost,stoppage.threshold=stoppage.threshold, set=TRUE)

      results[[1]] <- c(results[[1]], hyperoverlap.ij@entity1)
      results[[2]] <- c(results[[2]], hyperoverlap.ij@entity2)
      results[[3]] <- c(results[[3]], hyperoverlap.ij@shape)

      if (kernel=="linear") {
        results[[4]] <- c(results[[4]], hyperoverlap.ij@result)
        results[[5]] <- c(results[[5]], hyperoverlap.ij@number.of.points.misclassified)
      } else {
        results[[4]] <- c(results[[4]], hyperoverlap.ij@polynomial.order)
        results[[5]] <- c(results[[5]], hyperoverlap.ij@result)
        results[[6]] <- c(results[[6]], hyperoverlap.ij@number.of.points.misclassified)
      }

      if (write.to.file==TRUE) saveRDS(hyperoverlap.ij, file=paste0(path,hyperoverlap.ij@entity1,"_",hyperoverlap.ij@entity2,".rds"))
    }
  }
}
  results <-  data.frame(results)
  #results2 <- cbind(results$entity2, results$entity1, results[,3:ncol(results)])
  #colnames(results2) <- colnames(results)
  #results <- rbind(results,results2)
  #results <- unique(results)
  results$entity1 <- factor(results$entity1,levels=entities)
  results$entity2 <-  factor(results$entity2,levels=entities)
  results <- results[order(results$entity1,results$entity2),]

  return(results)

}
