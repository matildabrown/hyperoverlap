#' Overlap detection in n-dimensional space using support vector machines (SVMs)
#'
#' Given a matrix containing the ecological data (x) and labels (y) for two entities, a support vector machine is trained and the predicted label of each point is evaluated. If every point has been classified correctly, the entities can be separated and they do not overlap.
#'
#' @usage hyperoverlap_detect(x, y, kernel = "polynomial", kernel.degree = 3, cost = 50,
#'        stoppage.threshold = 0.2)
#'
#' @param x A matrix or data.frame containing the variables of interest for both entities.
#' @param y A vector of labels.
#' @param kernel Character. Either "linear" or "polynomial" (default = "polynomial").
#' @param kernel.degree Parameter needed for \code{kernel = polynomial} (default = 3).
#' @param cost Specifies the SVM margin 'hardness'. Default value is 50, but can be increased for improved accuracy (although this increases runtimes and memory usage).
#' @param stoppage.threshold Numeric. If the number of points misclassified using a linear hyperplane exceeds this proportion of the number of observations, non-linear separation is not attempted. Must be between 0 and 1 (default = 0.2).
#'
#' @return A \code{\link{hyperoverlap-class}} object
#' @export
#' @examples
#' \dontrun{
#' data = iris[which(iris$Species!=("versicolor")),]
#' x = hyperoverlap_detect(data[,1:3],data$Species, kernel="linear")
#' }
#' @details Input data should be preprocessed so that all variables are comparable (e.g. same order of magnitude). Polynomial kernels allow curvilinear decision boundaries to be found between entities (see \url{https://www.cs.cmu.edu/~ggordon/SVMs/new-svms-and-kernels.pdf}). Smaller values of \code{kernel.degree} permit less complex decision boundaries; biological significance is likely to be lost at values > 5.

hyperoverlap_detect <- function(x, y, kernel="polynomial",kernel.degree=3, cost=50,stoppage.threshold=0.2){

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

  if (length(unique(y))>2){
    stop("More than two entities detected. For multiple pairs, see hyperoverlap_set(). ")
  }

  if (nrow(x)<=(ncol(x)+1)){
    stop("Total number of points too small and will always produce a boundary. Decrease dimensionality or exclude this pair from analysis ")
  }

  if (stoppage.threshold<0){
    stop("Invalid stoppage threshold. Value must be between 0 and 1.")
  }

  if (stoppage.threshold>1){
    stop("Invalid stoppage threshold. Value must be between 0 and 1.")
  }

  n <- ncol(x)
  m <-  nrow(x)
  y <- factor(y)
  entity1 <-  as.character(levels(factor(y))[1])
  entity2 <-  as.character(levels(factor(y))[2])
  dimensions <-  colnames(x)
  occurrences <-  cbind(y,x)

  newdat.list <-  lapply(x, function(x) seq(min(x), max(x), len=50))
  newdat      <-  expand.grid(newdat.list)

  #fit the SVM, evaluate how well it separated the data
  fit <-  e1071::svm(formula = y~.,data=occurrences, scale=FALSE, kernel ="linear", cost=cost)
  pred <-  stats::predict(fit,x)
  accuracy <- table(pred,y)
  misclass <-  accuracy[2]+accuracy[3]


  # if the svm found a hyperplane that perfectly separates the entity data, return "linear non-overlap" and save the coordinates of the decision boundary
  if(misclass==0){

    #find the decision boundary coordinates
    newdat.pred <-  stats::predict(fit, newdata=newdat, decision.values=T)
    newdat.dv   <-  attr(newdat.pred, 'decision.values')
    newdat.dv   <-  array(newdat.dv, dim=rep(50, n))

    return(methods::new("Hyperoverlap",
               entity1=entity1,
               entity2=entity2,
               dimensions=dimensions,
               occurrences = as.matrix(occurrences),
               shape="linear",
               polynomial.order = 1,
               result = "non-overlap",
               accuracy = accuracy,
               number.of.points.misclassified = 0,
               decisionboundary=list(newdat.dv,newdat.list)
    ))

  } else {
    if(kernel=="linear"){
      #find the decision boundary coordinates

      newdat.pred <-  stats::predict(fit, newdata=newdat, decision.values=T)
      newdat.dv   <-  attr(newdat.pred, 'decision.values')
      newdat.dv   <-  array(newdat.dv, dim=rep(50, n))

      return(methods::new("Hyperoverlap",
                 entity1=entity1,
                 entity2=entity2,
                 dimensions=dimensions,
                 occurrences = as.matrix(occurrences),
                 shape="linear",
                 polynomial.order = 1,
                 result = "overlap",
                 accuracy = accuracy,
                 number.of.points.misclassified = misclass,
                 decisionboundary=list(newdat.dv,newdat.list)
      ))
    }
    if (kernel =="polynomial"){
      if (misclass > stoppage.threshold*nrow(x)) {
        stop("Misclassified points above stoppage threshold; curvilinear hyperoverlap not attempted.")
      }
      deg <- 2
      bestdeg <- deg
      bestfit <- fit
      while (deg <= kernel.degree){

        fit <-  e1071::svm(formula = y~.,data=occurrences, scale=FALSE, kernel ="polynomial",degree=deg, cost=cost)
        pred <-  stats::predict(fit,x)
        accuracy2 <- table(pred,y)
        misclass2 <-  accuracy[2]+accuracy[3]


        if(misclass2<misclass){
          misclass <- misclass2
          accuracy <- accuracy2
          bestdeg <- deg
          bestfit <- fit
        }

        if (misclass == 0) {

          #find the decision boundary coordinates
          newdat.pred <-  stats::predict(fit, newdata=newdat, decision.values=T)
          newdat.dv   <-  attr(newdat.pred, 'decision.values')
          newdat.dv   <-  array(newdat.dv, dim=rep(50, n))

          return(methods::new("Hyperoverlap",
                     entity1=entity1,
                     entity2=entity2,
                     dimensions=dimensions,
                     occurrences = as.matrix(occurrences),
                     shape="curvilinear",
                     polynomial.order = deg,
                     result = "non-overlap",
                     accuracy = accuracy,
                     number.of.points.misclassified = "0",
                     decisionboundary=list(newdat.dv,newdat.list)
          ))
          break()
        } else{
          deg <- deg+1
          if(deg>kernel.degree){

            #find the decision boundary coordinates
            newdat.pred <-  stats::predict(bestfit, newdata=newdat, decision.values=T)
            newdat.dv   <-  attr(newdat.pred, 'decision.values')
            newdat.dv   <-  array(newdat.dv, dim=rep(50, n))

            return(methods::new("Hyperoverlap",
                       entity1=entity1,
                       entity2=entity2,
                       dimensions=dimensions,
                       occurrences = as.matrix(occurrences),
                       shape="curvilinear",
                       polynomial.order = bestdeg,
                       result = "overlap",
                       accuracy=accuracy,
                       number.of.points.misclassified = misclass,
                       decisionboundary=list(newdat.dv,newdat.list)
            ))

          }
        }
      }
    }
  }

}
