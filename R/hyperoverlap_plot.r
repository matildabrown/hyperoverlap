#' Overlap plotting for low-dimensional spaces
#'
#' Plot the optimal separating hyperplane found by hyperoverlap_detect() in 3D .
#'
#' @param x An \code{\link{hyperoverlap-class}} object.
#' @seealso \code{\link{hyperoverlap_detect}} , \code{\link{hyperoverlap_lda}}
#'
#' @examples
#' \dontrun{
#' hyperoverlap_plot(x)
#' }
#' @export



hyperoverlap_plot = function(x){
  if(methods::is(x,"Hyperoverlap")==FALSE){
    stop("hyperoverlap_result must be a Hyperoverlap object")
  }

  if (length(x@dimensions)>=4){
    stop("Plotting of decision boundary not supported in >3 dimensions. Use hyperoverlap_lda()")
  }


  grid.length = 50
  newdat.list <-  lapply(x@occurrences[,-1], function(z) seq(min(z), max(z), len=grid.length))
  newdat      <-  expand.grid(newdat.list)

  #find the decision boundary coordinates
  newdat.pred <-  stats::predict(x@model, newdata=newdat, decision.values=T)
  newdat.dv   <-  attr(newdat.pred, 'decision.values')
  newdat.dv   <-  array(newdat.dv, dim=rep(50, n))

  dat = x@occurrences
  if (x@number.of.points.misclassified == length(which(dat[,1]==x@entity1))) stop("No decision boundary found.")
  if (x@number.of.points.misclassified == length(which(dat[,1]==x@entity2))) stop("No decision boundary found.")

  if (length(x@dimensions)==2){
    graphics::plot(dat[,-1], col=c("red","blue")[as.factor(dat[,1])], main = paste0("Hyperoverlap: ", x@entity1, " and ", x@entity2))
    graphics::contour(newdat.dv,level = 0, add=TRUE, x =newdat.list[[1]], y=newdat.list[[2]])
  }

  if (length(x@dimensions)==3){
  rgl::plot3d(dat[,-1], col=c("red","blue")[as.factor(dat[,1])], alpha=0.3, size = 7, main = paste0("Hyperoverlap: ", x@entity1, " and ", x@entity2))
  misc3d::contour3d(newdat.dv,level = 0, add=TRUE, alpha=0.7, x =newdat.list[[1]], y=newdat.list[[2]], z=newdat.list[[3]])
  }


}
