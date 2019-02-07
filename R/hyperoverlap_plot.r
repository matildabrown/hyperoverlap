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
  dat = x@occurrences
  rgl::plot3d(dat[,-1], col=c("red","blue")[as.factor(dat[,1])], alpha=0.3, size = 10, main = paste0("Hyperoverlap: ", x@entity1, " and ", x@entity2))
  if (x@number.of.points.misclassified == length(which(dat[,1]==x@entity1))) stop("No decision boundary found.")
  if (x@number.of.points.misclassified == length(which(dat[,1]==x@entity2))) stop("No decision boundary found.")

  misc3d::contour3d(x@decisionboundary[[1]],level = 0, add=TRUE, x =x@decisionboundary[[2]][[1]], y=x@decisionboundary[[2]][[2]], z=x@decisionboundary[[2]][[3]])


}
