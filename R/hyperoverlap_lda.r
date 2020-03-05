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


hyperoverlap_lda <- function (x, return.plot = TRUE, visualise3d = FALSE) {
  occ <- x@occurrences
  n=length(x@dimensions)


  #set up reference points for biplot
  ref = matrix(nrow=n, ncol=n, data=0)
  diag(ref)=1


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

  orth <- matlib::GramSchmidt(ord, normalize = TRUE) #force orthogonality

  while (length(which(orth==0))>0){ #not allow duplicate axes, zero axes
    for (i in 2:n){
      ord[,i] <- stats::runif(n, -1, 1)
    }
    orth <- matlib::GramSchmidt(ord, normalize = TRUE)
  }

  #set up transformation matrix
  tran <- occ
  cols <- c(colnames(occ)[1], "LDA1")
  for (i in 2:n) {
    cols <- c(cols, paste0("O", i - 1))
  }
  tran[,2:ncol(tran)] <- NA
  colnames(tran) <- cols

  #transform data and ref matrix
  tran[,2:ncol(tran)] <- as.matrix(occ[,-1])%*%orth
  ref1 = ref%*%orth



    if (visualise3d == FALSE) {

  ########2d VISUALISATION######
    pca <- stats::princomp(tran[, 3:(n + 1)])

    #transform ref matrix by pca
    ref2=ref1[,-1]%*%pca$loadings
    ref2 = cbind(ref1[,1],ref2)


    #set up segments for variable importance plot
    seg = matrix(nrow=n,ncol=4)
    colnames(seg)=c("x1","y1","x2","y2")
      
       tran2 <- tran[1:3]
    colnames(tran2)[3] <- "residualPCA"

    #add second axis to lda scores
    tran2[, 3] <- pca$scores[,1]
      
    midx = min(tran2[,2])+(max(tran2[,2]-min(tran2[,2]))/2)
    midy = min(tran2[,3])+(max(tran2[,3]-min(tran2[,3]))/2)
    seg[,1] = rep(midx,n)
    seg[,2] = rep(midy,n)
    for (i in 1:n){
      seg[i,3]=midx+ref2[i,1]
      seg[i,4]=midy+ref2[i,2]
    }

   

    if (return.plot == TRUE) {
    #plotting
      par(mfrow=c(1,2)) #plot data, variable transformation side-by-side
      graphics::plot(tran2[, 2:3], col = c("red", "blue")[as.factor(tran2$Entity)], main = paste0("Hyperoverlap: ", x@entity1, " and ", x@entity2))
      graphics::plot(tran2[, 2:3], col = "grey")
      segments(seg[,1],seg[,2],seg[,3],seg[,4])
      text(seg[,3],seg[,4],colnames(occ)[2:n])
      plotrix::draw.circle(seg[2,1],seg[1,2],1)
      segments(rep(0,3),rep(0,3),ref2[,1],ref2[,2])
      text(ref2[,1],ref2[,2],colnames(data[,1:4]))

    }
  }
    if (visualise3d == TRUE) {
      #########3D visualisation########
      pca <- stats::princomp(tran[, 3:(n + 1)])
      ref2=ref1[,-1]%*%pca$loadings
      ref2 = cbind(ref1[,1],ref2)

      tran2 <- tran[1:4]

      colnames(tran2)[3:4] <- c("residualPCA1", "residualPCA2")
      tran2[, 3:4] <- pca$scores[,1:2]

      #set up segments for variable importance plot
      seg = matrix(nrow=n,ncol=6)
      colnames(seg)=c("x1","y1","z1","x2","y2","z2")
      midx = min(tran2[,2])+(max(tran2[,2]-min(tran2[,2]))/2)
      midy = min(tran2[,3])+(max(tran2[,3]-min(tran2[,3]))/2)
      midz = min(tran2[,4])+(max(tran2[,4]-min(tran2[,4]))/2)
      seg[,1] = rep(midx,n)
      seg[,2] = rep(midy,n)
      seg[,3] <- rep(midz,n)

      for (i in 1:n){
        seg[i,4]=midx+ref2[i,1]
        seg[i,5]=midy+ref2[i,2]
        seg[i,6]=midz+ref2[i,3]
      }

      if (return.plot == TRUE){
        #####plotting
        rgl::mfrow3d(1, 1, sharedMouse = TRUE)
        rgl::plot3d(tran2[, 2:4], col = c("red", "blue")[as.factor(tran2$Entity)], alpha=0.3, size = 7, main = paste0("Hyperoverlap: ", x@entity1, " and ", x@entity2))
        aspect3d(1)
        rgl::plot3d(tran2[, 2:4], col = "grey", alpha=0.2,asp=1, size = 7)
        aspect3d(1)
        segments3d(x=as.vector(t(seg[,c(1,4)])),
                   y=as.vector(t(seg[,c(2,5)])),
                   z=as.vector(t(seg[,c(3,6)])))
        text3d(x=seg[,4],y=seg[,5],z=seg[,6],colnames(data[,1:4]))
    }
     }
  return(tran2)
}
