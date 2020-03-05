function (x, return.plot = TRUE, visualise3d = FALSE, showlegend = TRUE) 
{
  occ <- x@occurrences
  n <- length(x@dimensions)
  ref <- matrix(nrow = n, ncol = n, data = 0)
  diag(ref) <- 1
  if (is.data.frame(occ) == FALSE) 
    occ <- as.data.frame(occ)
  n <- ncol(occ) - 1
  colnames(occ)[1] <- "Entity"
  occ$Entity <- factor(occ$Entity)
  for (i in 2:ncol(occ)) {
    occ[, i] <- as.numeric(occ[, i])
  }
  lda1 <- MASS::lda(formula = Entity ~ ., data = occ)
  ord <- matrix(nrow = n, ncol = n)
  ord[, 1] <- lda1$scaling
  for (i in 2:n) {
    ord[, i] <- stats::runif(n, -1, 1)
  }
  orth <- matlib::GramSchmidt(ord, normalize = TRUE)
  while (length(which(orth == 0)) > 0) {
    for (i in 2:n) {
      ord[, i] <- stats::runif(n, -1, 1)
    }
    orth <- matlib::GramSchmidt(ord, normalize = TRUE)
  }
  tran <- occ
  cols <- c(colnames(occ)[1], "LDA1")
  for (i in 2:n) {
    cols <- c(cols, paste0("O", i - 1))
  }
  tran[, 2:ncol(tran)] <- NA
  colnames(tran) <- cols
  tran[, 2:ncol(tran)] <- as.matrix(occ[, -1]) %*% orth
  ref1 <- ref %*% orth
  if (visualise3d == FALSE) {
    pca <- stats::princomp(tran[, 3:(n + 1)])
    ref2 <- ref1[, -1] %*% pca$loadings
    ref2 <- cbind(ref1[, 1], ref2)
    seg <- matrix(nrow = n, ncol = 4)
    colnames(seg) <- c("x1", "y1", "x2", "y2")
    tran2 <- tran[1:3]
    colnames(tran2)[3] <- "residualPCA"
    tran2[, 3] <- pca$scores[, 1]
    midx <- min(tran2[, 2]) + (max(tran2[, 2] - min(tran2[, 
                                                         2]))/2)
    rangex <- max(tran2[, 2]) - min(tran2[, 2])
    midy <- min(tran2[, 3]) + (max(tran2[, 3] - min(tran2[, 
                                                         3]))/2)
    rangey <- max(tran2[, 3]) - min(tran2[, 3])
                 
    seg[, 1] <- rep(max(tran2[, 2]), n)
    seg[, 2] <- rep(midy+rangey, n)
    
    plot_factor <- min(c(rangex, rangey))*0.7
      
    for (i in 1:n) {
      seg[i, 3] <- max(tran2[, 2]) + plot_factor*ref2[i, 1]
      seg[i, 4] <- midy+rangey + plot_factor*ref2[i, 2]
    }
    if (return.plot == TRUE) {
      par(mfrow = c(1, 2))
      graphics::plot(tran2[, 2:3], col = c("red", "blue")[as.factor(tran2$Entity)]) 
          
      graphics::plot(tran2[, 2:3], col = "grey", 
                     xlim = c(min(tran2$LDA1), min(tran2$LDA1)+2*rangex),
                     ylim = c(min(tran2$residualPCA), min(tran2$residualPCA)+2*rangey))
      
      #plot the 'circle' 
      asp <- (par("pin")[1]/diff(par("usr")[1:2]))/(par("pin")[2]/diff(par("usr")[3:4]))
      plotrix::draw.ellipse(seg[2, 1], seg[1, 2], plot_factor/2,plot_factor/2, border="gray", lwd=1.5)
      
      segments(seg[, 1], seg[, 2], seg[, 3], seg[, 4], lwd=2, col="gray40")
      points(seg[, 3], seg[, 4], col='white', cex=2.7, pch=20)
      text(seg[, 3], seg[, 4], 1:n, cex=1.2)
      
    
      
      
      legend("bottomright", legend=colnames(occ)[2:(n+1)],
            pch = paste(seq(1:n)), cex=0.8)
      
      mtext(paste0("Hyperoverlap: ", x@entity1, " and ", 
                   x@entity2), side = 3, line = -3, outer = TRUE, cex=1.5, font=2)
      
    }
  }
  if (visualise3d == TRUE) {
    pca <- stats::princomp(tran[, 3:(n + 1)])
    ref2 <- ref1[, -1] %*% pca$loadings
    ref2 <- cbind(ref1[, 1], ref2)
    tran2 <- tran[1:4]
    colnames(tran2)[3:4] <- c("residualPCA1", "residualPCA2")
    tran2[, 3:4] <- pca$scores[, 1:2]
    seg <- matrix(nrow = n, ncol = 6)
    colnames(seg) <- c("x1", "y1", "z1", "x2", "y2", "z2")
    midx <- min(tran2[, 2]) + (max(tran2[, 2] - min(tran2[, 
                                                         2]))/2)
    midy <- min(tran2[, 3]) + (max(tran2[, 3] - min(tran2[, 
                                                         3]))/2)
    midz <- min(tran2[, 4]) + (max(tran2[, 4] - min(tran2[, 
                                                         4]))/2)
    seg[, 1] <- rep(midx, n)
    seg[, 2] <- rep(midy, n)
    seg[, 3] <- rep(midz, n)
    for (i in 1:n) {
      seg[i, 4] <- midx + ref2[i, 1]
      seg[i, 5] <- midy + ref2[i, 2]
      seg[i, 6] <- midz + ref2[i, 3]
    }
    if (return.plot == TRUE) {
      par3d(windowRect = c(20, 30, 1200, 600))
      rgl::mfrow3d(1, 2, sharedMouse = TRUE)
      rgl::plot3d(tran2[, 2:4], col = c("red", "blue")[as.factor(tran2$Entity)], 
                  alpha = 0.3, size = 7,)
      aspect3d(1)
      rgl::plot3d(tran2[, 2:4], col = "grey", alpha = 0.4, 
                  asp = 1, size = 7, xlab=NULL, ylab=NULL, zlab=NULL)
      aspect3d(1)
      segments3d(x = as.vector(t(seg[, c(1, 4)])), y = as.vector(t(seg[, 
                        c(2, 5)])), z = as.vector(t(seg[, c(3, 6)])), lwd=2, col="gray40")
     
      text3d(x = seg[, 4], y = seg[, 5], z = seg[, 6], 
             1:n, font=2, cex=1.5)
      
      if(showlegend==TRUE){
        legend3d("left", legend=colnames(occ)[2:(n+1)],
             pch = paste(seq(1:n)), cex=1)
      }
    }
  }
  return(tran2)
}
