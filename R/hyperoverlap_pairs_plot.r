#' Overlap heatmap plotting for analysis of multiple entities
#'
#' This function plots a matrix of overlap.
#'
#' @param x A matrix of the form produced by produced by \code{hyperoverlap_set()} (see Details).
#' @param cols A vector of colours (default: \code{c("red","blue)}).
#' @details Input matrix must contain columns named "entity1", "entity2" and "result"
#'
#' @return A \code{ggplot} object
#'
#'@examples
#'\dontrun{
#' hyperoverlap_pairs_plot(hyperoverlap.set.iris)
#'}
#' @export

hyperoverlap_pairs_plot <-  function(x, cols = pal){

  entity1 <- entity2 <- NULL
  pal = c("red","blue","lightgrey")
  overlap.plot = ggplot2::ggplot(data = x, aes(entity1, entity2, fill = as.factor(x$result)))+
                 geom_tile(color = "white", na.rm=TRUE, lwd=3)+
                 theme_void()+
                 scale_fill_manual(values=cols)+
                 scale_x_discrete(position = "top")+
                 theme(axis.text.x = element_text(angle = 90, vjust = 0,size = 8, hjust = 0),
                       axis.text.y = element_text(vjust = 0, size = 8, hjust = 0))+
                 coord_fixed()

  overlap.plot$labels$fill <- "Result"
  return(overlap.plot)

}



