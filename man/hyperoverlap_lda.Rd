% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyperoverlap_lda.r
\name{hyperoverlap_lda}
\alias{hyperoverlap_lda}
\title{Hyperoverlap visualisation using linear discriminant analysis (LDA)}
\usage{
hyperoverlap_lda(x, return.plot=TRUE, visualise3d=FALSE, showlegend=TRUE)
}
\arguments{
\item{x}{An \code{\link{hyperoverlap-class}} object.}

\item{return.plot}{Logical. If TRUE, data are plotted using \code{plot()}.}

\item{visualise3d}{Logical. If FALSE, data are projected onto two axes (LDA1, residualPCA1). If TRUE, data are projected onto three axes (LDA1, residualPCA1, residualPCA2)}

\item{showlegend}{Logical. Used for 3D plots.}
}
\value{
Returns a dataframe with columns "Entity", "LDA1", "residualPCA1",
"residualPCA2" (if \code{visualise3d = TRUE})
}
\description{
Hyperoverlap visualisation using linear discriminant analysis (LDA)
}
\details{
This function provides a way to visualise overlap (or non-overlap)
between classes of high dimensional data. For inspection, it is useful to
use the base graphics package (implemented by return.plot=TRUE). The
transformed coordinates of each point are also returned as a dataframe,
which can be plotted with user-defined parameters.
}
\examples{
#using iris dataset reduced to two species
data = iris[which(iris$Species!=("versicolor")),]
x = hyperoverlap_detect(data[1:4], data$Species)
hyperoverlap_lda(x)

}
\seealso{
\code{\link{hyperoverlap_detect}}
}
