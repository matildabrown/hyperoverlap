#' Storage class for the description of hyperoverlaps
#'
#' @slot entity1 A length-one character vector
#' @slot entity2 A length-one character vector
#' @slot dimensions A length n character vector containing the variables used to define the space
#' @slot occurrences A matrix containing the labelled input data
#' @slot shape shape of the decision boundary; either "linear" or "curvilinear"
#' @slot polynomial.order a length-one numeric vector showing the polynomial order of the most accurate kernel function. "0" if linear kernel.
#' @slot result a length-one character vector, either "overlap" or "non-overlap"
#' @slot accuracy a 2x2 table with the true (y) and predicted (pred) labels
#' @slot number.of.points.misclassified a length-one numeric vector
#' @slot decision.boundary list of 2; used for plotting decision boundary
#'
#' @name hyperoverlap-class
#' @rdname hyperoverlap-class

methods::setClass(Class="Hyperoverlap",
                  methods::representation(
                    entity1="character",
                    entity2="character",
                    dimensions="character",
                    occurrences = "data.frame",
                    shape="character",
                    polynomial.order = "numeric",
                    result = "character",
                    accuracy = "table",
                    number.of.points.misclassified = "numeric",
                    decisionboundary="list"
                  )
)
