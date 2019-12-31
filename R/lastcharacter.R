#' Test equality of last character to a character value
#'
#' Find the last character of each string and test if it matches
#' the character provided in the \code{y} argument.
#' 
#' @param x Vector of strings
#' @param y A single character to test for equality to
#'
lastcharacter <- function(x, y) {
  splits <- strsplit(x, split = "")
  ends <- lapply(lapply(splits, rev), "[", 1)
  return(ends == y)
}
