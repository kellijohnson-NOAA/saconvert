read.surveys=function (filen)
{
  lin <- readLines(filen, warn = FALSE)[-c(1:2)]
  empty <- which(lapply(lapply(strsplit(lin, split = "[[:space:]]+"),
                               paste, collapse = ""), nchar) == 0)
  if (length(empty) > 0) {
    lin <- lin[-empty]
  }
  lin <- sub("^\\s+", "", lin)
  idx1 <- grep("^[A-Z#]", lin, ignore.case = TRUE)
  idx2 <- c(idx1[-1] - 1, length(lin))
  names <- lin[idx1]
  years <- matrix(as.numeric(unlist(strsplit(lin[idx1 + 1],
                                             "[[:space:]]+"))), ncol = 2, byrow = TRUE)
  times <- matrix(as.numeric(unlist(strsplit(lin[idx1 + 2],
                                             "[[:space:]]+"))), ncol = 4, byrow = TRUE)[, 3:4, drop = FALSE]
  ages <- matrix(as.numeric(unlist(lapply(strsplit(lin[idx1 +
                                                         3], "[[:space:]]+"), function(x) x[1:2]))), ncol = 2,
                 byrow = TRUE)
  for (i in 1:length(names)) {
    if (!is.whole.positive.number(years[i, 1])) {
      stop(paste("In file", filen, ": Minimum year is expected to be a positive integer number for fleet number",
                 i))
    }
    if (!is.whole.positive.number(years[i, 2])) {
      stop(paste("In file", filen, ": Maximum year is expected to be a positive integer number for fleet number",
                 i))
    }
    if (years[i, 1] > years[i, 2]) {
      stop(paste("In file", filen, ": Maximum year is expected to be greater than minimum year for fleet number",
                 i))
    }
    if (ages[i, 1] > ages[i, 2]) {
      stop(paste("In file", filen, ": Maximum age is expected to be greater than minimum age for fleet number",
                 i))
    }
    if ((times[i, 1] < 0) | (times[i, 1] > 1)) {
      stop(paste("In file", filen, ": Minimum survey time is expected to be within [0,1] for fleet number",
                 i))
    }
    if ((times[i, 2] < 0) | (times[i, 2] > 1)) {
      stop(paste("In file", filen, ": Maximum survey time is expected to be within [0,1] for fleet number",
                 i))
    }
    if (times[i, 2] < times[i, 1]) {
      stop(paste("In file", filen, ": Maximum survey time is expected greater than minimum survey time for fleet number",
                 i))
    }
  }
  as.num <- function(x, na.strings = "NA") {
    stopifnot(is.character(x))
    na = x %in% na.strings
    x[na] = 0
    x = as.numeric(x)
    x[na] = NA_real_
    x
  }
  onemat <- function(i) {
    lin.local <- gsub("^[[:blank:]]*", "", lin[(idx1[i] +
                                                  4):idx2[i]])
    nr <- idx2[i] - idx1[i] - 3
    ret <- matrix(as.num(unlist((strsplit(lin.local, "[[:space:]]+")))),
                  nrow = nr, byrow = TRUE)[, , drop = FALSE]
    if (nrow(ret) != (years[i, 2] - years[i, 1] + 1)) {
      stop(paste("In file", filen, ": Year range specified does not match number of rows for survey fleet number",
                 i))
    }
    if ((ncol(ret) - 1) < (ages[i, 2] - ages[i, 1] + 1)) {
      stop(paste("In file", filen, ": Fewer columns than indicated by age range for survey fleet number",
                 i))
    }
    if (!is.numeric(ret)) {
      stop(paste("In file", filen, ": Non numeric data values detected for survey fleet number",
                 i))
    }
    ret <- as.matrix(ret[, -1]/ret[, 1])
    rownames(ret) <- years[i, 1]:years[i, 2]
    ret <- ret[, 1:length(ages[i, 1]:ages[i, 2]), drop = FALSE]
    colnames(ret) <- ages[i, 1]:ages[i, 2]
    attr(ret, "time") <- times[i, ]
    ret[ret < 0] <- NA
    ret
  }
  obs <- lapply(1:length(names), onemat)
  names(obs) <- names
  obs
}
