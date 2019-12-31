# code to read in ICES file structure and convert to ASAP
#  uses SAM read.ices fn modified by Dan Hennen,
read.ices=function (filen)
{
  if (grepl("^[0-9]", scan(filen, skip = 2, n = 1, quiet = TRUE,
                           what = ""))) {
    head <- scan(filen, skip = 2, n = 5, quiet = TRUE)
    minY <- head[1]
    maxY <- head[2]
    minA <- head[3]
    maxA <- head[4]
    datatype <- head[5]
    if (!is.whole.positive.number(minY)) {
      stop(paste("In file", filen, ": Minimum year is expected to be a positive integer number"))
    }
    if (!is.whole.positive.number(maxY)) {
      stop(paste("In file", filen, ": Maximum year is expected to be a positive integer number"))
    }
    if (!is.whole.positive.number(minA)) {
      stop(paste("In file", filen, ": Minimum age is expected to be a positive integer number"))
    }
    if (!is.whole.positive.number(maxA)) {
      stop(paste("In file", filen, ": Maximum age is expected to be a positive integer number"))
    }
    if (!(datatype %in% c(1, 2, 3, 5))) {
      stop(paste("In file", filen, ": Datatype code is expected to be one of the numbers 1, 2, 3, or 5"))
    }
    if (minY > maxY) {
      stop(paste("In file", filen, ": Minimum year is expected to be less than maximum year"))
    }
    if (minA > maxA) {
      stop(paste("In file", filen, ": Minimum age is expected to be less than maximum age"))
    }
    C <- as.matrix(read.table.nowarn(filen, skip = 5, header = FALSE))
    if (datatype == 1) {
      if ((maxY - minY + 1) != nrow(C)) {
        stop(paste("In file", filen, ": Number of rows does not match the year range given"))
      }
      if ((maxA - minA + 1) > ncol(C)) {
        stop(paste("In file", filen, ": Fewer columns than the age range given"))
      }
    }
    if (datatype == 2) {
      C <- as.matrix(read.table.nowarn(filen, skip = 5,
                                       header = FALSE))
      if (1 != nrow(C)) {
        stop(paste("In file", filen, ": For datatype 2 only one row of data is expected"))
      }
      if ((maxA - minA + 1) > ncol(C)) {
        stop(paste("In file", filen, ": Fewer columns than the age range given"))
      }
      C <- C[rep(1, maxY - minY + 1), ]
    }
    if (datatype == 3) {
      C <- as.matrix(read.table.nowarn(filen, skip = 5,
                                       header = FALSE))
      if (1 != nrow(C)) {
        stop(paste("In file", filen, ": For datatype 3 only one row of data is expected"))
      }
      if (1 != ncol(C)) {
        stop(paste("In file", filen, ": For datatype 3 only one column of data is expected"))
      }
      C <- C[rep(1, maxY - minY + 1), rep(1, maxA - minA +
                                            1)]
    }
    if (datatype == 5) {
      C <- as.matrix(read.table.nowarn(filen, skip = 5,
                                       header = FALSE))
      if ((maxY - minY + 1) != nrow(C)) {
        stop(paste("In file", filen, ": Number of rows does not match the year range given"))
      }
      if (1 != ncol(C)) {
        stop(paste("In file", filen, ": For datatype 5 only one column of data is expected"))
      }
      C <- C[, rep(1, maxA - minA + 1)]
    }
    rownames(C) <- minY:maxY
    C <- C[, 1:length(minA:maxA)]
    colnames(C) <- minA:maxA
    if (!is.numeric(C)) {
      stop(paste("In file", filen, ": Non numeric data values detected (could for instance be comma used as decimal operator)"))
    }
    return(C)
  }
  else {
    return(read.surveys(filen))
  }
}
