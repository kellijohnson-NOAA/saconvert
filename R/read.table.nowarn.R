read.table.nowarn=function (...)
{
  tryCatch.W.E <- function(expr) {
    W <- NULL
    w.handler <- function(w) {
      if (!grepl("incomplete final line", w))
        W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler), warning = W)
  }
  lis <- tryCatch.W.E(utils::read.table(...))
  if (!is.null(lis$warning))
    warning(lis$warning)
  lis$value
}
is.whole.positive.number=function (x, tol = .Machine$double.eps^0.5)
{
  (abs(x - round(x)) < tol) & (x >= 0)
}
