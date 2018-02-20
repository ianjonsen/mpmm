##' @method print mpmm
##' @export
print.mpmm <-
  function(x, digits = max(3, getOption("digits") - 3),
          ...)
  {
    print(x$call); cat("\n")

    invisible(x)
  }
