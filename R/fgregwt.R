
fcalibrate <- function(DT, target, cols, eps = 0.005, max_iter = 100) {
  .Call("C_fcalibrate_gregwt",
        DT,
        target,
        cols,
        eps,
        max_iter,
        PACKAGE = packageName())
}

fbcalibrate <- function(DT, target, cols, eps = 0.005, max_iter = 100) {
  .Call("C_bit_calibrate_gregwt",
        DT,
        target,
        cols,
        eps,
        max_iter,
        PACKAGE = packageName())
}
