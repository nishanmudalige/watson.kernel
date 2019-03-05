library(hypergeo)
library(fAsianOptions)

# Modified by examining code in the function vmfkde.tune in the Directional Package
# Credit to the following authors of the Directional package: 
# Michail Tsagris, Giorgos Athineou, Anamul Sajib, Eli Amson, Micah J. Waldstein

mwatson.tune = function (x, low = 0.1, up = 1) {
  p <- dim(x)[2]
  n <- dim(x)[1]
  d <- tcrossprod(x) * tcrossprod(x)
  
  # Square here to avoid NAN warning
  diag(d) <- NA
  con <- 2*(pi)^(p/2)
  
  funa <- function(h) {
    A <- d/(h^2)
    mpk <- gamma(p/2)/( con * Re(fAsianOptions::kummerM(1/h^2, 1/2, p/2)) )
    f <- rowSums(exp(A + log(mpk)), na.rm = TRUE)/(n - 1)
    mean(log(f))
  }
  
  a <- optimize(funa, c(low, up), maximum = TRUE)
  res <- c(a$maximum, a$objective)
  names(res) <- c("Optimal h", "cv")
  res
}


# Modified by examining code in the function vmf.kde in the Directional Package
# Credit to the following authors of the Directional package: 
# Michail Tsagris, Giorgos Athineou, Anamul Sajib, Eli Amson, Micah J. Waldstein

w.kde = function (x, h = NULL) {
  p <- dim(x)[2]
  n <- dim(x)[1]
  
  if (is.null(h)) {
    h <- as.numeric(mwatson.tune(x, low = 0.1, up = 1)[1])
  } else {
    h <- h
  }
  
  d <- (tcrossprod(x) * tcrossprod(x))/h^2
  
  mpk <- gamma(p/2)/( 2*(pi)^(p/2) * Re(fAsianOptions::kummerM(1/h^2, 1/2, p/2)) )
  f <- Rfast::rowmeans(exp(d + log(mpk) ))
  list(h = h, f = f)
}
