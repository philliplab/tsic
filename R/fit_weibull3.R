#' Fits a 3 parameter Weibull distribution to percentiles
#'
#' The Delaney paper reports the 25th, 50th, 75th and 99th percentiles of the window periods of various assays. This function will fit a three parameter Weibull distribution to these four points so that it can be included in the library.
#'
#' @param x The percentiles. The actual value that correspond to the y'th percentile.
#' @param y The percentile values (between 0 and 1)
#' @export

fit_weibull3 <- function(x, y){
  if (FALSE){
    # APTIMA
    x = c(11.5, 5.3,   29.1,  33)
    y = c(0.5,  0.025, 0.975, 0.99)
    # $par
    # [1] 4.764813 1.306748 8.916144

    # ARCHITECT
    x = c(13.6, 17.9, 23.1, 41.6)
    y = c(0.25,  0.5, 0.75, 0.99)
    # $par
    # [1]  7.209124  1.724618 13.185190

    # Geenius Fully Reactive
    x = c(28.2, 32.9, 38.6, 57.7)
    y = c(0.25,  0.5, 0.75, 0.99)
    #$par
    #[1] 21.151093  1.733457 14.483445
  }
  mse <- function(params){
    location = params[1]
    shape = params[2]
    scale = params[3]
    return (sum((y-(1 - exp(-((x - location)/scale)^shape)))^2))^(1/2)
  }
  mse(c(4.8, 1.35, 9))
  
  nlminb(c(4.8, 1.35, 9), mse)
}
