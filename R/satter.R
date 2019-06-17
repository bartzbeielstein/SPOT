###################################################################################################
#' Satterthwaite Function 
#' 
#' The Satterthwaite function can be used to estimate the magnitude of the variance component (sigma_beta)^2, 
#' when the random factor has significant main effects.
#' 
#' Note, the output from the \code{satter()} procedure is sigma_beta. 
#' 
#' @param MScoeff coefficients c_1, c_2
#' @param MSi mean squared values
#' @param dfi degrees of freedom
#' @param alpha error probability
#' 
#' @return vector with 1. estimate of variance 2. degrees of freedom, 3. lower value of 1-alpha confint
#' 4. upper value of 1-alpha confint
#'
#' @export
#' @examples
#' res <- satter(MScoeff= c(1/4, -1/4)
#'              , MSi = c(394.9, 73.3)
#'              , dfi = c(4,3)
#'              , alpha = 0.1)
###################################################################################################
satter <- function (MScoeff, MSi, dfi, alpha = 0.05)
{
    Lterm = MScoeff * MSi
    Lsum = sum(MScoeff * MSi)
    dff <- ((Lsum^2)/sum((Lterm^2)/dfi))
    dff2 <- round(dff)
    if (dff2 == 0) 
        dff2 = 1
    lower <- (dff2 * Lsum)/(qchisq(1 - alpha/2, dff2))
    upper <- (dff2 * Lsum)/(qchisq(alpha/2, dff2))
    return(cbind(estimate = Lsum, 
                 df = round(dff), 
                 lower = lower,
                 upper = upper))
}

