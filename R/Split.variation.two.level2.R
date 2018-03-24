Split.variation.two.level2 <-
function (X, factor1, factor2, sample)
{
    if (is.factor(sample)) {
        sample = as.numeric(sample)
        warning("the vector sample was converted into a numeric vector",
        call. = FALSE)
    }
    Xmi <- colMeans(X)
    Xm <- matrix(Xmi, nrow = nrow(X), ncol = ncol(X), byrow = T)
    indX <- cbind(sample, X)
    Xb <- apply(indX, MARGIN = 1, FUN = function(x, indX) {
        indice <- which(indX[, 1] == x[1])
        res <- colMeans(indX[indice, ])
        return(res[-1])
    }, indX = indX)
    Xs <- t(Xb)
    Xb <- t(Xb) - Xm
    xbfactor1 <- X
    for (i in levels(factor(factor1))) {
        indice <- which(factor1 == i)
        indXX <- indX[indice, ]
        res1 <- apply(indXX, MARGIN = 1, FUN = function(x, indXX) {
            indice <- which(indXX[, 1] == x[1])
            if (length(indice) == 1) {
                res <- colMeans(matrix(indXX[indice, ], nrow = 1,
                ncol = dim(indXX)[2]))
            }
            else {
                res <- colMeans(indXX[indice, ])
            }
            return(res[-1])
        }, indXX = indXX)
        xbfactor1[indice, ] <- t(res1)
    }
    xbfactor2 <- X
    for (i in levels(factor(factor2))) {
        indice <- which(factor2 == i)
        indXX <- indX[indice, ]
        res1 <- apply(indXX, MARGIN = 1, FUN = function(x, indXX) {
            indice <- which(indXX[, 1] == x[1])
            if (length(indice) == 1) {
                res <- colMeans(matrix(indXX[indice, ], nrow = 1,
                ncol = dim(indXX)[2]))
            }
            else {
                res <- colMeans(indXX[indice, ])
            }
            return(res[-1])
        }, indXX = indXX)
        xbfactor2[indice, ] <- t(res1)
    }
    matfactor1 <- matrix(factor1, nrow = 1, ncol = length(factor1))
    XFACTOR1 <- apply(matfactor1, MARGIN = 2, FUN = function(x,
    matfactor1) {
        indice <- which(matfactor1 == x[1])
        res <- colMeans(X[indice, ])
        return(res)
    }, matfactor1 = matfactor1)
    matfactor2 <- matrix(factor2, nrow = 1, ncol = length(factor2))
    XFACTOR2 <- apply(matfactor2, MARGIN = 2, FUN = function(x,
    matfactor2) {
        indice <- which(matfactor2 == x[1])
        res <- colMeans(X[indice, ])
        return(res)
    }, matfactor2 = matfactor2)
    XCS <- xbfactor1 - Xs + Xm - t(XFACTOR1)
    XTS <- xbfactor2 - Xs + Xm - t(XFACTOR2)
    Xw <- X - Xb - Xm - XCS - XTS
    res <- list(Xw = Xw, Xb = Xb, Xm = Xm, XCS = XCS, XTS = XTS)
}
