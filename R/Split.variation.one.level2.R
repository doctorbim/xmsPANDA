Split.variation.one.level2 <-
function (X, Y, sample)
{
    X = as.matrix(X)
    if (is.factor(sample)) {
        sample = as.numeric(sample)
        warning("the vector sample was converted into a numeric vector",
        call. = FALSE)
    }
    Xmi <- colMeans(X)
    Xm <- matrix(Xmi, nrow = nrow(X), ncol = ncol(X), byrow = T)
    indX <- cbind(sample, X)
    indsample <- unique(sample)
    n.sample <- length(indsample)
    Xbi <- t(apply(matrix(indsample, ncol = 1, nrow = n.sample),
    MARGIN = 1, FUN = function(x, indX) {
        indice <- which(indX[, 1] == x[1])
        res <- colMeans(indX[indice, ])[-1]
        return(c(x, res))
    }, indX = indX))
    Xb <- apply(matrix(sample, ncol = 1, nrow = length(sample)),
    MARGIN = 1, FUN = function(x, Xbi) {
        Xbi[which(Xbi[, 1] == x), -1]
    }, Xbi = Xbi)
    Xb <- t(Xb) - Xm
    Xw <- X - Xm - Xb
    res <- list(Xw = Xw, Xb = Xb, Xm = Xm)
}
