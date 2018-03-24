multilevel.spls2 <-
function (X, Y, cond, sample, ncomp = 2, keepX = rep(ncol(X),
ncomp), keepY = rep(ncol(Y), ncomp), tab.prob.gene = NULL,
max.iter = 500, tol = 1e-06, ...)
{
    Y = as.matrix(Y)
    if (length(dim(Y)) != 2 || !is.numeric(Y))
    stop("'Y' must be a numeric matrix.")
    mode = "canonical"
    Xw <- Split.variation.one.level(X, Y = cond, sample)$Xw
    Yw <- Split.variation.one.level(X = Y, Y = cond, sample)$Xw
    res <- spls(Xw, Yw, mode = mode, ncomp = ncomp, keepY = keepY,
    keepX = keepX)
    result <- c(res, list(Xw = Xw, Yw = Yw, sample = sample,
    name.condition = factor(cond), tab.prob.gene = tab.prob.gene))
    if (!is.null(tab.prob.gene)) {
        probeX <- result$names$X
        geneX <- result$tab.prob.gene[match(probeX, result$tab.prob.gene[,
        1]), 2]
        result$names$X <- as.character(geneX)
        colnames(result$X) <- as.character(geneX)
    }
    class(result) <- c("splslevel", "spls")
    return(invisible(result))
}
