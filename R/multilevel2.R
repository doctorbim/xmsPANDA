multilevel2 <-
function (X, Y = NULL, cond = NULL, sample = NULL, ncomp = 1,
keepX = rep(ncol(X), ncomp), keepY = NULL, method = NULL,
tab.prob.gene = NULL, max.iter = 500, tol = 1e-06, ...)
{
    check.one.level2(X, Y, cond, sample, ncomp, keepX, keepY,
    method, tab.prob.gene, max.iter, tol, ...)
    if (is.factor(sample)) {
        sample = as.numeric(sample)
        warning("the vector sample was converted into a numeric vector",
        call. = FALSE)
    }
    if (method == "splsda") {
        result = multilevel.splsda2(X = X, cond = cond, sample = sample,
        ncomp = ncomp, keepX = keepX, tab.prob.gene = tab.prob.gene)
    }
    else {
        result = multilevel.spls2(X = X, Y = Y, cond = cond, sample = sample,
        ncomp = ncomp, keepX = keepX, keepY = keepY, tab.prob.gene = tab.prob.gene)
    }
    return(result)
}
