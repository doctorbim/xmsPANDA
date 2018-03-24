multilevel.splsda2 <-
function (X, cond = NULL, sample = NULL, ncomp = 2, keepX = rep(ncol(X),
ncomp), tab.prob.gene = NULL, max.iter = 500, tol = 1e-06,...)
{
    factor1 = FALSE
    factor2 = FALSE
    if (is.null(dim(cond))) {
        factor1 = TRUE
        if (!is.factor(cond)) {
            cond = as.factor(cond)
            warning("cond was set as a factor", call. = FALSE)
        }
    }
    else {
        if (ncol(cond) == 2) {
            factor2 = TRUE
            if (!is.factor(cond[, 1])) {
                warning("First cond response was set as a factor",
                call. = FALSE)
            }
            if (!is.factor(cond[, 2])) {
                warning("Second cond response was set as a factor",
                call. = FALSE)
            }
            cond1 = as.factor(cond[, 1])
            cond2 = as.factor(cond[, 2])
        }
        else {
            stop("'cond' must be a matrix with max. 2 columns for the 2 factor analysis.")
        }
    }
    if (factor1 == TRUE) {
        if (nrow(X) != length(cond))
        stop("X and cond should have the same number of subjects")
    }
    if (factor2 == TRUE) {
        if (nrow(X) != nrow(cond))
        stop("X and cond should have the same number of subjects")
    }
    if (factor1 == TRUE) {
        Xw <- Split.variation.one.level2(X, cond, sample)$Xw
    }
    else {
        Xw <- Split.variation.two.level2(X, cond1, cond2, sample)$Xw
        cond <- as.factor(paste(cond1, cond2, sep = "."))
    }
    res <- splsda(Xw, cond, ncomp = ncomp, keepX, max.iter, tol,
    ...)
    class(res) <- "list"
    if (factor1) {
        result <- c(res, list(Xw = Xw, sample = sample, name.condition = factor(cond),
        tab.prob.gene = tab.prob.gene))
        class(result) <- c("splsda1fact", "splsda")
    }
    else {
        result <- c(res, list(Xw = Xw, sample = sample, name.condition = cond1,
        tab.prob.gene = tab.prob.gene, name.time = cond2))
        class(result) <- c("splsda2fact", "splsda")
    }
    return(invisible(result))
}
