check.one.level2 <-
function(X, Y, cond, sample, ncomp, keepX, keepY, method, tab.prob.gene,
max.iter, tol, ...)
{
    X = as.matrix(X)
    if (length(dim(X)) != 2 || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
    if (is.null(cond))
    stop("Vector cond is missing", call. = FALSE)
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of components, 'ncomp'.")
    if (is.null(sample))
    stop("Vector sample is missing", call. = FALSE)
    if (!is.null(dim(sample)))
    stop("sample should be a vector indicating the repeated measurements")
    if (length(sample) != nrow(X))
    stop("X and the vector sample should have the same number of subjects")
    if (length(summary(as.factor(sample))) == nrow(X))
    stop("Check that the vector sample reflects the repeated measurements")
    if (!any(names(summary(as.factor(sample))) == "1")) {
        cat("The vector sample includes the values: ", as.vector(names(summary(as.factor(sample)))),
        "\n")
        stop("sample vector", call. = FALSE)
    }
    if (is.null(method))
    stop("Input method missing, should be set to splsda or spls",
    call. = FALSE)
}
