score_pathway<-function (x, m, ranks, calcerr = FALSE, thresh = 5e-04, maxit = 200, 
    start, logfile = "") 
{
    x <- x[, apply(x, 2, sd) > 0.001]
    k <- dim(x)[2]
    if (k < 3) {
        c <- NULL
        cat(file = logfile, append = TRUE, "scoring failed (k=", 
            k, ").\n")
    }
    else {
        d <- matrix(0, 1, m)
        if (start == "by pca") {
            start <- NULL
        }
        else if (start == "by ranks") {
            start <- aggregate(x, by = list(ranks), FUN = mean)
            start <- as.matrix(start[, -1])
        }
        c <- principal_curve(x, start = start, thresh = thresh, 
            maxit = maxit)
    }
    if (!is.null(c)) {
        d[c$ord[1]] = 0
        for (j in 2:m) {
            d[c$ord[j]] <- d[c$ord[j - 1]] + dist(c$s[c$ord[(j - 
                1):j], ])
        }
        d = d/d[c$ord[m]]
        if (calcerr) {
            e <- matrix(0, 1, k)
            for (i in 1:k) {
                e[i] <- mean((c$s[, i] - x[, i])^2)
            }
        }
        else {
            e <- FALSE
        }
        list(score = d, error = e, thecurve = c)
    }
    else {
        cat(file = logfile, append = TRUE, "scoring failed.\n")
        NULL
    }
}
