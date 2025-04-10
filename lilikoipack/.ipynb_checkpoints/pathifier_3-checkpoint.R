PDSchanged<-function (data, allgenes, syms, 
    pathwaynames, normals = NULL, ranks = NULL, attempts = 100, 
    maximize_stability = TRUE, logfile = "", samplings = NULL, 
    min_exp = 4, min_std = 0.4) 
{
    cat(file = logfile, append = FALSE, "robust_score_bydist. min_exp=", 
        min_exp, ", min_std=", min_std, "\n")
    data[data < min_exp] = min_exp
    n <- ncol(data)
    if (is.null(normals)) {
        normals <- rep(TRUE, n)
        start <- "by pca"
    }
    else {
        start <- "by ranks"
    }
    if (is.null(ranks)) 
        ranks <- !normals
    ranks <- rank(ranks)
    if ((length(normals) != n) || (length(ranks) != n)) {
        stop("invalid dimentions")
    }
    l <- length(syms)
    nn <- sum(normals)
    m <- floor(0.8 * (n - nn)) + nn
    if (is.null(samplings)) {
        samplings <- matrix(0, attempts, m)
        w <- which(!normals)
        for (a in 1:attempts) {
            samplings[a, ] <- sort(c(w[sample(n - nn, m - nn)], 
                which(normals)))
        }
    }
    s <- NULL
    ind <- NULL
    for (i in 1:l) {
        pathway <- syms[[i]]
        pathwayindata <- getpathway(pathway, allgenes, data)
        k1 = sum(pathwayindata$isin)
        # Change the tolerance of overlapping metabolites
	if (k1 < 5) {
            si <- NULL
            cat(file = logfile, append = TRUE, "skipping pathway ", 
                i, " k1=", k1, "\n")
        }
        else {
            x <- pathwayindata$x
            pathway <- pathway[pathwayindata$isin]
            xm <- colMeans(x[normals, ])
            xs <- apply(x[normals, ], 2, sd)
            xs[xs < min_std] = min_std
            if (0 %in% xs) {
                si <- NULL
                cat(file = logfile, append = TRUE, "skipping pathway ", 
                  i, " (0 in xs)\n")
            }
            else {
                z <- (x - matrix(rep(xm, each = n), nrow = n))/(matrix(rep(xs, 
                  each = n), nrow = n))
                t <- prcomp(z)
                k2 = max(sum(t$sdev > 1.1), 4)
                k2 = min(k2, k1, 0.75 * dim(x)[1], sum(t$sdev > 
                  0.25))
                if (k2 < 3) {
                  si <- NULL
                  cat(file = logfile, append = TRUE, "skipping pathway ", 
                    i, " k2=", k2, "\n")
                }
                else {
                  pca <- t$x[, 1:k2]
                  res <- score_all_pathways_helper(pca, ranks, 
                    samplings, i, attempts, maximize_stability, 
                    logfile, start = start)
                  if (is.null(res)) {
                    si <- NULL
                    cat(file = logfile, append = TRUE, "skipping pathway ", 
                      i, "\n")
                  }
                  else {
                    ind <- c(ind, i)
                    si <- list(res$score, pathway, res$sig, res$origsig, 
                      res$k, res$thecurve$s, res$thecurve$ord, 
                      res$z, res$isin, xm, xs, t$center, t$rotation, 
                      k2)
                  }
                }
            }
        }
        s <- rbind(s, si)
    }
    cat(file = logfile, append = TRUE, length(ind), "pathways processed with start=", 
        start, "\n")
    rownames(s) <- pathwaynames[ind]
    list(scores = s[, 1], genesinpathway = s[, 2], newmeanstd = s[, 
        3], origmeanstd = s[, 4], pathwaysize = s[, 5], curves = s[, 
        6], curves_order = s[, 7], z = s[, 8], compin = s[, 9], 
        xm = s[, 10], xs = s[, 11], center = s[, 12], rot = s[, 
            13], pctaken = s[, 14], samplings = samplings, sucess = ind, 
        logfile = logfile)
}
