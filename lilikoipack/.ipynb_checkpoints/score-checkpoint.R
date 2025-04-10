score_all_pathways_helper<-function (z, ranks, samplings, i, attempts, maximize_stability, 
    logfile = "", start) 
{
    n <- dim(z)[1]
    k <- dim(z)[2]
    m <- dim(samplings)[2]
    mincheck <- 5
    kmin = max(floor(0.8 * k), mincheck + 1)
    mindelta = min(0.009, max(0.002, 1.5/k))
    sig <- matrix(0, 1, k)
    res <- score_pathway(z, n, ranks, calcerr = TRUE, start = start, 
        logfile = logfile)
    if (is.null(res)) {
        cat(file = logfile, append = TRUE, "pathway ", i, "> scoring failed 1.\n")
    }
    else {
        sig <- samplings_stdev(m, n, attempts, z, ranks, samplings, 
            start = start)
        if (sig > 10000) {
            cat(file = logfile, append = TRUE, "pathway ", i, 
                "> scoring failed 2 (sig:", sig, ").\n")
            res <- NULL
        }
        else {
            origsig <- sig
            cat(file = logfile, append = TRUE, "pathway ", i, 
                "> sig:", sig, "\n")
            isin <- 1:k
            if (maximize_stability) {
                testsig <- max(mincheck, floor(0.1 * k))
                newsig <- rep(0, testsig)
                while ((k >= kmin) & (sig > 0.05)) {
                  se <- sort(res$error, index.return = TRUE, 
                    decreasing = TRUE)
                  for (j in 1:testsig) {
                    newsig[j] <- samplings_stdev(m, n, attempts, 
                      z[, -se$ix[j]], ranks, samplings, start = start)
                  }
                  wj <- which.min(newsig)
                  cat(file = logfile, append = TRUE, "pathway ", 
                    i, " k=", k, "(", ncol(res$thecurve$s), ") wj=", 
                    wj, ">new sig:", newsig[wj])
                  if (sig - newsig[wj] < mindelta) {
                    cat(file = logfile, append = TRUE, " x rejected\n")
                    break
                  }
                  cat(file = logfile, append = TRUE, " | accepted!\n")
                  sig <- newsig[wj]
                  isin <- isin[-se$ix[wj]]
                  z <- z[, -se$ix[wj]]
                  k <- k - 1
                  res <- score_pathway(z, n, ranks, calcerr = TRUE, 
                    start = start, logfile = logfile)
                  if (is.null(res)) {
                    cat(file = logfile, append = TRUE, "pathway ", 
                      i, "> scoring failed 3.\n")
                    break
                  }
                }
            }
        }
    }
    if (is.null(res)) {
        NULL
    }
    else {
        list(score = res$score, thecurve = res$thecurve, z = z, 
            isin = isin, sig = sig, origsig = origsig, k = k)
    }
}
