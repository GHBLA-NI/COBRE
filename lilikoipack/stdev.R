samplings_stdev<-function (m, n, attempts, z, ranks, samplings, start, logfile = "") 
{
    dall <- array(dim = c(attempts, n))
    skip <- 0
    for (a in 1:attempts) {
        res <- score_pathway(z[samplings[a, ], ], m, ranks[samplings[a, 
            ]], start = start, logfile = logfile)
        if (!is.null(res)) {
            dall[a, samplings[a, ]] <- res$score
        }
        else {
            skip <- skip + 1
        }
    }
    if (skip < attempts/2) {
        mean(apply(dall, 2, sd, na.rm = TRUE), na.rm = TRUE)
    }
    else {
        Inf
    }
}
