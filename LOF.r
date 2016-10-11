# ===========================================================================
#
# Faster LOF than DMwR::lofactor, Rlof::lof (a single thread)
#
# ===========================================================================

library(nabor)  # for knn

LOF <- function(data, k) {
    knns <- knn(data, k = NROW(data))

    lrds <- numeric(NROW(data))
    knnIndsList <- vector("list", NROW(data))

    for (i in 1:NROW(data)) {
	# Count the number of k-th nearest neighbor
	numneigh <- sum(knns$nn.dists[i, -1] <= knns$nn.dist[i, k + 1])
	# knn's k-nearest neighbors
	knnIndsList[[i]] <- knns$nn.idx[i, 1:numneigh + 1]  
	# rbind("k-dist", "dist")
	tmp <- rbind(knns$nn.dists[knnIndsList[[i]], k + 1], knns$nn.dists[i, 1:numneigh + 1])
	# Calculate "Reachability Distance"
	reachDists <- apply(tmp, 2, max)
	lrds[i] <- 1 / (sum(reachDists) / numneigh)
    }

    lofs <- numeric(NROW(data))
    for (i in 1:NROW(data)) {
	lofs[i] <- sum(lrds[knnIndsList[[i]]] / lrds[i]) / length(knnIndsList[[i]])
    }

    return(lofs)
}


## example
if (0) {
    library(DMwR)
    library(Rlof)

    n <- 300
    d <- 300
    set.seed(1)
    dat <- matrix(rnorm(n * d), n, d)
    dat <- rbind(dat[1, ], dat[1, ], dat[4, ], dat)

    k <- 50

    # Check!
    print(system.time(lof1 <- LOF(dat, k)))
    print(system.time(lof2 <- lofactor(dat, k)))
    print(system.time(lof3 <- lof(dat, k, cores = 1)))
    print(identical(lof1, lof2))
    print(identical(lof1, lof3))
}
