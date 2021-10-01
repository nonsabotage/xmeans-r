# ***********************************************************
# xmeans_sub.r
# ***********************************************************
xmeans_sub <- setRefClass (
	Class = "xmeans_sub",

	methods = list (
		addPackages = function (CRAN = "http://cran.ism.ac.jp/",
								more = NULL)
		{
			addPackages <- c ("dplyr", "matrixStats", "mvtnorm", more)
			toInstPacks <- setdiff (addPackages, .packages(all.available = TRUE))
			for (pkg in toInstPacks)
			    install.packages (pkg, repos = CRAN)
			for (pkg in addPackages)
			    suppressPackageStartupMessages (library (pkg, character.only = TRUE))
		},
		xBIC = function (dat_, ignore.covar = TRUE)
		{
			dat <- dat_[, sapply (dat_, is.numeric)]
			nc  <- ncol (dat)
			nr  <- nrow (dat)
			if (ignore.covar) {
				covar <- dat %>%
							 summarise_each (funs (var)) %>%
							 diag ()
				df    <- nc * 2
			} else {
				covar <- cov (dat)
				df    <- nc * (nc + 3) * .5
			}
			center <- colMeans (dat)
			bic    <- -2 * sum (dmvnorm (dat, center, covar, log = TRUE)) + df * log (nr)
			return (bic)
		},
		xBICp = function (dat_, cluster_, ignore.covar = TRUE)
		{
			dat     <- dat_[, sapply (dat_, is.numeric)]
			cluster <- unlist (cluster_)
			nc      <- ncol (dat)
			nr      <- nrow (dat)
			if (ignore.covar) {
				covar_eachCluster <- cbind (dat, cc = cluster) %>%
									     group_by (cc) %>%
										 summarise_each (funs (var)) %>%
										 dplyr::select (., -cc) %>%
										 rowwise () %>%
										 do (rtn = diag (.)) %>%
										 as.list () %>%
										 "$"(rtn)
				sumDeterminant    <- sum (sapply (covar_eachCluster, function (mat) {prod (diag (mat))}))
				df                <- 4 * nc
			}
			else {
				covar_eachCluster <- cbind (dat, cc = cluster) %>%
										  group_by (cc) %>%
										  do (rtn = cov (dplyr::select (., -cc))) %>%
										  as.list () %>%
										  "$"(rtn)
				sumDeterminant    <- sum (sapply (covar_eachCluster, function (mat) {det (mat)}))
				df                <- nc * (nc + 3)
			}
			center_eachCluster <- cbind (dat, cc = cluster) %>%
									  group_by (cc) %>%
									  do (rtn = colMeans (dplyr::select (., -cc))) %>%
							          as.list () %>%
							          "$"(rtn)
			beta     <- sqrt (c (crossprod (center_eachCluster [[1]] - center_eachCluster [[2]])) / sumDeterminant)
			logAlpha <- - log (2) - pnorm (beta, log.p = TRUE)

			logLikelihood_eachCluster <- sapply (seq_along (unique (cluster)),
									 		     function (k) {
									 				 center_ <- center_eachCluster[[k]]
									 				 covar_  <- covar_eachCluster[[k]]
									 				 dat_    <- dat[cluster == k, ]
									 				 sum (dmvnorm (dat_, center_, covar_, log = TRUE))
												 })
			bicp <- - 2 * (nr * logAlpha + sum (logLikelihood_eachCluster)) + df * log (nr)
			return (bicp)
		},
		split2Cluster = function (dat_, ignore.covar = TRUE)
		{
			skm      <- s_kmeans$new ()
			cluster  <- skm$simple_kmeans (dat_, K_ = 2)
			ccounts  <- table (cluster)
			if (any (ccounts < 3)) {
				return (list (doSPLIT = FALSE))
			}
			bicp <- .self$xBICp (dat_, cluster = cluster)
			bic  <- .self$xBIC  (dat_)
			if (bic > bicp) {
				return (list (doSPLIT = TRUE,
							  d1      = dat_[cluster == 1, ],
							  d2      = dat_[cluster == 2, ]))
			}
			else {
				return (list (doSPLIT = FALSE))
			}
		},
		merge2Cluster = function (dat_i, dat_j)
		{
			par     <- data.frame (len = c (nrow (dat_i), nrow (dat_j)),
							       val = 1:2)
			cluster <- par %>%
						   rowwise () %>%
						   do (rval = rep (.$val, each = .$len)) %>%
						   "$"(rval) %>%
						   unlist ()
			dat_ij  <- rbind (dat_i, dat_j)
			bicp    <- .self$xBICp (dat_ij, cluster)
			bic     <- .self$xBIC (dat_ij)
			if (bic < bicp)
				return (list (doMERGE = TRUE,
							  D       = dat_ij))
			else
				return (list (doMERGE = FALSE))
		}
	)
)
