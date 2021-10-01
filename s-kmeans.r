# ***********************************************************
# s-kmeans.r
# ***********************************************************
s_kmeans <- setRefClass (

	Class  = "s_kmeans",

	methods = list (
		initCluster_random = function (dat, NUM_CLUSTER)
		{
			index        <- 1:nrow (dat)
			center_index <- sample (index, NUM_CLUSTER, replace = FALSE)
			cur_center   <- list ()
			for (i in seq_len (NUM_CLUSTER))
				cur_center[[i]] <- unlist (dat[center_index[i],])
			dist <- sapply (seq_len (NUM_CLUSTER),
							function (k) {
								center_ <- cur_center[[k]]
								rowSums ((sweep (dat, 2, center_, "-")) **2)
							})
			cur_clusters <- max.col (-dist)
			return (list (cur_center   = cur_center,
						  cur_clusters = cur_clusters))
		},
		initCluster_kmeanspp = function (dat, NUM_CLUSTER)
		{
			index        <- 1:nrow (dat)
			center_index <- sample (index, 1)
			center_      <- unlist (dat[center_index, ])
			min_dist     <- rowSums ((sweep (dat, 2, center_, "-")) ** 2)
			while (length (center_index) < NUM_CLUSTER) {
				prob         <- min_dist[-center_index] / sum (min_dist[-center_index])
				new_center   <- sample (index[-center_index], 1, prob = prob)
				center_      <- unlist (dat[new_center, ])
				new_dist     <- rowSums ((sweep (dat, 2, center_, "-")) ** 2)
				min_dist     <- pmin (min_dist, new_dist)
				center_index <- c (center_index, new_center)
			}
			cur_center <- list ()
			for (i in seq_len (NUM_CLUSTER))
				cur_center[[i]] <- unlist (dat[center_index[i],])
			dist <- sapply (seq_len (NUM_CLUSTER),
							function (k) {
								center_ <- cur_center[[k]]
								rowSums ((sweep (dat, 2, center_, "-")) **2)
							})
			cur_clusters <- max.col (-dist)
			return (list (cur_center   = cur_center,
						  cur_clusters = cur_clusters))
		},
		updateCenter = function (dat, cur_clusters)
		{
			cur_center <-
			 	cbind (dat, cc = cur_clusters) %>%
				group_by (cc) %>%
				do (rtn = colMeans (dplyr::select (., -cc))) %>%
				as.list () %>%
				"$"(rtn)
			return (cur_center)
		},
		simple_kmeans = function (dat_, K_ = 2, kmeanspp = TRUE)
		{
			dat          <- na.omit (dat_[, sapply (dat_, is.numeric)])
			NUM_CLUSTER  <- as.integer (K_)
			SIZE         <- as.integer (nrow (dat))
			if (kmeanspp)
				tmp <- .self$initCluster_kmeanspp (dat, NUM_CLUSTER)
			else
				tmp <- .self$initCluster_random (dat, NUM_CLUSTER)
			cur_clusters <- tmp$cur_clusters
			cur_center   <- tmp$cur_center
			old_clusters <- cur_clusters
			convergence  <- FALSE
			iter         <- 1L
			while (!convergence)
			{
				iter <- iter + 1L
				if (iter > 100L) {
					cat ("iter reached MAX_ITER\n")
					break
				}
				distances <-
					sapply (seq_len (NUM_CLUSTER),
							function (k) {
								center_ <- cur_center[[k]]
								d_      <- sweep (dat, 2, center_, "-")
								return (rowSums (d_ * d_))
							})
				cur_clusters  <- max.col (- distances)
				print (cur_clusters)
				if (length (unique (cur_clusters)) != NUM_CLUSTER)
					if (kmeanspp)
						cur_clusters <- .self$initCluster_kmeanspp (dat, NUM_CLUSTER)$cur_clusters
					else
						cur_clusters <- .self$initCluster_random (NUM_CLUSTER, SIZE)$cur_clusters
				convergence   <- identical (old_clusters, cur_clusters)
				old_clusters  <- cur_clusters
				cur_center    <- .self$updateCenter (dat, cur_clusters)
			}
			return (as.integer (cur_clusters))
		}
	)
)