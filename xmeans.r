# ***********************************************************
# xmeans.r
# ***********************************************************
xmeans <- setRefClass (
	Class = "xmeans",

	contains = c ("s_kmeans", "xmeans_sub"),

	fields = c ("num_split_cluster" = "integer",
				"stack"             = "list",
				"ignore.covar"      = "logical",
				"ISHIOKA"           = "logical",
				"Cluster"           = "integer"),

	methods = list (
		initialize = function (ignore.covar_ = TRUE,
							   ISHIOKA_      = TRUE)
		{
			addPackages ()
			initFields (num_split_cluster = 0L,
						stack             = list (),
						ignore.covar      = ignore.covar_,
						ISHIOKA           = ISHIOKA_)
		},
		ishioka_xmeans = function (DATA_)
		{
			.self$splitData (.self$cleanData (DATA_))
			if (ISHIOKA)
				.self$mergeSplitedData ()
			.self$assignCluster ()
		},
		splitData = function (data)
		{
			# simple x-means
			split_data <- split2Cluster (data, ignore.covar)
			if (split_data$doSPLIT) {
				splitData (split_data$d1)
				splitData (split_data$d2)
			}
			else {
				num_split_cluster          <<- num_split_cluster + 1L
				stack[[num_split_cluster]] <<- data
			}
		},
		mergeSplitedData = function () {
			ind_order           <- .self$getDataSize_order (.self$getDataSize_eachCluster ())
			num_cluster         <- length (ind_order)
			isMerge_eachCluster <- numeric (num_cluster)
			for (i in seq_len (num_cluster - 1)) {
				i_ <- ind_order[sprintf ("%d", i)]
				if (isMerge_eachCluster[i_] == 1)
					next
				for (j in (i+1):num_cluster) {
					j_ <- ind_order[sprintf ("%d", j)]
			 		if (isMerge_eachCluster[j_] == 0) {
						merge_data <- merge2Cluster (stack[[i_]], stack[[j_]])
						if (merge_data$doMERGE) {
							stack[[i_]] <<- merge_data$D
							isMerge_eachCluster[j_] <- 1
							break
						}
					}
				}
			}
			stack             <<- stack [isMerge_eachCluster == 0]
			num_split_cluster <<- length (stack)
		},
		cleanData = function (data_)
		{
			out            <- na.omit (data_ [, sapply (data_, is.numeric)])
			rownames (out) <- 1:nrow (out)
			return (out)
		},
		getDataSize_eachCluster  = function () {
			out <- sapply (stack, nrow)
			return (out)
		},
		getDataSize_order = function (DataSize_eachCluster)
		{
			out <- seq_len (length (DataSize_eachCluster))
			out <- setNames (out, order (DataSize_eachCluster))
			return (out)
		},
		assignCluster = function ()
		{
			out <- numeric (sum (.self$getDataSize_eachCluster ()))
			ind <- lapply (stack, function (d) as.integer (rownames(d)))
			for (i in seq_along (ind))
				out[ind[[i]]] <- i
			Cluster <<- as.integer (out)
		}
	)
)

