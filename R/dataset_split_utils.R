# Helpers for combining datasets and creating deterministic train/test splits.

split_sample_indices <- function(samp_info,
                                 prop_train,
                                 seed,
                                 stratify_by = NULL) {
  stopifnot(is.data.frame(samp_info), nrow(samp_info) > 0)
  stopifnot(is.numeric(prop_train), prop_train > 0, prop_train < 1)

  strata <- if (is.null(stratify_by)) {
    rep("all", nrow(samp_info))
  } else {
    stopifnot(
      is.character(stratify_by),
      length(stratify_by) == 1,
      stratify_by %in% colnames(samp_info)
    )
    samp_info[[stratify_by]]
  }

  set.seed(seed)
  train_idx <- vector("integer")
  seq_all <- seq_len(nrow(samp_info))

  for (level in unique(strata)) {
    level_idx <- which(strata == level)
    if (length(level_idx) == 1L) {
      # Nothing to stratify, assign the lone sample to training set.
      train_idx <- c(train_idx, level_idx)
      next
    }
    n_train_level <- max(1L, round(length(level_idx) * prop_train))
    n_train_level <- min(n_train_level, length(level_idx) - 1L)
    if (n_train_level <= 0L) {
      n_train_level <- 1L
    }
    selected <- sample(level_idx, n_train_level)
    train_idx <- c(train_idx, selected)
  }

  train_idx <- sort(unique(train_idx))
  val_idx <- setdiff(seq_all, train_idx)

  list(train = train_idx, val = val_idx)
}

subset_data_by_indices <- function(data, sample_indices, name_suffix = NULL) {
  stopifnot(is.list(data), !is.null(data$ex), !is.null(data$sampInfo))
  stopifnot(length(sample_indices) >= 1L)

  new_data <- list()
  new_data$ex <- data$ex[, sample_indices, drop = FALSE]
  new_data$featInfo <- data$featInfo
  new_data$sampInfo <- droplevels(data$sampInfo[sample_indices, , drop = FALSE])
  new_data$dataname <- if (is.null(name_suffix)) {
    data$dataname
  } else {
    paste(data$dataname, name_suffix, sep = "_")
  }
  new_data$samp_keeps <- which(new_data$sampInfo$keep == 1L)
  new_data
}

split_data_bundle <- function(data,
                              prop_train,
                              seed,
                              stratify_by = "dataset",
                              labels = c(train = "train", val = "val")) {
  indices <- split_sample_indices(
    samp_info = data$sampInfo,
    prop_train = prop_train,
    seed = seed,
    stratify_by = stratify_by
  )

  list(
    train = subset_data_by_indices(data, indices$train, labels[["train"]]),
    val = subset_data_by_indices(data, indices$val, labels[["val"]]),
    split = indices
  )
}
