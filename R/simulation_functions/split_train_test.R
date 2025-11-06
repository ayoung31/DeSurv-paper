split_train_test <- function(n, train_frac = 0.7) {
  idx <- sample(seq_len(n))
  n_train <- floor(train_frac * n)
  list(
    train = idx[seq_len(n_train)],
    test  = idx[(n_train+1):n]
  )
}