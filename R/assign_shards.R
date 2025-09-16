assign_shards <- function(grid, n_shards) {
  stopifnot(nrow(grid) > 0, n_shards >= 1)
  shard_of <- function(k, lambda, eta, lambdaW, lambdaH) {
    h <- digest::digest2int(paste(k, lambda, eta, lambdaW, lambdaH, sep = "|"))
    h %% n_shards
  }
  grid %>% mutate(shard = purrr::pmap_int(cur_data_all(), ~shard_of(..., n_shards = n_shards)))
}