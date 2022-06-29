compute_diff_history_norm <- function(ranks, max_rank, p = 2) {
  rank_t <- rep(0, max_rank + 1)
  history_norm <- rep(NA_real_, length(ranks))
  uniform_cdf <- seq(0, 1, length.out = max_rank + 2)[2:(max_rank + 2)]
  for(i in 1:length(ranks)) {
    rank_t[ranks[i] + 1] <- rank_t[ranks[i] + 1] + 1
    diff <- uniform_cdf - cumsum(rank_t) / i
    history_norm[i] <- sum((abs(diff)^p)/(max_rank + 1))^(1/p)
  }
  history_norm
}

get_precomputed_gamma_thresholds <- function(K, min_sims = 1, max_sims = 1000) {
  if(!dir.exists("./_SBC_cache")) {
    dir.create("./_SBC_cache")
  }

  gamma_file <- paste0("./_SBC_cache/gamma_thresholds_",K,"_",min_sims,"_", max_sims, ".rds")
  if(!file.exists(gamma_file)) {
    message("Precomputing gamma thresholds")
    # Computing time differs by N, randomly permute to get optimal division of work
    N_sims <- sample(min_sims:max_sims)
    thresholds <- future.apply::future_sapply(N_sims, FUN = function(n_sims) {
      SBC:::adjust_gamma(N = n_sims, L = 1, K = K)
    })

    thresholds_df <-
      dplyr::arrange(data.frame(N_sims = N_sims, log_gamma_threshold = log(thresholds)), N_sims)

    saveRDS(thresholds_df, gamma_file)
  } else {
    thresholds_df <- readRDS(gamma_file)
  }
  thresholds_df
}

compute_log_gamma_history <- function(ranks, max_rank) {
  rank_t <- rep(0, max_rank + 1)
  log_gamma <- rep(NA_real_, length(ranks))
  dummy <- rep(NA_real_, length(ranks))

  K <- max_rank + 1
  z <- (1:(K - 1)) / K

  for(i in 1:length(ranks)) {
    rank_t[ranks[i] + 1] <- rank_t[ranks[i] + 1] + 1
    scaled_ecdf <- cumsum(rank_t[1:max_rank])

    log_gamma[i] <- log(2) + min(
      pbinom(scaled_ecdf, i, z, log = TRUE),
      pbinom(scaled_ecdf - 1, i, z, lower.tail = FALSE, log = TRUE)
    )
    #dummy[i] <- log(2) + pbinom(scaled_ecdf[max_rank + 1] - 1, i, z[max_rank + 1], lower.tail = FALSE, log = TRUE)
  }
  #print(dummy)
  log_gamma
}

compute_ks_test_history <- function(ranks, max_rank) {
  ranks_cont <- (ranks + runif(length(ranks))) / (max_rank + 1)
  ks_p <- rep(NA_real_, length(ranks))
  for(i in 1:length(ranks)) {
    ks_p[i] <- ks.test(ranks_cont[1:i], "punif")$p.value
  }
  ks_p
}


plot_log_gamma_history <- function(results, min_sim_id = 0, max_sim_id = Inf) {
  unique_max_rank <- unique(results$stats$max_rank)
  if(length(unique_max_rank) > 1) {
    stop("Requires all max_rank to be equal")
  }
  gamma_thresholds_df <- get_precomputed_gamma_thresholds(K = unique_max_rank + 1)
  results$stats %>%
    group_by(variable) %>%
    mutate(log_gamma = compute_log_gamma_history(rank, unique_max_rank)) %>%
    filter(sim_id >= min_sim_id, sim_id <= max_sim_id) %>%
    inner_join(gamma_thresholds_df, by = c("sim_id" = "N_sims")) %>%
    ggplot(aes(x = sim_id, y = log_gamma - log_gamma_threshold, color = variable)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_line() +
    scale_y_continuous("Log Gamma - Threshold")
}

plot_ks_test_history <- function(results, min_sim_id = 0, max_sim_id = Inf) {
  results$stats %>%
    group_by(variable) %>%
    mutate(ks_p = compute_ks_test_history(rank, unique(max_rank))) %>%
    filter(sim_id >= min_sim_id, sim_id <= max_sim_id) %>%
    ggplot(aes(x = sim_id, y = ks_p, color = variable)) +
    geom_hline(yintercept = 0.05, color = "black") +
    geom_line() +
    scale_y_log10("P - value (KS test)")
}
