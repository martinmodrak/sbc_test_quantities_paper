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
    if(max_sims > 1000) {
      sims_above_1000 <- round(exp(seq(log(1001), log(max_sims), by = 0.02)))
      #Ensure max_sims is present
      sims_above_1000 <- unique(c(sims_above_1000, max_sims))
      if(min_sims <= 1000) {
        sims_to_compute <- c(min_sims:1000, sims_above_1000)
      } else {
        sims_to_compute <- sims_above_1000
      }
      interpolate <- TRUE
    } else {
      sims_to_compute <- min_sims:max_sims
      interpolate <- FALSE
    }
    # Computing time differs by N, randomly permute to get optimal division of work
    sims_to_compute <- sample(sims_to_compute)
    thresholds <- future.apply::future_sapply(sims_to_compute, FUN = function(n_sims) {
      SBC:::adjust_gamma(N = n_sims, L = 1, K = K)
    })


    if(interpolate) {
      thres_approx <- approx(log(sims_to_compute), log(thresholds), xout = log(min_sims:max_sims))
      thresholds_df <- data.frame(N_sims = min_sims:max_sims, log_gamma_threshold = thres_approx$y)
    } else {
      thresholds_df <-
        dplyr::arrange(data.frame(N_sims = sims_to_compute, log_gamma_threshold = log(thresholds)), N_sims)
    }

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


compute_ks_test_history_dgof <- function(ranks, max_rank) {
  ks_p <- rep(NA_real_, length(ranks))
  reference <- ecdf(0:max_rank)
  for(i in 1:length(ranks)) {
    ks_p[i] <- dgof::ks.test(ranks[1:i], reference)$p.value
  }
  ks_p
}


plot_log_gamma_history <- function(results, min_sim_id = 0, max_sim_id = Inf, wrap_cols = 4, variables_regex = NULL, ylim = NULL) {
  unique_max_rank <- unique(results$stats$max_rank)
  if(length(unique_max_rank) > 1) {
    stop("Requires all max_rank to be equal")
  }

  max_sim_id_to_show <- min(max_sim_id, max(results$stats$sim_id))

  gamma_thresholds_df <- get_precomputed_gamma_thresholds(K = unique_max_rank + 1,
                                                          min_sims = 1, max_sims = max(max_sim_id_to_show, 1000))

  stats <- results$stats
  if(!is.null(variables_regex)) {
     stats <- stats %>% filter(grepl(variables_regex, variable))
  }

  stats %>%
    filter(sim_id <= max_sim_id) %>%
    group_by(variable) %>%
    mutate(log_gamma = compute_log_gamma_history(rank, unique_max_rank)) %>%
    filter(sim_id >= min_sim_id) %>%
    inner_join(gamma_thresholds_df, by = c("sim_id" = "N_sims")) %>%
    ggplot(aes(x = sim_id, y = log_gamma - log_gamma_threshold)) +
    geom_hline(yintercept = 0, color = "lightblue") +
    geom_line() +
    scale_y_continuous("Log Gamma - Threshold", limits = ylim) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}

plot_ks_test_history <- function(results, min_sim_id = 0, max_sim_id = Inf, wrap_cols = 4, ylim = NULL) {
  results$stats %>%
    filter(sim_id <= max_sim_id) %>%
    group_by(variable) %>%
    mutate(ks_p = compute_ks_test_history(rank, unique(max_rank))) %>%
    filter(sim_id >= min_sim_id) %>%
    ggplot(aes(x = sim_id, y = ks_p)) +
    geom_hline(yintercept = 0.05,  color = "lightblue") +
    geom_line() +
    scale_y_log10("P - value (KS test)", limits = ylim) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}

plot_ks_test_history_dgof <- function(results, min_sim_id = 0, max_sim_id = Inf, wrap_cols = 4, ylim = NULL) {
  results$stats %>%
    filter(sim_id <= max_sim_id) %>%
    group_by(variable) %>%
    mutate(ks_p = compute_ks_test_history_dgof(rank, unique(max_rank))) %>%
    filter(sim_id >= min_sim_id) %>%
    ggplot(aes(x = sim_id, y = ks_p)) +
    geom_hline(yintercept = 0.05,  color = "lightblue") +
    geom_line() +
    scale_y_log10("P - value (KS test)", limits = ylim) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}

