data {
  int K;
  array[K] int<lower=0> observed;
  vector<lower=0>[K] prior_alpha;
}


parameters {
  positive_ordered[K - 1] y;
}

transformed parameters {
  simplex[K] x = softmax(append_row(0, y));
}

model {
  x ~ dirichlet(prior_alpha);
  // Jacobian
  target += sum(y) - (K - 1) * log_sum_exp(append_row(0, y));
  observed ~ multinomial(x);
}
