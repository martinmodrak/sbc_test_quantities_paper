data {
  int K;
  array[K] int<lower=0> observed;
  vector<lower=0>[K] prior_alpha;
}


parameters {
  positive_ordered[K - 1] y;
}

transformed parameters {
  simplex[K] x = (1 + append_row(0, y)) / (K + sum(y));
}

model {
  x ~ dirichlet(prior_alpha);
  // Jacobian
  target += - K * log(sum(y));
  observed ~ multinomial(x);
}
