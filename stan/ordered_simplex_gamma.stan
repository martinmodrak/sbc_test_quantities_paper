data {
  int K;
  array[K] int<lower=0> observed;
  vector<lower=0>[K] prior_alpha;
}

parameters {
  positive_ordered[K] lambda;
}

transformed parameters {
  simplex[K] x = lambda / sum(lambda);
}

model {
  lambda ~ gamma(prior_alpha, 1);
  observed ~ multinomial(x);
}
