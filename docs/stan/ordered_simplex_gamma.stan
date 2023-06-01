data {
  int K;
  array[K] int<lower=0> observed;
  vector<lower=0>[K] prior_alpha;
}

parameters {
  positive_ordered[K] w;
}

transformed parameters {
  simplex[K] x = w / sum(w);
}

model {
  w ~ gamma(prior_alpha, 1);
  observed ~ multinomial(x);
}
