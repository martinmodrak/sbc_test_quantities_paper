functions {
  vector ordered_simplex_constrain_softmax_lp(vector v) {
     int K = size(v) + 1;
     vector[K] v0 = append_row(0, v);
     // Jacobian
     target += sum(v) - K * log_sum_exp(v0);
     return softmax(v0);
  }
}

data {
  int K;
  array[K] int<lower=0> observed;
  vector<lower=0>[K] prior_alpha;
}


parameters {
  positive_ordered[K - 1] v;
}

transformed parameters {
  simplex[K] x =  ordered_simplex_constrain_softmax_lp(v);
}

model {
  x ~ dirichlet(prior_alpha);
  observed ~ multinomial(x);
}
