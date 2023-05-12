functions {
 //Input: ordered vector
 vector ordered_simplex_constrain_lp(vector y) {
    int Km1 = rows(y);
    vector[Km1 + 1] x;
    real stick_len = 1.0;

    // use only negative y's
    // start from the smallest -y, which is the last value in the vector
    for (k in 1:Km1) {
      real adj_y_k = -y[Km1 - k + 1] - log(Km1 - k + 1);
      real z_k = inv_logit(adj_y_k);
      x[k] = stick_len * z_k;
      target += log(stick_len) - log1p_exp(-adj_y_k) - log1p_exp(adj_y_k);
      stick_len -= x[k];
    }
    // this is new,
    // instead of adding the stick_len to the last element
    // distribute it to all the K - 1 values
    // could also just distribute evenly
    // comment out the x* = 1 + stick_len and just
    // return x + stick_len / Km1;
    x[Km1 + 1] = 0;
    //x *= 1 + stick_len;
    //return x + (1 - sum(x) ) / Km1;
    return x + stick_len / (Km1 + 1);
  }
}
data {
  int K;
  array[K] int<lower=0> observed;
  vector<lower=0>[K] prior_alpha;
}


parameters {
  positive_ordered[K - 1] y;
}

transformed parameters {
  simplex[K] x = ordered_simplex_constrain_lp(y);
}

model {
  x ~ dirichlet(prior_alpha);
  observed ~ multinomial(x);
}
