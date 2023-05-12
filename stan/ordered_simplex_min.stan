functions {
 //Input: vector of numbers constrained to [0,1]
 vector ordered_simplex_constrain_min_lp(vector y) {
    int Km1 = rows(y);
    vector[Km1 + 1] x;
    real remaining = 1; // Remaining amount to be distributed
    real base = 0; // The minimum for the next element
    for(i in 1:Km1) {
      if(y[i] <= 0 || y[i] >= 1) {
        reject("All elements of y have to be in [0,1]");
      }
      int K_prime = Km1 + 2 - i; // Number of remaining elements
      //First constrain to [0; 1 / K_prime]
      real x_cons = inv(K_prime) * y[i];
      // Jacobian for the constraint
      target += -log(K_prime);

      // Add the lowest element log density
      target += log(K_prime - 1) +  log(K_prime) + (K_prime - 2) * log1m(K_prime*x_cons);

      x[i] = base + remaining * x_cons;
      base = x[i];
      //We added  remaining * x_cons to each of the K_prime elements yet to be processed
      remaining -= remaining * x_cons * K_prime;
    }
    x[Km1 + 1] = base + remaining;

    return x;
 }
}
data {
  int K;
  array[K] int<lower=0> observed;
  vector<lower=0>[K] prior_alpha;
}


parameters {
  vector<lower=0, upper=1>[K - 1] y;
}

transformed parameters {
  simplex[K] x = ordered_simplex_constrain_min_lp(y);
}

model {
  x ~ dirichlet(prior_alpha);
  observed ~ multinomial(x);
}
