data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0.1> a;
  real<lower=0.1> b;
  real c;
  real intercept;
  real<lower=0> prec;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  a ~ uniform(0.1, 5);
  b ~ uniform(0.1, 4);
  c ~ normal(0, sqrt(1000));
  intercept ~ normal(0, sqrt(1000));
  prec ~ gamma(1, pow(1, -5));
  y ~ normal(intercept + c*log(1 + pow((x/a), b)), 1/sqrt(prec));
}

