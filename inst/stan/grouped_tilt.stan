// Discrete exponential tilt with ties as interval-censored ranks.
//
//   P(R = r | theta) = exp(theta * r / G) / Z(theta),   r = 1, ..., G
//
// Gene i observed in [lo[i], hi[i]] contributes the marginal of that block.
// theta ~ normal(0, prior_sd) is deliberately wide, so the mode matches the
// maximum likelihood when the likelihood has curvature, and stays finite when
// the likelihood is flat.
functions {
  real geom_log_expm1(real x) {
    if (x > 18) {
      return x;
    }
    return log(expm1(x));
  }

  // log sum_{r=lo}^{hi} exp(rate * r). Closed form of the geometric series.
  real geom_logsum(int lo, int hi, real rate) {
    int t = hi - lo + 1;
    if (abs(rate) < 1e-10) {
      return log(t) + rate * (lo + hi) / 2.0;
    } else if (rate > 0) {
      return rate * lo + geom_log_expm1(rate * t) - geom_log_expm1(rate);
    } else {
      return rate * hi + geom_log_expm1(-rate * t) - geom_log_expm1(-rate);
    }
  }
}
data {
  int<lower=0> S;
  int<lower=1> G;
  array[S] int<lower=1, upper=G> lo;
  array[S] int<lower=1, upper=G> hi;
  real<lower=0> prior_sd;
  real<lower=0> theta_max;
}
parameters {
  real<lower=-theta_max, upper=theta_max> theta;
}
model {
  real rate = theta / G;
  theta ~ normal(0, prior_sd);
  for (i in 1:S) {
    target += geom_logsum(lo[i], hi[i], rate);
  }
  target += -S * geom_logsum(1, G, rate);
}
