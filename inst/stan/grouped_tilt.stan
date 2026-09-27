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

  // log_sum_exp(rate * lo, rate * (lo + 1), ..., rate * hi): the log of the
  // total unnormalised mass on one rank interval. The arguments are an
  // arithmetic sequence, so this is a geometric series and the rank loop
  // collapses. With t = hi - lo + 1 and rate > 0,
  //
  //   log sum_r exp(rate * r)
  //     = rate * lo + log(exp(rate * t) - 1) - log(exp(rate) - 1).
  //
  // theta is one parameter; the sign split is not a second model. log(e^x - 1)
  // is only real for x > 0, and rate * t is negative when rate is, so the
  // rate < 0 line factors the same sum from the top of the interval:
  //
  //   rate * hi + log(exp(-rate * t) - 1) - log(exp(-rate) - 1).
  //
  // Both lines equal rate * lo + log|e^{rate*t} - 1| - log|e^{rate} - 1|.
  real geom_logsum(int lo, int hi, real rate) {
    int t = hi - lo + 1;
    if (abs(rate) < 1e-10) {
      return log(t) + rate * (lo + hi) / 2.0;
    } else if (rate > 0) {
      return rate * lo + geom_log_expm1(rate * t) - geom_log_expm1(rate);
    } 
    // The flipped sign only keeps the argument of geom_log_expm1 positive.
    // Near rate = 0 those two logs cancel, and the sum is
    // log(t) + rate * (lo + hi) / 2. The loop over signature genes is in the
    // model block; each gene adds one call, then the shared normaliser on
    // 1:G is subtracted once per gene.
  else {
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
