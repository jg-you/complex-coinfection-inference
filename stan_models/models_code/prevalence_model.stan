#include /all_functions.stan
data {
  /* Time series data */
  int<lower=1> T;                 // Number of observed states
  real ts[T];                     // Times of observation
  real<lower=0, upper=1> Y[T];    // Observed densities
  /* Hyperparameters */
  real<lower=0> scale_sigma;      // Control Cauchy prior on variance
  real<lower=0> scale_xi_mean;    // Control scale of prior on infection rate
  real<lower=0> scale_xi_spread;  // Control scale of prior on infection rate
  real<lower=0> scale_gamma;      // Control scale of prior on recovery rate
  int<lower=0> N;                 // Degree of Bernstein form
  /* Misc. */
  real max_I0;                    // Upper bound on the init condition
  real overshoot;                 // Overshoot on recovery rate domain
  int num_steps_beta;             // Output resolution (beta function)
  int num_steps_y;                // Output resolution (y_tilde)
  int max_iter;                   // Hard cut-off on the RK45 integrator.
}
transformed data {
  /* Declaration */
  real t0;
  real args_real[2];
  int  args_int[1];
  real tf[num_steps_y];
  real max_Y = max(Y);
  /* ODE args */
  t0 = 0;
  args_real[1] = -max_Y * overshoot;
  args_real[2] =  max_Y * (1 + overshoot);
  args_int[1] = N;
  /* Misc. */
  for (t in 1:num_steps_y) {tf[t] = ts[T] * t / num_steps_y;}  // for output
}
parameters {
  real<lower = 0>             gamma_rate;
  real<lower = 0>             mean_xi;               // Mean of Bernstein coeff
  real<lower = -mean_xi>      deviation_xi[N + 1];   // Deviation to the mean
  real<lower=0, upper=max_I0> I0;                    // Init. density of infected
  real<lower=0, upper=1-I0>   S0;                    // Init. density of susceptible
  real<lower=0>               sigma;                 // Variance parameter in rescaled space
}
transformed parameters {
  real<lower=0> all_rates[N + 2];  /* all_rates[1] :  gamma
                                    * all_rates[2:] : xi, the N + 1 coeffs of
                                    *             the Bernstein polynomials.
                                    */
  real<lower=0, upper=1> y0[2];    // RK45 needs init. conditions as an array
  y0[1] = S0;
  y0[2] = I0;
  all_rates[1] = gamma_rate;
  for (i in 1:N+1) {
    // Non-centered cauchy
    all_rates[i + 1] = mean_xi + scale_xi_spread * deviation_xi[i];
  }
}
model {
  real y_tilde[T, 2];
  // Prior
  mean_xi ~ normal(0, scale_xi_mean);
  sigma ~ cauchy(0, scale_sigma);
  gamma_rate ~ normal(0, scale_gamma);
  deviation_xi ~ cauchy(0, 1);
  // Integrate ODEs
  y_tilde = integrate_ode_rk45(complex_SIR, y0, t0, ts, all_rates, args_real, args_int, 1e-10, 1e-10, max_iter);
  // Add noise
  Y ~ normal(y_tilde[:, 2], sigma);
}
generated quantities {
  real beta_func[num_steps_beta];
  real y_tilde[num_steps_y, 2];
  y_tilde = integrate_ode_rk45(complex_SIR, y0, t0, tf, all_rates, args_real, args_int);
  beta_func = get_beta_func(all_rates, args_real[1], args_real[2], num_steps_beta, N);
}
