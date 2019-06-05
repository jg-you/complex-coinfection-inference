functions {
  /*
   * Computes the value of a Bernstein polynomial of degree N at a < t < b,
   * using De Casteljau's algorithm.
   *
   * @param t      Point in [a, b] where the the polynomial will be evaluated.
   * @param a      Lower bound of the domain of the function.
   * @param b      Upper bound of the domain of the function.
   * @param xi     Vector of the N + 1 coefficients of the Bernstein polynomial.
   * @param N      Degree N of the polynomial.
   *
   * @return Value of the polynomial at point t.
   */
  real bernstein(real t,
                 real a,
                 real b,
                 vector xi,
                 int N)
  {
    real t_prime = t / (b - a);
    vector[N + 1] curr_xi;
    vector[N + 1] next_xi;
    if (t < a) {
      t_prime = a;
      // reject("Outside of bertnstein range (t=", t," too low)");
    }
    if (t > b) {
      t_prime = b;
      // reject("Outside of bertnstein range (t=", t," too high)");
    }
    for (i in 1:N + 1) {
        curr_xi[i] = xi[i];
    }
    for (j in 1:N + 1) {
      for (i in 1:N + 1 - j) {
        next_xi[i] = curr_xi[i] * (1 - t_prime) + curr_xi[i + 1] * t_prime;
      }
      curr_xi = next_xi;
    }
    return next_xi[1];
  }
  /*
   * ODE system for the mean-field complex SIR dynamics on a well-mixed
   * population.
   *
   * Note: We do not actually need to track the R variable, so the system is
   *       effectively two dimensional.
   *
   * @param t           Time at which the equations will be evaluated.
   * @param y           State variables (y[1], y[2]) = (S, I). These
   *                    variables represent densities of susceptible and
   *                    infected individuals in the population.
   * @param all_rates   Rate parameters.
   *                       all_rates[1]:  gamma, the recovery rate.
   *                       all_rates[2:]: xi, N+1 coefficients of the Bernstein
                                          polynomials.
   * @param args_real   Array of size 2. Bounds of the domain of beta.
   * @param args_int    Array of size 1. Contains the degree of the polynomial.
   *
   * @return Values of the ODE system at time t.
   */
  real[] complex_SIR(real t,
                     real[] y,
                     real[] all_rates,
                     real[] args_real,
                     int[] args_int)
  {
    // Declarations
    real dydt[2];
    real beta_rate;
    real gamma_rate;
    vector[args_int[1] + 1] xi;

    // Get recovery rate and transmission rate for current state
    gamma_rate = all_rates[1];
    for (i in 1:args_int[1] + 1)
    { // Place Bernstein coefficients in a vector
      xi[i] = all_rates[i + 1];
    }
    beta_rate =  bernstein(y[2], args_real[1], args_real[2], xi, args_int[1]);

    // ODE system
    dydt[1] = -beta_rate * y[1] * y[2];
    dydt[2] = beta_rate * y[1] * y[2]  - gamma_rate * y[2];
    return dydt;
  }
  // =========================================================================
  // Output functions
  //
  // These functions' purpose is to clean-up the generated quantities block.
  // =========================================================================
  /* Generates the beta function on [a, b], at `num_steps` points.
   *
   * @param all_rates   Rate parameters.
   *                       all_rates[1]:  gamma, the recovery rate.
   *                       all_rates[2:]: xi, N+1 coefficients of the Bernstein
                                          polynomials.
   * @param a      Lower bound of the domain of the function.
   * @param b      Upper bound of the domain of the function.
   * @param num_steps    Resolution of the calculated xi function.
   * @param N      Degree N of the polynomial.
   */
  real[] get_beta_func(real[] all_rates,
                       real a,
                       real b,
                       int num_steps,
                       int N)
  {
    real beta_func[num_steps];
    vector[N + 1] xi;
    real t;
    real delta_t;
    // Place Bernstein coefficients in a vector
    for (i in 1:N + 1)
    {
      xi[i] = all_rates[i + 1];
    }

    delta_t = (b - a) / num_steps;
    for (i in 1:num_steps)
    {
      beta_func[i] =  bernstein(a + i * delta_t, a, b, xi, N);
    }
    return beta_func;
  }
  /* Generate an artificial time series of time series,
   *
   * @param T           Length of the time series.
   * @param y0          Initial conditions for (S0, I0). R0 is implicit.
   * @param all_rates   Rate parameters.
   *                       all_rates[1]:  gamma, the recovery rate.
   *                       all_rates[2:]: xi, N+1 coefficients of the Bernstein
                                          polynomials.
   * @param ts          Time at which the equations will be evaluated.
   * @param args_real   Array of size 2. Bounds of the domain of beta.
   * @param args_int    Array of size 1. Contains the degree of the polynomial.
   * @param population  Size of the population.
   */
  real[] get_z_tilde(int T,
                     real[] y0,
                     real[] all_rates,
                     real[] ts,
                     real[] args_real,
                     int[] args_int,
                     int population)
  {

    real y_tilde[T, 2];
    real z_tilde[T];
    y_tilde = integrate_ode_rk45(complex_SIR, y0, 0, ts, all_rates, args_real, args_int);
    z_tilde[1] = population * (y0[1] - y_tilde[1, 1]);
    for (t in 2:T)
    {
      z_tilde[t] = population * (y_tilde[t - 1, 1] - y_tilde[t, 1]);
    }
    return z_tilde;
  }
}