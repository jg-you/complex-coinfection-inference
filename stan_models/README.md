# STAN models


`prevalence_model.stan` implements a model that infers β(I) from a time series of the fraction of infected individuals in a population.
The model will work out of the box, but it will provide better results if the scale of the prior on the noise parameter, the recovery rate and the infection rate function are tuned to the data.
As such, the prevalence model should be preferred whenever this data is available.

`incidence_model.stan` implements a model that infers β(I) from a time series of the number of newly reported cases.
Because prevalence time series contain little information on the recovery rate, the model is likely to fail if the prior on this rate is poorly chosen.
Furthermore, one must specify the population size and the expected prevalence range **as parameters**, since the model explains the incidence time-series with an latent time-series of prevalences.
The prevalence can in theory span the full range [0, 1], even though it never does in practice.
If upon inspecting the posterior distribution, one notices that the latent prevalence times-series never exit some range [a, b] where 0 <= a < b <= 1, then one can tighten the range of allowed prevalences to reduce the number of rejected transitions and improve sampling.
