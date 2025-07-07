
var_labels <- c(
  "mu_long" = "<i>μ</i><sub>L</sub> (long-term boost)",
  "mu_short" = "<i>μ</i><sub>S</sub> (short-term boost)",
  "sigma_long" = "<i>σ</i><sub>L</sub> (long-term cross<br>reactivity)", 
  "sigma_short" = "<i>σ</i><sub>S</sub> (short-term cross<br>reactivity)",
  "tau" = "<i>τ</i> (seniority)",
  "omega" = "<i>ω</i> (waning)",
  "obs_sd" = "<i>σ</i><sub>obs</sub> (obs sd.)",
  "total_inf" = "Total infections",
  
  "beta" = "<i>β</i> (linear age effect)",
  
  "mixis_lp" = "Mixture-IS LOO score",
  "lp" = "Log-probability"
)


summaries_previous <- tribble(
  ~name, ~variable, ~run_name,  ~median, ~q95_lower, ~q95_upper,
  "kucharski_2018", "mu_long", "hanam_2018", 2.02, 1.96, 2.08,
  "kucharski_2018", "mu_short", "hanam_2018", 2.69, 2.50, 2.88,
  "kucharski_2018", "sigma_long", "hanam_2018", 0.130, 0.128, 0.132,
  "kucharski_2018", "sigma_short", "hanam_2018", 0.031, 0.026, 0.035,
  "kucharski_2018", "obs_sd", "hanam_2018", 1.29, 1.27, 1.31,
  "kucharski_2018", "tau", "hanam_2018", 0.039, 0.035, 0.042, 
  "kucharski_2018", "omega", "hanam_2018", 0.79,0.74, 0.84,
  
  "kucharski_2018", "mu_long", "fluscape_2009_HI", 0.97, 0.85, 1.13,
  "kucharski_2018", "sigma_long", "fluscape_2009_HI", 0.099, 0.084, 0.114,
  "kucharski_2018", "obs_sd", "fluscape_2009_HI", 1.50, 1.41, 1.59,
  "kucharski_2018", "tau", "fluscape_2009_HI", 0.016, 0.012, 0.021,
  
  "kucharski_2018", "mu_long", "fluscape_2009_neuts", 1.38, 1.14, 1.66,
  "kucharski_2018", "sigma_long", "fluscape_2009_neuts", 0.130, 0.106, 0.146,
  "kucharski_2018", "obs_sd", "fluscape_2009_neuts", 1.69, 1.56, 1.82,
  "kucharski_2018", "tau", "fluscape_2009_neuts", 0.020, 0.015, 0.027,
  
  
  
  "hay_2024", "mu_long", "hanam_2018", 1.96, 1.89, 2.03,
  "hay_2024", "mu_short", "hanam_2018", 2.65, 2.43, 2.92,
  "hay_2024", "sigma_long", "hanam_2018", 0.115, 0.112, 0.118,
  "hay_2024", "sigma_short", "hanam_2018", 0.0267, 0.02222, 0.0301,
  "hay_2024", "tau", "hanam_2018", 0.0427, 0.0391, 0.047,
  "hay_2024", "omega", "hanam_2018",  2.01 / 2.65, 1.74 / 2.65, 2.27 / 2.65, # Approximate
  "hay_2024", "obs_sd", "hanam_2018", 1.29, 1.27, 1.31
)

