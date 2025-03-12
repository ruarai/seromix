


theta_0 = (mu_long = 3.0, mu_short = 5.0, omega = 0.05, sigma_long = 0.2, sigma_short = 0.2, tau = 0.3)


theta = theta_0


logjoint(
    model,
    theta
)