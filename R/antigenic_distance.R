generate_antigenic_map <- function(antigenic_distances, buckets = 1, year_min = 1968, year_max = 2016) {
  ## Convert strains to correct time dimensions
  antigenic_distances$strain_year <- antigenic_distances$strain_year * buckets
  ## Fit spline through antigenic coordinates
  fit <- smooth.spline(antigenic_distances$X, antigenic_distances$Y, nknots = 10)
  
  ## Work out relationship between strain circulation time and x coordinate
  x_line <- lm(data = antigenic_distances, X ~ strain_year)
  
  ## Enumerate all strains that could circulate
  strain_year <- seq(year_min * buckets, year_max * buckets - 1, by = 1)
  
  ## Predict x and y coordinates for each possible strain from this spline
  x_predict <- predict(x_line, data.frame(strain_year))
  y_predict <- predict(fit, x = x_predict)
  
  fit_data <- tibble(
    x = y_predict$x, y = y_predict$y,
    strain_year = strain_year
  )
  
  return(fit_data)
}

generate_antigenic_distances <- function(df) {
  distance_matrix <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  
  # Fill the distance matrix with Euclidean distances
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(df)) {
      distance_matrix[i, j] <- sqrt((df$x[i] - df$x[j])^2 + (df$y[i] - df$y[j])^2)
    }
  }
  
  return(distance_matrix)
}
