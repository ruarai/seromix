
# Code taken from serosolver
# Hay, J. A., Minter, A., Ainslie, K. E. C., et al.
# An open source tool to infer epidemiological and
# immunological dynamics from serological data: serosolver. (2020).
# https://github.com/seroanalytics/serosolver


generate_antigenic_map <- function(antigenic_distances, modelled_years) {
  
  fit <- smooth.spline(antigenic_distances$X, antigenic_distances$Y)
  
  ## Predict x and y coordinates for each possible strain from this spline
  # x_predict <- predict(x_line, data.frame(strain_year))
  x_predict <- scalemap(modelled_years, modelled_years)
  y_predict <- predict(fit, x = x_predict)
  
  fit_data <- tibble(
    x = y_predict$x, y = y_predict$y,
    strain_year = modelled_years
  )
  
  return(fit_data)
}

scalemap<-function(xx, inf_years) {
  if(max(inf_years)>2012){stop("need infection range to be inside antigenic map")}
  map.range <- c(1968:2012)
  alen <- c(333.83,370.28); alenA <- (alen[2]-alen[1])/(max(xx)-min(xx)); alenB <- alen[1]-alenA*min(xx)
  s1 <- alenA*xx+alenB-alen[1]; s1*length(inf_years)/length(map.range) +alen[1]
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
