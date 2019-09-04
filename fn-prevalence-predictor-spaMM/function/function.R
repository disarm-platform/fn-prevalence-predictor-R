suppressPackageStartupMessages(library(spaMM))
suppressPackageStartupMessages(library(sf))
library(geojsonio)
library(httr)

get_posterior_metrics <- dget("function/helpers.R")

function(params) {

  # Read into memory
  point_data <- st_read(as.json(params$point_data), quiet = TRUE)
  
  # Send to covariate_extractor
  cov_ext_input_data_list <- list(points = geojson_list(point_data),
                                  layer_names = params$layer_names)
  
  response <-
    httr::POST(
      url = "https://faas.srv.disarm.io/function/fn-covariate-extractor",
      body = as.json(cov_ext_input_data_list),
      content_type_json(),
      timeout(90)
    )
  
  # Get contents of the response
  response_content <- content(response)
  
  points_sf <- st_read(as.json(response_content$result), quiet = TRUE)
  
  # Annoyingly values come back as factors, so change back to numeric
  points_sf$n_trials <- as.numeric(as.character(points_sf$n_trials))
  points_sf$n_positive <- as.numeric(as.character(points_sf$n_positive))
  #points_sf$id <- 1:nrow(points_sf)
  
  # Now pass into fn-cv-ml
  cv_hal_input_data_list <- list(points = geojson_list(points_sf),
                                 layer_names = params$layer_names)
  
  response_cv_ml <-
    httr::POST(
      url = "https://faas.srv.disarm.io/function/fn-cv-ml",
      body = as.json(cv_hal_input_data_list),
      content_type_json(),
      timeout(90)
    )
  
  response_content_cv_ml <- content(response_cv_ml)
  
  # Clean up data ready to pass into model
  mod_data_sf <- st_read(as.json(response_content_cv_ml$result$points), quiet=TRUE)
  mod_data <- as.data.frame(mod_data_sf)
  mod_data <- cbind(mod_data, st_coordinates(mod_data_sf))
  mod_data$cv_predictions <- as.numeric(as.character(mod_data$cv_predictions))
  mod_data$fitted_predictions <- as.numeric(as.character(mod_data$fitted_predictions))
  mod_data$n_trials <- as.numeric(as.character(mod_data$n_trials))
  mod_data$n_positive <- as.numeric(as.character(mod_data$n_positive))
  mod_data$n_neg <- mod_data$n_trials - mod_data$n_positive
  train_data <- mod_data[!is.na(mod_data$n_trials),]
  pred_data <- mod_data[is.na(mod_data$n_trials),]
  
  # Fit model
  spaMM_mod <- fitme(cbind(n_positive, n_neg) ~
                       cv_predictions +
                       Matern(1|X+Y),
                     data=train_data,
                     family=binomial())
  
  pred_data$cv_predictions <- pred_data$fitted_predictions
  
  # Get posterior metrics
  mod_data$cv_predictions <- mod_data$fitted_predictions
  
  posterior_metrics <- get_posterior_metrics(spaMM_mod,
                        mod_data,
                        200,
                        params$exceedance_threshold)
  
  # Bind to point_data
  for(i in names(posterior_metrics)){
    point_data[[i]] <- posterior_metrics[[i]]
  }
  return(point_data)
}
