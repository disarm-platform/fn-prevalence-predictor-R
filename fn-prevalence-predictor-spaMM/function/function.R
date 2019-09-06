suppressPackageStartupMessages(library(spaMM))
suppressPackageStartupMessages(library(sf))
library(geojsonio)
library(httr)

get_posterior_metrics <- dget("function/helpers.R")


choose_batch <- function(XY, entropy, rho, nu, batch_size) {
  # XY : candidate locations where the batch is choosen from
  # entropy: entropy computed at the locations XY
  # rho, nu: Matern covariance parameters of the model
  # batch_size: number of observations in the batch

  XY_ix <- seq(entropy) # Index of potential locations
  ix_max <- which.max(entropy)  # Location with highest entropy
  XY_ix_lo <- XY_ix[-c(ix_max)] # Locations after removing max entropy
  ix_batch <- c(ix_max) # Batch of locations to survey next
  step <- 2
  while (length(ix_batch) < batch_size) {
    tradeoff <- c()
    for (ix_ in XY_ix_lo) {
      new_batch <- c(ix_batch, ix_)
      new_XY <- XY[new_batch, ]
      Kxx <- spaMM::MaternCorr(d = as.matrix(proxy::dist(new_XY)), rho = rho, nu = nu)
      hlogD <- .5 * determinant(Kxx, logarithm = TRUE)$modulus
      tradeoff <- c(tradeoff, entropy[ix_] + sqrt(log(step)) * hlogD)
    }
    relative_max <- which.max(tradeoff)  # Location with highest entropy
    ix_batch <- c(ix_batch, XY_ix_lo[relative_max])
    XY_ix_lo <- XY_ix_lo[-c(relative_max)]
    step <- step + 1
  }

  # If we want as output the indices of locations in XY return this
  #return(ix_batch)

  # If we want as output the locations themselves, return this one
  return(XY[ix_batch, ])
}


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
  
  # If there are >250 points, sample 200 of them to estimate 
  # covariance parameters and then update model with fixed cov pars 
  if(nrow(train_data)>=250){
    set.seed(1981)
    train_data_samp <- train_data[sample(1:nrow(train_data), 200),]
    
    spaMM_mod_samp <- fitme(cbind(n_positive, n_neg) ~
                         cv_predictions +
                         Matern(1|X+Y),
                       data=train_data_samp,
                       family=binomial())
    
    spaMM_mod <- fitme(cbind(n_positive, n_neg) ~
                         cv_predictions +
                         Matern(1|X+Y),
                       data=train_data,
                       fixed = list(nu = spaMM_mod_samp$CorrEst_and_RanFix$corrPars$'1'$nu,
                                    rho = spaMM_mod_samp$CorrEst_and_RanFix$corrPars$'1'$rho),
                       family=binomial())
    
  }else{
  
  # Fit model
  spaMM_mod <- fitme(cbind(n_positive, n_neg) ~
                       cv_predictions +
                       Matern(1|X+Y),
                     data=train_data,
                     family=binomial())
  }
  
  pred_data$cv_predictions <- pred_data$fitted_predictions
  
  # Get posterior metrics
  mod_data$cv_predictions <- mod_data$fitted_predictions
  
  posterior_metrics <- get_posterior_metrics(spaMM_mod,
                        mod_data,
                        200,
                        params$exceedance_threshold)
  
  ### Wrapper for adaptive sampling ###
  batch_size_param <- 10
  new_batch <- choose_batch(XY = train_data[, c('X', 'Y')],
                            entropy = posterior_metrics$entropy,
                            rho = spaMM_mod$corrPars[[1]]$rho,
                            nu = spaMM_mod$corrPars[[1]]$nu,
                            batch_size = batch_size_param)
  ### That's all :) ###
  
  # Bind to point_data
  for(i in names(posterior_metrics)){
    point_data[[i]] <- posterior_metrics[[i]]
  }
  return(geojson_list(point_data))
}
