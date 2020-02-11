function(params) {
  # Individual check for each parameter
  if (is.null(params[['point_data']])) {
    stop('Missing `point_data` parameter')
  }
}