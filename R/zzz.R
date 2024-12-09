#Run 'tidymodels_prefer()' when the package is loaded

.onLoad <- function(libname, pkgname) {
  # Safely call tidymodels_prefer() to avoid conflicts
  if (requireNamespace("tidymodels", quietly = TRUE)) {
    tidymodels::tidymodels_prefer()
  } else {
    warning("The 'tidymodels' package is not installed, so tidymodels_prefer() was not called.")
  }
}