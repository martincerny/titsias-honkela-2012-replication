getNumDetailedTime <- function(numTime, numIntegrationPoints) {
  return((numTime - 1) * numIntegrationPoints + 1);
}