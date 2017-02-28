normalizeGenes <- function(data) {
  #based on .getProcessedData() function from tigre package
  times <- data$modeltime
  experiments <- data$experiments
  
  y <- assayDataElement(data, 'exprs')
  if ('var.exprs' %in% assayDataElementNames(data))
    yvar <- assayDataElement(data, 'var.exprs')
  else
    yvar <- 0 * y
  
  scale <- sqrt(rowMeans(y^2))
  scale[scale==0] <- 1
  scaleMat <- scale %*% array(1, dim = c(1, ncol(data)))
  y <- y / scaleMat
  yvar <- yvar / scaleMat^2
  
  
  expids <- unique(experiments)
  times <- times[experiments==expids[1]]
  
  genes <- featureNames(data)
  
  newData <- list(y = y, yvar = yvar, genes = genes, times = times, experiments = expids)
  return (newData)
}


