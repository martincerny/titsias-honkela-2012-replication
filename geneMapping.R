loadMapping <- function() {
  mapping = read.table("array_data/CG_AffyOligo_Gadfly3_01_13_03", sep="\t",col.names = c("Gadfly", "Symbol", "FBgn", "Oligo"), stringsAsFactors = FALSE, quote = "")
  mapping$Oligo = sapply(mapping$Oligo, function(oligo) {unique(lapply(strsplit(oligo,";"), function(x) { sub("_[0-9]*$", "", x) })[[1]])}, USE.NAMES=FALSE)
  return(mapping);
}

spotsFromFBgn <- function(mapping, flybase_id)
{
  matchingIndices = mapping$FBgn == flybase_id;
  if(any(matchingIndices)) {
    return(mapping[matchingIndices,"Oligo"][[1]]);
  }
  else
  {
    return(character(0));
  }
}

spotsFromSymbol <- function(mapping, symbol)
{
  matchingIndices = mapping$Symbol == symbol;
  if(any(matchingIndices)) {
    return(mapping[matchingIndices,"Oligo"][[1]]);
  }
  else
  {
    return(character(0));
  }
}


rowIDsFromSpotNames <- function(normalizedData, spotNames) {
  return(match(spotNames, normalizedData$genes));
}