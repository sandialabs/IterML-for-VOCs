###################################################Scoring helper functions#########################

rmse = function(x){
  return(sqrt(mean(x^2)))
}

rmsr = function(x, e){
  return(sqrt(mean((x/e)^2)))
}

wmse = function(x, w){
  #w = 1/e or 1/e^2
  x = x[is.finite(w)]
  w = w[is.finite(w)]
  return(sqrt(sum(w*(x)^2)/sum(w)))
}

q2 = function(x, p){
  #x are actual endpoint values, p are the predictions
  return(1 - sum((x - p)^2)/sum((x - mean(x))^2) )
}