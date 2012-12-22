rmixchisq <-
function(n,s,q){
  if(q>0){
    mixprobs <- numeric(q+1)
    for(k in (s:(q+s))){
      mixprobs[k-s+1] <- choose(q,k-s)*2^(-q)
    }
    mix <- n*mixprobs
  }else{
    mix=n
  }
  
  sample <- NULL
  for(k in (s:(q+s))){
    sample <- c(sample, rchisq(mix[k-s+1],df=k))
  }
  
  return(sample) 
}
