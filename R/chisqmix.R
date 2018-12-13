#'Chi-Squared Mixtures Distribution
#'
#'Density, distribution function, quantile function and random generation
#'for mixtures of chi-squared distributions that corresponds to the null 
#'distribution of the Likelihood Ratio between 2 nested mixed models.
#'
#'The approximate null distribution of a likelihood ratio for 2 nested mixed
#'models, where both fixed and random effects are tested simultaneously, is a
#'very specific mixture of \eqn{\chi^2}{chi-square} distributions [\cite{Self & Liang
#'(1987), Stram & Lee (1994) and Stram & Lee (1995)}].  It depends on both the
#'number of random effects and the number of fixed effects to be tested
#'simultaneously: 
#'\deqn{LRT_{H_0}\sim\sum_{k=q}^{q+r}{{r}\choose{k-q}}2^{-r}\chi^2_{(k)}}{LRT_H0~\sum k=q..q+r combination(r,k-q) 2^(-r) \chi^2 (k)}
#'
#'@param n number of observations.
#'@param p a probability.
#'@param x,quant a quantile.
#'@param s number of fixed effects to be tested.
#'@param q number of random effects to be tested.
#'@param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#'
#'@return A vector of random independent observations of the \eqn{\chi^2}{chi-square} mixture
#'identified by the values of \code{s} and \code{q}.
#'
#'@aliases chisqmix rchisqmix dchisqmix qchisqmix pchisqmix
#'@rdname chisqmix
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link{pval_simu}}
#'
#'@references Self, S. G. and Liang, K., 1987, Asymptotic properties of maximum
#'likelihood estimators and likelihood ratio tests under nonstandard
#'conditions, \emph{Journal of the American Statistical Association} 82:
#'605--610.
#'
#'Stram, D. O. and Lee, J. W., 1994, Variance components testing in the
#'longitudinal mixed effects model, \emph{Biometrics} 50: 1171--1177.
#'
#'Stram, D. O. and Lee, J. W., 1995, Corrections to "Variance components
#'testing in the longitudinal mixed effects model" by Stram, D. O. and Lee, J.
#'W.; 50: 1171--1177 (1994), \emph{Biometrics} 51: 1196.
#'
#'
#'@importFrom stats rchisq
#'
#'@export
#'
#'@examples
#'library(graphics)
#'library(stats)
#'
#'sample_mixt <- rchisqmix(n=1000, s=3, q=3)
#'plot(density(sample_mixt))
#'
#'
rchisqmix <- function(n,s,q){
  if(q>0){
    mixprobs <- numeric(q+1)
    for(k in (s:(q+s))){
      mixprobs[k-s+1] <- choose(q,k-s)*2^(-q)
    }
    mix <- n*mixprobs
    #s <- s + q*(q-1)/2 #conservative corrections for the covariances terms
  }else{
    mix=n
  }
  
  sample <- numeric(sum(mix))
  mix0 <- cumsum(c(0, mix))
  for(k in (s:(q+s))){
    sample[(mix0[k-s+1] + 1):mix0[k-s+2]] <- stats::rchisq(mix[k-s+1],df=k)
  }
  
  return(sample) 
}

#'@rdname chisqmix
#'@export
dchisqmix <- function(x, s, q){
	if(q>0){
		mixprobs <- numeric(q+1)
		for(k in (s:(q+s))){
			mixprobs[k-s+1] <- choose(q,k-s)*2^(-q)
		}
		mix <- mixprobs
		#s <- s + q*(q-1)/2 #conservative corrections for the covariances terms
	}else{
		mix=1
	}
	
	res <- numeric(length(mix))
	for(k in (s:(q+s))){
		res[k-s+1] <- mix[k-s+1]*stats::dchisq(x,df=k)
	}
	return(sum(res))
}

#'@rdname chisqmix
#'@export
qchisqmix <- function(p, s, q){
	if(q>0){
		mixprobs <- numeric(q+1)
		for(k in (s:(q+s))){
			mixprobs[k-s+1] <- choose(q,k-s)*2^(-q)
		}
		mix <- mixprobs
		#s <- s + q*(q-1)/2 #conservative corrections for the covariances terms
	}else{
		mix=1
	}
	
	res <- numeric(length(mix))
	for(k in (s:(q+s))){
		res[k-s+1] <- mix[k-s+1]*stats::qchisq(p, df=k)
	}
	return(sum(res)/sum(mix))
}

#'@rdname chisqmix
#'@export
pchisqmix <- function(quant, s, q, lower.tail=TRUE){
	if(q>0){
		mixprobs <- numeric(q+1)
		for(k in s:(q+s)){
			mixprobs[k-s+1] <- choose(q, k-s)*2^(-q)
		}
		mix <- mixprobs
		#s <- s + q*(q-1)/2 #conservative corrections for the covariances terms
	}else{
		mix=1
	}
	
	res <- numeric(length(mix))
	for(k in (s:(q+s))){
		res[k-s+1] <- mix[k-s+1]*stats::pchisq(quant, df=k, lower.tail = lower.tail)
	}
	return(sum(res)/sum(mix))
}