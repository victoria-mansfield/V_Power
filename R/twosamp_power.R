#'Power Simulations for Two-Sample Comparisons Tests (T-test or Wilcox Rank-Sum)
#'
#'This function takes inputs for two distributions of data (as vectors) and runs either
#'a parametric two-sample T-Test or non-parametric Wilcox rank-sum test. Each simulation
#'samples from the candidate distributions with replacement for each iteration. Returns
#'data frame of simulated power for each candidate N. Will add functionality for differing 
#'sample sizes.
#'
#'@param N Vector representing sample size per group
#'@param n.sims Number of simulations to run for each N
#'@param alpha Significance level for Wilcox Test
#'@param seed Set your seed
#'@param par 'T' for parametric T-test, 'F' for non-parametric Wilcox Rank-Sum test
#'@param dist1 Vector containing values (either real or simulated) for first distribution being compared
#'@param dist2 Vector containing values (either real or simulated) for second distribution being compared
#'
#'@export


twosamp_power <- function(N = 1000,
                         n.sims = 1000,
                         alpha=.05,
                         seed = 12345,
                         par = F,
                         dist1 = a,
                         dist2 = b){

  ## initialize vector
  wil_power <- rep(0, length(N))
  ## set seed
  set.seed(seed)

  for(k in 1:length(N)){

    wil_sig <- vector(length=n.sims)

    for(i in 1:n.sims){

      ## simulate distributions based on inputs
      ## sample with replacement from dist1 and dist2
      dist1_sim <- sample(dist1, N[k], replace = T)
      dist2_sim <- sample(dist2, N[k], replace = T)

      if(par){
        wil_sig[i] <- (t.test(dist1_sim, dist2_sim)$p.value < alpha)
      }else{
        wil_sig[i] <- (wilcox.test(dist1_sim, dist2_sim, exact = F)$p.value < alpha)
      }
    }
    wil_power[k] <- sum(wil_sig)/length(wil_sig)
  }
  return(data.frame(TwoSamp_Power = wil_power, n=N))
}
