#'Power Simulations for McNemar's Test for Sensitivity/Specificity
#'
#'This function takes values for sensitivity and specificity of
#'two tests for comparison and calculates the power, given n,
#'needed to observe significant differences in the performance
#'between the two. Power is returned for both sensitivity and
#'specificity by running McNemar's test on simulated true positives
#'and true negatives respectively. Returns a matrix with simulated
#'power given sample size.
#'
#' @param N vector of candidate sample sizes
#' @param n.sims number of simulations to run for each N
#' @param alpha significance level for McNemar
#' @param seed set yo seed!
#' @param sens1  sensitivity of first test
#' @param sens2  sensitivity of second test
#' @param spec1  specificity of first test
#' @param spec2  specificity of second test
#' @param prev  true positive prevalence
#'
#' @export
#'
mcnemar_power <- function(N = 1000,
                          n.sims = 1000,
                          alpha=.05,
                          seed = 12345,
                          prev = p,
                          sens1 = sn1,
                          spec1 = sp1,
                          sens2 = sn2,
                          spec2 = sp2){

  ## initialize vectors
  mc_sens_power <- vector(length=length(N))
  mc_spec_power <- vector(length=length(N))

  set.seed(seed)

  for(k in 1:length(N)){

    mc_pos <- vector(length=n.sims)
    mc_neg <- vector(length=n.sims)

    for(i in 1:n.sims){

      ## simulate group of true positives based on pred

      truepos <- sum(rbinom( n=N[k], size=1, prob = prev))

      ## now simulate test result in participants that are positive and negative

      ## test 1
      posresult_truepos1 <- rbinom( n=truepos, size=1, prob = sens1)
      posresult_trueneg1 <- rbinom (n= (N[k]-truepos), size=1, prob=(1-spec1))

      ## test 2
      posresult_truepos2 <- rbinom( n=truepos, size=1, prob = sens2)
      posresult_trueneg2 <- rbinom (n= (N[k]-truepos), size=1, prob=(1-spec2))

      ## combine test results
      res_tab <- data.frame(group = c(rep("tp", truepos), rep("tn", N[k]-truepos)),
                            test1 = c(posresult_truepos1, posresult_trueneg1),
                            test2 = c(posresult_truepos2, posresult_trueneg2))

      ## mcnemar test will cause an error if top left or bottom right corner in contingency tables are zeroes
      ## check to make sure table is 2x2 with sum and dim functions

      if(sum(dim(with(res_tab[res_tab$group=="tp",], table(test1, test2)))==2)==2){
        m1 <- mcnemar.test(with(res_tab[res_tab$group=="tp",], table(test1, test2)))
        if(is.na(m1$p.value)){
          mc_pos[i] <- FALSE ## to say 'not signifcant' in cases where tests line up perfectly
        }else{
          mc_pos[i] <- (m1$p.value < alpha)
        }
      }else{
        mc_pos[i] <- FALSE
      }

      if(sum(dim(with(res_tab[res_tab$group=="tn",], table(test1, test2)))==2)==2){
        m2 <- mcnemar.test(with(res_tab[res_tab$group=="tn",], table(test1, test2)))
        if(is.na(m2$p.value)){
          mc_neg[i] <- FALSE ## to say 'not signifcant' in cases where tests line up perfectly
        }else{
          mc_neg[i] <- (m2$p.value < alpha)
        }
      }else{
        mc_neg[i] <- FALSE
      }

    }
    mc_sens_power[k] <- sum(mc_pos)/length(mc_pos)
    mc_spec_power[k] <- sum(mc_neg)/length(mc_neg)

  }

  return(data.frame(Sens_Power = mc_sens_power, Spec_Power = mc_spec_power, n=N))
}

