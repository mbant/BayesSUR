#' BayesSUR
#' @title uniform random numbers
#' @description
#' Subfunction for 
#' 
#' @param x response matrix
#' @export
randVec <- function(x){
   
     rb <- randBeta(1,1)
     re <- randExponential(1)
     rbin <- randBinomial(10,.5)
     rmulti <- randMultinomial(10, c(.4,.6))
     rnor <- randNormal(1., 2.)
     rndt <- randT(1.)
     rwish <- randWishart(3., matrix(c(3.,2.,1.,2.,4.,1.,1.,1.,5.),nrow=3,ncol=3))
     rber <- randBernoulli(.4)
     ru01 <- randU01()
     rlogu01 <- randLogU01()
     rintunif <- randIntUniform(1, 5)
     rg <- randGamma(.5,.5)
     rig <- randIGamma(.5,.5)
     rvecexp <- randVecExponential(5, 1.5)
     rvecnor <- randVecNormal(3)
     rvect <- randVecT(4, 1.3)
     rmt <- randMvT(1., c(.6,.8), matrix(c(3.,1.,1.,3.),nrow=2,ncol=2))
     rmvnor <- randMvNormal(c(.6,.8), matrix(c(3.,1.,1.,3.),nrow=2,ncol=2))
     riwish <- randIWishart(3., matrix(c(3.,2.,1.,2.,4.,1.,1.,1.,5.),nrow=3,ncol=3))
     rmn <- randMN(matrix(c(3.,1.,.5,2.),nrow=2,ncol=2), matrix(c(3.,1.,1.,3.),nrow=2,ncol=2), matrix(c(3.,1.,1.,3.),nrow=2,ncol=2))
     rtruncnorm <- randTruncNorm(.5, 1, .1, .9)
     rvecsamp <- randVecSampleWithoutReplacement(5, c(1,0,1,1,0), 5)
     rsamp <- randSampleWithoutReplacement(10, c(1,2), 2)
     rvecwsamp <- randVecWeightedSampleWithoutReplacement(2, c(.4,.6), 2, c(0,1))
     rwsamp <- randWeightedSampleWithoutReplacement(2, c(.4,.6), c(0,1))
     rwisamp <- randWeightedIndexSampleWithoutReplacement(2, c(.4,.6), 2)
     risamp <- randIndexSampleWithoutReplacement(10, 2)
     ruwisamp <- randWeiIndexSampleWithoutReplacement(10, c(.4,.6))
     
  return(list(rb, re, rbin, rmulti, rnor, rndt, rwish, rber, ru01, rlogu01, rintunif, rg, rig, rvecexp, rvecnor, rvect,rmt,
              rmvnor, riwish, rmn, rtruncnorm, rvecsamp, rsamp, rvecwsamp, rwsamp, rwisamp, risamp, ruwisamp))
}
