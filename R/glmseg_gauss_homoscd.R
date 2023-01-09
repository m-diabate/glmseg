#' simuData_gauss_homoscd
#'
#' @param X a n x p matrix of covariates. The first column of X is filled with 1 and always corresponds to the intercept. Thus, if p = 1, then there are no covariates
#' @param par a vector of p*K+1 parameters: p parameters for each of the K sub-segments + 1 parameter corresponding to the standard deviation common to all sub-segments
#' @param bp a vector of K-1 fixed breakpoints
#'
#' @return a vector of n (homoscedastic) gaussian observations
#' @export
#' @importFrom stats rnorm
#' @examples
#'n = 1000 ; X = matrix(1,n,1)
#'bp = c(553)
#'par = c(10,12,3)
#'y <- simuData_gauss_homoscd(X,par,bp)
#'y
simuData_gauss_homoscd=function(X,par,bp=NULL) {
  n=nrow(X)
  p=ncol(X)
  K=length(bp)+1
  if (length(par)!=(p*K+1)) stop("Error !") # par[(p*K)+1] = sigma (constant pour tous les segments)
  beta=matrix(par[1:(p*K)],nrow=p)
  Xbeta=X%*%beta
  bp0=c(0,bp,n)
  y=rep(NA,n)
  for (k in 1:K) {
    idx=(bp0[k]+1):(bp0[k+1])
    y[idx] <- rnorm(n=idx, mean=Xbeta[idx,k], sd = par[p*K+1])
  }
  return(y)
}

#' estim_gauss_homoscd
#'
#' @param y a vector of n observations (supposed gaussian and homoscedastic)
#' @param X a n x p matrix of covariates. The first column of X is filled with 1 and always corresponds to the intercept. Thus, if p = 1, then there are no covariates
#' @param bp a vector of K-1 fixed breakpoints
#'
#' @return a list containing: a vector beta of estimated regression parameters for each sub-segment, the estimated standard deviation sigma common to all sub-segments, and the estimated logLikelihood loglik corresponding to the sum of the log-likelihoods on all sub-segments
#' @export
#' @importFrom stats coefficients lm logLik sigma
#' @examples
#'n = 1000 ; X = matrix(1,n,1)
#'bp = c(553)
#'par = c(10,12,3)
#'y <- simuData_gauss_homoscd(X,par,bp)
#'res_estimated <- estim_gauss_homoscd(y,X,bp)
#'res_estimated
estim_gauss_homoscd=function(y,X,bp=NULL) {
  n=nrow(X)
  K=length(bp)+1
  bp0=c(0,bp,n)
  fit=vector(K,mode="list")
  for (k in 1:K) {
    idx=(bp0[k]+1):(bp0[k+1])
    fit[[k]]=glm(y~.-1,family=gaussian,data=data.frame(y=y,X=X)[idx,])
  }
  return(list(beta=lapply(fit,coefficients),sigma=lapply(fit,sigma),loglik=sum(sapply(fit,logLik))))
}


#' ComputeCritL_gauss_homoscd
#'
#' @param y a vector of n observations (supposed gaussian and homoscedastic)
#' @param X a n x p matrix of covariates. The first column of X is filled with 1 and always corresponds to the intercept. Thus, if p = 1, then there are no covariates
#' @param K a parameter specifying that we want K-1 breakpoints
#' @param bp_ProbRange a range of values including a given breakpoint
#'
#' @return a list containing respectively: the log-likelihood values and the BIC values for each of the points in bp_ProbRange
#' @export
#' @examples
#' n = 1000 ; X = matrix(1,n,1)
#' bp = c(553)
#' par = c(10,12,3)
#' y <- simuData_gauss_homoscd(X,par,bp)
#' K <- 2
#' largeRange <- round(unname(quantile(1:n, c(0.05,0.95))))
#' bp_ProbRange = largeRange[1]:largeRange[2]
#' ComputeCritL_gauss_homoscd(y,X,K,bp_ProbRange)
ComputeCritL_gauss_homoscd <- function(y,X,K,bp_ProbRange){
  critL=rep(NA,length(bp_ProbRange))
  # K = length(bp)+1
  # #K = dim(lpi_jk)[1]
  #---- critBIC : pour le critere BIC (pour plus tard) -----
  critBIC=rep(NA,length(bp_ProbRange))
  n=nrow(X)
  p=ncol(X)
  for (i in 1:length(bp_ProbRange)) {
    critL[i]=estim_gauss_homoscd(y,X,bp=bp_ProbRange[i])$loglik
    critBIC[i]=-2*critL[i] + p*K*log(n)
  }
  return(list(critL = critL, critBIC = critBIC))
}

#' Title
#'
#' @param y a vector of n observations (supposed gaussian and homoscedastic)
#' @param X a n x p matrix of covariates. The first column of X is filled with 1 and always corresponds to the intercept. Thus, if p = 1, then there are no covariates
#' @param K a parameter specifying that we want K-1 breakpoints
#'
#' @return a list containing: the regression results and the computed breakpoint maximizing the log-likelihood
#' @export
#' @importFrom stats quantile
#' @examples
#' n = 1000 ; X = matrix(1,n,1)
#' bp = c(553)
#' par = c(10,12,3)
#' y <- simuData_gauss_homoscd(X,par,bp)
#' K <- 2
#' res_bruteForce_gauss_homoscd =  bruteForce_gauss_homoscd(y,X,K)
#' res_bruteForce_gauss_homoscd
#' res_bruteForce_gauss_homoscd$bpForceBrute
bruteForce_gauss_homoscd <- function(y,X,K){ # bruteForce qui marche que pour un bp
  n = length(y)
  largeRange <- round(unname(quantile(1:n, c(0.05,0.95))))
  bp_ProbRange = largeRange[1]:largeRange[2]

  resCrit <- ComputeCritL_gauss_homoscd(y,X,K,bp_ProbRange)
  critL = resCrit$critL
  fit2=estim_gauss_homoscd(y,X,bp_ProbRange[which.max(critL)])
  bpForceBrute = c(bp_ProbRange[which.max(critL)])
  return(list(resForceBrute = fit2, bpForceBrute = bpForceBrute))
}


#' max_crit_gauss_homoscd
#'
#' @param y a vector of n observations (supposed gaussian and homoscedastic)
#' @param X a n x p matrix of covariates. The first column of X is filled with 1 and always corresponds to the intercept. Thus, if p = 1, then there are no covariates
#' @param par a vector of p*K+1 parameters: p parameters for each of the K sub-segments + 1 parameter corresponding to the standard deviation common to all sub-segments
#' @param lpi_jk a KxK matrix specifying that we want K-1 breakpoints (to be simplified !!!)
#'
#' @return log-likelihood
#' @export
#' @importFrom stats dnorm
#' @examples
#' n = 1000 ; X = matrix(1,n,1)
#' bp = c(553)
#' par = c(10,12,3)
#' y <- simuData_gauss_homoscd(X,par,bp)
#' pi_jk = matrix(c(1,0,1,1),nrow = 2) ; lpi_jk = log(pi_jk)
#' max_crit_gauss_homoscd(y,X,par,lpi_jk)
max_crit_gauss_homoscd=function(y,X,par,lpi_jk){
  p = ncol(X)
  n = nrow(X)
  K = dim(lpi_jk)[1]
  beta=matrix(par[1:(K*p)],nrow=p)
  Xbeta=X%*%beta
  le = matrix(dnorm(y, mean = Xbeta, sd = par[K*p+1],log = TRUE),n,K)
  res = Fw_factLog_Niv2(exp(lpi_jk), le, max)
  return(res)
}

#' map_crit_gauss_homoscd
#'
#' @param y a vector of n observations (supposed gaussian and homoscedastic)
#' @param X a n x p matrix of covariates. The first column of X is filled with 1 and always corresponds to the intercept. Thus, if p = 1, then there are no covariates
#' @param par a vector of p*K+1 parameters: p parameters for each of the K sub-segments + 1 parameter corresponding to the standard deviation common to all sub-segments
#' @param lpi_jk a KxK matrix specifying that we want K-1 breakpoints (to be simplified !!!)
#'
#' @return log-likelihood, segmentation and corresponding breakpoint(s)
#' @export
#' @importFrom stats dnorm
#' @examples
#' n = 1000 ; X = matrix(1,n,1)
#' bp = c(553)
#' par = c(10,12,3)
#' y <- simuData_gauss_homoscd(X,par,bp)
#' pi_jk = matrix(c(1,0,1,1),nrow = 2) ; lpi_jk = log(pi_jk)
#' map_crit_gauss_homoscd(y,X,par,lpi_jk)
map_crit_gauss_homoscd=function(y,X,par,lpi_jk){
  p = ncol(X)
  n = nrow(X)
  K = dim(lpi_jk)[1]
  beta=matrix(par[1:(K*p)],nrow=p)
  Xbeta=X%*%beta
  le = matrix(dnorm(y, mean = Xbeta, sd = par[K*p+1],log = TRUE),n,K) #log(e)
  resM = FwBw_factLog_Niv2(exp(lpi_jk), le, max,2)
  res1 = resM$res
  res2 = resM$Rseg
  bp_opt=which(diff(res2)==1)
  names(bp_opt)=NULL
  return(list(res=res1,seg=res2,bp_opt=bp_opt))
}


#' glmseg_gauss_homoscd
#'
#' @param y a vector of n observations (supposed gaussian and homoscedastic)
#' @param X a n x p matrix of covariates. The first column of X is filled with 1 and always corresponds to the intercept. Thus, if p = 1, then there are no covariates
#' @param K a parameter specifying that we want K-1 breakpoints
#' @param ... optional glm parameters (to be done completely !!!)
#'
#' @return list containing MaxEM results (log-likelihood and estimated parameters) and computed breakpiint(s)
#' @export
#' @importFrom stats glm sd gaussian
#' @importFrom glmnet glmnet
#' @importFrom gtools combinations
#' @importFrom bazar is.empty
#' @examples
#' n = 1000 ; X = matrix(1,n,1)
#' bp = c(553)
#' par = c(10,12,3)
#' y <- simuData_gauss_homoscd(X,par,bp)
#' K = 2
#' glmseg_gauss_homoscd(y,X,K)
glmseg_gauss_homoscd <- function(y,X,K, ...){
  glmOptionalParams <- list(...)
  if(!is.null(glmOptionalParams$weights)){ # corriger l'acces a l'element (glmOptionalParams$weights donne null ...)
    weights = glmOptionalParams$weights
    print("@@@@@@ weights @@@@@@")
    print(length(weights))
  }

  n = nrow(X)
  p = ncol(X)

  lpi_jk = matPijk(K)
  #================== COMBINAISON FUSED LASSO ==================
  # Fused-Lasso
  D = matrix(1,n,n)
  D[upper.tri(D)] = 0

  res.fused0 <- glmnet(x = D, y = y)

  ChoiceLambdaGoodDim <- apply(res.fused0$beta,2, function(x) sum(abs(x)>0))
  break.fused1 <- which(abs(res.fused0$beta[,min(which(ChoiceLambdaGoodDim >= 5*(K-1)))])> 0)

  break.fused <- break.fused1

  print("------------- candidats fused-lasso avant selection -------------")
  print(break.fused)

  s_tol <- 50
  size_bp_lasso_prov <- length(break.fused)
  iterwh <- 1

  while(size_bp_lasso_prov > (3*(K-1))/2 && iterwh < s_tol){
    break.fused <- almost_unique_tol(break.fused, tol = iterwh)
    size_bp_lasso_prov <- length(break.fused)
    iterwh <- iterwh + 1
  }

  break.fused <- break.fused[break.fused!= 1]
  break.fused <- break.fused[break.fused!= 2]
  break.fused <- break.fused[break.fused!= n-1]
  break.fused <- break.fused[break.fused!= n]

  print("------------- candidats fused-lasso apres selection -------------")
  print(break.fused)
  size_bp_lasso <- length(break.fused)

  Comb_bp_using_lasso = combinations(size_bp_lasso, K-1, v=break.fused, set=TRUE, repeats.allowed=FALSE)
  nb_l <- dim(Comb_bp_using_lasso)[1]
  nb_c <- dim(Comb_bp_using_lasso)[2]
  vectLikInit <- NULL
  print(paste("@@@@@@ ", nb_l,"INITIALISATIONS (VIA FUSED LASSO) @@@@@@"))

  for(ind_init in 1:nb_l){
    bpFused0=c(0,Comb_bp_using_lasso[ind_init,],n)
    fit=vector(nb_c+1,mode="list")
    for(k_fl in 1:(nb_c+1)){
      idx=(bpFused0[k_fl]+1):(bpFused0[k_fl+1])
      #---------- avec glm --------------
      fit[[k_fl]]= glm(y~.-1,family=gaussian,data=data.frame(y=y,X=X)[idx,])
    }

    par_init=sapply(fit,coefficients)
    sigFused = sd(y)
    par_init = c(par_init, sigFused)
    logLvalue = max_crit_gauss_homoscd(y,X,par_init,lpi_jk) # (estimation dans un contexte homoscedastique)
    #logLvalue = sum(sapply(fit,logLik)) # (estimation dans un contexte heteroscedastique, pour faire des tests)
    bp_MaxEM = c(map_crit_gauss_homoscd(y,X,par_init,lpi_jk)$bp_opt)

    bpFused0While=c(0,bp_MaxEM,n)
    segProb = which(diff(bpFused0While) < p) # < max(3,p)
    if(is.empty(segProb) && !is.empty(bp_MaxEM)){
      logLvalue_prec = logLvalue - 999

      prec_bp = Comb_bp_using_lasso[ind_init,]
      diff_bp = FALSE

      ItWhil <- 1
      ItWhilMax = 100

      tol = 1e-5
      #while(abs(logLvalue_prec-logLvalue) > tol && !diff_bp){
      while(abs(logLvalue_prec-logLvalue) > tol && ItWhil <= ItWhilMax && !diff_bp && is.empty(segProb) && !is.empty(bp_MaxEM)){

        bpFused1=c(0,bp_MaxEM,n)
        fit_whil=vector(nb_c+1,mode="list")
        for(k_fl in 1:(nb_c+1)){
          idxWhil=(bpFused1[k_fl]+1):(bpFused1[k_fl+1])
          #---------- avec glm --------------
          fit_whil[[k_fl]]= glm(y~.-1,family=gaussian,data=data.frame(y=y,X=X)[idxWhil,])
        }

        logLvalue_prec = logLvalue
        prec_bp = bp_MaxEM

        par_init=sapply(fit_whil,coefficients)
        par_init = c(par_init, sigFused)

        logLvalue = max_crit_gauss_homoscd(y,X,par_init,lpi_jk)
        #logLvalue = sum(sapply(fit_whil,logLik))
        bp_MaxEM = c(map_crit_gauss_homoscd(y,X,par_init,lpi_jk)$bp_opt)

        bpFused0While=c(0,bp_MaxEM,n)
        segProb = which(diff(bpFused0While) < p) # < max(3,p)

        diff_bp = length(bp_MaxEM) == length(prec_bp) && (bp_MaxEM == prec_bp)
        if(!diff_bp){
          if(length(bp_MaxEM) != length(prec_bp)){
            break
          }
        }
        ItWhil = ItWhil+1
      }

      if(!is.infinite(logLvalue)){
        opt2 <- list(par = par_init, value = logLvalue)

        print(paste("@@@ INITIALISATION NUMERO :",ind_init,"@@@"))
        print("logVraisemblance :")
        print(logLvalue)
        print("bpMaxEM :")
        print(bp_MaxEM)

        print("nbIter :")
        print(ItWhil)

        vectLikInit[[ind_init]] = opt2
      }
    }
  }
  vectLikInit <- vectLikInit[!sapply(vectLikInit,is.null)]

  indexInitOpti <- which.max(sapply(vectLikInit, function(xx) xx$value))
  opt2 <- vectLikInit[[indexInitOpti]]
  bp_MaxEM = c(map_crit_gauss_homoscd(y,X,opt2$par,lpi_jk)$bp_opt)

  return(list(resMaxEM = opt2, bp_MaxEM = bp_MaxEM))
}
