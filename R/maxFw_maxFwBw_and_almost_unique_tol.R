#' Fw_factLog_Niv2
#'
#' @param pi_jk a KxK matrix specifying that we want K-1 breakpoints (to be simplified !!!)
#' @param le a n*K matrix corresponding to the logarithm of the emission probabilities
#' @param FUN an operator allowing to choose between Forward (FUN = sum) and MaxForward (FUN = max) algorithms. We use MaxForward algorithm here
#'
#' @return log-likelihood
#' @export
Fw_factLog_Niv2 = function(pi_jk,le,FUN){
  dim = dim(le)
  n=dim[1]
  K=dim[2]
  Fw =matrix(NA,n,K)
  L = rep(NA,n)

  l=apply(le,1,max)
  e_fL=exp(le-l)

  Fw[1,1]=e_fL[1,1]; Fw[1,-1]=0

  tmp=max(Fw[1,])
  L[1]=l[1]+log(tmp)
  Fw[1,]=Fw[1,]/tmp

  for(i in 2:n){
    for(k in 1:K){
      Fw[i,k]=FUN(Fw[i-1,]*pi_jk[,k]*e_fL[i,k])
    }
    tmp=max(Fw[i,])
    L[i]=L[i-1]+l[i]+log(tmp)
    Fw[i,]=Fw[i,]/tmp
  }

  logRes = L+log(Fw); return(logRes[n,K]) # loglikelihood
}



#' FwBw_factLog_Niv2
#'
#' @param pi_jk a KxK matrix specifying that we want K-1 breakpoints (to be simplified !!!)
#' @param le a n*K matrix corresponding to the logarithm of the emission probabilities
#' @param FUN an operator allowing to choose between Forward (FUN = sum) and MaxForward (FUN = max) algorithms. We use MaxForward algorithm here
#' @param algo a parameter allowing to choose between Forward-Backward (algo = 1) and MaxForward-MaxBackward (algo = 2) algorithms. We use MaxForward-MaxBackward algorithm here
#'
#' @return log-likelihood, marginal probabilities and segmentation
#' @export
FwBw_factLog_Niv2 = function(pi_jk, le, FUN, algo){ # algo = 1 pour FwBw et 2 pour MaxFwBw
  dim = dim(le)
  n=dim[1]
  K=dim[2]
  Fw_fL =matrix(NA,n,K)
  L = rep(NA,n)

  l=apply(le,1,max)
  e_fL=exp(le-l)

  Fw_fL[1,1]=e_fL[1,1]; Fw_fL[1,-1]=0
  tmp = max(Fw_fL[1,])
  L[1] = l[1]+log(tmp)
  Fw_fL[1,]= Fw_fL[1,]/tmp
  for(i in 2:n){
    for(k in 1:K){
      Fw_fL[i,k]  = FUN(Fw_fL[i-1,]*pi_jk[,k]*e_fL[i,k])
    }
    tmp = max(Fw_fL[i,])
    L[i] = L[i-1] + l[i] + log(tmp)
    Fw_fL[i,]= Fw_fL[i,]/tmp
  }

  Bw_fL=matrix(NA,n,K)
  M = rep(NA,n)
  Bw_fL[n,]=1
  tmp = max(Bw_fL[n,])
  M[n] = log(tmp)
  Bw_fL[n,]= Bw_fL[n,]/tmp
  for(i in n:2){
    for(j in 1:K){
      Bw_fL[i-1,j]=FUN(pi_jk[j,]*e_fL[i,]*Bw_fL[i,])
    }
    tmp = max(Bw_fL[i-1,])
    M[i-1] = M[i] + l[i]+ log(tmp)
    Bw_fL[i-1,]= Bw_fL[i-1,]/tmp
  }

  res_fL = Fw_fL*Bw_fL
  res2_fL = exp(L+M)*Fw_fL*Bw_fL
  w = res_fL/apply(res_fL,1,FUN)
  if(algo == 1){
    return(list(res = res2_fL[n,K], w = w)) # On renvoie les probabilités marginales
  }else{
    R = apply(res_fL,1,which.max)
    return(list(res = res2_fL[n,K], w = w, Rseg = R)) # On renvoie surtout la segmentation optimale
  }
}

#' almost_unique_tol
#'
#' @param vect a vector
#' @param tol tolerance
#'
#' @return a vector
#' @export
#'
#' @examples
#' vect <- c(69, 70, 77, 81, 85, 110)
#' tol = 2
#' almost_unique_tol(vect, tol)
almost_unique_tol <- function(vect, tol){
  l = length(vect)
  d = diff(vect)
  ind  = which(d <= tol)
  l_ind = length(ind)
  if(l_ind == 0){
    return(vect)
  }else{
    v_res <- NULL
    i = 1
    v_res <- c(v_res,vect[i])
    while(i < l){
      depart = vect[i]
      nv <- which(vect - depart > tol)[1]
      if(is.na(nv)){
        break
      }else{
        v_res <- c(v_res, vect[which(vect - depart > tol)[1]])
        i = nv
      }
    }
    return(v_res)
  }
}


#' matPijk
#'
#' @param K represents the number of segments (K-1 breakpoints)
#'
#' @return a KxK ("transition") matrix (in logarithmic scale)
#' @export
#'
#' @examples
#' matPijk(4)
matPijk <- function(K){
  pi_jk <- matrix(rep(0,K*K),nrow = K)
  diag(pi_jk) <- rep(1,K)
  for(i in 1:(K-1)){
    pi_jk[i,i+1] <- 1
  }
  lpi_jk = log(pi_jk)
  return(lpi_jk)
}
