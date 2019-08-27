#' Fits the stochastic block model using maximum pseudolikelihood estimation
#'
#' Fits the stochastic block model using maximum
#' pseudolikelihood estimation, as proposed by Amini et. al. (2013). 
#' This function implements the conditional pseudolikelihood algorithm 
#' from Amini et al. (2013). 
#'
#' @param X n x n adjacency matrix
#' @param K number of communities; by default, chosen using the method of Le and Levina (2015)
#' @param max.iter maximum number of iterations for the EM algorithm 
#' @param tol the EM algorithm stops when the relative tolerance is less than this value
#' @param parallel An optional argument allowing for parallel computing using the 
#' doParallel package 
#'
#' @import randnet
#' @import irlba
#' @import doParallel
#' @import matrixStats
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom stats kmeans
#' @export
#' @return 
#' A list containing the following components:
#' \item{eta}{Estimate of eta, a K x K matrix defined in Amini et. al. (2013)}
#' \item{pi}{Estimate of the community membership probabilities}
#' \item{ploglik}{The maximum of the pseudolikelihood function}
#' \item{logphi}{n x K matrix, where (i, k)th entry contains the log p.m.f. of a multinomial
#' random variable with probability vector eta_k (the kth row of eta), evaluated at b_i, 
#' which is the ith row of the block compression matrix defined in Amini et. al. (2013)}
#' \item{responsibilities}{n x K matrix containing the responsibilities/soft cluster 
#' memberships for the n nodes}
#' \item{class}{A vector containing n (hard) cluster memberships for the n nodes}
#' \item{converged}{whether the algorithm converged to the desired tolerance}
#' 
#' @examples
#' # 50 draws from a stochastic block model for two network data views 
#' # where the communities are dependent
#' n <- 50
#' Pi <- diag(c(0.5, 0.5))
#' theta1 <- rbind(c(0.5, 0.1), c(0.1, 0.5))
#' theta2 <- cbind(c(0.1, 0.5), c(0.5, 0.1))
#' 
#' dat <- mv_sbm_gen(n, Pi, theta1, theta2)
#' 
#' # Fit SBM to view 1
#' results <- pl_est_com(X=dat$data$view1, K = 2)
#' table(results$class, dat$communities$view1)
#' 
#' @references
#' Amini, A. A., Chen, A., Bickel, P. J., & Levina, E. (2013). 
#' Pseudo-likelihood methods for community detection in large sparse networks. 
#' The Annals of Statistics, 41(4), 2097-2122.
#' 
#' Le, C. M., & Levina, E. (2015). Estimating the number of communities 
#' in networks by spectral methods. arXiv preprint arXiv:1507.00827.

pl_est_com <- function(X, K=NULL, max.iter = 1000, tol = 1e-8, parallel=FALSE) {
  input_check <- function(X, K, max.iter, tol, parallel) { 
    if(!is.numeric(X) | !is.numeric(max.iter) | !is.numeric(tol)) { 
      stop("All inputs must be numeric")  
    }
    
    if(length(max.iter) > 1 | length(tol) > 1) { 
      stop("max.iter and tol should be scalars")  
    }
    
    if(is.vector(X)) stop("X is not in the right format; should be a n x n matrix") 
    if(nrow(X) != ncol(X)) stop("X is not in the right format; should be a n x n matrix") 
    if(!isSymmetric(X)) stop("X must be a symmetric n x n matrix")
    if(!is.null(K)) { 
      if(!is.numeric(K)) stop("Number of communities must be a positive integer")
      if(length(K) > 1) stop("K should be scalar") 
      if(K %% 1 != 0 | K < 0) stop('Number of communities must be a positive integer')
      if(K < 1) stop('Number of communities must be greater than 1')
    }
    if(max.iter %% 1 != 0 | max.iter <= 0) stop("max number of iterations must be a positive integer")
    if(tol <= 0) stop("tolerance for convergence criteria must be positive")
    if(!is.logical(parallel)) stop("parallel should be T or F")
  }
  
  input_check(X, K, max.iter, tol, parallel)
  
  # if K is null, estimate using spectral method
  if(is.null(K)) K <- max(randnet::BHMC.estimate(X)$K, 2)
  
  # Zero degree nodes
  zd <- rowSums(X) == 0
  if(K > sum(!zd)) stop("Number of communities should be greater than the number of non-zero nodes.")
  
  # Save original data 
  orig.X <- X 
  orig.n <- nrow(X)
  
  # Work on data with zero degree nodes removed - we will handle the 
  # zero degree nodes separately at the end
  X <- X[!zd, !zd]
  n <- nrow(X) 
  
  # A parallelizable implementation of spectral clustering with perturbations
  reg.SP.parallel <- function(A, K, tau=0.25, lap=TRUE, parallel=TRUE) { 
    avg.d <- mean(colSums(A))
    A.tau <- A + tau * avg.d/nrow(A)
    if (!lap) {
      SVD <- irlba::irlba(A.tau, nu = K, nv = K)
    } else {
      d.tau <- colSums(A.tau)
      L.tau <- diag(1/sqrt(d.tau)) %*% A.tau %*% diag(1/sqrt(d.tau))
      SVD <- irlba::irlba(L.tau, nu = K, nv = K)
    }
    
    one_kmeans <- function() { 
      km <- kmeans(SVD$v[, 1:K], centers=K, nstart=1, iter.max=30)  
      return(list(cluster=km$cluster, loss=km$tot.withinss))
    }
    
    if(parallel) { 
      comb <- function(x, ...) {
        lapply(seq_along(x),
               function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
      }
      
      km.res <- foreach::foreach(i = 1:30, .combine='comb', .multicombine=TRUE,
                      .init=list(list(), list())) %dopar% one_kmeans()
      best.km.res <- list(cluster=unlist(km.res[[1]][which.min(unlist(km.res[[2]]))]), 
                          loss=min(unlist(km.res[[2]])))
    } else { 
      km.res <- replicate(30, one_kmeans())
      best.km.res <- list(cluster=km.res[1, which.min(km.res[2, ])]$cluster, 
                          loss=km.res[2, which.min(km.res[2, ])]$loss)
    }
    
    return(list(cluster=best.km.res$cluster, loss=best.km.res$loss))
  }
  
  # Initialize the pseudolikelihood algorithm
  if(K < n) { 
    Z.hat <- reg.SP.parallel(X, K, parallel=parallel)$cluster 
  } else { 
    Z.hat <- seq(1, n)  
  } 
  
  M.hat <- matrix(0, n, K)
  M.hat[cbind(1:n, Z.hat)] <- 1
  b.hat <- crossprod(X, M.hat) 
  
  # initialize EM algorithm 
  count.h <- tabulate(Z.hat)
  pi.h <- count.h/n
  
  thet.h <- matrix(0, K, K)
  for(k1 in 1:K) { 
    for(k2 in 1:K) { 
      thet.h[k1, k2] <- mean(X[Z.hat == k1, Z.hat == k2, drop=F])  
    }  
  }
  
  lam.h <- n*t(thet.h*pi.h)
  eta.h <- lam.h/rowSums2(lam.h)
  eta.h[eta.h == 0] <- 1e-16
  
  converged <- F
  iter <- 1
  while((iter <= max.iter) && !converged) { 
    # Compute objective function 
    log.phi.h <- tcrossprod(b.hat, log(eta.h))
    log.phi.pi.h <- t(log.phi.h) + log(pi.h) 
    cpl <- sum(colLogSumExps(log.phi.pi.h))
    
    # Check tolerance
    if(iter > 1) { 
      diffcpl <- abs((cplOld - cpl)/cplOld)
      converged <- (diffcpl < tol)
    }
    
    # E-step 
    log.phi.pi.h.norm <- t(log.phi.pi.h) - colMaxs(log.phi.pi.h)
    gam.h <- exp(log.phi.pi.h.norm)/rowSums2(exp(log.phi.pi.h.norm))
    
    # M-step 
    pi.h <- colMeans2(gam.h)
    eta.h <- t(gam.h)%*%b.hat
    eta.h <- eta.h/rowSums2(eta.h)
    eta.h[eta.h == 0] <- 1e-16
    
    cplOld <- cpl
    iter <- iter + 1
  }
  
  # Estimates from EM
  Z.hat <- max.col(gam.h)
  log.phi.h <- tcrossprod(b.hat, log(eta.h))
  log.phi.pi.h <- t(log.phi.h) + log(pi.h) 
  cpl.est <- sum(colLogSumExps(log.phi.pi.h))
  
  # Deal with any zero degree nodes
  if(sum(zd) > 0) { 
    log.phi.h.final <- matrix(NA, orig.n, K)
    log.phi.h.final[!zd, ] <- log.phi.h
    log.phi.h.final[zd, ] <- 0
    
    log.phi.pi.h <- t(log.phi.h) + log(pi.h) 
    log.phi.pi.h.norm <- t(log.phi.pi.h) - colMaxs(log.phi.pi.h)
    
    gam.h.final <- matrix(NA, orig.n, K) 
    gam.h.final[!zd, ] <- exp(log.phi.pi.h.norm)/rowSums2(exp(log.phi.pi.h.norm))
    gam.h.final[zd, ] <- pi.h
    
    Z.hat.final <- rep(NA, orig.n)
    Z.hat.final[!zd] <- Z.hat
    # Randomly assign hard clusters for the zero degree nodes 
    Z.hat.final[zd] <- sample(1:K, sum(zd), replace=T, prob=pi.h)
    
    Z.hat <- Z.hat.final
    log.phi.h <- log.phi.h.final
    gam.h <- gam.h.final
    cpl.est <-  sum(colLogSumExps(t(log.phi.h.final) + log(pi.h) ))
  }
  
  return(list(logphi=log.phi.h, pi=pi.h, eta=eta.h, class=Z.hat, 
              responsibilities=gam.h, ploglik=cpl.est, converged=converged))
}
  