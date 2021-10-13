#' Pseudo pseudolikelihood ratio test for dependent communities
#'
#'
#' Implements the pseudo pseudolikelihood ratio test described in Section 3 of
#' Gao et. al. (2019) "Testing for Association in Multi-View Network Data"
#' for testing for dependence between communities in two network data views.
#' Fits stochastic block models in each view.
#'
#' @param X Multi-view data with two views; a list of two n x n adjacency matrices.
#' @param nperm An integer specifying the number of permutations to
#' use for the permutation procedure. The default number is 200.
#' @param step A numeric value containing the fixed step size to be used in the optimization
#' algorithm for estimating Pi. The default step size is 0.001. 
#' @param maxiter A numeric value containing the maximum number of iterations to run in
#'  the optimization algorithm. The default maximum is 1000.
#' @param K1 An optional argument containing the number of communities in View 1.
#'  If left out, then the number of communities is chosen with the method of 
#'  Le and Levina (2015).
#' @param K2 An optional argument containing the number of communities in View 2.
#'  If left out, then the number of communities is chosen with the method of 
#'  Le and Levina (2015). 
#' @param parallel An optional argument; if true, do parallel computing using the 
#' doParallel package 
#' @import randnet
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @export
#' @return
#' A list containing the following output components:
#' \item{K1}{The number of communities in view 1}
#' \item{K2}{The number of communities in view 2}
#' \item{Pi.est}{The estimated Pi matrix}
#' \item{P2LRstat}{The pseudo likelihood ratio test statistic}
#' \item{pval}{The p-value}
#' \item{modelfit1}{The parameter estimates and community assignment estimates from View 1.}
#' \item{modelfit2}{The parameter estimates and community assignment estimates from View 2.}
#' @examples
#' set.seed(1)
#' n <- 50
#' Pi <- diag(c(0.5, 0.5))
#' theta1 <- rbind(c(0.5, 0.1), c(0.1, 0.5))
#' theta2 <- cbind(c(0.1, 0.5), c(0.5, 0.1))
#' 
#' # 50 draws from a multi-view SBM with perfectly dependent communities
#' dat <- mv_sbm_gen(n, Pi, theta1, theta2)
#' 
#' # Test H0: communities are independent
#' # Data was generated under the alternative hypothesis 
#' results <- test_indep_com(dat$data, nperm=25)
#' results$pval
#'
#' @references
#' Amini, A. A., Chen, A., Bickel, P. J., & Levina, E. (2013). Pseudo-likelihood methods 
#' for community detection in large sparse networks. The Annals of Statistics, 41(4), 2097-2122.
#'
#' Gao, L.L., Witten, D., Bien, J. Testing for Association in Multi-View Network Data, preprint. 
#'
#' Le, C. M., & Levina, E. (2015). Estimating the number of communities 
#' in networks by spectral methods. arXiv preprint arXiv:1507.00827.
test_indep_com <- function(X, K1=NULL, K2=NULL, nperm=200, step=0.001, maxiter=1000, 
                           parallel=FALSE) {
  input_check <- function(X, K1, K2, nperm, step, maxiter, parallel) { 
    if(typeof(X) != "list") {
      stop("X is not in the right format; should be a list of two matrices")
    } 
    
    if(length(X) != 2) {
      stop("X is not in the right format; should be a list of two matrices")
    }
   
    if(!is.matrix(X[[1]]) | !is.numeric(X[[1]]) | !is.matrix(X[[2]]) | !is.numeric(X[[2]])) {
      stop("X is not in the right format; should be a list of two matrices")
    } 
    
    if(nrow(X[[1]]) != ncol(X[[1]]) | nrow(X[[2]]) != ncol(X[[2]]) | nrow(X[[1]]) != ncol(X[[2]])) { 
      stop("X is not in the right format; should be a list of two n x n matrices")  
    } 
    
    if(!isSymmetric(X[[1]]) | !isSymmetric(X[[2]])) { 
      stop("Networks should be undirected, so adjacency matrices should be symmetric")        
    }
    
    if(any(is.na(X[[1]])) || any(is.na(X[[2]])))
      stop(paste("Views cannot contain any missing values"))
    
    if(step <= 0) stop('Step size for optimization algorithm should be > 0')
    if(nperm < 1) stop('Number of permutations should be > 0')
    if(nperm <= 50) warning('Number of permutation iterations is specified to be small;
                     p-value approximation will likely be imprecise')
    if(maxiter %% 1 != 0|| maxiter <= 0) stop('Maximum number of iterations should be
                                                a positive integer')
    
    if(!is.null(K1)) { 
      if(!is.numeric(K1)) stop('Number of communities must be a positive integer in View 1')
      if(K1 %% 1 != 0 || K1 < 0) stop('Number of communities must be a positive integer in View 1')
      if(K1 < 1) stop('Number of communities must be greater than 1 in View 1')
    }
    
    if(!is.null(K2)) { 
      if(!is.numeric(K2)) stop('Number of communities must be a positive integer in View 2')
      if(K2 %% 1 != 0 || K2 < 0) stop('Number of communities must be a positive integer in View 2')
      if(K2 < 1) stop('Number of communities must be greater than 1 in View 2')
    }
    
    if(!is.logical(parallel)) stop("parallel should be TRUE or FALSE")
  } 
  
  input_check(X, K1, K2, nperm, step, maxiter, parallel)
  
  if(is.null(K1)) K1 <- max(randnet::BHMC.estimate(X[[1]])$K, 2)  
  if(is.null(K2)) K2 <- max(randnet::BHMC.estimate(X[[2]])$K, 2)  
  
  cluster.results1 <- pl_est_com(X[[1]], K1, parallel=parallel)
  cluster.results2 <- pl_est_com(X[[2]], K2, parallel=parallel)
  
  logphi1 <- cluster.results1$logphi
  logphi2 <- cluster.results2$logphi
  pi1.h <- cluster.results1$pi
  pi2.h <- cluster.results2$pi 

  # Computes diag(pi1.h)*C.h*diag(pi2.h)
  Pi.hat <- optimize_over_Pi(logphi1, logphi2, pi1.h, pi2.h,
                                 max.iter=maxiter, stepsz=step)

  # Compute P2LRT
  null.pll <- cluster.results1$ploglik + cluster.results2$ploglik 
  p2lr.stat <- (-Pi.hat$obj) - null.pll
  n <- nrow(X[[1]])
  
  do_one_perm <- function() { 
    perm <- sample(1:n, replace=F)
    Pi.hat.perm <- optimize_over_Pi(logphi1, logphi2[perm, ], 
                                     pi1.h, pi2.h, max.iter=maxiter, stepsz=step)
    -Pi.hat.perm$obj - null.pll
  }
  
  if(nperm <= 1000) { 
    if(parallel) { 
      perms <- foreach(i = 1:nperm, .combine=c) %dopar% do_one_perm()
    } else { 
      perms <- replicate(nperm, do_one_perm())  
      perms <- perms[!is.na(perms)] # remove any NA's
      pval <- (sum(perms >= p2lr.stat) + 1)/(length(perms) + 1)
      return(list(K1=K1, K2=K2, Pi.est=Pi.hat$Pi, P2LRstat=p2lr.stat, 
                  pval=pval, modelfit1=cluster.results1, modelfit2=cluster.results2))
    }
  } else { 
    if(parallel) { 
      interim.perm <- foreach::foreach(i = 1:1000, .combine=c) %dopar% do_one_perm()
    } else { 
      interim.perm <- replicate(1000, do_one_perm())  
    }
    
    interim.perm <- interim.perm[!is.na(interim.perm)] # Remove any NA's 
    interim.pval <- (sum(interim.perm >= p2lr.stat) + 1)/(length(interim.perm) + 1) 
  
    if(interim.pval > 0.1) { 
      cat("Stopped at 1000 permutations (p-value > 0.1)")
      return(list(K1=K1, K2=K2, Pi.est=Pi.hat$Pi, P2LRstat=p2lr.stat, 
                  pval=interim.pval, modelfit1=cluster.results1, modelfit2=cluster.results2))
    } 
  
    if(parallel) { 
      rest.perm <- foreach::foreach(i = 1001:nperm, .combine=c) %dopar% do_one_perm()
    } else { 
      rest.perm <- replicate(nperm - 1000 + 1, do_one_perm())  
    }
    
    final.perm <- c(interim.perm, rest.perm)
    final.perm <- final.perm[!is.na(final.perm)] # Remove any NA's 
    pval <- (sum(final.perm >= p2lr.stat) + 1)/(length(final.perm) + 1)
    
    return(list(K1=K1, K2=K2, Pi.est=Pi.hat$Pi, P2LRstat=p2lr.stat, 
                pval=pval, modelfit1=cluster.results1, modelfit2=cluster.results2))
  }
}