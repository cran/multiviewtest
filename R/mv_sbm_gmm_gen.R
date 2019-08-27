#' Generates data from a stochastic block model for a network view and a multivariate view
#'
#' Generates data from a stochastic block model for a network view and a multivariate view
#' with n observations. The data for the multivariate view is drawn from a Gaussian
#' mixture model. 
#'
#' @param n number of observations
#' @param Pi K1 x K2 matrix where the (k, k')th entry contains the probability of an observation 
#' belonging to community k in View 1 and cluster k' in View 2
#' @param theta1 K1 x K1 matrix containing the between-community edge probabilities for View 1
#' @param mu2 mu2 p2 x K2 matrix where the columns contain the K2 cluster means in View 2
#' @param Sigma2 p2 x p2 matrix containing the covariance matrix for View 2
#' @param sparse If true, return matrix views in sparseMatrix format
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats rnorm
#' @export
#' @return 
#' A list containing the following components:
#' \item{data}{A list with two items: the view 1 n x n adjacency matrix and 
#' the view 2 n x p multivariate data set }
#' \item{communities}{A list with two items: the view 1 community memberships and 
#' the view 2 cluster memberships}
#' 
#' @examples
#' # 50 draws from a stochastic block model for a network view and a multivariate view
#' # where the communities and the clusters are independent
#' n <- 50
#' Pi <- tcrossprod(c(0.5, 0.5), c(0.5, 0.5))
#' theta1 <- rbind(c(0.5, 0.1), c(0.1, 0.5))
#' mu2 <- cbind(c(2, 2), c(-2, 2))
#' Sigma2 <- diag(rep(0.5, 2))
#' 
#' mv_sbm_gmm_gen(n, Pi, theta1, mu2, Sigma2)
#' 
#' @references
#' Gao, L.L., Witten, D., Bien, J. Testing for Association in Multi-View Network Data, preprint. 
mv_sbm_gmm_gen <- function(n, Pi, theta1, mu2, Sigma2, sparse=FALSE) {
  input_check <- function(n, Pi, theta1, mu2, Sigma2) { 
    if(!is.matrix(Pi) | !is.numeric(Pi)) {
      stop("Pi must be a numeric matrix")  
    } else { 
      if(any(Pi < 0) | (abs(sum(Pi) - 1) > 1e-15)) { 
        stop("Pi must be non-negative and sum to 1")  
      }  
    }
    
    if(!is.numeric(n)) { 
      stop("n must be a positive integer")
    } else { 
      if((n %% 1 != 0) | n <= 0) { 
        stop("n must be a positive integer")  
      }  
    }
    
    if(!is.matrix(theta1) | !is.numeric(theta1) | !is.numeric(mu2) | !is.matrix(mu2)) { 
      stop("theta1 and mu2 must be numeric matrices")
    } else { 
      if((ncol(theta1) != nrow(Pi)) | nrow(theta1) != nrow(Pi) | 
         (ncol(mu2) != ncol(Pi))) { 
        stop("theta1 must be K1 x K1 and mu2 must be p2 x K2, where Pi is K1 x K2")  
      } else { 
        if(any(theta1 < 0) | any(theta1 > 1)) { 
          stop("theta1 must have non-negative elements that are < 1")  
        } else { 
          if(!isSymmetric(theta1)) { 
            stop("theta1 must be a symmetric matrix")  
          }    
        }
      }
    }
    
    if(!is.matrix(Sigma2) | !is.numeric(Sigma2)) { 
      stop("Sigma2 must be a numeric matrix")  
    } else { 
      if((nrow(Sigma2) != ncol(Sigma2))) { 
        stop("Sigma2 must be PSD")  
      } else { 
        if(any(eigen(Sigma2)$values < 0)) { 
          stop("Sigma2 must be PSD")  
        }
      }
    }
    
  }
  
  # output: random adjacency matrix X with E[X] = pmat 
  generate_random_graph <- function(pmat) { 
    n <- nrow(pmat)
    upper.ind <- which(upper.tri(pmat), arr.ind=T)
    upper.unif <- stats::runif(n = nrow(upper.ind))
    edge.ind <- upper.ind[which(upper.unif < pmat[upper.ind]), ]
    Matrix::sparseMatrix(edge.ind[, 1], edge.ind[, 2], dims=c(n, n), symmetric=T)
  }
  
  generate_sv_sbm <- function(theta, com) {
    n <- length(com)
    K <- nrow(theta)
    
    # expected value of adjacency matrix: 
    # n x n matrix with ijth entry = theta[com[i], com[j]]
    outer <- matrix(0, n, K)
    outer[cbind(1:n, com)] <- 1
    outer <- outer
    pmat <- outer%*%tcrossprod(theta, outer)
    
    generate_random_graph(pmat)
  }
  
  generate_sv_gaussian <- function(mu, Sig, cl) {
    n <- length(cl)
    p <- nrow(Sig)
    
    matrix(stats::rnorm(n * p), n, p)%*%chol(Sig) + t(mu[, cl])
  }
  
  input_check(n, Pi, theta1, mu2, Sigma2)
  
  subgrp <- mv_memberships_gen(n, Pi)
  X1 <- generate_sv_sbm(theta1, subgrp[, 1])
  X2 <- generate_sv_gaussian(mu2, Sigma2, subgrp[, 2])
  
  if(!sparse) { 
    X1 <- as.matrix(X1)*1  
  }
  
  return(list(data=list(view1=X1, view2=X2), subgroups=list(view1=subgrp[, 1], view2=subgrp[, 2])))
}