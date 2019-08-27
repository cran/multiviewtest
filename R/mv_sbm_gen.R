#' Generates data from a stochastic block model for multiple network data views
#'
#' Generates data from a stochastic block model for multiple network data views 
#' with n observations and two views.  
#'
#' @param n number of observations
#' @param Pi K1 x K2 matrix where the (k, k')th entry contains the probability of an observation 
#' belonging to community k in View 1 and community k' in View 2
#' @param theta1 K1 x K1 matrix containing the between-community edge probabilities for View 1
#' @param theta2 K2 x K2 matrix containing the between-community edge probabilities for View 2
#' @param sparse If true, return data views in sparseMatrix format
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats runif
#' @export
#' @return 
#' A list containing the following components:
#' \item{data}{A list with two items: the n x n view 1 adjacency matrix and 
#' the n x n view 2 adjacency matrix}
#' \item{communities}{A list with two items: the view 1 community memberships and 
#' the view 2 community memberships}
#' 
#' @examples
#' # 50 draws from a stochastic block model for two network data views 
#' # where the communities are dependent
#' n <- 50
#' Pi <- diag(c(0.5, 0.5))
#' theta1 <- rbind(c(0.5, 0.1), c(0.1, 0.5))
#' theta2 <- cbind(c(0.1, 0.5), c(0.5, 0.1))
#' 
#' mv_sbm_gen(n, Pi, theta1, theta2)
#' 
#' @references
#' Gao, L.L., Witten, D., Bien, J. Testing for Association in Multi-View Network Data, preprint. 
mv_sbm_gen <- function(n, Pi, theta1, theta2, sparse=FALSE) {
  input_check <- function(n, Pi, theta1, theta2) { 
    if(!is.matrix(Pi) | !is.numeric(Pi)) {
      stop("Pi must be a numeric matrix")  
    } 
    
    if(any(Pi < 0) | (abs(sum(Pi) - 1) > 1e-15)) { 
      stop("Pi must be non-negative and sum to 1")  
    }
    
    if(!is.numeric(n)) { 
      stop("n must be a positive integer")
    } 
    
    if((n %% 1 != 0) | n <= 0) { 
        stop("n must be a positive integer")  
    }
    
    if(!is.matrix(theta1) | !is.numeric(theta1) | !is.matrix(theta2) | !is.numeric(theta2)) { 
      stop("theta1 and theta2 must be numeric matrices")
    } 
    
    if((ncol(theta1) != nrow(Pi)) | nrow(theta1) != nrow(Pi) | 
         nrow(theta2) != ncol(Pi) | (ncol(theta2) != ncol(Pi))) { 
        stop("theta1 must be K1 x K1 and theta2 must be K2 x K2, where Pi is K1 x K2")  
    }
    
    if(any(theta1 < 0) | any(theta1 > 1) | any(theta2 < 0) | any(theta2 > 1)) { 
          stop("theta1 and theta2 must have non-negative elements that are < 1")  
    } 
    
    if(!isSymmetric(theta1) | !isSymmetric(theta2)) { 
            stop("theta1 and theta2 must be symmetric matrices")  
    }   
  }
  
  # outputs random adjacency matrix X with E[X] = pmat 
  generate_random_graph <- function(pmat) { 
    n <- nrow(pmat)
    upper.ind <- which(upper.tri(pmat), arr.ind=T)
    upper.unif <- runif(n = nrow(upper.ind))
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
  
  input_check(n, Pi, theta1, theta2)
  
  com <- mv_memberships_gen(n, Pi)
  X1 <- generate_sv_sbm(theta1, com[, 1])
  X2 <- generate_sv_sbm(theta2, com[, 2])
  
  if(!sparse) { 
    X1 <- as.matrix(X1)*1 # convert to regular matrix object
    X2 <- as.matrix(X2)*1
  }
  
  return(list(data=list(view1=X1, view2=X2), communities=list(view1=com[, 1], view2=com[, 2])))
}
