#' Generates data from a multi-view Gaussian mixture model
#'
#' Generates data from a multi-view Gaussian mixture model with n observations 
#' and two views. 
#'
#' @param n number of observations
#' @param Pi K1 x K2 matrix where the (k, k')th entry contains the probability of an observation 
#' belonging to cluster k in View 1 and cluster k' in View 2
#' @param mu1 p1 x K1 matrix where the columns contain the K1 cluster means in View 1
#' @param mu2 p2 x K2 matrix where the columns contain the K2 cluster means in View 2
#' @param Sigma1 p1 x p1 matrix containing the covariance matrix for View 1
#' @param Sigma2 p2 x p2 matrix containing the covariance matrix for View 2
#'
#' @importFrom stats rnorm
#' @export
#' @return 
#' A list containing the following components:
#' \item{data}{A list with two items: the view 1 n x p1 multivariate data set and 
#' the view 2 n x p2 multivariate data set }
#' \item{clusters}{A list with two items: the view 1 cluster memberships and 
#' the view 2 cluster memberships}
#' 
#' @examples
#' # 25 draws from a two-view Gaussian mixture model where the clusters are independent
#' n <- 25
#' Pi <-  tcrossprod(c(0.5, 0.5), c(0.25, 0.25, 0.5))
#' mu1 <- cbind(c(2, 2), c(-2, 2))
#' mu2 <- cbind(c(0, 1), c(1, 0), c(-1, 0))
#' Sigma1 <- diag(rep(1, 2))
#' Sigma2 <- diag(rep(0.5, 2))
#' 
#' mv_gmm_gen(n, Pi, mu1, mu2, Sigma1, Sigma2)
#' 
#' @references
#' Gao, L.L., Bien, J., Witten, D. (2019) Are Clusterings of Multiple Data Views Independent?
#' Biostatistics,  <DOI:10.1093/biostatistics/kxz001>
#'
#' Gao, L.L., Witten, D., Bien, J. Testing for Association in Multi-View Network Data, preprint. 
mv_gmm_gen <- function(n, Pi, mu1, mu2, Sigma1, Sigma2) {
  input_check <- function(n, Pi, mu1, mu2, Sigma1, Sigma2) { 
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
    
    if(!is.matrix(mu1) | !is.numeric(mu1) | !is.matrix(mu2) | !is.numeric(mu2)) { 
      stop("mu1 and mu2 must be numeric matrices")
    } else { 
      if((ncol(mu1) != nrow(Pi)) | (ncol(mu2) != ncol(Pi))) { 
        stop("mu1 must be n x K1 and mu2 must be n x K2, where Pi is K1 x K2")  
      }  
    }
    
    if(!is.matrix(Sigma1) | !is.numeric(Sigma1) | !is.matrix(Sigma2) | !is.numeric(Sigma2)) { 
      stop("Sigma1 and Sigma2 must be numeric matrices")  
    } 
    
    if((nrow(Sigma1) != ncol(Sigma1)) | (nrow(Sigma2) != ncol(Sigma2))) { 
        stop("Sigma1 and Sigma2 must be PSD")  
    } 
    
    if(any(eigen(Sigma1)$values < 0) | any(eigen(Sigma2)$values < 0)) { 
          stop("Sigma1 and Sigma2 must be PSD")  
    }
  }
  
  generate_sv_gaussian <- function(mu, Sig, cl) {
    n <- length(cl)
    p <- nrow(Sig)
    
    matrix(stats::rnorm(n * p), n, p)%*%chol(Sig) + t(mu[, cl])
  }
  
  input_check(n, Pi, mu1, mu2, Sigma1, Sigma2)
  
  cl <- mv_memberships_gen(n, Pi)
  X1 <- generate_sv_gaussian(mu1, Sigma1, cl[, 1])
  X2 <- generate_sv_gaussian(mu2, Sigma2, cl[, 2])
  
  return(list(data=list(view1=X1, view2=X2), clusters=list(view1=cl[, 1], view2=cl[, 2])))
}

