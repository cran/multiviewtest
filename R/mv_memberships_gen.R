#' Generates multi-view subgroup memberships
#'
#' Generates subgroup membership pairs (Z1, Z2) for a multi-view data set with 
#' n observations and two views. 
#'
#' @param n number of observations
#' @param Pi K1 x K2 matrix where the (k, k')th entry contains the probability of an observation 
#' belonging to subgroup k in View 1 and subgroup k' in View 2
#'
#' @export
#' @return n x 2 matrix where the first column contains subgroup memberships for View 1 (Z1)
#' and the second column contains subgroup memberships for View 2 (Z2). 
#' 
#' @examples
#' Pi <- tcrossprod(c(0.5, 0.5), c(0.25, 0.25, 0.5))
#' n <- 25
#' mv_memberships_gen(n, Pi)
#' 
#' @references
#' Gao, L.L., Bien, J., Witten, D. (2019) Are Clusterings of Multiple Data Views Independent?
#' to appear in Biostatistics, <DOI:10.1093/biostatistics/kxz001>
#'
#' Gao, L.L., Witten, D., Bien, J. Testing for Association in Multi-View Network Data, preprint. 
mv_memberships_gen <- function(n, Pi) {
  sanity_check <- function(n, Pi) { 
    if(!is.matrix(Pi)) {
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
  }
  
  sanity_check(n, Pi)
  
  ii <- sample(length(Pi), n, replace = TRUE, prob = Pi)
  K1 <- nrow(Pi)
  K2 <- ncol(Pi)
  cpairs <- cbind(rep(1:K1, K2), rep(1:K2, each = K1))
  memberships <- cpairs[ii, ]
  
  return(memberships)
}