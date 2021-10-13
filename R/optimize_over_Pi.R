#' Exponentiated gradient descent for estimating Pi
#'
#'
#' Implements the optimization algorithm for solving equation (2.8)
#' in Section 2.3.2 of Gao et al. (2019) "Are Clusterings of Multiple
#' Data Views Independent?" Derivation of the algorithm is given in Appendix B.
#'
#' @param logphi1 log(phi1), where the n x K1 matrix phi1 is defined in equation (2.9)
#' @param logphi2 log(phi1), where the n x K2 matrix phi2 is defined in equation (2.9)
#' @param row K1-vector containing the estimated View 1 mixture component probabilities
#' @param col K2-vector containing the estimated View 2 mixture component probabilities
#' @param max.iter Maximum number of iterations to be run.
#' @param stepsz Fixed step size to be used in the optimization; see Appendix B for details.
#' @param Pi.init Initializes the optimization algorithm with Pi.init. (Optional)
#'
#' @import matrixStats
#' @export
#' @return
#' List of:
#' \item{Pi.est}{Estimate of Pi; maximizes the log-likelihood function of X1 and X2.}
#' \item{obj}{The log-likelihood function evaluated at Pi.est.}
#'
#' @examples
#' # Generate two-view Gaussian mixture model data with Pi = I/3
#' set.seed(1)
#' n <- 100
#' K <- 3
#' p <- 2
#' Sigma <- diag(rep(0.5^2, p))
#' Pi <- diag(rep(1, K))/K
#' mu1 <- cbind(c(2, 0), c(0, 2),  c(2, -2))
#' mu2 <- cbind(c(-2, 0), c(0, -2), c(-2, 2))
#' dat <- mv_gmm_gen(n, Pi, mu1, mu2, Sigma, Sigma)
#' view1dat <- dat$data$view1
#' view2dat <- dat$data$view2

# Cluster each view
#' library(mclust)
#' EM.View1 <- Mclust(view1dat, G=K, modelNames=c("EII"))
#' EM.View2 <- Mclust(view2dat, G=K, modelNames=c("EII"))

# Constructs inputs to algorithm
#' logphi1 <-  cdens(modelName="EII", data=view1dat, logarithm=TRUE, parameters=EM.View1$parameters)
#' logphi2 <-  cdens(modelName="EII", data=view2dat, logarithm=TRUE, parameters=EM.View2$parameters)
#' pi1.est <- EM.View1$parameters$pro
#' pi2.est <- EM.View2$parameters$pro

# Estimates Pi
#' Pi.est <- optimize_over_Pi(logphi1, logphi2, pi1.est, pi2.est)
#' Pi.est$Pi
#' 
#' @references
#' Gao, L.L., Bien, J., Witten, D. (2019) Are Clusterings of Multiple Data Views Independent?
#' to appear in Biostatistics, <DOI:10.1093/biostatistics/kxz001>
#' 
#' Gao, L.L., Witten, D., Bien, J. Testing for Association in Multi-View Network Data, preprint. 

optimize_over_Pi <- function(logphi1, logphi2, row, col,
                              max.iter = 1000, stepsz = 0.001, Pi.init=NULL) {
  
  input_check <- function(logphi1, logphi2, row, col, max.iter, stepsz) { 
    # Sanity checks for input
    if(!is.numeric(logphi1) | !is.numeric(logphi2) | !is.numeric(row) |
       !is.numeric(col) | !is.numeric(max.iter) | !is.numeric(stepsz)) {
      stop("all inputs must be numeric")
    }
    
    if(is.null(ncol(logphi1)) | is.null(ncol(logphi2))) {
      stop("logphi1 and logphi2 must be matrices")
    }
    
    if(length(row) == 1 | length(col) == 1) {
      stop("number of clusters in the 2 views must be > 1")
    }
    
    if(any(row <= 0) | any(col <= 0)) {
      stop("cluster membership probabilities must be positive")
    }
    
    if(abs(sum(row) - 1) > 1e-10 | abs(sum(col) - 1) > 1e-10) {
      stop("cluster membership probabilities in each view must sum to 1")
    }
    
    if(max.iter %% 1 != 0) {
      stop("max number of iterations must be an integer")
    }
    
    if(max.iter <= 0) {
      stop("max number of iterations must be > 0")
    }
    
    if(stepsz <= 0) {
      stop("step size must be positive")
    }
    
    if(nrow(logphi1) != nrow(logphi2)) {
      stop("logphi1 and logphi2 should have the same number of rows")
    }
  }
  
  input_check(logphi1, logphi2, row, col, max.iter, stepsz)
  
  K1 <- length(row)
  K2 <- length(col)
  
  if(ncol(logphi1) != K1 || ncol(logphi2) != K2) {
    stop(paste("number of columns in logphi1 and length of row must agree &",
               "number of columns in logphi2 and length of col must agree"))
  }
  
  # Initialize
  if(is.null(Pi.init)) { 
    Pi <- matrix(1 / (K1 * K2), K1, K2)
  } else { 
    Pi <- Pi.init  
  }
  obj <- rep(NA, max.iter)
  
  # For "vectorizing" the problem: (makes computation easier)
  # Construct n x (K1 + K2) matrices where each row contains
  # the combination of the kth mixture component density in View 1
  # and the k'th mixture component density in View 2
  log.phis <- logphi1[, rep(seq(K1), K2)] +
    logphi2[, rep.int(seq(K2), rep.int(K1, K2))]
  tlog.phis <- t(log.phis)
  
  for(iter in seq(max.iter))  {
    # Compute joint log-likelihood function
    log.phis.Pi <- tlog.phis + as.vector(log(Pi)) # K1*K2 x n
    obj.vect <- matrixStats::colLogSumExps(log.phis.Pi)
    obj[iter] <- -sum(obj.vect)
    
    # If the relative tolerance is less than 1e-10, terminate
    if(abs((obj[iter] - obj[iter-1])/obj[iter-1]) < 1e-8 && iter > 1) {
      break
    }
    
    logG <- matrix(matrixStats::colLogSumExps(log.phis - obj.vect), K1, K2)
    
    logM <- log(Pi) - 1 + exp(log(stepsz) + logG)
    M <- exp(logM)
    
    # Sinkhorn-Knopp algorithm from Cuturi (2013) to rescale M
    v <- rep(1, K2)
    u <- rep(1, K1)
    
    scale.M <- M/row
    M.t <- t(M)
    for(i in seq(100)) {
      uprev <- u
      u <- 1/scale.M%*%(col/M.t%*%u)
      
      if(max(abs(u/uprev - 1)) < 1e-4 && i > 1) {
        break
      }
    }
    
    v <- as.numeric(col/(M.t%*%u))
    
    # Update Pi
    Pi <- as.numeric(u)*M%*%diag(v)
  }
  
  list(Pi = Pi, obj = obj[iter])
}

