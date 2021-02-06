Gibbs_sampling <- function(R, X_U = NULL, X_V = NULL, D = 32, totalepoch = 300, burnin = 200, sn_max = 10, sn_init = 1, method){
  #-----------------------parameters----------------------#
  #     R                   # matrix that should be imputed
  #     X_U = NULL          # U feature matrix
  #     X_V = NULL          # V feature matrix
  #     D = 32              # num_latent
  #     alpha = 0.7         # precision
  #     totalepoch = 300    # total number of epoch
  #     burnin = 200        # number of burn-in epoch
  #     sn_max = 10         # for adaptive precision
  #     sn_init = 1         # for adaptive precision
  #     0             # corresponding

  # default hyper-parameters #
  mu_0 = matrix(rep(0,D))
  nu_0 = D
  W_0 = diag(rep(1,D))
  beta_0 = 2
  mu = 1
  nu = 1

  #-----------------------Initialize----------------------#
  N = nrow(R)
  M = ncol(R)

  if(is.null(X_U)){
    X_U = t(rep(0,N))
  }
  if(is.null(X_V)){
    X_V = t(rep(0,M))
  }

  #-------------------------------------------------------#

  U = 0.1*matrix(rnorm(D*N),nrow = D)     # initialize latent matrix for genes
  V = 0.1*matrix(rnorm(D*M),nrow = D)     # initialize latent matrix for cells

  F_U = nrow(X_U)
  F_V = nrow(X_V)

  beta_U = 0.1*matrix(rnorm(F_U*D),nrow = F_U)
  beta_V = 0.1*matrix(rnorm(F_V*D),nrow = F_V)

  lambda_beta_U = abs(rnorm(1))
  lambda_beta_V = abs(rnorm(1))

  R_calculate = matrix(rep(0, M*N), nrow = N)   # store the result

  #--------------Calculate initial precision--------------#
  # var_all = mean(R[setB].^2)
  ff = (R != 0)                 ####
  var_all = sd(R[ff])^2
  alpha = (sn_init+1)/var_all
  alpha_max = (sn_max)/var_all

  #-------------------Impute using macau------------------#
  for(epoch in 1:totalepoch){
    list_ml = list()
    # Sampling #
    list_ml = sample_mu_U_Lambda_U(U, X_U, beta_U, lambda_beta_U, mu_0, nu_0, W_0, beta_0)
    mu_U = list_ml[[1]]
    Lambda_U = list_ml[[2]]

    list_ml = sample_mu_U_Lambda_U(V, X_V, beta_V, lambda_beta_V, mu_0, nu_0, W_0, beta_0)
    mu_V = list_ml[[1]]
    Lambda_V = list_ml[[2]]

    U = sample_U0(U, R, V, X_U, beta_U, Lambda_U, mu_U, alpha)                 ####
    V = sample_U0(V, t(R), U, X_V, beta_V, Lambda_V, mu_V, alpha)                 ####

    U = sample_U0(U, R, V, X_U, beta_U, Lambda_U, mu_U, alpha)                 ####
    V = sample_U0(V, t(R), U, X_V, beta_V, Lambda_V, mu_V, alpha)                 ####

    lambda_beta_U = sample_lambda_beta_U(beta_U, Lambda_U, mu, nu)
    lambda_beta_V = sample_lambda_beta_U(beta_V, Lambda_V, mu, nu)

    beta_U = sample_beta_U(mu_U, Lambda_U, U, X_U, lambda_beta_U)
    beta_V = sample_beta_U(mu_V, Lambda_V, V, X_V, lambda_beta_V)

    Rsample1 = t(U) %*% V

    # Update precision alpha
    alpha = sample_alpha0(R,Rsample1,alpha_max)                 ####
    # print(paste0('Prec:',alpha))

    # Assume we reach the stationary distribution of a Markov chain after
    # burn-in and gain a bunch of samples. Then consider the expectation as
    # the result
    if(epoch > burnin & epoch > 1){
      R_calculate = (R_calculate*(epoch-burnin-1)+Rsample1)/(epoch-burnin)
    }
    else{
      R_calculate = Rsample1
    }
  }
  if(method == 1){
    R_calculate[ff] = R[ff]
  }
  return(R_calculate)
}
