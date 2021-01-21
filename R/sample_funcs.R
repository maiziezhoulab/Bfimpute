sample_mu_U_Lambda_U <- function(U, X_U, beta_U, lambda_beta_U, mu_0, nu_0, W_0, beta_0){
  D = dim(U)[1]
  N = dim(U)[2]
  F_U = dim(X_U)[1]
  U_2 = U - t(beta_U) %*% X_U
  U_bar = matrix(rowMeans(U_2))
  
  S_bar = matrix(rep(0,D*D),nrow = D)
  for(i in 1:N){
    S_bar = S_bar + U_2[,i] %*% t(U_2[,i])
  }
  S_bar = S_bar/N
  
  mu_0_2 = (beta_0 * mu_0 + N * U_bar) / (beta_0+N)
  beta_0_2 = beta_0+N
  nu_0_2 = nu_0+N+F_U
  
  W_0_2 = pd.solve(pd.solve(W_0) + N * S_bar + beta_0 * (mu_0 %*% t(mu_0)) - beta_0_2 * (mu_0_2 %*% t(mu_0_2)) + lambda_beta_U * (t(beta_U) %*% beta_U))
  
  Lambda_U = rWishart(n = 1, Sigma = W_0_2, df = nu_0_2)[,,1]
  
  L = t(chol(pd.solve(beta_0_2 * Lambda_U)))
  mu_U = L %*% rnorm(D)+mu_0_2
  
  return(list(mu_U,Lambda_U))
}


sample_U0 <- function(U, R, V, X_U, beta_U, Lambda_U, mu_U, alpha){
  D = dim(U)[1]
  N = dim(U)[2]
  U = sapply(1:N, function(i){
    ff = which(R[i,] != 0)                 ####
    V_2 = V[,ff]
    Lambda_U_2 = Lambda_U + alpha * (V_2 %*% t(V_2))
    Lambda_U_2_Inv = pd.solve(Lambda_U_2)
    
    L = t(chol(Lambda_U_2_Inv))

    mu_U_2 = Lambda_U_2_Inv %*% (Lambda_U %*% (mu_U+t(beta_U) %*% X_U[,i]) + alpha * V_2 %*% matrix(R[i,ff]))
    U_i = L %*% rnorm(D)+mu_U_2
    return(U_i)
  })
  return(U)
}


sample_lambda_beta_U <- function(beta_U, Lambda_U, mu, nu){
  # the parameters of gamma distribution may be wrong in matlab Macau
  F_U = dim(beta_U)[1]
  D = dim(beta_U)[2]
  nu_bar = F_U * D + nu
  mu_bar = nu_bar * mu / (nu + mu * sum(diag(t(beta_U) %*% beta_U %*% Lambda_U)))
  lambda_beta_U = rgamma(n = 1, shape = nu_bar / 2, scale = mu_bar * 2 / nu_bar)
  return(lambda_beta_U)
}


sample_beta_U <- function(mu_U, Lambda_U, U, X_U, lambda_beta_U){
  D = dim(U)[1]
  N = dim(U)[2]
  F_U = dim(X_U)[1]
  
  E_1_temp = matrix(rnorm(N * D), nrow = N)
  E_2_temp = matrix(rnorm(F_U * D), nrow = F_U)
  
  L = t(chol(pd.solve(Lambda_U)))
  E_1 = t(L %*% t(E_1_temp))
  E_2 = t(L %*% t(E_2_temp))
  
  U_bar = t(U-matrix(rep(mu_U,N),ncol = N))
  A = X_U %*% t(X_U) + lambda_beta_U * diag(rep(1,F_U))
  B = X_U %*% (U_bar + E_1) + sqrt(lambda_beta_U) * E_2
  
  beta_U = solve(A,B)
  return(beta_U)
}


sample_alpha0 <- function(R,Rsample1,alpha_max){
  ff = (R != 0)                 ####
  numnozero = sum(ff)
  
  # var_all = mean(R(ff).^2,'all')
  var_all = sd(R[ff])^2
  sumsq = sum((Rsample1[ff]-R[ff])^2)
  
  # (a0, b0) correspond to a prior of 1 sample of noise with full variance
  a0 = 0.5
  b0 = 0.5 * var_all
  aN = a0 + numnozero/2
  bN = b0 + sumsq/2
  alpha = rgamma(1, shape = aN, scale = 1/bN)
  if(alpha > alpha_max){
    alpha = alpha_max
  }
  return(alpha)
}




