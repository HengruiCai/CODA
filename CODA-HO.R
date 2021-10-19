# This is an official implementation for CODA with Homogeneous Baseline Covariates (CODA-HO)
# for multiple covariates and mediators.

#Load R Packages
library(foreach)
library(doParallel)
library(policytree)
library(randomForest)
library(MASS)
#install.packages('devtools')
#devtools::install_github("grf-labs/policytree", subdir = "r-package/policytree") #install.packages('policytree')

#simulation
#pre-setting
repnum = 500

info = matrix(0, nrow=repnum, ncol=42)
# Real physical cores in the computer
cores <- detectCores(logical=F)
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl, cores=cores)

info <- foreach(ith=1:repnum, .combine='rbind', .packages=c("policytree", "MASS", "truncnorm", "randomForest")) %dopar%
{ 
  set.seed(2333)
  seeds <- ceiling(runif(repnum, 10000, 1e+09))
  set.seed(seeds[ith])
  cat("Reps", ith, "\n") 
  n_e.list = c(500, 1000) # experimental sample
  n_u = 2000 # different size for the auxillary sample
  K = length(n_e.list)
    
  r = 10 # dimension of X
  s = 2 # dimension of M 
    
  V_hat = V_b_hat = V_e_hat = V_b_e_hat = sigma = sigma_e = rho_hat = sigma_m_hat = W_diff = V_e_best_rule = matrix(0, nrow=1, ncol=K)
  
  for(k in 1:length(n_e.list)){
       
    n_e = n_e.list[k] 
    n = n_u + n_e

    # data generation
    x = matrix(runif(n*r, -2, 2), nrow=n, ncol=r)
    a = rep(0, n)
    for(id in 1:n){
        prob_i = exp(0.4 + 0.2 * x[id, 1] - 0.2 * x[id, 2]) / (1 + exp (0.4 + 0.2 * x[id, 1] - 0.2 * x[id, 2]))
        a[id] = rbinom(1, 1, prob_i)
    }
    m = matrix(0, nrow=n, ncol=s)
    m[, 1] = x[, 1] + 2 * x[, 2] + a * (x[, 1] * x[, 2])
    m[, 2] = x[, 1]^2/2 + 2 * x[, 2] + a * (x[, 1] * x[, 2]) 

    my_rho = 0.7
    my_s1 = 2
    my_s2 = 1.5
    my_mu <- c(0, 0) # Mean
    my_sigma <- matrix(c(my_s1^2, my_s1*my_s2*my_rho, my_s1*my_s2*my_rho, my_s2^2), 2) # Covariance matrix
  
    noises = mvrnorm(n_e * s, mu = my_mu, Sigma = my_sigma)
    e_m = noises[,1] # from MASS package
    e_y = noises[1:n_e,2] 
             
    x_e = x[1:n_e, ]
    a_e = a[1:n_e]
    m_e = m[1:n_e, ] + cbind(e_m[1:n_e],e_m[(n_e+1):(s*n_e)]) 
    y_e = 2 * cos(x_e[, 1]) + x_e[, 2] + 2 * a_e * (x_e[, 2] * x_e[, 1]) + e_y  
      
    df_e = cbind(x_e, a_e, m_e, y_e) # experimental data

    x_u = x[(n_e + 1):n, ]
    a_u = a[(n_e + 1):n]
    m_u = m[(n_e + 1):n, ] + runif(n_u * s, -1, 1)  
    df_u = cbind(x_u, a_u, m_u) # auxillary data
      
    m = rbind(m_e, m_u) 

    #build the regression model on the long-term outcome    
    par_mu_y = lm(y_e ~ x_e + a_e : x_e)$coefficients
    mu_y <- function(rule, df_x){
        as.matrix(cbind(1, df_x, rule * df_x)) %*% as.matrix(par_mu_y)
        } 
      
    #build the regression model on the intermediate outcomes    
    par_mu_m_1 = lm(m[, 1] ~ x + a : x)$coefficients
    mu_m_1 <- function(rule, df_x){
        as.matrix(cbind(1, df_x, rule * df_x)) %*% as.matrix(par_mu_m_1)
        }
    par_mu_m_2 = lm(m[, 2] ~ x + a : x)$coefficients
    mu_m_2 <- function(rule, df_x){
        as.matrix(cbind(1, df_x, rule * df_x)) %*% as.matrix(par_mu_m_2)
        }
      
    ps_fit_e <- glm(a_e ~ x_e, family = binomial)
    pi_e <- predict(ps_fit_e, type = "response")

    ps_fit_u <- glm(a_u ~ x_u, family = binomial)
    pi_u <- predict(ps_fit_u, type = "response")

    #value functions
    V_e_b <- function(rule){ 
        (a_e == rule) * (y_e - mu_y(rule, x_e)) / (pi_e * a_e + (1 - pi_e) * (1 - a_e)) + mu_y(rule, x_e) 
        }
    W_e_b_1 <- function(rule){ 
        (a_e == rule) * (m_e[, 1] - mu_m_1(rule, x_e)) / (pi_e * a_e + (1 - pi_e) * (1 - a_e)) + mu_m_1(rule, x_e)
        }
    W_u_b_1 <- function(rule){
        (a_u == rule) * (m_u[, 1] - mu_m_1(rule, x_u)) / (pi_u * a_u + (1 - pi_u) * (1 - a_u)) + mu_m_1(rule, x_u)
        }
    W_e_b_2 <- function(rule){ 
        (a_e == rule) * (m_e[, 2] - mu_m_2(rule, x_e)) / (pi_e * a_e + (1 - pi_e) * (1 - a_e)) + mu_m_2(rule, x_e)
        }
    W_u_b_2 <- function(rule){
        (a_u == rule) * (m_u[, 2] - mu_m_2(rule, x_u)) / (pi_u * a_u + (1 - pi_u) * (1 - a_u)) + mu_m_2(rule, x_u)
        }

    # Calibrated ODR
    sample_ratio = n_e / n_u

    Gamma_V_e = Gamma_W_e_1 = Gamma_W_e_2 = matrix(0, nrow=n_e, ncol=2)
    Gamma_W_u_1 = Gamma_W_u_2 = matrix(0, nrow=n_u, ncol=2)
    Gamma_V_e[,1] = V_e_b(0)
    Gamma_V_e[,2] = V_e_b(1) 
    Gamma_W_e_1[,1] = W_e_b_1(0)
    Gamma_W_e_1[,2] = W_e_b_1(1) 
    Gamma_W_u_1[,1] = W_u_b_1(0)
    Gamma_W_u_1[,2] = W_u_b_1(1) 
    Gamma_W_e_2[,1] = W_e_b_2(0)
    Gamma_W_e_2[,2] = W_e_b_2(1) 
    Gamma_W_u_2[,1] = W_u_b_2(0)
    Gamma_W_u_2[,2] = W_u_b_2(1)       
      
    #value functions
    V_e <- function(b){
        rule = predict(b, x_e) 
        mean((rule == 1) * Gamma_V_e[, 1] + (rule == 2) * Gamma_V_e[, 2])
        }
    W_e_1 <- function(b){
        rule = predict(b, x_e) 
        mean((rule == 1) * Gamma_W_e_1[, 1] + (rule == 2) * Gamma_W_e_1[, 2])
        }
    W_u_1 <- function(b){
        rule = predict(b, x_u) 
        mean((rule == 1) * Gamma_W_u_1[, 1] + (rule == 2) * Gamma_W_u_1[, 2])
        }
    W_e_2 <- function(b){
        rule = predict(b, x_e) 
        mean((rule == 1) * Gamma_W_e_2[, 1] + (rule == 2) * Gamma_W_e_2[, 2])
        }
    W_u_2 <- function(b){
        rule = predict(b, x_u) 
        mean((rule == 1) * Gamma_W_u_2[, 1] + (rule == 2) * Gamma_W_u_2[, 2])
        }
      
    #value functions
    V_e_i <- function(b){
        rule = predict(b, x_e)
        (rule == 1) * Gamma_V_e[, 1] + (rule == 2) * Gamma_V_e[, 2]
        }
    W_e_i_1 <- function(b){
        rule = predict(b, x_e)
        (rule == 1) * Gamma_W_e_1[, 1] + (rule == 2) * Gamma_W_e_1[, 2]
        }
    W_u_i_1 <- function(b){
        rule = predict(b, x_u)
        (rule == 1) * Gamma_W_u_1[, 1] + (rule == 2) * Gamma_W_u_1[, 2]
        }

    W_e_i_2 <- function(b){
        rule = predict(b, x_e)
        (rule == 1) * Gamma_W_e_2[, 1] + (rule == 2) * Gamma_W_e_2[, 2]
        }
    W_u_i_2 <- function(b){
        rule = predict(b, x_u)
        (rule == 1) * Gamma_W_u_2[, 1] + (rule == 2) * Gamma_W_u_2[, 2]
        }

      
    rho <- function(b){ 
        c(mean((V_e_i(b) - V_e(b)) * (W_e_i_1(b) - W_e_1(b))), 
          mean((V_e_i(b) - V_e(b)) * (W_e_i_2(b) - W_e_2(b))))
        }

    Sigma_m <- function(b){
        S_m = matrix(0, nrow=s, ncol=s)
        S_m[1,1] = mean((W_e_i_1(b) - W_e_1(b)) * (W_e_i_1(b) - W_e_1(b))) + sample_ratio * mean((W_u_i_1(b) - W_u_1(b)) * (W_u_i_1(b) - W_u_1(b)))
        S_m[1,2] = mean((W_e_i_1(b) - W_e_1(b)) * (W_e_i_2(b) - W_e_2(b))) + sample_ratio * mean((W_u_i_1(b) - W_u_1(b)) * (W_u_i_2(b) - W_u_2(b)))
        S_m[2,1] = mean((W_e_i_2(b) - W_e_2(b)) * (W_e_i_1(b) - W_e_1(b))) + sample_ratio * mean((W_u_i_2(b) - W_u_2(b)) * (W_u_i_1(b) - W_u_1(b)))
        S_m[2,2] = mean((W_e_i_2(b) - W_e_2(b)) * (W_e_i_2(b) - W_e_2(b))) + sample_ratio * mean((W_u_i_2(b) - W_u_2(b)) * (W_u_i_2(b) - W_u_2(b)))
        S_m
        }
      

    #find the ODR with experimental sample only
    tree_e <- policy_tree(x_e, Gamma_V_e, depth = 2)
    V_e_hat[k] = V_e(tree_e) 
      
    # value under true
    N = 5000
    X_N = matrix(runif(N*r, -2 , 2), nrow=N, ncol=r)
    V_b_e_hat[k] = mean(2 * cos(X_N[, 1]) + X_N[, 2] + (predict(tree_e, X_N) - 1) * 2 * (X_N[, 2] * X_N[, 1])) 
      
    # variance for V_e under regular ODR
    sigma_e[k] = sd(V_e_i(tree_e)) / sqrt(n_e)   
      
    #find CODA with two samples 
    Gamma_V = matrix(0, nrow=n_e, ncol=2)
    Gamma_V[,1] = Gamma_V_e[,1] - (rho(tree_e) %*% solve(Sigma_m(tree_e))) %*% t(as.matrix(cbind(Gamma_W_e_1[,1] - W_u_1(tree_e), Gamma_W_e_2[,1] - W_u_2(tree_e))))
    Gamma_V[,2] = Gamma_V_e[,2] - (rho(tree_e) %*% solve(Sigma_m(tree_e))) %*% t(as.matrix(cbind(Gamma_W_e_1[,2] - W_u_1(tree_e), Gamma_W_e_2[,2] - W_u_2(tree_e)))) 
      
    best_tree <- policy_tree(x_e, Gamma_V, depth = 2)

    best_rule = predict(best_tree, x_e)
      
    Gamma_V[,1] = Gamma_V_e[,1] - (rho(best_tree) %*% solve(Sigma_m(best_tree))) %*% t(as.matrix(cbind(Gamma_W_e_1[,1] - W_u_1(best_tree), Gamma_W_e_2[,1] - W_u_2(best_tree))))
    Gamma_V[,2] = Gamma_V_e[,2] - (rho(best_tree) %*% solve(Sigma_m(best_tree))) %*% t(as.matrix(cbind(Gamma_W_e_1[,2] - W_u_1(best_tree), Gamma_W_e_2[,2] - W_u_2(best_tree))))
      
    V_hat[k] = mean((best_rule == 1) * Gamma_V[, 1] + (best_rule == 2) * Gamma_V[, 2])
 
    # variance for V under CODA
    sigma[k] = sqrt(var(V_e_i(best_tree)) - (rho(best_tree) %*% solve(Sigma_m(best_tree))) %*% rho(best_tree)) / sqrt(n_e)
      
    # value under true
    V_b_hat[k] = mean(2 * cos(X_N[, 1]) + X_N[, 2] + (predict(best_tree, X_N) - 1) * 2 * (X_N[, 2] * X_N[, 1])) 
      
  }
  cat("Repe",ith,"\n") 
  c(V_hat,V_b_hat,V_e_hat,V_b_e_hat,sigma,sigma_e)
      
}
stopImplicitCluster()
stopCluster(cl) 

save(info,file="Res.RData")
