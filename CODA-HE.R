# This is an official implementation for CODA with Heterogeneous Baseline Covariates (CODA-HE).

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
    
  r = 2 # dimension of X
  s = 1 # dimension of M 
    
  V_hat = V_b_hat = V_e_hat = V_b_e_hat = sigma = sigma_e = rho_hat = sigma_m_hat = W_diff_hat = V_e_best_rule = matrix(0, nrow=1, ncol=K)
    
  for(k in 1:length(n_e.list)){     
       
    n_e = n_e.list[k] 
    n = n_u + n_e

    # data generation 
    x_e = matrix(runif(n_e*r, -2, 2), nrow=n_e, ncol=r)
    x_u = matrix(runif(n_u*r, -1, 1.5), nrow=n_u, ncol=r) 
    x = rbind(x_e, x_u)

    a_e = rep(0, n_e)
    for(id in 1:n_e){
        prob_i = exp(0.4 + 0.2 * x_e[id, 1] - 0.2 * x_e[id, 2]) / (1 + exp (0.4 + 0.2 * x_e[id, 1] - 0.2 * x_e[id, 2]))
        a_e[id] = rbinom(1, 1, prob_i)
    }
    a_u = rep(0, n_u)
    for(id in 1:n_u){
        prob_i = exp(0.4 + 0.2 * x_u[id, 1] - 0.2 * x_u[id, 2]) / (1 + exp (0.4 + 0.2 * x_u[id, 1] - 0.2 * x_u[id, 2]))
        a_u[id] = rbinom(1, 1, prob_i)
    }   
    a = c(a_e, a_u)

    my_rho = 0.7
    my_s1 = 2
    my_s2 = 1.5
    my_mu <- c(0, 0) # Mean
    my_sigma <- matrix(c(my_s1^2, my_s1*my_s2*my_rho, my_s1*my_s2*my_rho, my_s2^2), 2) # Covariance matrix
  
    noises = mvrnorm(n_e, mu = my_mu, Sigma = my_sigma)
    e_m = noises[,1] # from MASS package
    e_y = noises[,2] 
      
    m_e = x_e[, 1] + 2 * x_e[, 2] + a_e * (x_e[, 1] - x_e[, 2]) + e_m 
    y_e = 2 * x_e[, 1] + x_e[, 2] + 2 * a_e * (x_e[, 2] - x_e[, 1]) + e_y 
    df_e = cbind(x_e, a_e, m_e, y_e) # experimental data
 
    m_u = x_u[, 1] + 2 * x_u[, 2] + a_u * (x_u[, 1] - x_u[, 2]) + runif(n_u, -1, 1) 
    df_u = cbind(x_u, a_u, m_u) # auxillary data

    m = c(m_e, m_u)

    #build the regression model on the long-term outcome    
    par_mu_y = lm(y_e ~ x_e + a_e : x_e)$coefficients
    mu_y <- function(rule, df_x){
        as.matrix(cbind(1, df_x, rule * df_x)) %*% as.matrix(par_mu_y)
        } 
      
    #build the regression model on the intermediate outcomes    
    par_mu_m = lm(c(m_e, m_u) ~ x + a : x)$coefficients
    mu_m <- function(rule, df_x){
        as.matrix(cbind(1, df_x, rule * df_x)) %*% as.matrix(par_mu_m)
        }

    ps_fit <- glm(a ~ x, family = binomial)
    pi <- predict(ps_fit, type = "response")

    ps_fit_e <- glm(a_e ~ x_e, family = binomial)
    pi_e <- predict(ps_fit_e, type = "response")
      
    sampling_fit_data = as.matrix(x) 
    sample_index = c(rep(1,n_e), rep(0,n_u))
    colnames(sampling_fit_data) = NULL
    rf = randomForest(sampling_fit_data, as.factor(sample_index)) 

    #value functions 
    V_e_b <- function(rule){  
        (a_e == rule) * (y_e - mu_y(rule, x_e)) / (pi_e * a_e + (1 - pi_e) * (1 - a_e)) + mu_y(rule, x_e) 
        }
      
    W_diff_b <- function(rule){ 
        sampling_new_data = as.matrix(x)   
        colnames(sampling_new_data) = NULL 
        sampling_prob = predict(rf, sampling_new_data, type ='prob')[,2] 
        c(((sample_index / sampling_prob) * (a == rule) * (m - mu_m(rule, x)) / (pi * a + (1 - pi) * (1 - a)))[1:n_e], -(((1 - sample_index) / (1 - sampling_prob)) * (a == rule) * (m - mu_m(rule, x)) / (pi * a + (1 - pi) * (1 - a)))[(1+n_e):n])
        } 
      
      
    Gamma_V_e = matrix(0, nrow=n_e, ncol=2) 
    Gamma_V_e[,1] = V_e_b(0)
    Gamma_V_e[,2] = V_e_b(1) 
    
    Gamma_W_diff = matrix(0, nrow=n, ncol=2) 
    Gamma_W_diff[,1] = W_diff_b(0)
    Gamma_W_diff[,2] = W_diff_b(1)

    #value functions
    V_e <- function(b){
        rule = predict(b, x_e)  
        mean((rule == 1) * Gamma_V_e[, 1] + (rule == 2) * Gamma_V_e[, 2])
        } 

    W_diff <- function(b){
        rule = predict(b, x)  
        mean((rule == 1) * Gamma_W_diff[, 1] + (rule == 2) * Gamma_W_diff[, 2])
        }

    #value functions
    V_e_i <- function(b){
        rule = predict(b, x_e) 
        (rule == 1) * Gamma_V_e[, 1] + (rule == 2) * Gamma_V_e[, 2]
        } 
    W_diff_i <- function(b){
        rule = predict(b, x) 
        (rule == 1) * Gamma_W_diff[, 1] + (rule == 2) * Gamma_W_diff[, 2]
        }

    rho <- function(b){
        sum((V_e_i(b) - V_e(b)) * W_diff_i(b)[1:n_e]) / sqrt(n * n_e)
        }

    Sigma_m <- function(b){
        mean((W_diff_i(b) - W_diff(b)) * (W_diff_i(b) - W_diff(b))) 
    }

    #find the ODR with experimental sample only
    tree_e <- policy_tree(x_e, Gamma_V_e, depth = 2)
    V_e_hat[k] = V_e(tree_e) 
      
    # value under true
    N = 5000
    X_N = matrix(runif(N*r, -2 , 2), nrow=N, ncol=r)
    V_b_e_hat[k] = mean(2 * X_N[, 1] + X_N[, 2] + (predict(tree_e, X_N) - 1) * 2 * (X_N[, 2] - X_N[, 1])) 

    # variance for V_e under regular ODR
    sigma_e[k] = sd(V_e_i(tree_e)) / sqrt(n_e)  
       
    #find the CODA with two samples 
    Gamma_V = Gamma_V_e - sqrt(n/n_e) * ((n_e/n) * Gamma_W_diff[1:n_e,] + sum(W_diff_i(tree_e)[(1+n_e):n])/n) * (rho(tree_e) / (Sigma_m(tree_e)))
      
    best_tree <- policy_tree(x_e, Gamma_V, depth = 2)

    best_rule = predict(best_tree, x_e)
    
    Gamma_V = Gamma_V_e - sqrt(n/n_e) * ((n_e/n) * Gamma_W_diff[1:n_e,] + sum(W_diff_i(best_tree)[(1+n_e):n])/n) * (rho(best_tree) / (Sigma_m(best_tree)))

    V_hat[k] = mean((best_rule == 1) * Gamma_V[, 1] + (best_rule == 2) * Gamma_V[, 2])
 
    # variance for V under CODA
    sigma[k] = sqrt(var(V_e_i(best_tree)) - (rho(best_tree) %*% solve(Sigma_m(best_tree))) %*% rho(best_tree)) / sqrt(n_e)
      
    # value under true
    V_b_hat[k] = mean(2 * X_N[, 1] + X_N[, 2] + (predict(best_tree, X_N) - 1) * 2 * (X_N[, 2] - X_N[, 1])) 
     
    rho_hat[k] = rho(best_tree)
    sigma_m_hat[k] = Sigma_m(best_tree)
    W_diff_hat[k] = W_diff(best_tree)
    V_e_best_rule[k] = V_e(best_tree)
      
  }
  cat("Repe",ith,"\n") 
  c(V_hat,V_b_hat,V_e_hat,V_b_e_hat,sigma,sigma_e,rho_hat,sigma_m_hat, W_diff_hat,V_e_best_rule)
      
}
stopImplicitCluster()
stopCluster(cl) 

save(info,file="Res.RData")
