
# g(): the competing risks function used in the 'incidence' part of the mixture model
# Input:
#   - l1: the rate parameter for progression to CIN2+
#   - l2: the rate parameter for viral clearance
#   - t: time
# Output:
#   - the cumulative risk of CIN2+ among HPV-positive women wihtout CIN2+ at baseline
g <- function(l1, l2, t) {
  # calculate a, b from l1, l2 at the start of the g function
  a <- l1 / (l1 + l2)
  b <- l1 + l2
  return(ifelse(t == Inf, 1,  a * (1 - exp(-b * t))))
}

# create.par(): function to create the parameter names which will be used in the mstep expressions
# Input:
#   - n1: number of covariates used for the progression rate parameter + 1
#   - n2: number of covariates used for the clearance rate parameter + 1
#   - n3: number of covariates used in the parameter for the probability of prevalent disease + 1
# Output:
#   - vector of names of the parameters with (1) g0, g1... for l1, (2) w0, w1... for l2 and (3) p0, p1... for pi
create.par <- function(n1, n2, n3) {
  par1 <- paste0(c("g"), seq(0, n1 - 1))
  par2 <- paste0(c("w"), seq(0, n2 - 1))
  par3 <- paste0(c("p"), seq(0, n3 - 1))
  return(c(par1, par2, par3))
}

# create.covariate.data(): creates matrices for each parameter of observed independent variables (like the 'X' matrices in a simple regression model)
# Input:
#   - l1_x: vector containing names of covariates used in the lambda_1 (progression rate) parameter (must be column names in data input)
#   - l2_x: vector containing names of covariates used in the lambda_2 (clearance rate) parameter (must be column names in data input)
#   - pi_x: vector containing names of covariates used in the pi parameter (must be column names in data input)
#   - data: first three columns must be (1) left interval, (2) right interval and (3) z indicator for prevalent/incident disease,
#           following columns must match the column names given in l1_x, l2_x and pi_x
# Output:
#   - list of 3 matrices corresponding to relevant covariates for each of the parameters lambda_1, lambda_2 and pi
create.covariate.data <- function(l1_x, l2_x, pi_x, data){
  n1 <- length(l1_x) + 1 # number of parameters for lambda_1 (add 1 for intercept)
  n2 <- length(l2_x) + 1 # number of parameters for lambda_2 (add 1 for intercept)
  n3 <- length(pi_x) + 1 # number of parameters for pi (add 1 for intercept)

  data1 <- matrix(c(rep(1, dim(data)[1]*n1)), ncol=n1) # create empty matrix
  if (n1 != 1){
    for (i in 1:(n1-1)){
      data1[,(i+1)] <- data[[l1_x[i]]] # use relevant column from the main data
    }
  }
  data2 <- matrix(c(rep(1, dim(data)[1]*n2)), ncol=n2) # create empty matrix
  if (n2 != 1){
    for (i in 1:(n2-1)){
      data2[,(i+1)] <- data[[l2_x[i]]] # use relevant column from the main data
    }
  }
  data3 <- matrix(c(rep(1, dim(data)[1]*n3)), ncol=n3) # create empty matrix
  if (n3 != 1){
    for (i in 1:(n3-1)){
      data3[,(i+1)] <- data[[pi_x[i]]] # use relevant column from the main data
    }
  }
  return(list(data1, data2, data3))
}

estep.h <- function(current_par, n1, n2, n3, data1, data2, data3, data_int){
  right <- data_int[[2]] # right intervals from the data
  z <- data_int[[3]] # indicator variable for prevalent or incident disease

  h <- current_par[1]
  current_par <- current_par[-1]
  # multiply data by corresponding parameters to calculate l1, l2 and pi
  l1 <- exp(data1 %*% current_par[1:(n1)])
  l2 <- exp(data2 %*% current_par[(n1+1):(n1+n2)])
  if (n3 == 1){
    # if pi has no covariates then the value of p is the same for all women
    p <- data3 %*% current_par[n1+n2+n3]
  }else if (n3 > 1){
    p <- exp(data3 %*% current_par[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% current_par[(n1+n2+1):(n1+n2+n3)]))
  }

  # the expected value of z only depends on the right interval since the left is zero when z is unknown
  est <- p/(p + (1 - p)*(1-exp(-exp(h)*right) + exp(-exp(h)*right)*g(l1, l2, right)))
  est[which(right==Inf)] <- p[which(right==Inf)]
  # the expected value of the latent variable is the estimate above if z is unknown
  # otherwise it is the true value of z
  new_z <- ifelse(is.na(z), est, z)
  return(new_z)
}

build.q.function.h <- function(n1, n2, n3, model_par, data1, data2, data3){
  # function to 'build' the expected log-likelihood function as an expression which
  # allows it to be differentiated w.r.t. different parameters in the mstep() function

  # we loop through the names of covariate parameters (e.g., g0, g1, w0, p0, p1, p2) for each of l1, l2 and pi
  # and paste them together with the corresponding data matrix, this string is then used to create the
  # expression which will be differentiated w.r.t each parameter in the mstep() function
  model_par <- model_par[-1]
  # expression for lambda_1 (progression rate parameter)
  expr1 <- model_par[1]
  if (n1 != 1){
    for (i in 1:(n1-1)){
      expr1 <- paste0(expr1, "+", model_par[1+i], "*", "data1", i)
    }}

  # expression for lambda_2 (clearance rate parameter)
  expr2 <- model_par[n1+1]
  if (n2!=1){
    for (i in 1:(n2-1)){
      expr2 <- paste0(expr2, "+", model_par[n1+1+i], "*", "data2", i)
    }}

  # expression for pi (prevalent probability parameter)
  expr3 <- model_par[n1+n2+1]
  if (n3 > 1){
    for (i in 1:(n3-1)){
      expr3 <- paste0(expr3, "+", model_par[n1+n2+1+i], "*", "data3", i)
    }

    #  expected complete log-likelihood for Ri!=Inf
    q1 <- paste0("z*(", expr3, ") - log(1 + exp(",expr3, ")) + (1-z)*log(((1-exp(-exp(h)*right)) + exp(-exp(h)*right)*(exp(",
                 expr1, ")/(exp(", expr1, ") + exp(", expr2, ")))*(1-exp(-(exp(", expr1, ") + exp(", expr2,
                 "))*right))) - ((1-exp(-exp(h)*left)) + exp(-exp(h)*left)*(exp(", expr1, ")/(exp(", expr1, ") + exp(",
                 expr2,")))*(1-exp(-(exp(", expr1, ") + exp(", expr2, "))*left))))" )


    # expected log-likelihood for Ri==Inf because g(Inf)=1
    q2 <- paste0("z*(", expr3, ") - log(1 + exp(",expr3, ")) + (1-z)*log(1 - (1-exp(-exp(h)*left) + exp(-exp(h)*left)*(exp(",
                 expr1, ")/(exp(", expr1, ") + exp(", expr2,")))*(1-exp(-(exp(", expr1, ") + exp(", expr2, "))*left))))" )


  }else if (n3 == 1){ # if there are no covariates for pi we use a different expression (don't need to use logit function)
    #  expected complete log-likelihood for Ri!=Inf for when there is only intercept for pi parameter
    q1 <- paste0("z*log(", expr3, ") + (1-z)*(log(1-", expr3, ") + log((1-exp(-exp(h)*right)) + exp(-exp(h)*right)*(exp(",
                 expr1, ")/(exp(", expr1, ") + exp(", expr2, ")))*(1-exp(-(exp(", expr1, ") + exp(", expr2,
                 "))*right)) - ((1-exp(-exp(h)*left)) + exp(-exp(h)*left)*(exp(", expr1, ")/(exp(", expr1, ") + exp(",
                 expr2,")))*(1-exp(-(exp(", expr1, ") + exp(", expr2, "))*left)))))" )


    # expected log-likelihood for Ri==Inf (because g(Inf)=1) for when there is only intercept for pi parameter
    q2 <- paste0("z*log(", expr3, ") + (1-z)*(log(1-", expr3, ") + log(1 - (1-exp(-exp(h)*left) + exp(-exp(h)*left)*(exp(",
                 expr1, ")/(exp(", expr1, ") + exp(", expr2,")))*(1-exp(-(exp(", expr1, ") + exp(", expr2, "))*left)))))" )
  }

  # convert from string to expression to be used in the differentiate function in the mstep()
  q1 <- str2expression(q1)
  q2 <- str2expression(q2)
  return(list(q1, q2))
}

mstep.h <- function(current_par, expected_z, n1, n2, n3, data1, data2, data3, data_int){
  left <- data_int[[1]] # left intervals from the data
  right <- data_int[[2]] # right intervals from the data
  z <- expected_z # expected value of z = output from the estep.h() function
  pars <- c('h', create.par(n1, n2, n3)) # creates a vector of the names of parameters that need to be estimated

  #current_par[1] <- log(current_par[1])
  # create temp data for each covariate used for each parameter (necessary for calculating the hessian)
  # each covariate as a vector rather than a matrix or relevant covariates for each parameter because
  # the differentiation function cannot differentiate indexed data frames (e.g., data[,1])
  if (n1 !=1){
    for (i in 1:(n1-1)){
      assign(paste0("data1", i), data1[,i+1], envir=environment())
    }}
  if (n2 !=1){
    for (i in 1:(n2-1)){
      assign(paste0("data2", i), data2[,i+1], envir=environment())
    }}
  if (n3!=1){
    for (i in 1:(n3-1)){
      assign(paste0("data3", i), data3[,i+1], envir=environment())
    }}

  # assign the values of current_par to vectors with names from the 'pars' list
  for (i in 1:length(pars)){
    assign(pars[i], current_par[i], envir=environment())
  }
  # expected log-likelihoods:
  exp_log_lik <- build.q.function.h(n1, n2, n3, pars, data1, data2, data3)
  q1 <- exp_log_lik[[1]] # q1 is the expected log-likelihood for when Right!=Inf
  q2 <- exp_log_lik[[2]] # q2 is the expected log-likelihood for when Right=Inf, so only depends on left interval

  if (n3 > 1){
    hess <- matrix(rep(0, (n1+n2+n3+1)^2), nrow=(n1+n2+n3+1)) # empty hessian matrix
    grad <- rep(0, (n1+n2+n3+1)) # empty gradient matrix
  }else if(n3 == 1){
    hess <- matrix(rep(0, (n1+n2+1)^2), nrow=(n1+n2+1)) # empty hessian matrix
    grad <- rep(0, (n1+n2+1)) # empty gradient matrix
  }
  # parameters for l1 and l2 depend on left/right so they need a separate Q expression when R=Inf
  # we can't have both right=0 and left=0 because then the code will try evaluate log(exp(0) - exp(0)) = log(0) = -Inf
  # we also exclude cases with right = Inf and left = 0 because they aren't informative and also make an error
  for (i in 1:(n1+n2 + 1)){
    # we loop through the parameter list and calculate the first derivatives for gradient vector for the l1 and l2 parameters
    u1 <- sum(eval(D(q1, pars[i]))[which(right<Inf & right>0)])
    u2 <- sum(eval(D(q2, pars[i]))[which(right==Inf & left!=0)])
    grad[i] <- u1 + u2
    for (j in 1:i){
      # calculate second derivatives for the Hessian matrix for the l1 and l2 parameters
      h1 <- sum(eval(D(D(q1, pars[i]), pars[j]))[which(right<Inf & right>0)])
      h2 <- sum(eval(D(D(q2, pars[i]), pars[j]))[which(right==Inf & left!=0)])
      # the hessian is a symmetric matrix so we can fill both hess[i,j] and hess[j,i] together to save time
      hess[i,j] <- hess[j,i] <- h1 + h2
    }
  }
  if (n3 > 1){
    # parameters for pi (for example p0, p1, p2, p3, p4) do not depend on left/right so we can use q1 only
    for (i in (n1+n2+2):(n1+n2+n3+1)){
      # calculate first derivatives for gradient vector for the pi parameters
      grad[i] <- sum(eval(D(q1, pars[i])))
      for (j in (n1+n2+2):i){
        # calculate second derivatives for Hessian matrix for the pi parameters
        # the hessian is a symmetric matrix so we can fill both hess[i,j] and hess[j,i] together to save time
        hess[i,j] <- hess[j,i] <- sum(eval(D(D(q1, pars[i]), pars[j])))
      }
    }

    new_par <-  current_par - solve(hess)%*%grad
  }else if (n3 == 1) {
    new_pi <- mean(z)
    new_par <- current_par[1:(n1+n2+1)] - solve(hess)%*%grad
    new_par <- c(new_par, new_pi)
  }
  #new_par[1] <- exp(new_par[1])
  # single Newton step to calculate updated parameters
  return(list(as.vector(new_par), hess)) # return the new parameters
}

# log.likelihood(): function for the log-likelihood
# Input:
#   - current_par: the current value of the parameters
#   - n1: number of covariates used for the progression rate parameter + 1
#   - n2: number of covariates used for the clearance rate parameter + 1
#   - n3: number of covariates used in the parameter for the probability of prevalent disease + 1
#   - data1: matrix with columns corresponding to the covariates used in the progression rate parameter
#   - data2: matrix with columns corresponding to the covariates used in the clearance rate parameter
#   - data3: matrix with columns corresponding to the covariates used in the prevalent probability parameter
#   - data_int: data frame with columns (1) left, (2) right, and (3) z, representing the left and right interval
#               and the prevalent/incident indicator variable
#   -h: value of background risk for all women
# Output:
#   - value of the log-likelihood for the current parameter values
log.likelihood.h <- function(current_par,  data1, data2, data3, data_int){
  left <- data_int[[1]] # left intervals from the data
  right <- data_int[[2]] # right intervals from the data
  z <- data_int[[3]]

  h <- exp(current_par[1])
  current_par <- current_par[-1]

  n1 <- ncol(data1)
  n2 <- ncol(data2)
  n3 <- ncol(data3)

  l1 <- exp(data1 %*% current_par[1:n1])
  l2 <- exp(data2 %*% current_par[(n1+1):(n1+n2)])
  if (n3 > 1) {
    p <- exp(data3 %*% current_par[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% current_par[(n1+n2+1):(n1+n2+n3)]))
  }else if (n3 == 1) {
    p <- data3 %*% current_par[n1+n2+n3]
  }

  llk <- rep(0, length(right))
  llk[which(right==0)] <- I(log(p))[which(right==0)] # (same as z=1)
  llk[which(z==0 & right<Inf)] <- I(log((1-p)*((1- exp(-h*right) + exp(-h*right)*g(l1, l2, right)) -
                                                 (1-exp(-h*left) + exp(-h*left)*g(l1, l2, left)))))[which(z==0 & right <Inf)]
  llk[which(z==0 & right==Inf)] <- I(log((1-p)*(1 - (1- exp(-h*left) + exp(-h*left)*g(l1, l2, left)))))[which(z==0 & right==Inf)]

  llk[which(is.na(z) & right<Inf)] <- I(log(p + (1-p)*(1 - exp(-h*right) + exp(-h*right)*g(l1, l2, right))))[which(is.na(z) & right<Inf)]
  llk[which(is.na(z) & right==Inf)] <- 0
  return(sum(llk))
}

# em_function(): combining the E- and M-step and repeating until the parameter estimates converge
# Input:
#   - current_theta: initial value for the EM algorithm
#   - l1_x: vector containing names of covariates used in the lambda_1 (progression rate) parameter (must be column names in data input)
#   - l2_x: vector containing names of covariates used in the lambda_2 (clearance rate) parameter (must be column names in data input)
#   - pi_x: vector containing names of covariates used in the pi parameter (must be column names in data input)
#   - data: first three columns must be (1) left interval, (2) right interval and (3) z indicator for prevalent/incident disease,
#           following columns must match the column names given in l1_x, l2_x and pi_x
#   - h: background risk value
# Output:
#   - inital.values: inital values, usualy calculated with the find.init() function
#   - theta.hat: optimum parameter values estimated by the EM algorithm
#   - num.iterations: number of iterations until algorithm converged
#   - log.likelihood: value of log.likelihood at optimum parameter values
#   - hess: hessian matrix (can be used to calculated parameter variances)
#   - summary: data frame with estimate, std.dev, and 95% CI for each parameter (to be used in data set comparisons)
em.function.h <- function(init, l1_x, l2_x, pi_x, data, epsilon = 1e-06, silent=T){

  n1 <- length(l1_x) + 1
  n2 <- length(l2_x) + 1
  n3 <- length(pi_x) + 1

  covariate_data <- create.covariate.data(l1_x, l2_x, pi_x, data)
  data1 <- covariate_data[[1]]
  data2 <- covariate_data[[2]]
  data3 <- covariate_data[[3]]
  current_theta <- init
  old_llk <- 0
  new_theta <- as.vector(mstep.h(current_theta, estep.h(current_theta, n1, n2, n3, data1, data2, data3, data),
                                 n1, n2, n3, data1, data2, data3, data)[[1]])
  new_llk <- log.likelihood.h(new_theta, data1, data2, data3, data)
  if (!silent) print(round(c(new_theta, new_llk), 4))
  iter <- 1
  while (abs(new_llk - old_llk) > epsilon){
    if (iter>=1) {
      current_theta <- new_theta
      old_llk <- new_llk}
    new_theta <- as.vector(mstep.h(current_theta, estep.h(current_theta, n1, n2, n3, data1, data2, data3, data),
                                   n1, n2, n3, data1, data2, data3, data)[[1]])
    new_llk <- log.likelihood.h(new_theta, data1, data2, data3, data)
    if (!silent) print(round(c(new_theta, new_llk), 4))
    iter <- iter + 1
  }
  hess <- mstep.h(new_theta, estep.h(new_theta, n1, n2, n3, data1, data2, data3, data),
                  n1, n2, n3, data1, data2, data3, data)[[2]]
  names(new_theta) <- c('h', create.par(n1, n2, n3))
  std.dev <- round(sqrt(-diag(solve(hess))),4)

  if (n3 > 1){
    temp_theta <- new_theta
  }else if (n3==1){
    temp_theta <- new_theta[-length(new_theta)] # remove the last one cos we dont have a std. dev estimate for pi
  }

  summary <- round(data.frame(theta.hat = temp_theta, std.dev, lower = temp_theta - 1.96*std.dev,
                              upper = temp_theta + 1.96*std.dev),4)
  rownames(summary) <- names(temp_theta)
  return(list(initial.values = init, theta.hat = new_theta, num.iterations=iter,
              log.likelihood = new_llk, hess=hess, std.dev=std.dev, summary=summary))
}

short.em.h <- function(l1_x, l2_x, pi_x, data, short.epsilon=1e-1, silent=T){
  n1 <- length(l1_x) + 1
  n2 <- length(l2_x) + 1
  n3 <- length(pi_x) + 1

  covariate_data <- create.covariate.data(l1_x, l2_x, pi_x, data)
  data1 <- covariate_data[[1]]
  data2 <- covariate_data[[2]]
  data3 <- covariate_data[[3]]

  current_theta <- log(runif(n1+n2+n3))
  h.init <- log(runif(1, 0, 0.2))
  current_theta <- c(h.init, current_theta)
  current_theta[n1+n2+n3+1] <- ifelse(n3==1, exp(current_theta[n1+n2+n3+1]), current_theta[n1+n2+n3+1])

  old_llk <- 0
  new_theta <- as.vector(mstep.h(current_theta, estep.h(current_theta, n1, n2, n3, data1, data2, data3, data),
                                 n1, n2, n3, data1, data2, data3, data)[[1]])
  new_llk <- log.likelihood.h(new_theta, data1, data2, data3, data)

  if (abs(new_llk) == Inf | is.nan(new_llk) | is.na(new_llk)){
    return(short.em.h(l1_x, l2_x, pi_x, data, short.epsilon))
  }
  while (abs(new_llk - old_llk) > short.epsilon){
    current_theta <- new_theta
    old_llk <- new_llk
    new_theta <- try(as.vector(mstep.h(current_theta, estep.h(current_theta, n1, n2, n3, data1, data2, data3, data),
                                       n1, n2, n3, data1, data2, data3, data)[[1]]), silent=T)

    if (class(new_theta) == 'try-error'){
      return(short.em.h(l1_x, l2_x, pi_x, data, short.epsilon))
    }
    new_llk <- log.likelihood.h(new_theta, data1, data2, data3, data)
    if (! silent) print(round(c(new_theta, new_llk), 4))

    if (abs(new_llk) == Inf | is.nan(new_llk) | is.na(new_llk)){
      return(short.em.h(l1_x, l2_x, pi_x, data, short.epsilon))
    }
  }
  return(list(theta.hat = new_theta, log.likelihood = new_llk))
}


#' Simulator
#'
#' Cervical cancer screening data simulator.
#'
#' @import Rlab
#'
#' @param n Number of women in the simulated data set.
#' @param l1_x A vector containing the names of covariates used in the \ifelse{html}{\out{\eqn{\lambda}<sub>1</sub>}}{ \eqn{\lambda_1}} (progression rate) parameter. Options are "age", "HPV16" and "cytology".
#' @param l2_x A vector containing the names of covariates used in the \eqn{\lambda_2} (clearance rate) parameter. Options are "age", "HPV16" and "cytology".
#' @param pi_x A vector containing the names of covariates used in the \eqn{\pi} parameter (probability of prevalent disease). Options are "age", "HPV16" and "cytology".
#' @param params Numerical vector the parameter values to be used in the data simulation (first value is background risk, then l1, l2, pi)
#' @param show_prob A value representing the probability of a woman showing up for the screening visit. Defaults to 0.9.
#' @param i A value representing the interval between screening rounds (in years). Defaults to 5.
#'
#' @return A data frame containing the left and right interval of CIN2+ detection, the indicator of prevalent disease, age (younger or older than 39, 1=younger),
#' HPV genotype (HPV16 or other, 1=HPV16), and cytology result (normal or abnormal, 1=abnormal).

#' @author Kelsi Kroon \email{k.kroon@amsterdamumc.nl}
#' @export
#'
simulator <- function(n, l1_x, l2_x, pi_x, params, show_prob = 0.9, i=5){
  require(Rlab)
  # Function to simulate data with baseline characteristics
  # Inputs:
  # n = sample size,
  # params = c(h, g0, g1, w0, p0, p1, p2, p3) i.e., model parameters,
  # show_prob = the probability the patient shows up for a screening round,
  # i = screening interval that observation times are scattered around

  # observation process - 5 screens
  # if the subject didn't show for screening then we make the corresponding screening time equal NA
  # if they do show up then their time is scattered around each screening round
  d <- 1
  screening_times <- data.frame(x1 = ifelse(rbern(n, show_prob), 0, NA),
                                x2 = ifelse(rbern(n, show_prob), rnorm(n, i*1, sqrt(d)), NA),
                                x3 = ifelse(rbern(n, show_prob), rnorm(n, i*2, sqrt(d)), NA),
                                x4 = ifelse(rbern(n, show_prob), rnorm(n, i*3, sqrt(d)), NA),
                                x5 = ifelse(rbern(n, show_prob), rnorm(n, i*4, sqrt(d)), NA))
  #x6 = ifelse(rbern(n, show_prob), rnorm(n, i*5, sqrt(d)), NA))
  # baseline characteristics
  # Age at baseline between 30 and 60
  age <- rbern(n, 1/3)
  # Cytology Results: this is an indicator variable so 1 means abnormal cytology and 0 means not abnormal (or unknown for z=NA)
  # if they did not show up for screening at time 0 then their cytology result is 0 because it is unknown
  cytology <- ifelse(is.na(screening_times$x1), 1, rbern(n, 0.5))

  # HPV genotype (HPV 16 or other) - this is an indicator variable so 1 means they have HPV16 and 0 means other HPV type
  hpv <- rbern(n, 1/3)

  covariate_data <- create.covariate.data(l1_x, l2_x, pi_x, data=data.frame(age=age, hpv=hpv, cyt=cytology))
  data1 <- covariate_data[[1]]
  data2 <- covariate_data[[2]]
  data3 <- covariate_data[[3]]

  n1 <- length(l1_x) + 1
  n2 <- length(l2_x) + 1
  n3 <- length(pi_x) + 1

  h <- params[1]
  params <- params[-1]
  l1 <- exp(data1 %*% params[1:(n1)])
  l2 <- exp(data2 %*% params[(n1+1):(n1+n2)])
  if (n3 == 1){
    # if pi has no covariates then the value of p is the same for all women
    p <- params[n1+n2+n3]
  }else if (n3 > 1){
    p <- exp(data3 %*% params[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% params[(n1+n2+1):(n1+n2+n3)]))
  }

  # disease process
  t1 <- ifelse(rbern(n, p)==1, 0, Inf) # prevalent CIN2/3
  t2 <- ifelse(rbern(n, l1/(l1+l2))==1, rexp(n, rate=(l1+l2)), Inf) # due to the HPV infection at baseline
  t3 <- rexp(n, rate=h) # due to the background risk (simulate value for everyone)
  t <- pmin(t1, t2, t3) # keep the minimum value as the actual event time

  screening_times$actual <- t
  # create a list of the screening times with NA values removed for each subject
  # the code makes the following change: 0   NA    t2    t3    NA  ---> 0    t2    t3
  # this allows us to find the interval where CIN2/3 was detected
  screens <- apply(screening_times, 1, function(x) c(na.omit(unlist(x, use.names=FALSE))))

  z <- rep(0, n) # create the indicator variable Z

  # create left intervals by finding the last value in the list of screens that is smaller than the actual event time ß
  left <- vapply(screens, function(x) x[Position(function(z) z <= x[length(x)], x[c(1:(length(x))-1)], right=TRUE)], 1)

  # if left interval is NA then disease was unknown at baseline because it was not checked
  z[is.na(left)] <- NA

  left[is.na(left)] <- 0 # set unknown left intervals to 0 because CIN2/3 could have been there at baseline

  # create a list of right intervals by finding the first value in the
  # list of screen times that is greater than the actual event time
  right <- vapply(screens, function(x) x[Position(function(z) z > x[length(x)], x[c(1:(length(x))-1)])], 1)

  # if the actual event time t=0 and left interval l=0 and the indicator is not unknown
  # (meaning disease was checked at baseline), then the right interval is also zero
  right[left==0 & t==0 & !is.na(z)] <-  0

  z[which(right==0)] <- 1 # right is only zero when disease is prevalent (defined above)

  # if the actual time of CIN2/3 development is after the last screening time, then the set the time to Inf
  last_screening_time <- vapply(screens, function(x) tail(x, 2)[1], 1)
  right[screening_times$actual > last_screening_time] <- Inf

  # if the right interval is NA then set it to infinity - this happens if all screening
  # rounds were NA (very rare, this is just to avoid errors in case it happens)
  right[is.na(right)] <- Inf

  return(data.frame(left, right, z = z, age = age, hpv = hpv, cytology))
}


#' Competing Cause Prevalence-Incidence Mixture Models
#'
#' This function fits Competing Cause Prevalence-Incidence mixture models to interval-censored cervical cancer screening data and obtains parameter estimates.
#' It is possible for the user to select the covariates that will be used for each parameter.
#' @import survival Rlab ggplot2
#' @param l1_x A vector containing the names of covariates used in the \ifelse{html}{\out{&lambda<sub>1</sub>}}{ \eqn{\lambda_1}} (progression rate) parameter (must match column name(s) in the input data)
#' @param l2_x A vector containing the names of covariates used in the \ifelse{html}{\out{&lambda<sub>2</sub>}}{ \eqn{\lambda_2}} (clearance rate) parameter (must match column name(s) in the input data)
#' @param pi_x A vector containing the names of covariates used in the \eqn{\pi} parameter (probability of prevalent disease) (must match column name(s) in the input data)
#' @param data Data used to fit the model containing columns for each term in l1_x, l2_x and pi_x. The first three columns must be (1)
#' left interval, (2) right interval and (3) z indicator for prevalent/incident disease.The data must contain observed prevalent and incident cases, however some cases may be unknown (`NA`).
#' @param num.runs Number of runs of the 'Short' EM algorithm used to determine initial values for the EM algorithm. Defaults to 30.
#' @param short.epsilon Convergence criteria used in the 'Short' EM algorithm to determine initial values. Defaults to 0.1.
#' @param epsilon Convergence criteria for the change in log-likelihood value used for stopping the EM algorithm. Defaults to 1e-08.

#' @return The output is a list containing the following elements:
#' \itemize{
#' \item initial.values - initial values determined by the Short EM algorithm process
#' \item theta.hat - optimum parameter values estimated by the EM algorithm
#' \item num.iterations - number of iterations until the EM algorithm converged
#' \item log.likelihood - value of the log.likelihood at the optimum parameter values
#' \item hess - hessian matrix
#' \item std.dev - standard deviation of parameter estimates
#' \item summary - data frame with estimate, std.dev, and 95\% CI for each parameter (useful for data set comparisons)
#' }
#'
#' @author Kelsi Kroon \email{k.kroon@amsterdamumc.nl}
#' @export
#'
#' @examples
#' sim_dat <- simulator(2000, c("hpv"), c(), c(),
#'                      c(exp(-5), -3, 2, -2, 0.25), show_prob = 0.9, i=5)
#' CCmixture.fit(c("hpv"), c(), c(), sim_dat, silent=F)

CCmixture.fit <- function(l1_x, l2_x, pi_x, data, num.runs=30, short.epsilon=1e-1, epsilon=1e-8, silent=T){
  short.inits <- list()
  for (i in 1:num.runs){
    short.inits[[i]] <- short.em.h(l1_x, l2_x, pi_x, data,  short.epsilon=1e-1, silent)[c("theta.hat", "log.likelihood")]
  }
  short.inits.mat <- matrix(unlist(short.inits), nrow=length(short.inits), byrow=T)
  # find the parameter values that results in the maximum likelihood
  init <- short.inits.mat[which.max(short.inits.mat[,dim(short.inits.mat)[2]]),1:(dim(short.inits.mat)[2]-1)]
  final.res <- em.function.h(init, l1_x, l2_x, pi_x, data, epsilon=epsilon, silent=silent)
  return(final.res)
}


#' Competing Cause Prevalence-Incidence Mixture Model Predictions
#'
#' @param l1_x A vector containing the names of covariates used in the \ifelse{html}{\out{\eqn{\lambda}<sub>1</sub>}}{ \eqn{\lambda_1}} (progression rate) parameter (must match column name(s) in the input data)
#' @param l2_x A vector containing the names of covariates used in the \eqn{\lambda_2} (clearance rate) parameter (must match column name(s) in the input data)
#' @param pi_x A vector containing the names of covariates used in the \eqn{\pi} parameter (probability of prevalent disease) (must match column name(s) in the input data)
#' @param data Data set of covariates from which to make the predictions.
#' @param time.points Numeric vector of time points used to make cumulative risk predictions
#' @param theta.hat Parameter estimates for the model to be used (output from CCmixture.fit)
#'
#' @export
CCmixture.predict <- function(l1_x, l2_x, pi_x, data, time.points, theta.hat){
  h <- exp(theta.hat[1])
  theta.hat <- theta.hat[-1]

  covariate_data <- create.covariate.data(l1_x, l2_x, pi_x, data=data)
  data1 <- covariate_data[[1]]
  data2 <- covariate_data[[2]]
  data3 <- covariate_data[[3]]

  n1 <- length(l1_x) + 1
  n2 <- length(l2_x) + 1
  n3 <- length(pi_x) + 1

  l1 <- as.numeric(exp(data1 %*% theta.hat[1:(n1)]))
  l2 <- as.numeric(exp(data2 %*% theta.hat[(n1+1):(n1+n2)]))
  if (n3 == 1){
    # if pi has no covariates then the value of p is the same for all women
    p <- theta.hat[n1+n2+n3]
  }else if (n3 > 1){
    p <- as.numeric(exp(data3 %*% theta.hat[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% theta.hat[(n1+n2+1):(n1+n2+n3)])))
  }

  prob <- p + (1-p)*(1 - exp(-h*time.points) + exp(-h*time.points)*g(l1, l2, time.points))
  return(prob)
}
