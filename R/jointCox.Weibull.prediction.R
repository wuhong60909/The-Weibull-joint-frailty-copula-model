jointCox.Weibull.prediction = function(time, p, Z1, Z2,
                                       object = NULL,
                                       scale1 = NULL, shape1 = NULL,
                                       scale2 = NULL, shape2 = NULL,
                                       eta = NULL,
                                       theta = NULL,
                                       beta1 = NULL,
                                       beta2 = NULL,
                                       alpha = NULL,
                                       V = NULL,
                                       method = NULL,
                                       No.of.Monte.Carlo.simulations = 1000, ...) {
  if(class(object) == "jointCox.Weibull.reg") {
    scale1 = object$g[1, 1]
    shape1 = object$g[2, 1]
    scale2 = object$h[1, 1]
    shape2 = object$h[2, 1]
    eta = object$eta[1]
    theta = object$theta[1]
    beta1 = object$beta1[, 1]
    beta2 = object$beta2[, 1]
    alpha = object$alpha
    V = object$est_var
  }
  N.MC = No.of.Monte.Carlo.simulations
  meanX = meanD = mrlX = mrlD = varX = varD = NULL
  survX = survD = denX = denD = hazX = hazD = NULL
  quantX = quantD = NULL

  if(eta >= shape1/2) {
    if(eta >= shape1) {
      warning("eta >= shape1, the mean and MRL of X do not exist!")
    }
    warning("eta >= shape1/2, the variance of X does not exist!")
  }
  if(eta >= shape2/(2 * alpha)) {
    if(eta >= shape2/alpha) {
      warning("eta >= shape2/alpha, the mean and MRL of D do not exist!")
    }
    warning("eta >= shape2/(2 * alpha), the variance of D does not exist!")
  }
  if(alpha != 0 & alpha != 1){
    warning("When alpha != 0 and alpha != 1, the numerical solutions are required.
            In addition, MRL of D only provides the case when alpha is equal to 0 or 1.")
  }

  Estimate = MeanX(Z1 = Z1,
                   scale1 = scale1, shape1 = shape1,
                   scale2 = scale2, shape2 = shape2,
                   eta = eta, theta = theta,
                   beta1 = beta1, beta2 = beta2,
                   alpha = alpha,
                   V = V,
                   method = method,
                   No.of.Monte.Carlo.simulations = N.MC)
  meanX = data.frame(Z1, Estimate, row.names = NULL)

  Estimate = MeanD(Z2 = Z2,
                   scale1 = scale1, shape1 = shape1,
                   scale2 = scale2, shape2 = shape2,
                   eta = eta, theta = theta,
                   beta1 = beta1, beta2 = beta2,
                   alpha = alpha,
                   V = V,
                   method = method,
                   No.of.Monte.Carlo.simulations = N.MC)
  meanD = data.frame(Z2, Estimate, row.names = NULL)

  Estimate = VarX(Z1 = Z1,
                  scale1 = scale1, shape1 = shape1,
                  scale2 = scale2, shape2 = shape2,
                  eta = eta, theta = theta,
                  beta1 = beta1, beta2 = beta2,
                  alpha = alpha,
                  V = V,
                  method = method,
                  No.of.Monte.Carlo.simulations = N.MC)
  varX = data.frame(Z1, Estimate, row.names = NULL)

  Estimate = VarD(Z2 = Z2,
                  scale1 = scale1, shape1 = shape1,
                  scale2 = scale2, shape2 = shape2,
                  eta = eta, theta = theta,
                  beta1 = beta1, beta2 = beta2,
                  alpha = alpha,
                  V = V,
                  method = method,
                  No.of.Monte.Carlo.simulations = N.MC)
  varD = data.frame(Z2, Estimate, row.names = NULL)

  Estimate = MRLX(time = time, Z1 = Z1,
                  scale1 = scale1, shape1 = shape1,
                  scale2 = scale2, shape2 = shape2,
                  eta = eta, theta = theta,
                  beta1 = beta1, beta2 = beta2,
                  alpha = alpha,
                  V = V,
                  method = method,
                  No.of.Monte.Carlo.simulations = N.MC)
  mrlX = data.frame(Z1, t = time, Estimate, row.names = NULL)

  Estimate = MRLD(time = time, Z2 = Z2,
                  scale1 = scale1, shape1 = shape1,
                  scale2 = scale2, shape2 = shape2,
                  eta = eta, theta = theta,
                  beta1 = beta1, beta2 = beta2,
                  alpha = alpha,
                  V = V,
                  method = method,
                  No.of.Monte.Carlo.simulations = N.MC)
  mrlD = data.frame(Z2, t = time, Estimate, row.names = NULL)

  Estimate = SurvX(time = time, Z1 = Z1,
                   scale1 = scale1, shape1 = shape1,
                   scale2 = scale2, shape2 = shape2,
                   eta = eta, theta = theta,
                   beta1 = beta1, beta2 = beta2,
                   alpha = alpha,
                   V = V,
                   method = method,
                   No.of.Monte.Carlo.simulations = N.MC)
  survX = data.frame(Z1, t = time, Estimate, row.names = NULL)

  Estimate = SurvD(time = time, Z2 = Z2,
                   scale1 = scale1, shape1 = shape1,
                   scale2 = scale2, shape2 = shape2,
                   eta = eta, theta = theta,
                   beta1 = beta1, beta2 = beta2,
                   alpha = alpha,
                   V = V,
                   method = method,
                   No.of.Monte.Carlo.simulations = N.MC)
  survD = data.frame(Z2, t = time, Estimate, row.names = NULL)

  Estimate = DensityX(time = time, Z1 = Z1,
                      scale1 = scale1, shape1 = shape1,
                      scale2 = scale2, shape2 = shape2,
                      eta = eta, theta = theta,
                      beta1 = beta1, beta2 = beta2,
                      alpha = alpha,
                      V = V,
                      method = method,
                      No.of.Monte.Carlo.simulations = N.MC)
  denX = data.frame(Z1, t = time, Estimate, row.names = NULL)

  Estimate = DensityD(time = time, Z2 = Z2,
                      scale1 = scale1, shape1 = shape1,
                      scale2 = scale2, shape2 = shape2,
                      eta = eta, theta = theta,
                      beta1 = beta1, beta2 = beta2,
                      alpha = alpha,
                      V = V,
                      method = method,
                      No.of.Monte.Carlo.simulations = N.MC)
  denD = data.frame(Z2, t = time, Estimate, row.names = NULL)

  Estimate = HazardX(time = time, Z1 = Z1,
                     scale1 = scale1, shape1 = shape1,
                     scale2 = scale2, shape2 = shape2,
                     eta = eta, theta = theta,
                     beta1 = beta1, beta2 = beta2,
                     alpha = alpha,
                     V = V,
                     method = method,
                     No.of.Monte.Carlo.simulations = N.MC)
  hazX = data.frame(Z1, t = time, Estimate, row.names = NULL)

  Estimate = HazardD(time = time, Z2 = Z2,
                     scale1 = scale1, shape1 = shape1,
                     scale2 = scale2, shape2 = shape2,
                     eta = eta, theta = theta,
                     beta1 = beta1, beta2 = beta2,
                     alpha = alpha,
                     V = V,
                     method = method,
                     No.of.Monte.Carlo.simulations = N.MC)
  hazD = data.frame(Z2, t = time, Estimate, row.names = NULL)

  Estimate = QuantileX(p = p, Z1 = Z1,
                       scale1 = scale1, shape1 = shape1,
                       scale2 = scale2, shape2 = shape2,
                       eta = eta, theta = theta,
                       beta1 = beta1, beta2 = beta2,
                       alpha = alpha,
                       V = V,
                       method = method,
                       No.of.Monte.Carlo.simulations = N.MC)
  quantX = data.frame(Z1, p = p, Estimate, row.names = NULL)

  Estimate = QuantileD(p = p, Z2 = Z2,
                       scale1 = scale1, shape1 = shape1,
                       scale2 = scale2, shape2 = shape2,
                       eta = eta, theta = theta,
                       beta1 = beta1, beta2 = beta2,
                       alpha = alpha,
                       V = V,
                       method = method,
                       No.of.Monte.Carlo.simulations = N.MC)
  quantD = data.frame(Z2, p = p, Estimate, row.names = NULL)

  ss = list(Mean.of.X = meanX,
            Mean.of.D = meanD,
            Variance.of.X = varX,
            Variance.of.D = varD,
            Survival.of.X = survX,
            Survival.of.D = survD,
            Density.of.X = denX,
            Density.of.D = denD,
            Hazard.of.X = hazX,
            Hazard.of.D = hazD,
            MRL.of.X = mrlX,
            MRL.of.D = mrlD,
            Quantile.of.X = quantX,
            Quantile.of.D = quantD)
  ss
}

## 1. Mean of X ##
MeanX = function(Z1, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                 method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z1 = as.matrix(Z1)
  N = nrow(Z1)
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z1 %*% beta1)
    if(eta < shape1) {
      a = exp(lgamma(1/eta - 1/shape1) + (1/shape1) * log(1/eta) - lgamma(1/eta)) *
        (scale1 * exp(bz))^(-1/shape1) * gamma(1 + 1/shape1)
    }
    else {
      a = rep(NA, N)
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 2. Mean of D ##
MeanD = function(Z2, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                 method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z2 = as.matrix(Z2)
  N = nrow(Z2)
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z2 %*% beta2)
    if(eta < shape2/alpha) {
      a = exp(lgamma(1/eta - alpha/shape2) + (alpha/shape2) * log(1/eta) - lgamma(1/eta)) *
        (scale2 * exp(bz))^(-1/shape2) * gamma(1 + 1/shape2)
    }
    else {
      a = rep(NA, N)
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 3. Variance of X ##
VarX = function(Z1, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z1 = as.matrix(Z1)
  N = nrow(Z1)
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z1 %*% beta1)
    if(eta < shape1/2) {
      a1 = lgamma(1/eta - 2/shape1) + (2/shape1) * log(1/eta) - lgamma(1/eta)
      a2 = lgamma(1/eta - 1/shape1) + (1/shape1) * log(1/eta) - lgamma(1/eta)
      a = (exp(a1) * gamma(1 + 2/shape1) - exp(2 * a2) * (gamma(1 + 1/shape1))^2) *
        (scale1 * exp(bz))^(-2/shape1)
    }
    else {
      a = rep(NA, N)
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 4. Variance of D ##
VarD = function(Z2, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z2 = as.matrix(Z2)
  N = nrow(Z2)
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z2 %*% beta2)
    if(eta < shape2/(2 * alpha)) {
      a1 = lgamma(1/eta - 2 * alpha/shape2) + (2 * alpha/shape2) * log(1/eta) - lgamma(1/eta)
      a2 = lgamma(1/eta - alpha/shape2) + (alpha/shape2) * log(1/eta) - lgamma(1/eta)
      a = (exp(a1) * gamma(1 + 2/shape2) - exp(2 * a2) * (gamma(1 + 1/shape2))^2) *
        (scale2 * exp(bz))^(-2/shape2)
    }
    else {
      a = rep(NA, N)
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 5. Survival function of X ##
SurvX = function(time, Z1, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                 method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z1 = as.matrix(Z1)
  if(nrow(Z1) > 1 & length(time) > 1) {
    stop("Please fix one of the arguments time or Z1 first")
  }
  else if(nrow(Z1) == 1) {
    N = length(time)
  }
  else if(length(time) == 1) {
    N = nrow(Z1)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z1 %*% beta1)
    a1 = eta * scale1 * exp(bz) * time^shape1
    a = (1 + a1)^(-1/eta)
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 6. Survival function of D ##
SurvD = function(time, Z2, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                 method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z2 = as.matrix(Z2)
  if(nrow(Z2) > 1 & length(time) > 1) {
    stop("Please fix one of the arguments time or Z2 first")
  }
  else if(nrow(Z2) == 1) {
    N = length(time)
  }
  else if(length(time) == 1) {
    N = nrow(Z2)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z2 %*% beta2)
    ## weibull-gamma ##
    if(alpha == 1) {
      a1 = eta * scale2 * exp(bz) * time^shape2
      a = (1 + a1)^(-1/eta)
    }
    ## weibull ##
    else if(alpha == 0) {
      a = exp(-scale2 * exp(bz) * time^shape2)
    }
    ## numerical solution ##
    else {
      a1 = scale2 * exp(bz) * time^shape2
      Func = function(a) {
        InnerFunc = function(u) {
          exp(-u^alpha * a) * dgamma(u, shape = 1/eta, scale = eta)
        }
        Int = integrate(InnerFunc, lower = 0, upper = Inf, stop.on.error = FALSE)
        Int$value
      }
      s = sapply(a1, Func)
      a = s
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate, SE = SE, Lower = Lower, Upper = Upper, row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 7. Density function of X ##
DensityX = function(time, Z1, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                    method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z1 = as.matrix(Z1)
  if(nrow(Z1) > 1 & length(time) > 1) {
    stop("Please fix one of the arguments time or Z1 first")
  }
  else if(nrow(Z1) == 1) {
    N = length(time)
  }
  else if(length(time) == 1) {
    N = nrow(Z1)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z1 %*% beta1)
    a1 = scale1 * exp(bz) * shape1 * time^(shape1 - 1)
    a2 = (1 + eta * scale1 * exp(bz) * time^shape1)^(1/eta + 1)
    a = a1 / a2
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 8. Density function of D ##
DensityD = function(time, Z2, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                    method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z2 = as.matrix(Z2)
  if(nrow(Z2) > 1 & length(time) > 1) {
    stop("Please fix one of the arguments time or Z2 first")
  }
  else if(nrow(Z2) == 1) {
    N = length(time)
  }
  else if(length(time) == 1) {
    N = nrow(Z2)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z2 %*% beta2)
    ## weibull-gamma ##
    if(alpha == 1) {
      a1 = scale2 * exp(bz) * shape2 * time^(shape2 - 1)
      a2 = (1 + eta * scale2 * exp(bz) * time^shape2)^(1/eta + 1)
      a = a1 / a2
    }
    ## weibull ##
    else if(alpha == 0) {
      a1 = scale2 * exp(bz) * shape2 * time^(shape2 - 1)
      a2 = exp(-scale2 * exp(bz) * time^(shape2))
      a = a1 * a2
    }
    ## numerical solution ##
    else {
      a1 = scale2 * exp(bz) * shape2 * time^(shape2 - 1)
      a2 = scale2 * exp(bz) * time^shape2
      Func = function(a) {
        InnerFunc = function(u) {
          #u^alpha * exp(-u^alpha * a) * dgamma(u, shape = 1/eta, scale = eta)
          exp(alpha * log(u) - u^alpha * a + log(dgamma(u, shape = 1/eta, scale = eta)))
        }
        Int = integrate(InnerFunc, lower = 0, upper = Inf, stop.on.error = FALSE)
        Int$value
      }
      f = a1 * sapply(a2, Func)
      a = f
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 9. Hazard function of X ##
HazardX = function(time, Z1, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                   method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z1 = as.matrix(Z1)
  if(nrow(Z1) > 1 & length(time) > 1) {
    stop("Please fix one of the arguments time or Z1 first")
  }
  else if(nrow(Z1) == 1) {
    N = length(time)
  }
  else if(length(time) == 1) {
    N = nrow(Z1)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z1 %*% beta1)
    a1 = scale1 * exp(bz) * shape1 * time^(shape1 - 1)
    a2 = 1 + eta * scale1 * exp(bz) * time^shape1
    a = a1 / a2
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 10. Hazard function of D ##
HazardD = function(time, Z2, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                   method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z2 = as.matrix(Z2)
  if(nrow(Z2) > 1 & length(time) > 1) {
    stop("Please fix one of the arguments time or Z2 first")
  }
  else if(nrow(Z2) == 1) {
    N = length(time)
  }
  else if(length(time) == 1) {
    N = nrow(Z2)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z2 %*% beta2)
    ## weibull-gamma ##
    if(alpha == 1) {
      a1 = scale2 * exp(bz) * shape2 * time^(shape2 - 1)
      a2 = 1 + eta * scale2 * exp(bz) * time^shape2
      a = a1 / a2
    }
    ## weibull ##
    else if(alpha == 0) {
      a = scale2 * exp(bz) * shape2 * time^(shape2 - 1)
    }
    ## numerical solution ##
    else {
      a1 = scale2 * exp(bz) * shape2 * time^(shape2 - 1)
      a2 = scale2 * exp(bz) * time^shape2
      Func1 = function(a) {
        InnerFunc = function(u) {
          #exp(-u^alpha * a) * dgamma(u, shape = 1/eta, scale = eta)
          exp(-u^alpha * a + log(dgamma(u, shape = 1/eta, scale = eta)))
        }
        Int = integrate(InnerFunc, lower = 0, upper = Inf, stop.on.error = FALSE)
        Int$value
      }
      Func2 = function(a) {
        InnerFunc = function(u) {
          #u^alpha * exp(-u^alpha * a) * dgamma(u, shape = 1/eta, scale = eta)
          exp(alpha * log(u) - u^alpha * a + log(dgamma(u, shape = 1/eta, scale = eta)))
        }
        Int = integrate(InnerFunc, lower = 0, upper = Inf, stop.on.error = FALSE)
        Int$value
      }
      s = sapply(a2, Func1)
      f = a1 * sapply(a2, Func2)
      h = f / s
      a = h
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 11. MRL of X ##
MRLX = function(time = 0, Z1, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z1 = as.matrix(Z1)
  if(nrow(Z1) > 1 & length(time) > 1) {
    stop("Please fix one of the arguments time or Z1 first")
  }
  else if(nrow(Z1) == 1) {
    N = length(time)
  }
  else if(length(time) == 1) {
    N = nrow(Z1)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z1 %*% beta1)
    if(eta >= shape1) {
      a = rep(NA, N)
    }
    else {
      SurvX_func = function(x) {
        (1 + eta * scale1 * exp(bz) * x^(shape1))^(-1/eta)
      }
      IBeta = function(x, par.a, par.b) {
        l.func = function(w) {
          w^(par.a - 1) * (1 - w)^(par.b - 1)
        }
        Int = integrate(l.func, lower = 0, upper = x, stop.on.error = F)
        Int$value
      }
      a1 = 1/eta
      a2 = 1/shape1
      a3 = 1/(1 + eta * scale1 * exp(bz) * time^(shape1))
      l1 = a1^a2 * (scale1 * exp(bz))^(-a2) * a2 / SurvX_func(time)
      l2 = sapply(a3, FUN = IBeta, par.a = a1 - a2, par.b = a2)
      a = l1 * l2
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 12. MRL of D ##
MRLD = function(time = 0, Z2, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z2 = as.matrix(Z2)
  if(nrow(Z2) > 1 & length(time) > 1) {
    stop("Please fix one of the arguments time or Z2 first")
  }
  else if(nrow(Z2) == 1) {
    N = length(time)
  }
  else if(length(time) == 1) {
    N = nrow(Z2)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z2 %*% beta2)
    if(eta >= shape2/alpha) {
      a = rep(NA, N)
    }
    else {
      ## weibull-gamma ##
      if(alpha == 1) {
        SurvD_func = function(y) {
          (1 + eta * scale2 * exp(bz) * y^(shape2))^(-1/eta)
        }
        IBeta = function(x, par.a, par.b) {
          l.func = function(w) {
            w^(par.a - 1) * (1 - w)^(par.b - 1)
          }
          Int = integrate(l.func, lower = 0, upper = x, stop.on.error = F)
          Int$value
        }
        a1 = 1/eta
        a2 = 1/shape2
        a3 = 1/(1 + eta * scale2 * exp(bz) * time^(shape2))
        l1 = a1^a2 * (scale2 * exp(bz))^(-a2) * a2 / SurvD_func(time)
        l2 = sapply(a3, FUN = IBeta, par.a = a1 - a2, par.b = a2)
        a = l1 * l2
      }
      ## weibull ##
      else if(alpha == 0) {
        SurvD_func = function(y) {
          exp(-scale2 * exp(bz) * y^shape2)
        }
        IGamma = function(x, par.a) {
          l.func = function(w) {
            w^(par.a - 1) * exp(-w)
          }
          Int = integrate(l.func, lower = x, upper = Inf, stop.on.error = F)
          Int$value
        }
        a1 = 1/shape2
        a2 = scale2 * exp(bz) * time^shape2
        l1 = (scale2 * exp(bz))^(-a1) * a1 / SurvD_func(time)
        l2 = sapply(a2, FUN = IGamma, par.a = a1)
        a = l1 * l2
      }
      ## numerical solution ..........##
      else {
        a = rep(NA, N)
      }
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 13. Quantile of X ##
QuantileX = function(p, Z1, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                     method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z1 = as.matrix(Z1)
  if(nrow(Z1) > 1 & length(p) > 1) {
    stop("Please fix one of the arguments p or Z1 first")
  }
  else if(nrow(Z1) == 1) {
    N = length(p)
  }
  else if(length(p) == 1) {
    N = nrow(Z1)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z1 %*% beta1)
    a1 = (1 - (1 - p)^eta) / (eta * (1 - p)^eta)
    a2 = scale1 * exp(bz)
    a = a1^(1/shape1) * a2^(-1/shape1)
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}

## 14. Quantile of D ##
QuantileD = function(p, Z2, scale1, shape1, scale2, shape2, eta, theta, beta1, beta2, alpha, V = NULL,
                     method = NULL, No.of.Monte.Carlo.simulations = 1000, ...) {
  p1 = length(beta1)
  p2 = length(beta2)
  est = c(scale1, shape1, scale2, shape2, eta, theta, beta1, beta2)
  Z2 = as.matrix(Z2)
  if(nrow(Z2) > 1 & length(p) > 1) {
    stop("Please fix one of the arguments p or Z2 first")
  }
  else if(nrow(Z2) == 1) {
    N = length(p)
  }
  else if(length(p) == 1) {
    N = nrow(Z2)
  }
  est_func = function(phi) {
    scale1 = phi[1]
    shape1 = phi[2]
    scale2 = phi[3]
    shape2 = phi[4]
    eta = phi[5]
    theta = phi[6]
    beta1 = phi[(6 + 1):(6 + p1)]
    beta2 = phi[(6 + p1 + 1):(6 + p1 + p2)]
    bz = as.vector(Z2 %*% beta2)
    ## weibull-gamma ##
    if(alpha == 1) {
      a1 = (1 - (1 - p)^eta) / (eta * (1 - p)^eta)
      a2 = scale2 * exp(bz)
      a = a1^(1/shape2) * a2^(-1/shape2)
    }
    ## weibull ##
    else if(alpha == 0) {
      a1 = -log(1 - p)
      a2 = scale2 * exp(bz)
      a = a1^(1/shape2) * a2^(-1/shape2)
    }
    ## numerical solution ##
    else {
      ## inverse function ##
      inverse = function(f, interval = c(-100, 100)) {
        function(y) {
          uniroot((function(x) f(x) - y), interval = interval)$root
        }
      }
      f1 = function(p, bz) {
        ## Survival function of D ##
        SurvD_func = function(t) {
          InnerFunc2 = function(y) {
            InnerFunc1 = function(u) {
              exp(-u^alpha * scale2 * exp(bz) * y^shape2) * dgamma(u, shape = 1/eta, scale = eta)
            }
            Int = integrate(InnerFunc1, lower = 0, upper = Inf, stop.on.error = FALSE)
            Int$value
          }
          sapply(t, InnerFunc2)
        }
        inv.SurvD_func = inverse(SurvD_func, interval = c(0, exp(500)))
        inv.SurvD_func(1 - p)
      }
      if(length(bz) == 1) {
        a = sapply(p, f1, bz = bz)
      }
      else if(length(p) == 1) {
        a = sapply(bz, f1, p = p)
      }
    }
    a
  }
  if(is.null(method) | is.null(V)) {
    result = est_func(est)
  }
  else if(method == "Delta Method") {
    if(min(eigen(V)$values) > 0) {
      J = jacobian(est_func, est)
      est_var = diag(J %*% V %*% t(J))
      Estimate = est_func(est)
      SE = sqrt(est_var)
      Lower = est_func(est) - 1.96 * sqrt(est_var)
      Upper = est_func(est) + 1.96 * sqrt(est_var)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else if(method == "Monte Carlo") {
    if(min(eigen(V)$values) > 0) {
      est_tilde = est
      est_tilde[1:6] = log(est_tilde[1:6])
      V_tilde = diag(c(1/est[1:6], rep(1, p1 + p2))) %*% V %*% diag(c(1/est[1:6], rep(1, p1 + p2)))
      est_tilde_samples = mvrnorm(n = No.of.Monte.Carlo.simulations, mu = est_tilde, Sigma = V_tilde)
      est_tilde_samples[, 1:6] = exp(est_tilde_samples[, 1:6])
      est_func_hat = apply(est_tilde_samples, MARGIN = 1, FUN = est_func)
      est_func_hat = matrix(est_func_hat, nrow = N)
      CI = apply(est_func_hat, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      Lower = CI[1, ]
      Upper = CI[2, ]
      SE = apply(est_func_hat, MARGIN = 1, FUN = sd, na.rm = TRUE)
      Estimate = est_func(est)
      result = data.frame(Estimate = Estimate,
                          SE = SE,
                          Lower = Lower,
                          Upper = Upper,
                          row.names = NULL)
    }
    else {
      result = data.frame(Estimate = est_func(est),
                          SE = NA,
                          Lower = NA,
                          Upper = NA,
                          row.names = NULL)
    }
  }
  else {
    result = est_func(est)
  }
  result
}
