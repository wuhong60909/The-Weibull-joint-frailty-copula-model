WeibullCox.reg = function(t.event, event, Z, Randomize_num = 10) {
  d = event
  N = length(t.event)
  T_mean = mean(t.event)
  Name = names(Z)
  Z = as.matrix(Z)
  p = ncol(Z)

  ## Likelihood function ##
  l.func = function(phi) {
    g = exp(pmin(phi[1:2], 500))
    b = phi[(2 + 1):(2 + p)]
    l = 0
    ## WEibull baseline hazard ##
    weibull.h = function(x, scale0, shape0) {
      scale0 * shape0 * x^(shape0 - 1)
    }
    weibull.H = function(x, scale0, shape0) {
      scale0 * x^(shape0)
    }
    r = as.vector(weibull.h(t.event, scale0 = g[1], shape0 = g[2]))
    R = as.vector(weibull.H(t.event, scale0 = g[1], shape0 = g[2]))
    bZ = as.numeric(Z %*% b)
    l = l + sum(d * (log(r) + bZ))
    l = l - sum(pmin(exp(bZ) * R, exp(500)))
    -l
  }
  p0 = rep(0, 2 + p)
  p0[1] = p0[1] - log(T_mean)
  res = nlm(l.func, p = p0, hessian = TRUE)
  MPL = -res$minimum
  R_num = 0
  repeat {
    if ((min(eigen(res$hessian)$values) > 0) & (res$code == 1)) {
      break
    }
    if (R_num >= Randomize_num) {
      break
    }
    R_num = R_num + 1
    p0_Rand = runif(2 + p, -1, 1)
    p0_Rand[1] = p0_Rand[1] - log(T_mean)
    res_Rand = nlm(l.func, p = p0_Rand, hessian = TRUE)
    MPL_Rand = -res_Rand$minimum
    if (MPL_Rand > MPL) {
      res = res_Rand
      MPL = -res$minimum
    }
  }
  H_PL = -res$hessian
  DF_upper = 18 + p
  temp = (det(H_PL) == 0) | is.na(det(H_PL))
  if (temp) {
    V = solve(-H_PL + diag(rep(1e-04, 2 + p)), tol = 10^(-50))
  }
  else {
    V = solve(-H_PL, tol = 10^(-50))
  }
  D_PL = diag(c(1/exp(res$estimate[1:2]), rep(1, p)))
  H_PL = D_PL %*% H_PL %*% D_PL
  H = H_PL
  if (is.na(det(H_PL)) | det(H_PL) == 0) {
    DF = DF_upper
  }
  else {
    DF = min(max(sum(diag(solve(H_PL, tol = 10^(-50)) %*% H)), p), DF_upper)
  }
  LCV = -l.func(res$estimate) - DF
  AIC = 2 * DF - 2 * (-l.func(res$estimate))
  BIC = log(N) * DF - 2 * (-l.func(res$estimate))
  convergence_res = c(MPL = MPL, DF = DF, LCV = LCV, AIC = AIC, BIC = BIC,
                      code = res$code, No.of.iterations = res$iterations,
                      No.of.randomizations = R_num)

  est = c(exp(res$est[1:2]), res$est[(2 + 1):(2 + p)])
  est_var = diag(c(est[1:2], rep(1, p))) %*% V %*% diag(c(est[1:2], rep(1, p)))

  beta_est = res$est[(2 + 1):(2 + p)]
  beta_se = sqrt(diag(V)[(2 + 1):(2 + p)])
  z_score = beta_est/beta_se
  p_value = 1 - pnorm(abs(z_score))
  g_est = exp(res$est[1:2])
  g_var = diag(g_est) %*% V[1:2, 1:2] %*% diag(g_est)
  g_se = sqrt(diag(g_var))

  beta_res = data.frame(Estimate = beta_est, SE = beta_se,
                        z_score = z_score, p_value = p_value,
                        row.names = Name)
  g_res = data.frame(Estimate = g_est, SE = g_se,
                     Lower = g_est * exp(-1.96 * sqrt(diag(V)[1: 2])),
                     Upper = g_est * exp(1.96 * sqrt(diag(V)[1: 2])),
                     row.names = c("scale", "shape"))

  ss = list(est = est, est_var = est_var,
            beta = beta_res, g = g_res, g_var = g_var,
            convergence = convergence_res)
  class(ss) = "WeibullCox.reg"
  ss
}
