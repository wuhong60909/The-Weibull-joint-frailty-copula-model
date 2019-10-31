## Example 1 ##
time = c(0, 1, 3, 5)
p = c(0.3, 0.5, 0.7, 0.9)
jointCox.Weibull.prediction(time = time, p = p, Z1 = 1, Z2 = 1,
                            scale1 = 0.4, shape1 = 1.5, scale2 = 0.4, shape2 = 1.5,
                            eta = 0.5, theta = 2, beta1 = 1, beta2 = 1, alpha = 1)

jointCox.Weibull.prediction(time = time, p = p, Z1 = 1, Z2 = 1,
                            scale1 = 0.4, shape1 = 0.6, scale2 = 0.4, shape2 = 0.6,
                            eta = 0.5, theta = 2, beta1 = 1, beta2 = 1, alpha = 0.5)

## Example 2 ##
G = 30
N = 20
theta_true = 2
beta1_true = 1
beta2_true = 1
scale1_true = 5^(1 - 1.5)
scale2_true = 5^(1 - 2.5)
shape1_true = 1.5
shape2_true = 2.5
eta_true = 0.5
alpha_true = 1

t.event = t.death = event = death = Z1 = group = NULL
ij = 0
for(i in 1:G){
  u_i = rgamma(1, shape = 1/eta_true, scale = eta_true)
  for(j in 1:N){
    ij = ij + 1
    group[ij] = i
    Z1[ij] = runif(1)
    V1 = runif(1)
    V2 = runif(1)
    w = (1 + (V2^(-theta_true / (theta_true + 1)) - 1) * V1^(-theta_true))^(-1 / theta_true)
    X_ij = (-log(V1) / (scale1_true * u_i * exp(beta1_true * Z1[ij])))^(1 / shape1_true)
    D_ij = (-log(w) / (scale2_true * u_i^alpha_true * exp(beta2_true*Z1[ij])))^(1 / shape2_true)
    C_ij = runif(1, min = 0, max = 5)
    t.event[ij] = min(X_ij, D_ij, C_ij)
    t.death[ij] = min(D_ij, C_ij)
    event[ij] = as.numeric( t.event[ij] == X_ij )
    death[ij] = as.numeric( t.death[ij] == D_ij )
  }
}
obj = jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                     Z1 = Z1, Z2 = Z1, group = group, Adj = 500, alpha = alpha_true)

obj

time = c(0, 1, 3, 5)
p = c(0.3, 0.5, 0.7, 0.9)
jointCox.Weibull.prediction(object = obj, time = time, p = p, Z1 = 1, Z2 = 1)
#jointCox.Weibull.prediction(time = time, p = p, Z1 = 1, Z2 = 1,
#                            scale1 = obj$g[1, 1], shape1 = obj$g[2, 1], scale2 = obj$h[1, 1], shape2 = obj$h[2, 1],
#                            eta = obj$eta[1], theta = obj$theta[1], beta1 = obj$beta1[, 1], beta2 = obj$beta2[, 1],
#                            alpha = obj$alpha)

