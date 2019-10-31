library(numDeriv)
library(MASS)
getwd()
weibull.h = function(x, scale0, shape0)
{
  scale0 * shape0 * x^(shape0 - 1)
}

weibull.H = function(x, scale0, shape0)
{
  scale0 * x^(shape0)
}

########## scenario I, G = 5, Ni = 100 or 200 ##########
G = 5; N = 100; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 700;
G = 5; N = 100; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 700;
G = 5; N = 100; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 700;
G = 5; N = 100; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 700;

G = 5; N = 100; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 700;
G = 5; N = 100; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 700;
G = 5; N = 100; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 700;
G = 5; N = 100; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 700;

G = 5; N = 100; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 700;
G = 5; N = 100; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 700;
G = 5; N = 100; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 700;
G = 5; N = 100; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 600;
G = 5; N = 200; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 700;


########## scenario II, G = 30, Ni = 10 or 20 ##########
G = 30; N = 10; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 10; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;
G = 30; N = 10; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 10; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 0.8); shape1_true = 0.8; scale2_true = 5^(1 - 0.8); shape2_true = 0.8; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;

G = 30; N = 10; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 10; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;
G = 30; N = 10; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 10; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 1.5); shape1_true = 1.5; scale2_true = 5^(1 - 1.5); shape2_true = 1.5; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;

G = 30; N = 10; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 2; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 10; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 2; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;
G = 30; N = 10; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 6; beta1_true = 1; beta2_true = 1; alpha_true = 1; Adj = 500;
G = 30; N = 10; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;
G = 30; N = 20; scale1_true = 5^(1 - 2.5); shape1_true = 2.5; scale2_true = 5^(1 - 2.5); shape2_true = 2.5; eta_true = 0.5; theta_true = 6; beta1_true = -1; beta2_true = -1; alpha_true = 1; Adj = 500;


## Weibull prediction ##
z = 0.5
p = 0.5
time = 3

A0 = jointCox.Weibull.prediction(time = time, p = p, Z1 = z, Z2 = z,
                                 scale1 = scale1_true, shape1 = shape1_true,
                                 scale2 = scale2_true, shape2 = shape2_true,
                                 eta = eta_true, theta = theta_true,
                                 beta1 = beta1_true, beta2 = beta2_true,
                                 alpha = alpha_true)

## true median ##
MedianX_true = A0$Quantile.of.X$Estimate
MedianD_true = A0$Quantile.of.D$Estimate

## true mean ##
MeanX_true = A0$Mean.of.X$Estimate
MeanD_true = A0$Mean.of.D$Estimate

## true MRL at time = 3 ##
MRLX_true = A0$MRL.of.X$Estimate
MRLD_true = A0$MRL.of.D$Estimate


# Weibull simulation =================================================
R = 10
## Estimate ## =======================================================
beta1_hat = c()
beta2_hat = c()
eta_hat = c()
theta_hat = c()
scale1_hat = c()
shape1_hat = c()
scale2_hat = c()
shape2_hat = c()
## Delta Method ##
MedianX_hat.Delta.Method = c()
MedianD_hat.Delta.Method = c()
MeanX_hat.Delta.Method = c()
MeanD_hat.Delta.Method = c()
MRLX_hat.Delta.Method = c()
MRLD_hat.Delta.Method = c()
## Monte Carlo ##
MedianX_hat.Monte.Carlo = c()
MedianD_hat.Monte.Carlo = c()
MeanX_hat.Monte.Carlo = c()
MeanD_hat.Monte.Carlo = c()
MRLX_hat.Monte.Carlo = c()
MRLD_hat.Monte.Carlo = c()
## Estimate ## =======================================================

## SE ## =============================================================
beta1_SE = c()
beta2_SE = c()
eta_SE = c()
theta_SE = c()
scale1_SE = c()
shape1_SE = c()
scale2_SE = c()
shape2_SE = c()
## Delta Method ##
MedianX_SE.Delta.Method = c()
MedianD_SE.Delta.Method = c()
MeanX_SE.Delta.Method = c()
MeanD_SE.Delta.Method = c()
MRLX_SE.Delta.Method = c()
MRLD_SE.Delta.Method = c()
## Monte Carlo ##
MedianX_SE.Monte.Carlo = c()
MedianD_SE.Monte.Carlo = c()
MeanX_SE.Monte.Carlo = c()
MeanD_SE.Monte.Carlo = c()
MRLX_SE.Monte.Carlo = c()
MRLD_SE.Monte.Carlo = c()
## SE ## =============================================================

## Lower ## =============================================================
beta1_low = c()
beta2_low = c()
eta_low = c()
theta_low = c()
scale1_low = c()
shape1_low = c()
scale2_low = c()
shape2_low = c()
## Delta Method ##
MedianX_low.Delta.Method = c()
MedianD_low.Delta.Method = c()
MeanX_low.Delta.Method = c()
MeanD_low.Delta.Method = c()
MRLX_low.Delta.Method = c()
MRLD_low.Delta.Method = c()
## Monte Carlo ##
MedianX_low.Monte.Carlo = c()
MedianD_low.Monte.Carlo = c()
MeanX_low.Monte.Carlo = c()
MeanD_low.Monte.Carlo = c()
MRLX_low.Monte.Carlo = c()
MRLD_low.Monte.Carlo = c()
## Lower ## =============================================================

## Upper ## =============================================================
beta1_up = c()
beta2_up = c()
eta_up = c()
theta_up = c()
scale1_up = c()
shape1_up = c()
scale2_up = c()
shape2_up = c()
## Delta Method ##
MedianX_up.Delta.Method = c()
MedianD_up.Delta.Method = c()
MeanX_up.Delta.Method = c()
MeanD_up.Delta.Method = c()
MRLX_up.Delta.Method = c()
MRLD_up.Delta.Method = c()
## Monte Carlo ##
MedianX_up.Monte.Carlo = c()
MedianD_up.Monte.Carlo = c()
MeanX_up.Monte.Carlo = c()
MeanD_up.Monte.Carlo = c()
MRLX_up.Monte.Carlo = c()
MRLD_up.Monte.Carlo = c()
## Upper ## =============================================================

CEN = c()
MPL = c()
DF = c()
LCV = c()
AIC = c()
BIC = c()
code = c()
No.of.iterations = c()
No.of.randomizations = c()
ADJ = c()

err = 1:R
ptm = proc.time()
for(ii in err)
{
  set.seed(ii)
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
  CEN[ii] = mean((event == 0) & (death == 0))

  repeat{
    check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                     Z1 = Z1, Z2 = Z1, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
    if(class(check) != "try-error" | Adj < 0){break}
    Adj = Adj - 50
  }
  joint = check
  ADJ[ii] = Adj

  MPL[ii] = joint$convergence["MPL"]
  DF[ii] = joint$convergence["DF"]
  LCV[ii] = joint$convergence["LCV"]
  AIC[ii] = joint$convergence["AIC"]
  BIC[ii] = joint$convergence["BIC"]
  code[ii] = joint$convergence["code"]
  No.of.iterations[ii] = joint$convergence["No.of.iterations"]
  No.of.randomizations[ii] = joint$convergence["No.of.randomizations"]

  beta1_hat[ii] = joint$beta1[1, "Estimate"]
  beta2_hat[ii] = joint$beta2[1, "Estimate"]
  eta_hat[ii] = joint$eta["Estimate"]
  theta_hat[ii] = joint$theta["Estimate"]
  scale1_hat[ii] = joint$g[1, "Estimate"]
  shape1_hat[ii] = joint$g[2, "Estimate"]
  scale2_hat[ii] = joint$h[1, "Estimate"]
  shape2_hat[ii] = joint$h[2, "Estimate"]

  beta1_SE[ii] = joint$beta1[1, "SE"]
  beta2_SE[ii] = joint$beta2[1, "SE"]
  eta_SE[ii] = joint$eta["SE"]
  theta_SE[ii] = joint$theta["SE"]
  scale1_SE[ii] = joint$g[1, "SE"]
  shape1_SE[ii] = joint$g[2, "SE"]
  scale2_SE[ii] = joint$h[1, "SE"]
  shape2_SE[ii] = joint$h[2, "SE"]

  beta1_low[ii] = joint$beta1[1, "Lower"]
  beta2_low[ii] = joint$beta2[1, "Lower"]
  eta_low[ii] = joint$eta["Lower"]
  theta_low[ii] = joint$theta["Lower"]
  scale1_low[ii] = joint$g[1, "Lower"]
  shape1_low[ii] = joint$g[2, "Lower"]
  scale2_low[ii] = joint$h[1, "Lower"]
  shape2_low[ii] = joint$h[2, "Lower"]

  beta1_up[ii] = joint$beta1[1, "Upper"]
  beta2_up[ii] = joint$beta2[1, "Upper"]
  eta_up[ii] = joint$eta["Upper"]
  theta_up[ii] = joint$theta["Upper"]
  scale1_up[ii] = joint$g[1, "Upper"]
  shape1_up[ii] = joint$g[2, "Upper"]
  scale2_up[ii] = joint$h[1, "Upper"]
  shape2_up[ii] = joint$h[2, "Upper"]

  ## Delta Method ## =========================================================
  A1 = jointCox.Weibull.prediction(object = joint,
                                   time = time, p = p, Z1 = z, Z2 = z,
                                   method = "Delta Method")

  QX_D = A1$Quantile.of.X
  QD_D = A1$Quantile.of.D
  MeanX_D = A1$Mean.of.X
  MeanD_D = A1$Mean.of.D
  MRLX_D = A1$MRL.of.X
  MRLD_D = A1$MRL.of.D

  ## Estimate ##
  MedianX_hat.Delta.Method[ii] = QX_D$Estimate
  MedianD_hat.Delta.Method[ii] = QD_D$Estimate
  MeanX_hat.Delta.Method[ii] = MeanX_D$Estimate
  MeanD_hat.Delta.Method[ii] = MeanD_D$Estimate
  MRLX_hat.Delta.Method[ii] = MRLX_D$Estimate
  MRLD_hat.Delta.Method[ii] = MRLD_D$Estimate

  ## SE ##
  MedianX_SE.Delta.Method[ii] = QX_D$SE
  MedianD_SE.Delta.Method[ii] = QD_D$SE
  MeanX_SE.Delta.Method[ii] = MeanX_D$SE
  MeanD_SE.Delta.Method[ii] = MeanD_D$SE
  MRLX_SE.Delta.Method[ii] = MRLX_D$SE
  MRLD_SE.Delta.Method[ii] = MRLD_D$SE

  ## Lower ##
  MedianX_low.Delta.Method[ii] = QX_D$Lower
  MedianD_low.Delta.Method[ii] = QD_D$Lower
  MeanX_low.Delta.Method[ii] = MeanX_D$Lower
  MeanD_low.Delta.Method[ii] = MeanD_D$Lower
  MRLX_low.Delta.Method[ii] = MRLX_D$Lower
  MRLD_low.Delta.Method[ii] = MRLD_D$Lower

  ## Upper ##
  MedianX_up.Delta.Method[ii] = QX_D$Upper
  MedianD_up.Delta.Method[ii] = QD_D$Upper
  MeanX_up.Delta.Method[ii] = MeanX_D$Upper
  MeanD_up.Delta.Method[ii] = MeanD_D$Upper
  MRLX_up.Delta.Method[ii] = MRLX_D$Upper
  MRLD_up.Delta.Method[ii] = MRLD_D$Upper
  ## Delta Method ## =========================================================


  ## Monte Carlo ## ==========================================================
  A2 = jointCox.Weibull.prediction(object = joint,
                                   time = time, p = p, Z1 = z, Z2 = z,
                                   method = "Monte Carlo")

  QX_M = A2$Quantile.of.X
  QD_M = A2$Quantile.of.D
  MeanX_M = A2$Mean.of.X
  MeanD_M = A2$Mean.of.D
  MRLX_M = A2$MRL.of.X
  MRLD_M = A2$MRL.of.D

  ## Estimate ##
  MedianX_hat.Monte.Carlo[ii] = QX_M$Estimate
  MedianD_hat.Monte.Carlo[ii] = QD_M$Estimate
  MeanX_hat.Monte.Carlo[ii] = MeanX_M$Estimate
  MeanD_hat.Monte.Carlo[ii] = MeanD_M$Estimate
  MRLX_hat.Monte.Carlo[ii] = MRLX_M$Estimate
  MRLD_hat.Monte.Carlo[ii] = MRLD_M$Estimate

  ## SE ##
  MedianX_SE.Monte.Carlo[ii] = QX_M$SE
  MedianD_SE.Monte.Carlo[ii] = QD_M$SE
  MeanX_SE.Monte.Carlo[ii] = MeanX_M$SE
  MeanD_SE.Monte.Carlo[ii] = MeanD_M$SE
  MRLX_SE.Monte.Carlo[ii] = MRLX_M$SE
  MRLD_SE.Monte.Carlo[ii] = MRLD_M$SE

  ## Lower ##
  MedianX_low.Monte.Carlo[ii] = QX_M$Lower
  MedianD_low.Monte.Carlo[ii] = QD_M$Lower
  MeanX_low.Monte.Carlo[ii] = MeanX_M$Lower
  MeanD_low.Monte.Carlo[ii] = MeanD_M$Lower
  MRLX_low.Monte.Carlo[ii] = MRLX_M$Lower
  MRLD_low.Monte.Carlo[ii] = MRLD_M$Lower

  ## Upper ##
  MedianX_up.Monte.Carlo[ii] = QX_M$Upper
  MedianD_up.Monte.Carlo[ii] = QD_M$Upper
  MeanX_up.Monte.Carlo[ii] = MeanX_M$Upper
  MeanD_up.Monte.Carlo[ii] = MeanD_M$Upper
  MRLX_up.Monte.Carlo[ii] = MRLX_M$Upper
  MRLD_up.Monte.Carlo[ii] = MRLD_M$Upper
  ## Monte Carlo ## ==========================================================

  ## plot baseline hazard for TTP ##
  if(ii == 1)
  {
    curve(weibull.h(x, scale0 = scale1_true, shape0 = shape1_true), from = 0, to = 5, ylim = c(0, 5),
          xlab = "time", ylab = "Baseline hazard for TTP", type = "l", lwd = 2, col = 1)
    curve(weibull.h(x, scale0 = joint$g[1, "Estimate"], shape0 = joint$g[2, "Estimate"]),
          from = 0, to = 5, col = 2, add = T)
  }
  if(ii > 1 & ii <= R)
  {
    curve(weibull.h(x, scale0 = joint$g[1, "Estimate"], shape0 = joint$g[2, "Estimate"]),
          from = 0, to = 5, col = 2, add = T)
  }
  legend("topleft", paste(ii))
}
take_time = proc.time() - ptm
curve(weibull.h(x, scale0 = scale1_true, shape0 = shape1_true), from = 0, to = 5, ylim = c(0, 5),
      xlab = "time", ylab = "Baseline hazard for TTP", type = "l", lwd = 2, col = 1, add = T)
legend("topright", c("True", "Estimated"), lwd = c(2, 1), col = c(1, 2))

result = cbind(ADJ, CEN, MPL, DF, LCV, AIC, BIC, code, No.of.iterations, No.of.randomizations,
               beta1_hat, beta2_hat, eta_hat, theta_hat, scale1_hat, shape1_hat, scale2_hat, shape2_hat,
               beta1_SE, beta2_SE, eta_SE, theta_SE, scale1_SE, shape1_SE, scale2_SE, shape2_SE,
               beta1_low, beta2_low, eta_low, theta_low, scale1_low, shape1_low, scale2_low, shape2_low,
               beta1_up, beta2_up, eta_up, theta_up, scale1_up, shape1_up, scale2_up, shape2_up,
               MedianX_hat.Delta.Method, MedianX_SE.Delta.Method, MedianX_low.Delta.Method, MedianX_up.Delta.Method,
               MedianD_hat.Delta.Method, MedianD_SE.Delta.Method, MedianD_low.Delta.Method, MedianD_up.Delta.Method,
               MeanX_hat.Delta.Method, MeanX_SE.Delta.Method, MeanX_low.Delta.Method, MeanX_up.Delta.Method,
               MeanD_hat.Delta.Method, MeanD_SE.Delta.Method, MeanD_low.Delta.Method, MeanD_up.Delta.Method,
               MedianX_hat.Monte.Carlo, MedianX_SE.Monte.Carlo, MedianX_low.Monte.Carlo, MedianX_up.Monte.Carlo,
               MedianD_hat.Monte.Carlo, MedianD_SE.Monte.Carlo, MedianD_low.Monte.Carlo, MedianD_up.Monte.Carlo,
               MeanX_hat.Monte.Carlo, MeanX_SE.Monte.Carlo, MeanX_low.Monte.Carlo, MeanX_up.Monte.Carlo,
               MeanD_hat.Monte.Carlo, MeanD_SE.Monte.Carlo, MeanD_low.Monte.Carlo, MeanD_up.Monte.Carlo)

True = c(beta1 = beta1_true, beta2 = beta2_true,
         eta = eta_true, theta = theta_true,
         scale1 = scale1_true, shape1 = shape1_true,
         scale2 = scale2_true, shape2 = shape2_true,
         MedianX.D = MedianX_true,
         MedianD.D = MedianD_true,
         MeanX.D = MeanX_true,
         MeanD.D = MeanD_true,
         MRLX.D = MRLX_true,
         MRLD.D = MRLD_true,
         MedianX.M = MedianX_true,
         MedianD.M = MedianD_true,
         MeanX.M = MeanX_true,
         MeanD.M = MeanD_true,
         MRLX.M = MRLX_true,
         MRLD.M = MRLD_true
)

Mean = c(beta1 = mean(beta1_hat), beta2 = mean(beta2_hat),
         eta = mean(eta_hat), theta = mean(theta_hat),
         scale1 = mean(scale1_hat), shape1 = mean(shape1_hat),
         scale2 = mean(scale2_hat), shape2 = mean(shape2_hat),
         MedianX.D = mean(MedianX_hat.Delta.Method, na.rm = TRUE),
         MedianD.D = mean(MedianD_hat.Delta.Method, na.rm = TRUE),
         MeanX.D = mean(MeanX_hat.Delta.Method, na.rm = TRUE),
         MeanD.D = mean(MeanD_hat.Delta.Method, na.rm = TRUE),
         MRLX.D = mean(MRLX_hat.Delta.Method, na.rm = TRUE),
         MRLD.D = mean(MRLD_hat.Delta.Method, na.rm = TRUE),
         MedianX.M = mean(MedianX_hat.Monte.Carlo, na.rm = TRUE),
         MedianD.M = mean(MedianD_hat.Monte.Carlo, na.rm = TRUE),
         MeanX.M = mean(MeanX_hat.Monte.Carlo, na.rm = TRUE),
         MeanD.M = mean(MeanD_hat.Monte.Carlo, na.rm = TRUE),
         MRLX.M = mean(MRLX_hat.Monte.Carlo, na.rm = TRUE),
         MRLD.M = mean(MRLD_hat.Monte.Carlo, na.rm = TRUE)
)

SD = c(beta1 = sd(beta1_hat), beta2 = sd(beta2_hat),
       eta = sd(eta_hat), theta = sd(theta_hat),
       scale1 = sd(scale1_hat), shape1 = sd(shape1_hat),
       scale2 = sd(scale2_hat), shape2 = sd(shape2_hat),
       MedianX.D = sd(MedianX_hat.Delta.Method, na.rm = TRUE),
       MedianD.D = sd(MedianD_hat.Delta.Method, na.rm = TRUE),
       MeanX.D = sd(MeanX_hat.Delta.Method, na.rm = TRUE),
       MeanD.D = sd(MeanD_hat.Delta.Method, na.rm = TRUE),
       MRLX.D = sd(MRLX_hat.Delta.Method, na.rm = TRUE),
       MRLD.D = sd(MRLD_hat.Delta.Method, na.rm = TRUE),
       MedianX.M = sd(MedianX_hat.Monte.Carlo, na.rm = TRUE),
       MedianD.M = sd(MedianD_hat.Monte.Carlo, na.rm = TRUE),
       MeanX.M = sd(MeanX_hat.Monte.Carlo, na.rm = TRUE),
       MeanD.M = sd(MeanD_hat.Monte.Carlo, na.rm = TRUE),
       MRLX.M = sd(MRLX_hat.Monte.Carlo, na.rm = TRUE),
       MRLD.M = sd(MRLD_hat.Monte.Carlo, na.rm = TRUE)
)

SE = c(beta1 = mean(beta1_SE, na.rm = TRUE), beta2 = mean(beta2_SE, na.rm = TRUE),
       eta = mean(eta_SE, na.rm = TRUE), theta = mean(theta_SE, na.rm = TRUE),
       scale1 = mean(scale1_SE, na.rm = TRUE), shape1 = mean(shape1_SE, na.rm = TRUE),
       scale2 = mean(scale2_SE, na.rm = TRUE), shape2 = mean(shape2_SE, na.rm = TRUE),
       MedianX.D = mean(MedianX_SE.Delta.Method, na.rm = TRUE),
       MedianD.D = mean(MedianD_SE.Delta.Method, na.rm = TRUE),
       MeanX.D = mean(MeanX_SE.Delta.Method, na.rm = TRUE),
       MeanD.D = mean(MeanD_SE.Delta.Method, na.rm = TRUE),
       MRLX.D = mean(MRLX_SE.Delta.Method, na.rm = TRUE),
       MRLD.D = mean(MRLD_SE.Delta.Method, na.rm = TRUE),
       MedianX.M = mean(MedianX_SE.Monte.Carlo, na.rm = TRUE),
       MedianD.M = mean(MedianD_SE.Monte.Carlo, na.rm = TRUE),
       MeanX.M = mean(MeanX_SE.Monte.Carlo, na.rm = TRUE),
       MeanD.M = mean(MeanD_SE.Monte.Carlo, na.rm = TRUE),
       MRLX.M = mean(MRLX_SE.Monte.Carlo, na.rm = TRUE),
       MRLD.M = mean(MRLD_SE.Monte.Carlo, na.rm = TRUE)
)

CP = c(beta1 = mean((beta1_low < beta1_true) & (beta1_true < beta1_up), na.rm = TRUE),
       beta2 = mean((beta2_low < beta2_true) & (beta2_true < beta2_up), na.rm = TRUE),
       eta = mean((eta_low < eta_true) & (eta_true < eta_up), na.rm = TRUE),
       theta = mean((theta_low < theta_true) & (theta_low < theta_up), na.rm = TRUE),
       scale1 = mean((scale1_low < scale1_true) & (scale1_true < scale1_up), na.rm = TRUE),
       shape1 = mean((shape1_low < shape1_true) & (shape1_true < shape1_up), na.rm = TRUE),
       scale2 = mean((scale2_low < scale2_true) & (scale2_true < scale2_up), na.rm = TRUE),
       shape2 = mean((shape2_low < shape2_true) & (shape2_true < shape2_up), na.rm = TRUE),
       MedianX.D = mean((MedianX_low.Delta.Method < MedianX_true) & (MedianX_true < MedianX_up.Delta.Method), na.rm = TRUE),
       MedianD.D = mean((MedianD_low.Delta.Method < MedianD_true) & (MedianD_true < MedianD_up.Delta.Method), na.rm = TRUE),
       MeanX.D = mean((MeanX_low.Delta.Method < MeanX_true) & (MeanX_true < MeanX_up.Delta.Method), na.rm = TRUE),
       MeanD.D = mean((MeanD_low.Delta.Method < MeanD_true) & (MeanD_true < MeanD_up.Delta.Method), na.rm = TRUE),
       MRLX.D = mean((MRLX_low.Delta.Method < MRLX_true) & (MRLX_true < MRLX_up.Delta.Method), na.rm = TRUE),
       MRLD.D = mean((MRLD_low.Delta.Method < MRLD_true) & (MRLD_true < MRLD_up.Delta.Method), na.rm = TRUE),
       MedianX.M = mean((MedianX_low.Monte.Carlo < MedianX_true) & (MedianX_true < MedianX_up.Monte.Carlo), na.rm = TRUE),
       MedianD.M = mean((MedianD_low.Monte.Carlo < MedianD_true) & (MedianD_true < MedianD_up.Monte.Carlo), na.rm = TRUE),
       MeanX.M = mean((MeanX_low.Monte.Carlo < MeanX_true) & (MeanX_true < MeanX_up.Monte.Carlo), na.rm = TRUE),
       MeanD.M = mean((MeanD_low.Monte.Carlo < MeanD_true) & (MeanD_true < MeanD_up.Monte.Carlo), na.rm = TRUE),
       MRLX.M = mean((MRLX_low.Monte.Carlo < MRLX_true) & (MRLX_true < MRLX_up.Monte.Carlo), na.rm = TRUE),
       MRLD.M = mean((MRLD_low.Monte.Carlo < MRLD_true) & (MRLD_true < MRLD_up.Monte.Carlo), na.rm = TRUE)
)

ERR = c(beta1 = mean(is.na(beta1_SE)), beta2 = mean(is.na(beta2_SE)),
        eta = mean(is.na(eta_SE)), theta = mean(is.na(theta_SE)),
        scale1 = mean(is.na(scale1_SE)), shape1 = mean(is.na(shape1_SE)),
        scale2 = mean(is.na(scale2_SE)), shape2 = mean(is.na(shape2_SE)),
        MedianX.D = mean(is.na(MedianX_SE.Delta.Method)),
        MedianD.D = mean(is.na(MedianD_SE.Delta.Method)),
        MeanX.D = mean(is.na(MeanX_SE.Delta.Method)),
        MeanD.D = mean(is.na(MeanD_SE.Delta.Method)),
        MRLX.D = mean(is.na(MRLX_SE.Delta.Method)),
        MRLD.D = mean(is.na(MRLD_SE.Delta.Method)),
        MedianX.M = mean(is.na(MedianX_SE.Monte.Carlo)),
        MedianD.M = mean(is.na(MedianD_SE.Monte.Carlo)),
        MeanX.M = mean(is.na(MeanX_SE.Monte.Carlo)),
        MeanD.M = mean(is.na(MeanD_SE.Monte.Carlo)),
        MRLX.M = mean(is.na(MRLX_SE.Monte.Carlo)),
        MRLD.M = mean(is.na(MRLD_SE.Monte.Carlo))
)
Parameter = c("beta1", "beta2", "eta", "theta", "scale1", "shape1", "scale2", "shape2",
              "MedianX.D", "MedianD.D", "MeanX.D", "MeanD.D", "MRLX.D", "MRLD.D",
              "MedianX.M", "MedianD.M", "MeanX.M", "MeanD.M", "MRLX.M", "MRLD.M")
voilate1 = sum(eta_hat > shape1_hat)
voilate2 = sum(eta_hat > shape2_hat/alpha_true)
data_structure = data.frame(G = G, N = N, CEN = mean(CEN), code = mean(code == 1),
                            voilate1 = voilate1, voilate2 = voilate2)
joint_data = data.frame(Parameter = Parameter, True = round(True, 3), Mean = round(Mean, 3), SD = round(SD, 3), SE = round(SE, 3), CP = round(CP, 2), ERR = ERR)
data_structure
joint_data
mean(code == 1)
take_time

## output ##
write.table(joint_data, file = "joint.weibull.csv", sep = ",", row.names = F)
write.table(data_structure, file = "data_structure.csv", sep = ",", row.names = F)
write.table(result, file = "result.csv", sep = ",", row.names = F)

# plot1 ==============================================================
png("Baseline_hazard_TTP.png", width = 600, height = 450)
m = seq(1, R)
for(ii in m)
{
  if(ii == m[1])
  {
    curve(weibull.h(x, scale0 = scale1_hat[ii], shape0 = shape1_hat[ii]), from = 0, to = 5, ylim = c(0, 5),
          xlab = "time", ylab = "Baseline hazard for TTP", type = "l", col = 2)
  }
  if(ii > m[1])
  {
    curve(weibull.h(x, scale0 = scale1_hat[ii], shape0 = shape1_hat[ii]),
          from = 0, to = 5, col = 2, add = T)
  }
}
curve(weibull.h(x, scale0 = scale1_true, shape0 = shape1_true),
      from = 0, to = 5, lwd = 2, col = 1, add = T)
legend("topright", c("True", "Estimated"), lwd = c(2, 1), col = c(1, 2))
dev.off()
# plot1 ==============================================================


# plot2 ==============================================================
png("Baseline_hazard_OS.png", width = 600, height = 450)
m = seq(1, R)
for(ii in m)
{
  if(ii == m[1])
  {
    curve(weibull.h(x, scale0 = scale2_hat[ii], shape0 = shape2_hat[ii]), from = 0, to = 5, ylim = c(0, 5),
          xlab = "time", ylab = "Baseline hazard for OS", type = "l", col = 2)
  }
  if(ii > m[1])
  {
    curve(weibull.h(x, scale0 = scale2_hat[ii], shape0 = shape2_hat[ii]),
          from = 0, to = 5, col = 2, add = T)
  }
}
curve(weibull.h(x, scale0 = scale2_true, shape0 = shape2_true),
      from = 0, to = 5, lwd = 2, col = 1, add = T)
legend("topright", c("True", "Estimated"), lwd = c(2, 1), col = c(1, 2))
dev.off()
# plot2 ==============================================================


## Find out the abnormal results manually, replace different Adj and do it again.
err = which(scale1_hat > 10)










