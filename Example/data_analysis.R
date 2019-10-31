library(survival)

## data analysis ## ===================================
ovarian_cancer = read.table(file = "ovarian_cancer.csv", header = T, sep = ",")
index = which(is.na(ovarian_cancer$tumorstage))
ovarian_cancer = ovarian_cancer[-index, ]

t.event = ovarian_cancer$days_to_tumor_recurrence
event = ovarian_cancer$recurrence_status
t.death = ovarian_cancer$days_to_death
death = ovarian_cancer$vital_status
group = ovarian_cancer$group

t.event[t.event == 0] = 1 ## data can not be zero ##
t.death[t.death == 0] = 1 ## data can not be zero ##

CEN = mean(event == 0 & death == 0)
CEN

## median follow-up ##
death2 = death
death2[death == 1] = 0
death2[death == 0] = 1
s = Surv(time = t.death, event = death2)
ss = survfit(s ~ group)
ss
length(death)

## descriptive statistics ##
Z = ovarian_cancer[, c("CXCL12", "tumorstage", "summarystage")]
Z$tumorstage = factor(Z$tumorstage)
Z$summarystage = factor(Z$summarystage, levels = c("early", "late"))
summary(Z)


## tumorstage ## ===========================================================================
## TTP ##
km.for.X = survfit(Surv(time = t.event, event = event) ~ Z$tumorstage, conf.type = "log-log")
survdiff.for.X = survdiff(Surv(time = t.event, event = event) ~ Z$tumorstage)
p_value.for.X = 1 - pchisq(survdiff.for.X$chisq, length(survdiff.for.X$n) - 1)
p_value.for.X

plot(km.for.X, mark.time = F, xlim = c(0, 6420), col = c(4, 3, 2, 1), lwd = 2,
     xlab = "time (days)", ylab = "Survival probability for TTP",
     main = NULL, conf.int = F)
text(4000, 0.4, 'log-rank P = ')
text(5100, 0.4, round(p_value.for.X, digits = 9))
legend("topright", c("stage I (n = 39)", "stage II (n = 44)", "stage III (n = 788)", "stage IV(n = 128)"),
       col = c(4, 3, 2, 1), lty = c(1, 1, 1, 1), lwd = c(2, 2, 2, 2))

## OS ##
km.for.D = survfit(Surv(time = t.death, event = death) ~ Z$tumorstage, conf.type = "log-log")
survdiff.for.D = survdiff(Surv(time = t.death, event = death) ~ Z$tumorstage)
p_value.for.D = 1 - pchisq(survdiff.for.D$chisq, length(survdiff.for.D$n) - 1)
p_value.for.D

plot(km.for.D, mark.time = F, xlim = c(0, 6420), col = c(4, 3, 2, 1), lwd = 2,
     xlab = "time (days)", ylab = "Survival probability for OS",
     main = NULL, conf.int = F)
text(4000, 0.4, 'log-rank P = ')
text(5200, 0.4, round(p_value.for.D, digits = 9))
legend("topright", c("stage I (n = 39)", "stage II (n = 44)", "stage III (n = 788)", "stage IV(n = 128)"),
       col = c(4, 3, 2, 1), lty = c(1, 1, 1, 1), lwd = c(2, 2, 2, 2))
## tumorstage ## ===========================================================================

## summarystage ## =========================================================================
## TTP ##
km.for.X = survfit(Surv(time = t.event, event = event) ~ Z$summarystage, conf.type = "log-log")
survdiff.for.X = survdiff(Surv(time = t.event, event = event) ~ Z$summarystage)
p_value.for.X = 1 - pchisq(survdiff.for.X$chisq, length(survdiff.for.D$n) - 1)
p_value.for.X

plot(km.for.X, mark.time = F, xlim = c(0, 6420), col = c(4, 2), lwd = 2,
     xlab = "time (days)", ylab = "Survival probability for TTP",
     main = NULL, conf.int = F)

text(4000, 0.5, 'log-rank P = ')
text(5200, 0.5, round(p_value.for.X, digits = 9))
legend("topright", c("stage : early (n = 83)", "stage : late (n = 916)"),
       col = c(4, 2), lty = c(1, 1), lwd = c(2, 2))

## OS ##
km.for.D = survfit(Surv(time = t.death, event = death) ~ Z$summarystage, conf.type = "log-log")
survdiff.for.D = survdiff(Surv(time = t.death, event = death) ~ Z$summarystage)
p_value.for.D = 1 - pchisq(survdiff.for.D$chisq, length(survdiff.for.D$n) - 1)
p_value.for.D

plot(km.for.D, mark.time = F, xlim = c(0, 6420), col = c(4, 2), lwd = 2,
     xlab = "time (days)", ylab = "Survival probability for OS",
     main = NULL, conf.int = F)
text(3500, 0.5, 'log-rank P = ', round(p_value.for.D, digits = 3))
text(5500, 0.5, round(p_value.for.D, digits = 9))
legend("topright", c("stage : early (n = 83)", "stage : late (n = 916)"),
       col = c(4, 2), lty = c(1, 1), lwd = c(2, 2))
## summarystage ## =========================================================================



## data transformation ##
trst = data.frame(CXCL12 = Z$CXCL12, tumorstage.2 = 0, tumorstage.3 = 0, tumorstage.4 = 0, summarystage.late = 0)
trst$tumorstage.2[Z$tumorstage == 2] = 1
trst$tumorstage.3[Z$tumorstage == 3] = 1
trst$tumorstage.4[Z$tumorstage == 4] = 1
trst$summarystage.late[Z$summarystage == "late"] = 1
Z = trst

set.seed(1)
## case 1_1, Z1 = CXCL12, Z2 = CXCL12 ## =============================================
Z1 = Z$CXCL12
Z2 = Z$CXCL12
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint1_1 = check
Adj
take.time
## case 1_1, Z1 = CXCL12, Z2 = CXCL12 ## =============================================

set.seed(1)
## case 1_2, Z1 = CXCL12, Z2 = tumorstage ## =============================================
Z1 = Z$CXCL12
Z2 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]

alpha_true = 0
Adj = 600
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint1_2 = check
Adj
take.time
## case 1_2, Z1 = CXCL12, Z2 = tumorstage ## =============================================

set.seed(1)
## case 1_3, Z1 = CXCL12, Z2 = summarystage ## =============================================
Z1 = Z$CXCL12
Z2 = Z$summarystage.late
alpha_true = 0
Adj = 600
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint1_3 = check
Adj
take.time
## case 1_3, Z1 = CXCL12, Z2 = summarystage ## =============================================

set.seed(1)
## case 2_1, Z1 = tumorstage, Z2 = CXCL12 ## =============================================
Z1 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z$CXCL12
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint2_1 = check
Adj
take.time
## case 2_1, Z1 = tumorstage, Z2 = CXCL12 ## =============================================

set.seed(1)
## case 2_2, Z1 = tumorstage, Z2 = tumorstage ## =============================================
Z1 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint2_2 = check
Adj
take.time
## case 2_2, Z1 = tumorstage, Z2 = tumorstage ## =============================================

set.seed(1)
## case 2_3, Z1 = tumorstage, Z2 = summarystage ## =============================================
Z1 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z$summarystage.late
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint2_3 = check
Adj
take.time
## case 2_3, Z1 = tumorstage, Z2 = summarystage ## =============================================

set.seed(1)
## case 3_1, Z1 = summarystage, Z2 = CXCL12 ## =============================================
Z1 = Z$summarystage.late
Z2 = Z$CXCL12
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint3_1 = check
Adj
take.time
## case 3_1, Z1 = summarystage, Z2 = CXCL12 ## =============================================

set.seed(1)
## case 3_2, Z1 = summarystage, Z2 = tumorstage ## =============================================
Z1 = Z$summarystage.late
Z2 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]
alpha_true = 0
Adj = 600
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint3_2 = check
Adj
take.time
## case 3_2, Z1 = summarystage, Z2 = tumorstage ## =============================================

set.seed(1)
## case 3_3, Z1 = summarystage, Z2 = summarystage ## =============================================
Z1 = Z$summarystage.late
Z2 = Z$summarystage.late
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint3_3 = check
Adj
take.time
## case 3_3, Z1 = summarystage, Z2 = summarystage ## =============================================

set.seed(1)
## case 1_12, Z1 = CXCL12, Z2 = CXCL12 + tumorstage ## =============================================
Z1 = Z$CXCL12
Z2 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint1_12 = check
Adj
take.time
## case 1_12, Z1 = CXCL12, Z2 = CXCL12 + tumorstage ## =============================================

set.seed(1)
## case 1_13, Z1 = CXCL12, Z2 = CXCL12 + summarystage ## =============================================
Z1 = Z$CXCL12
Z2 = Z[, c("CXCL12", "summarystage.late")]
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint1_13 = check
Adj
take.time
## case 1_13, Z1 = CXCL12, Z2 = CXCL12 + summarystage ## =============================================

set.seed(1)
## case 2_12, Z1 = tumorstage, Z2 = CXCL12 + tumorstage ## =============================================
Z1 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
alpha_true = 0
Adj = 600
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint2_12 = check
Adj
take.time
## case 2_12, Z1 = tumorstage, Z2 = CXCL12 + tumorstage ## =============================================

set.seed(1)
## case 2_13, Z1 = tumorstage, Z2 = CXCL12 + summarystage ## =============================================
Z1 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z[, c("CXCL12", "summarystage.late")]
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint2_13 = check
Adj
take.time
## case 2_13, Z1 = tumorstage, Z2 = CXCL12 + summarystage ## =============================================

set.seed(1)
## case 3_12, Z1 = summarystage, Z2 = CXCL12 + tumorstage ## =============================================
Z1 = Z$summarystage.late
Z2 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint3_12 = check
Adj
take.time
## case 3_12, Z1 = summarystage, Z2 = CXCL12 + tumorstage ## =============================================

set.seed(1)
## case 3_13, Z1 = summarystage, Z2 = CXCL12 + summarystage ## =============================================
Z1 = Z$summarystage.late
Z2 = Z[, c("CXCL12", "summarystage.late")]
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint3_13 = check
Adj
take.time
## case 3_13, Z1 = summarystage, Z2 = CXCL12 + summarystage ## =============================================

set.seed(1)
## case 12_1, Z1 = CXCL12 + tumorstage, Z2 = CXCL12 ## =============================================
Z1 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z$CXCL12
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint12_1 = check
Adj
take.time
## case 12_1, Z1 = CXCL12 + tumorstage, Z2 = CXCL12 ## =============================================

set.seed(1)
## case 13_1, Z1 = CXCL12 + summarystage, Z2 = CXCL12 ## =============================================
Z1 = Z[, c("CXCL12", "summarystage.late")]
Z2 = Z$CXCL12
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint13_1 = check
Adj
take.time
## case 13_1, Z1 = CXCL12 + summarystage, Z2 = CXCL12 ## =============================================

set.seed(1)
## case 12_2, Z1 = CXCL12 + tumorstage, Z2 = tumorstage ## =============================================
Z1 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]
alpha_true = 0
Adj = 700
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint12_2 = check
Adj
take.time
## case 12_2, Z1 = CXCL12 + tumorstage, Z2 = tumorstage ## =============================================

set.seed(1)
## case 13_2, Z1 = CXCL12 + summarystage, Z2 = tumorstage ## =============================================
Z1 = Z[, c("CXCL12", "summarystage.late")]
Z2 = Z[, c("tumorstage.2", "tumorstage.3", "tumorstage.4")]
alpha_true = 0
Adj = 650
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint13_2 = check
Adj
take.time
## case 13_2, Z1 = CXCL12 + summarystage, Z2 = tumorstage ## =============================================

set.seed(1)
## case 12_3, Z1 = CXCL12 + tumorstage, Z2 = summarystage ## =============================================
Z1 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z$summarystage.late
alpha_true = 0
Adj = 400
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint12_3 = check
Adj
take.time
## case 12_3, Z1 = CXCL12 + tumorstage, Z2 = summarystage ## =============================================

set.seed(1)
## case 13_3, Z1 = CXCL12 + summarystage, Z2 = summarystage ## =============================================
Z1 = Z[, c("CXCL12", "summarystage.late")]
Z2 = Z$summarystage.late
alpha_true = 0
Adj = 650
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint13_3 = check
Adj
take.time
## case 13_3, Z1 = CXCL12 + summarystage, Z2 = summarystage ## =============================================

set.seed(1)
## case 12_12, Z1 = CXCL12 + tumorstage, Z2 = CXCL12 + tumorstage ## =============================================
Z1 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
alpha_true = 0
Adj = 600
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint12_12 = check
Adj
take.time
## case 12_12, Z1 = CXCL12 + tumorstage, Z2 = CXCL12 + tumorstage ## =============================================

set.seed(1)
## case 12_13, Z1 = CXCL12 + tumorstage, Z2 = CXCL12 + summarystage ## =============================================
Z1 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
Z2 = Z[, c("CXCL12", "summarystage.late")]
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint12_13 = check
Adj
take.time
## case 12_13, Z1 = CXCL12 + tumorstage, Z2 = CXCL12 + summarystage ## =============================================

set.seed(1)
## case 13_12, Z1 = CXCL12 + summarystage, Z2 = CXCL12 + tumorstage ## =============================================
Z1 = Z[, c("CXCL12", "summarystage.late")]
Z2 = Z[, c("CXCL12", "tumorstage.2", "tumorstage.3", "tumorstage.4")]
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint13_12 = check
Adj
take.time
## case 13_12, Z1 = CXCL12 + summarystage, Z2 = CXCL12 + tumorstage ## =============================================

set.seed(1)
## case 13_13, Z1 = CXCL12 + summarystage, Z2 = CXCL12 + summarystage ## =============================================
Z1 = Z[, c("CXCL12", "summarystage.late")]
Z2 = Z[, c("CXCL12", "summarystage.late")]
alpha_true = 0
Adj = 550
ptm = proc.time()
repeat{
  check = try(jointCox.Weibull.reg(t.event = t.event, event = event, t.death = t.death, death = death,
                                   Z1 = Z1, Z2 = Z2, group = group, Randomize_num = 10, Adj = Adj, alpha = alpha_true))
  if(class(check) != "try-error" | Adj <= 0){break}
  Adj = Adj - 50
}
take.time = proc.time() - ptm
joint13_13 = check
Adj
take.time
## case 13_13, Z1 = CXCL12 + summarystage, Z2 = CXCL12 + summarystage ## =============================================


## eta ##
joint1_1$eta
joint1_2$eta
joint1_3$eta
joint2_1$eta
joint2_2$eta
joint2_3$eta
joint3_1$eta
joint3_2$eta
joint3_3$eta

joint1_12$eta
joint1_13$eta
joint2_12$eta
joint2_13$eta
joint3_12$eta
joint3_13$eta

joint12_1$eta
joint13_1$eta
joint12_2$eta
joint13_2$eta
joint12_3$eta
joint13_3$eta

joint12_12$eta
joint12_13$eta
joint13_12$eta
joint13_13$eta


## AIC and BIC ##
AIC = c(
  joint1_1$convergence["AIC"],
  joint1_2$convergence["AIC"],
  joint1_3$convergence["AIC"],
  joint2_1$convergence["AIC"],
  joint2_2$convergence["AIC"],
  joint2_3$convergence["AIC"],
  joint3_1$convergence["AIC"],
  joint3_2$convergence["AIC"],
  joint3_3$convergence["AIC"],

  joint1_12$convergence["AIC"],
  joint1_13$convergence["AIC"],
  joint2_12$convergence["AIC"],
  joint2_13$convergence["AIC"],
  joint3_12$convergence["AIC"],
  joint3_13$convergence["AIC"],

  joint12_1$convergence["AIC"],
  joint13_1$convergence["AIC"],
  joint12_2$convergence["AIC"],
  joint13_2$convergence["AIC"],
  joint12_3$convergence["AIC"],
  joint13_3$convergence["AIC"],

  joint12_12$convergence["AIC"],
  joint12_13$convergence["AIC"],
  joint13_12$convergence["AIC"],
  joint13_13$convergence["AIC"]
)

names(AIC) = c(
  "joint1_1",
  "joint1_2",
  "joint1_3",
  "joint2_1",
  "joint2_2",
  "joint2_3",
  "joint3_1",
  "joint3_2",
  "joint3_3",

  "joint1_12",
  "joint1_13",
  "joint2_12",
  "joint2_13",
  "joint3_12",
  "joint3_13",

  "joint12_1",
  "joint13_1",
  "joint12_2",
  "joint13_2",
  "joint12_3",
  "joint13_3",

  "joint12_12",
  "joint12_13",
  "joint13_12",
  "joint13_13"
)

BIC = c(
  joint1_1$convergence["BIC"],
  joint1_2$convergence["BIC"],
  joint1_3$convergence["BIC"],
  joint2_1$convergence["BIC"],
  joint2_2$convergence["BIC"],
  joint2_3$convergence["BIC"],
  joint3_1$convergence["BIC"],
  joint3_2$convergence["BIC"],
  joint3_3$convergence["BIC"],

  joint1_12$convergence["BIC"],
  joint1_13$convergence["BIC"],
  joint2_12$convergence["BIC"],
  joint2_13$convergence["BIC"],
  joint3_12$convergence["BIC"],
  joint3_13$convergence["BIC"],

  joint12_1$convergence["BIC"],
  joint13_1$convergence["BIC"],
  joint12_2$convergence["BIC"],
  joint13_2$convergence["BIC"],
  joint12_3$convergence["BIC"],
  joint13_3$convergence["BIC"],

  joint12_12$convergence["BIC"],
  joint12_13$convergence["BIC"],
  joint13_12$convergence["BIC"],
  joint13_13$convergence["BIC"]
)

names(BIC) = c(
  "joint1_1",
  "joint1_2",
  "joint1_3",
  "joint2_1",
  "joint2_2",
  "joint2_3",
  "joint3_1",
  "joint3_2",
  "joint3_3",

  "joint1_12",
  "joint1_13",
  "joint2_12",
  "joint2_13",
  "joint3_12",
  "joint3_13",

  "joint12_1",
  "joint13_1",
  "joint12_2",
  "joint13_2",
  "joint12_3",
  "joint13_3",

  "joint12_12",
  "joint12_13",
  "joint13_12",
  "joint13_13"
)

DF = c(
  joint1_1$convergence["DF"],
  joint1_2$convergence["DF"],
  joint1_3$convergence["DF"],
  joint2_1$convergence["DF"],
  joint2_2$convergence["DF"],
  joint2_3$convergence["DF"],
  joint3_1$convergence["DF"],
  joint3_2$convergence["DF"],
  joint3_3$convergence["DF"],

  joint1_12$convergence["DF"],
  joint1_13$convergence["DF"],
  joint2_12$convergence["DF"],
  joint2_13$convergence["DF"],
  joint3_12$convergence["DF"],
  joint3_13$convergence["DF"],

  joint12_1$convergence["DF"],
  joint13_1$convergence["DF"],
  joint12_2$convergence["DF"],
  joint13_2$convergence["DF"],
  joint12_3$convergence["DF"],
  joint13_3$convergence["DF"],

  joint12_12$convergence["DF"],
  joint12_13$convergence["DF"],
  joint13_12$convergence["DF"],
  joint13_13$convergence["DF"]
)

names(DF) = c(
  "joint1_1",
  "joint1_2",
  "joint1_3",
  "joint2_1",
  "joint2_2",
  "joint2_3",
  "joint3_1",
  "joint3_2",
  "joint3_3",

  "joint1_12",
  "joint1_13",
  "joint2_12",
  "joint2_13",
  "joint3_12",
  "joint3_13",

  "joint12_1",
  "joint13_1",
  "joint12_2",
  "joint13_2",
  "joint12_3",
  "joint13_3",

  "joint12_12",
  "joint12_13",
  "joint13_12",
  "joint13_13"
)

res1 = data.frame(AIC = names(AIC)[order(rank(AIC))], DF = DF[order(rank(AIC))], row.names = NULL)
res1

res2 = data.frame(BIC = names(BIC)[order(rank(BIC))], DF = DF[order(rank(BIC))], row.names = NULL)
res2

joint12_12$eta
joint12_13$eta
joint12_3$eta



## joint12_12 ## ===========================================
## beta1 ##
round(exp(joint12_12$beta1[, c(1, 3, 4)]), digits = 2)
## beta2 ##
round(exp(joint12_12$beta2[, c(1, 3, 4)]), digits = 2)
## shape1 ##
round(joint12_12$g[2, c(1, 3, 4)], digits = 3)
## shape2 ##
round(joint12_12$h[2, c(1, 3, 4)], digits = 3)
## eta ##
round(joint12_12$eta[c(1, 3, 4)], digits = 3)
## theta ##
round(joint12_12$theta[c(1, 3, 4)], digits = 3)
## tau ##
round(joint12_12$tau[c(1, 3, 4)], digits = 3)
## alpha ##
joint12_12$alpha

## AIC ##
joint12_12$convergence["AIC"]
## BIC ##
joint12_12$convergence["BIC"]
## DF $$
joint12_12$convergence["DF"]
## joint12_12 ## ===========================================


## joint12_13 ## ===========================================
## beta1 ##
round(exp(joint12_13$beta1[, c(1, 3, 4)]), digits = 2)
## beta2 ##
round(exp(joint12_13$beta2[, c(1, 3, 4)]), digits = 2)
## shape1 ##
round(joint12_13$g[2, c(1, 3, 4)], digits = 3)
## shape2 ##
round(joint12_13$h[2, c(1, 3, 4)], digits = 3)
## eta ##
round(joint12_13$eta[c(1, 3, 4)], digits = 3)
## theta ##
round(joint12_13$theta[c(1, 3, 4)], digits = 3)
## tau ##
round(joint12_13$tau[c(1, 3, 4)], digits = 3)
## alpha ##
joint12_13$alpha

## AIC ##
joint12_13$convergence["AIC"]
## BIC ##
joint12_13$convergence["BIC"]
## DF ##
joint12_12$convergence["DF"]
## joint12_13 ## ===========================================




## Table 13 ~ 16 ##
joint = joint12_13
time = 1825
p = 0.5

## CXCL12 = -2 -1 0 1 2
ZZ = data.frame(CXCL12 = c(rep(-2, 4), rep(-1, 4), rep(0, 4), rep(1, 4), rep(2, 4)),
                stage = rep(c(1, 2, 3, 4), 5), summarystage = 0)

ZZ$summarystage[ZZ$stage <= 2] = "early"
ZZ$summarystage[ZZ$stage > 2] = "late"

## Z1 ##
Z1 = data.frame(CXCL12 = ZZ$CXCL12, tumorstage.2 = 0, tumorstage.3 = 0, tumorstage.4 = 0)
Z1$tumorstage.2[ZZ$stage == 2] = 1
Z1$tumorstage.3[ZZ$stage == 3] = 1
Z1$tumorstage.4[ZZ$stage == 4] = 1


## Z2 ##
Z2 = data.frame(CXCL12 = ZZ$CXCL12, summarystage.late = 0)
Z2$summarystage.late[ZZ$summarystage == "late"] = 1



A = jointCox.Weibull.prediction(object = joint, time = time, p = p,
                                Z1 = Z1, Z2 = Z2,
                                method = "Monte Carlo")

## TTP ##
median.X = data.frame(CXCL12 = ZZ$CXCL12, stage = ZZ$stage, round(A$Quantile.of.X[, c("Estimate", "Lower", "Upper")]))
mean.X = data.frame(CXCL12 = ZZ$CXCL12, stage = ZZ$stage, round(A$Mean.of.X[, c("Estimate", "Lower", "Upper")]))
mrl.X = data.frame(CXCL12 = ZZ$CXCL12, stage = ZZ$stage, round(A$MRL.of.X[, c("Estimate", "Lower", "Upper")]))
surv.X = data.frame(CXCL12 = ZZ$CXCL12, stage = ZZ$stage, round(A$Survival.of.X[, c("Estimate", "Lower", "Upper")], digits = 3))
median.X
mean.X
mrl.X
surv.X

## OS ##
median.D = data.frame(CXCL12 = ZZ$CXCL12, summarystage = ZZ$summarystage, round(A$Quantile.of.D[, c("Estimate", "Lower", "Upper")]))
mean.D = data.frame(CXCL12 = ZZ$CXCL12, summarystage = ZZ$summarystage, round(A$Mean.of.D[, c("Estimate", "Lower", "Upper")]))
mrl.D = data.frame(CXCL12 = ZZ$CXCL12, summarystage = ZZ$summarystage, round(A$MRL.of.D[, c("Estimate", "Lower", "Upper")]))
surv.D = data.frame(CXCL12 = ZZ$CXCL12, summarystage = ZZ$summarystage, round(A$Survival.of.D[, c("Estimate", "Lower", "Upper")], digits = 3))
unique(median.D)
unique(mean.D)
unique(mrl.D)
unique(surv.D)












## Figure 11, 12 ## =====================================================================
joint = joint12_13
time = seq(0, 6420, length.out = 100)
p = 0.5

## CXCL12 = -2, 2
CXCL12 = -2
col0 = c(4, 3, 2, 1)
lty0 = c(1, 2, 4, 5)

## TTP ## =============================================================================
## MRLX ##
plot(0, 0, type = "l", xlim = c(0, 6420), ylim = c(0, 10000),
     xlab = "time (days)", ylab = "MRL for TTP", main = paste("CXCL12 = ", CXCL12))
for (ii in 1:4) {
  stage = ii
  if(stage <= 2) {
    summarystage = "early"
  }else {
    summarystage = "late"
  }

  Z1 = data.frame(CXCL12 = CXCL12, tumorstage.2 = 0, tumorstage.3 = 0, tumorstage.4 = 0)
  if(stage == 2) { Z1$tumorstage.2 = 1 }
  if(stage == 3) { Z1$tumorstage.3 = 1 }
  if(stage == 4) { Z1$tumorstage.4 = 1 }

  Z2 = data.frame(CXCL12 = CXCL12, summarystage.late = 0)
  if(summarystage == "late") { Z2$summarystage.late = 1 }

  A = jointCox.Weibull.prediction(object = joint, time = time, p = p,
                                  Z1 = Z1, Z2 = Z2)
  res = A$MRL.of.X
  lines(res$t, res$Estimate, lwd = 2, lty = lty0[ii], col = col0[ii])
}
legend("topright", legend = c("stage I", "stage II", "stage III", "stage IV"),
       lwd = 2, lty = lty0, col = col0)


## SurvX ##
plot(0, 0, type = "l", xlim = c(0, 6420), ylim = c(0, 1),
     xlab = "time (days)", ylab = "Survival probability for TTP", main = paste("CXCL12 = ", CXCL12))
for (ii in 1:4) {
  stage = ii
  if(stage <= 2) {
    summarystage = "early"
  }else {
    summarystage = "late"
  }

  Z1 = data.frame(CXCL12 = CXCL12, tumorstage.2 = 0, tumorstage.3 = 0, tumorstage.4 = 0)
  if(stage == 2) { Z1$tumorstage.2 = 1 }
  if(stage == 3) { Z1$tumorstage.3 = 1 }
  if(stage == 4) { Z1$tumorstage.4 = 1 }

  Z2 = data.frame(CXCL12 = CXCL12, summarystage.late = 0)
  if(summarystage == "late") { Z2$summarystage.late = 1 }

  A = jointCox.Weibull.prediction(object = joint, time = time, p = p,
                                  Z1 = Z1, Z2 = Z2)
  res = A$Survival.of.X
  lines(res$t, res$Estimate, lwd = 2, lty = lty0[ii], col = col0[ii])
}
legend("topright", legend = c("stage I", "stage II", "stage III", "stage IV"),
       lwd = 2, lty = lty0, col = col0)


## HazardX ##
plot(0, 0, type = "l", xlim = c(0, 6420), ylim = c(0, 0.002),
     xlab = "time (days)", ylab = "Hazard rate for TTP", main = paste("CXCL12 = ", CXCL12))
for (ii in 1:4) {
  stage = ii
  if(stage <= 2) {
    summarystage = "early"
  }else {
    summarystage = "late"
  }

  Z1 = data.frame(CXCL12 = CXCL12, tumorstage.2 = 0, tumorstage.3 = 0, tumorstage.4 = 0)
  if(stage == 2) { Z1$tumorstage.2 = 1 }
  if(stage == 3) { Z1$tumorstage.3 = 1 }
  if(stage == 4) { Z1$tumorstage.4 = 1 }

  Z2 = data.frame(CXCL12 = CXCL12, summarystage.late = 0)
  if(summarystage == "late") { Z2$summarystage.late = 1 }

  A = jointCox.Weibull.prediction(object = joint, time = time, p = p,
                                  Z1 = Z1, Z2 = Z2)
  res = A$Hazard.of.X
  lines(res$t, res$Estimate, lwd = 2, lty = lty0[ii], col = col0[ii])
}
legend("topright", legend = c("stage I", "stage II", "stage III", "stage IV"),
       lwd = 2, lty = lty0, col = col0)


## OS ## =============================================================================
## MRLD ##
plot(0, 0, type = "l", xlim = c(0, 6420), ylim = c(0, 10000),
     xlab = "time (days)", ylab = "MRL for OS", main = paste("CXCL12 = ", CXCL12))
for (ii in c(1, 3)) {
  stage = ii
  if(stage <= 2) {
    summarystage = "early"
  }else {
    summarystage = "late"
  }

  Z1 = data.frame(CXCL12 = CXCL12, tumorstage.2 = 0, tumorstage.3 = 0, tumorstage.4 = 0)
  if(stage == 2) { Z1$tumorstage.2 = 1 }
  if(stage == 3) { Z1$tumorstage.3 = 1 }
  if(stage == 4) { Z1$tumorstage.4 = 1 }

  Z2 = data.frame(CXCL12 = CXCL12, summarystage.late = 0)
  if(summarystage == "late") { Z2$summarystage.late = 1 }

  A = jointCox.Weibull.prediction(object = joint, time = time, p = p,
                                  Z1 = Z1, Z2 = Z2)
  res = A$MRL.of.D
  if(ii == 1) {
    lines(res$t, res$Estimate, lwd = 2, lty = 1, col = 4)
  }
  if(ii == 3) {
    lines(res$t, res$Estimate, lwd = 2, lty = 2, col = 2)
  }
}
legend("topright", legend = c("stage early", "stage late"),
       lwd = 2, lty = c(1, 2), col = c(4, 2))



## SurvD ##
plot(0, 0, type = "l", xlim = c(0, 6420), ylim = c(0, 1),
     xlab = "time (days)", ylab = "Survival probability for OS", main = paste("CXCL12 = ", CXCL12))
for (ii in c(1, 3)) {
  stage = ii
  if(stage <= 2) {
    summarystage = "early"
  }else {
    summarystage = "late"
  }

  Z1 = data.frame(CXCL12 = CXCL12, tumorstage.2 = 0, tumorstage.3 = 0, tumorstage.4 = 0)
  if(stage == 2) { Z1$tumorstage.2 = 1 }
  if(stage == 3) { Z1$tumorstage.3 = 1 }
  if(stage == 4) { Z1$tumorstage.4 = 1 }

  Z2 = data.frame(CXCL12 = CXCL12, summarystage.late = 0)
  if(summarystage == "late") { Z2$summarystage.late = 1 }

  A = jointCox.Weibull.prediction(object = joint, time = time, p = p,
                                  Z1 = Z1, Z2 = Z2)
  res = A$Survival.of.D
  if(ii == 1) {
    lines(res$t, res$Estimate, lwd = 2, lty = 1, col = 4)
  }
  if(ii == 3) {
    lines(res$t, res$Estimate, lwd = 2, lty = 2, col = 2)
  }
}
legend("topright", legend = c("stage early", "stage late"),
       lwd = 2, lty = c(1, 2), col = c(4, 2))


## HazardD ##
plot(0, 0, type = "l", xlim = c(0, 6420), ylim = c(0, 0.002),
     xlab = "time (days)", ylab = "Hazard rate for OS", main = paste("CXCL12 = ", CXCL12))
for (ii in c(1, 3)) {
  stage = ii
  if(stage <= 2) {
    summarystage = "early"
  }else {
    summarystage = "late"
  }

  Z1 = data.frame(CXCL12 = CXCL12, tumorstage.2 = 0, tumorstage.3 = 0, tumorstage.4 = 0)
  if(stage == 2) { Z1$tumorstage.2 = 1 }
  if(stage == 3) { Z1$tumorstage.3 = 1 }
  if(stage == 4) { Z1$tumorstage.4 = 1 }

  Z2 = data.frame(CXCL12 = CXCL12, summarystage.late = 0)
  if(summarystage == "late") { Z2$summarystage.late = 1 }

  A = jointCox.Weibull.prediction(object = joint, time = time, p = p,
                                  Z1 = Z1, Z2 = Z2)
  res = A$Hazard.of.D
  if(ii == 1) {
    lines(res$t, res$Estimate, lwd = 2, lty = 1, col = 4)
  }
  if(ii == 3) {
    lines(res$t, res$Estimate, lwd = 2, lty = 2, col = 2)
  }
}
legend("topright", legend = c("stage early", "stage late"),
       lwd = 2, lty = c(1, 2), col = c(4, 2))
## Figure 11, 12 ## =====================================================================










