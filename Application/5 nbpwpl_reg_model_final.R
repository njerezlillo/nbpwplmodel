library(dplyr)
library(maxLik)
library(xtable)
library(stringr)
library(survival)
library(survminer)
source("nbpwpl.R")

# Dataset -----------------------------------------------------------------

load("./Application/preliminar2.RData")
n <- nrow(df)

# Final estimation --------------------------------------------------------

target <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df, p)
}

fit <- optim(fn = target, par = rep(0, 9), method = "CG", hessian = T)
round(2 * 9 + 2 * fit$value, 2) # AIC
round(log(nrow(df)) * 9 + 2 * fit$value, 2) # BIC
H <- solve(fit$hessian)
se <- H %>% diag %>% sqrt

# Theta -------------------------------------------------------------------

coef_final <- rep(0, 16) 
coef_final[!zero_position] <- fit$par

theta_matrix(coef_final, df, p)

results <- data.frame(est = fit$par,
                      se = se,
                      ci_l = fit$par - 1.96 * se,
                      ci_u = fit$par + 1.96 * se,
                      pval = 2 * pnorm(abs(fit$par/se), lower.tail = F))

print(xtable(t(results)), include.rownames = F)

# Alpha -------------------------------------------------------------------

Alpha <- alpha_matrix(coef_final, df, p)

# Group 1
se_1 <- deltamethod(coef_final, c(1, 0, 0, 0), H, p)
results_1 <- data.frame(est = Alpha[1, ],
                        se = se_1,
                        ci_l = Alpha[1, ] - 1.96 * se_1,
                        ci_u = Alpha[1, ] + 1.96 * se_1)
results_1

# Group 2
se_2 <- deltamethod(coef_final, c(1, 0, 0, 1), H, p)
results_2 <- data.frame(est = Alpha[4, ],
                        se = se_2,
                        ci_l = Alpha[4, ] - 1.96 * se_2,
                        ci_u = Alpha[4, ] + 1.96 * se_2)
results_2

# Group 3
se_3 <- deltamethod(coef_final, c(1, 1, 0, 1), H, p)
results_3 <- data.frame(est = Alpha[6, ],
                        se = se_3,
                        ci_l = Alpha[6, ] - 1.96 * se_3,
                        ci_u = Alpha[6, ] + 1.96 * se_3)
results_3

print(xtable(rbind(results_1, results_3, results_2)), include.rownames = F)

# Figure ------------------------------------------------------------------

s_1 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha[1, 4], Alpha[1, 1:3], p))
s_2 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha[4, 4], Alpha[4, 1:3], p))
s_3 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha[6, 4], Alpha[6, 1:3], p))

fit_km <- 
  survfit(Surv(time, status) ~ disease + age + gender,
          data = df[-c(39, 40, 50, 172),] %>% 
            mutate(time = time - 1))

ggsurvplot(fit_km, legend = "none", palette = rep("gray90", 9),
           ggtheme = theme_bw(), size = 0.5, censor = F)$plot +
  labs(x = "Duration (years)", y = "Survival Rate") +
  geom_function(fun = s_1, aes(colour = "Group I")) +
  geom_function(fun = s_3, aes(colour = "Group II")) +
  geom_function(fun = s_2, aes(colour = "Group III")) +
  theme(legend.position = "top") + 
  scale_color_manual(NULL, values = 
                       c("Group I" = "#E34A33", "Group III" = "#2C7FB8",
                         "Group II" = "#FEB24C"))

ggsave("./Application/Figure3.eps", width = 8, height = 5,
       dpi = 600, device = "eps", units = "in")
ggsave("./Application/Figure3.tif", width = 8, height = 5,
       dpi = 600, device = "tiff")

# RQR ---------------------------------------------------------------------

temp1 <- s_1(df$time[df$age == 0 & df$disease == 0] - 1)
temp2 <- s_2(df$time[df$age == 0 & df$disease == 1] - 1)
temp3 <- s_3(df$time[df$age == 1 & df$disease == 1] - 1)

rqr <- qnorm(c(temp1))

ggplot(data.frame(x = rqr), aes(sample = x)) +
  stat_qq() + labs(y = "Weight") +
  stat_qq_line() +
  theme_bw()

# Cox-snell ---------------------------------------------------------------

ei <- -c(log(temp1), log(temp2), log(temp3))
km_ei <- survfit(Surv(ei, c(df$status[df$age == 0 & df$disease == 0],
                            df$status[df$age == 0 & df$disease == 1],
                            df$status[df$age == 1 & df$disease == 1])) ~ 1)

ggplot() + aes(x = km_ei$surv, y = exp(-km_ei$time)) + 
  geom_point(pch = 1) + theme_bw() + 
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1.2) +
  labs(x = "S(ei): Kaplan-Meier", y = "S(ei): Standard Exponential")

ggsave("./Application/FigureF3.eps", width = 8, height = 4,
       dpi = 600, device = "eps", units = "in")
ggsave("./Application/FigureF3.tif", width = 8, height = 4,
       dpi = 600, device = "tiff")

save(df, p, q, zero_position, coef_final, H, Alpha, n,
     file = "./Application/preliminar3.RData")
save.image("./Application/run_script_5.RData")
