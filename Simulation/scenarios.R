library(ggplot2)
source("nbpwpl.R")
load("./Application/preliminar3.RData")

# Figure 1 ----------------------------------------------------------------

case_1 <- c(-0.74,  0.00,  0.00,
            0.67, 0.00, 0.00,
            -1.81, 0.00, 0.00,
            -5.71, 0.00, 4.81)

Alpha_1 <- alpha_matrix(case_1, df[,-6], p)

c_1_s_1 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_1[1, 4], Alpha_1[1, 1:3], p))
c_1_s_2 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_1[2, 4], Alpha_1[2, 1:3], p))
c_1_s_3 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_1[3, 4], Alpha_1[3, 1:3], p))
c_1_s_4 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_1[4, 4], Alpha_1[4, 1:3], p))

case_2 <- c(-0.74,  0.00,  0.00,
            0.67, 0.00, 0.00,
            -1.81, 0.00, 0.00,
            -5.71, -0.55, 4.81)

Alpha_2 <- alpha_matrix(case_2, df[,-6], p)

c_2_s_1 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_2[1, 4], Alpha_2[1, 1:3], p))
c_2_s_2 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_2[2, 4], Alpha_2[2, 1:3], p))
c_2_s_3 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_2[3, 4], Alpha_2[3, 1:3], p))
c_2_s_4 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_2[4, 4], Alpha_2[4, 1:3], p))

fig_1 <- 
  ggplot() +
  stat_function(fun = c_1_s_1, xlim = c(0, 12)) + 
  stat_function(fun = c_1_s_2, xlim = c(0, 12)) + 
  stat_function(fun = c_1_s_3, xlim = c(0, 12)) + 
  stat_function(fun = c_1_s_4, xlim = c(0, 12)) + 
  theme_bw() + labs(x = "Time", y = "Survival")
fig_2 <- 
  ggplot() +
  stat_function(fun = c_2_s_1, xlim = c(0, 12)) + 
  stat_function(fun = c_2_s_2, xlim = c(0, 12)) + 
  stat_function(fun = c_2_s_3, xlim = c(0, 12)) + 
  stat_function(fun = c_2_s_4, xlim = c(0, 12)) + 
  theme_bw() + labs(x = "Time", y = "Survival")

plot_scenarios<- gridExtra::grid.arrange(fig_1, fig_2, ncol = 2)
ggsave(plot = plot_scenarios, "./Simulation/scenarios.pdf", height = 5, width = 10)

