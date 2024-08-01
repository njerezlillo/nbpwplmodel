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

case_3 <- c(-0.74, 0.61,  0.00,
            0.67, 0.00, 0.00,
            -1.81, 0.00, 0.00,
            -5.71, 0.00, 4.81)

Alpha_3 <- alpha_matrix(case_3, df[,-6], p)

c_3_s_1 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_3[1, 4], Alpha_3[1, 1:3], p))
c_3_s_2 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_3[2, 4], Alpha_3[2, 1:3], p))
c_3_s_3 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_3[3, 4], Alpha_3[3, 1:3], p))
c_3_s_4 <- Vectorize(function(x) snbpwpl(x + 1, q, Alpha_3[4, 4], Alpha_3[4, 1:3], p))


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

fig_3 <- 
  ggplot() +
  stat_function(fun = c_3_s_1, xlim = c(0, 12)) + 
  stat_function(fun = c_3_s_2, xlim = c(0, 12)) + 
  stat_function(fun = c_3_s_3, xlim = c(0, 12)) + 
  stat_function(fun = c_3_s_4, xlim = c(0, 12)) + 
  theme_bw() + labs(x = "Time", y = "Survival")

plot_scenarios<- gridExtra::grid.arrange(fig_1, fig_2, fig_3, ncol = 3)
ggsave(plot = plot_scenarios, "./Simulation/FigureB1.eps", 
       height = 5, width = 12,
       dpi = 600, device = "eps", units = "in")
ggsave(plot = plot_scenarios, "./Simulation/FigureB1.tif",
       width = 12, height = 5,
       dpi = 600, device = "tiff")


