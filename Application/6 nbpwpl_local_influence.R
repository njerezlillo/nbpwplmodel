library(stringr)
library(dplyr)
library(maxLik)
library(survival)
library(survminer)
library(xtable)
source("nbpwpl.R")

# Dataset -----------------------------------------------------------------

load("./Application/preliminar3.RData")
df <- arrange(df, time)

# Influence ---------------------------------------------------------------

I_cw <- Influence(df, coef_final, p, H, "cw")
I_rs <- Influence(df, coef_final, p, H, "rs")

# case-weight

fig_Bi_cw <- 
  ggplot(data.frame(id = df$time - 1, Bi = I_cw$Bi, di = df$status)) + 
  geom_point(aes(x = id, y = Bi, pch = di)) + theme_bw() +
  labs(x = "Time (i)", y = "Conformal curvature (Bi)") +
  scale_shape_manual("", values = c(1, 4),
                     labels = c(expression(delta[i] == 0), 
                                expression(delta[i] == 1))) +
  geom_segment(aes(x = 0, y = mean(I_cw$Bi) + 3.5 * sd(I_cw$Bi), 
                   xend = max(df$time) - 1, yend = mean(I_cw$Bi) + 3.5 * sd(I_cw$Bi)),
              col = "gray70", lty = 2) +
  annotate("text", label = "cut-off = 0.00587",
    x = 13.5, y = mean(I_cw$Bi) + 3.8 * sd(I_cw$Bi), 
    size = 4, colour = "black") +
  theme(legend.text = element_text(size = 12),
        legend.position = "top")

fig_dmax_cw <- 
  ggplot(data.frame(id = df$time - 1, Bi = I_cw$dmax)) + 
  geom_point(aes(x = id, y = Bi), pch = 1) + theme_bw() +
  labs(x = "Time (i)", y = "dmax (Bi)") +
  geom_segment(aes(x = p[2] - 1, y = 0,
                   xend = p[2] - 1, yend = max(I_cw$dmax)*1.01), 
               col = "gray70", lty = 2) +
  geom_segment(aes(x = p[3] - 1, y = 0,
                   xend = p[3] - 1, yend = max(I_cw$dmax)*1.01), 
               col = "gray70", lty = 2)

fig_Ui_cw <- 
  ggplot(data.frame(id = df$time - 1, Bi = I_cw$Ui)) + 
  geom_line(aes(x = id, y = Bi)) + theme_bw() +
  labs(x = "Time (i)", y = "(Ui)")

# response

fig_Bi_rs <- 
  ggplot(data.frame(id = df$time - 1, Bi = I_rs$Bi, di = df$status)) + 
  geom_point(aes(x = id, y = Bi, pch = di)) + theme_bw() +
  labs(x = "Time (i)", y = "Conformal curvature (Bi)") +
  scale_shape_manual("", values = c(1, 4),
                     labels = c(expression(delta[i] == 0), 
                                expression(delta[i] == 1))) +
  geom_segment(aes(x = 0, y = mean(I_rs$Bi) + 3.5 * sd(I_rs$Bi), 
                   xend = max(df$time) - 1, 
                   yend = mean(I_rs$Bi) + 3.5 * sd(I_rs$Bi)),
               col = "gray70", lty = 2) +
  annotate("text", label = "cut-off = 0.00647",
           x = 13.5, y = mean(I_rs$Bi) + 3.9 * sd(I_rs$Bi),
           size = 4, colour = "black") +
  theme(legend.text = element_text(size = 12),
        legend.position = "top")

fig_dmax_rs <- 
  ggplot(data.frame(id = df$time - 1, Bi = I_rs$dmax)) + 
  geom_point(aes(x = id, y = Bi), pch = 1) + theme_bw() +
  labs(x = "Time (i)", y = "dmax (Bi)") +
  geom_segment(aes(x = p[2] - 1, y = 0,
                   xend = p[2] - 1, yend = max(I_rs$dmax)*1.01), 
               col = "gray70", lty = 2) +
  geom_segment(aes(x = p[3] - 1, y = 0,
                   xend = p[3] - 1, yend = max(I_rs$dmax)*1.01), 
               col = "gray70", lty = 2)

fig_Ui_rs <- 
  ggplot(data.frame(id = df$time - 1, Bi = I_rs$Ui)) + 
  geom_line(aes(x = id, y = Bi)) + theme_bw() +
  labs(x = "Time (Year)", y = "(Ui)")

# out
plot_influence <- gridExtra::grid.arrange(fig_Bi_cw, fig_Bi_rs, ncol = 2)
ggsave(plot = plot_influence, "./Application/bi_influence.pdf",
       height = 5, width = 10)

# CW: post-deletion analysis ----------------------------------------------

index_cw <- which(I_cw$Bi > mean(I_cw$Bi) + 3.5 * sd(I_cw$Bi)) # 51
df_cw <- df[index_cw, ]

df_cw$risk_interval <- rep(c("R1", "R2", "R3"), n_each_interval(df_cw$time, p))

df_cw$group[df_cw$age == 0 & df_cw$disease == 0] <-  "Group I"
df_cw$group[df_cw$age == 0 & df_cw$disease == 1] <-  "Group II"
df_cw$group[df_cw$age == 1 & df_cw$disease == 1] <-  "Group III"

table(df_cw$risk_interval, df_cw$group, useNA = "ifany")

n_each_interval(df_cw$time[df_cw$group == "Group I"], p)
index_cw[which(df_cw$group == "Group I")]

index_each_interval(df_cw$time[df_cw$group == "Group II"], p)
index_cw[which(df_cw$group == "Group II")]

index_each_interval(df_cw$time[df_cw$group == "Group III"], p)
index_cw[which(df_cw$group == "Group III")]

# RS: post-deletion analysis ----------------------------------------------

index_rs <- which(I_rs$Bi > mean(I_rs$Bi) + 3.5 * sd(I_rs$Bi)) # 74
df_rs <- df[index_rs, ]

df_rs$risk_interval <- rep(c("R1", "R2", "R3"), n_each_interval(df_rs$time, p))

df_rs$group[df_rs$age == 0 & df_rs$disease == 0] <-  "Group I"
df_rs$group[df_rs$age == 0 & df_rs$disease == 1] <-  "Group II"
df_rs$group[df_rs$age == 1 & df_rs$disease == 1] <-  "Group III"

table(df_rs$risk_interval, df_rs$group, useNA = "ifany")

n_each_interval(df_rs$time[df_rs$group == "Group I"], p)
index_rs[which(df_rs$group == "Group I")]

index_each_interval(df_rs$time[df_rs$group == "Group II"], p)
index_rs[which(df_rs$group == "Group II")]

index_each_interval(df_rs$time[df_rs$group == "Group III"], p)
index_rs[which(df_rs$group == "Group III")]

####

intersect(index_rs, index_cw)

save.image("./Application/run_script_6.RData")