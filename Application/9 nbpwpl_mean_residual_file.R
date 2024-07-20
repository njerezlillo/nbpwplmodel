library(stringr)
library(dplyr)
library(maxLik)
library(xtable)
library(survival)
library(survminer)
source("nbpwpl.R")

# Dataset -----------------------------------------------------------------

load("./Application/preliminar3.RData")

Alpha <- alpha_matrix(coef_final, df, p)

s_1 <- Vectorize(function(x) snbpwpl(x + 1, q, p0_1, Alpha[1, 1:3], p))
s_2 <- Vectorize(function(x) snbpwpl(x + 1, q, p0_2, Alpha[4, 1:3], p))
s_3 <- Vectorize(function(x) snbpwpl(x + 1, q, p0_3, Alpha[6, 1:3], p))

ggplot() +
  geom_function(fun = s_1, aes(colour = "Group I")) +
  geom_function(fun = s_3, aes(colour = "Group II")) +
  geom_function(fun = s_2, aes(colour = "Group III")) +
  xlim(c(0, 15)) +
  scale_colour_manual(name = "Groups", 
                      values = c("Group I" = "#E34A33", 
                                 "Group II" = "#FEB24C", 
                                 "Group III" = "#2C7FB8")) +
  theme(legend.position = "top")

p0_1 <- Alpha[1, 4]
p0_2 <- Alpha[4, 4]
p0_3 <- Alpha[6, 4]

s_1 <- Vectorize(function(x) (snbpwpl(x + 1, q, p0_1, Alpha[1, 1:3], p) - p0_1) / (1 - p0_1))
s_2 <- Vectorize(function(x) (snbpwpl(x + 1, q, p0_2, Alpha[4, 1:3], p) - p0_2) / (1 - p0_2))
s_3 <- Vectorize(function(x) (snbpwpl(x + 1, q, p0_3, Alpha[6, 1:3], p) - p0_3) / (1 - p0_3))

ggplot() +
  geom_function(fun = s_1, aes(colour = "Group I")) +
  geom_function(fun = s_3, aes(colour = "Group II")) +
  geom_function(fun = s_2, aes(colour = "Group III")) +
  xlim(c(0, 15)) +
  scale_colour_manual(name = "Groups", 
                      values = c("Group I" = "#E34A33", 
                                 "Group II" = "#FEB24C", 
                                 "Group III" = "#2C7FB8")) +
  theme(legend.position = "top")

###

mrf_1 <- function (u) {
  integrate(function(x) s_1(x), u, Inf)$value / s_1(u)
}

mrf_2 <- function (u) {
  integrate(function(x) s_2(x), u, Inf)$value / s_2(u)
}

mrf_3 <- function (u) {
  integrate(function(x) s_3(x), u, Inf)$value / s_3(u)
}

mrf_1 <- Vectorize(mrf_1)
mrf_2 <- Vectorize(mrf_2)
mrf_3 <- Vectorize(mrf_3)

ggplot() +
  geom_function(fun = mrf_1, aes(colour = "Group I")) +
  geom_function(fun = mrf_3, aes(colour = "Group II")) +
  geom_function(fun = mrf_2, aes(colour = "Group III")) +
  xlim(c(0, 7)) +
  scale_colour_manual(
    name = "Groups",
    values = c(
      "Group I" = "#E34A33",
      "Group II" = "#FEB24C",
      "Group III" = "#2C7FB8"
    )
  ) +
  theme(legend.position = "top")

seq_time <- seq(0, 6, by = 1)

xtable(t(
  data.frame(Time = seq_time,
             "Group I" = mrf_1(seq_time),
             "Group II" = mrf_3(seq_time),
             "Group III" = mrf_2(seq_time))
))

f1 <- function(x) s_1(x)
f2 <- function(x) s_2(x)
f3 <- function(x) s_3(x)

f4 <- function(x) x * s_1(x)
f5 <- function(x) x * s_2(x)
f6 <- function(x) x * s_3(x)

sqrt(2 * integrate(Vectorize(f4), 0, Inf)$value -
       integrate(Vectorize(f1), 0, Inf)$value ^ 2)
sqrt(2 * integrate(Vectorize(f6), 0, Inf)$value -
       integrate(Vectorize(f3), 0, Inf)$value ^ 2)
sqrt(2 * integrate(Vectorize(f5), 0, Inf)$value -
       integrate(Vectorize(f2), 0, Inf)$value ^ 2)
