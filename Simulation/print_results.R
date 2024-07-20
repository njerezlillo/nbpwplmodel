require(dplyr)
require(xtable)
require(numDeriv)
source("nbpwpl.R")
load("./Application/preliminar2.RData")

# Case 1 ------------------------------------------------------------------

load("./Simulation/sim_1.RData")

# Estimates

Table1 <- data.frame(estimates[[1]]$Beta[, 1:3],
                     estimates[[2]]$Beta[, 1:3],
                     estimates[[3]]$Beta[, 1:3])

Table2 <- data.frame(estimates[[1]]$Alpha[, 1:3],
                     estimates[[2]]$Alpha[, 1:3],
                     estimates[[3]]$Alpha[, 1:3])

Table1
Table2

xtable(Table1, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

xtable(Table2, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

# Monte Carlo

Table1 <- data.frame(estimates[[1]]$Beta[, 4:6],
                     estimates[[2]]$Beta[, 4:6],
                     estimates[[3]]$Beta[, 4:6])

Table2 <- data.frame(estimates[[1]]$Alpha[, 4:6],
                     estimates[[2]]$Alpha[, 4:6],
                     estimates[[3]]$Alpha[, 4:6])

Table1
Table2

xtable(Table1, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

xtable(Table2, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

# Case 2 ------------------------------------------------------------------

load("./Simulation/sim_2.RData")

# Estimates

Table1 <- data.frame(estimates[[1]]$Beta[, 1:3],
                     estimates[[2]]$Beta[, 1:3],
                     estimates[[3]]$Beta[, 1:3])

Table2 <- data.frame(estimates[[1]]$Alpha[, 1:3],
                     estimates[[2]]$Alpha[, 1:3],
                     estimates[[3]]$Alpha[, 1:3])

Table1
Table2

xtable(Table1, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

xtable(Table2, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

# Monte Carlo

Table1 <- data.frame(estimates[[1]]$Beta[, 4:6],
                     estimates[[2]]$Beta[, 4:6],
                     estimates[[3]]$Beta[, 4:6])

Table2 <- data.frame(estimates[[1]]$Alpha[, 4:6],
                     estimates[[2]]$Alpha[, 4:6],
                     estimates[[3]]$Alpha[, 4:6])

Table1
Table2

xtable(Table1, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

xtable(Table2, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

# Case 3 ------------------------------------------------------------------

load("./Simulation/sim_3.RData")

# Estimates

Table1 <- data.frame(estimates[[1]]$Beta[, 1:3],
                     estimates[[2]]$Beta[, 1:3],
                     estimates[[3]]$Beta[, 1:3])

Table2 <- data.frame(estimates[[1]]$Alpha[, 1:3],
                     estimates[[2]]$Alpha[, 1:3],
                     estimates[[3]]$Alpha[, 1:3])

Table1
Table2

xtable(Table1, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

xtable(Table2, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

# Monte Carlo

Table1 <- data.frame(estimates[[1]]$Beta[, 4:6],
                     estimates[[2]]$Beta[, 4:6],
                     estimates[[3]]$Beta[, 4:6])

Table2 <- data.frame(estimates[[1]]$Alpha[, 4:6],
                     estimates[[2]]$Alpha[, 4:6],
                     estimates[[3]]$Alpha[, 4:6])

Table1
Table2

xtable(Table1, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})

xtable(Table2, digits = 3) %>% 
  print(include.rownames = F, sanitize.text.function = function (x) {x})
