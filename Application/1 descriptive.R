library(dplyr)
library(survival)
library(survminer)
library(cuRe)
source("nbpwpl.R")
source("nbcfwei.R")

# Dataset -----------------------------------------------------------------

df <- read.csv("./Application/kidney.csv")[,-1] %>% na.omit() %>% 
  filter(Disease.Type %in% c("Adenomas and Adenocarcinomas",
                             "Complex Mixed and Stromal Neoplasms")) %>% 
  mutate(gender = case_when(Gender == "Female" ~ 1, Gender == "Male" ~ 0),
         disease = case_when(Disease.Type == "Adenomas and Adenocarcinomas" ~ 0,
                             Disease.Type == "Complex Mixed and Stromal Neoplasms" ~ 1),
         age = as.numeric(Age.at.diagnosis > mean(Age.at.diagnosis)),
         intercept = 1) %>% 
  select(time, status, intercept, age, gender, disease)

fit_km_disease <-
  survfit(Surv(time, status) ~ disease, data = df %>% mutate(time = time - 1))

q <- -1

# Figure 1 ----------------------------------------------------------------

ggsurvplot(fit_km_disease, palette = c("#E7B800", "#2E9FDF"),
           legend.labs = c("A&A",
                           "CM&SN."), 
           conf.int = T, ggtheme = theme_bw(),
           legend.title = "Disease", size = 0.5, censor = F) +
  labs(x = "Duration (years)", y = "Survival Rate")

ggsave("./Application/km_disease.pdf", width = 8, height = 5)

# Cure fraction weibull model ---------------------------------------------

l_wei <- function(w) loglik_nbpwwei(w, df[, -c(4, 5)], q)

fit_wei <-
  maxBFGS(
    l_wei,
    start = c(2, 2, 1, 1),
    constraints = list(
      ineqA = cbind(diag(c(1, 1)), matrix(0, ncol = 2, nrow = 2)),
      ineqB = c(0.01, 0.01)
    )
  )

theta <- fit_wei$estimate

s_AA <- Vectorize(function(x) snbpwwei(x + 1, q, plogis(theta[3]), theta[1:2]), "x")
s_CMSN <- Vectorize(function(x) snbpwwei(x + 1, q, plogis(sum(theta[3:4])), theta[1:2]), "x")

ggsurvplot(fit_km_disease, palette = rep("gray90", 2),
           conf.int = F, ggtheme = theme_bw(), size = 0.5,
           legend.title = "Disease", censor = F)$plot +
  labs(x = "Duration (years)", y = "Survival Rate") + 
  geom_function(fun = s_AA, aes(colour = "A&A")) +
  geom_function(fun = s_CMSN, aes(colour = "CM&SN")) +
  theme(legend.position = "top") + 
  scale_colour_manual("Disease", values = 
                       c("A&A" = "#FEB24C",
                         "CM&SN" = "#2C7FB8"))
ggsave("./Application/wei_cure_disease.pdf", width = 8, height = 5)

save.image("./Application/run_script_1.RData")
