library(dplyr)

# Diagnosis ---------------------------------------------------------------

t <- 36 
dir_temp <- paste0('./Application/resources/', paste0("temp", 1:t, ".tsv"))

df_diagnosis <- NULL

for (i in 1:t) {
  df_temp <- read.delim(dir_temp[i], sep = '\t') %>% 
    rename(submitter_id = Case.ID) %>% select(-Vital.Status)
  df_diagnosis <- df_diagnosis %>% rbind(df_temp) 
}
  
# Gender ------------------------------------------------------------------

df_female <- 
  read.delim('./Application/resources/kidney_female.tsv', sep = '\t') %>% 
  mutate(time = time + 1, 
         status = 
           case_when(censored == "true" ~ F,
                     censored == "false" ~ T)
  )

df_male <- 
  read.delim('./Application/resources/kidney_male.tsv', sep = '\t') %>% 
  mutate(time = time + 1, 
         status = 
           case_when(censored == "true" ~ F,
                     censored == "false" ~ T)
  ) 

df_gender <- rbind(df_female, df_male) %>% 
  select(time, status, submitter_id)

df <- left_join(df_gender, df_diagnosis) %>% 
  mutate(Age.at.diagnosis = gsub( " .*$", "", Age.at.diagnosis),
         Age.at.diagnosis = as.numeric(Age.at.diagnosis),
         Age.at.diagnosis = ifelse(Age.at.diagnosis > 90, 0, Age.at.diagnosis))

# Out ---------------------------------------------------------------------

write.csv(df, "./Application/kidney.csv")

# Analysis ----------------------------------------------------------------

df$Gender %>% table

df$Disease.Type %>% table

df$Age.at.diagnosis %>% quantile(c(0.25, 0.5, 0.75), na.rm = T)

df$Primary.Diagnosis %>% table %>% sort 

df$Ethnicity %>% table

df$Race %>% table

df2 <- df %>% filter(Disease.Type %in% c("Adenomas and Adenocarcinomas",
                                  "Complex Mixed and Stromal Neoplasms"))

survfit(Surv(time, status) ~ Primary.Diagnosis + Disease.Type, 
        data = df2 %>% filter(Primary.Diagnosis %in% c("Clear cell adenocarcinoma, NOS", "Papillary adenocarcinoma, NOS"),
                             )) %>% plot



