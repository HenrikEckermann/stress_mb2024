library(readxl)
library(car)
library(here)

# Read data #####
path <- here("data/SMILEY_maternal_stress_all.xlsx")
df_maternal_stress <- read_xlsx(path)

# Define missingness ####
df_maternal_stress[df_maternal_stress == 99999] <- NA

# Recode data ####
# STAI does not need to be reversed (recoding included in surveys answer options)
# EPDS needs reverse coding for item 3,5,6,7,8,9,10
df_maternal_stress$EPDS_18w_3_r <- recode(df_maternal_stress$EPDS_18w_3, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_18w_5_r <- recode(df_maternal_stress$EPDS_18w_5, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_18w_6_r <- recode(df_maternal_stress$EPDS_18w_6, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_18w_7_r <- recode(df_maternal_stress$EPDS_18w_7, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_18w_8_r <- recode(df_maternal_stress$EPDS_18w_8, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_18w_9_r <- recode(df_maternal_stress$EPDS_18w_9, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_18w_10_r <- recode(df_maternal_stress$EPDS_18w_10, '0=3;1=2;2=1;3=0')

df_maternal_stress$EPDS_32w_3_r <- recode(df_maternal_stress$EPDS_32w_3, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_32w_5_r <- recode(df_maternal_stress$EPDS_32w_5, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_32w_6_r <- recode(df_maternal_stress$EPDS_32w_6, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_32w_7_r <- recode(df_maternal_stress$EPDS_32w_7, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_32w_8_r <- recode(df_maternal_stress$EPDS_32w_8, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_32w_9_r <- recode(df_maternal_stress$EPDS_32w_9, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_32w_10_r <- recode(df_maternal_stress$EPDS_32w_10, '0=3;1=2;2=1;3=0')

df_maternal_stress$EPDS_8mo_3_r <- recode(df_maternal_stress$EPDS_8mo_3, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_8mo_5_r <- recode(df_maternal_stress$EPDS_8mo_5, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_8mo_6_r <- recode(df_maternal_stress$EPDS_8mo_6, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_8mo_7_r <- recode(df_maternal_stress$EPDS_8mo_7, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_8mo_8_r <- recode(df_maternal_stress$EPDS_8mo_8, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_8mo_9_r <- recode(df_maternal_stress$EPDS_8mo_9, '0=3;1=2;2=1;3=0')
df_maternal_stress$EPDS_8mo_10_r <- recode(df_maternal_stress$EPDS_8mo_10, '0=3;1=2;2=1;3=0')

#PSS-10 needs reverse coding for item 4,5,7,8.
df_maternal_stress$PSS_18w_4_r <- recode(df_maternal_stress$PSS_18w_4, '0=4;1=3;2=2;3=1;4=0')
df_maternal_stress$PSS_18w_5_r <- recode(df_maternal_stress$PSS_18w_5, '0=4;1=3;2=2;3=1;4=0')
df_maternal_stress$PSS_18w_7_r <- recode(df_maternal_stress$PSS_18w_7, '0=4;1=3;2=2;3=1;4=0')
df_maternal_stress$PSS_18w_8_r <- recode(df_maternal_stress$PSS_18w_8, '0=4;1=3;2=2;3=1;4=0')

df_maternal_stress$PSS_32w_4_r <- recode(df_maternal_stress$PSS_32w_4, '0=4;1=3;2=2;3=1;4=0')
df_maternal_stress$PSS_32w_5_r <- recode(df_maternal_stress$PSS_32w_5, '0=4;1=3;2=2;3=1;4=0')
df_maternal_stress$PSS_32w_7_r <- recode(df_maternal_stress$PSS_32w_7, '0=4;1=3;2=2;3=1;4=0')
df_maternal_stress$PSS_32w_8_r <- recode(df_maternal_stress$PSS_32w_8, '0=4;1=3;2=2;3=1;4=0')

df_maternal_stress$PSS_8mo_4_r <- recode(df_maternal_stress$PSS_8mo_4, '0=4;1=3;2=2;3=1;4=0')
df_maternal_stress$PSS_8mo_5_r <- recode(df_maternal_stress$PSS_8mo_5, '0=4;1=3;2=2;3=1;4=0')
df_maternal_stress$PSS_8mo_7_r <- recode(df_maternal_stress$PSS_8mo_7, '0=4;1=3;2=2;3=1;4=0')
df_maternal_stress$PSS_8mo_8_r <- recode(df_maternal_stress$PSS_8mo_8, '0=4;1=3;2=2;3=1;4=0')

#PES-BRIEF needs to be recoded when 'non applicable' is answered 
df_maternal_stress$PES_18w_11_r <- recode(df_maternal_stress$PES_18w_11, '4=0')
df_maternal_stress$PES_18w_12_r <- recode(df_maternal_stress$PES_18w_12, '4=0')
df_maternal_stress$PES_18w_13_r <- recode(df_maternal_stress$PES_18w_13, '4=0')
df_maternal_stress$PES_18w_14_r <- recode(df_maternal_stress$PES_18w_14, '4=0')
df_maternal_stress$PES_18w_15_r <- recode(df_maternal_stress$PES_18w_15, '4=0')
df_maternal_stress$PES_18w_16_r <- recode(df_maternal_stress$PES_18w_16, '4=0')
df_maternal_stress$PES_18w_17_r <- recode(df_maternal_stress$PES_18w_17, '4=0')
df_maternal_stress$PES_18w_18_r <- recode(df_maternal_stress$PES_18w_18, '4=0')
df_maternal_stress$PES_18w_19_r <- recode(df_maternal_stress$PES_18w_19, '4=0')
df_maternal_stress$PES_18w_20_r <- recode(df_maternal_stress$PES_18w_20, '4=0')

df_maternal_stress$PES_32w_11_r <- recode(df_maternal_stress$PES_32w_11, '4=0')
df_maternal_stress$PES_32w_12_r <- recode(df_maternal_stress$PES_32w_12, '4=0')
df_maternal_stress$PES_32w_13_r <- recode(df_maternal_stress$PES_32w_13, '4=0')
df_maternal_stress$PES_32w_14_r <- recode(df_maternal_stress$PES_32w_14, '4=0')
df_maternal_stress$PES_32w_15_r <- recode(df_maternal_stress$PES_32w_15, '4=0')
df_maternal_stress$PES_32w_16_r <- recode(df_maternal_stress$PES_32w_16, '4=0')
df_maternal_stress$PES_32w_17_r <- recode(df_maternal_stress$PES_32w_17, '4=0')
df_maternal_stress$PES_32w_18_r <- recode(df_maternal_stress$PES_32w_18, '4=0')
df_maternal_stress$PES_32w_19_r <- recode(df_maternal_stress$PES_32w_19, '4=0')
df_maternal_stress$PES_32w_20_r <- recode(df_maternal_stress$PES_32w_20, '4=0')

#PRAQ does not need to be reversed
#PSAS and PSAS-RSF-C does not need to be reversed (recoding included in surveys 
#answer options). Answer '5' (=not applicable) has been recoded to '0' manually in excel

# Sumscores ####
#STAI sumscores
df_maternal_stress$STAI_18w_sum <- df_maternal_stress$STAI_18w_1+df_maternal_stress$STAI_18w_2+
                                   df_maternal_stress$STAI_18w_3+df_maternal_stress$STAI_18w_4+
                                   df_maternal_stress$STAI_18w_5+df_maternal_stress$STAI_18w_6+
                                   df_maternal_stress$STAI_18w_7+df_maternal_stress$STAI_18w_8+
                                   df_maternal_stress$STAI_18w_9+df_maternal_stress$STAI_18w_10+
                                   df_maternal_stress$STAI_18w_11+df_maternal_stress$STAI_18w_12+
                                   df_maternal_stress$STAI_18w_13+df_maternal_stress$STAI_18w_14+
                                   df_maternal_stress$STAI_18w_15+df_maternal_stress$STAI_18w_16+
                                   df_maternal_stress$STAI_18w_17+df_maternal_stress$STAI_18w_18+
                                   df_maternal_stress$STAI_18w_19+df_maternal_stress$STAI_18w_20

df_maternal_stress$STAI_32w_sum <- df_maternal_stress$STAI_32w_1+df_maternal_stress$STAI_32w_2+
                                   df_maternal_stress$STAI_32w_3+df_maternal_stress$STAI_32w_4+
                                   df_maternal_stress$STAI_32w_5+df_maternal_stress$STAI_32w_6+
                                   df_maternal_stress$STAI_32w_7+df_maternal_stress$STAI_32w_8+
                                   df_maternal_stress$STAI_32w_9+df_maternal_stress$STAI_32w_10+
                                   df_maternal_stress$STAI_32w_11+df_maternal_stress$STAI_32w_12+
                                   df_maternal_stress$STAI_32w_13+df_maternal_stress$STAI_32w_14+
                                   df_maternal_stress$STAI_32w_15+df_maternal_stress$STAI_32w_16+
                                   df_maternal_stress$STAI_32w_17+df_maternal_stress$STAI_32w_18+
                                   df_maternal_stress$STAI_32w_19+df_maternal_stress$STAI_32w_20

df_maternal_stress$STAI_8mo_sum <- df_maternal_stress$STAI_8mo_1+df_maternal_stress$STAI_8mo_2+
                                   df_maternal_stress$STAI_8mo_3+df_maternal_stress$STAI_8mo_4+
                                   df_maternal_stress$STAI_8mo_5+df_maternal_stress$STAI_8mo_6+
                                   df_maternal_stress$STAI_8mo_7+df_maternal_stress$STAI_8mo_8+
                                   df_maternal_stress$STAI_8mo_9+df_maternal_stress$STAI_8mo_10+
                                   df_maternal_stress$STAI_8mo_11+df_maternal_stress$STAI_8mo_12+
                                   df_maternal_stress$STAI_8mo_13+df_maternal_stress$STAI_8mo_14+
                                   df_maternal_stress$STAI_8mo_15+df_maternal_stress$STAI_8mo_16+
                                   df_maternal_stress$STAI_8mo_17+df_maternal_stress$STAI_8mo_18+
                                   df_maternal_stress$STAI_8mo_19+df_maternal_stress$STAI_8mo_20

#EPDS sumscores
df_maternal_stress$EPDS_18w_sum <- df_maternal_stress$EPDS_18w_1+df_maternal_stress$EPDS_18w_2+
                                   df_maternal_stress$EPDS_18w_3_r+df_maternal_stress$EPDS_18w_4+
                                   df_maternal_stress$EPDS_18w_5_r+df_maternal_stress$EPDS_18w_6_r+
                                   df_maternal_stress$EPDS_18w_7_r+df_maternal_stress$EPDS_18w_8_r+
                                   df_maternal_stress$EPDS_18w_9_r+df_maternal_stress$EPDS_18w_10_r

df_maternal_stress$EPDS_32w_sum <- df_maternal_stress$EPDS_32w_1+df_maternal_stress$EPDS_32w_2+
                                   df_maternal_stress$EPDS_32w_3_r+df_maternal_stress$EPDS_32w_4+
                                   df_maternal_stress$EPDS_32w_5_r+df_maternal_stress$EPDS_32w_6_r+
                                   df_maternal_stress$EPDS_32w_7_r+df_maternal_stress$EPDS_32w_8_r+
                                   df_maternal_stress$EPDS_32w_9_r+df_maternal_stress$EPDS_32w_10_r

df_maternal_stress$EPDS_8mo_sum <- df_maternal_stress$EPDS_8mo_1+df_maternal_stress$EPDS_8mo_2+
                                   df_maternal_stress$EPDS_8mo_3_r+df_maternal_stress$EPDS_8mo_4+
                                   df_maternal_stress$EPDS_8mo_5_r+df_maternal_stress$EPDS_8mo_6_r+
                                   df_maternal_stress$EPDS_8mo_7_r+df_maternal_stress$EPDS_8mo_8_r+
                                   df_maternal_stress$EPDS_8mo_9_r+df_maternal_stress$EPDS_8mo_10_r

#PSS-10 sumscores
df_maternal_stress$PSS_18w_sum <- df_maternal_stress$PSS_18w_1+df_maternal_stress$PSS_18w_2+
                                  df_maternal_stress$PSS_18w_3+df_maternal_stress$PSS_18w_4_r+
                                  df_maternal_stress$PSS_18w_5_r+df_maternal_stress$PSS_18w_6+
                                  df_maternal_stress$PSS_18w_7_r+df_maternal_stress$PSS_18w_8_r+
                                  df_maternal_stress$PSS_18w_9+df_maternal_stress$PSS_18w_10

df_maternal_stress$PSS_32w_sum <- df_maternal_stress$PSS_32w_1+df_maternal_stress$PSS_32w_2+
                                  df_maternal_stress$PSS_32w_3+df_maternal_stress$PSS_32w_4_r+
                                  df_maternal_stress$PSS_32w_5_r+df_maternal_stress$PSS_32w_6+
                                  df_maternal_stress$PSS_32w_7_r+df_maternal_stress$PSS_32w_8_r+
                                  df_maternal_stress$PSS_32w_9+df_maternal_stress$PSS_32w_10

df_maternal_stress$PSS_8mo_sum <- df_maternal_stress$PSS_8mo_1+df_maternal_stress$PSS_8mo_2+
                                  df_maternal_stress$PSS_8mo_3+df_maternal_stress$PSS_8mo_4_r+
                                  df_maternal_stress$PSS_8mo_5_r+df_maternal_stress$PSS_8mo_6+
                                  df_maternal_stress$PSS_8mo_7_r+df_maternal_stress$PSS_8mo_8_r+
                                  df_maternal_stress$PSS_8mo_9+df_maternal_stress$PSS_8mo_10

# PES-BRIEF Intensity score (sum of intensity of experience of each hassle)
df_maternal_stress$PESbrief_18w_sum <- df_maternal_stress$PES_18w_11_r+df_maternal_stress$PES_18w_12_r+
                                  df_maternal_stress$PES_18w_13_r+df_maternal_stress$PES_18w_14_r+
                                  df_maternal_stress$PES_18w_15_r+df_maternal_stress$PES_18w_16_r+
                                  df_maternal_stress$PES_18w_17_r+df_maternal_stress$PES_18w_18_r+
                                  df_maternal_stress$PES_18w_19_r+df_maternal_stress$PES_18w_20_r

df_maternal_stress$PESbrief_32w_sum <- df_maternal_stress$PES_32w_11_r+df_maternal_stress$PES_32w_12_r+
                                  df_maternal_stress$PES_32w_13_r+df_maternal_stress$PES_32w_14_r+
                                  df_maternal_stress$PES_32w_15_r+df_maternal_stress$PES_32w_16_r+
                                  df_maternal_stress$PES_32w_17_r+df_maternal_stress$PES_32w_18_r+
                                  df_maternal_stress$PES_32w_19_r+df_maternal_stress$PES_32w_20_r

# PRAQR Sumscore of subscales: fear of handicapped child + fear of giving birth
df_maternal_stress$PRAQR_handicap_18w_sum <- df_maternal_stress$PRAQR_18w_1+df_maternal_stress$PRAQR_18w_2+
                                 df_maternal_stress$PRAQR_18w_3+df_maternal_stress$PRAQR_18w_5

df_maternal_stress$PRAQR_birth_18w_sum <- df_maternal_stress$PRAQR_18w_4+df_maternal_stress$PRAQR_18w_6+
                                 df_maternal_stress$PRAQR_18w_7

df_maternal_stress$PRAQR_handicap_32w_sum <- df_maternal_stress$PRAQR_32w_1+df_maternal_stress$PRAQR_32w_2+
                                 df_maternal_stress$PRAQR_32w_3+df_maternal_stress$PRAQR_32w_5

df_maternal_stress$PRAQR_birth_32w_sum <- df_maternal_stress$PRAQR_32w_4+df_maternal_stress$PRAQR_32w_6+
                                df_maternal_stress$PRAQR_32w_7

# PSAS and PSAS-RSF sumscores
df_maternal_stress$PSAS_RSFC_2w_sum <- df_maternal_stress$PSAS_2w_RSFC1+df_maternal_stress$PSAS_2w_RSFC2+
  df_maternal_stress$PSAS_2w_RSFC3+df_maternal_stress$PSAS_2w_RSFC4+df_maternal_stress$PSAS_2w_RSFC5+
  df_maternal_stress$PSAS_2w_RSFC6+df_maternal_stress$PSAS_2w_RSFC7+df_maternal_stress$PSAS_2w_RSFC8+
  df_maternal_stress$PSAS_2w_RSFC9+df_maternal_stress$PSAS_2w_RSFC10+df_maternal_stress$PSAS_2w_RSFC11+
  df_maternal_stress$PSAS_2w_RSFC12

df_maternal_stress$PSAS_6w_sum <- rowSums(df_maternal_stress[,c("PSAS_6w_1","PSAS_6w_2","PSAS_6w_3","PSAS_6w_4","PSAS_6w_5","PSAS_6w_6","PSAS_6w_7","PSAS_6w_8","PSAS_6w_9","PSAS_6w_10",
                                      "PSAS_6w_11","PSAS_6w_12","PSAS_6w_13","PSAS_6w_14","PSAS_6w_15","PSAS_6w_16","PSAS_6w_17","PSAS_6w_18","PSAS_6w_19","PSAS_6w_20",
                                      "PSAS_6w_21","PSAS_6w_22","PSAS_6w_23","PSAS_6w_24","PSAS_6w_25","PSAS_6w_26","PSAS_6w_27","PSAS_6w_28","PSAS_6w_29","PSAS_6w_30",
                                      "PSAS_6w_31","PSAS_6w_32","PSAS_6w_33","PSAS_6w_34","PSAS_6w_35","PSAS_6w_36","PSAS_6w_37","PSAS_6w_38","PSAS_6w_39","PSAS_6w_40",
                                      "PSAS_6w_41","PSAS_6w_42","PSAS_6w_43","PSAS_6w_44","PSAS_6w_45","PSAS_6w_46","PSAS_6w_47","PSAS_6w_48","PSAS_6w_49","PSAS_6w_50",
                                      "PSAS_6w_51")], na.rm = F)

df_maternal_stress$PSAS_RSFC_12w_sum <- df_maternal_stress$PSAS_12w_RSFC1+df_maternal_stress$PSAS_12w_RSFC2+
  df_maternal_stress$PSAS_12w_RSFC3+df_maternal_stress$PSAS_12w_RSFC4+df_maternal_stress$PSAS_12w_RSFC5+
  df_maternal_stress$PSAS_12w_RSFC6+df_maternal_stress$PSAS_12w_RSFC7+df_maternal_stress$PSAS_12w_RSFC8+
  df_maternal_stress$PSAS_12w_RSFC9+df_maternal_stress$PSAS_12w_RSFC10+df_maternal_stress$PSAS_12w_RSFC11+
  df_maternal_stress$PSAS_12w_RSFC12

# Check sumscores
hist(df_maternal_stress$EPDS_18w_sum)
hist(df_maternal_stress$EPDS_32w_sum)
hist(df_maternal_stress$EPDS_8mo_sum)

hist(df_maternal_stress$STAI_18w_sum)
hist(df_maternal_stress$STAI_32w_sum)
hist(df_maternal_stress$STAI_8mo_sum)

hist(df_maternal_stress$PSS_18w_sum)
hist(df_maternal_stress$PSS_32w_sum)
hist(df_maternal_stress$PSS_8mo_sum)

hist(df_maternal_stress$PESbrief_18w_sum)
hist(df_maternal_stress$PESbrief_32w_sum)

hist(df_maternal_stress$PRAQR_handicap_18w_sum)
hist(df_maternal_stress$PRAQR_handicap_32w_sum)
hist(df_maternal_stress$PRAQR_birth_18w_sum)
hist(df_maternal_stress$PRAQR_birth_32w_sum)

summary(df_maternal_stress$PESbrief_18w_sum)
summary(df_maternal_stress$PESbrief_32w_sum)

hist(df_maternal_stress$PSAS_RSFC_2w_sum)
hist(df_maternal_stress$PSAS_6w_sum)
hist(df_maternal_stress$PSAS_RSFC_12w_sum)

summary(df_maternal_stress$PSAS_RSFC_2w_sum)
summary(df_maternal_stress$PSAS_6w_sum)
summary(df_maternal_stress$PSAS_RSFC_12w_sum)

# export to henriks project
library(tidyverse)
ms <- select(df_maternal_stress, ID, contains("_sum")) %>%
  select(-contains("PSAS")) %>%
  pivot_longer(-ID, names_to = c("Q", "week"), names_pattern = "(\\w+)_(\\d+)[wm]o?_sum") %>%
  pivot_wider(names_from = Q, values_from = value) %>%
  mutate(
    pre = ifelse(week == "8", FALSE, TRUE),
    week = ifelse(week == "8", 32, as.numeric(week)))
colnames(ms) <- str_to_lower(colnames(ms))

# store raw data before going on
path <- here("data/ms_raw.csv")
write_csv(ms, file = path)
mlr::summarizeColumns(select(ms, stai, epds, pss))
colnames(ms)
ms <-   mutate(
  ms, 
  stai = scale(stai)[, 1],
  epds = scale(epds)[, 1],
  pss = scale(pss)[, 1],
  ms = stai + epds + pss,
  pesbrief = scale(pesbrief)[, 1],
  praqr_handicap = scale(praqr_handicap)[, 1],
  praqr_birth = scale(praqr_birth)[, 1]
)



psas <- select(
    df_maternal_stress, ID, 
    psas_2w_sum = PSAS_RSFC_2w_sum, 
    psas_6w_sum = PSAS_6w_sum, 
    psas_12w_sum = PSAS_RSFC_12w_sum
    ) %>%
  pivot_longer(-ID, names_to = c("Q", "week"), names_pattern = "(\\w+)_(\\d+)[wm]o?_sum") %>%
  pivot_wider(names_from = Q, values_from = value) %>%
  mutate(
    pre = FALSE,
    week =  as.numeric(week))

mlr::summarizeColumns(select(ms, ms))
path <- here("data/ms.csv")
write_csv(ms, file = path)
path <- here("data/psas.csv")
write_csv(psas, file = path)

# PSAS at 2, 6 & 12 weeks , PRAQ at 18 and 32 weeks
