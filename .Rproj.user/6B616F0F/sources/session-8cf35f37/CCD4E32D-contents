# for imputation
set.seed(1)

library(tidyverse)
library(readxl)
library(glue)
library(here)
library(mice)

path <- here("data")
files <- list.files(
  path, 
  full.names = TRUE, 
  include.dirs = FALSE, 
  recursive = TRUE,
  pattern = ".xlsx$"
  )


# FFQ
files <- files <- list.files(
  here("data/ffq"), 
  full.names = TRUE, 
  include.dirs = FALSE, 
  recursive = TRUE,
  pattern = ".xlsx$"
)

files
ffq <- map_dfr(files, function(file) {
  week <- str_extract(file, "(\\d\\d)wks.xlsx$", group = 1)
  read_excel(file) %>%
    mutate(
      ID = as.numeric(ID),
      week = as.numeric(week),
      pre = ifelse(
        week == 12, FALSE, ifelse(
          week == 18, TRUE, ifelse(
            week == 32, TRUE, NA))),
      week = ifelse(week == 12, 32, week) # for joining 
      )
})
colnames(ffq) <- str_to_lower(colnames(ffq))
ffq <- select(ffq, id, week, pre, contains("dhd"))


# STAI, EPDS & PSS

file <- here("data/ms.csv")
ms <- read_csv(file)


# parity 
parity <- read_excel(here("data/SMILEY_partner_parity.xlsx")) %>%
  select(id = ID, parity = firstchild)

# antibiotics mother
# files <- c(
#   "data/SMILEY_18wks_Stoolquestionnaire_Rawdata.xlsx",
#   "data/SMILEY_32wks_Stoolquestionnaire_rawdata.xlsx",
#   "data/SMILEY_8mo_FU_Stoolquestionnaire_Rawdata.xlsx"
# )
# sheets <- c("18wkn", "32wkn", "8MOFU_mother_NAME")
# abx <- map2_dfr(files, sheets, function(file, sheet) {
#   read_excel(file, sheet = sheet, na = "99999") %>%
#     select(id = ID, abx = Stool_6) %>% 
#     mutate(
#       id = as.numeric(id),
#       pre = ifelse(str_detect(sheet, "18"), TRUE, ifelse(str_detect(sheet, "32"), TRUE, FALSE)),
#       abx = ifelse(abx == 2, 0, ifelse(abx == 1, 1, ifelse(abx == 88888, 0, NA))),
#       week = ifelse(str_detect(sheet, "18"), 18, ifelse(str_detect(sheet, "32"), 32, 32))) %>%
#     filter(id != "EXAMPLE")
# })

file <- here("data/Maternal_antibiotics.xlsx")
cols <- c(
  "antibiotics_18wks_pregnancy",
  "antibiotics_32wks_pregnancy",
  "antibiotics_32wks_postpartum"
)

abx <- map_dfr(cols, function(col) {
  read_excel(file, na = "") %>%
    select(id = ID, abx = glue("{col}")) %>%
    mutate(
      id = as.numeric(id),
      pre = ifelse(str_detect(col, "pregnancy"), TRUE, ifelse(str_detect(col, "postpartum"), FALSE, NA)),
      abx = ifelse(abx == 2, 0, ifelse(abx == 1, 1, ifelse(abx == 88888, NA, NA))),
      week = ifelse(str_detect(col, "18"), 18, ifelse(str_detect(col, "32"), 32, NA))) 
})



count(abx, week, pre, abx)

# gestational week 
ga <- read_excel(here("data/SMILEY_gestational_age.xlsx")) %>%
  select(
    id = ID, 
    #ga_end = Gestational_age_birth_days,
    ga_18wks = Gestational_age_prenatal_18wks_weeks,
    ga_32wks = Gestational_age_prenatal_32wks_weeks) %>%
  pivot_longer(contains("ga"), names_to = "week", values_to = "ga", names_pattern = "(\\d\\d)") %>%
  mutate(week = as.numeric(week))

# maternal age
age <- read_excel(here("data/SMILEY_age_mother.xlsx")) %>%
  select(id = ID, age = Age_mother_prenatal_18wks) 
age


# pre-pregnancy BMI 
file <- list.files(here("data"), pattern = "BMI", full.names = TRUE)[[1]]
bmi <- read_excel(file, na = "99999", col_types = "numeric")
warnings()
colnames(bmi) <- c("id", "length", "weight")
bmi$pbmi <- bmi$weight / bmi$length^2


# Physical activity
files <- list.files(here("data/"), pattern = "PPAQ", full.names = TRUE)
activity <- map_dfr(files, function(file) {
  read_excel(file, na = "99999") %>%
    select(id = ID, PPAQ_20, PPAQ_21,PPAQ_22, PPAQ_23, PPAQ_24, PPAQ_25, 
           PPAQ_26, PPAQ_27a, PPAQ_27b, PPAQ_28a, PPAQ_28b) %>%
    mutate(
      id = as.numeric(id),
      week = str_extract(file, "(\\d+)wks", group = 1),
      PPAQ_27 = ifelse(!is.na(PPAQ_27a), PPAQ_27b, 1),
      PPAQ_28 = ifelse(!is.na(PPAQ_28a), PPAQ_28b, 1),
      PPAQ_27 = ifelse(is.na(PPAQ_20), NA, PPAQ_27)
      ) %>%
    select(-PPAQ_27a, -PPAQ_27b, -PPAQ_28a, -PPAQ_28b) %>%
    mutate_all(function(x) as.numeric(x)) %>%
    mutate(across(contains("PPAQ"), function(x) x-1)) %>%
    mutate(activity = PPAQ_20 + PPAQ_21 + PPAQ_22 + PPAQ_23 + PPAQ_24 + PPAQ_25
           + PPAQ_26 + PPAQ_27 + PPAQ_28) %>%
    select(id, week, activity)
})
ggplot(activity, aes(as.factor(week), activity)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.1)


# maternal education & delivery mode
edu <- read_excel(here("data/SMILEY_extracovariates/SMILEY_extracovariates_data.docx.xlsx")) %>%
  select(id = ID, edu = Demographics_3, csection = Delivery_mode)  %>%
  mutate(
    csection = ifelse(csection == "1", 0, ifelse(csection == "2", 1, NA)),
    edu = as.numeric(edu),
    edu = ifelse(edu == 99999, NA, ifelse(edu >= 6, 1, 0))
  )

file <- here("data/SMILEY_hair_preg32wks/SMILEY_32wks_hair_cortisol_RawClean.xlsx")
cort <- read_excel(file) %>% 
  mutate(across(everything(), function(x) ifelse(x == 99999, NA, x))) %>%
  select(-contains("mass"), -contains("comments")) %>%
  mutate(across(everything(), function(x) as.character(x))) %>%
  pivot_longer(-ID, names_to = "var", values_to = "values") %>%
  mutate(
    values = as.numeric(values),
    var = str_remove(var, "SMILEY_32wk_"),
    t = str_extract(var, "\\d"),
    t = ifelse(t == 1, "t2", ifelse(t == 2, "t1", NA)),
    var = str_extract(var, "(cortiso\\w+)\\d", group = 1),
    ID = as.double(ID)
    ) %>%
  pivot_wider(names_from = var, values_from = values) %>%
  rename(id = ID) %>%
  filter(!is.na(t))
head(cort)
file <- here("data/SMILEY_hair_preg32wks/SMILEY_32wks_hair_questionnaire_RawClean.xlsx")
# i checked descriptive and by far most have either one or two hair washings or three or four
# therefore I make two groups
cort2 <- read_excel(file) %>% 
  select(id = ID, lessthanone = "4a", oneortwo = "4b", threeorfour = "4c", morethanfour = "4d") %>%
  pivot_longer(-id, names_to = "frequency", values_to = "values") %>%
  mutate(values = ifelse(values == 99999, NA, ifelse(values == 88888, NA, values))) %>%
  # group_by(frequency) %>%
  # summarise(sum = sum(values, na.rm = TRUE))
  pivot_wider(names_from = "frequency", values_from = "values") %>%
  mutate(hw = ifelse(lessthanone, 0, ifelse(oneortwo, 0, ifelse(threeorfour, 1,  ifelse(morethanfour, 1, NA)))),
         hw = ifelse(is.na(hw), 0, hw)
         ) %>%
  select(id, hw)

cort <- full_join(cort, cort2, by = "id") %>%
  mutate(ratio = log(cortisol/cortisone), cortisol = ifelse(cortisol >= 50, NA, log(cortisol)), cortisone = log(cortisone))

mlr::summarizeColumns(cort)
plot(cort$cortisol)



dm <- full_join(activity, bmi, by = "id") %>%
  full_join(age, by = "id")%>%
  full_join(ga, by = c("id", "week"))%>%
  full_join(abx, by = c("id", "week"))%>%
  full_join(parity, by = "id") %>%
  full_join(ms, by = c("id", "week", "pre")) %>%
  full_join(ffq, by = c("id", "week", "pre")) %>%
  full_join(edu, by = c("id")) %>%
  mutate(t = ifelse(week == 18, "t1", ifelse(pre, "t2", ifelse(!pre, "t3", NA)))) %>%
  full_join(cort, by = c("id", "t"))
dim(dm)
write_csv(dm, file = here("data/dm.csv"))

# impute predictors

cortisol <- dm$cortisol
cortisone <- dm$cortisone
ratio <- dm$ratio
dmimp <- mice(dm, m = 5)
dmimp <- map(1:5, function(m) {
  out <- complete(dmimp, action = m) %>%
    group_by(id) %>%
    mutate(
      pbmi = min(pbmi), 
      csection = min(csection),
      length = min(length),
      weight = min(weight),
      age = min(age),
      parity = min(parity),
      edu = max(edu)
      ) %>%
    select(-ga) %>%
    ungroup()
  out$cortison <- cortisone
  out$cortisol <- cortisol
  out$ratio <- ratio
  out
})
test <- dmimp[[1]]

# the psas needs to be analyzed separately for the mother models 
# get psas
psas <- read_csv(here("data/psas.csv")) %>%
  group_by(ID) %>%
  summarise(psas = mean(psas, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    #id = as.character(ID), 
    psas = scale(psas)[, 1],
    stress_group = ifelse(psas <= -1, "low", ifelse(psas < 1, "medium", ifelse(psas >= 1, "high", NA)))
  ) %>%
  select(id = ID, psas)

dmimp[[1]]
psas_dmimp <- map(dmimp, function(dimp) {
  dtemp <- filter(dimp, t == "t3") %>%
    full_join(psas, by = c("id"))
  dtemp_imp <- mice(dtemp, m = 1)
  complete(dtemp_imp) %>% as.data.frame()
})

save(dm, dmimp, psas, psas_dmimp, file = here::here("data/rdata/dm.Rds"))



