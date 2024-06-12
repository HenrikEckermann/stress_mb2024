set.seed(1)
library(mia)
library(tidyverse)
library(glue)
library(vegan) 
library(here)
library(brms)
library(bayesplot)
library(posterior)
library(tidybayes)


load(file = here::here("data/rdata/tse.Rds"))
tse <- transformAssay(tse_s, method = "clr", name = "clr", pseudocount = 0.000001)
# we calculate volatility for mothers and infants separately

# mothers 
tse_m <- tse[, colData(tse)$origin == "m"]

asv <- t(assay(tse_m, "clr"))
ait <- vegdist(asv, method = "euclidean")
ait <- as.matrix(ait)
ids <- sort(unique(colData(tse_m)$id))

vol <- map_dfr(ids, function(id) {
  
  # each subject has maximal 3 samples:
  s1 <- glue("{id}mTRUE18")
  s2 <- glue("{id}mTRUE32")
  s3 <- glue("{id}mFALSE32")
  
  # we can only obtain vol if there is no missing sample in a pair 
  vol1 <- ifelse(s1 %in% rownames(ait) & s2 %in% rownames(ait), ait[s1, s2], NA)
  vol2 <- ifelse(s2 %in% rownames(ait) & s3 %in% rownames(ait), ait[s2, s3], NA)
  
  tibble(
    id,
    time = c("18-32", "pre-post"),
    vol = c(vol1, vol2)
  )
})


vol_by_comp <- group_by(vol, time) %>% nest()
# visualize distributions of volatility per time point pair 
voldist <- map(vol_by_comp[[2]], function(df) {
  ggplot(df, aes(vol)) +
    geom_density()
})

mvol1 <- vol_by_comp[[2]][[1]]
mvol2 <- vol_by_comp[[2]][[2]]
colnames(mvol1) <- c("id", "mvol1")
colnames(mvol2) <- c("id", "mvol2")




####################### determine model structure ############################

# as we did for alpha diversity, we must determine which structure performs better
# out of sample predictions: mlm or single models. we start with the mlm, so
# first we need data in longformat

load(here("data/rdata/dm.Rds"))

d <- select(dm, -week) %>%
  select(-matches("dhd15_\\d+")) %>%
  pivot_wider(names_from = t, values_from = c(
    activity, ga, abx, pre, stai, epds, pss, ms, dhd15_total,
    contains("praqr"), pesbrief, contains("cort"), hw
  )) 
dvol <- select(d, id, pbmi, parity, edu, age) 
dvol$id <- as.character(dvol$id)

dvol$activity1 <- d$activity_t1
dvol$activity2 <- d$activity_t2 
dvol$ms1 <- d$ms_t2
dvol$ms2 <- d$ms_t3
dvol$stai1 <- d$stai_t2
dvol$stai2 <- d$stai_t3
dvol$epds1 <- d$epds_t2
dvol$epds2 <- d$epds_t3
dvol$pss1 <- d$pss_t2
dvol$pss2 <- d$pss_t3
dvol$praqr_handicap1 <- d$praqr_handicap_t1
dvol$praqr_handicap2 <- d$praqr_handicap_t2
dvol$praqr_birth1 <- d$praqr_birth_t1
dvol$praqr_birth2 <- d$praqr_birth_t2
dvol$pesbrief1 <- d$pesbrief_t1
dvol$pesbrief2 <- d$pesbrief_t2
dvol$cortisol1 <- d$cortisol_t1
dvol$cortisol2 <- d$cortisol_t2
dvol$cortisone1 <- d$cortisone_t1
dvol$cortisone2 <- d$cortisone_t2
dvol$hw1 <- d$hw_t1
dvol$hw2 <- d$hw_t2
colnames(dvol)
dvol <- full_join(dvol, mvol1, by = "id") %>%
  full_join(mvol2, by = "id") %>%
  mutate(edu = as.factor(edu))
dlong <- pivot_longer(dvol, matches(".*\\d$"), names_to = c("var", "t"), names_pattern = "(\\w+)(\\d)") %>%
  pivot_wider(names_from = var, values_from = value) %>%
  mutate(
    hw = as.factor(hw),
    across(where(is.numeric), function(x) scale(x)[, 1]))
dvol <- mutate(dvol, hw1 = as.factor(hw1), hw2 = as.factor(hw2), across(where(is.numeric), function(x) scale(x)[, 1]))


# use 10 fold crossvalidation
folds <- caret::createFolds(dvol$id, k = 10)

# first the multilevel model
m1_mm_msqe <- map2_dbl(folds, 1:10, function(fold, k) {
  
  # splits
  train <- dlong[-fold, ]
  test <- dlong[fold, ]
  
  # fit model 
  m1_mm <- brm(
    family = student(),
    formula = mvol ~ activity + age + edu + parity + t + (1|id),
    data = train,
    file = here(glue("data/rdata/models/mvol_m1_multilevel_{k}"))
  )
  
  pred <- predict(m1_mm, newdata = test, re_formula = NA) %>%
    as.data.frame() %>%
    .$Estimate
  actual <- test$mvol
  # mse
  squared_diff <- (actual - pred)^2
  mse <- mean(squared_diff, na.rm = TRUE)
  mse
})
mean(m1_mm_msqe)

# now the separate models
m2_mm_msqe <- map2_dbl(folds, 1:10, function(fold, k) {
  
  # splits
  train <- dvol[-fold, ]
  test <- dvol[fold, ]
  
  # fit models
  m1 <- brm(
    family = student(),
    formula = mvol1 ~ activity1 + age + edu + parity,
    data = train,
    file = here(glue("data/rdata/models/mvol_m1_multilevel_t1_{k}"))
  )
  
  m2 <- brm(
    family = student(),
    formula = mvol2 ~ activity2 + age + edu + parity,
    data = train,
    file = here(glue("data/rdata/models/mvol_m1_multilevel_t2_{k}"))
  )
  

  
  # get mses
  pred1 <- predict(m1, newdata = test, re_formula = NA) %>%
    as.data.frame() %>%
    .$Estimate
  actual1 <- test %>% .$mvol1
  squared_diff1 <- (actual1 - pred1)^2
  mse1 <- mean(squared_diff1, na.rm = TRUE)
  
  pred2 <- predict(m2, newdata = test, re_formula = NA) %>%
    as.data.frame() %>%
    .$Estimate
  actual2 <- test %>% .$mvol2
  squared_diff2 <- (actual2 - pred2)^2
  mse2 <- mean(squared_diff2, na.rm = TRUE)
  
  
  mean(c(squared_diff1, squared_diff2), na.rm = TRUE)
  
})

# again, the multilevel model beats the single models
mean(m1_mm_msqe)
sd(m1_mm_msqe)
mean(m2_mm_msqe)
sd(m2_mm_msqe)


load(here("data/rdata/dm.Rds"))
dsvol <- map(dmimp, function(dimps) {
  dtemp <- select(dimps, -week, -matches("dhd15_\\d+")) %>%
    pivot_wider(id_cols = "id", names_from = t, values_from = c(
      activity, abx, pre, stai, epds, pss, ms, dhd15_total,
      contains("praqr"), pesbrief, contains("cort"), contains("hw")
    )) 
  
  # i am not sure why bmi disappears but due to being time efficient i just add it manually with a full join
  dtemp$pbmi <- filter(dimps, week == 18) %>% .$pbmi
  dtemp$parity <- filter(dimps, week == 18) %>% .$parity
  dtemp$edu <- filter(dimps, week == 18) %>% .$edu
  dtemp$age <- filter(dimps, week == 18) %>% .$age
  dtemp$hw <- filter(dimps, week == 18) %>% .$hw
  

  dvol <- select(dtemp, id, pbmi, parity, edu, age, hw) 
  dvol$id <- as.character(dvol$id)
  
  dvol$activity1 <- dtemp$activity_t1
  dvol$activity2 <- dtemp$activity_t2 
  dvol$ms1 <- dtemp$ms_t2
  dvol$ms2 <- dtemp$ms_t3
  dvol$stai1 <- dtemp$stai_t2
  dvol$stai2 <- dtemp$stai_t3
  dvol$epds1 <- dtemp$epds_t2
  dvol$epds2 <- dtemp$epds_t3
  dvol$pss1 <- dtemp$pss_t2
  dvol$pss2 <- dtemp$pss_t3
  dvol$praqr_handicap1 <- dtemp$praqr_handicap_t1 
  dvol$praqr_handicap2 <- dtemp$praqr_handicap_t1 
  dvol$praqr_birth1 <- dtemp$praqr_birth_t1 
  dvol$praqr_birth2 <- dtemp$praqr_birth_t1
  dvol$pesbrief1 <- dtemp$pesbrief_t1
  dvol$pesbrief2 <- dtemp$pesbrief_t1 
  dvol$cortisol1 <- dtemp$cortisol_t1
  dvol$cortisol2 <- dtemp$cortisol_t2
  dvol$cortisone1 <- dtemp$cortisone_t1
  dvol$cortisone2 <- dtemp$cortisone_t2

  
  dvol <- full_join(dvol, mvol1, by = "id") %>%
    full_join(mvol2, by = "id") %>%
    mutate(edu = as.factor(edu))
  dlong <- pivot_longer(dvol, matches(".*\\d$"), names_to = c("var", "t"), names_pattern = "(\\w+)(\\d)") %>%
    pivot_wider(names_from = var, values_from = value) %>%
    mutate(
      hw = as.factor(hw),
      across(where(is.numeric), function(x) scale(x)[, 1]))
  filter(dlong, !is.na(mvol))
}) 

psas_dsvol <- map(psas_dmimp, function(dimps) {
  dvol <- dimps %>% mutate(id = as.character(id)) %>% select(-t)
  dvol <- full_join(dvol, mvol1, by = "id") %>%
    full_join(mvol2, by = "id") %>%
    mutate(edu = as.factor(edu))
  dlong <- pivot_longer(dvol, contains("mvol"), names_to = "t", names_pattern = "\\w+(\\d)", values_to = "mvol") %>%
    #pivot_wider(names_from = var, values_from = value) %>%
    mutate(across(where(is.numeric), function(x) scale(x)[, 1]))
  filter(dlong, !is.na(mvol))
}) 

# fit final multilevel model
mmvol <- brm_multiple(
    family = student(),
    formula = mvol ~ activity + age + edu + parity + t * ms + (1|id),
    data = dsvol,
    file = here(glue("data/rdata/models/mvol_multilevel"))
  )

source(here("R/helper_functions.R"))
b <- summarise_posterior(mmvol, parameters = c("b_ms", "b_t2"))
save(b, file = here("data/rdata/bmvol.Rds"))

summary(mmvol)

# check individual questionnaires
qs <- c("ms", "stai", "epds", "pss", "praqr_birth", "praqr_handicap", "pesbrief", "cortisol", "cortisone", "ratio")
models <- map_dfr(qs, function(stress) {
  if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
    f <- as.formula(glue("mvol ~ activity + age + edu + parity + t * {stress} + hw + (1|id)"))
  } else {
    f <- as.formula(glue("mvol ~ activity + age + edu + parity + t * {stress} + (1|id)"))
  }
  
  # fit final multilevel model
  mmvol <- brm_multiple(
    family = student(),
    formula = f,
    data = dsvol,
    file = here(glue("data/rdata/models/mvol_multilevel_{stress}"))
  )
  
  summarise_posterior(mmvol, parameters = c(glue("b_{stress}"), glue("b_t2:{stress}"))) %>%
    mutate(stress = stress)
})
models


# for exploratory purpose we still look at the individual time points but given
# the interaction, that should not differ.
msvol <- map_dfr(qs, function(stress) {
  map_dfr(c(1, 2), function(tp) {
    if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
      f <- as.formula(glue("mvol ~ activity + age + edu + parity + {stress} + hw"))
    } else {
      f <- as.formula(glue("mvol ~ activity + age + edu + parity + {stress}"))
    }
    # fit model for individual time point
    mmvol <- brm_multiple(
      family = student(),
      formula = f,
      data = map(dsvol, ~filter(.x, t == tp)),
      file = here(glue("data/rdata/models/mvol_{stress}_{tp}"))
    )
    
    summarise_posterior(mmvol, parameters = glue("b_{stress}")) %>%
      mutate(stress = stress, t = tp) %>%
      filter(parameter == glue("b_{stress}"))
  })
})

msvol





# the results are in line: maternal stress shows no clear/strong relationship
# with volatility. 
# check if the PSAS is related to vol between 32 weeks and 8 months pp 

mmvol <- brm_multiple(
  family = student(),
  formula = mvol ~ activity + age + edu + parity + t * psas + (1|id),
  data = psas_dsvol,
  file = here(glue("data/rdata/models/mvol_multilevel_psas"))
)

summary(mmvol)


