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
# infants
tse_i <- tse[, colData(tse)$origin == "i"]

asv <- t(assay(tse_i, "clr"))
ait <- vegdist(asv, method = "euclidean")
ait <- as.matrix(ait)

ids <- sort(unique(colData(tse_i)$id))
colData(tse_i)
vol <- map_dfr(ids, function(id) {
  
  # each subject has maximal 4 samples:
  s1 <- glue("{id}iFALSE2")
  s2 <- glue("{id}iFALSE6")
  s3 <- glue("{id}iFALSE12")
  s4 <- glue("{id}iFALSE32")
  
  # we can only obtain vol if there is no missing sample in a pair 
  vol1 <- ifelse(s1 %in% rownames(ait) & s2 %in% rownames(ait), ait[s1, s2], NA)
  vol2 <- ifelse(s2 %in% rownames(ait) & s3 %in% rownames(ait), ait[s2, s3], NA)
  vol3 <- ifelse(s3 %in% rownames(ait) & s4 %in% rownames(ait), ait[s3, s4], NA)
  
  tibble(
    id,
    time = c("2-6", "6-12", "12-32"),
    vol = c(vol1, vol2, vol3)
  )
})
rownames(ait)
head(vol)
tail(vol)

vol_by_comp <- group_by(vol, time) %>% nest()
# visualize distributions of volatility per time point pair 
voldist <- map(vol_by_comp[[2]], function(df) {
  ggplot(df, aes(vol)) +
    geom_density()
})
voldist[[3]]
ivol1 <- vol_by_comp[[2]][[1]]
ivol2 <- vol_by_comp[[2]][[2]]
ivol3 <- vol_by_comp[[2]][[3]]
colnames(ivol1) <- c("id", "ivol1")
colnames(ivol2) <- c("id", "ivol2")
colnames(ivol3) <- c("id", "ivol3")



####################### determine model structure ############################

# as we did for alpha diversity, we must determine which structure performs better
# out of sample predictions: mlm or single models. we start with the mlm, so
# first we need data in longformat

load(here("data/rdata/di.Rds"))

dvol <- select(
  dlong,
  id,
  week,
  dmode,
  parity,
  feeding,
  gage,
  contains("cort"),
  contains("ratio"),
  contains("ms"),
  stai,
  epds,
  hw,
  pss,
  praq_fear_birth_g,
  praq_worries_handicap_g,
  psas) %>%
  filter(week != 32)

dvol$id <- as.character(dvol$id)
vol <- full_join(ivol1, ivol2, by = "id") %>%
  full_join(ivol3, by = "id") %>%
  pivot_longer(contains("vol"), names_to = "week", values_to = "ivol") %>%
  mutate(week = ifelse(week == "ivol1", "2", ifelse(
    week == "ivol2", "6", ifelse(
      week == "ivol3", "12", NA))))

dvol <- full_join(dvol, vol, by = c("id", "week")) %>%
  mutate(across(where(is.numeric), function(x) scale(x)[, 1]))



# use 10 fold crossvalidation
folds <- caret::createFolds(dvol$id, k = 10)

# first the multilevel model
i1_msqe <- map2_dbl(folds, 1:10, function(fold, k) {
  
  # splits
  train <- dvol[-fold, ]
  test <- dvol[fold, ]
  
  # fit model 
  i1 <- brm(
    family = student(),
    formula = ivol ~ dmode + gage + feeding + parity + week + (1|id),
    data = train,
    file = here(glue("data/rdata/models/ivol_i1_multilevel_{k}"))
  )
  
  pred <- predict(i1, newdata = test, re_formula = NA) %>%
    as.data.frame() %>%
    .$Estimate
  actual <- test$ivol
  # mse
  squared_diff <- (actual - pred)^2
  mse <- mean(squared_diff, na.rm = TRUE)
  mse
})
mean(i1_msqe)
dvol
# now the separate models



# splits
train <- dvol[-folds[[1]], ]
test <- dvol[folds[[1]], ]


i2_msqe <- map2_dbl(folds, 1:10, function(fold, k) {
  
  # splits
  train <- dvol[-fold, ]
  test <- dvol[fold, ]
  
  # fit models
  i1 <- brm(
    family = student(),
    formula = ivol ~ dmode + gage + feeding + parity,
    data = filter(train, week == "2"),
    file = here(glue("data/rdata/models/ivol_i1_multilevel_t1_{k}"))
  )
  
  i2 <- brm(
    family = student(),
    formula = ivol ~ dmode + gage + feeding + parity,
    data = filter(train, week == "6"),
    file = here(glue("data/rdata/models/ivol_i1_multilevel_t2_{k}"))
  )
  
  i3 <- brm(
    family = student(),
    formula = ivol ~ dmode + gage + feeding + parity,
    data = filter(train, week == "12"),
    file = here(glue("data/rdata/models/ivol_i1_multilevel_t3_{k}"))
  )
  

  
  # get mses
  pred1 <- predict(i1, newdata = filter(test, week == "6"), re_formula = NA) %>%
    as.data.frame() %>%
    .$Estimate
  actual1 <- filter(test, week == "6") %>% .$ivol
  squared_diff1 <- (actual1 - pred1)^2
  mse1 <- mean(squared_diff1, na.rm = TRUE)
  
  pred2 <- predict(i2, newdata = filter(test, week == "12"), re_formula = NA) %>%
    as.data.frame() %>%
    .$Estimate
  actual2 <- filter(test, week == "12") %>% .$ivol
  squared_diff2 <- (actual2 - pred2)^2
  mse2 <- mean(squared_diff2, na.rm = TRUE)
  
  pred3 <- predict(i3, newdata = filter(test, week == "12"), re_formula = NA) %>%
    as.data.frame() %>%
    .$Estimate
  actual3 <- filter(test, week == "12") %>% .$ivol
  squared_diff3 <- (actual3 - pred3)^2
  mse3 <- mean(squared_diff3, na.rm = TRUE)
  
  
  mean(c(squared_diff1, squared_diff2, squared_diff3), na.rm = TRUE)
  
})

mean(i1_msqe)
sd(i1_msqe)
mean(i2_msqe)
sd(i2_msqe)


# I will need to get the imputations for the covariates with missing values:
dsvol <- map(dlongiimp, function(dimps) {
  dvol <- select(
    dimps,
    id,
    week,
    dmode,
    hw,
    parity,
    feeding,
    gage,
    contains("cort"),
    contains("ratio"),
    contains("ms"),
    stai,
    epds,
    pss,
    praq_fear_birth_g,
    praq_worries_handicap_g,
    psas) %>%
    filter(week != 32)
  
  dvol$id <- as.character(dvol$id)
  
  
  dvol <- full_join(dvol, vol, by = c("id", "week")) %>%
    mutate(across(where(is.numeric), function(x) scale(x)[, 1]))
  dvol

}) 

dsvol[[1]] %>% head()

# fit a multilevel model
mivol <- brm_multiple(
    family = student(),
    formula = ivol ~ week * ms + parity + feeding + dmode + gage + (1|id),
    data = dsvol,
    file = here(glue("data/rdata/models/ivol_multilevel2"))
  )

summary(mivol)


# for exploratory purpose we also look at the individual time points
msvol <-map(c("2", "6", "12"), function(tp) {
  # fit model for individual time point
  mivol <- brm_multiple(
    family = student(),
    formula = ivol ~ ms + parity + feeding + dmode + gage,
    data = map(dsvol, ~filter(.x, week == tp)),
    file = here(glue("data/rdata/models/ivol_multilevel_{tp}"))
  )
})

map(msvol, ~summary(.x))

source(here("R/helper_functions.R"))
# now for all stress scores
stress_measures <- c("ms", "stai", "epds", "pss", "psas", 
                     "praq_fear_birth_g", "praq_worries_handicap_g",
                     "cortisol_pp", "cortisone_pp", "ratio_pp", 
                     "cortisol_g1", "cortisone_g1", "ratio_g1",
                     "cortisol_g2", "cortisone_g2", "ratio_g2",
                     "cortisol_g3", "cortisone_g3", "ratio_g3")
b <- map_dfr(stress_measures, function(stress) {
  map_dfr(c("2", "6", "12"), function(tp) {
    if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
      f <- as.formula(glue("ivol ~ {stress} + hw + parity + feeding + dmode + gage"))
    } else {
      f <- as.formula(glue("ivol ~ {stress} + parity + feeding + dmode + gage"))
    }
    
    # fit model for individual time point
    mivoltemp <- brm_multiple(
      family = student(),
      formula = f,
      data = map(dsvol, ~filter(.x, week == tp)),
      file = here(glue("data/rdata/models/ivol_{stress}_{tp}"))
    )
    
    param <- glue("b_{stress}")
    summarise_posterior(mivoltemp, parameters = param) %>%
      mutate(week = tp, stress = stress)
  })
})
print(b, n = 57)
filter(b, p <= 0.05 | p >= 0.95)


# make proper Bayesian plots 
sig_stress_measures <- c("ms", "pss","ratio_g2")
sig_weeks <- c("2", "6", "12")
models <- map2(sig_stress_measures, sig_weeks, function(stress, tp) {
  if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
    f <- as.formula(glue("ivol ~ {stress} + hw + parity + feeding + dmode + gage"))
  } else {
    f <- as.formula(glue("ivol ~ {stress} + parity + feeding + dmode + gage"))
  }
  
  # fit model for individual time point
  mivoltemp <- brm_multiple(
    family = student(),
    formula = f,
    data = map(dsvol, ~filter(.x, week == tp)),
    file = here(glue("data/rdata/models/ivol_{stress}_{tp}"))
  )
  mivoltemp
})




# first ms week 2
nd <- expand_grid(
  ms = seq(min(dvol$ms[dvol$week == "2"], na.rm = TRUE), max(dvol$ms[dvol$week == "2"], na.rm = TRUE), 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  hw = 1,
  parity = 1,
  feeding = 1,
  dmode = 1,
  gage = median(dsvol[[1]]$gage, na.rm = TRUE)
)
test <- posterior_linpred(models[[1]], newdata = nd, re.form = NA)

test_sum <- test %>%
  as.data.frame()
nd$i <- colnames(test_sum)
test_sum <- pivot_longer(test_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")


lbl <- c(
  `2` = "Week 2 - week 6 (I1 - I2)"
)



p1 <- dvol %>% filter(week == "2") %>%
  ggplot(aes(x = ms, y = ivol)) +
  geom_smooth(
    data = test_sum, 
    aes(x = ms, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#edf8e9', color = "black", alpha = 4/5, linewidth = 1/4
    ) +
  geom_point(alpha = 0.5, size = 3) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Volatility") + xlab("Maternal Stress Composite") +
  ylim(c(-4, 4)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p1
# then pss week 6 
nd <- expand_grid(
  #ms = seq(min(dvol$ms[dvol$week == "2"], na.rm = TRUE), max(dvol$ms[dvol$week == "2"], na.rm = TRUE), 0.1),
  pss = seq(min(dvol$pss[dvol$week == "6"], na.rm = TRUE), max(dvol$pss[dvol$week == "6"], na.rm = TRUE), 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  hw = 1,
  parity = 1,
  feeding = 1,
  dmode = 1,
  gage = median(dsvol[[1]]$gage, na.rm = TRUE)
)
test <- posterior_linpred(models[[2]], newdata = nd, re.form = NA)

test_sum <- test %>%
  as.data.frame()
nd$i <- colnames(test_sum)
test_sum <- pivot_longer(test_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")


lbl <- c(
  `6` = "Week 6 - week 12 (I2 - I3)"
)

p2 <- dvol %>% filter(week == "6") %>%
  ggplot(aes(x = pss, y = ivol)) +
  geom_smooth(
    data = test_sum, 
    aes(x = pss, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#bae4b3', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_point(alpha = 0.5, size = 3) +
  ylab("") + xlab("Perceived Stress Scale") +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylim(c(-4, 4)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )
p2

# lastly ratio_g2
nd <- expand_grid(
  #ms = seq(min(dvol$ms[dvol$week == "2"], na.rm = TRUE), max(dvol$ms[dvol$week == "2"], na.rm = TRUE), 0.1),
  #pss = seq(min(dvol$ms[dvol$week == "6"], na.rm = TRUE), max(dvol$ms[dvol$week == "6"], na.rm = TRUE), 0.1),
  ratio_g2 = seq(min(dvol$ratio_g2[(dvol$week == "12" & !is.na(dvol$ivol))], na.rm = TRUE), max(dvol$ratio_g2[(dvol$week == "12" & !is.na(dvol$ivol))], na.rm = TRUE), 0.1),
  id = NA,
  hw = 1,
  parity = 1,
  feeding = 1,
  dmode = 1,
  gage = median(dsvol[[1]]$gage, na.rm = TRUE)
)
test <- posterior_linpred(models[[3]], newdata = nd, re.form = NA)

test_sum <- test %>%
  as.data.frame()
nd$i <- colnames(test_sum)
test_sum <- pivot_longer(test_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")

dvol %>% filter(week == "12") %>%
  summarise(max(ratio_g2, na.rm = TRUE))

lbl <- c(
  `12` = "Week 12 - months 8 (I3 - I4)"
)

p3 <- dvol %>% filter(week == "12") %>%
  select(ratio_g2, ivol, week) %>%
  na.omit() %>%
  ggplot(aes(x = ratio_g2, y = ivol)) +
  geom_smooth(
    data = test_sum, 
    aes(x = ratio_g2, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#74c476', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_point(alpha = 0.5, size = 3) +
  ylab("") + xlab("Log-ratio of cortisol/cortisone (HCR2)") +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylim(c(-4, 4)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )
p3
save(b, p1, p2, p3, file = here("data/rdata/postnatal_vol_i.Rds"))

