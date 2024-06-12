set.seed(1)
library(mia)
library(tidyverse)
library(vegan)
library(here)
library(brms)
library(glue)
library(HDInterval)
library(dagitty)
source(here("R/helper_functions.R"))

load(file = here::here("data/rdata/tse.Rds"))
load(here("data/rdata/di.Rds"))
tse <- transformAssay(tse_s, method = "clr", name = "clr", pseudocount = 0.000001)

# prepare data for analisys of infants
tse_i <- tse[, colData(tse)$origin == "i"]
colData(tse_i) %>% colnames()
dlong$week <- as.numeric(dlong$week)
dlong$id <- as.character(dlong$id)
# the 2 year samples were not ready on time for this project, so I must filter them
tse_i <- tse_i[, tse_i$week != 104]
# now leftjoin meta data with tse
colData(tse_i) <- colData(tse_i) %>%
                    as.data.frame() %>%
                    rownames_to_column("rn") %>%
                    mutate(
                      #week = factor(week, levels = c("2", "6", "12", "32")),
                      t = ifelse(week == "2", "t1", ifelse(week == "6", "t2", ifelse(week == "12", "t3",ifelse(week == "32", "t4",NA)))),
                      t = factor(t, levels = c("t1", "t2", "t3", "t4"))
                      ) %>%
                    left_join(dlong, by = c("id", "week")) %>%
                    column_to_rownames("rn") %>%
                    DataFrame()


# alpha diversity
tse_i <- estimateDiversity(
    tse_i, 
    assay.type = "counts", 
    index = "shannon")

# visualize 
iad_plots <- map(c("shannon", "faith"), function(alpha) {
  ylabtitle <- ifelse(alpha == "shannon", glue("{str_to_title(alpha)} index"),
                      str_to_title(alpha))
  colData(tse_i) %>%
    as.data.frame() %>%
    ggplot(aes_string("t", alpha, fill = "t")) +
    geom_boxplot(outlier.alpha = 1) +
    geom_jitter(width = 0.1, alpha = 0.6, size = 3) +
    ylab(ylabtitle) + xlab("") +
    scale_fill_manual(values = c('#edf8e9','#bae4b3','#74c476','#238b45')) +
    theme_bw(base_size = 15) +
    theme(legend.position = "none")
})
iad_plots[[2]]

save(iad_plots, file = here("data/rdata/iad_plots.Rds"))

# store relevant data for statistical analysis in df and prepare
d <- colData(tse_i) %>% as.data.frame() %>%
    mutate(across(
      all_of(c("faith", "bmi_pp", "gage", "psas", "birthweight", "age", 
               "shannon", "stai", "epds", "pss")),
            function(x) scale(x)[, 1])) %>%
    mutate(
      ms = scale(epds + stai + pss)[, 1], 
      across(all_of(c(contains("cort"), contains("praq"))), function(x) scale(x)[, 1]),
      stress_group = ifelse(ms <= -1, "low", ifelse(ms < 1, "medium", ifelse(ms >= 1, "high", NA))))

cnames <- colnames(dlong)[!colnames(dlong) %in% c("id", "t")]
# also long format version including imputed covariates
dlongiimp <- map(dlongiimp, function(dimp) {
    dtemp <- dimp
    dtemp$id <- as.character(dtemp$id)
    dtemp <- mutate(
      dtemp, 
      t = ifelse(week == "2", "t1", ifelse(week == "6", "t2", ifelse(week == "12", "t3",ifelse(week == "32", "t4",NA)))),
      t = factor(t, levels = c("t1", "t2", "t3", "t4")))
    # store relevant data in df
    dout <- colData(tse_i) %>% as.data.frame() %>%
        select(-all_of(c(cnames, "week", "pre"))) %>%
        left_join(dtemp, by = c("id", "t")) %>%
        mutate(across(all_of(c("faith", "bmi_pp", "gage", "psas", "birthweight", "age", 
                               "shannon", "stai", "epds", "pss")), 
            function(x) scale(x)[, 1])) %>%
        mutate(
          ms = scale(epds + stai + pss)[, 1], 
          across(all_of(c(contains("cort"), contains("praq"), contains("ratio"))), function(x) scale(x)[, 1]),
          stress_group = ifelse(ms <= -1, "low", ifelse(ms < 1, "medium", ifelse(ms >= 1, "high", NA))))
    dout
})



# first we need to define the optimal structure suitable to answer
# our research question. which ones
# are the covariates that we MUST include?
colnames(d)
g <- dagitty('dag {
  bb="0,0,1,1"
  b_abx [pos="0.518,0.397"]
  dmode [pos="0.211,0.245"]
  edu [pos="0.125,0.177"]
  feeding [pos="0.444,0.134"]
  gage [pos="0.252,0.114"]
  ms [pos="0.332,0.292"]
  parity [pos="0.111,0.385"]
  sex [pos="0.288,0.405"]
  shannon [pos="0.535,0.252"]
  b_abx -> shannon
  dmode -> ms
  dmode -> shannon
  edu -> dmode
  edu -> feeding
  edu -> shannon
  feeding -> shannon
  gage -> dmode
  gage -> feeding
  gage -> ms
  gage -> shannon
  ms -> feeding
  ms -> shannon
  parity -> edu
  parity -> feeding
  parity -> gage
  parity -> ms
  parity -> shannon
  sex -> dmode
  sex -> shannon
}
')

# we must include dmode, gage and parity minimally 
adjustmentSets(g, exposure = "ms", outcome = "shannon", type = "minimal")
# which ones are optional
adjustmentSets(g, exposure = "ms", outcome = "shannon", type = "all")
# if this DAG is correct, we can include all these:
# { b_abx, dmode, edu, feeding, gage, parity, sex } based on model fit . But first I check CI
# to potentially modify the DAG:
impliedConditionalIndependencies(g)


# abx _||_ dmod
summary(glm(b_abx ~dmode, d, family = binomial(link = "logit")))
# abx _||_ edu
summary(glm(b_abx ~edu, d, family = binomial(link = "logit")))
# abx _||_ fdng
summary(glm(b_abx ~feeding, d, family = binomial(link = "logit")))
# abx _||_ gage
summary(glm(b_abx ~gage, d, family = binomial(link = "logit")))
# abx _||_ ms
summary(glm(b_abx ~ms, d, family = binomial(link = "logit")))
# abx _||_ prty
summary(glm(b_abx ~parity, d, family = binomial(link = "logit")))
# abx _||_ sex
summary(glm(b_abx ~sex, d, family = binomial(link = "logit")))
# dmod _||_ fdng | edu
summary(glm(dmode ~feeding + edu, d, family = binomial(link = "logit")))
# dmod _||_ prty | edu
summary(glm(dmode ~parity + edu, d, family = binomial(link = "logit")))
# edu _||_ gage
summary(glm(edu ~gage, d, family = binomial(link = "logit")))
# edu _||_ ms | dmod, gage, prty
summary(glm(edu ~ms + dmode + gage + parity, d, family = binomial(link = "logit")))
# edu _||_ sex
summary(glm(edu ~sex, d, family = binomial(link = "logit")))
# fdng _||_ ms | dmod, gage, prty
summary(glm(feeding ~ms + dmode + gage + parity, d, family = binomial(link = "logit")))
# fdng _||_ ms | edu, gage, prty
summary(glm(feeding ~ms  + gage + parity + edu, d, family = binomial(link = "logit")))
# fdng _||_ sex
summary(glm(feeding ~sex, d, family = binomial(link = "logit")))
# gage _||_ sex
summary(lm(gage ~ sex, d))
# ms _||_ sex | dmod, edu, gage
summary(lm(ms ~ sex + dmode + gage, d))
# ms _||_ sex | dmod, gage, prty
summary(lm(ms ~ sex + dmode + gage + parity, d))
# prty _||_ sex
summary(glm(parity ~sex, d, family = binomial(link = "logit")))
# the DAG has been adjusted based on testing the conditional independences. 

# is it better to treat the microbiota at the 
# different time points as different system (for each time point fit a model) or should we fit a 
# multi level model. I will try both and depending on which model makes better predictions (msqe)
# we will use that for inference. I am aware that this CV check is not optimal as the predicitons
# cannot benefit use the random intercept for prediction. However, the ML model would still be
# expected to outperform because it should better capture the actual signal in the predictors.

# use 10 fold crossvalidation
folds <- caret::createFolds(d$sid, k = 10)
# first the multilevel model
i1_mm_msqe <- map2_dbl(folds, 1:10, function(fold, k) {

    # splits
    train <- d[-fold, ]
    test <- d[fold, ]

    # fit model 
    i1_mm <- brm(
        family = student(),
        formula = shannon ~ dmode + gage  + parity + t + (1|id),
        data = train,
        file = here(glue("data/rdata/models/i1_multilevel_{k}"))
    )

    pred <- predict(i1_mm, newdata = test, re_formula = NA) %>%
      as.data.frame() %>%
      .$Estimate
    actual <- test$shannon
    # mse
    squared_diff <- (actual - pred)^2
    mse <- mean(squared_diff, na.rm = TRUE)
    mse
})
mean(i1_mm_msqe)

# now the separate models
i2_mm_msqe <- map2_dbl(folds, 1:10, function(fold, k) {

    # splits
    train <- d[-fold, ]
    test <- d[fold, ]

    # fit models
    i1 <- brm(
        family = student(),
        formula = shannon ~ dmode + gage + parity,
        data = filter(train, t == "t1"),
        file = here(glue("data/rdata/models/i1_multilevel_t1_{k}"))
    )

    i2 <- brm(
        family = student(),
        formula = shannon ~ dmode + gage + parity,
        data = filter(train, t == "t2"),
        file = here(glue("data/rdata/models/i2_multilevel_t2_{k}"))
    )

    i3 <- brm(
        family = student(),
        formula = shannon ~ dmode + gage + parity,
        data = filter(train, t == "t3"),
        file = here(glue("data/rdata/models/i3_multilevel_t3_{k}"))
    )
    
    i4 <- brm(
      family = student(),
      formula = shannon ~ dmode + gage + parity,
      data = filter(train, t == "t4"),
      file = here(glue("data/rdata/models/i4_multilevel_t4_{k}"))
    )

    # get mses
    pred1 <- predict(i1, newdata = filter(test, t == "t1"), re_formula = NA) %>%
        as.data.frame() %>%
        .$Estimate
    actual1 <- filter(test, t == "t1") %>% .$shannon
    squared_diff1 <- (actual1 - pred1)^2
    mse1 <- mean(squared_diff1, na.rm = TRUE)

    pred2 <- predict(i2, newdata = filter(test, t == "t2"), re_formula = NA) %>%
        as.data.frame() %>%
        .$Estimate
    actual2 <- filter(test, t == "t2") %>% .$shannon
    squared_diff2 <- (actual2 - pred2)^2
    mse2 <- mean(squared_diff2, na.rm = TRUE)

    pred3 <- predict(i3, newdata = filter(test, t == "t3"), re_formula = NA) %>%
        as.data.frame() %>%
        .$Estimate
    actual3 <- filter(test, t == "t3") %>% .$shannon
    squared_diff3 <- (actual3 - pred3)^2
    mse3 <- mean(squared_diff3, na.rm = TRUE)
    
    pred4 <- predict(i4, newdata = filter(test, t == "t4"), re_formula = NA) %>%
      as.data.frame() %>%
      .$Estimate
    actual4 <- filter(test, t == "t4") %>% .$shannon
    squared_diff4 <- (actual4 - pred4)^2
    mse4 <- mean(squared_diff4, na.rm = TRUE)
    
    

    mean(c(squared_diff1, squared_diff2, squared_diff3, squared_diff4), na.rm = TRUE)

})


# they perform equally well
mean(i1_mm_msqe)
sd(i1_mm_msqe)
mean(i2_mm_msqe)
sd(i2_mm_msqe)




# last step before we continue is: 3) which of the optional covariates improve model fit 
# and should therefore be included. These are obligated: { dmode, gage, parity }
# optional: { b_abx, edu, feeding, sex } 

loo_comp <- map_dfr(1:5, function(m) {
  map_dfr(unique(d$t), function(time) {
    m_base <- brm(
      formula = shannon ~ parity + dmode + gage,
      data = filter(dlongiimp[[m]], t == time),
      file = glue(here("data/rdata/models/i_base_{m}_{time}"))
    )
    loo_mbase <- add_criterion(
      m_base,
      "loo",
      file = glue(here("data/rdata/models/modelsloo_ibase_{m}_{time}")),
      moment_match = FALSE
    )
    
    if (time == "t4") {
      f <- bf(shannon ~ parity + dmode + gage)
    } else {
      f <- bf(shannon ~ parity + dmode + gage  + b_abx)
    }
    m_1 <- brm(
      formula = f,
      data = filter(dlongiimp[[m]], t == time),
      file = glue(here("data/rdata/models/i_1_{m}_{time}"))
    )
    loo_1 <- add_criterion(
      m_1,
      "loo",
      file = glue(here("data/rdata/models/modelsloo_i1_{m}_{time}")),
      moment_match = FALSE
    )
    
    
    m_2 <- brm(
      formula = shannon ~ parity + dmode + gage + edu,
      data = filter(dlongiimp[[m]], t == time),
      file = glue(here("data/rdata/models/i_2_{m}_{time}"))
    )
    loo_2 <- add_criterion(
      m_2,
      "loo",
      file = glue(here("data/rdata/models/modelsloo_i2_{m}_{time}")),
      moment_match = FALSE
    )
    
    m_3 <- brm(
      formula = shannon ~ parity + dmode + gage + feeding,
      data = filter(dlongiimp[[m]], t == time),
      file = glue(here("data/rdata/models/i_3_{m}_{time}"))
    )
    loo_3 <- add_criterion(
      m_3,
      "loo",
      file = glue(here("data/rdata/models/modelsloo_i3_{m}_{time}")),
      moment_match = FALSE
    )
    
    m_4 <- brm(
      formula = shannon ~ parity + dmode + gage + sex,
      data = filter(dlongiimp[[m]], t == time),
      file = glue(here("data/rdata/models/i_4_{m}_{time}"))
    )
    loo_4 <- add_criterion(
      m_4,
      "loo",
      file = glue(here("data/rdata/models/modelsloo_i4_{m}_{time}")),
      moment_match = FALSE
    )
    
    
    
    tbl <- loo_compare(loo_mbase, loo_1, loo_2, loo_3, loo_4)
    tbl %>% as.data.frame() %>%
      mutate(order = 1:dim(as.data.frame(tbl))) %>%
      rownames_to_column("model") %>%
      mutate(t = time)
  })
})

count(loo_comp, t, order, model) %>%
  arrange(t, order, desc(n))


dlongiimp[[1]]
map(dlongiimp, function(dtemp) {
  filter(dtemp, t == "t4")
})



# lets see if postnatal maternal stress has an influence on alpha diversity over time?
stress_measures <- c("ms", "stai", "epds", "pss", "psas", "cortisol_pp", "cortisone_pp", "ratio_pp")
mad_ms <- map(stress_measures, function(stress) {
    map(c("shannon", "faith"), function(alpha) {
      map(unique(dlongiimp[[1]]$t), function(time) {
        
        if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
          form <- as.formula(glue("{alpha} ~ {stress} + parity + dmode + gage + hw"))
        } else {
          form <- as.formula(glue("{alpha} ~ {stress} + parity + dmode + gage"))
        }
        
        brm_multiple(
          family = student(),
          formula = form,
          data = map(dlongiimp, ~filter(.x, t == time)),
          file = here(glue("data/rdata/models/iad_{stress}_{alpha}_{time}"))
        )
      })
    })
})

# from the results printed below that if there is a linear relationship, it is likely a positive
# association and it is strongest at T3. However, we cannot state any of that with confidence. 
# furthermore, the effects are similar between stress measurements

b <- map2_dfr(mad_ms, stress_measures, function(l, stress) {
    map2_dfr(l, c("shannon", "faith"), function(models, ad) {
      map2_dfr(models, unique(d$t), function(x, time) {
        params <- colnames(x$rhats)[str_detect(colnames(x$rhats), "b_")]
        params <- str_replace(params, "\\.", ":")
        summarise_posterior(x, parameters = params) %>%
          mutate(alpha = ad, stress = stress, t = time) %>%
          filter(str_detect(parameter, stress))
      })
  })
})
print(b, n = 64)
# b <- b[map_lgl(b$parameter, function(param) mean(str_detect(param, stress_measures)) > 0 ), ] %>%
#   arrange(t, alpha, parameter)
save(b, file = here("data/rdata/b_stress_ad_i.Rds"))




# we explore the same now using low and high stress groups

iad_ms_group <- map(c("shannon", "faith"), function(alpha) {
  map(unique(d$t), function(time) {
    form <- as.formula(glue("{alpha} ~ stress_group + parity + dmode + gage"))
    brm_multiple(
      family = student(),
      formula = form,
      data = map(dlongiimp, ~filter(.x, t == time)),
      file = here(glue("data/rdata/models/iad_ms_{alpha}_{time}_group"))
    )
  })
})


b <- map_dfr(c("shannon", "faith"), function(alpha) {
  map(seq_along(unique(d$t)), function(t_d) {
      nd <- expand.grid(
          dmode = c("1", "2"),
          parity = c(0, 1), 
          gage = median(d$gage, na.rm = TRUE),
          stress_group = c("low", "medium", "high"),
          id = NA
      )
      temp_join <- select(nd, parity, dmode, gage, stress_group) %>% mutate(V = glue("V{1:dim(nd)[1]}"))
      epred <- posterior_epred(
          iad_ms_group[[ifelse(alpha == "shannon", 1, 2)]][[t_d]],
          re.formula = NA,
          newdata = nd) %>% as.data.frame() %>%
          pivot_longer(everything(), names_to = "V", values_to = "mu") %>%
          left_join(temp_join, by = "V")
      # now we can calculate any contrast we like:
      # main effects stress group
      l <- filter(epred, stress_group == "low") %>% .$mu
      m <- filter(epred, stress_group == "medium") %>% .$mu
      h <- filter(epred, stress_group == "high") %>% .$mu
  
      # calculate contrasts
      lm <- l - m
      lh <- l - h 
      mh <- m - h 
  
      map2_dfr(
        list(lm, lh, mh), 
          c("low - medium", "low - high", "medium - high"), 
          function(contrast, name) {
            m <- median(contrast)
            lower <- hdi(contrast)[1]
            upper <- hdi(contrast)[2]
            p <- mean(contrast > 0)
            tibble(time = unique(d$t)[t_d], alpha = alpha, contrast = name, m, lower, upper, p) %>%
                mutate(across(where(is.numeric), round, 3))
      })
  })
})
print(arrange(b, time, alpha, contrast), n = 24)
save(b, file = here("data/rdata/b_stress_group_i.Rds"))



# make plots

lbl <- c(
  `t4` = "MS 8 months postpartum"
)
p1 <- filter(d, t == "t4", !is.na(ms)) %>%
  mutate(stress_group = factor(stress_group, levels = c("low", "medium", "high"))) %>% 
  ggplot(aes(stress_group, shannon)) +
  geom_boxplot(outlier.alpha = 0, fill = '#238b45') +
  facet_wrap(~t, labeller = as_labeller(lbl), strip.position = "top") +
  geom_jitter(width = 0.2, size = 3) +
  #facet_wrap(~t, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index 8 months postpartum") + xlab("Maternal Stress Composite") +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )


p1

# proper bayesian plots to show the regression coefs
f <- as.formula(glue("shannon ~ epds + parity + dmode + gage"))
m <- brm_multiple(
  family = student(),
  formula = form,
  data = map(dlongiimp, ~filter(.x, t == time)),
  file = here(glue("data/rdata/models/iad_epds_shannon_t4"))
)

nd <- expand_grid(
  epds = seq(min(d$epds[d$week == "32"], na.rm = TRUE), max(d$epds[d$week == "32"], na.rm = TRUE), 0.1),
  id = NA,
  gage = median(dlongiimp[[1]]$gage, na.rm = TRUE),
  dmode = 1,
  parity = 1
)

post <- posterior_linpred(m, newdata = nd, re.form = NA)

post_sum <- post %>%
  as.data.frame()
nd$i <- colnames(post_sum)
post_sum <- pivot_longer(post_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")

lbl <- c(
  `t4` = "EPDS 8 months postpartum"
)
p2 <- d %>% filter(t == "t4") %>%
  ggplot(aes(x = epds, y = shannon)) +
  geom_smooth(
    data = post_sum, 
    aes(x = epds, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#238b45', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~t, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index 8 months postpartum") + xlab("EPDS") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p2


f <- as.formula(glue("shannon ~ pss + parity + dmode + gage"))
m <- brm_multiple(
  family = student(),
  formula = form,
  data = map(dlongiimp, ~filter(.x, t == time)),
  file = here(glue("data/rdata/models/iad_pss_shannon_t4"))
)

nd <- expand_grid(
  pss = seq(min(d$pss[d$week == "32"], na.rm = TRUE), max(d$pss[d$week == "32"], na.rm = TRUE), 0.1),
  id = NA,
  gage = median(dlongiimp[[1]]$gage, na.rm = TRUE),
  dmode = 1,
  parity = 1
)

post <- posterior_linpred(m, newdata = nd, re.form = NA)

post_sum <- post %>%
  as.data.frame()
nd$i <- colnames(post_sum)
post_sum <- pivot_longer(post_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")

lbl <- c(
  `t4` = "PSS-10 8 months postpartum"
)
p3 <- d %>% filter(t == "t4") %>%
  ggplot(aes(x = pss, y = shannon)) +
  geom_smooth(
    data = post_sum, 
    aes(x = pss, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#238b45', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~t, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index 8 months postpartum") + xlab("PSS-10") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p3

save(p1, p2, p3, file = here("data/rdata/ad_postnatal_stress_p.Rds"))

# we found that at t3, the low group had lower shannon than the other groups
# same effect direction but not significant for Faith. This probably responsible
# for the slight positive trend we found in the continuous models
# now that we checked postnatal stress and alpha diversity, we will check associations
# with prenatal stress

# get the data also for mothers 
load(here("data/rdata/dm.Rds"))

stress_measures_i <- c("ms", "stai", "epds", "pss")
stress_measures_m <- c("ms", "stai", "epds", "pss", "praqr_birth", "praqr_handicap")
stress_measures <- c("ms", "stai", "epds", "pss", "praqr_birth", "praqr_handicap",
                     "cortisol_g1", "cortisol_g2", "cortisol_g3",
                     "cortisone_g1", "cortisone_g2", "cortisone_g3",
                     "ratio_g1", "ratio_g2", "ratio_g3",
                     "cortisol_pp", "cortisone_pp", "ratio_pp"
                     )
indices <- c("shannon", "faith")
mother_times <- c("t1", "t2")
infant_times <- c(2, 6, 12, 32)


b <- map_dfr(stress_measures, function(stress) {
  map_dfr(indices, function(alpha) {
    map_dfr(mother_times, function(mt) {
      map_dfr(infant_times, function(it) {
        
        if ((str_detect(stress, "cort") | str_detect(stress, "ratio") & mt == "t2")) {
          return(NULL)
        }
        
        # select and filter data for this RQ
        md <- map(dmimp, function(mdtemp) {
          filter(mdtemp, t == mt) %>%
            select(id, all_of(stress_measures_m)) %>%
            mutate(id = as.character(id))
        })
        id <- map(dlongiimp, function(idtemp) {
          filter(idtemp, week == it) %>%
            select(id, parity, feeding, gage, dmode, all_of(alpha),
                   -all_of(stress_measures_i), contains("cort"), contains("ratio"), hw)
        })
        
        # combine the data for modeling
        dtemp <- map2(md, id, function(mdtemp, idtemp) {
          full_join(mdtemp, idtemp, by = "id")
        })
        
        # model
        if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
          f <- as.formula(glue("{alpha} ~ {stress} + parity + dmode + gage + hw"))
        } else {
          f <- as.formula(glue("{alpha} ~ {stress} + parity + dmode + gage"))
        }
        
        m <- brm_multiple(
            family = student(),
            formula = f,
            data = dtemp,
            file = here(glue("data/rdata/models/iad_{stress}_{alpha}_{mt}_{it}"))
          )
        param <- glue("b_{stress}")
        summarise_posterior(m, parameters = param) %>%
          mutate(
            infant_time = it,
            mother_time = mt,
            stress = stress,
            alpha = alpha
          )
      })
    })
  })
})

b

# b1 <- filter(b, (str_detect(parameter, "cort") | str_detect(parameter, "ratio"))) %>%
#   filter(mother_time == "t1") %>%
#   mutate(mother_time = "")
# b2 <- filter(b, !(str_detect(parameter, "cort") | str_detect(parameter, "ratio")))
# b <- bind_rows(b1, b2) 


save(b, file = here("data/rdata/b_stress_prenatal_i.Rds"))

b %>%
  arrange(infant_time)  %>%
  filter(p <= 0.05 | p >= 0.95) %>%
  print(n = 80) 

# obtain models for plotting
stress_measures <- c("ms", "stai", "epds", "pss", "praqr_birth", "praqr_handicap")
stress_measures_m <- c("ms", "stai", "epds", "pss", "praqr_birth", "praqr_handicap")
indices <- c("shannon", "faith")
mother_times <- c("t1", "t2")
infant_times <- c(2, 6, 12, 32)
mlist <- map(stress_measures, function(stress) {
  map(indices, function(alpha) {
    map(mother_times, function(mt) {
      map(infant_times, function(it) {
        
        if ((str_detect(stress, "cort") | str_detect(stress, "ratio") & mt == "t2")) {
          return(NULL)
        }
        
        # select and filter data for this RQ
        md <- map(dmimp, function(mdtemp) {
          filter(mdtemp, t == mt) %>%
            select(id, all_of(stress_measures_m)) %>%
            mutate(id = as.character(id))
        })
        id <- map(dlongiimp, function(idtemp) {
          filter(idtemp, week == it) %>%
            select(id, parity, feeding, gage, dmode, all_of(alpha),
                   -all_of(stress_measures_i), contains("cort"))
        })
        
        # combine the data for modeling
        dtemp <- map2(md, id, function(mdtemp, idtemp) {
          full_join(mdtemp, idtemp, by = "id")
        })
        
        # model
        if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
          f <- as.formula(glue("{alpha} ~ {stress} + parity + dmode + gage + hw"))
        } else {
          f <- as.formula(glue("{alpha} ~ {stress} + parity + dmode + gage"))
        }
        
        m <- brm_multiple(
          family = student(),
          formula = f,
          data = dtemp,
          file = here(glue("data/rdata/models/iad_{stress}_{alpha}_{mt}_{it}"))
        )
        m
      })
    })
  })
})



# make proper bayesian plot for these relations

# stai, shannon t1, 2

nd <- expand_grid(
  stai = seq(min(dm$stai[dm$week == "18"], na.rm = TRUE), max(dm$stai[dm$week == "18"], na.rm = TRUE), 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  #hw = 1,
  gage = median(dlongiimp[[1]]$gage, na.rm = TRUE),
  dmode = 1,
  parity = 1,
  
)

post <- posterior_linpred(mlist[[2]][[1]][[1]][[1]], newdata = nd, re.form = NA)

post_sum <- post %>%
  as.data.frame()
nd$i <- colnames(post_sum)
post_sum <- pivot_longer(post_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")


lbl <- c(
  `2` = "STAI prenatal week 18"
)
p1 <- d %>% filter(week == "2") %>%
  ggplot(aes(x = stai, y = shannon)) +
  geom_smooth(
    data = post_sum, 
    aes(x = stai, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#edf8e9', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index 2 weeks postpartum (I1)") + xlab("STAI") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p1


# handicap, shannon, t2, 2

nd <- expand_grid(
  praqr_handicap = seq(min(dm$praqr_handicap[dm$week == "32"], na.rm = TRUE), max(dm$praqr_handicap[dm$week == "32"], na.rm = TRUE), 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  #hw = 1,
  gage = median(dlongiimp[[1]]$gage, na.rm = TRUE),
  dmode = 1,
  parity = 1,
  
)

post <- posterior_linpred(mlist[[6]][[1]][[2]][[1]], newdata = nd, re.form = NA)

post_sum <- post %>%
  as.data.frame()
nd$i <- colnames(post_sum)
post_sum <- pivot_longer(post_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")


lbl <- c(
  `2` = "PRAQR2-H prenatal week 32"
)
p2 <- d %>% filter(week == "2") %>%
  ggplot(aes(x = praq_worries_handicap_g, y = shannon)) +
  geom_smooth(
    data = post_sum, 
    aes(x = praqr_handicap, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#edf8e9', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index 2 weeks postpartum (I1)") + xlab("PRAQR2-H (fear of a handicapped child)") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p2


# ms, faith, t2, 6

nd <- expand_grid(
  ms = seq(min(scale(dm$ms[dm$week == "32"])[, 1], na.rm = TRUE), max(scale(dm$ms[dm$week == "32"])[, 1], na.rm = TRUE), 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  #hw = 1,
  gage = median(dlongiimp[[1]]$gage, na.rm = TRUE),
  dmode = 1,
  parity = 1,
  
)

post <- posterior_linpred(mlist[[1]][[2]][[2]][[2]], newdata = nd, re.form = NA)

post_sum <- post %>%
  as.data.frame()
nd$i <- colnames(post_sum)
post_sum <- pivot_longer(post_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")


lbl <- c(
  `6` = "MS prenatal week 32"
)

p3 <- d %>% filter(week == "6") %>%
  ggplot(aes(x = ms, y = faith)) +
  geom_smooth(
    data = post_sum, 
    aes(x = ms, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#bae4b3', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Faith 6 weeks postpartum (I2)") + xlab("MS") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p3



# epds, faith, t2, 6

nd <- expand_grid(
  epds = seq(min(dm$epds[dm$week == "32"], na.rm = TRUE), 4, 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  #hw = 1,
  gage = median(dlongiimp[[1]]$gage, na.rm = TRUE),
  dmode = 1,
  parity = 1,
  
)

post <- posterior_linpred(mlist[[3]][[2]][[2]][[2]], newdata = nd, re.form = NA)

post_sum <- post %>%
  as.data.frame()
nd$i <- colnames(post_sum)
post_sum <- pivot_longer(post_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")


lbl <- c(
  `6` = "EPDS prenatal week 32"
)

p4 <- d %>% filter(week == "6") %>%
  ggplot(aes(x = epds, y = faith)) +
  geom_smooth(
    data = post_sum, 
    aes(x = epds, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#bae4b3', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Faith 6 weeks postpartum (I2)") + xlab("EPDS") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p4


# birth, shannon, t1/t2, 12
nd <- expand_grid(
  praqr_birth = seq(min(dm$praqr_birth[dm$week == "18"], na.rm = TRUE), max(dm$praqr_birth[dm$week == "18"], na.rm = TRUE), 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  #hw = 1,
  gage = median(dlongiimp[[1]]$gage, na.rm = TRUE),
  dmode = 1,
  parity = 1,
  
)

post <- posterior_linpred(mlist[[5]][[1]][[1]][[3]], newdata = nd, re.form = NA)

post_sum <- post %>%
  as.data.frame()
nd$i <- colnames(post_sum)
post_sum <- pivot_longer(post_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")


lbl <- c(
  `12` = "PRAQR2-B prenatal week 18"
)

p5 <- d %>% filter(week == "12") %>%
  ggplot(aes(x = praq_fear_birth_g, y = shannon)) +
  geom_smooth(
    data = post_sum, 
    aes(x = praqr_birth, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#74c476', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index 12 weeks postpartum (I3)") + xlab("PRAQR2-B (fear of giving birth)") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p5

# again but onluy cort vars 
# obtain models for plotting
stress_measures <- c("cortisone_g3", "ratio_g2")
stress_measures_m <- c("ms", "stai", "epds", "pss", "praqr_birth", "praqr_handicap")
indices <- c("shannon", "faith")
mlist <- map(stress_measures, function(stress) {
  map(indices, function(alpha) {
    
    # if ((str_detect(stress, "cort") | str_detect(stress, "ratio") & mt == "t2")) {
    #   return(NULL)
    # }
    
    # select and filter data for this RQ
    md <- map(dmimp, function(mdtemp) {
      filter(mdtemp, t == "t1") %>%
        select(id, all_of(stress_measures_m)) %>%
        mutate(id = as.character(id))
    })
    id <- map(dlongiimp, function(idtemp) {
      filter(idtemp, week == 32) %>%
        select(id, parity, feeding, gage, dmode, all_of(alpha),
               -all_of(stress_measures_i), contains("cort"))
    })
    
    # combine the data for modeling
    dtemp <- map2(md, id, function(mdtemp, idtemp) {
      full_join(mdtemp, idtemp, by = "id")
    })
    
    # model
    if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
      f <- as.formula(glue("{alpha} ~ {stress} + parity + dmode + gage + hw"))
    } else {
      f <- as.formula(glue("{alpha} ~ {stress} + parity + dmode + gage"))
    }
    
    m <- brm_multiple(
      family = student(),
      formula = f,
      data = dtemp,
      file = here(glue("data/rdata/models/iad_{stress}_{alpha}_{'t1'}_{32}"))
    )
    m

  })
})

# HCN1, faith, t1, 32

nd <- expand_grid(
  cortisone_g3 = seq(min(dlongiimp[[1]]$cortisone_g3, na.rm = TRUE), max(dlongiimp[[1]]$cortisone_g3, na.rm = TRUE), 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  #hw = 1,
  gage = median(dlongiimp[[1]]$gage, na.rm = TRUE),
  dmode = 1,
  parity = 1,
  hw = 1
  
)

post <- posterior_linpred(mlist[[1]][[2]], newdata = nd, re.form = NA)

post_sum <- post %>%
  as.data.frame()
nd$i <- colnames(post_sum)
post_sum <- pivot_longer(post_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")


lbl <- c(
  `32` = "Hair Cortisone (HCN1)"
)

p6 <- d %>% filter(week == "32") %>%
  ggplot(aes(x = cortisone_g3, y = shannon)) +
  geom_smooth(
    data = post_sum, 
    aes(x = cortisone_g3, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#238b45', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index 8 month postpartum (I4)") + xlab("Hair Cortisone (HCN1)") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p6


# hcr2, shannon, t1, 32

nd <- expand_grid(
  ratio_g2 = seq(-2.77, 0.1, 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  #hw = 1,
  gage = median(dlongiimp[[1]]$gage, na.rm = TRUE),
  dmode = 1,
  parity = 1,
  hw = 1
  
)

post <- posterior_linpred(mlist[[2]][[1]], newdata = nd, re.form = NA)

post_sum <- post %>%
  as.data.frame()
nd$i <- colnames(post_sum)
post_sum <- pivot_longer(post_sum, everything(), names_to = "i", values_to = "Estimate") %>%
  group_by(i) %>%
  summarise(
    mu = median(Estimate),
    lower = quantile(Estimate, 0.025),
    upper = quantile(Estimate, 0.975)
  ) %>%
  ungroup() %>%
  full_join(nd, by = "i")


lbl <- c(
  `32` = "Log-ratio of hair cortisol/cortisone (HCR2)"
)

p7 <- d %>% filter(week == "32") %>%
  ggplot(aes(x = ratio_g2, y = shannon)) +
  geom_smooth(
    data = post_sum, 
    aes(x = ratio_g2, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#238b45', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index 8 month postpartum (I4)") + xlab("Log-ratio of hair cortisol/cortisone (HCR2)") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

p7

save(p6, p7, file = here("data/rdata/shannon_stress_i_2.Rds"))

save(p1, p2, p3, p4, p5, file = here("data/rdata/shannon_stress_i.Rds"))

p1