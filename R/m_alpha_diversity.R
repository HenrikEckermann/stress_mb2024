set.seed(1)
library(mia)
library(tidyverse)
library(vegan)
library(here)
library(brms)
library(glue)
library(HDInterval)
library(dagitty)

load(file = here::here("data/rdata/tse.Rds"))
load(here("data/rdata/dm.Rds"))
tse <- transformAssay(tse_s, method = "clr", name = "clr", pseudocount = 0.000001)
# mothers 
tse_m <- tse[, colData(tse)$origin == "m"]
colData(tse_m) %>% colnames()
dm$id <- as.character(dm$id)
dm <- select(dm, -week, -pre)
colData(tse_m) <- colData(tse_m) %>%
                    as.data.frame() %>%
                    left_join(dm, by = c("id", "t")) %>%
                    DataFrame()


# alpha diversity
tse_m <- estimateDiversity(
    tse_m, 
    assay.type = "counts", 
    index = "shannon")

mad_plots <- map(c("shannon", "faith"), function(alpha) {
    ylabtitle <- ifelse(alpha == "shannon", glue("{str_to_title(alpha)} index"),
                        str_to_title(alpha))
    colData(tse_m) %>%
        as.data.frame() %>%
        mutate(
          t = ifelse(t == "t3", "M3", ifelse(
            t == "t2", "M2", ifelse(
              t == "t1", "M1", NA))),
          t = factor(t, levels = c("M1", "M2", "M3"))
          ) %>%
        ggplot(aes_string("t", alpha, fill = "t")) +
        geom_boxplot(outlier.alpha = 1) +
        geom_jitter(width = 0.1, alpha = 0.6, size = 3) +
        ylab(ylabtitle) + xlab("") +
        scale_fill_manual(values = c('#fee0d2','#fc9272','#de2d26')) + 
        labs(fill = "Sample") +
        theme_bw(base_size = 15) 
})
mad_plots[[1]]
save(mad_plots, file = here("data/rdata/mad_plots.Rds"))

# store relevant data in df
d <- colData(tse_m) %>% as.data.frame() %>%
    mutate(across(all_of(c("faith", "activity", "pbmi", "age", "shannon", "ms", "stai", "epds", "pss")), 
            function(x) scale(x)[, 1])) %>%
    mutate(stress_group = ifelse(ms <= -1, "low", ifelse(ms < 1, "medium", ifelse(ms >= 1, "high", NA))))

cnames <- colnames(dm)[!colnames(dm) %in% c("id", "t")]
dmimp <- map(dmimp, function(dimp) {
    dtemp <- dimp
    dtemp$id <- as.character(dtemp$id)
    # store relevant data in df
    dout <- colData(tse_m) %>% as.data.frame() %>%
        select(-all_of(c(cnames, "week", "pre"))) %>%
        left_join(dtemp, by = c("id", "t")) %>%
        mutate(across(all_of(c("faith", "activity", "pbmi", "age", "shannon", "ms", "stai", "epds", "pss", "cortisol", "cortisone")), 
            function(x) scale(x)[, 1])) %>%
        mutate(stress_group = ifelse(ms <= -1, "low", ifelse(ms < 1, "medium", ifelse(ms >= 1, "high", NA))))
    dout
})



# looks like alpha diversity 8 months postpartum is lower compared to within preg
# can we confirm this statistically?
# Before we dive into that, we first need to answer several questions:

# 1) Our DAG includes the following covariates: parity, pBMI, edu, activity,
# gestational age (time point), maternal stress, antibiotics & diet.
# Some of these covariates are interesting to look at because of previous
# research (parity and pBMI) while a certain selection must be included/excluded
# in order to focus on our research question MS --> alpha diversity, which ones
# are the covariates that we MUST include?
colnames(d)
g <- dagitty('dag {
  bb="0,0,1,1"
  abx [pos="0.362,0.450"]
  activity [pos="0.101,0.164"]
  age [pos="0.033,0.061"]
  dhd15_total [pos="0.556,0.146"]
  edu [pos="0.123,0.353"]
  ms [pos="0.341,0.204"]
  parity [pos="0.063,0.291"]
  pbmi [pos="0.358,0.373"]
  shannon [pos="0.729,0.178"]
  t [pos="0.345,0.293"]
  abx -> shannon
  activity -> abx
  activity -> ms
  activity -> shannon
  age -> activity
  age -> edu
  age -> ms
  age -> parity
  age -> shannon
  dhd15_total -> shannon
  edu -> abx
  edu -> activity
  edu -> dhd15_total
  edu -> ms
  edu -> pbmi
  ms -> dhd15_total
  ms -> shannon
  parity -> activity
  parity -> ms
  parity -> pbmi
  parity -> shannon
  pbmi -> shannon
  t -> dhd15_total
  t -> shannon
  }
'
)

# we must include activity, education and parity minimally
adjustmentSets(g, exposure = "ms", outcome = "shannon", type = "minimal")
# which ones are optional
adjustmentSets(g, exposure = "ms", outcome = "shannon", type = "all")
# if this DAG is correct, we can include all these:
# { abx, activity, age, edu, parity, pbmi, t } based on model fit . But first I check CI:
impliedConditionalIndependencies(g)
# abx _||_ ms | actv, age, edu, prty
summary(glm(abx ~ ms + activity + age + edu + parity, data = d))
# abx _||_ t
summary(glm(abx ~ t, data = d))
# actv _||_ t
summary(lm(activity ~ t, data = d))
# age _||_ pbmi | edu, prty
summary(lm(age ~ pbmi + edu + parity, data = d))
# age _||_ t
# edu _||_ prty | age
summary(glm(edu ~ parity + age, data = mutate(d, edu = ifelse(edu == "high", 1, 0))))
# edu _||_ shnn | abx, actv, age, d15_, ms, prty, pbmi, t
summary(glm(edu ~ shannon + abx + activity + dhd15_total + ms + parity + age + pbmi, data = mutate(d, edu = ifelse(edu == "high", 1, 0))))
# edu _||_ t
summary(glm(edu ~ t, data = mutate(d, edu = ifelse(edu == "high", 1, 0))))
# ms _||_ pbmi | actv, age, edu, prty
summary(lm(ms ~ pbmi + edu + parity + age + activity, data = d))
# ms _||_ t
summary(lm(ms ~ t, data = d))
# prty _||_ t
summary(glm(parity ~ t, data = mutate(d, parity = ifelse(parity == "1", 1, 0))))
# pbmi _||_ t
summary(lm(pbmi ~ t, data = d))

# the DAG has been adjusted based on testing the conditional independences. 


# question 2 & 3 have to do with model structure. 2) is it better to treat the microbiota at the 
# different time points as different system (for each time point fit a model) or should we fit a 
# multi level model. I will try both and depending on which model makes better predictions (msqe)
# we will use that for inference. 


# use 10 fold crossvalidation
folds <- caret::createFolds(d$sid, k = 10)

# first the multilevel model
m1_mm_msqe <- map2_dbl(folds, 1:10, function(fold, k) {

    # splits
    train <- d[-fold, ]
    test <- d[fold, ]

    # fit model 
    m1_mm <- brm(
        family = student(),
        formula = shannon ~ activity + age + edu + parity + t + (1|id),
        data = train,
        file = here(glue("data/rdata/models/m1_multilevel_{k}"))
    )

    pred <- predict(m1_mm, newdata = test, re_formula = NA) %>%
      as.data.frame() %>%
      .$Estimate
    actual <- test$shannon
    # mse
    squared_diff <- (actual - pred)^2
    mse <- mean(squared_diff, na.rm = TRUE)
    mse
})
mean(m1_mm_msqe)

# now the separate models
m2_mm_msqe <- map2_dbl(folds, 1:10, function(fold, k) {

    # splits
    train <- d[-fold, ]
    test <- d[fold, ]

    # fit models
    m1 <- brm(
        family = student(),
        formula = shannon ~ activity + age + edu + parity,
        data = filter(train, t == "t1"),
        file = here(glue("data/rdata/models/m1_multilevel_t1_{k}"))
    )

    m2 <- brm(
        family = student(),
        formula = shannon ~ activity + age + edu + parity,
        data = filter(train, t == "t2"),
        file = here(glue("data/rdata/models/m1_multilevel_t2_{k}"))
    )

    m3 <- brm(
        family = student(),
        formula = shannon ~ age + edu + parity,
        data = filter(train, t == "t3"),
        file = here(glue("data/rdata/models/m1_multilevel_t3_{k}"))
    )

    # get mses
    pred1 <- predict(m1, newdata = filter(test, t == "t1"), re_formula = NA) %>%
        as.data.frame() %>%
        .$Estimate
    actual1 <- filter(test, t == "t1") %>% .$shannon
    squared_diff1 <- (actual1 - pred1)^2
    mse1 <- mean(squared_diff1, na.rm = TRUE)

    pred2 <- predict(m2, newdata = filter(test, t == "t2"), re_formula = NA) %>%
        as.data.frame() %>%
        .$Estimate
    actual2 <- filter(test, t == "t2") %>% .$shannon
    squared_diff2 <- (actual2 - pred2)^2
    mse2 <- mean(squared_diff2, na.rm = TRUE)

    pred3 <- predict(m3, newdata = filter(test, t == "t3"), re_formula = NA) %>%
        as.data.frame() %>%
        .$Estimate
    actual3 <- filter(test, t == "t3") %>% .$shannon
    squared_diff3 <- (actual3 - pred3)^2
    mse3 <- mean(squared_diff3, na.rm = TRUE)

    mean(c(squared_diff1, squared_diff2, squared_diff3), na.rm = TRUE)

})

# the multilevel model beats the single models
mean(m1_mm_msqe)
sd(m1_mm_msqe)
mean(m2_mm_msqe)
sd(m2_mm_msqe)


# last step before we continue is: 3) which of the optional covariates improve model fit 
# and should therefore be included. These are obligated: { activity, age, edu, parity }
# optional: { abx, pbmi}

loo_comp <- map_dfr(1:5, function(m) {
    m_base <- brm(
        formula = shannon ~ t + parity + activity + edu + age + (1|id),
        data = dmimp[[m]],
        file = glue(here("data/rdata/models/m_base_{m}"))
    )
    loo_mbase <- add_criterion(
        m_base,
        "loo",
        file = glue(here("data/rdata/models/modelsloo_mbase_{m}")),
        moment_match = FALSE
    )

    m_1 <- brm(
        formula = shannon ~ t + parity + activity + edu + age + pbmi + (1|id),
        data = dmimp[[m]],
        file = glue(here("data/rdata/models/m_1_{m}"))
    )
    loo_1 <- add_criterion(
        m_1,
        "loo",
        file = glue(here("data/rdata/models/modelsloo_m1_{m}")),
        moment_match = FALSE
    )


    m_2 <- brm(
        formula = shannon ~ t + parity + activity + edu + age + abx + (1|id),
        data = dmimp[[m]],
        file = glue(here("data/rdata/models/m_2_{m}"))
    )
    loo_2 <- add_criterion(
        m_2,
        "loo",
        file = glue(here("data/rdata/models/modelsloo_m2_{m}")),
        moment_match = FALSE
    )

    m_3 <- brm(
        formula = shannon ~ t * parity * pbmi + activity + edu + age + (1|id),
        data = dmimp[[m]],
        file = glue(here("data/rdata/models/m_3_{m}"))
    )
    loo_3 <- add_criterion(
        m_3,
        "loo",
        file = glue(here("data/rdata/models/modelsloo_m3_{m}")),
        moment_match = FALSE
    )

    m_4 <- brm(
        formula = shannon ~ t * parity * pbmi + activity + edu + age + abx + (1|id),
        data = dmimp[[m]],
        file = glue(here("data/rdata/models/m_4_{m}"))
    )
    loo_4 <- add_criterion(
        m_4,
        "loo",
        file = glue(here("data/rdata/models/modelsloo_m4_{m}")),
        moment_match = FALSE
    )



    tbl <- loo_compare(loo_mbase, loo_1, loo_2, loo_3, loo_4)
    tbl %>% as.data.frame() %>%
        mutate(order = 1:dim(as.data.frame(tbl))) %>%
        rownames_to_column("model")
})

count(loo_comp, order, model)
# the results suggest to use model 1 or the base model. Even though we will use therefore model 1
# to evaluate our research questions regarding ms, we still perform sensitivity analyses and exploratory
# analyses to see if our data support what previous studies found.



# first lets just look at change over time to add data to kennedy et al
mad1 <- map(c("shannon", "faith"), function(alpha) {
    form <- as.formula(glue("{alpha} ~ t + parity + activity + edu + age + pbmi + (1|id)"))
    brm_multiple(
        family = student(),
        formula = form,
        data = dmimp,
        file = here(glue("data/rdata/models/mad1_{alpha}_t"))
    )
})

map(mad1, ~summary(.x))

mad3 <- map(c("shannon", "faith"), function(alpha) {
    form <- as.formula(glue("{alpha} ~ t * parity * pbmi + activity + edu + age + (1|id)"))
    brm_multiple(
        family = student(),
        formula = form,
        data = dmimp,
        file = here(glue("data/rdata/models/mad3_{alpha}_t"))
    )
})

map(mad3, ~summary(.x))

# what is the average change over time across different levels of covariates?
b <- map_dfr(c("shannon", "faith"), function(alpha) {
    nd <- expand.grid(
        t = c("t1", "t2", "t3"),
        parity = c("0", "1"), 
        edu = c("low", "high"),
        age = median(d$age, na.rm = TRUE),
        pbmi = median(d$pbmi, na.rm = TRUE),
        activity = median(d$activity, na.rm = TRUE),
        id = NA
    )
    temp_join <- select(nd, t, parity, edu) %>% mutate(V = glue("V{1:12}"))
    temp_join
    epred <- posterior_epred(
        mad1[[ifelse(alpha == "shannon", 1, 2)]],
        re.formula = NA,
        newdata = nd) %>% as.data.frame() %>%
        pivot_longer(everything(), names_to = "V", values_to = "mu") %>%
        left_join(temp_join, by = "V")
    tail(epred, 20)
    # now we can calculate any contrast we like:
    # t2- t1
    t1 <- filter(epred, t == "t1") %>% .$mu
    t2 <- filter(epred, t == "t2") %>% .$mu
    t3 <- filter(epred, t == "t3") %>% .$mu

    t2_t1 <- t2 - t1 
    t3_t1 <- t3 - t1 
    t3_t2 <- t3 - t2 

    map2_dfr(list(t2_t1, t3_t1, t3_t2), c("t2 - t1", "t3 - t1", "t3 - t2"), function(contrast, name) {
        m <- median(contrast)
        lower <- hdi(contrast)[1]
        upper <- hdi(contrast)[2]
        p <- mean(contrast > 0)
        tibble(alpha = alpha, contrast = name, m, lower, upper, p) %>%
            mutate(across(where(is.numeric), round, 3))
    })

})
b
save(b, file = here("data/rdata/b_change.Rds"))

# even if the interaction of pbmi and parity does not lead to better model fit, do we see a similar 
# trend as kennedy et al when looking at the effect size?
# what is the average change over time across different levels of covariates?
b <- map_dfr(c("shannon", "faith"), function(alpha) {
    nd <- expand.grid(
        t = c("t1", "t2", "t3"),
        parity = c("0", "1"), 
        edu = c("low", "high"),
        age = median(d$age, na.rm = TRUE),
        pbmi = c(quantile(d$pbmi, 0.1, na.rm = TRUE), 
                median(d$pbmi, na.rm = TRUE), 
                quantile(d$pbmi, 0.9, na.rm = TRUE)
                ),
        activity = median(d$activity, na.rm = TRUE),
        id = NA
    )

    temp_join <- select(nd, t, parity, edu, pbmi) %>% mutate(V = glue("V{1:dim(nd)[1]}"))
    epred <- posterior_epred(
        mad3[[ifelse(alpha == "shannon", 1, 2)]],
        re.formula = NA,
        newdata = nd) %>% as.data.frame() %>%
        pivot_longer(everything(), names_to = "V", values_to = "mu") %>%
        left_join(temp_join, by = "V")
    print(epred, n =  36)


    # now we can calculate any contrast we like:
    # is there an interaction that ad changes over time, but only in the low bmi group?
    t1_lbmi <- filter(epred, t == "t1", pbmi == quantile(d$pbmi, 0.1, na.rm = TRUE)) %>% .$mu
    t2_lbmi <- filter(epred, t == "t2", pbmi == quantile(d$pbmi, 0.1, na.rm = TRUE)) %>% .$mu
    t3_lbmi <- filter(epred, t == "t3", pbmi == quantile(d$pbmi, 0.1, na.rm = TRUE)) %>% .$mu

    t1_hbmi <- filter(epred, t == "t1", pbmi == quantile(d$pbmi, 0.9, na.rm = TRUE)) %>% .$mu
    t2_hbmi <- filter(epred, t == "t2", pbmi == quantile(d$pbmi, 0.9, na.rm = TRUE)) %>% .$mu
    t3_hbmi <- filter(epred, t == "t3", pbmi == quantile(d$pbmi, 0.9, na.rm = TRUE)) %>% .$mu

    t2_t1_lbmi <- t2_lbmi - t1_lbmi 
    t3_t1_lbmi <- t3_lbmi - t1_lbmi 
    t3_t2_lbmi <- t3_lbmi - t2_lbmi 

    t2_t1_hbmi <- t2_hbmi - t1_hbmi 
    t3_t1_hbmi <- t3_hbmi - t1_hbmi 
    t3_t2_hbmi <- t3_hbmi - t2_hbmi 

    t1_lbmi_pairty0 <- filter(epred, t == "t1", pbmi == quantile(d$pbmi, 0.1, na.rm = TRUE), parity == "0") %>% .$mu
    t2_lbmi_pairty0 <- filter(epred, t == "t2", pbmi == quantile(d$pbmi, 0.1, na.rm = TRUE), parity == "0") %>% .$mu
    t3_lbmi_pairty0 <- filter(epred, t == "t3", pbmi == quantile(d$pbmi, 0.1, na.rm = TRUE), parity == "0") %>% .$mu

    t1_lbmi_pairty1 <- filter(epred, t == "t1", pbmi == quantile(d$pbmi, 0.1, na.rm = TRUE), parity == "1") %>% .$mu
    t2_lbmi_pairty1 <- filter(epred, t == "t2", pbmi == quantile(d$pbmi, 0.1, na.rm = TRUE), parity == "1") %>% .$mu
    t3_lbmi_pairty1 <- filter(epred, t == "t3", pbmi == quantile(d$pbmi, 0.1, na.rm = TRUE), parity == "1") %>% .$mu

    t1_hbmi_pairty0 <- filter(epred, t == "t1", pbmi == quantile(d$pbmi, 0.9, na.rm = TRUE), parity == "0") %>% .$mu
    t2_hbmi_pairty0 <- filter(epred, t == "t2", pbmi == quantile(d$pbmi, 0.9, na.rm = TRUE), parity == "0") %>% .$mu
    t3_hbmi_pairty0 <- filter(epred, t == "t3", pbmi == quantile(d$pbmi, 0.9, na.rm = TRUE), parity == "0") %>% .$mu

    t1_hbmi_pairty1 <- filter(epred, t == "t1", pbmi == quantile(d$pbmi, 0.9, na.rm = TRUE), parity == "1") %>% .$mu
    t2_hbmi_pairty1 <- filter(epred, t == "t2", pbmi == quantile(d$pbmi, 0.9, na.rm = TRUE), parity == "1") %>% .$mu
    t3_hbmi_pairty1 <- filter(epred, t == "t3", pbmi == quantile(d$pbmi, 0.9, na.rm = TRUE), parity == "1") %>% .$mu

    t2_t1_lbmi_pairty0 <- t2_lbmi_pairty0 - t1_lbmi_pairty0 
    t3_t1_lbmi_pairty0 <- t3_lbmi_pairty0 - t1_lbmi_pairty0 
    t3_t2_lbmi_pairty0 <- t3_lbmi_pairty0 - t2_lbmi_pairty0 

    t2_t1_lbmi_pairty1 <- t2_lbmi_pairty1 - t1_lbmi_pairty1 
    t3_t1_lbmi_pairty1 <- t3_lbmi_pairty1 - t1_lbmi_pairty1 
    t3_t2_lbmi_pairty1 <- t3_lbmi_pairty1 - t2_lbmi_pairty1 

    t2_t1_hbmi_pairty0 <- t2_hbmi_pairty0 - t1_hbmi_pairty0 
    t3_t1_hbmi_pairty0 <- t3_hbmi_pairty0 - t1_hbmi_pairty0 
    t3_t2_hbmi_pairty0 <- t3_hbmi_pairty0 - t2_hbmi_pairty0 

    t2_t1_hbmi_pairty1 <- t2_hbmi_pairty1 - t1_hbmi_pairty1 
    t3_t1_hbmi_pairty1 <- t3_hbmi_pairty1 - t1_hbmi_pairty1 
    t3_t2_hbmi_pairty1 <- t3_hbmi_pairty1 - t2_hbmi_pairty1 

    # does pmbi have different effects depending on the timepoint?
    t1_lbmi_hbmi <- t1_lbmi - t1_hbmi 
    t2_lbmi_hbmi <- t2_lbmi - t2_hbmi 
    t3_lbmi_hbmi <- t3_lbmi - t3_hbmi 


    contrasts <- map2_dfr(list(
        t2_t1_lbmi, t3_t1_lbmi, t3_t2_lbmi, 
        t2_t1_hbmi, t3_t1_hbmi, t3_t2_hbmi,

        t2_t1_lbmi_pairty0,
        t3_t1_lbmi_pairty0,
        t3_t2_lbmi_pairty0,
        t2_t1_lbmi_pairty1,
        t3_t1_lbmi_pairty1,
        t3_t2_lbmi_pairty1,
        t2_t1_hbmi_pairty0,
        t3_t1_hbmi_pairty0,
        t3_t2_hbmi_pairty0,
        t2_t1_hbmi_pairty1,
        t3_t1_hbmi_pairty1,
        t3_t2_hbmi_pairty1,

        t1_lbmi_hbmi <- t1_lbmi - t1_hbmi,
        t2_lbmi_hbmi <- t2_lbmi - t2_hbmi,
        t3_lbmi_hbmi <- t3_lbmi - t3_hbmi 
        ), 
        c(
            "t2 - t1 low bmi", "t3 - t1 low bmi", "t3 - t2 low bmi",
            "t2 - t1 high bmi", "t3 - t1 high bmi", "t3 - t2 high bmi",

            "t2 - t1 low bmi parity0", "t3 - t1 low bmi parity0", "t3 - t2 low bmi parity0",
            "t2 - t1 low bmi parity1", "t3 - t1 low bmi parity1", "t3 - t2 low bmi parity1",

            "t2 - t1 high bmi parity0", "t3 - t1 high bmi parity0", "t3 - t2 high bmi parity0",
            "t2 - t1 high bmi parity1", "t3 - t1 high bmi parity1", "t3 - t2 high bmi parity1",

            "t1 low bmi - high bmi", "t2 low bmi - high bmi", "t3 low bmi - high bmi"
            ), function(contrast, name) {
        m <- median(contrast)
        lower <- hdi(contrast)[1]
        upper <- hdi(contrast)[2]
        p <- mean(contrast > 0)
        tibble(alpha = alpha, contrast = name, m, lower, upper, p) %>%
            mutate(across(where(is.numeric), round, 3))
    })

})
b
print(b, n = 42)
save(b, file = here("data/rdata/b_kennedy.Rds"))
# exploratory contrasts indidcate that the change in diversity is highest in mothers that had low bmi and
# it was their first child. But unlike kennedy we talk about the shift to postpartum. There was no 
# change from t1 to t2 really across any contrasts I investigated. 

# we did not observe a correlation between pbmi and alpha diversity and our effect sizes point more 
# to a relation in the other direction, but keep in mind that we did not add GWG to the model, which
# would change the quesiton to: keeping GWG constant, does pbmi relate to ad. We ask, does pbmi related
# to ad independent of GWG?






# lets see if maternal stress has an influence on alpha diversity over time?



mad_ms <- map(c("ms", "stai", "epds", "pss", "praqr_handicap", "praqr_birth"), function(stress) {
    map(c("shannon", "faith"), function(alpha) {
      if (stress %in% c("praqr_handicap", "praqr_birth")) {
        dimptemp <- map(dmimp, ~filter(.x, t != "t3"))
      } else {
        dimptemp <- dmimp
      }
      
        form <- as.formula(glue("{alpha} ~ t * {stress} + parity + activity + edu + age + pbmi + (1|id)"))
        brm_multiple(
            family = student(),
            formula = form,
            data = dimptemp,
            file = here(glue("data/rdata/models/mad_{stress}_{alpha}_t"))
        )
    })
})

# from the results printed below that if there is a linear relationship, it is likely a positive
# association and it is strongest at T3. However, we cannot state any of that with confidence. 
# furthermore, the effects are similar between stress measurements
source(here("R/helper_functions.R"))

b <- map2_dfr(mad_ms, c("ms", "stai", "epds", "pss", "praqr_handicap", "praqr_birth"), function(l, stress) {
    map2_dfr(l, c("shannon", "faith"), function(x, ad) {
      params <- colnames(x$rhats)[str_detect(colnames(x$rhats), "b_")]
      params <- str_replace(params, "\\.", ":")
      params <- params[str_detect(params, stress)]
      
      post1 <- as_draws_df(x) %>%
        select(all_of(params))
      temp1 <- glue("b_{stress}")
      temp2 <- glue("b_tt2:{stress}")
      temp3 <- glue("b_tt3:{stress}")
      post2 <- tibble(
        post1[[temp1]],
        post1[[temp1]] + post1[[temp2]],
        post1[[temp1]] + post1[[temp3]]
      )
      colnames(post2) <- colnames(post1)

      summarise_posterior2(post2, parameters = params) %>%
        mutate(alpha = ad, stress = stress) %>%
        filter(str_detect(parameter, stress))
  })
})

print(b, n = 32)
# i need an extra copy for the table later
btemp <- b
save(b, file = here("data/rdata/b_stress_ad_m.Rds"))



# we explore the same now using low and high stress groups
mad_ms_group <- map(c("shannon", "faith"), function(alpha) {
    form <- as.formula(glue("{alpha} ~ t * stress_group + parity + activity + edu + age + pbmi + (1|id)"))
    brm_multiple(
        family = student(),
        formula = form,
        data = dmimp,
        file = here(glue("data/rdata/models/mad_stress_group_{alpha}_t"))
    )
})


b <- map_dfr(c("shannon", "faith"), function(alpha) {
    nd <- expand.grid(
        t = c("t1", "t2", "t3"),
        parity = c("0", "1"), 
        edu = c(0, 1),
        age = median(d$age, na.rm = TRUE),
        pbmi = median(d$pbmi, na.rm = TRUE),
        activity = median(d$activity, na.rm = TRUE),
        stress_group = c("low", "medium", "high"),
        id = NA
    )
    temp_join <- select(nd, t, parity, edu, stress_group) %>% mutate(V = glue("V{1:dim(nd)[1]}"))
    epred <- posterior_epred(
        mad_ms_group[[ifelse(alpha == "shannon", 1, 2)]],
        re.formula = NA,
        newdata = nd) %>% as.data.frame() %>%
        pivot_longer(everything(), names_to = "V", values_to = "mu") %>%
        left_join(temp_join, by = "V")
    # now we can calculate any contrast we like:
    # main effects stress group
    l <- filter(epred, stress_group == "low") %>% .$mu
    m <- filter(epred, stress_group == "medium") %>% .$mu
    h <- filter(epred, stress_group == "high") %>% .$mu

    # per time point 
    lt1 <- filter(epred, t == "t1", stress_group == "low") %>% .$mu
    lt2 <- filter(epred, t == "t2", stress_group == "low") %>% .$mu
    lt3 <- filter(epred, t == "t3", stress_group == "low") %>% .$mu

    mt1 <- filter(epred, t == "t1", stress_group == "medium") %>% .$mu
    mt2 <- filter(epred, t == "t2", stress_group == "medium") %>% .$mu
    mt3 <- filter(epred, t == "t3", stress_group == "medium") %>% .$mu 

    ht1 <- filter(epred, t == "t1", stress_group == "high") %>% .$mu
    ht2 <- filter(epred, t == "t2", stress_group == "high") %>% .$mu
    ht3 <- filter(epred, t == "t3", stress_group == "high") %>% .$mu

    # calculate contrasts
    lm <- l - m
    lh <- l - h 
    mh <- m - h 

    lmt1 <- lt1 - mt1 
    lht1 <- lt1 - ht1 
    mht1 <- mt1 - ht1 

    lmt2 <- lt2 - mt2 
    lht2 <- lt2 - ht2 
    mht2 <- mt2 - ht2 

    lmt3 <- lt3 - mt3 
    lht3 <- lt3 - ht3 
    mht3 <- mt3 - ht3 

    map2_dfr(
        list(
                lm,
                lh, 
                mh, 

                lmt1, 
                lht1, 
                mht1, 

                lmt2, 
                lht2, 
                mht2, 

                lmt3, 
                lht3, 
                mht3 
            ), 
        c(
            "low - medium", "low - high", "medium - high",
            "low - medium t1", "low - high t1", "medium - high t1",
            "low - medium t2", "low - high t2", "medium - high t2",
            "low - medium t3", "low - high t3", "medium - high t3"
            ), 
        function(contrast, name) {
        m <- median(contrast)
        lower <- hdi(contrast)[1]
        upper <- hdi(contrast)[2]
        p <- mean(contrast > 0)
        tibble(alpha = alpha, contrast = name, m, lower, upper, p) %>%
            mutate(across(where(is.numeric), round, 3))
    })

})
print(b, n = 24)
save(b, file = here("data/rdata/b_stress_group_m.Rds"))

# we found that at t3, the low group had lower shannon than the other groups
# same effect direction but not significant for Faith. This probably responsible
# for the slight positive trend we found in the continuous models

# PRAQ, PES and PSAS (i left this code in but I actually added praq above for the
# supplementary table in the end)


mad_ms <- map(c("praqr_handicap", "praqr_birth"), function(stress) {
  map(c("shannon", "faith"), function(alpha) {
    form <- as.formula(glue("{alpha} ~ t * {stress} + parity + activity + edu + age + pbmi + (1|id)"))
    brm_multiple(
      family = student(),
      formula = form,
      data = map(dmimp, ~filter(.x, t != "t3")),
      file = here(glue("data/rdata/models/mad_{stress}_{alpha}_t"))
    )
  })
})


map(mad_ms, function(l) {
  map(l, ~summary(.x))
})

b_praq <- map2_dfr(mad_ms, c("praqr_handicap", "praqr_birth"), function(mlist, stress) {
  map2_dfr(mlist, c("shannon", "faith"), function(model, alpha) {
    summarise_posterior(model, parameters = glue("b_{stress}")) %>%
      mutate(alpha = alpha)
  })
})
b_praq
# no effect of praw handicap, but for praqr_birth, the model 
# is certain that alpha diversity is decreased in high stress mothers but only
# at 18 weeks of pregnancy. The effect is visible for shannong and faith

save(b_praq, file = here("data/rdata/b_praq_m.Rds"))

# make proper bayesian plot for this relation
nd <- expand_grid(
  praqr_birth = seq(min(dmimp[[1]]$praqr_birth[dmimp[[1]]$week != "32"], na.rm = TRUE), max(dmimp[[1]]$praqr_birth[dmimp[[1]]$week != "32"], na.rm = TRUE), 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  #hw = 1,
  parity = 1,
  activity = median(dmimp[[1]]$activity, na.rm = TRUE),
  edu = median(dmimp[[1]]$edu, na.rm = TRUE),
  age = median(dmimp[[1]]$age, na.rm = TRUE),
  pbmi = median(dmimp[[1]]$pbmi, na.rm = TRUE),
  t = "t1"
)

post <- posterior_linpred(mad_ms[[2]][[1]], newdata = nd, re.form = NA)

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
  `18` = "M1"
)



shannon_m1 <- d %>% filter(week != "32") %>%
  ggplot(aes(x = praqr_birth, y = shannon)) +
  geom_smooth(
    data = post_sum, 
    aes(x = praqr_birth, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#fee0d2', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index") + xlab("PRAQR2-B (fear of giving birth)") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

shannon_m1



mad_ms <- map(c("shannon", "faith"), function(alpha) {
  form <- as.formula(glue("{alpha} ~ psas + parity + activity + edu + age + pbmi"))
  brm_multiple(
    family = student(),
    formula = form,
    data = psas_dmimp,
    file = here(glue("data/rdata/models/mad_psas_{alpha}_t"))
  )
})

map(mad_ms, ~summary(.x))



b_psas <- map2_dfr(mad_ms, c("shannon", "faith"), function(model, alpha) {
  summarise_posterior(model, parameters = glue("b_psas")) %>%
    mutate(alpha = alpha)
})
b_psas
save(b_psas, file = here("data/rdata/b_psas_m.Rds"))







mad_ms_group <- map(c("shannon", "faith"), function(alpha) {
  form <- as.formula(glue("{alpha} ~ stress_group + parity + activity + edu + age + pbmi"))
  brm_multiple(
    family = student(),
    formula = form,
    data = psas_dmimp,
    file = here(glue("data/rdata/models/mad_psas_group_{alpha}_t"))
  )
})






b <- map_dfr(c("shannon", "faith"), function(alpha) {
  nd <- expand.grid(
    parity = c("0", "1"), 
    edu = c(0, 1),
    age = median(d$age, na.rm = TRUE),
    pbmi = median(d$pbmi, na.rm = TRUE),
    activity = median(d$activity, na.rm = TRUE),
    stress_group = c("low", "medium", "high"),
    id = NA
  )
  temp_join <- select(nd, parity, edu, stress_group) %>% mutate(V = glue("V{1:dim(nd)[1]}"))
  epred <- posterior_epred(
    mad_ms_group[[ifelse(alpha == "shannon", 1, 2)]],
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
      tibble(alpha = alpha, contrast = name, m, lower, upper, p) %>%
        mutate(across(where(is.numeric), round, 3))
    })
  
})
b
save(b, file = here("data/rdata/b_stress_group_m.Rds"))





# lastly, cortisol and cortisone 
mad_ms <- map(c("cortisol", "cortisone", "ratio"), function(stress) {
  map(c("shannon", "faith"), function(alpha) {
    form <- as.formula(glue("{alpha} ~ t * {stress} + parity + activity + edu + age + pbmi + hw + (1|id)"))
    brm_multiple(
      family = student(),
      formula = form,
      data = map(dmimp, ~filter(.x, t != "t3")),
      file = here(glue("data/rdata/models/mad_{stress}_{alpha}_t"))
    )
  })
})


b <- map2_dfr(mad_ms, c("cortisol", "cortisone", "ratio"), function(l, stress) {
  map2_dfr(l, c("shannon", "faith"), function(x, ad) {
    params <- colnames(x$rhats)[str_detect(colnames(x$rhats), "b_")]
    params <- str_replace(params, "\\.", ":")
    params <- params[str_detect(params, stress)]
    
    post1 <- as_draws_df(x) %>%
      select(all_of(params))
    temp1 <- glue("b_{stress}")
    temp2 <- glue("b_tt2:{stress}")
    post2 <- tibble(
      post1[[temp1]],
      post1[[temp1]] + post1[[temp2]]
    )
    colnames(post2) <- colnames(post1)
    
    summarise_posterior2(post2, parameters = params) %>%
      mutate(alpha = ad, stress = stress) %>%
      filter(str_detect(parameter, stress))
  })
})



# make supplementary table from here:
bm_supp <- bind_rows(btemp, b, b_psas, ) 
save(bm_supp, file = here("data/rdata/bm_supp.Rds"))


# indication that only at t1 hair cortisol is associated, cortisone potentially negative
b
save(b, file = here("data/rdata/b_cort_ad_m.Rds"))


lbl <- c(
  `t3` = "M3"
)

shannon_m2 <- filter(d, t == "t3", !is.na(ms)) %>%
  mutate(stress_group = ifelse(ms <= -1, "low", ifelse(ms < 1, "medium", ifelse(ms >= 1, "high", NA)))) %>%
  mutate(stress_group = factor(stress_group, levels = c("low", "medium", "high"))) %>%
  ggplot(aes(stress_group, shannon)) +
  geom_boxplot(outlier.alpha = 0, fill = '#de2d26') +
  geom_jitter(width = 0.2, size = 3) +
  facet_wrap(~t, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index") + xlab("Maternal Stress Composite") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

shannon_m2



# make proper bayesian plot for log ratio relations
nd <- expand_grid(
  ratio = seq(min(dmimp[[1]]$ratio[dmimp[[1]]$week == 18], na.rm = TRUE), max(dmimp[[1]]$ratio[dmimp[[1]]$week == 18], na.rm = TRUE), 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  hw = 1,
  parity = 1,
  activity = median(dmimp[[1]]$activity, na.rm = TRUE),
  edu = median(dmimp[[1]]$edu, na.rm = TRUE),
  age = median(dmimp[[1]]$age, na.rm = TRUE),
  pbmi = median(dmimp[[1]]$pbmi, na.rm = TRUE),
  t = "t1"
)


post <- posterior_linpred(mad_ms[[3]][[1]], newdata = nd, re.form = NA)

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
  `18` = "M1"
)



shannon_m3 <- d %>% filter(week == "18", !is.na(ratio)) %>%
  ggplot(aes(x = ratio, y = shannon)) +
  geom_smooth(
    data = post_sum, 
    aes(x = ratio, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#fee0d2', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Shannon index") + xlab("Log-ratio of cortisol/cortisone (HCR2)") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )

shannon_m3


nd <- expand_grid(
  ratio = seq(min(dmimp[[1]]$ratio[dmimp[[1]]$week == 18], na.rm = TRUE), max(dmimp[[1]]$ratio[dmimp[[1]]$week == 18], na.rm = TRUE), 0.1),
  #pss = seq(-2, 2, 0.1),
  #ratio_g2 = seq(-2, 2, 0.1),
  id = NA,
  hw = 1,
  parity = 1,
  activity = median(dmimp[[1]]$activity, na.rm = TRUE),
  edu = median(dmimp[[1]]$edu, na.rm = TRUE),
  age = median(dmimp[[1]]$age, na.rm = TRUE),
  pbmi = median(dmimp[[1]]$pbmi, na.rm = TRUE),
  t = "t1"
)


post <- posterior_linpred(mad_ms[[3]][[2]], newdata = nd, re.form = NA)

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
  `18` = "M1"
)



shannon_m4 <- d %>% filter(week == "18", !is.na(ratio)) %>%
  ggplot(aes(x = ratio, y = faith)) +
  geom_smooth(
    data = post_sum, 
    aes(x = ratio, y = mu, ymin = lower, ymax = upper), 
    stat = "identity",
    fill = '#fee0d2', color = "black", alpha = 4/5, linewidth = 1/4
  ) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1) +
  facet_wrap(~week, labeller = as_labeller(lbl), strip.position = "top") +
  ylab("Faith") + xlab("Log-ratio of cortisol/cortisone (HCR2)") +
  ylim(c(-3.5, 3.5)) +
  theme_bw(base_size = 15) +
  theme(
    strip.placement = "inside",
    strip.text.x = element_text(size = 14)
  )
shannon_m3
shannon_m4

save(shannon_m1, shannon_m2, shannon_m3, shannon_m4, file = here("data/rdata/shannon_m_p.Rds"))
