set.seed(1)
library(mia)
library(tidyverse)
library(vegan)
library(here)
library(brms)
library(glue)
library(HDInterval)
library(dagitty)
library(Maaslin2)

load(file = here::here("data/rdata/tse.Rds"))
load(here("data/rdata/dm.Rds"))
# for analyses we apply prevalence fitlering
#tse <- agglomerateByRank(tse, rank = "genus")
tse <- subsetByPrevalentFeatures(tse_s, detection = 0, prevalence = 0.1)


# mothers 
tse_m <- tse[, colData(tse)$origin == "m"]
dm$id <- as.character(dm$id)
dm <- select(dm, -week, -pre)
colData(tse_m) %>% rownames()
colData(tse_m) <- colData(tse_m) %>%
                    as.data.frame() %>%
                    rownames_to_column("sid2") %>%
                    left_join(dm, by = c("id", "t")) %>%
                    column_to_rownames("sid2") %>%
                    DataFrame()


stress_indices <- c("ms", "stai", "epds", "pss", "praqr_handicap", "praqr_birth")

# store relevant data in df
d <- colData(tse_m) %>% as.data.frame() %>%
    mutate(across(all_of(c("activity", "pbmi", "age", stress_indices)), 
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
        mutate(across(all_of(c("activity", "pbmi", "age", stress_indices)), 
            function(x) scale(x)[, 1])) %>%
        mutate(stress_group = ifelse(ms <= -1, "low", ifelse(ms < 1, "medium", ifelse(ms >= 1, "high", NA))))
    dout
})




# lets start to check which species may drive the change we see between t2/t3
asv_tab <- t(assay(tse_m))
meta <- colData(tse_m) %>% as.data.frame() %>%
  filter(t != "t1")
asv_tab <- asv_tab[meta$sid,]

# you can specifiy different GLMs/normalizations/transforms. We used similar
# settings as in Nearing et al. (2021) here:
fit_data <- Maaslin2(
  asv_tab,
  meta,
  output = here::here("data/maaslin/1"),
  transform = "AST",
  fixed_effects = c("t"),
  random_effects = "id", 
  reference = "t2",  
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)

daa_tbl_1 <- filter(fit_data$results, qval <= 0.2, metadata == "t") %>%
  select(feature, coef, qval) %>%
  mutate(
    feature = str_extract(feature, "g__.*"),
    across(where(is.numeric), round, 3)
    ) %>%
  arrange(desc(abs(coef))) %>%
  rename(Feature = feature, Coefficient = coef, q = qval)
daa_tbl_1
save(daa_tbl_1, file = here("data/rdata/daa_tbl_1.Rds"))

# even though we did not find differences in beta diversity, lets check if
# individual bacteria change during pregnancy
asv_tab <- t(assay(tse_m))
meta <- colData(tse_m) %>% as.data.frame() %>%
  filter(t != "t3")
asv_tab <- asv_tab[meta$sid,]

# you can specify different GLMs/normalizations/transforms. We used similar
# settings as in Nearing et al. (2021) here:
fit_data <- Maaslin2(
  asv_tab,
  meta,
  output = here::here("data/maaslin/1"),
  transform = "AST",
  fixed_effects = c("t"),
  random_effects = "id", 
  reference = "t1",  
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)

daa_tbl_1_2 <- filter(fit_data$results, qval <= 0.2, metadata == "t") %>%
  select(feature, coef, qval) %>%
  mutate(
    feature = str_extract(feature, "g__.*"),
    across(where(is.numeric), round, 3)
  ) %>%
  arrange(desc(abs(coef))) %>%
  rename(Feature = feature, Coefficient = coef, q = qval)
daa_tbl_1_2



# now we look at stress across all time points

stress_indices <- c("ms", "stai", "epds", "pss", "praqr_handicap", "praqr_birth", "cortisol", "cortisone")
daa_tbl_2 <- map_dfr(stress_indices, function(stress) {
  coefs <- c(stress, "parity", "activity", "edu", "age", "pbmi", "t")
  if (str_detect(stress, "praqr")) {
    asv_tab <- t(assay(tse_m))
    meta <- colData(tse_m) %>% as.data.frame() %>%
      filter(t != "t3")
    asv_tab <- asv_tab[meta$sid,]
  } else if (str_detect(stress, "cort")) {
    asv_tab <- t(assay(tse_m))
    meta <- colData(tse_m) %>% as.data.frame() %>%
      filter(t != "t3")
    asv_tab <- asv_tab[meta$sid,]
    coefs <- c(stress, "parity", "activity", "edu", "age", "pbmi", "t", "hw")
  } else {
    asv_tab <- t(assay(tse_m))
    meta <- colData(tse_m) %>% as.data.frame()
    asv_tab <- asv_tab[meta$sid,]
  }
  fit_data <- Maaslin2(
    asv_tab,
    meta,
    output = here::here("data/maaslin/1"),
    transform = "AST",
    fixed_effects = coefs,
    random_effects = "id", 
    reference = c("t", "t1"),  
    normalization = "TSS",
    standardize = FALSE,
    min_prevalence = 0 # prev filterin already done
  )
  filter(fit_data$results, qval <= 0.2, metadata == stress) %>%
    mutate(stress = stress)
  }) %>%
  select(stress, feature, coef, qval) %>%
  mutate(
    feature = str_extract(feature, "g__.*"),
    across(where(is.numeric), round, 3)
  ) %>%
  arrange(desc(abs(coef))) %>%
  rename("Stress Variable" = stress, Feature = feature, Coefficient = coef, q = qval)



save(daa_tbl_2, file = here("data/rdata/daa_tbl_2.Rds"))


# do these results differ if we analyze per time point? Note that interactions
# are not possible in Maaslin2 and therefore we must explore this anyways
stress_indices <- c("ms", "stai", "epds", "pss", "praqr_handicap", "praqr_birth", "psas", "cortisol", "cortisone")
tps <- c("t1", "t2", "t3")
daa_tbl_3 <- map_dfr(stress_indices, function(stress) {
  map_dfr(tps, function(time) {
    
    coefs <- c(stress, "parity", "activity", "edu", "age", "pbmi")
    asv_tab <- t(assay(tse_m))
    meta <- colData(tse_m) %>% as.data.frame() %>%
      filter(t == time)
    if (time == "t3" & str_detect(stress, "praq")) {
      return(NULL)
    } else if(stress == "psas" & time %in% c("t1", "t2")) {
      return(NULL)
    } else if(stress == "psas" & time == "t3") {
      psas$id <- as.character(psas$id)
      meta <- colData(tse_m) %>% as.data.frame() %>%
        filter(t == time) %>%
        rownames_to_column("rn") %>%
        left_join(psas, by = "id") %>%
        column_to_rownames("rn")
    } else if(str_detect(stress, "cort") & time == "t3") {
      return(NULL)
    }
    if (str_detect(stress, "cort")) {
      coefs <- c(stress, "parity", "activity", "edu", "age", "pbmi", "hw")
    }
    asv_tab <- asv_tab[meta$sid,]
    fit_data <- Maaslin2(
      asv_tab,
      meta,
      output = here::here("data/maaslin/1"),
      transform = "AST",
      fixed_effects = coefs,
      #random_effects = "id", 
      #reference = c("t", "t1"),  
      normalization = "TSS",
      standardize = FALSE,
      min_prevalence = 0 # prev filterin already done
    )
    filter(fit_data$results, qval <= 0.2, metadata == stress) %>%
      mutate(stress = stress, t = time)
  })
}) %>%
  select(stress, t, feature, coef, qval) %>%
  mutate(
    feature = str_extract(feature, "g__.*"),
    across(where(is.numeric), round, 3)
  ) %>%
  arrange(desc(abs(coef))) %>%
  rename("Stress Variable" = stress, Feature = feature, Coefficient = coef, q = qval)

daa_tbl_2
daa_tbl_3

# in that case we do not find any associations, probably due to lack of power
save(daa_tbl_3, file = here("data/rdata/daa_tbl_3.Rds"))


