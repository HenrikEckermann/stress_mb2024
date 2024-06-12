set.seed(1)
library(mia)
library(tidyverse)
library(vegan)
library(here)
library(glue)
library(Maaslin2)

load(file = here::here("data/rdata/tse.Rds"))
load(here("data/rdata/di.Rds"))
# for analyses we apply prevalence filtering
#tse <- agglomerateByRank(tse, rank = "genus")
tse <- subsetByPrevalentFeatures(tse_s, detection = 0, prevalence = 0.1)
tse_i <- tse[, colData(tse)$origin == "i"]
# join all relevant data
dlong$week <- as.numeric(dlong$week)
dlong$id <- as.character(dlong$id)

# the year samples didnt arrive on time for this project
tse_i <- tse_i[, tse_i$week != 104]

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


# store relevant data in df
df <- colData(tse_i) %>% as.data.frame() %>%
  mutate(across(
    all_of(c("gage", "psas", 
             "stai", "epds", "pss")),
    function(x) scale(x)[, 1])) %>%
  mutate(
    ms = scale(epds + stai + pss)[, 1], 
    across(all_of(c(contains("cort"), contains("praq"), contains("ratio"))), function(x) scale(x)[, 1]),
    stress_group = ifelse(ms <= -1, "low", ifelse(ms < 1, "medium", ifelse(ms >= 1, "high", NA))))

cnames <- colnames(dlong)[!colnames(dlong) %in% c("id", "t")]
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
    mutate(across(all_of(c("gage", "psas", 
                           "stai", "epds", "pss")), 
                  function(x) scale(x)[, 1])) %>%
    mutate(
      ms = scale(epds + stai + pss)[, 1], 
      across(all_of(c(contains("cort"), contains("praq"), contains("ratio"))), function(x) scale(x)[, 1]),
      stress_group = ifelse(ms <= -1, "low", ifelse(ms < 1, "medium", ifelse(ms >= 1, "high", NA))))
  dout
})




# here we look at stress across all time points
stress_indices <- c("ms", "stai", "epds", "pss", "praq_worries_handicap_g", "praq_fear_birth_g",
                    "cortisol_g1", "cortisol_g2", "cortisol_g3", "cortisone_g1",
                    "cortisone_g2", "cortisone_g3", "ratio_g1", "ratio_g2", "ratio_g3",
                    "cortisol_pp", "cortisone_pp", "ratio_pp")


if (!file.exists(here("data/rdata/i_daa_tbl.Rds"))){
  daa_tbl <- map_dfr(stress_indices, function(stress) {
  
    asv_tab <- t(assay(tse_i))
    meta <- colData(tse_i) %>% as.data.frame() %>% 
      select(all_of(stress), dmode, parity, gage, hw,  everything())
    asv_tab <- asv_tab[meta$sid,]
    if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
      fe <- c("t", "dmode", "gage", "parity", "hw", stress)
    } else {
      fe <- c("t", "dmode", "gage", "parity", stress)
    }
    fit_data <- Maaslin2(
      asv_tab,
      meta,
      output = here::here("data/maaslin/1"),
      transform = "AST",
      fixed_effects = fe,
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
  save(daa_tbl, file = here("data/rdata/i_daa_tbl.Rds"))
} else {
  load(here("data/rdata/i_daa_tbl.Rds"))
}

daa_tbl



# do these results differ if we analyze per time point? Note that interactions
# are not possible in Maaslin2 and therefore we must explore this anyways.
tps <- c("t1", "t2", "t3", "t4")
if (!file.exists(here("data/rdata/i_daa_tbl2.Rds"))) {
  daa_tbl_2 <- map_dfr(stress_indices, function(stress) {
    map_dfr(tps, function(time) {
      if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
        fe <- c("dmode", "gage", "parity", "hw", stress)
      } else {
        fe <- c("dmode", "gage", "parity", stress)
      }
      asv_tab <- t(assay(tse_i))
      meta <- colData(tse_i) %>% as.data.frame() %>%
        select(all_of(stress), dmode, parity, gage, hw,  everything(), t) %>%
        filter(t == time)
  
      asv_tab <- asv_tab[meta$sid,]
      fit_data <- Maaslin2(
        asv_tab,
        meta,
        output = here::here("data/maaslin/1"),
        transform = "AST",
        fixed_effects = fe,
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
    select(t, stress, feature, coef, qval) %>%
    mutate(
      feature = str_extract(feature, "g__.*"),
      across(where(is.numeric), round, 3)
    ) %>%
    arrange(desc(abs(coef))) %>%
    rename("Time Point" = t, "Stress Variable" = stress, Feature = feature, Coefficient = coef, q = qval)
  save(daa_tbl_2, file = here("data/rdata/i_daa_tbl2.Rds"))
} else {
  load(here("data/rdata/i_daa_tbl2.Rds"))
}

daa_tbl_2

# lastly we check prenatal stress --> postnatal infant microbiota
load(here("data/rdata/dm.Rds"))
itps <- c("t1", "t2", "t3", "t4")
mtps <- c("t1", "t2", "t3")
stress_indices_m <- c("ms", "stai", "epds", "pss", "praqr_handicap", "praqr_birth",
                    "cortisol", "cortisol", "cortisol", "cortisone", 
                    "cortisone", "cortisone", "ratio", "ratio", "ratio")
stress_indices_i <- c("ms", "stai", "epds", "pss", "praq_worries_handicap_g", "praq_fear_birth_g",
                      "cortisol_g1", "cortisol_g2", "cortisol_g3", "cortisone_g1", 
                    "cortisone_g2", "cortisone_g3", "ratio_g1", "ratio_g2", "ratio_g3")
colData(tse_i) %>% colnames()

if (!file.exists(here("data/rdata/i_daa_tbl3.Rds"))) {
  daa_tbl_3 <- map2_dfr(stress_indices_m, stress_indices_i, function(stress_m, stress_i) {
    map_dfr(tps, function(time) {
      map_dfr(mtps, function(mt) {
        
        if (!(str_detect(stress_i, "cort") | str_detect(stress_i, "ratio"))) {
          stress <- stress_m
        } else {
          stress <- stress_i
        }
        
          if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
            fe <- c("dmode", "gage", "parity", "hw", stress)
          } else {
            fe <- c("dmode", "gage", "parity", stress)
          }
          asv_tab <- t(assay(tse_i))
          meta <- colData(tse_i) %>% as.data.frame() %>%
            filter(t == time)
          if (!(str_detect(stress_i, "cort") | str_detect(stress_i, "ratio"))) {
            meta <- select(meta, -all_of(stress_i))
          }
            
          temp <- filter(dm, t == mt, pre) %>% 
            select(id, all_of(stress_m)) %>%
            mutate(id = as.character(id))
          
          meta <- meta %>%
            rownames_to_column("temp") %>%
            left_join(temp, by = "id") %>%
            column_to_rownames("temp")
          
          asv_tab <- asv_tab[meta$sid,]
          fit_data <- Maaslin2(
            asv_tab,
            meta,
            output = here::here("data/maaslin/1"),
            transform = "AST",
            fixed_effects = fe,
            #random_effects = "id", 
            #reference = c("t", "t1"),  
            normalization = "TSS",
            standardize = FALSE,
            min_prevalence = 0 # prev filterin already done
          )
          filter(fit_data$results, qval <= 0.2, metadata == stress) %>%
            mutate(stress = stress, mt = mt, it = time)
        })
      })
    }) %>%
    
    select(mt, it, stress, feature, coef, qval) %>%
    mutate(
      feature = str_extract(feature, "g__.*"),
      across(where(is.numeric), round, 3),
      mt = ifelse(
        (str_detect(stress, "cort") | str_detect(stress, "ratio")), "32", mt)
    ) %>%
    arrange(desc(abs(coef))) %>%
    rename("Stress Variable" = stress, Feature = feature, Coefficient = coef, q = qval)
  save(daa_tbl_3, file = here("data/rdata/i_daa_tbl3.Rds"))
} else {
  load(here("data/rdata/i_daa_tbl3.Rds"))
}






arrange(daa_tbl, `Stress Variable`, Coefficient, q)
arrange(daa_tbl_2, `Time Point`, `Stress Variable`, Coefficient, q)
arrange(daa_tbl_3, mt, it, `Stress Variable`, Coefficient, q) %>%
  distinct()




# make plot for b fragilis
asv_tab <- t(assay(tse_i))
meta <- colData(tse_i) %>% as.data.frame() %>%
  filter(t == "t2") 

asv_tab <- asv_tab[meta$sid,]
test <- as.data.frame(asv_tab) %>% 
  rownames_to_column("sid") %>%
  select(sid, contains("fragilis"))
colnames(test) <- c("sid", "B.fragilis")
dfp <- full_join(meta, test, by = "sid")
ggplot(dfp, aes(ratio_g1, B.fragilis)) +
  geom_point() +
  geom_smooth()








# check the same at phylum level 
tse_phylum <- agglomerateByRank(tse_i, rank = "Phylum")


# now we look at stress across all time points

stress_indices <- c("ms", "stai", "epds", "pss", "praq_worries_handicap_g", "praq_fear_birth_g",
                    "cortisol_g1", "cortisol_g2", "cortisol_g3", "cortisone_g1",
                    "cortisone_g2", "cortisone_g3", "ratio_g1", "ratio_g2", "ratio_g3",
                    "cortisol_pp", "cortisone_pp", "ratio_pp")


if (!file.exists(here("data/rdata/i_daa_tbl_p.Rds"))){
  daa_tbl <- map_dfr(stress_indices, function(stress) {
    
    asv_tab <- t(assay(tse_phylum))
    meta <- colData(tse_phylum) %>% as.data.frame() %>% 
      select(all_of(stress), dmode, parity, gage, hw,  everything())
    asv_tab <- asv_tab[meta$sid,]
    if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
      fe <- c("t", "dmode", "gage", "parity", "hw", stress)
    } else {
      fe <- c("t", "dmode", "gage", "parity", stress)
    }
    fit_data <- Maaslin2(
      asv_tab,
      meta,
      output = here::here("data/maaslin/1"),
      transform = "AST",
      fixed_effects = fe,
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
      #feature = str_extract(feature, "g__.*"),
      across(where(is.numeric), round, 3)
    ) %>%
    arrange(desc(abs(coef))) %>%
    rename("Stress Variable" = stress, Feature = feature, Coefficient = coef, q = qval)
  save(daa_tbl, file = here("data/rdata/i_daa_tbl_p.Rds"))
} else {
  load(here("data/rdata/i_daa_tbl_p.Rds"))
}

daa_tbl



# do these results differ if we analyze per time point? Note that interactions
# are not possible in Maaslin2 and therefore we must explore this.
tps <- c("t1", "t2", "t3", "t4")
if (!file.exists(here("data/rdata/i_daa_tbl2_p.Rds"))) {
  daa_tbl_2 <- map_dfr(stress_indices, function(stress) {
    map_dfr(tps, function(time) {
      if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
        fe <- c("dmode", "gage", "parity", "hw", stress)
      } else {
        fe <- c("dmode", "gage", "parity", stress)
      }
      asv_tab <- t(assay(tse_phylum))
      meta <- colData(tse_phylum) %>% as.data.frame() %>%
        select(all_of(stress), dmode, parity, gage, hw,  everything(), t) %>%
        filter(t == time)
      
      asv_tab <- asv_tab[meta$sid,]
      fit_data <- Maaslin2(
        asv_tab,
        meta,
        output = here::here("data/maaslin/1"),
        transform = "AST",
        fixed_effects = fe,
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
    select(t, stress, feature, coef, qval) %>%
    mutate(
      #feature = str_extract(feature, "g__.*"),
      across(where(is.numeric), round, 3)
    ) %>%
    arrange(desc(abs(coef))) %>%
    rename("Time Point" = t, "Stress Variable" = stress, Feature = feature, Coefficient = coef, q = qval)
  save(daa_tbl_2, file = here("data/rdata/i_daa_tbl2_p.Rds"))
} else {
  load(here("data/rdata/i_daa_tbl2_p.Rds"))
}

daa_tbl_2

# lastly we check prenatal stress --> postnatal infant microbiota

load(here("data/rdata/dm.Rds"))
itps <- c("t1", "t2", "t3", "t4")
mtps <- c("t1", "t2", "t3")
stress_indices_m <- c("ms", "stai", "epds", "pss", "praqr_handicap", "praqr_birth",
                      "cortisol", "cortisol", "cortisol", "cortisone", 
                      "cortisone", "cortisone", "ratio", "ratio", "ratio")
stress_indices_i <- c("ms", "stai", "epds", "pss", "praq_worries_handicap_g", "praq_fear_birth_g",
                      "cortisol_g1", "cortisol_g2", "cortisol_g3", "cortisone_g1", 
                      "cortisone_g2", "cortisone_g3", "ratio_g1", "ratio_g2", "ratio_g3")
colData(tse_i) %>% colnames()

if (!file.exists(here("data/rdata/i_daa_tbl3_p.Rds"))) {
  daa_tbl_3 <- map2_dfr(stress_indices_m, stress_indices_i, function(stress_m, stress_i) {
    map_dfr(tps, function(time) {
      map_dfr(mtps, function(mt) {
        
        if (!(str_detect(stress_i, "cort") | str_detect(stress_i, "ratio"))) {
          stress <- stress_m
        } else {
          stress <- stress_i
        }
        
        if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
          fe <- c("dmode", "gage", "parity", "hw", stress)
        } else {
          fe <- c("dmode", "gage", "parity", stress)
        }
        asv_tab <- t(assay(tse_phylum))
        meta <- colData(tse_phylum) %>% as.data.frame() %>%
          filter(t == time)
        if (!(str_detect(stress_i, "cort") | str_detect(stress_i, "ratio"))) {
          meta <- select(meta, -all_of(stress_i))
        }
        
        temp <- filter(dm, t == mt, pre) %>% 
          select(id, all_of(stress_m)) %>%
          mutate(id = as.character(id))
        
        meta <- meta %>%
          rownames_to_column("temp") %>%
          left_join(temp, by = "id") %>%
          column_to_rownames("temp")
        
        asv_tab <- asv_tab[meta$sid,]
        fit_data <- Maaslin2(
          asv_tab,
          meta,
          output = here::here("data/maaslin/1"),
          transform = "AST",
          fixed_effects = fe,
          #random_effects = "id", 
          #reference = c("t", "t1"),  
          normalization = "TSS",
          standardize = FALSE,
          min_prevalence = 0 # prev filterin already done
        )
        filter(fit_data$results, qval <= 0.2, metadata == stress) %>%
          mutate(stress = stress, mt = mt, it = time)
      })
    })
  }) %>%
    
    select(mt, it, stress, feature, coef, qval) %>%
    mutate(
      #feature = str_extract(feature, "g__.*"),
      across(where(is.numeric), round, 3),
      mt = ifelse(
        (str_detect(stress, "cort") | str_detect(stress, "ratio")), "32", mt)
    ) %>%
    arrange(desc(abs(coef))) %>%
    rename("Stress Variable" = stress, Feature = feature, Coefficient = coef, q = qval)
  save(daa_tbl_3, file = here("data/rdata/i_daa_tbl3_p.Rds"))
} else {
  load(here("data/rdata/i_daa_tbl3_p.Rds"))
}


filter(daa_tbl, q <= 0.2) %>% arrange(`Stress Variable`, Coefficient, q)
filter(daa_tbl_2, q <= 0.2) %>% arrange(`Time Point`, `Stress Variable`, Coefficient, q)
filter(daa_tbl_3, q <= 0.2) %>% arrange(mt, it, `Stress Variable`, Coefficient, q) %>%
  distinct()

