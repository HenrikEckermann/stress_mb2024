library(mia)
library(tidyverse)
library(here)
library(glue)
library(table1)


load(file = here::here("data/rdata/tse.Rds"))
load(here("data/rdata/dm.Rds"))
load(here("data/rdata/di.Rds"))


ga <- readxl::read_excel(here("data/SMILEY_gestational_age.xlsx"), na = "99999") %>%
  select(id = ID, ga_18 = "Gestational_age_prenatal_18wks_weeks", ga_32 = "Gestational_age_prenatal_32wks_weeks") %>%
  mutate(id = as.integer(id))
iage <- readxl::read_excel("data/infant_ages.xlsx", na = "99999") %>%
  select(id = ID, ia_8m = "Infant_age_8MO_days") %>%
  mutate(id = as.integer(id))

d <- left_join(d, ga, by = "id") %>%
  left_join(iage, by = "id")
colnames(d)


# table 1

# get ids for all with at least one stool sample
stool_ids <- colData(tse_s) %>%
  as.data.frame() %>%
  select(id) %>%
  distinct() %>%
  .$id
stool_ids %>% length() 

assay(tse_s)[, str_detect(colnames(assay(tse_s)), "501")]
stress_non_scaled <- read_csv(here("data/ms_raw.csv"))
cnames <- select(stress_non_scaled, -id, -week, -pre) %>% colnames()

dm <- select(dm, -all_of(cnames)) %>%
  left_join(stress_non_scaled, by = c("id", "week", "pre"))

count(dm, week, pre, abx)
colnames(dm)
factor(dm$abx, levels = c("0", "1"), labels = c("No", "Yes"))
dmwide <- mutate(dm, week = ifelse(pre, glue("{week}_pre"), glue("{week}_post"))) %>%
  select(id, week, activity, pbmi, age, parity, stai, epds, pss,
       contains("praq"), dhd15_total, edu, csection, abx) %>%
  mutate(
    parity = factor(parity, levels = c("0", "1"), labels = c("No", "Yes")),
    csection = factor(csection, levels = c("0", "1"), labels = c("vaginal", "c-section")),
    abx = factor(abx, levels = c("0", "1"), labels = c("No", "Yes")),
    edu = factor(edu, levels = c("0", "1"), labels = c("low", "high"))
  ) %>%
  pivot_wider(id_cols = c(id, pbmi, age, parity, edu, csection), names_from = week, values_from = c(activity, stai, epds, pss, praqr_handicap, praqr_birth, dhd15_total, abx)) %>%
  filter(id %in% stool_ids)





dinf <- select(d, id, sex, birthweight, gage, pet, contains("cort"), contains("ratio"), 
       contains("abx"), contains("age_wk"), contains("feeding_wk"), ga_18, ga_32, ia_8m) %>%
  filter(id %in% stool_ids)
min(dinf$age_wk12, na.rm = TRUE)
dinf$age_wk12
colnames(psas)
psas <- read_csv(here("data/psas.csv")) %>%
  select(id = ID, week, psas) %>%
  pivot_wider(names_from = week, values_from = psas, names_prefix = "psas_") %>%
  filter(id %in% stool_ids)

joined <- full_join(dmwide, dinf, by = "id") %>%
  full_join(psas, by = "id") %>%
  select(
    -"praqr_handicap_32_post",
    -"praqr_birth_32_post",
    -"dhd15_total_32_post",
    -contains("m_abx")
  )


f <- select_if(joined, is.factor) %>% 
  select(-parity, -csection, -edu, -contains("abx")) %>%
  colnames()

for (i in seq_along(f)) {
  if (f[i] == "sex") {
    joined[[f[i]]] <- factor(
      joined[[f[i]]], 
      levels = c("1", "2"),
      labels = c("Boy", "Girl")
    )
  } else if (str_detect(f[i], "feeding")) {
    joined[[f[i]]] <- factor(
      joined[[f[i]]], 
      levels = c("1", "2", "3"),
      labels = c("Breastmilk", "Formula", "Mixed")
    )
  } else if (f[i] == "csection") {
    joined[[f[i]]] <- factor(
      joined[[f[i]]], 
      levels = c("0", "1"),
      labels = c("vaginal", "c-section")
    )
  } else {
    joined[[f[i]]] <- factor(
      joined[[f[i]]], 
      levels = c("0", "1"),
      labels = c("No", "Yes")
    )
  }
}

colnames(joined)
count(joined, abx_18_pre)

mid <- select(joined, id, sex) %>%
  filter(is.na(sex)) %>%
  .$id
mid %in% stool_ids
stool_ids
# change labels
l <- c(
  "id",
  
  "Pre-pregnancy BMI",
  "Maternal age (in years)",
  "Firstborn",
  "Maternal education",
  "Delivery mode",
  
  
  "Physical activity gestational week 18",
  "Physical activity postpartum week 32",
  "Physical activity gestational week 32",
  
  "STAI gestational week 18", 
  "STAI postpartum week 32", 
  "STAI gestational week 32", 
  
  "EPDS gestational week 18", 
  "EPDS postpartum week 32", 
  "EPDS gestational week 32", 
  
  "PSS-10 gestational week 18", 
  "PSS-10 postpartum week 32", 
  "PSS-10 gestational week 32", 
  
  "PRAQR2-H gestational week 18",
  "PRAQR2-H gestational week 32",
  "PRAQR2-Bgestational week 18",
  "PRAQR2-B gestational week 32",
  
  "Maternal Dutch Health Diet Index gestational week 18",
  "Maternal Dutch Health Diet Index gestational week 32",
  
  "Antibiotics gestational week 18",
  "Antibiotics gestational week 32",
  "Antibiotics postpartum week 32",
  
  "Infant sex",
  "Birthweight",
  "Gestational age",
  "Pet",
  
  "Cortisol (log) postpartum",
  "Cortisone (log) postpartum",
  "Cortisol (log) gestational week 23-32",
  "Cortisol (log) gestational week 15-23",
  "Cortisol (log) gestational week 6-15",
  
  "Cortisone (log) gestational week 23-32",
  "Cortisone (log) gestational week 15-23",
  "Cortisone (log) gestational week 6-15",
  
  "Log-ratio Cortisol/Cortisone postpartum",
  "Log-ratio Cortisol/Cortisone gestational week 23-32",
  "Log-ratio Cortisol/Cortisone gestational week 15-23",
  "Log-ratio Cortisol/Cortisone gestational week 6-15",
  
  "Infant antibiotics week 12",
  "Infant antibiotics week 2",
  "Infant antibiotics week 6",
  "Infant antibiotics week 32",
  
  "Infant age week 2",
  "Infant age week 6",
  "Infant age week 12",
  
  "Infant feeding week 12",
  "Infant feeding week 2",
  "Infant feeding week 6",
  "Infant feeding week 32",
  
  "Maternal gestational age gestational week 18",
  "Maternal gestational age gestational week 32",
  "Infant age week 32",
  
  "PSAS (short) postpartum week 2",
  "PSAS (long) postpartum week 6",
  "PSAS (short) postpartum week 12"
)

for (i in seq_along(colnames(joined))) {
  label(joined[[colnames(joined)[i]]]) <- l[i]
}

count(joined, abx_32_post)

tbl1 <- table1(~ pbmi + age + parity + edu + csection +
                 activity_18_pre + activity_32_pre +
                 cortisol_g1 + cortisol_g2 +
                 cortisol_g3 + cortisol_pp + cortisone_g1 + cortisone_g2 +
                 cortisone_g3 + cortisone_pp + ratio_g1 + ratio_g2 +
                 ratio_g3 + ratio_pp + dhd15_total_18_pre + dhd15_total_32_pre +
                 
                 abx_18_pre + abx_32_post + abx_32_pre +
                 
                 gage + epds_18_pre + epds_32_pre + epds_32_post + praqr_birth_18_pre +
                 praqr_birth_32_pre + praqr_handicap_18_pre + praqr_handicap_32_pre +
                 psas_2 + psas_6 + psas_12 +
                 pss_18_pre + pss_32_pre + pss_32_post + stai_18_pre + stai_32_pre +
                 stai_32_post +
                 
                 
                 age_wk2 + age_wk6 + age_wk12 + b_abx_wk2 + b_abx_wk6 + b_abx_wk12 +
                 b_abx_wk32 + birthweight + feeding_wk2 + feeding_wk6 + feeding_wk12 +
                 feeding_wk32 + ga_18 + ga_32 + ia_8m + pet + sex, 
               
               data = joined,
               caption = "Demographics Table",
               transpose = FALSE,
               overall = ""
)


tbl1

colData(tse_s) %>%
  as.data.frame() %>%
  group_by(origin, pre, week) %>%
  summarise(n = n())

# table 2
colData(tse_s) %>%
  as.data.frame() %>%
  filter(week != 104) %>%
  group_by(origin, id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(origin, n) %>%
  summarise(n())


tbl2 <- colData(tse_s) %>%
  as.data.frame() %>%
  group_by(origin, pre, week) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(
    origin = ifelse(origin == "i", "Infant", "Mother"),
    Stage = ifelse(pre, "Prenatal", "Postnatal"),
    Abbreviation = ifelse(week == 2, "I1", ifelse(week == 6, "I2", ifelse(
      week == 12, "I3", ifelse((week == 32 & origin == "Infant"), "I4", ifelse(
        (week == 32 & Stage == "Prenatal"), "M2", ifelse(
          (week == 32 & Stage == "Postnatal" & origin == "Mother"), "M3", ifelse(
            week == 18, "M1", NA))))))),
    Abbreviation = factor(Abbreviation, levels = c(
      "M1", "M2", "M3", "I1", "I2", "I3", "I4"
    ))
  ) %>%
  select(Origin = origin, Stage, Week = week, n, Abbreviation) %>%
  filter(!is.na(Abbreviation)) %>%

  arrange(Abbreviation)

tbl2



hc_vars <- c("cortisol_g3", "cortisol_g2", "cortisol_g1", "cortisol_pp")
ns <- map_dbl(hc_vars, function(hc) {
  rows <- select(d, all_of(hc)) %>%
    na.omit() %>%
    dim()
  rows[1]
})

tbl3 <- tibble(
  Stage = c("Prenatal", "Prenatal", "Prenatal", "Postpartum"),
  "Collection week" = c("32", "32", "32", "8"),
  "Actual week (mean \u00b1 SD)" = c(
    "32.28 \u00b1 1.52",
    "32.28 \u00b1 1.52",
    "32.28 \u00b1 1.52",
    "8.19 \u00b1 1.17"
    ),
  "Measured period (weeks)" = c(
    "6 - 15",
    "15 - 23",
    "23 - 32",
    "4 - 8"
  ),
  n = ns,
  "Cortisol" = c("HCS1", "HCS2", "HCS3",  "HCS4"),
  "Cortisone" = c("HCN1", "HCN2", "HCN3", "HCN4"),
  "Log-ratio" = c("HCR1", "HCR2", "HCR3", "HCR4")
)
tbl3

save(tbl1, tbl2, tbl3, file = here("data/rdata/tables.Rds"))

# correlation table

mtimes <- list(
  c(18, 1),
  c(32, 1),
  c(32, 0)
)
not_all_na <- function(x) any(!is.na(x))
cort <- select(
  d,
  id,
  "HCS3" = cortisol_g1,
  "HCS2" = cortisol_g2,
  "HCS1" = cortisol_g3,
  "HCN3" = cortisone_g1,
  "HCN2" = cortisone_g2,
  "HCN1" = cortisone_g3,
  "HCS4" = cortisol_pp,
  "HCN4" = cortisone_pp
  )

cors <- map(mtimes, function(vec) {
  weekt <- vec[1]
  pret <- vec[2]
  
  temp <- select(
    dm, id, pre,  week,  STAI = stai, EPDS = epds, 
    `PSS-10` = pss, PRAQ1 = praqr_handicap, PRAQ2 = praqr_birth,
    MS = ms) %>%
    filter(week == weekt, pre == pret) %>%
    left_join(cort, by = "id") %>%
    select(-id, -week, -pre)
  
  cor(temp, use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    select(where(not_all_na)) %>%
    na.omit() %>%
    rownames_to_column("var2") %>%
    pivot_longer(-var2, names_to = "var1", values_to = "r") 
  
  # ggplot(aes(var2, var1, fill = r, label = round(r, 2))) +
  #   geom_tile() +
  #   geom_text() +
  #   scale_fill_gradient2(mid="#FBFEF9", low="#0C6291", high="#A63446", limits=c(-1,1)) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
})
cors[[1]]

cors2 <- map2(mtimes, cors, function(vec, dtemp) {
  weekt <- vec[1]
  pret <- vec[2]
  
  temp <- select(
    dm, id, pre,  week,  STAI = stai, EPDS = epds, 
    `PSS-10` = pss, PRAQ1 = praqr_handicap, PRAQ2 = praqr_birth,
    MS = ms) %>%
    filter(week == weekt, pre == pret) %>%
    left_join(cort, by = "id") %>%
    select(-id, -week, -pre)
  
  df_r <- map2_dfr(dtemp$var1, dtemp$var2, function(var1, var2) {
      test <- cor.test(temp[[var1]], temp[[var2]])
      tibble(var1 = var1, var2 = var2, r = test$estimate, p = test$p.value)
    })
  
  mutate(df_r, r_if_sig = ifelse(p < 0.05, r, NA)) %>%
    mutate(
      var1 = ifelse(var1 == "PRAQ1", "PRAQR2-H", ifelse(var1 == "PRAQ2", "PRAQR2-B", var1)),
      var2 = ifelse(var2 == "PRAQ1", "PRAQR2-H", ifelse(var2 == "PRAQ2", "PRAQR2-B", var2))
      ) %>%
    
  ggplot(aes(var2, var1, fill = r, label = round(r_if_sig, 2))) +
    geom_tile() +
    geom_text() +
    scale_fill_gradient2(mid="#FBFEF9", low="#0C6291", high="#A63446", limits=c(-1,1)) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
})

cors2[[3]]


itimes <- c("2", "6", "12", "32")

icors <- map(itimes, function(it) {
  
  temp <- select(
    dlong, id, week,  STAI = stai, EPDS = epds, 
    `PSS-10` = pss, PRAQ1 = praq_worries_handicap_g, PRAQ2 = praq_fear_birth_g,
    MS = ms) %>%
    filter(week == it) %>%
    left_join(cort, by = "id") %>%
    select(-id, -week)
  
  cor(temp, use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    select(where(not_all_na)) %>%
    na.omit() %>%
    rownames_to_column("var2") %>%
    pivot_longer(-var2, names_to = "var1", values_to = "r") 
  
  # ggplot(aes(var2, var1, fill = r, label = round(r, 2))) +
  #   geom_tile() +
  #   geom_text() +
  #   scale_fill_gradient2(mid="#FBFEF9", low="#0C6291", high="#A63446", limits=c(-1,1)) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
})


icors2 <- map2(itimes, icors, function(it, dtemp) {

  temp <- select(
    dlong, id, week,  STAI = stai, EPDS = epds, 
    `PSS-10` = pss, PRAQ1 = praq_worries_handicap_g, PRAQ2 = praq_fear_birth_g,
    MS = ms) %>%
    filter(week == it) %>%
    left_join(cort, by = "id") %>%
    select(-id, -week)
  
  df_r <- map2_dfr(dtemp$var1, dtemp$var2, function(var1, var2) {
    test <- cor.test(temp[[var1]], temp[[var2]])
    tibble(var1 = var1, var2 = var2, r = test$estimate, p = test$p.value)
  })
  
  mutate(df_r, r_if_sig = ifelse(p < 0.05, r, NA)) %>%
    mutate(
      var1 = ifelse(var1 == "PRAQ1", "PRAQR2-H", ifelse(var1 == "PRAQ2", "PRAQR2-B", var1)),
      var2 = ifelse(var2 == "PRAQ1", "PRAQR2-H", ifelse(var2 == "PRAQ2", "PRAQR2-B", var2))
    ) %>%
    
    ggplot(aes(var2, var1, fill = r, label = round(r_if_sig, 2))) +
    geom_tile() +
    geom_text() +
    scale_fill_gradient2(mid="#FBFEF9", low="#0C6291", high="#A63446", limits=c(-1,1)) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
    
})
icors2[[3]]

save(cors2, icors2, file = here("data/rdata/corplots.Rds"))

