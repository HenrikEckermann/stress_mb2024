set.seed(1)
library(tidyverse)
library(readxl)
library(glue)
library(mice)


path <- here::here("data")
files <- list.files(
  path, 
  full.names = TRUE, 
  include.dirs = FALSE, 
  recursive = TRUE,
  pattern = ".xlsx$"
  )

files

chosen_files <- c(
  "/SMILEY_extracovariates_data.docx.xlsx",
  "/SMILEY_6wks_hair_cortisol_RawClean.xlsx",
  "/SMILEY_32wks_hair_cortisol_RawClean.xlsx",
  "/SMILEY_12wks_Feedings_RawClean.xlsx",
  "/SMILEY_12wks_Health-mother-baby_RawClean.xlsx",
  "/SMILEY_12wks_MSSS_Clean.xlsx",
  "/SMILEY_12wks_PSAS-RSF-C_Clean.xlsx",
  "/SMILEY_12wks_Stress.xlsx",
  "/SMILEY_2wks_Feedings_RawClean.xlsx",
  "/SMILEY_2wks_Health-mother-baby_RawClean.xlsx",
  "/SMILEY_2wks_PSAS-RSF-C_Clean.xlsx",
  "/SMILEY_2wks_stress.xlsx",
  "/SMILEY_6wks_Feedings_RawClean.xlsx",
  "/SMILEY_6wks_Health-mother-baby_RawClean.xlsx",
  "/SMILEY_6wks_PSAS_Clean.xlsx",
  "/SMILEY_6wks_Stress.xlsx",
  "/SMILEY_32wks_Baby-PAWS_Clean.xlsx",
  "/SMILEY_32wks_MSSS_Clean.xlsx",
  "/SMILEY_32wks_PES-brief_RawClean.xlsx",
  "/SMILEY_32wks_PRAQ-R2_Clean.xlsx",
  "/SMILEY_32wks_Stress.xlsx",
  "/SMILEY_32wks_TPDS_Clean.xlsx",
  "/Smiley_12wk_Stool_RawClean.xlsx",
  "/SMILEY_2wks_Stool_RawClean.xlsx",
  "/SMILEY_6wk_Stool_RawClean.xlsx",
  "/SMILEY_child_age.xlsx",
  "/SMILEY_partner_parity.xlsx",
  "/SMILEY_8mo_InfantStool_Rawdata.xlsx",
  "/SMILEY_8mo_Feedings_RawClean"
)





# all data will be stored in a list of dataframes
all_d <- map2(chosen_files, seq_along(chosen_files), function(file, i) {
  fp <- files[str_detect(files, file)]
  fp
  if (i %in% c(23, 24, 25, 26)) {
    tmp1 <- read_excel(fp, na = "99999")
    tmp <- read_excel(fp, col_types = c("guess", "date", rep("guess", dim(tmp1)[2] - 2)), na = "99999")
  } else {
    tmp <- read_excel(fp, na = "99999")
  }
  colnames(tmp) <- str_to_lower(colnames(tmp))
  tmp
})




# does data have same number of participants as expected?
map(all_d, ~dim(.x))
# what kind of data do we have?
map(all_d, ~colnames(.x))
map(all_d, ~head(.x))

# for each dataframe, I will make a selection of variables
tmp1 <- select(
  all_d[[1]],
  id,
  sex = baby_sex,
  birthweight = baby_birthweight,
  dmode = delivery_mode, 
  gage = gestational_age_birth_weeks,
  edu = demographics_3,
  hh_n = demographics_7,
  hh_who = demographic_8,
  pet = demographics_9,
  pet_which = demographics_9a,
  bmi_g = bmi_p32,
  bmi_pp = bmi_6wks 
)



tmp2 <- select(
  all_d[[2]],
  id,
  cortisol_pp = smiley_6wk_cortisol1,
  cortisone_pp = smiley_6wk_cortisone1) %>%
  mutate(
    cortisol_pp = as.numeric(cortisol_pp),
    cortisone_pp = as.numeric(cortisone_pp),
    ratio_pp = cortisol_pp/cortisone_pp
    ) %>%
  mutate(across(contains("_pp"), function(x) log(x)))

# get the hair covs:
file <- here("data/SMILEY_hair_6wks/SMILEY_6wks_hair_questionnaire_RawClean.xlsx")
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
         hw = ifelse(is.na(hw), 0, hw),
         id = as.character(id)
  ) %>%
  select(id, hw)

tmp2 <- left_join(tmp2, cort2, by = "id")

tmp3 <- select(
  all_d[[3]],
  id, 
  cortisol_g1 = "smiley_32wk_cortisol1\r\n",
  cortisol_g2 = "smiley_32wk_cortisol2\r\n",  
  cortisol_g3 = "smiley_32wk_cortisol3\r\n",  
  cortisone_g1 = "smiley_32wk_cortisone1\r\n", 
  cortisone_g2 = "smiley_32wk_cortisone2\r\n",
  cortisone_g3 = "smiley_32wk_cortisone3\r\n") %>%
  mutate(
    across(contains("cort"), function(x) as.numeric(x)), 
    ratio_g1 = cortisol_g1 / cortisone_g1,
    ratio_g2 = cortisol_g2 / cortisone_g2,
    ratio_g3 = cortisol_g3 / cortisone_g3
  ) %>%
  mutate(across(contains("_g"), function(x) log(x)))
  

# add hw
file <- here("data/SMILEY_hair_preg32wks/SMILEY_32wks_hair_questionnaire_RawClean.xlsx")

cort2 <- read_excel(file) %>% 
  select(id = ID, lessthanone = "4a", oneortwo = "4b", threeorfour = "4c", morethanfour = "4d") %>%
  pivot_longer(-id, names_to = "frequency", values_to = "values") %>%
  mutate(values = ifelse(values == 99999, NA, ifelse(values == 88888, NA, values))) %>%
  # group_by(frequency) %>%
  # summarise(sum = sum(values, na.rm = TRUE))
  pivot_wider(names_from = "frequency", values_from = "values") %>%
  mutate(hw_g = ifelse(lessthanone, 0, ifelse(oneortwo, 0, ifelse(threeorfour, 1,  ifelse(morethanfour, 1, NA)))),
         hw_g = ifelse(is.na(hw_g), 0, hw_g),
         id = as.character(id)
  ) %>%
  select(id, hw_g)

tmp3 <- left_join(tmp3, cort2, by = "id")

tmp4 <- select(
  all_d[[4]],
  id,
  # solid left out as nobody gave solids
  # solid_wk12 = feeding_solidfoods_7to12wks,
  feeding_wk12 = feeding_general_7to12wks,
)


tmp5 <- select(
  all_d[[5]],
  id,
  b_abx_wk12 = antibiotics_baby_6to12wks,
  m_abx_wk12 = antibiotics_mother_6to12wks,
  novac = "vaccin_12wks#mijn_kind_heeft_nog_geen_vaccinaties_gekregen"
)

tmp6 <- select(
  all_d[[6]],
  id,
  msss_wk12 = sum
  )

tmp7 <- select(
  all_d[[7]],
  id, 
  psas_wk12 = sum
)

tmp8 <- select(
  all_d[[8]],
  id,
  epds_wk12 = wks12_epds,
  pss_wk12 = `wks12_pss-10`,
  stai_wk12 = wks12_stai
)

tmp9 <- select(
  all_d[[9]],
  id,
  feeding_wk2 = feeding_general_birthtill2wks
)


tmp10 <- select(
  all_d[[10]],
  id,
  b_abx_wk2 = antibiotics_baby_2wks,
  m_abx_wk2 = antibiotic_mother_2wks
)



tmp11 <- select(
  all_d[[11]],
  id,
  psas_wk2 = sum
)


tmp12 <- select(
  all_d[[12]],
  id,
  epds_wk2 = wks2_epds,
  pss_wk2 = `wks2_pss-10`,
  stai_wk2 = wks2_stai
)


tmp13 <- select(
  all_d[[13]],
  id,
  feeding_wk6 = feeding_general_3to6wks_1
)

tmp14 <- select(
  all_d[[14]],
  id,
  b_abx_wk6 = antibiotics_baby_2to6wks,
  m_abx_wk6 = antibiotics_mother_2to6wks
)


tmp15 <- select(
  all_d[[15]],
  id,
  psas_wk6 = sum
)


tmp16 <- select(
  all_d[[16]],
  id,
  epds_wk6 = wks6_epds,
  pss_wk6 = `wks6_pss-10`,
  stai_wk6 = wks6_stai
)


tmp17 <- select(
  all_d[[17]],
  id,
  paws_g = sum_val
)


tmp18 <- select(
  all_d[[18]],
  id,
  msss_g = sum
)




tmp20 <- select(
  all_d[[20]],
  id,
  praq_fear_birth_g = fear_birth,
  praq_worries_handicap_g = worries_handicap
)


tmp21 <- select(
  all_d[[21]],
  id,
  epds_g = wks32_epds,
  pss_g = `wks32_pss-10`,
  stai_g = wks32_stai
)


tmp22 <- select(
  all_d[[22]],
  id,
  tpds_g = sum
)


tmp23 <- select(
  all_d[[23]],
  id,
  stool_date_wk12 = stool_1,
  stool_time_wk12 = stool_2
)


tmp24 <- select(
  all_d[[24]],
  id,
  stool_date_wk2 = stool_1,
  stool_time_wk2 = stool_2
)


tmp25 <- select(
  all_d[[25]],
  id,
  stool_date_wk6 = stool_1
  # stool_time_wk6 = stool_2
)

tmp26 <- select(
  all_d[[26]],
  id,
  birthdate = "child_date_of_birth"
)

tmp27 <- select(
  all_d[[27]],
  id,
  partner = partner_yesno,
  parity = firstchild
)

tmp28 <- select(
  all_d[[28]],
  id,
  t1 = stool_7,
  t2 = stool_8) %>%
  mutate(
    t1 = ifelse(t1 == 2, 0, ifelse(t1 == 1, 1, ifelse(is.na(t1), 0, 0))),
    t2 = ifelse(t2 == 2, 0, ifelse(t2 == 1, 1, ifelse(is.na(t2), 0, 0))),
    b_abx_wk32 = ifelse(as.numeric(t1) + as.numeric(t2) > 0, 1, 0),
    b_abx_wk32 = ifelse(is.na(b_abx_wk32), 0, b_abx_wk32)
  ) %>%
  filter(id != "EXAMPLE") %>%
  select(id, b_abx_wk32)


tmp29 <- select(
  all_d[[29]],
  id, 
  t = feeding_general_fu) %>%
  mutate(feeding_wk32 = ifelse(t == 1, 2, ifelse(t == 2, 2, ifelse(t == 3, 3, ifelse(t == 4, 1, NA))))) %>%
  select(id, feeding_wk32)

# get stress from mothers at 8 months pp
tmp30 <- read_csv(here("data/ms_raw.csv")) %>%
  filter(!pre, week == 32) %>%
  select(
    id, 
    stai_wk32 = stai,
    epds_wk32 = epds,
    pss_wk32 = pss
)





all_d <- map(list(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10,
                  tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp18,
                  tmp20, tmp21, tmp22, tmp23, tmp24, tmp25, tmp26, tmp27,
                  tmp28, tmp29, tmp30), function(d) {
                    d$id <- as.integer(d$id)
                    d
                  })

excl <- c("id", "hh_who", "pet_which", "stool_date_wk12", "stool_date_wk6", "stool_date_wk2", "stool_time_wk2", "stool_time_wk6", "stool_time_wk12", "birthdate")
d <- all_d %>% reduce(full_join, by = "id")
incl <- colnames(d)[!colnames(d) %in% excl]
d <- all_d %>% reduce(full_join, by = "id") %>%
  mutate(
    across(all_of(incl), function(x) as.numeric(x)),
    birthdate = as.POSIXct(birthdate), 
    stool_date_wk2 = as.POSIXct(stool_date_wk2),
    age_wk2 = as.numeric(stool_date_wk2 - birthdate),
    age_wk6 = as.numeric(stool_date_wk6 - birthdate),
    age_wk12 = as.numeric(stool_date_wk12 - birthdate),
    across(-birthdate, function(x) ifelse(x == 99999, NA, ifelse(x == 88888, NA, x))),
    edu = as.numeric(edu), sex = as.factor(sex), dmode = as.factor(dmode), 
    pet = as.factor(pet), novac = as.factor(novac), 
    edu = factor(ifelse(edu <=5, "low", ifelse(edu > 5, "high", edu))),
    hw = as.factor(hw), hw_g = as.factor(hw_g),
    parity = as.factor(parity)
    )

# ask Hellen if ever anyone coded the household member variable to siblings,
# parity etc.:
# hellen said we must generate our own siblings variable but parity will be added
# first I check the unique terms used to describe siblings
str_split(d$hh_who, ",") %>% unlist() %>% 
  str_split(" ") %>% unlist() %>%
  unique()

sibling_terms <- c(
  "zoon", 
  "Kind", 
  "kind", 
  "jongen",
  "meisje",
  "Dochter",
  "dochter",
  "Stiefdochter",
  "stiefzoon",
  "stiefdochter",
  "stiefkind",
  "pleegkind",
  "Bonuskind",
  "zus",
  "Broer",
  "Zoon"
)
term <- paste(sibling_terms, collapse = "|")
d <- mutate(
  d, 
  siblings = ifelse(str_detect(hh_who, term), 1, 0),
  # note that Hellen informed me that for id 589 "broer" means the brother
  # of the actual participant, therefore I change here manually
  siblings = ifelse(id == 589, 0, siblings),
  # also for some of the missing values of siblings, we can borrow from
  # firstchild as whoever generated that variable knows more than me
  siblings = ifelse(is.na(siblings), as.numeric(!parity), siblings),
  siblings = as.factor(siblings)
  )

# there should be no overlap between siblings and firstchild, except for 
# stepchildren, which I count as siblings for the microbiota project
filter(d, siblings == parity) %>%
  select(id, siblings, parity, hh_n, hh_who)


d <- select(d, -c(hh_n, hh_who)) 
str(d)
dlong <- select(d, -birthdate, -contains("stool_date_wk"), -contains("stool_time_wk")) %>%
  pivot_longer(
    matches("_wk\\d+"),
    names_to = c("variable", "week"),
    names_pattern = "(.+)_wk(\\d+)") %>%
  pivot_wider(
    names_from = variable,
    values_from = value) %>%
  mutate(across(contains("feeding"), function(x) as.factor(x)),
         across(contains("abx"), function(x) as.factor(x)))


d <- mutate(d, across(contains("feeding"), function(x) as.factor(x)),
            across(contains("abx"), function(x) as.factor(x)))
# explore the data to decide what approaches can/should bet taken

# how many have low edu?
count(d, edu)

# how many mothers took abx
count(dlong, week, m_abx)

# how many babies took abx
count(dlong, week, b_abx)

# how many have siblings
count(d, siblings)

# partners 
count(d, sex)

# there might be issues with the date columns because we have some weird ages
filter(dlong, age <7 | age > (12*7) + 50) %>%
  select(id, week, age)

# 507, 580, 618,  had impossible date combinations without knowing which date is incorrect
# fixed year also for 529, 253, 578, 586, 589, 593, 612, 613

filter(dlong, id %in% c(507, 580, 618)) %>%
  select(id, week, age)
# as of now we must filter 507 w2, 580 w2, 580 w6 and 618w12
dlong$age[dlong$week == 6 & dlong$id == 507] <- NA
dlong$age[dlong$week == 2 & dlong$id == 580] <- NA
dlong$age[dlong$week == 6 & dlong$id == 580] <- NA
dlong$age[dlong$week == 12 & dlong$id == 618] <- NA

# how many kids were not vaccinate at week 12?
count(d, novac)


# change age var also in the wide df
d$age_wk6[d$id == 507] <- NA
d$age_wk2[d$id == 580] <- NA
d$age_wk6[d$id == 580] <- NA
d$age_wk12[d$id == 618] <- NA
colnames(d)

# now that we have the df, create imps for all vars except cortisol
diimp <- mice(select(d, -pet_which, -contains("date")), m = 5)
dlongiimp <- mice(select(dlong, -pet_which, -contains("date")), m = 5)
diimp <- map(1:5, function(m) {
  dimp <- complete(diimp, action = 5) %>% 
    select(-contains("cort"), -contains("ratio"))
  dtemp <- select(d, id, contains("cort"), contains("ratio"))
  dimp <- left_join(dimp, dtemp, by = "id")
})
colnames(dlong)
dlongiimp <- map(1:5, function(m) {
  dimp <- complete(dlongiimp, action = 5) %>% 
    select(-contains("cort"), -contains("ratio")) %>%
    mutate(
      stai = scale(stai)[, 1],
      epds = scale(epds)[, 1],
      pss = scale(pss)[, 1],
      ms = scale(stai + epds + pss)[, 1],
    )
  dtemp <- select(dlong, id, week, contains("cort"), contains("ratio"))
  dimp <- left_join(dimp, dtemp, by = c("id", "week"))
})

dlong <-     mutate(
  dlong,
  stai = scale(stai)[, 1],
  epds = scale(epds)[, 1],
  pss = scale(pss)[, 1],
  ms = scale(stai + epds + pss)[, 1],
)




# store new dfs
write_csv(dlong, file = here::here("data/dlong.csv"))
write_csv(
    select(d, -c(
      "stool_date_wk12", 
      "stool_date_wk6", 
      "stool_date_wk2", 
      "stool_time_wk2", 
      # "stool_time_wk6", 
      "stool_time_wk12", 
      "birthdate")),
    file = here::here("data/dwide.csv")
  )


save(d, dlong, diimp, dlongiimp, file = here("data/rdata/di.Rds"))

