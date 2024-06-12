library(mia)
library(tidyverse)
library(vegan)
library(here)
library(glue)
library(scater)
library(patchwork)



load(file = here::here("data/rdata/tse.Rds"))
load(here("data/rdata/dm.Rds"))



# what counts are the highest and lowest
assay(tse_s) %>% as.data.frame() %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  filter(value != 0) %>%
  summarise(min = min(value), max = max(value))

tse_s <- tidySummarizedExperiment::filter(tse_s, week != "104")

tse_g <- agglomerateByRank(tse_s, rank = "Genus")
tse <- transformAssay(tse_s, method = "clr", name = "clr", pseudocount = 0.000001)
tse_g <- transformAssay(tse_g, method = "clr", name = "clr", pseudocount = 0.000001)
# we calculate volatility for mothers and infants separately

# mothers 
tse_m <- tse[, colData(tse)$origin == "m"]


colData(tse_m) %>% rownames()
sum(duplicated(colData(tse_m) %>% rownames()))
sum(duplicated(dm$sid))
dm$id <- as.character(dm$id)
dm$sid <- glue("{dm$id}m{dm$pre}{dm$week}")
unique(colData(tse_m)$t)
colData(tse_m) <- colData(tse_m) %>%
                    as.data.frame() %>%
                    #rownames_to_column("sid") %>%
                    select(-id, -pre, -week) %>%
                    left_join(select(dm, -t), by = "sid") %>%
                    column_to_rownames("sid") %>%
                    DataFrame()



# visualize the data

# Beta diversity metrics like Bray-Curtis are often applied to relabundances
tse_m <- transformAssay(tse_m,
                      assay.type = "counts",
                      method = "relabundance")

# Run PCoA on relabundance assay with Bray-Curtis distances
tse_m <- runMDS(tse_m,
              FUN = vegan::vegdist,
              method = "bray",
              assay.type = "relabundance",
              name = "MDS_bray")

# Create ggplot object
p_bray <- plotReducedDim(
  tse_m, "MDS_bray",
  colour_by = "t",
  point_alpha = 1,
  point_size = 3
  ) +
  scale_color_manual(values = c('#fee0d2','#fc9272','#de2d26')) +
  theme_bw(base_size = 15)

# Calculate explained variance
e <- attr(reducedDim(tse_m, "MDS_bray"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p_bray <- p_bray + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")) 

p_bray



tse_m <- runMDS(tse_m,
              FUN = mia::calculateUnifrac,
              name = "Unifrac",
              tree = rowTree(tse_m),
              ntop = nrow(tse_m),
              assay.type = "counts")

p_uni <- plotReducedDim(
  tse_m, "Unifrac", colour_by = "t",
  point_alpha = 1, 
  point_size = 3
  ) +
  scale_color_manual(values = c('#fee0d2','#fc9272','#de2d26')) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")


# plot all samples 


# colData(tse_s) <- colData(tse_s) %>% as.data.frame() %>%
#   rownames_to_column("rntemp") %>%
#   mutate(test, Sample = ifelse(pre & week == "18", "M1", ifelse(
#     pre & week == "32", "M2", ifelse(origin == "m" & week == "32", "M3", ifelse(
#       week == "2", "I1", ifelse(week == "6", "I2", ifelse(week == "12", "I3", ifelse(
#         origin == "i" & week == "32", "I4", NA)))))))) %>%
#   column_to_rownames("rntemp") %>%
#   DataFrame()


# Beta diversity metrics like Bray-Curtis are often applied to relabundances
tse_s <- transformAssay(tse_s,
                        assay.type = "counts",
                        method = "relabundance")

# Run PCoA on relabundance assay with Bray-Curtis distances
tse_s <- runMDS(tse_s,
                FUN = vegan::vegdist,
                method = "bray",
                assay.type = "relabundance",
                name = "MDS_bray")

# Create ggplot object
p <- plotReducedDim(
  tse_s, "MDS_bray",
  colour_by = "t",
  point_alpha = 1,
  point_size = 3
) +
  scale_color_manual(values = c('#fee0d2','#fc9272','#de2d26', '#edf8e9','#bae4b3','#74c476','#238b45')) +
  theme_bw(base_size = 15)

# Calculate explained variance
e <- attr(reducedDim(tse_s, "MDS_bray"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")) 

p



tse_s <- runMDS(tse_s,
                FUN = mia::calculateUnifrac,
                name = "Unifrac",
                tree = rowTree(tse_s),
                ntop = nrow(tse_s),
                assay.type = "counts")

plotReducedDim(
  tse_s, "Unifrac", colour_by = "t",
  point_alpha = 1, 
  point_size = 2
) +
  scale_color_manual(values = c('#fee0d2','#fc9272','#de2d26', '#edf8e9','#bae4b3','#74c476','#238b45')) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")



# import biplot function
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/mb_helper.R")

# i perform PCA on all data first, then only maternal samples
pseq <- makePhyloseqFromTreeSE(tse_s)
pseq_clr <- microbiome::transform(pseq, transform = "clr")
# sample_data(pseq_clr) <- sd_to_df(pseq_clr) %>%
#   mutate(
#     Infancy = ifelse(week == 52, "Late", "Early"),
#     Condition = ifelse(condition == 1, "SSC", ifelse(condition == 0, "CAU", NA))
#   ) %>%
#   df_to_sd()

sids <- colData(tse_s) %>% 
  as.data.frame() %>%
  filter(origin == "m") %>%
  .$sid

# all samples
bp <- biplot(
  pseq_clr, 
  color = "t", 
  point_size = 3, 
  otu_alpha = 0,
  #colors = c("#fc8d62", "#8da0cb"),
  colors = c('#fee0d2','#fc9272','#de2d26'),
  filter_samples = sids,
  #shape = "origin"
)


p <- bp[[1]] + theme_bw(base_size = 15) + 
  labs(color = "Sample") +
  scale_color_manual(
    values = c("t1" = '#fee0d2', "t2" = '#fc9272', "t3" = '#de2d26'),
    labels = c("M1", "M2", "M3"))
p

# only mothers
pseq <- makePhyloseqFromTreeSE(tse_m)
pseq_clr <- microbiome::transform(pseq, transform = "clr")
bp <- biplot(
  pseq_clr, 
  color = "t", 
  point_size = 3, 
  otu_alpha = 0,
  #colors = c("#fc8d62", "#8da0cb"),
  colors = c('#fee0d2','#fc9272','#de2d26'),
  #filter_samples = sids,
  #shape = "origin"
)


p <- bp[[1]] + theme_bw(base_size = 15) + 
  labs(color = "Sample") +
  scale_color_manual(
    values = c("t1" = '#fee0d2', "t2" = '#fc9272', "t3" = '#de2d26'),
    labels = c("M1", "M2", "M3"))
p

save(p, p_uni, p_bray, file = here("data/rdata/biplot_m.Rds"))
ggsave(
  filename = here("fig/fig2.png"),
  unit = "in",
  height = 10,
  width = 10
  )


# # PCA on maternal samples only
# pseq <- makePhyloseqFromTreeSE(tse_m)
# pseq_clr <- microbiome::transform(pseq, transform = "clr")



# show intra individual distances from t1-t2 and t2 to t3
sids <- colData(tse_s) %>% 
  as.data.frame() %>%
  filter(origin == "m", t != "t3") %>%
  .$sid
sample_data(pseq_clr) <- sd_to_df(pseq_clr) %>%
  rename(subject_id = id) %>%
  df_to_sd()
bp_t1_t2 <- biplot(
  pseq_clr, 
  #color = "t", 
  point_size = 3, 
  alpha = 0,
  otu_alpha = 0,
  #colors = c("#fc8d62", "#8da0cb"),
  #colors = c('#fee0d2','#fc9272','#de2d26'),
  filter_samples = sids,
  connect_series = "subject_id"
)
bp[[1]]

# show intra individual distances from t2 to t3
sids <- colData(tse_s) %>% 
  as.data.frame() %>%
  filter(origin == "m", t != "t1") %>%
  .$sid

bp_t2_t3 <- biplot(
  pseq_clr, 
  #color = "t", 
  point_size = 3, 
  alpha = 0,
  otu_alpha = 0,
  #colors = c("#fc8d62", "#8da0cb"),
  #colors = c('#fee0d2','#fc9272','#de2d26'),
  filter_samples = sids,
  connect_series = "subject_id"
)
p <- bp_t1_t2[[1]] + ylim(c(-45, 40)) + theme_bw(base_size = 15) + 
  bp_t2_t3[[1]] + ylim(c(-45, 40)) + ylab("") + theme_bw(base_size = 15)
p
ggsave(
  filename = here("fig/fig3.png"),
  unit = "in",
  height = 10,
  width = 20
)


# no visual separation based on time point indicating that at least no drastic
# shift occur between 18, 32 week prenatally and 8 months postnatally
# often we can detect smaller differences statistically, lets see:


# do we detect differences between t1 and t2?
pm_df <- map_dfr(unique(dm$t), function(time) {
# extract relevant meta data and omit na as adonis doesnt accept them.
  meta <- colData(tse_m) %>% as.data.frame() %>%
    rownames_to_column("sid") %>%
    select(id, sid, t, pbmi, parity, ms, stai, epds, pss, 
          edu, activity) %>%
    filter(t != time) %>%
    na.omit()
  # we need to account for non-independence of data in the infancy model 
  ids <- meta %>% mutate(id = as.factor(id)) %>%
    .$id
  h <- how(plots = Plots(strata = ids, type = "none"),
          nperm = 999)
  # according to omitted NAs I need to select stool samples 
  asv <- t(assay(tse_m, "clr"))
  asv <- asv[meta$sid, ]
  colnames(t(asv))
  head(asv)
  dim(asv)
  dim(meta)
  # fit and inspect model 
  permanova <- adonis2(asv ~ t,
                      # by = "margin", # each term analyzed individually
                      data = meta,
                      method = "euclidean",
                      # h does not work if trend is in data (therefore use 999), see Gavins post
                      permutations = h
  )

  df_temp <- as.data.frame(permanova) %>%
    rownames_to_column("parameter") %>%
    mutate(t = ifelse(time == "t1", "t2 vs t3", ifelse(time == "t2", "t1 vs t3", "t1 vs t2"))) %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    filter(parameter == "t") %>%
    select(-parameter)

  colnames(df_temp) <- c("df", "sqs", "R2", "F", "p", "t")
  df_temp

})
pm_df
save(pm_df, file = here("data/rdata/pm_df_t.Rds"))

# how do results differ when permuting all observations?
pm_df <- map_dfr(unique(dm$t), function(time) {
  # extract relevant meta data and omit na as adonis doesnt accept them.
  meta <- colData(tse_m) %>% as.data.frame() %>%
    rownames_to_column("sid") %>%
    select(id, sid, t, pbmi, parity, ms, stai, epds, pss, 
           edu, activity) %>%
    filter(t != time) %>%
    na.omit()
  # according to omitted NAs I need to select stool samples 
  asv <- t(assay(tse_m, "clr"))
  asv <- asv[meta$sid, ]
  colnames(t(asv))
  head(asv)
  dim(asv)
  dim(meta)
  # fit and inspect model 
  permanova <- adonis2(asv ~ t,
                       # by = "margin", # each term analyzed individually
                       data = meta,
                       method = "euclidean",
                       # h does not work if trend is in data (therefore use 999), see Gavins post
                       permutations = 999
  )
  
  df_temp <- as.data.frame(permanova) %>%
    rownames_to_column("parameter") %>%
    mutate(t = ifelse(time == "t1", "t2 vs t3", ifelse(time == "t2", "t1 vs t3", "t1 vs t2"))) %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    filter(parameter == "t") %>%
    select(-parameter)
  
  colnames(df_temp) <- c("df", "sqs", "R2", "F", "p", "t")
  df_temp
  
})
pm_df
save(pm_df, file = here("data/rdata/pm_df_t_noh.Rds"))

colData(tse_m) %>% as.data.frame() %>% colnames()

stress_indices <- c("ms", "stai", "epds", "pss", "praqr_handicap", "praqr_birth",
                    "cortisol", "cortisone", "ratio"
                    )
head(colData(tse_m) %>% as.data.frame())
# is stress related to microbiota at each time point? I test also at separate time points because it
# is not clear if PERMANOVA works well with dependence between samples, at least I need to double check results
pm_df <- map_dfr(unique(dm$t), function(time) {
  map_dfr(stress_indices, function(stress) {
  if (!(time == "t3" & str_detect(stress, "praqr"))) {
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_m) %>% as.data.frame() %>%
      rownames_to_column("sid") %>%
      select(id, sid, t, pbmi, parity, all_of(stress),
             edu, activity, age) %>%
      filter(t == time) %>%
      na.omit()
    
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_m, "clr"))
    asv <- asv[meta$sid, ]
    # fit and inspect model 
    form <- as.formula(glue("asv ~ activity + pbmi + {stress} + edu + parity + age"))
    permanova <- adonis2(form,
                         # by = "margin", # each term analyzed individually
                         data = meta,
                         method = "euclidean",
                         # h does not work if trend is in data (therefore use 999), see Gavins post
                         permutations = 999
    )
    
    df_temp <- as.data.frame(permanova) %>%
      rownames_to_column("paramter") %>%
      mutate(t = time, stress = stress)
    
    colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p", "t", "stress")
    out <- filter(df_temp, parameter == stress)
    return(out)
  }
  })
})

pm_df


# do results differ if we utilize maximum power?
pm_df <- map_dfr(stress_indices, function(stress) {
  if (str_detect(stress, "praqr")) {
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_m) %>% as.data.frame() %>%
      rownames_to_column("sid") %>%
      select(id, sid, t, pbmi, parity, all_of(stress),
             edu, activity, age) %>%
      filter(t != "t3") %>%
      na.omit()
  } else {
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_m) %>% as.data.frame() %>%
      rownames_to_column("sid") %>%
      select(id, sid, t, pbmi, parity, all_of(stress),
             edu, activity, age) %>%
      na.omit()
  }


  # according to omitted NAs I need to select stool samples 
  asv <- t(assay(tse_m, "clr"))
  asv <- asv[meta$sid, ]
  # fit and inspect model 
  form <- as.formula(glue("asv ~ activity + pbmi + t * {stress} + edu + parity + age"))
  permanova <- adonis2(form,
                      # by = "margin", # each term analyzed individually
                      data = meta,
                      method = "euclidean",
                      # h does not work if trend is in data (therefore use 999), see Gavins post
                      permutations = 999
  )

  df_temp <- as.data.frame(permanova) %>%
    rownames_to_column("paramter") %>%
    mutate(stress = stress)
  
  colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p", "stress")
  filter(df_temp, parameter == stress)
})

pm_df

# do results differ if we utilize maximum power? Here we account for non-independence
pm_df <- map_dfr(stress_indices, function(stress) {
  
  if (str_detect(stress, "praqr")) {
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_m) %>% as.data.frame() %>%
      rownames_to_column("sid") %>%
      select(id, sid, t, pbmi, parity, all_of(stress),
             edu, activity, age) %>%
      filter(t != "t3") %>%
      na.omit()
  } else {
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_m) %>% as.data.frame() %>%
      rownames_to_column("sid") %>%
      select(id, sid, t, pbmi, parity, all_of(stress),
             edu, activity, age) %>%
      na.omit()
  }

  # we need to account for non-independence of data in the infancy model 
  ids <- meta %>% mutate(id = as.factor(id)) %>%
    .$id
  h <- how(plots = Plots(strata = ids, type = "none"),
          nperm = 999)
  # according to omitted NAs I need to select stool samples 
  asv <- t(assay(tse_m, "clr"))
  asv <- asv[meta$sid, ]
  # fit and inspect model 
  form <- as.formula(glue("asv ~ activity + pbmi + t * {stress} + edu + parity + age"))
  permanova <- adonis2(form,
                      # by = "margin", # each term analyzed individually
                      data = meta,
                      method = "euclidean",
                      # h does not work if trend is in data (therefore use 999), see Gavins post
                      permutations = h
  )

  df_temp <- as.data.frame(permanova) %>%
    rownames_to_column("paramter") %>%
    mutate(stress = stress)
  
  colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p", "stress")
  #filter(df_temp, parameter == stress)
  df_temp
})



# now we need to explore psas
# do results differ if we utilize maximum power? Here we account for non-independence


# extract relevant meta data and omit na as adonis doesnt accept them.
psas$id <- as.character(psas$id)
meta <- colData(tse_m) %>% as.data.frame() %>%
  rownames_to_column("sid") %>%
  select(id, sid, t, pbmi, parity,
         edu, activity, age) %>%
  filter(t == "t3") %>%
  left_join(psas, by = "id") %>%
  na.omit()

  

# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_m, "clr"))
asv <- asv[meta$sid, ]
# fit and inspect model 
form <- as.formula(glue("asv ~ activity + pbmi + psas + edu + parity + age"))
permanova <- adonis2(form,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "euclidean",
                     # h does not work if trend is in data (therefore use 999), see Gavins post
                     permutations = 999
)

df_temp <- as.data.frame(permanova) %>%
  rownames_to_column("paramter") 

colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p")
#filter(df_temp, parameter == stress)
df_temp


# do we find an interaction between parity and pbmi?
  # extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_m) %>% as.data.frame() %>%
  rownames_to_column("sid") %>%
  select(id, sid, t, pbmi, parity,
        edu, activity, age) %>%
  na.omit()

# we need to account for non-independence of data in the infancy model 
ids <- meta %>% mutate(id = as.factor(id)) %>%
  .$id
h <- how(plots = Plots(strata = ids, type = "none"),
        nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_m, "clr"))
asv <- asv[meta$sid, ]
# fit and inspect model 
form <- as.formula(glue("asv ~ activity + pbmi * parity + t + edu + age"))
permanova <- adonis2(form,
                    # by = "margin", # each term analyzed individually
                    data = meta,
                    method = "euclidean",
                    # h does not work if trend is in data (therefore use 999), see Gavins post
                    permutations = h
)
permanova

# what if we adjust to permutations across?

permanova <- adonis2(form,
                     # by = "margin", # each term analyzed individually
                     data = meta,
                     method = "euclidean",
                     # h does not work if trend is in data (therefore use 999), see Gavins post
                     permutations = 999
)
permanova

# there seems to be an interaction between pbmi and parity across time such that
# the microbiota composition either of them only matters if the other is set to
# a certain value, each pbmi may only matter if parity is firstchild. 

# lets see if this is only one time points or all:

pm_df <- map_dfr(unique(dm$t), function(time) {
      # extract relevant meta data and omit na as adonis doesnt accept them.
      meta <- colData(tse_m) %>% as.data.frame() %>%
        rownames_to_column("sid") %>%
        select(id, sid, t, pbmi, parity,
               edu, activity, age) %>%
        filter(t == time) %>%
        na.omit()
      
      # according to omitted NAs I need to select stool samples 
      asv <- t(assay(tse_m, "clr"))
      asv <- asv[meta$sid, ]
      # fit and inspect model 
      form <- as.formula(glue("asv ~ activity + pbmi * parity + edu + age"))
      permanova <- adonis2(form,
                           # by = "margin", # each term analyzed individually
                           data = meta,
                           method = "euclidean",
                           # h does not work if trend is in data (therefore use 999), see Gavins post
                           permutations = 999
      )
      
      df_temp <- as.data.frame(permanova) %>%
        rownames_to_column("paramter") %>%
        mutate(t = time)
      
      colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p", "t")
      out <- df_temp
      return(out)
})
pm_df





# this section is used to explore the interaction
  # extract relevant meta data and omit na as adonis doesnt accept them.
meta <- colData(tse_m) %>% as.data.frame() %>%
  rownames_to_column("sid") %>%
  select(id, sid, t, pbmi, parity,
        edu, activity, age) %>%
  filter(parity == "0") %>%
  na.omit()

# we need to account for non-independence of data in the infancy model 
ids <- meta %>% mutate(id = as.factor(id)) %>%
  .$id
h <- how(plots = Plots(strata = ids, type = "none"),
        nperm = 999)
# according to omitted NAs I need to select stool samples 
asv <- t(assay(tse_m, "clr"))
asv <- asv[meta$sid, ]
dim(asv)
# fit and inspect model 
form <- as.formula(glue("asv ~ activity + pbmi  + t + edu + age"))
permanova <- adonis2(form,
                    # by = "margin", # each term analyzed individually
                    data = meta,
                    method = "euclidean",
                    # h does not work if trend is in data (therefore use 999), see Gavins post
                    permutations = h
)
permanova





# lastly cortisol and cortisone 
pm_df <- map_dfr(c("cortisol", "cortisone", "ratio"), function(stress) {
  
  # extract relevant meta data and omit na as adonis doesnt accept them.
  meta <- colData(tse_m) %>% as.data.frame() %>%
    rownames_to_column("sid") %>%
    select(id, sid, t, pbmi, parity, all_of(stress),
           edu, activity, age, hw) %>%
    filter(t != "t3") %>%
    na.omit()

  
  # we need to account for non-independence of data in the infancy model 
  ids <- meta %>% mutate(id = as.factor(id)) %>%
    .$id
  h <- how(plots = Plots(strata = ids, type = "none"),
           nperm = 999)
  # according to omitted NAs I need to select stool samples 
  asv <- t(assay(tse_m, "clr"))
  asv <- asv[meta$sid, ]
  # fit and inspect model 
  form <- as.formula(glue("asv ~ activity + pbmi + t * {stress} + edu + parity + age + hw"))
  permanova <- adonis2(form,
                       # by = "margin", # each term analyzed individually
                       data = meta,
                       method = "euclidean",
                       # h does not work if trend is in data (therefore use 999), see Gavins post
                       permutations = h
  )
  
  df_temp <- as.data.frame(permanova) %>%
    rownames_to_column("paramter") %>%
    mutate(stress = stress)
  
  colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p", "stress")
  #filter(df_temp, parameter == stress)
  df_temp
})

pm_df
# no significant association with cortisol or cortisone




# i perform PCA on all data first, then only maternal samples
pseq <- makePhyloseqFromTreeSE(tse_m)
pseq_clr <- microbiome::transform(pseq, transform = "clr")

sids <- colData(tse_m) %>% 
  as.data.frame() %>%
  rownames_to_column("sid") %>%
  mutate(
  lowbmi = ifelse(pbmi <= quantile(pbmi, 0.33, na.rm = TRUE), 1, ifelse(pbmi >= quantile(pbmi, 0.66, na.rm = TRUE), 0, NA))
  ) %>%
  filter(origin == "m", !is.na(lowbmi)) %>%
  .$sid

sample_data(pseq_clr) <- sd_to_df(pseq_clr) %>%
  mutate(
    lowbmi = ifelse(pbmi <= quantile(pbmi, 0.33, na.rm = TRUE), 1, ifelse(pbmi >= quantile(pbmi, 0.66, na.rm = TRUE), 0, NA))
  ) %>%
  df_to_sd()

# all samples
bp <- biplot(
  pseq_clr, 
  color = "lowbmi", 
  point_size = 5, 
  otu_alpha = 0,
  #colors = c("#fc8d62", "#8da0cb"),
  colors = c('#e41a1c','#377eb8','#4daf4a'),
  facet = "firstchild", 
  filter_samples = sids
  #shape = "firstchild",
  
)
bp[[1]]
sample_data(pseq_clr)

# we see maybe a little bit more separation between low and high bmi in the firstborn group 

