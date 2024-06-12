library(mia)
library(tidyverse)
library(vegan)
library(here)
library(glue)
library(scater)
library(patchwork)

load(file = here::here("data/rdata/tse.Rds"))
load(here("data/rdata/di.Rds"))

# what are lowest and highest counts in this dataset?
assay(tse_s) %>% as.data.frame() %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  filter(value != 0) %>%
  summarise(min = min(value), max = max(value))
tse <- transformAssay(tse_s, method = "clr", name = "clr", pseudocount = 0.000001)


# infants
tse_i <- tse[, colData(tse)$origin == "i"]
# exclude two year samples because they arrived to late for this project
tse_i <- tse_i[, colData(tse_i)$week != 104]
dlong$id <- as.character(dlong$id)
dlong$sid <- glue("{dlong$id}iFALSE{dlong$week}")
colData(tse_i) <- colData(tse_i) %>%
                    as.data.frame() %>%
                    rownames_to_column("sid2") %>%
                    select(-id, -t, -pre, -week) %>%
                    left_join(dlong, by = "sid") %>%
                    column_to_rownames("sid2") %>%
                    DataFrame()

colData(tse_i)$week <- factor(colData(tse_i)$week, levels = c("2", "6", "12", "32"))

# visualize the data
tse_i <- transformAssay(tse_i,
                      assay.type = "counts",
                      method = "relabundance")

# Run PCoA on relabundance assay with Bray-Curtis distances
tse_i <- runMDS(tse_i,
              FUN = vegan::vegdist,
              method = "bray",
              assay.type = "relabundance",
              name = "MDS_bray")

# Create ggplot object
p <- plotReducedDim(
  tse_i, "MDS_bray",
  colour_by = "week",
  point_alpha = 1,
  point_size = 2
  ) +
  scale_color_manual(values = c('#edf8e9','#bae4b3','#74c476','#238b45')) +
  theme_bw(base_size = 15)

# Calculate explained variance
e <- attr(reducedDim(tse_i, "MDS_bray"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")) 

p

tse_i <- runMDS(tse_i,
              FUN = mia::calculateUnifrac,
              name = "Unifrac",
              tree = rowTree(tse_i),
              ntop = nrow(tse_i),
              assay.type = "counts")

plotReducedDim(
  tse_i, "Unifrac", colour_by = "week",
  point_alpha = 1, 
  point_size = 2
  ) +
  scale_color_manual(values = c('#edf8e9','#bae4b3','#74c476','#238b45')) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")

# import biplot function
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/mb_helper.R")

# i perform PCA on all data first, then only maternal samples
pseq <- makePhyloseqFromTreeSE(tse_i)
pseq_clr <- microbiome::transform(pseq, transform = "clr")
sids <- colData(tse_i) %>% 
  as.data.frame() %>%
  filter(origin == "i") %>%
  .$sid

# all samples
bp <- biplot(
  pseq_clr, 
  color = "week", 
  point_size = 3, 
  otu_alpha = 0,
  #colors = c("#fc8d62", "#8da0cb"),
  colors = c('#edf8e9','#bae4b3','#74c476','#238b45'),
  filter_samples = sids,
  #shape = "origin"
)
p <- bp[[1]] + theme_bw(base_size = 15) + theme(legend.position = "none")
p

# perform permanova analyses
stress_indices <- c("ms", "stai", "epds", "pss", "praq_worries_handicap_g", 
                    "praq_fear_birth_g", "cortisol_g1", "cortisol_g2", 
                    "cortisol_g3", "cortisol_pp", "cortisone_g1", "cortisone_g2", 
                    "cortisone_g3", "cortisone_pp", "psas", "ratio_g1", "ratio_g2",
                    "ratio_g3")
# is stress related to microbiota at each time point? I test also at separate time points because it
# is not clear if PERMANOVA works well with dependence between samples, at least I need to double check results

if (!file.exists(here("data/rdata/stress_infant_pm.Rds"))) {
  pm_df <- map_dfr(unique(dlong$week), function(weekn) {
    map_dfr(stress_indices, function(stress) {
      if (!(weekn == "32" & stress == "psas")) {
        if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
          form <- as.formula(glue("asv ~ {stress} + parity + dmode + gage + hw"))
          # extract relevant meta data and omit na as adonis doesnt accept them.
          meta <- colData(tse_i) %>% as.data.frame() %>%
            rownames_to_column("sid2") %>%
            select(id, sid, week, parity, dmode, gage, all_of(stress), hw) %>%
            filter(week == weekn) %>%
            na.omit()
        } else {
          form <- as.formula(glue("asv ~ {stress} + parity + dmode + gage"))
          # extract relevant meta data and omit na as adonis doesnt accept them.
          meta <- colData(tse_i) %>% as.data.frame() %>%
            rownames_to_column("sid2") %>%
            select(id, sid, week, parity, dmode, gage, all_of(stress)) %>%
            filter(week == weekn) %>%
            na.omit()
        }

        
        # according to omitted NAs I need to select stool samples 
        asv <- t(assay(tse_i, "clr"))
        asv <- asv[meta$sid, ]
        # fit and inspect model 

        
        permanova <- adonis2(form,
                             # by = "margin", # each term analyzed individually
                             data = meta,
                             method = "euclidean",
                             # h does not work if trend is in data (therefore use 999), see Gavins post
                             permutations = 999
        )
        
        df_temp <- as.data.frame(permanova) %>%
          rownames_to_column("parameter") %>%
          mutate(week = weekn, stress = stress)
        
        colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p", "week", "stress")
        out <- filter(df_temp, parameter == stress)
        return(out)
      }
    })
  })
  
  save(pm_df, file = here("data/rdata/stress_infant_pm.Rds"))
} else {
  load(here("data/rdata/stress_infant_pm.Rds"))
}


filter(pm_df, p < 0.05)



if (!file.exists(here("data/rdata/stress_infant_allsamples_pm.Rds"))) {
# do results differ if we utilize maximum power?
pm_df <- map_dfr(stress_indices, function(stress) {
  
  
  if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
    form <- as.formula(glue("asv ~ week * {stress} + dmode + gage + parity + hw"))
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_i) %>% as.data.frame() %>%
      rownames_to_column("sid2") %>%
      select(id, sid, week, parity, gage, dmode, all_of(stress), hw) %>%
      na.omit()
  } else {
    form <- as.formula(glue("asv ~ week * {stress} + dmode + gage + parity"))
    # extract relevant meta data and omit na as adonis doesnt accept them.
    meta <- colData(tse_i) %>% as.data.frame() %>%
      rownames_to_column("sid2") %>%
      select(id, sid, week, parity, gage, dmode, all_of(stress)) %>%
      na.omit()
  }
  



  # according to omitted NAs I need to select stool samples 
  asv <- t(assay(tse_i, "clr"))
  asv <- asv[meta$sid, ]
  # fit and inspect model 

  
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
save(pm_df, file = here("data/rdata/stress_infant_allsamples_pm.Rds"))

} else {
  load(here("data/rdata/stress_infant_allsamples_pm.Rds"))
}
filter(pm_df, p < 0.05)

# do results differ if we utilize maximum power? Here we account for non-independence
if (!file.exists(here("data/rdata/stress_infant_allsamples_nonindependence_pm.Rds"))) {
  pm_df <- map_dfr(stress_indices, function(stress) {
    
    if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
      form <- as.formula(glue("asv ~week * {stress} + dmode + parity + gage + hw"))
      # extract relevant meta data and omit na as adonis doesnt accept them.
      meta <- colData(tse_i) %>% as.data.frame() %>%
        rownames_to_column("sid2") %>%
        select(id, sid, week, parity, gage, dmode, all_of(stress), hw) %>%
        na.omit()
    } else {
      form <- as.formula(glue("asv ~week * {stress} + dmode + parity + gage"))
      # extract relevant meta data and omit na as adonis doesnt accept them.
      meta <- colData(tse_i) %>% as.data.frame() %>%
        rownames_to_column("sid2") %>%
        select(id, sid, week, parity, gage, dmode, all_of(stress)) %>%
        na.omit()
    }
    
  

  
    # we need to account for non-independence of data in the infancy model 
    ids <- meta %>% mutate(id = as.factor(id)) %>%
      .$id
    h <- how(plots = Plots(strata = ids, type = "none"),
            nperm = 999)
    # according to omitted NAs I need to select stool samples 
    asv <- t(assay(tse_i, "clr"))
    asv <- asv[meta$sid, ]
    # fit and inspect model 

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
    filter(df_temp, parameter == stress)
  })
  save(pm_df, file = here("data/rdata/stress_infant_allsamples_nonindependence_pm.Rds"))
} else {
  load(here("data/rdata/stress_infant_allsamples_nonindependence_pm.Rds"))
}

filter(pm_df, p < 0.05)





# now we ask the same questions for prenatal stress

infant_time <- c("2", "6", "12", "32")
mother_time <- c("18", "32")
# i leave out cort because we covered that above
stress_indices <- c("ms", "stai", "epds", "pss", "praqr_handicap", 
                    "praqr_birth")
load(here("data/rdata/dm.Rds"))

if (!file.exists(here("data/rdata/pn_stress_infant_pm.Rds"))) {
  pm_df <- map_dfr(infant_time, function(it) {
    map_dfr(mother_time, function(mt) {
      map_dfr(stress_indices, function(stress) {
        if (!(it == "32" & stress == "psas")) {
          if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
            form <- as.formula(glue("asv ~ {stress} + parity + dmode + gage + hw"))
            # extract relevant meta data and omit na as adonis doesnt accept them.
            meta <- colData(tse_i) %>% as.data.frame() %>%
              rownames_to_column("sid2") %>%
              select(id, sid, week, parity, dmode, gage, hw) %>%
              filter(week == it) %>%
              na.omit()
            
            temp <- filter(dm, week == mt, pre) %>% 
              select(id, all_of(stress)) %>%
              mutate(id = as.character(id))
            
            
            meta <- left_join(meta, temp, by = "id") %>% na.omit()
          } else {
            form <- as.formula(glue("asv ~ {stress} + parity + dmode + gage"))
            # extract relevant meta data and omit na as adonis doesnt accept them.
            meta <- colData(tse_i) %>% as.data.frame() %>%
              rownames_to_column("sid2") %>%
              select(id, sid, week, parity, dmode, gage) %>%
              filter(week == it) %>%
              na.omit()
            
            temp <- filter(dm, week == mt, pre) %>% 
              select(id, all_of(stress)) %>%
              mutate(id = as.character(id))
            
            
            meta <- left_join(meta, temp, by = "id") %>% na.omit()
          }
  
          
          # according to omitted NAs I need to select stool samples 
          asv <- t(assay(tse_i, "clr"))
          asv <- asv[meta$sid, ]
          # fit and inspect model 
  
          
          permanova <- adonis2(form,
                               # by = "margin", # each term analyzed individually
                               data = meta,
                               method = "euclidean",
                               # h does not work if trend is in data (therefore use 999), see Gavins post
                               permutations = 999
          )
          
          df_temp <- as.data.frame(permanova) %>%
            rownames_to_column("parameter") %>%
            mutate(infant_time = it, mother_time = mt, stress = stress)
          
          colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p", "infant_time", "mother_time", "stress")
          out <- filter(df_temp, parameter == stress)
          return(out)
        }
      })
    })
  })
  save(pm_df, file = here("data/rdata/pn_stress_infant_pm.Rds"))
} else {
  load(here("data/rdata/pn_stress_infant_pm.Rds"))
}

filter(pm_df, p < 0.05)



if (!file.exists(here("data/rdata/pn_stress_infant_allsamples_pm.Rds"))) {
  # do results differ if we utilize maximum power?
  pm_df <- map_dfr(stress_indices, function(stress) {
    map_dfr(mother_time, function(mt) {
      if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
        form <- as.formula(glue("asv ~ week * {stress} + dmode + gage + parity + hw"))
        # extract relevant meta data and omit na as adonis doesnt accept them.
        meta <- colData(tse_i) %>% as.data.frame() %>%
          rownames_to_column("sid2") %>%
          select(id, sid, week, parity, gage, dmode, hw) %>%
          na.omit()
        
        temp <- filter(dm, week == mt, pre) %>% 
          select(id, all_of(stress)) %>%
          mutate(id = as.character(id))
        
        meta <- left_join(meta, temp, by = "id") %>% na.omit()
      } else {
        form <- as.formula(glue("asv ~ week * {stress} + dmode + gage + parity"))
        # extract relevant meta data and omit na as adonis doesnt accept them.
        meta <- colData(tse_i) %>% as.data.frame() %>%
          rownames_to_column("sid2") %>%
          select(id, sid, week, parity, gage, dmode) %>%
          na.omit()
        
        temp <- filter(dm, week == mt, pre) %>% 
          select(id, all_of(stress)) %>%
          mutate(id = as.character(id))
        
        meta <- left_join(meta, temp, by = "id") %>% na.omit()
      }
      

      
      
      # according to omitted NAs I need to select stool samples 
      asv <- t(assay(tse_i, "clr"))
      asv <- asv[meta$sid, ]
      # fit and inspect model 

      permanova <- adonis2(form,
                           # by = "margin", # each term analyzed individually
                           data = meta,
                           method = "euclidean",
                           # h does not work if trend is in data (therefore use 999), see Gavins post
                           permutations = 999
      )
      
      df_temp <- as.data.frame(permanova) %>%
        rownames_to_column("paramter") %>%
        mutate(mother_time = mt, stress = stress)
      
      colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p", "mother_time", "stress")
      filter(df_temp, parameter == stress)
    })
  })
  save(pm_df, file = here("data/rdata/pn_stress_infant_allsamples_pm.Rds"))
} else {
  load(here("data/rdata/pn_stress_infant_allsamples_pm.Rds"))
}

filter(pm_df, p <= 0.05)

# do results differ if we utilize maximum power? Here we account for non-independence
if (!file.exists(here("data/rdata/stress_infant_allsamples_nonindependence_pm.Rds"))) {
  pm_df <- map_dfr(stress_indices, function(stress) {
    map_dfr(mother_time, function(mt) {
    
      # extract relevant meta data and omit na as adonis doesnt accept them.
      meta <- colData(tse_i) %>% as.data.frame() %>%
        rownames_to_column("sid2") %>%
        select(id, sid, week, parity, gage, dmode) %>%
        na.omit()
      
      temp <- filter(dm, week == mt, pre) %>% 
        select(id, all_of(stress)) %>%
        mutate(id = as.character(id))
      
      meta <- left_join(meta, temp, by = "id") %>% na.omit()
      
      # we need to account for non-independence of data in the infancy model 
      ids <- meta %>% mutate(id = as.factor(id)) %>%
        .$id
      h <- how(plots = Plots(strata = ids, type = "none"),
               nperm = 999)
      # according to omitted NAs I need to select stool samples 
      asv <- t(assay(tse_i, "clr"))
      asv <- asv[meta$sid, ]
      # fit and inspect model 
      if ((str_detect(stress, "cort") | str_detect(stress, "ratio"))) {
        form <- as.formula(glue("asv ~week * {stress} + dmode + parity + gage + hw"))
      } else {
        form <- as.formula(glue("asv ~week * {stress} + dmode + parity + gage"))
      }
      
      permanova <- adonis2(form,
                           # by = "margin", # each term analyzed individually
                           data = meta,
                           method = "euclidean",
                           # h does not work if trend is in data (therefore use 999), see Gavins post
                           permutations = h
      )
      
      df_temp <- as.data.frame(permanova) %>%
        rownames_to_column("paramter") %>%
        mutate(mother_time = mt, stress = stress)
      
      colnames(df_temp) <- c("parameter", "df", "sqs", "R2", "F", "p", "mother_time", "stress")
      filter(df_temp, parameter == stress)
    })
  })
  
  save(pm_df, file = here("data/rdata/stress_infant_allsamples_nonindependence_pm.Rds"))
} else {
  load(here("data/rdata/stress_infant_allsamples_nonindependence_pm.Rds"))
}

filter(pm_df, p <= 0.05)







# i perform PCA on all data first, then only maternal samples
pseq <- makePhyloseqFromTreeSE(tse_i)
pseq_clr <- microbiome::transform(pseq, transform = "clr")

sids <- colData(tse_i) %>% 
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

