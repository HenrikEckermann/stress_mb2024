library(mia)
library(tidyverse)
library(tidySummarizedExperiment)
library(ape)
library(here)



### taxonomic abundances
path <- here("data/processed_by_turku/joined_metaphlan_table.tsv")
tse <- loadFromMetaphlan(path)

# Reassigning rownames to match with tree labels
rownames(tse) <- rowData(tse)$clade_name
tree <- read.tree(here("data/processed_by_turku/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk"))
# the tree tip labels are the SGB ids. In order to map we need the following files:
compressed_file_path <- here("data/processed_by_turku/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2")
species_names <- read_delim(compressed_file_path, col_names = FALSE)
colnames(species_names) <- c("sgb_id", "name")
filter(species_names, str_detect(sgb_id, "SGB225"))
species_names <- mutate(
  species_names, 
  sgb_id = str_remove(sgb_id, "SGB"),
  sgb_id = str_remove(sgb_id, "_group")
  ) %>%
  filter(sgb_id %in% tree$tip.label)
# check if we have all names available
length(tree$tip.label) - length(species_names$sgb_id)
# before, if I did not remove _group as well it did not work
tip_labels <- map_chr(tree$tip.label, function(tl) {
  species_names$name[species_names$sgb_id == tl]
})
rowTree(tse) <- tree


# derive metadata from ids and create new sample ids that are easier to read

# the subject ids are in the end of the current id names
rx <- "[MC]FE[p1682][2M31w][wO28]?(\\d{3})_"
id <- str_match(colnames(tse), rx)[, 2]
# for the sampling time
week <- ifelse(str_detect(colnames(tse), "E2w"), 2, ifelse(
  str_detect(colnames(tse), "E12w"), 12, ifelse(
    str_detect(colnames(tse), "Ep32"), 32, ifelse(
      str_detect(colnames(tse), "Ep18"), 18, ifelse(
        str_detect(colnames(tse), "E6w"), 6, ifelse(
          str_detect(colnames(tse), "E8MO"), 32, NA))))))
# maternal vs infant sample
origin <- ifelse(str_starts(colnames(tse), pattern = "M"), "m", "i")
# prenatal vs postnatal
prepost <- ifelse(str_detect(colnames(tse), "p[31][28]"), "pre", "post")

md <- tibble(
  id,
  origin,
  prepost,
  week,
  old_sample_ids = colnames(tse)
  ) %>%
    mutate(sid = glue::glue('{id}{origin}{prepost}{week}'))


colData(tse) <- colData(tse) %>%
                  as.data.frame() %>%
                  rownames_to_column("old_sample_ids") %>%
                  full_join(md, by = "old_sample_ids") %>%
                  column_to_rownames("old_sample_ids") %>%
                  DataFrame()

colnames(tse) <- colData(tse)$sid

# how many samples per time point
colData(tse) %>%
  as.data.frame() %>%
  count(origin, prepost, week)

# if there are duplicate samples we must disregard one
sids <- colData(tse) %>%
  as.data.frame() %>%
  rownames_to_column("sid_tse") %>%
  group_by(sid) %>%
  summarise(sid_tse = sid_tse, n = n()) %>%
  filter(n > 1) %>%
  .$sid_tse
# there are 2 duplicate samples. one sample of a mother, another of an infant
# in order to see if one of the is an outlier lets plot them

# to proceed we must have unique sample ids,
# we can for now take the unique ids created automatically by this operation
usid <- colData(tse) %>%
  as.data.frame() %>%
  rownames_to_column("sid_tse") %>%
  .$sid_tse

table(colnames(tse)) %>% as.data.frame() %>% filter(Freq > 1)
colnames(tse) <- usid

source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/mb_helper.R")
tse_g <- mergeFeaturesByRank(tse, rank = "Genus")
pseq <- makePhyloseqFromTreeSE(tse_g)
pseq_clr <- microbiome::transform(pseq, transform = "clr")
sids
# all samples
bp <- biplot(
  pseq_clr, 
  filter_samples = sids,
  color = "sid",
  point_size = 5, 
  otu_alpha = 0,
  text = TRUE,
  colors = c("#fc8d62", "#8da0cb")
)
bp[[1]]
# the samples cluster together closely indicating that we can use either sample
colData(tse)$sid <- colnames(tse)
dim(tse)

tse <- filter(tse, !sid %in% c("554mpre18.1", "636ipost6.1"))
dim(tse)
# now the double samples are filtered out and we can proceed creating a tse 


save(tse, file = here::here("data/rdata/tse.Rds"))


# for Eva project specifically

tse_i <- tse[, colData(tse)$origin == "i"]


# combine with metadata
md2 <- read_csv(here::here("data/dlong.csv")) %>%
  mutate(sample_id = glue::glue('{id}ipost{week}')) %>%
  select(-id, -week)

colData(tse_i) <- colData(tse_i) %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  select(-sid) %>%
  left_join(md2, by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  DataFrame()


colData(tse_i) %>%
  as.data.frame() %>%
  count(origin, prepost, week)

# for the mbage model we need one at genus level and prevalence filtered
tse_smiley <- agglomerateByRank(tse_i, rank = "Genus")
tse_smiley <- subsetByPrevalentFeatures(tse_smiley, rank = "Genus", prevalence = 5 / 100)
dim(tse_smiley)
save(tse_smiley, file = here::here("data/rdata/tse_smiley.Rds"))



### functional abundance
path_f <- here::here("data/processed_by_turku/joined_pathabundance.tsv")
tse_f <- loadFromHumann(path_f)
af <- assay(tse_f)
rownames(af)
colnames(af)



