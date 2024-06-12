library(mia)
library(tidySummarizedExperiment)
library(ape)
library(here)
library(tidyverse)

### taxonomic abundances
path1 <- here("data/processed_by_turku/joined_metaphlan_table.tsv")
tse1 <- loadFromMetaphlan(path1)

path2 <- here("data/processed_by_turku/merged_metaphlan_feb26.tsv")
tse2 <- loadFromMetaphlan(path2)


tse <- mergeSEs(tse1, tse2, missing_values = 0)
tse_s <- agglomerateByRank(tse, rank = "Species")
rn <- as.data.frame(rowData(tse_s))$clade_name %>% 
  map_chr(~unname(.x)) %>%
  str_remove("\\|t__.*$")

rownames(tse_s) <- rn
tree <- read.tree(here("data/processed_by_turku/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk"))
# the tree tip labels are the SGB ids. In order to map we need the following file:
file <- here("data/processed_by_turku/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt")
species_names <- read_delim(file, col_names = FALSE)
colnames(species_names) <- c("sgb_id", "name")

# these are the 5 taxa that remained left out after joining the tree
s <- c(
  "k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Debaryomycetaceae|g__Candida|s__Candida_parapsilosis",                                       
  "k__Eukaryota|p__Basidiomycota|c__Malasseziomycetes|o__Malasseziales|f__Malasseziaceae|g__Malassezia|s__Malassezia_globosa",                                        
  "k__Eukaryota|p__Basidiomycota|c__Malasseziomycetes|o__Malasseziales|f__Malasseziaceae|g__Malassezia|s__Malassezia_restricta",                                      
  "k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Debaryomycetaceae|g__Candida|s__Candida_albicans",                                           
  "k__Eukaryota|p__Eukaryota_unclassified|c__Eukaryota_unclassified|o__Eukaryota_unclassified|f__Eukaryota_unclassified|g__Blastocystis|s__Blastocystis_sp_subtype_1"
)
filter(species_names, name %in% s)


species_names <- mutate(
  species_names, 
  sgb_id = str_remove(sgb_id, "SGB"),
  sgb_id = str_remove(sgb_id, "_group"),
  name = str_split(name, ",") %>% map_chr(1) %>% map_chr(~unname(.x))
  ) %>%
  filter(sgb_id %in% tree$tip.label)


# check if we have all names available (before, if I did not remove _group this test failed)
length(tree$tip.label) - length(species_names$sgb_id)
# map names to tip labels in the order of the tip labels
tip_labels <- map_chr(tree$tip.label, function(tl) {
  species_names$name[species_names$sgb_id == tl]
})

tip_labels <- make.unique(tip_labels)
tree$tip.label <- tip_labels
rowTree(tse_s) <- tree

tree$tip.label[str_detect(tree$tip.label, "2154")]
rownames(tse_s)[str_detect(rownames(tse_s), "2154")]


mean(tip_labels %in% rownames(tse_s))
length(rn) - length(rownames(tse_s))
rn[!rn %in% rownames(tse_s)]



colnames(tse_s)
# the subject ids are in the end of the current id names
rx <- "[MC][_F2][2YE][YFE]?[_p1682E][F2M31w]?[EwO28]?_?(\\d{3})_"
id <- str_match(colnames(tse_s), rx)[, 2]
test <- tibble(id)
# for the sampling time
week <- ifelse(str_detect(colnames(tse_s), "E2w"), 2, ifelse(
  str_detect(colnames(tse_s), "E12w"), 12, ifelse(
    str_detect(colnames(tse_s), "Ep32"), 32, ifelse(
      str_detect(colnames(tse_s), "Ep18"), 18, ifelse(
        str_detect(colnames(tse_s), "E6w"), 6, ifelse(
          str_detect(colnames(tse_s), "E8MO"), 32, ifelse(
            str_detect(colnames(tse_s), "2Y"), 104, NA)))))))
# maternal vs infant sample
origin <- ifelse(str_starts(colnames(tse_s), pattern = "M"), "m", "i")
# prenatal vs postnatal
pre <- ifelse(str_detect(colnames(tse_s), "p[31][28]"), TRUE, FALSE)

md <- tibble(
  id,
  origin,
  pre,
  week,
  old_sample_ids = colnames(tse_s)
  ) %>%
    mutate(sid = glue::glue('{id}{origin}{pre}{week}'))


colData(tse_s) <- colData(tse_s) %>%
                  as.data.frame() %>%
                  rownames_to_column("old_sample_ids") %>%
                  full_join(md, by = "old_sample_ids") %>%
                  column_to_rownames("old_sample_ids") %>%
                  DataFrame()

colnames(tse_s) <- colData(tse_s)$sid

# how many samples per time point
colData(tse_s) %>%
  as.data.frame() %>%
  count(origin, pre, week)

# if there are duplicate samples we must disregard one
sids <- colData(tse_s) %>%
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
usid <- colData(tse_s) %>%
  as.data.frame() %>%
  rownames_to_column("sid_tse") %>%
  .$sid_tse

table(colnames(tse_s)) %>% as.data.frame() %>% filter(Freq > 1)
colnames(tse_s) <- usid

rownames(tse_s)
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/mb_helper.R")


tse_g <- mergeFeaturesByRank(tse_s, rank = "Genus")
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
colData(tse_s)$sid <- colnames(tse_s)
dim(tse_s)

tse_s <- tidySummarizedExperiment::filter(tse_s, !sid %in% c("554mpre18.1", "636ipost6.1")) 
colData(tse_s) <- colData(tse_s) %>% 
  as.data.frame() %>%
  mutate(t = ifelse(origin == "i", NA, ifelse(week == 18, "t1", ifelse(pre, "t2", ifelse(!pre, "t3", NA))))) %>%
  DataFrame()
dim(tse_s)
# now the double samples are filtered out and we can proceed creating a tse 

# calculate phylogenetic metrics
tse_s <- estimateFaith(tse_s, assay.type = "counts")

tse_i <- tse_s[, colData(tse_s)$origin == "i"]
tse_m <- tse_s[, colData(tse_s)$origin == "m"]

# combine with metadata (infancy)
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
  count(origin, pre, week)


# combine with metadata (mothers)
pp <- c("post", "pre")
md2 <- read_csv(here::here("data/d.csv")) %>%
  mutate(sample_id = glue::glue('{id}m{pp[1 + pre]}{week}')) %>%
  select(-id, -week, -pre, -t) %>%
  select(sample_id, everything())
head(md2)
dim(md2)

colData(tse_m) <- colData(tse_m) %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  select(-sid) %>%
  left_join(md2, by = "sample_id") %>%
  column_to_rownames("sample_id") %>%
  DataFrame()

colnames(colData(tse_m))
colData(tse_m) %>%
  as.data.frame() %>%
  count(origin, pre, week)
colData(tse_m) %>% colnames()

save(tse_s, file = here::here("data/rdata/tse.Rds"))

# for the mbage model we need one at genus level and prevalence filtered
tse_smiley <- agglomerateByRank(tse_i, rank = "Genus")
tse_smiley <- subsetByPrevalentFeatures(tse_smiley, rank = "Genus", prevalence = 5 / 100)
dim(tse_smiley)
save(tse_smiley, file = here::here("data/rdata/tse_smiley.Rds"))



### functional abundance (2y samples not yet added here)
path_f <- here::here("data/processed_by_turku/joined_pathabundance.tsv")
tse_f <- loadFromHumann(path_f)
af <- assay(tse_f)
rownames(af)
colnames(af)




