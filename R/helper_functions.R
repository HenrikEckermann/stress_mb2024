library(HDInterval)
library(tidyverse)


# calculate contrasts for brms models
calc_contrast <- function(tb, str_vs_str = c("con0", "con1"), ctname = str_vs_str[2]) {
  group1 <- select(tb, contains(str_vs_str[1])) %>%
              pivot_longer(
              everything(), 
              names_to = "subgroups", 
              values_to = "value"
            ) 
  group2 <- select(tb, contains(str_vs_str[2])) %>%
              pivot_longer(
              everything(), 
              names_to = "subgroups", 
              values_to = "value"
            ) 
  contrast <- group2$value - group1$value 
  tb <- tibble(
      name = ctname,
      beta = median(contrast),
      lower = hdi(contrast)[1],
      upper = hdi(contrast)[2]) %>%
    mutate(across(where(is.numeric), round, 2))
  return(tb)
}


summarise_posterior <- function(model, parameters, n_digits = 3) {
  as_draws_df(model) %>%
    select(all_of(parameters)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    group_by(parameter) %>%
    summarise(
      m = median(value),
      sd = sd(value),
      lower = hdi(value)[1],
      upper = hdi(value)[2],
      p = mean(value > 0)
    ) %>%
    mutate(across(where(is.numeric), round, n_digits))
}

# custom function to summarize interactions
summarise_posterior2 <- function(drawsdf, parameters, n_digits = 3) {
  pivot_longer(drawsdf, everything(), names_to = "parameter", values_to = "value") %>%
    group_by(parameter) %>%
    summarise(
      m = median(value),
      sd = sd(value),
      lower = hdi(value)[1],
      upper = hdi(value)[2],
      p = mean(value > 0)
    ) %>%
    mutate(across(where(is.numeric), round, n_digits))
}