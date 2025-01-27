## function from christian to get population standard rates ##

load_std_population_sizes <- function() {
  std_popsizes <- tibble::tribble(
    ~age_group,    ~std_popsize,
    "00 years",    13818,
    "01-04 years", 55317,
    "05-09 years", 72533,
    "10-14 years", 73032,
    "15-19 years", 72169,
    "20-24 years", 66478,
    "25-29 years", 64529,
    "30-34 years", 71044,
    "35-39 years", 80762,
    "40-44 years", 81851,
    "45-49 years", 72118,
    "50-54 years", 62716,
    "55-59 years", 48454,
    "60-64 years", 38793,
    "65-69 years", 34264,
    "70-74 years", 31773,
    "75-79 years", 26999,
    "80-84 years", 17842,
    "85+ years", 15508
  )
  
  std_popsizes_aggregated <- std_popsizes %>%
    mutate(
      age_group = case_when(
        age_group == "00 years"    ~ "0-24",
        age_group == "01-04 years" ~ "0-24",
        age_group == "05-09 years" ~ "0-24",
        age_group == "10-14 years" ~ "0-24",
        age_group == "15-19 years" ~ "0-24",
        age_group == "20-24 years" ~ "0-24",
        age_group == "25-29 years" ~ "25-44",
        age_group == "30-34 years" ~ "25-44",
        age_group == "35-39 years" ~ "25-44",
        age_group == "40-44 years" ~ "25-44",
        age_group == "45-49 years" ~ "45-64",
        age_group == "50-54 years" ~ "45-64",
        age_group == "55-59 years" ~ "45-64",
        age_group == "60-64 years" ~ "45-64",
        age_group == "65-69 years" ~ "65-74",
        age_group == "70-74 years" ~ "65-74",
        age_group == "75-79 years" ~ "75+",
        age_group == "80-84 years" ~ "75+",
        age_group == "85+ years"   ~ "75+")) %>%
    group_by(age_group) %>%
    summarize(std_popsize = sum(std_popsize)) %>%
    ungroup()
  
  std_popsizes_aggregated %<>% mutate(std_popsize_proportion = std_popsize / sum(std_popsize))
  
  # sourced from table 1
  #
  # Provisional COVID-19 Age-Adjusted Death Rates, by Race and Ethnicity —
  # United States, 2020–2021
  #
  # https://www.cdc.gov/mmwr/volumes/71/wr/mm7117e2.htm
  # April 29, 2022 / 71(17);601-605
  covid19_deaths_by_age <-
    tibble::tribble(
      ~age_group, ~popsize_estimate, ~covid19_deaths,
      "0-24",         102849110,            1595,
      "25-44",          88205838,           21550,
      "45-64",          82769810,          108838,
      "65-74",          32549398,          101408,
      "75+",          23109967,          178072
    )
  
  # calculate rate per capita
  covid19_deaths_by_age %<>% mutate(rate = covid19_deaths / popsize_estimate)
  
  # join in covid19 2021 provisional rate data
  std_popsizes_aggregated %<>% left_join(covid19_deaths_by_age %>% select(age_group, rate))
  
  # calculate expected deaths for the US 2000 million
  std_popsizes_aggregated %<>% mutate(expected_deaths = rate * std_popsize)
  
  return(std_popsizes_aggregated)
}
