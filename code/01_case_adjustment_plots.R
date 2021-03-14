#' model original and adjusted covid case positivity rates using data from https://covidtracking.com/
#' https://covidtracking.com/ was discontinued 7th March 2021
#' the positivity adjustment formula is courtesy of http://freerangestats.info see following links
#' http://freerangestats.info/blog/2020/05/09/covid-population-incidence
#' http://freerangestats.info/blog/2020/05/17/covid-texas-incidence
#' https://github.com/ellisp/blog-source/blob/master/_working/0178-covid-prevalence-inference.R

# install.packages(c("tidyverse", "janitor", "scales", "Cairo", "mgcv", "EpiEstim", "patchwork", "lubridate", "ggrepl", "ggseas"), dependencies = TRUE)
# devtools::install_github("ellisp/frs-r-package/pkg")
# remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"))

library(here)
library(tidyverse)
library(janitor)
library(scales)
library(Cairo)
library(mgcv)
library(EpiEstim)
library(patchwork)
library(frs)
library(tabulizer)
library(lubridate)
library(ggrepel)
library(ggseas)
# CairoWin()

source(here("code/utils.R"))

#' import data

#' covid case data by state
covid_tracking_state_data_raw = read_csv("https://covidtracking.com/api/v1/states/daily.csv") 
covid_tracking_state_info = read_csv("https://covidtracking.com/api/v1/states/info.csv")

glimpse(covid_tracking_state_data_raw)
glimpse(covid_tracking_state_info)

names(covid_tracking_state_data_raw)
names(covid_tracking_state_info)

covid_tracking_state_data_raw = covid_tracking_state_data_raw %>% rename("state_abb"="state")
covid_tracking_state_info = covid_tracking_state_info %>% rename("state_abb"="state")

unique(covid_tracking_state_data_raw$state_abb)
unique(covid_tracking_state_info$state_abb)

#' population by state
census_2019_state_pop_url = "https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/state/detail/SCPRC-EST2019-18+POP-RES.csv"
census_2019_state_pop = read_csv(census_2019_state_pop_url)

state.name_dc = c(state.name, "District of Columbia")
state.abb_dc = c(state.abb, "DC")

glimpse(census_2019_state_pop)
names(census_2019_state_pop)

census_2019_state_pop = census_2019_state_pop %>% filter(NAME %in% state.name_dc) %>% 
  rename('state_name'='NAME', 'pop_estimate_2019'="POPESTIMATE2019") %>% 
  mutate(state_abb = state.abb_dc[match(state_name, state.name_dc)], state_name) %>%
  select(state_name, state_abb,  pop_estimate_2019) %>% 
  arrange(desc(pop_estimate_2019))

top_12_populous_states = census_2019_state_pop %>% slice(1:12)

covid_tracking_state_data_raw = covid_tracking_state_data_raw %>% left_join(., census_2019_state_pop %>% select(-state_name), by="state_abb")

covid_tracking_state_data = covid_tracking_state_data_raw %>%
  filter(state_abb %in% state.abb_dc) %>%
  mutate(date = as.Date(as.character(date), format = "%Y%m%d")) %>%
  clean_names() %>%
  #' force total number of tests to be at least as many as the number of positives:
  mutate(total_test_results_increase = pmax(positive_increase, total_test_results_increase)) %>%
  mutate(pos_rate = positive_increase / total_test_results_increase) %>%
  arrange(date) %>%
  mutate(date_n = as.numeric(date))  %>%
  left_join(select(covid_tracking_state_info, state_abb, state_name = name), by = "state_abb")

sort(unique(covid_tracking_state_data$state_abb))
length(sort(unique(covid_tracking_state_data$state_abb)))
names(covid_tracking_state_data)

#' just the 12 biggest states
top_12_populous_states_case_data = covid_tracking_state_data %>%
  group_by(state_abb) %>%
  dplyr::summarise(max_pos = max(positive)) %>%
  arrange(desc(max_pos)) %>%
  filter(state_abb %in% top_12_populous_states$state_abb) %>%
  inner_join(covid_tracking_state_data, by = "state_abb") %>%
  #' state has to be a factor for use in mgcv::gam:
  mutate(state_name = fct_reorder(state_name, positive, .fun = sum),
         pos_rate = ifelse(pos_rate<0, 0, pos_rate)) %>%
  #' we want deaths in 7-14 days time as a crude indicator of cases now, for use later
  #' Tried various methods and 7 was best. Obviously, if doing this 'for real', 7 should
  #' be a parameter we estimate from the data
  group_by(state_abb) %>%
  arrange(date) %>%
  mutate(deaths_x_days_later = lead(death_increase, 7)) %>%
  ungroup()

sort(unique(top_12_populous_states_case_data$state_abb))

#' smooth the positive test rates
top_12_populous_states_mod = gam(pos_rate ~ state_name + s(date_n, by = state_name), 
                    data = top_12_populous_states_case_data, 
                    family = quasibinomial,
                    weights = total_test_results_increase)

top_12_populous_states_plot_data = top_12_populous_states_case_data %>%
  generate_smooth_data(., model=top_12_populous_states_mod)

texas_plot_data = top_12_populous_states_case_data %>%
  filter(state_name == "Texas") %>%
  generate_smooth_data(., model=top_12_populous_states_mod)

top_12_populous_states_case_adjustment_plot = top_12_populous_states_plot_data %>%
  ggplot(aes(x = date, y = value, colour = variable)) +
  facet_wrap(~state_name, scale = "free_y")  +
  stat_rollapplyr(index.ref = 60, width = 7) +
  the_theme +
  the_labs +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Trends in daily COVID-19 cases (rolling seven-day average, scale-free index).
           With and without case adjustment for test positivity rate")

#' plot suggests relatively more unknown cases in March and April

texas_case_adjustment_plot = texas_plot_data %>%
  ggplot(aes(x = date, y = value, colour = variable)) +
  stat_rollapplyr(index.ref = 70, width = 7) +
  the_theme +
  the_labs +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("Trends in daily COVID-19 cases in Texas (rolling seven-day average, scale-free index). 
          With and without case adjustment for test positivity rate")

#' after adjustment for test-positivity, plot suggests new cases are still accelerating

ifelse(!dir.exists(file.path(here("output"))), dir.create(file.path(here("output"))), FALSE)
ifelse(!dir.exists(file.path(here("rdata"))), dir.create(file.path(here("rdata"))), FALSE)

ggsave(here("output/top_12_populous_states_case_adjustment_plot.png"),
       top_12_populous_states_case_adjustment_plot, width=10.67, height=6, dpi=120)

ggsave(here("output/texas_case_adjustment_plot.png"),
       texas_case_adjustment_plot, width=10.67, height=6, dpi=120)
