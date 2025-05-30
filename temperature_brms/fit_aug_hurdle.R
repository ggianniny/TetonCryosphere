# august hurdle

# rescale data by dividing by the global max
# hurdle model for August
# this script just fits the hurdle model for august. 
# need to make a separate script for summarizing output, making figures, etc. 

# script name for running in the terminal
# source("temperature_brms/fit_aug_hurdle.R")


# libraries
library(tidyverse)
library(brms)

# set up directory path and model parameters
write_dir <- "Temperature_brms"
model_month <- 8
full_iter <- 2000


# read in data #### 

source_info <- read.csv("stream_temp_data/source_info.csv")%>%
  rename(site = stream) #rename for merge

# change names of groups to match Khatiwada
source_info <- source_info %>%
  mutate(source = case_when(
    site == "cloudveil" ~ "sub_ice",
    .default = source)) |>
  mutate(source = case_when(
    source == "sub_ice" ~ "rock_glacier",
    source == "snowmelt" ~ "snowfield",
    .default = source))


# read in Temperature data:

temp_clean <- read_csv("stream_temp_data/cleaned_full_datasets/temps_hourly.csv")

# sites to remove
# these sites are from other areas and are not included in the present analysis
rm_sites <- c("silver_run", "death_canyon", "peterson", "quad_cr", "schoolroom", "paintbrush", "windcave")

# modify data to have a month and year column
# remove sites and NA temperature values
temp_clean <- temp_clean |>
  mutate(month = month(date1),
         year = year(date1)) |>
  filter(!is.na(temp_c),
         !site %in% rm_sites)

# join the source info
temp_clean <- left_join(temp_clean,
                        source_info)

# keep only columns for model
# convert year into a z-score
temp <- temp_clean %>% 
  select(temp_c, temp_raw, year, month, source, site) %>%
  mutate(year_s = (year - mean(year)) / sd(year))

# make a table with the mean and sd of year
# this makes back calculating the posterior later easier
fit_data_year_summary <- temp |>
  ungroup() |>
  summarize(mean_year = mean(year), 
            sd_year = sd(year))

# save fit data year summary
saveRDS(fit_data_year_summary, 
        paste(write_dir, 
              "/fit_data_year_summary.rds",
              sep = ""))

# only keep the august data
# scale response variable (temp_c) to be between 0 and 1
# this speeds up model fitting
aug_dat <- temp |>
  filter(month == model_month) |>
  mutate(temp_1 = temp_c / max(temp_c))

# save model data as an RDS object
# This can be used later for plotting
mod_dat <- aug_dat
saveRDS(mod_dat, paste(write_dir, 
                       "/hurdle_fit_data.rds",
                       sep = ""))

# see what priors are needed for the model
get_prior(temp_1 ~ year_s * source +
            (1 + year_s |site) + (1|year_s),
          family = hurdle_gamma(),
          data = aug_dat)



my_prior <- c(
  # linear predictor, b
  #~ 95% of prior between -1 and 1
  # because of the log-link in a gamma model, you exp() these values to get the effect
  # so 95% of this prior allows temperature to change by a factor of exp(-1) = 0.367 to exp(1) = 2.718 from year to year. 
  prior(normal(0, 0.5), class = b), 
  # the response variable has been scaled to be from 0 to 1
  # the intercept is the mean temperature when x = 0
  # since year is centered, x = 0 in the middle of the data
  # We also need to exp() these values to estimate the log temperature value
  # so this distribution has the temperature from exp(-2) * 20 = 2.71c to exp(1) * 20 = 54.37
  prior(normal(-0.75, 0.5), class = Intercept),
  # this prior is for the random site intercept hyperparameter 
  # this allows the effects of site to vary between ~0 and ~3
  # exp(0) = 1, exp(3) = 20
  prior(exponential(2), class = sd))

###
# all other priors are left as the default values
###


hurdle_aug <- brm(
  bf(temp_1 ~ year_s + source + year_s:source + (1 + year_s |site)
     # muting hurdle formula for now
     # few 0's in august, so generic hu ~1 seems ok
    # hu ~ temp_1 ~ year_s + source + year_s:source + (1 + year_s |site) + (1|year_s)
    ),
  family = hurdle_gamma(),
  data = mod_dat,
  prior = my_prior,
  iter = 10)
# threads = 12
# backend = "cmdstanr"

hurdle_aug <- update(hurdle_aug, 
                     iter = full_iter, 
                     chains = 4, 
                     cores = 4, 
                     save_pars = save_pars(all = TRUE),
                     sample_prior = TRUE)
# save the output
saveRDS(hurdle_aug, paste(write_dir,
                    "/fit_model_aug_hurdle.rds", 
                    sep = ""))

