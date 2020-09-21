# Download and Danish weekly death and population counts
#
# Jonas Schöley
#
# 2020-09-21
#
# Merge population and death count data and summarize into
# province aggregation level.

# Init ------------------------------------------------------------

library(tidyverse)

dat <- list()
fig <- list()

# Constants -------------------------------------------------------

cnst <- list(
  # this is where we want our time-series to start
  # this is the origin of the "weeks_since_origin" variable
  origin_date = lubridate::as_date('2007-12-31'),
  week_epi_year_starts = 40
)

# Functions -------------------------------------------------------

#' Convert Week of Year to Date
#'
#' @param year Year integer.
#' @param week Week of year integer (1 to 53).
#' @param weekday Weekday integer (1, Monday to 7, Sunday).
#' @param offset Integer offset added to `week` before date calculation.
#'
#' @return A date object.
#' 
#' @source https://en.wikipedia.org/wiki/ISO_8601
#'
#' @author Jonas Schöley
#'
#' @examples
#' # the first Week of 2020 actually starts Monday, December 30th 2019
#' ISOWeekDate2Date(2020, 1, 1)
ISOWeekDate2Date <- function (year, week, weekday = 1, offset = 0) {
  require(ISOweek)
  isoweek_string <-
    paste0(
      year, '-W',
      formatC(
        week+offset,
        flag = '0',
        format = 'd',
        digits = 1
      ),
      '-', weekday
    )
  ISOweek2date(isoweek_string)
}

#' Calculate Weeks Since Some Origin Date
#'
#' @param date Date string.
#' @param origin_date Date string.
#' @param week_format Either 'integer' for completed weeks or
#' 'fractional' for completed fractional weeks.
#'
#' @return Time difference in weeks.
#'
#' @author Jonas Schöley
#'
#' @examples
#' # My age in completed weeks
#' WeeksSinceOrigin(Sys.Date(), '1987-07-03')
WeeksSinceOrigin <-
  function(date, origin_date, week_format = "integer") {
    require(ISOweek)
    fractional_weeks_since_origin <-
      as.double(difftime(
        as.Date(date),
        as.Date(origin_date),
        units = "weeks"
      ))
    switch(
      week_format,
      fractional = fractional_weeks_since_origin,
      integer = as.integer(fractional_weeks_since_origin)
    )
  }

# Codebook --------------------------------------------------------

codebook <- list()

codebook$age_group <-
  tribble(
    ~code,                ~label,     ~starting_age, ~width,
    '0-4 years',          '[00,05)',   0,             5,
    '5-9 years',          '[05,10)',   5,             5,
    '10-14 years',        '[10,15)',  10,             5,
    '15-19 years',        '[15,20)',  15,             5,
    '20-24 years',        '[20,25)',  20,             5,
    '25-29 years',        '[25,30)',  25,             5,
    '30-34 years',        '[30,35)',  30,             5,
    '35-39 years',        '[35,40)',  35,             5,
    '40-44 years',        '[40,45)',  40,             5,
    '45-49 years',        '[45,50)',  45,             5,
    '50-54 years',        '[50,55)',  50,             5,
    '55-59 years',        '[55,60)',  55,             5,
    '60-64 years',        '[60,65)',  60,             5,
    '65-69 years',        '[65,70)',  65,             5,
    '70-74 years',        '[70,75)',  70,             5,
    '75-79 years',        '[75,80)',  75,             5,
    '80-84 years',        '[80,85)',  80,             5,
    '85-89 years',        '[85,90)',  85,             5,
    '90-94 years',        '[90,95)',  90,             5,
    '95-99 years',        '[95,100)', 95,             5,
    '100 years and over', '100+',     100,            NA
  )


# lookup-table region, province, municipality
codebook$municipality <-
  read.csv(
    'code/dk_region_province_municipality.csv',
    fileEncoding = 'UTF-8'
  )

# Load raw data ---------------------------------------------------

load('data/raw/dodc2_and_folk1b')

# Prepare weekly death counts -------------------------------------

dat$deaths <-
  dat$dodc2 %>%
  # standardize variable names
  rename(
    province = 'OMRÅDE',
    sex = 'KØN',
    age_group = 'ALDER',
    year_week = 'TID',
    deaths = 'INDHOLD'
  ) %>%
  # separate year and week of year variables
  separate(
    year_week, into = c('year', 'week'), sep = 'U'
  ) %>%
  # subset to data of interest
  filter(
    sex != 'Total',
    age_group != 'Total',
    province != 'All Denmark'
  ) %>%
  # convert date variables to integer
  mutate_at(
    vars(year, week, deaths),
    as.integer
  ) %>%
  mutate(
    # age group as factor variable
    age_group =
      factor(
        age_group,
        levels = codebook$age_group$code,
        labels = codebook$age_group$label
      ),
    # convert ISO8601 week format to standard date format
    # week starts at 01 and contains 53 in case of leap-year
    date = ISOWeekDate2Date(year, week)
  ) %>%
  arrange(
    year, week
  ) %>%
  ungroup() %>%
  select(
    province, sex, age_group,
    year, iso_week = week, date,
    deaths
  )

# Prepare quarterly population counts -----------------------------

dat$population <-
  dat$folk1b %>%
  # standardize variable names
  rename(
    municipality = 'OMRÅDE',
    sex = 'KØN',
    age_group = 'ALDER',
    year_quarter = 'TID',
    population = 'INDHOLD'
  ) %>%
  # subset to data of interest
  filter(
    sex != 'Total',
    age_group != 'Total',
    municipality != 'All Denmark'
  ) %>%
  # separate year and quarter variables
  separate(
    year_quarter, into = c('year', 'quarter'), sep = 'Q'
  ) %>%
  # convert date variables to integer
  mutate_at(
    vars(year, quarter, population),
    as.integer
  ) %>%
  mutate(
    # age group as factor variable
    age_group =
      factor(
        age_group,
        levels = codebook$age_group$code,
        labels = codebook$age_group$label
      ),
    # convert year-quarter into date assuming quarters start
    # at week 1, 14, 27, 40 
    iso_week =
      quarter*13-13+1,
    date =
      ISOWeekDate2Date(
        year,
        week = iso_week
      ),
    # which province does this municipality belong to?
    province =
      factor(
        municipality,
        codebook$municipality$municipality,
        codebook$municipality$province
      )
  ) %>%
  # summarise population counts over municipality within each province
  group_by(
    year, iso_week, date, sex, province, age_group
  ) %>%
  summarise(
    population_quarterly = sum(population)
  ) %>%
  ungroup() %>%
  select(
    province, sex, age_group,
    year, iso_week, date,
    population_quarterly
  )

# Merge deaths and population -------------------------------------

# merge deaths and population data
# careful: this join may not work on Windows due
# to bad support for UTF-8. Check join for provinces
# with special Danish characters.
dat$deaths_population <-
  left_join(
    dat$deaths,
    dat$population
  ) %>%
  arrange(
    province, sex, age_group,
    year, date
  ) %>%
  filter(
    date >= cnst$origin_date
  ) %>%
  # fill in values between quarters
  group_by(
    year, sex, province, age_group
  ) %>%
  arrange(date) %>%
  fill(population_quarterly, .direction = 'down') %>%
  ungroup() %>%
  arrange(province, sex, age_group, date)

# Add additional variables of interest ----------------------------

dat$deaths_population <-
  dat$deaths_population %>%
  mutate(
    # date that corresponding epi-year started
    start_of_epi_year =
      ISOWeekDate2Date(
        ifelse(iso_week<cnst$week_epi_year_starts, year-1, year),
        cnst$week_epi_year_starts, 1
      ),
    # current epi-year
    epi_year =
      paste0(lubridate::year(start_of_epi_year),'/',
             lubridate::year(start_of_epi_year)+1),
    # weeks into epi-year (starting at 0)
    epi_week =
      WeeksSinceOrigin(date, start_of_epi_year),
    # age group start and width
    start_of_age_group =
      factor(
        age_group,
        levels = codebook$age_group$label,
        labels = codebook$age$starting_age
      ) %>%
      as.character() %>%
      as.numeric(),
    width_of_age_group =
      factor(
        age_group,
        levels = codebook$age$label,
        labels = codebook$age$width
      )
  ) %>%
  select(
    province, sex, age_group, start_of_age_group, width_of_age_group,
    date, year, start_of_epi_year, epi_year, iso_week, epi_week,
    deaths, population_quarterly
  )

# Visualize merged data -------------------------------------------

fig$population <-
  dat$deaths_population %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = population_quarterly, color = sex), size = 0.1) +
  facet_grid(age_group~province) +
  scale_x_date(
    date_minor_breaks = '1 year',
    date_breaks = '2 years',
    date_labels = '%y'
  ) +
  scale_y_continuous(
    'Quarterly population count (thousands)',
    labels = scales::number_format(scale = 0.001)
  )

fig$deaths <-
  dat$deaths_population %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = deaths, color = sex), size = 0.1) +
  facet_grid(age_group~province) +
  scale_x_date(
    date_minor_breaks = '1 year',
    date_breaks = '2 years',
    date_labels = '%y'
  ) +
  scale_y_continuous(
    'Weekly death count'
  )

# Export ----------------------------------------------------------

saveRDS(dat$deaths_population, 'out/deaths_population.RData')
