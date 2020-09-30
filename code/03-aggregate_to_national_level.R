# Aggregate data to national level
#
# Jonas Schöley
#
# 2020-09-30
#
# Aggregate weekly deaths and population data into national level.

# Init ------------------------------------------------------------

library(tidyverse)

fig <- list()

# Function --------------------------------------------------------

#' Export ggplot
#' 
#' @author Jonas Schöley
ExportFigure <-
  function(figure,
           path,
           filename,
           width = 170,
           height = 100,
           scale = 1,
           device = 'png',
           dpi = 300,
           add_date = FALSE) {
    require(ggplot2)
    
    if (missing(filename)) {
      filename <- tolower(gsub('\\.', '_', make.names(deparse(substitute(figure)))))
    }
    if (isTRUE(add_date)) {
      filename <- paste0(Sys.Date(), '-', filename)
    }
    
    arguments <-
      list(
        filename = paste0(filename, '.', device),
        plot = figure,
        path = path,
        width = width,
        height = height,
        units = "mm",
        scale = scale,
        dpi = dpi,
        device = device
      )
    if (device == 'pdf') {
      arguments$useDingbats <- FALSE 
    }
    
    do.call(ggsave, arguments)
  }

#' Export ggplots Stored in List
#' 
#' @author Jonas Schöley
ExportFiguresFromList <- function(lst, path, ...) {
  figure_names <- tolower(gsub('\\.+', '_', make.names(names(lst))))
  Fun <- function (figure, filename, ...) {
    ExportFigure(figure = figure, filename = filename, ...)
  }
  purrr::pwalk(
    list(lst, figure_names),
    Fun, path = path, ...
  )
}

# Load ------------------------------------------------------------

deaths_population <-
  readRDS('out/2020-09-30-deaths_population.RData')

dk_weekly_deaths_quarterly_population <-
  deaths_population %>%
  group_by(
    sex, age_group, start_of_age_group, width_of_age_group,
    date, year, start_of_epi_year, epi_year, iso_week, epi_week
  ) %>%
  summarise(
    population_quarterly = sum(population_quarterly),
    deaths = sum(deaths)
  ) %>%
  ungroup()

# Visualize aggregated data ---------------------------------------

fig$population <-
  dk_weekly_deaths_quarterly_population %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = population_quarterly, color = sex)) +
  facet_wrap(~age_group, scales = 'free_y') +
  scale_x_date(
    date_minor_breaks = '1 year',
    date_breaks = '3 years',
    date_labels = '%y'
  ) +
  scale_y_continuous(
    labels = scales::number_format(scale = 0.001)
  ) +
  scale_color_manual(values = c(Men = '#0f26bc', Women = '#bc1a0f')) +
  labs(
    y = 'Quarterly population count (thousands)',
    x = 'year'
  ) +
  theme_minimal() +
  theme(legend.position = c(0.8,0.05), legend.title = element_blank())

fig$deaths <-
  dk_weekly_deaths_quarterly_population %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = deaths, color = sex), size = 0.1) +
  facet_wrap(~age_group, scales = 'free_y') +
  scale_x_date(
    date_minor_breaks = '1 year',
    date_breaks = '2 years',
    date_labels = '%y'
  ) +
  scale_color_manual(values = c(Men = '#0f26bc', Women = '#bc1a0f')) +
  labs(
    y = 'Weekly death counts',
    x = 'year'
  ) +
  theme_minimal() +
  theme(legend.position = c(0.8,0.05), legend.title = element_blank())

# Export ----------------------------------------------------------

save(
  dk_weekly_deaths_quarterly_population,
  file = paste0('out/', Sys.Date(), '-dk_weekly_deaths_quarterly_population.RData')
)

ExportFiguresFromList(fig, 'out', add_date = TRUE)