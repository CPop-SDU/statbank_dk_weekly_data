# Download and Danish weekly death and population counts
#
# Jonas Schöley
#
# 2020-09-30
#
# Download data on weekly death counts and quarterly population
# sizes by sex and region from StatBank Denmark.

# Init ------------------------------------------------------------

dat <- list()

# Functions -------------------------------------------------------

#' Download a Table from Statistics Denmark StatBank
#'
#' The default is to request the table at the maximum level
#' of cross-classification among variable levels.
#'
#' @param table_id Single table id.
#' @param total_vars Vector of variable id's for which only totals over
#' all sub-groups should be returned.
#' @param chunk_var Single variable id which is used for chunked API
#' requests. The chunks are defined by the variable values.
#' 
#' @author Jonas Schöley
#' 
#' @note This function connects to the API at http://api.statbank.dk.
#' See https://www.dst.dk/en/Statistik/statistikbanken/api for the API
#' specification.
#' 
#' This function can request a lot of data from the api, therefore, keep
#' its use to a minumum or otherwise risk the IT department at
#' Statistics Denmark to further restrict access to its API.
#' 
#' @return A data frame.
DownloadStatBankTable <- function (
  table_id, total_vars = '', chunk_var = ''
) {
  
  require(danstat)
  require(purrr)
  
  # prepare API request
  meta <- get_table_metadata(table_id)
  variable_ids <- meta$variables$id
  variable_values <- meta$variables$values
  total_loc <- sapply(total_vars, function (x) which(x == variable_ids))
  chunk_values <- ''
  if (chunk_var != '') {
    chunk_values <- variable_values[[which(variable_ids == chunk_var)]][['id']]
  }
  
  # build chunked api requests
  # 1. api request chunk
  # 1.1 variable
  # 1.1$code
  # 1.1$values
  # API request for a given value of a chunk variable
  BuildAPIRequest <-
    function (chunk_value, chunk_var, var_ids, var_values, total_vars) {
      
      api_request <-
        map2(
          var_ids, var_values,
          ~ {
            # unless otherwise specified,
            # request all variable levels
            spec <- list(code = .x, values = NA)
            # if a variable is specified as total,
            # download only totals
            if (.x %in% total_vars) {
              spec <- list(
                code = .x,
                # variable values for totals are
                # always in first position
                values = .y[['id']][1]
              )
            }
            if (.x %in% chunk_var) {
              spec <- list(
                code = .x,
                values = chunk_value
              )
            }
            spec
          }
        )
      
      return(api_request)
      
    }
  
  chunked_api_request <-
    map(
      chunk_values,
      BuildAPIRequest,
      chunk_var = chunk_var,
      var_ids = variable_ids,
      var_values = variable_values,
      total_vars = total_vars
    )
  
  DownloadTableChunks <- function (chunked_api_request, table_id) {
    cat(paste0(unlist(chunked_api_request), collapse = ' '), sep = '\n')
    
    the_data <-
      get_data(
        table_id = table_id,
        variables = chunked_api_request,
        language = 'en'
      )
    
    # give the server some time to breathe
    Sys.sleep(1)
    
    return(the_data)
  }
  
  the_table <-
    map_dfr(
      chunked_api_request,
      DownloadTableChunks,
      table_id = table_id,
      .id = 'chunk'
    )
  
  
  return(the_table)
  
}

# Download StatBank Data ------------------------------------------

dat <- list()

# weekly death counts by sex, 5-year age group, and province
dat$dodc2 <-
  DownloadStatBankTable(
    'DODC2',
    # download chunked over age levels
    chunk_var = 'ALDER'
  )

# quarterly populations numbers by
# sex, 5-year age group, and municipality
dat$folk1b <-
  DownloadStatBankTable(
    'FOLK1B',
    # only select totals for citizenship
    total_vars = 'STATSB',
    chunk_var = 'ALDER'
  )

save(dat, file = 'out/dodc2_and_folk1b.RData')
