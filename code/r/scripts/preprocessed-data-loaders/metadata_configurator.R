source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_utils()

install_package_if_missing("googlesheets")
library(googlesheets)

load_metadata_configuration <- function(sheet_name = "child_study_data_description_metadata_filtered") {
  # Retrives the data processing configuration from a CSV on the web.
  # Any changes to variable data processing should be handled the CSV config file only.
  
  # Load configuration sheet from Googlesheets
  google_sheet <- gs_title(sheet_name)
  google_sheet.df <- as.data.frame(gs_read(google_sheet))
  variable.descriptions <- google_sheet.df
  
  # Retrieve whitelist items
  rownames(variable.descriptions) <- variable.descriptions$VARIABLE
  variable.descriptions.whitelist <- variable.descriptions[variable.descriptions$STAT_DATA_TYPE != "-", ]
  
  # Grab variables of the default data types
  nominal <- variable.descriptions.whitelist[variable.descriptions.whitelist$STAT_DATA_TYPE == "nominal", ]
  ordinal <- variable.descriptions.whitelist[variable.descriptions.whitelist$STAT_DATA_TYPE == "ordinal", ]
  continuous <- variable.descriptions.whitelist[variable.descriptions.whitelist$STAT_DATA_TYPE == "continuous", ]
  binary <- variable.descriptions.whitelist[variable.descriptions.whitelist$STAT_DATA_TYPE == "binary", ]
  
  CONTINUOUS_OR_DISCRETE_VARIABLES <- rownames(continuous)
  DICHOTOMOUS_VARIABLES <- rownames(binary)
  ORDINAL_VARIABLES <- rownames(ordinal)
  NOMINAL_VARIABLES <- rownames(nominal)
  
  return(list(CONTINUOUS_OR_DISCRETE_VARIABLES = CONTINUOUS_OR_DISCRETE_VARIABLES,
              DICHOTOMOUS_VARIABLES = DICHOTOMOUS_VARIABLES,
              ORDINAL_VARIABLES = ORDINAL_VARIABLES,
              NOMINAL_VARIABLES = NOMINAL_VARIABLES))
}