source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_utils()

## Load the metadata variable configurations from a spreadsheet configuration.
import.data_loading_manager("metadata_configurator.R")
metadata_config_vars = load_metadata_configuration()
CONTINUOUS_OR_DISCRETE_VARIABLES = metadata_config_vars[["CONTINUOUS_OR_DISCRETE_VARIABLES"]]
DICHOTOMOUS_VARIABLES = metadata_config_vars[["DICHOTOMOUS_VARIABLES"]]
NOMINAL_VARIABLES = metadata_config_vars[["NOMINAL_VARIABLES"]]

## This script handles:
# 1. Grabs the data processing configuration from a Google Sheets CSV
# 2. Merging the metadata from the microbiome metadata (n = 1000) with the external metadata file (n = 3264)
# 3. Converts to numeric variables according to the configuration instructions

######### Recommendations for use:
## Use .MERGED_METADATA if you want access to both the original variables AND the converted variables. It contains both.
##  (*) Scenarios where it's useful is when you're looking at single variables or subsets in univariate fashion and want access to:
##        1/ the original variable or 2/ the converted variable at different times.
##  (*) If the analysis package/function you're using can handle non-numeric variables, then it is be more valuable to use this.
##
## .MERGED_METADATA_NUMERIC contains all the original variables, processed and converted into all possible numeric data types. 
##  (*) Nominal variables are dummy-coded/one-hot encoded. This was an intentional decision. Rather than converting it to numerical categories, 
##    which will take on unintended ordering, dummy-coding makes it easier to interpret. In addition, dummy coding ensures spurious relationships 
##    do not result from the forced ordering that can result from converting natural categories to numeric categories.
##  (*) Date variables are converted from a string representation to a unix time representation, in seconds. This can also be changed to days for 
##    easier readability/interpretability. Unix time is calculated as the time elapsed since Jan 1, 1970 at 00:00:00. For dates beyond this time, 
##    using the default conversion should be ok. But for dates before this time then you will need to handle it differently when used in the models 
##    (will the time be negative?). CHILD study dates are all year 2008+ so no need to be concerned.
##
## Use .MERGED_METADATA_NUMERIC if you want quick access to just the numeric variables. It is subsetted for the numeric versions only.
##  (*) This is useful if you're doing analysis on all the variables at once (multivariate fashion). You will not require any knowledge of what the data looks 
##    like to use it. 
##  (*) It's a one stop shop for running analysis on all variables. 

convert_all_to_numeric_friendly <- function(metadata) {
  ## Converts all columns to a numeric friendly type.
  ## This will result in multiple variations and instances of a single variable. 
  ## Returns the dataframe with multiple extra versions for each variable. 
  
  continuous.subset <- metadata[, colnames(metadata) %in% CONTINUOUS_OR_DISCRETE_VARIABLES]
  binary.subset <- metadata[, colnames(metadata) %in% DICHOTOMOUS_VARIABLES]
  nominal.subset <- metadata[, colnames(metadata) %in% NOMINAL_VARIABLES]
  
  # Handle diseasestatus separately
  metadata$diseasestatus_3y <- ifelse(metadata$diseasestatus_3y == "Asthmatic", 1, 0)
  metadata$diseasestatus_5y <- ifelse(metadata$diseasestatus_5y == "Asthmatic", 1, 0)
  # Convert date string variables to continuous
  metadata <- convert_cols_to_unix_datetime(data = metadata, date_variables = DATE_VARIABLES, suffix = "_continuous", replace = FALSE)
  # Convert the rest of the continuous variables to continuous
  metadata <- convert_cols_to_numeric(metadata, which = CONTINUOUS_OR_DISCRETE_VARIABLES, suffix = "_continuous", replace = FALSE)
  # Convert binary variables to numeric
  metadata <- convert_cols_to_numeric(metadata, which = DICHOTOMOUS_VARIABLES, suffix = "_binary", replace = FALSE)
  # Convert binary variables to factors
  metadata <- convert_cols_to_numeric(metadata, which = DICHOTOMOUS_VARIABLES, suffix = "_fbinary", replace = FALSE, is_binary_factor = TRUE)
  # Convert nominal variables to dummy variables
  metadata <- convert_cols_to_dummy_vars(metadata, which = NOMINAL_VARIABLES, suffix = "_fdnominal", replace = FALSE, as_factor = TRUE)
  metadata <- convert_cols_to_dummy_vars(metadata, which = NOMINAL_VARIABLES, suffix = "_ndnominal", replace = FALSE, as_factor = FALSE)
  # Convert to factors, and add suffix to nominal variables; do not dummy-code them 
  metadata <- convert_cols_to_factors(metadata, which = NOMINAL_VARIABLES, suffix = "_nominal", replace = FALSE)
  return(metadata)
}

ext_metadata <- read.csv(concat_paths(DIR_PATH.PRIMARY_DATA_INPUT, "processed/metadata/data.cli_metadata.added_features.variables_processed.csv"), na.strings = c("", " ", "NA"), row.names = 1)
mb_metadata <- read.csv(concat_paths(DIR_PATH.PRIMARY_DATA_INPUT, "processed/metadata/metadata_filtered_two_timepoints.added_features.csv"), na.strings = c("", " ", "NA"), row.names = 1)

ext_metadata <- data_processor.select_percentage_complete_cols(ext_metadata)

mb_metadata.3m <- subset(mb_metadata, Visit == "3 month")
mb_metadata.1y <- subset(mb_metadata, Visit == "1 year")

ext_metadata.overlapping_subset.3m <- ext_metadata[ext_metadata$SubjectNumber %in% mb_metadata.3m$SubjectNumber, ]
ext_metadata.overlapping_subset.1y <- ext_metadata[ext_metadata$SubjectNumber %in% mb_metadata.1y$SubjectNumber, ]

# Sort them all
mb_metadata.3m <- mb_metadata.3m[order(mb_metadata.3m$SubjectNumber), ]
mb_metadata.1y <- mb_metadata.1y[order(mb_metadata.1y$SubjectNumber), ]
ext_metadata.overlapping_subset.3m <- ext_metadata.overlapping_subset.3m[order(ext_metadata.overlapping_subset.3m$SubjectNumber), ]
ext_metadata.overlapping_subset.1y <- ext_metadata.overlapping_subset.1y[order(ext_metadata.overlapping_subset.1y$SubjectNumber), ]

# Sanity check
message("Sanity check that subject IDs are ordered and in agreement.")
message(paste0("Subject IDs (3m) in disagreement: ", sum(mb_metadata.3m$SubjectNumber != ext_metadata.overlapping_subset.3m$SubjectNumber, na.rm = TRUE)))
message(paste0("Subject IDs (1y) in disagreement: ", sum(mb_metadata.1y$SubjectNumber != ext_metadata.overlapping_subset.1y$SubjectNumber, na.rm = TRUE)))

common_vars <- intersect(colnames(ext_metadata), colnames(mb_metadata))
common_vars <- common_vars[common_vars != "SubjectNumber"]

message("Testing data consistency for 3-month samples:")
for (var in common_vars) {
  num_diff <- sum(mb_metadata.3m[[var]] != ext_metadata.overlapping_subset.3m[[var]], na.rm = TRUE)
  if (num_diff != 0) {
    print(paste0(var, " || number of disagreements: ", as.character(num_diff)))
  }
}
message("Testing data consistency for 1-year samples:")
for (var in common_vars) {
  num_diff <- sum(mb_metadata.1y[[var]] != ext_metadata.overlapping_subset.1y[[var]], na.rm = TRUE)
  if (num_diff != 0) {
    print(paste0(var, " || number of disagreements: ", as.character(num_diff)))
  }
}

ext_metadata.overlapping_subset.3m.col_filtered <- ext_metadata.overlapping_subset.3m[, !(colnames(ext_metadata.overlapping_subset.3m) %in% common_vars)] # Use for debugging
ext_metadata.overlapping_subset.1y.col_filtered <- ext_metadata.overlapping_subset.1y[, !(colnames(ext_metadata.overlapping_subset.1y) %in% common_vars)] # Use for debugging

all_covariates.3m <- join(mb_metadata.3m, ext_metadata.overlapping_subset.3m.col_filtered, by = c("SubjectNumber"))
all_covariates.1y <- join(mb_metadata.1y, ext_metadata.overlapping_subset.1y.col_filtered, by = c("SubjectNumber"))

# Replace all "missing value codes" with NA
missing_val_codes <- c("8888", "999")
all_covariates.3m <- replace_val_with_na(all_covariates.3m, missing_val_codes)
all_covariates.1y <- replace_val_with_na(all_covariates.1y, missing_val_codes)

# Create multiple versions of the existing variables (factor, binary, dummy-coded)
.MERGED_METADATA.3m <- convert_all_to_numeric_friendly(all_covariates.3m)
.MERGED_METADATA.1y <- convert_all_to_numeric_friendly(all_covariates.1y)

# If a variable has a space in it (e.g. because it was converted to a dummy-coded variable, so the class had a space), then replace with underscore.
.MERGED_METADATA.3m <- data_processor.fill_col_name_whitespace(.MERGED_METADATA.3m)
.MERGED_METADATA.1y <- data_processor.fill_col_name_whitespace(.MERGED_METADATA.1y)

# Subset the dataframe to grab metadata variables of the following numeric data types. 
numeric_options <- c("_ndnominal", "_binary", "_continuous")
.MERGED_METADATA_NUMERIC.3m <- data_processor.select_metadata_by_data_processing_type(.MERGED_METADATA.3m, options = numeric_options)
.MERGED_METADATA_NUMERIC.1y <- data_processor.select_metadata_by_data_processing_type(.MERGED_METADATA.1y, options = numeric_options)

# Exclude variables with constant values as downstream models are not designed to handle
.MERGED_METADATA_NUMERIC.3m <- subset_out_constants(.MERGED_METADATA_NUMERIC.3m)
.MERGED_METADATA_NUMERIC.1y <- subset_out_constants(.MERGED_METADATA_NUMERIC.1y)

# Write this to CSV for external use. 
write.csv(x=.MERGED_METADATA_NUMERIC.3m, file=concat_paths(DIR_PATH.PRIMARY_DATA_INPUT, "merged-metadata/3m_numeric.csv"))
write.csv(x=.MERGED_METADATA_NUMERIC.1y, file=concat_paths(DIR_PATH.PRIMARY_DATA_INPUT, "merged-metadata/1y_numeric.csv"))

# Save R objects.
save(.MERGED_METADATA.3m, .MERGED_METADATA.1y, .MERGED_METADATA_NUMERIC.3m, .MERGED_METADATA_NUMERIC.1y, 
     list = c(".MERGED_METADATA.3m", ".MERGED_METADATA.1y", ".MERGED_METADATA_NUMERIC.3m", ".MERGED_METADATA_NUMERIC.1y"), 
     file = concat_paths(DIR_PATH.TMP_DATA_OUTPUT, "processed-data-inputs/merged_metadata_objects.RData"))
