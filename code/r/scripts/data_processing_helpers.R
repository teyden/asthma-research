source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_helper_script("core_utilities.R")

# TODO - consider creating a "commons" package for all these helpers for easier import
# NOTE - naming convention: underscore, to distinguish from R core utilities

# DATAFRAME MANAGEMENT (VALUE AND COLUMN NAME CONVERTERS, DATAFRAME SUBSETTERS)
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
replace_val_with_na <- function(df, values) {
  ## Replaces the provided values with NA. 
  ## Returns the updated dataframe. 
  for (value in values) {
    df[df == value] <- NA
  }
  return(df)
}

round_numeric_columns <- function(df, sig_digs) {
  ## Rounds all the numeric columns to the precision of the significant digits provided.
  ## Returns the updated dataframe. 
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[, nums] <- round(df[, nums], digits = sig_digs)
  return(df)
}

convert_missing_to_median <- function(data) {
  ## Performs imputation of missing values to be replaced with the median for the column.
  ## Returns a new data frame. 
  imputed <- data.table::copy(data)
  names <- c()
  i <- 1
  for(colname in colnames(imputed)){
    na_rows <- is.na(imputed[, c(colname)])
    if (sum(na_rows) > 0) {
      imputed[na_rows, c(colname)] <- median(imputed[, c(colname)], na.rm = TRUE)
      names[i] <- colname
      i <- i + 1
    }
  }
  message("Columns modified", names)
  return(imputed)
}

convert_cols_to_numeric <- function(data, which = "all", suffix = "_num", replace = FALSE, is_binary_factor = FALSE) {
  ## Convert categories to numeric. If it's a numerical value already, it'll preserve the number value.
  ## Returns the updated dataframe. 
  vars_to_convert <- which
  col_names <- colnames(data)
  if (length(which) == 1 && which == "all") {
    vars_to_convert <- col_names
  }
  for (colname in vars_to_convert) {
    if (!colname %in% col_names) next 
    name <- paste0(colname, suffix)
    if (replace == TRUE) {
      name <- colname
    }
    if (is_binary_factor == TRUE) {
      col <- data[!is.na(data[[colname]]), ][[colname]]
      if (length(unique(col)) > 2) {
        stop(paste0("Incorrect conversion. Attempt to convert a variable to a binary factor when there are more than 2 values found. Variable = ", colname))
      } else {
        data[[name]] <- as.factor(data[[colname]])
      }
    } else {
      data[[name]] <- as.numeric(data[[colname]])
    }
  }
  return(data)
}

convert_cols_to_factors <- function(data, which = "all", suffix = "_nominal", replace = FALSE) {
  ## Convert to factors. If it's a numerical value already, it'll preserve the number value.
  ## Returns the updated dataframe. 
  vars_to_convert <- which
  col_names <- colnames(data)
  if (length(which) == 1 && which == "all") {
    vars_to_convert <- col_names
  }
  for (colname in vars_to_convert) {
    if (!colname %in% col_names) next 
    name <- paste0(colname, suffix)
    if (replace == TRUE) {
      name <- colname
    }
    data[[name]] <- as.factor(data[[colname]])
  }
  return(data)
}

append_suffix_to_colname <- function(data, which, suffix, replace = FALSE) {
  ## Appends a suffix to all given column names in a dataframe.
  ## Returns the updated dataframe. 
  if (length(which) == 0) return(data)
  col_names <- colnames(data)
  for (colname in which) {
    if (!colname %in% col_names) next 
    col <- data[[colname]]
    name <- paste0(colname, suffix)
    if (replace == TRUE) {
      data[[colname]] <- NULL
    }
    data[[name]] <- col
  }
  return(data)
}

convert_cols_to_unix_datetime <- function(data, date_variables, date_format = "%Y-%m-%d", suffix = "_num", replace = FALSE, unit = "days") {
  ## Given a vector of date variables inside of the data dataframe, it converts all those dates to the unix time representation
  ## Purpose is to convert dates to a continuous variable that can be used for statistical tests and models.
  ## Returns the updated dataframe. 
  col_names <- colnames(data)
  for (date_variable in date_variables) {
    if (!date_variable %in% col_names) next 
    else {
      name <- paste0(date_variable, suffix)
      if (replace == TRUE) {
        name <- date_variable
      }
      if (unit == "days") {
        data[[name]] <- round(as.numeric(as.POSIXct(data[[date_variable]], format = date_format)) / (1000*60*60*24))
      } else {
        data[[name]] <- as.numeric(as.POSIXct(data[[date_variable]], format = date_format))
      }
    }
  }
  return(data)
}

convert_cols_to_dummy_vars <- function(data, which = "all", suffix = "_fdnominal", replace = FALSE, as_factor = TRUE) {
  ## Create dummy variables from a set of variables.
  ## Returns the updated dataframe. 
  install_package_if_missing("psych")
  library(psych)
  
  vars_to_dummy_code <- which
  col_names <- colnames(data)
  if (length(which) == 1 && which == "all") {
    vars_to_dummy_code <- col_names
  }
  for (orig_var in vars_to_dummy_code) {
    if (!orig_var %in% col_names) next
    dummy_vars <- as.data.frame(dummy.code(data[[orig_var]]))
    for (dummy_var in colnames(dummy_vars)) {
      if (!as_factor) {
        data[[paste0(orig_var, "___", dummy_var, suffix)]] <- dummy_vars[[dummy_var]]
      } else {
        data[[paste0(orig_var, "___", dummy_var, suffix)]] <- as.factor(dummy_vars[[dummy_var]])
      }
    }
    if (replace == TRUE) {
      data[[orig_var]] <- NULL
    }
  }
  return(data)
}

data_processor.convert_cols_to_dummy_vars <- function(data, which = "all", suffix = "_fdnominal", replace = FALSE) {
  return(convert_cols_to_dummy_vars(data, which = "all", suffix = "_fdnominal", replace = FALSE))
}

subset_out_constants <- function(df) {
  ## Removes the columns with constant value. 
  ## Returns the dataframe subset.
  constant_variables <- get_constant_cols(df)
  # Variables that are filtered out due to having constant values
  constant_variables
  df <- df[, !(names(df) %in% constant_variables)]
  return(df)
}

data_processor.remove_blank_rows <- function(df, cols_to_ignore = NULL) {
  ## Removes the blank rows (100% NA) from the dataframe. Ignores the columns provided.
  ## Returns the dataframe subset.
  if (is.null(cols_to_ignore)) {
    # Faster implementation if no columns to ignore
    return(df[rowSums(is.na(df) == 1), ])
  } else {
    rnames = rownames(df)
    for (rowname in rnames) {
      if (is_row_blank(df, rowname, cols_to_ignore)) {
        df <- df[which(rownames(df) != rowname), ]
      }
    }
    return(df)
  }
}

data_processor.select_complete_rows <- function(df) {
  ## Selects rows that are 100% complete (no NA's at all). 
  ## Returns the dataframe subset.
  return(df[rowSums(is.na(df)) == 0, ])
}

data_processor.select_percentage_complete_cols <- function(df, percent_complete_threshold = 70) {
  ## Selects columns with a minimum percentage complete (non-NA) level. Drops the columns from the dataframe that do not meet the threshold. 
  ## Returns the dataframe subset meeting the above requirement.
  N <- dim(df)[1]
  message(paste0("Removing columns with high percentage of blanks/NA's (minimum threshold = ", percent_complete_threshold, "%):"))
  for (col in colnames(df)) {
    total_avail <- N - count_na(df, col)
    percent_complete <- (total_avail / N)*100
    if (percent_complete < percent_complete_threshold) {
      df[[col]] <- NULL
      print(paste0(round(percent_complete), "% complete - ", col))
    }
  }
  return(df)
}

data_processor.fill_col_name_whitespace <- function(data, variable_prefixes = NULL) {
  ## Fills the whitespace in the dataframe column names.
  ## Columns with any prefixes given by variable_prefixes will have whitespace filled, otherwise all columns will be handled. 
  ## Returns the updated dataframe.
  if (is_vector_blank(variable_prefixes)) {
    target_columns <- colnames(data)
  } else {
    target_columns <- data_processor.get_col_names_by_substring(data, variable_prefixes)
  }
  for (orig_var in target_columns) {
    new_name <- gsub(" ", "_", orig_var)
    if (orig_var == new_name) {
      next
    }
    data[[new_name]] <- data[[orig_var]]
    data[[orig_var]] <- NULL
  }
  return(data)
}

data_processor.get_col_names_by_substring <- function(data, substrings) {
  target_col_names <- c()
  for (substring in substrings) {
    target_col_names <- c(target_col_names, names(data[, (grepl(substring, names(data)))]))
  }
  return(target_col_names)
}

data_processor.get_default_data_processing_types <- function() {
  ## Returns a vector of the prefixes for the default data types for the standard modelling methods. 
  return(c("_continuous", "_fbinary", "_nominal"))
}

data_processor.get_valid_data_processing_types <- function() {
  ## Returns a vector of the prefixes for all possible data types existent in the data processing utils.
  return(c("_continuous", "_binary", "_fbinary", "_nominal", "_ndnominal", "_fdnominal"))
}

data_processor.get_disease_and_phenotype_col_names <- function(data, variable_prefixes = c("diseasestatus", "atopicstatus")) {
  return(data_processor.get_col_names_by_substring(data, variable_prefixes))
}

data_processor.set_reference_level_for_response <- function(data, outcomes) {
  for (outcome in outcomes) {
    if (grepl("atopicstatus", outcome)) {
      data[[outcome]] <- relevel(data[[outcome]], ref="Control")
    }
    # For now, no need to bother with diseasestatus because it's binary (0, 1) in which 0 (healthy) is automatically the reference.
  }
  return(data)
}

data_processor.select_metadata_by_data_processing_type <- function(data, use_default_options = TRUE, options = NULL) {
  # Subsets the data matrix by columns, identified by a suffix that indicates how it was processed
  
  # The default option is just the original variables with the suffix appended, and converted to the correct R object type, 
  # without modifying the actual values, feature-engineering, or dummy-coding.
  
  # NOTE: the default options are suitable for the glm() functions. Has mixed data types (characters, int, double, etc).
  # dnominal = dummy coded nominal variables. e.g. one nominal variable with 4 categories will become 4 variables.
  # fbinary = factor coded/converted binary variables; this is so that they are not mistakenly handled as continuous variables
  
  valid_options <- data_processor.get_valid_data_processing_types()
  default_options <- data_processor.get_default_data_processing_types()
  if (use_default_options && is_vector_blank(options)) {
    options <- default_options
  } else if (use_default_options && !is_vector_blank(options)) {
    message("Ignoring default options as alternative options are provided.")
  } else if (!use_default_options && is_vector_blank(options)) {
    stop("Must provide valid options if not using the default options.")
  }
  
  is_valid = TRUE
  blacklist <- data
  for (option in options) {
    if (!(option %in% valid_options)) {
      message(paste0("Processing data type is not valid: [", option, "]. Use one of the valid options:"))
      print(valid_options)
      is_valid = FALSE
      break
    }
    # Deselect all the columns that are undesired until complete
    blacklist <- blacklist[, !(grepl(option, names(blacklist)))]
  }
  # Subset the original data matrix, by selecting all columns that are not in the blacklist subset
  output <- data[, !(names(data) %in% names(blacklist))]
  if (is_valid == TRUE) {
    return(output)
  } else {
    stop("Invalid options provided. Will not produce output unless all are valid.")
  }
}

data_processor.get_colnames_by_suffix <- function(df, suffix) {
  return(colnames(df[, grepl(suffix, colnames(df))]))
}
data_processor.get_continuous_colnames <- function(df, suffix = "_continuous") {
  return(data_processor.get_colnames_by_suffix(df, suffix))
}
data_processor.get_binary_colnames <- function(df, suffix = "_binary") {
  return(data_processor.get_colnames_by_suffix(df, suffix))
}
data_processor.get_nominal_colnames <- function(df, suffix = "_nominal") {
  return(data_processor.get_colnames_by_suffix(df, suffix))
}
data_processor.get_fdnominal_colnames <- function(df, suffix = "_fdnominal") {
  return(data_processor.get_colnames_by_suffix(df, suffix))
}
data_processor.get_ordinal_colnames <- function(df, suffix = "_ordinal") {
  return(data_processor.get_colnames_by_suffix(df, suffix))
}

data_processor.get_name_roots <- function(var_names, suffix) {
  ## Finds the root name of a variable given a suffix. 
  ## Returns the root names.
  output <- data.frame(var = 1:length(var_names))
  i <- 1
  for (var_name in var_names) {
    root <- gsub(suffix, "", var_name)
    output[i, c("var")] <- root
    i <- i + 1
  }
  return(output$var)
}

data_processor.find_variable_alternate_names <- function(full_df, curr_suffix, alt_suffix, variables_of_interest = NULL) {
  ## Given a suffix for a variable, finds all the other equivalent variables (based on variable root names) for an alternative sufix. 
  ## Returns the alternative variable names. 
  
  ## Example: variables present in the dataframe: gender__fbinary and gender__ndnominal
  ## Provided: curr_suffix = "__fbinary", alt_suffix = "__ndnominal"
  ## Function finds all variables with curr_suffix and finds all variable root name matches with alt_suffix.
  ## The function will return "gender__ndnominal"
  
  if (is_vector_blank(variables_of_interest)) {
    curr_names <- data_processor.get_colnames_by_suffix(full_df, curr_suffix)
  } else {
    curr_names <- variables_of_interest
  }
  name_roots <- data_processor.get_name_roots(curr_names, curr_suffix)
  alt_names <- data_processor.get_col_names_by_substring(full_df, name_roots)
  alt_names_subset <- full_df[, colnames(full_df) %in% alt_names]
  alt_names <- data_processor.get_col_names_by_substring(alt_names_subset, alt_suffix)
  return(alt_names)
}