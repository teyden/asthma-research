source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_helper_script("data_processing_helpers.R")

install_package_if_missing <- function(package_name = NULL, package_names = NULL) {
  if (is_vector_blank(package_names) && !is_string_blank(package_name)) {
    if (!require(package_name, character.only = TRUE)) install.packages(package_name)
  } else if (!is_vector_blank(package_names) && is_string_blank(package_name)) {
    for (package_name in package_names) {
      if (!require(package_name, character.only = TRUE)) install.packages(package_name)
    }
  } else {
    stop("Provide only one parameter.")
  }
}

print_process_message <- function(curr_file, message) {
  message(paste0(curr_file, ": ", message))
}

## ESSENTIAL VALIDATION UTILITIES
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
is_col_constant <- function(data, colname, remove.na = TRUE) {
  nonblank_subset <- data[!is_column_blank(data, colname), ]
  first_val <- nonblank_subset[1, c(colname)]
  is_constant <- all(data[, c(colname)] == first_val, na.rm = remove.na)
  return(is_constant)
}

is_constant <- function(vect, remove.na = TRUE) {
  nonblank_subset <- vect[!is_vector_el_blank(vect)]
  first_val <- nonblank_subset[1]
  is_constant <- all(vect == first_val, na.rm = remove.na)
  return(is_constant)
}

is_vector_el_blank <- function(vect) {
  ## Checks whether the elements of the vector are blank.
  ## Returns a logical vector.
  return(is.na(vect) | vect == "")
}

any_blank <- function(vect) {
  ## Checks if any values in the vector is blank.
  ## Returns a logical vector.
  return(any(is_vector_el_blank(vect)))
}

all_blank <- function(vect) {
  ## Checks if all values in the vector are blank.
  ## Returns a logical vector.
  return(all(is_vector_el_blank(vect)))
}

is_row_blank <- function(df, rowname, cols_to_ignore = NULL) {
  ## Checks if the row in a dataframe is blank. 
  ## Returns a logical vector.
  if (!is.null(cols_to_ignore) && !all_blank(cols_to_ignore)) {
    row <- df[rowname, !(names(df) %in% cols_to_ignore)]
  } else {
    row <- df[rowname, ]
  }
  return(all_blank(row))
}

get_constant_cols <- function(data, remove.na = TRUE) {
  ## Finds all the columns with constant values.
  ## Return vector of names of columns with constant values.
  constant_variables <- c()
  i = 1
  for (colname in colnames(data)) {
    is_constant <- is_col_constant(data, colname, remove.na = remove.na)
    if (!is.na(is_constant) && is_constant == TRUE) {
      constant_variables[i] <- colname
      i <- i + 1
    }
  } 
  return(constant_variables)
}

# Checks blank condition for a column in a dataframe
is_column_blank <- function(data, colname) {
  return(is.na(data[[colname]]) | data[[colname]] == "")
}

count_na <- function(data, col_name="") {
  if (col_name == "") {
    return(sum(is.na(data)))
  } else {
    return(sum(is.na(data[[col_name]])))
  }
}

is_dataframe_blank <- function(df) {
  if (typeof(df) != "NULL" && typeof(df) != "list") {
    stop("Incorrect object type. Provide list/dataframes only.")
  } else {
    return(is.null(df) || nrow(df) == 0)
  }
}

is_vector_blank <- function(vect) {
  ## Checks whether the input is blank. Does not check it element-wise. Returns a logical.
  ## Note: there isn't an explicit vector type. Character vectors are recognized just as "character" and doubles as "double". It's weird, but just avoiding using "list".
  if (typeof(vect) == "list") {
    stop("Incorrect object type [list]. Provide vectors only.")
  } else {
    return(is.null(vect) || length(vect) == 0)
  }
}

is_string_blank <- function(string) {
  if (typeof(string) != "NULL" && (typeof(string) != "character" || length(string) != 1)) {
    stop("Incorrect object type. Provide characters only. If a vector, then must have a length of 1, or provide the character string directly.")
  } else {
    # FYI: is.null(NULL) == is.null(c(NULL))
    return(is.na(string) || is.null(string) || nchar(string) == 0)
  }
}