concat_paths <- function(parent_path, ...) {
  ## Concatenates path partials
  to_concat <- c(parent_path)
  paths <- c(...)
  if (length(paths) > 0 || !is.null(paths)) {
    for (item in paths) {
      if (item == "") {
        stop("Empty string provided for path string concatenation.")
      }
    }
    to_concat <- c(parent_path, paths)
  }
  return(paste(to_concat, collapse = "/"))
}

DIR_PATH.PROJECT_PARENT <- "~/Projects/asthma"
DIR_PATH.PROJECT_ROOT <- concat_paths(DIR_PATH.PROJECT_PARENT, "asthma-research")
DIR_PATH.DOCS_HTML <- concat_paths(DIR_PATH.PROJECT_ROOT, "docs/analysis")
DIR_PATH.R_CODE <- concat_paths(DIR_PATH.PROJECT_ROOT, "code/r")
DIR_PATH.R_CODE.SCRIPTS <- concat_paths(DIR_PATH.R_CODE, "scripts")
DIR_PATH.R_CODE.RMD <- concat_paths(DIR_PATH.R_CODE, "rmarkdown")

# Path to the directory containing data output from analysis. It is labeled "tmp" because this data can 
# get replaced between runs.
DIR_PATH.TMP_DATA_OUTPUT <- concat_paths(DIR_PATH.PROJECT_PARENT, "data-objects")

# Path to the directory containing input data for analysis. Some files do change, but only during data
# processing. This folder should not contain any data analysis output files and should strictly be 
# outputted to DIR_PATH.TMP_DATA_OUTPUT.
DIR_PATH.PRIMARY_DATA_INPUT <- concat_paths(DIR_PATH.PROJECT_PARENT, "data/child-study-data-jan-2019")

import.abs_filepaths <- function(path_to_files, purpose_display = "files") {
  ## Imports files given by the vector path_to_files. If a file doesn't exist, an error is thrown and 
  ## import is halted. All files must be valid for this to run successfully.
  for (path_to_file in c(path_to_files)) {
    if (file.exists(path_to_file)) {
      source(path_to_file)
      print(paste("Project script imported: [", path_to_file, "]"))
    } else {
      stop(paste0("Importing ", purpose_display, ": can't import a file that doesn't exist [", path_to_file, "]"))
    } 
  }
}

import.data_loading_manager <- function(filename, file_directory_path = DIR_PATH.R_CODE.SCRIPTS) {
  import.abs_filepaths(concat_paths(DIR_PATH.R_CODE.SCRIPTS, "preprocessed-data-loaders", filename))
}

import.project_helper_script <- function(filename, file_directory_path = DIR_PATH.R_CODE.SCRIPTS) {
  import.abs_filepaths(concat_paths(DIR_PATH.R_CODE.SCRIPTS, filename))
}

import.data_modeling_utils <- function(filename, file_directory_path = DIR_PATH.R_CODE.SCRIPTS) {
  import.abs_filepaths(concat_paths(DIR_PATH.R_CODE.SCRIPTS, "data-modeling-utils", filename))
}

import.feature_selection_utils <- function(filename, file_directory_path = DIR_PATH.R_CODE.SCRIPTS) {
  import.abs_filepaths(concat_paths(DIR_PATH.R_CODE.SCRIPTS, "feature-selection-utils", filename))
}

import.project_utils <- function(filename = "", file_directory_path = DIR_PATH.R_CODE.SCRIPTS, file_purpose = "project utilities") {
  if (is.null(filename) || filename == "") {
    path_to_data_processing_helpers <- concat_paths(DIR_PATH.R_CODE.SCRIPTS, "data_processing_helpers.R")
    path_to_core_utils <- concat_paths(DIR_PATH.R_CODE.SCRIPTS, "core_utilities.R")
    path_to_files <- c(path_to_data_processing_helpers, path_to_core_utils)
  } else {
    path_to_files <- c(concat_paths(file_directory_path, filename))
  }
  import.abs_filepaths(path_to_files, file_purpose)
}

import.microbiome_data_processing <- function() {
  import.data_loading_manager("microbiome_data_processing.R")
}

import.mirkat_modeling_dependencies <- function() {
  ## Imports the essential scripts for modeling with the microbiome kernels utils
  ## Includes loading microbiome data, merging and converting the metadata, and matching samples.
  import.data_loading_manager("merged_metdata_loader.R")
  import.data_loading_manager("microbiome_data_loader.R")
  import.data_modeling_utils("mirkat_utilities.R")
}
