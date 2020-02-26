source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_utils()
import.project_helper_script("batch_run_rmarkdown_utilities.R")

MAIN_RMARKDOWN_RUN_FILE <- concat_paths(DIR_PATH.R_CODE.RMD, "covariate-association", "covariate_logistic_regression.Rmd")
RENDERED_OUTPUT_DIR <- concat_paths(DIR_PATH.DOCS_HTML, "covariate-association", "merged-metadata")

##############################################################################################
##############################################################################################
##############################################################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0 && (args[1] == "--s" || args[1] == "--simple")) {
  
  # Handle loading data. 
  if (length(args) > 1 && args[2] == "--skip-data-processing") {
    # Load from previously saved objects.
    load(concat_paths(DIR_PATH.TMP_DATA_OUTPUT, "merged_metadata_generic.RData"))
  } else {
    # Fresh run of data processing.
    import.data_loading_manager("merge_metadata.R")
  }
  
  ## Case 1: simple linear regression, no control for confounding. <disease/phenotype> ~ <single covariate>
  data <- data_processor.select_metadata_by_data_processing_type(.MERGED_METADATA.3m, use_default_options = TRUE)
  outcomes <- data_processor.get_disease_and_phenotype_col_names(data)
  data <- data_processor.set_reference_level_for_response(data, outcomes)
  output_filename_base <- "simple_logistic_regression"

  cases <- list(
    "1" = list(
      data = data,
      outcomes = data_processor.get_disease_and_phenotype_col_names(data), # Might need to use the dnominal atopicstatus variables instead
      predictors = colnames(data[, !(names(data) %in% names(outcomes))]),  # This gets handled by the regression wrapper functions anyway, but good to be explicit.
      output_filename = generate_html_filename(output_filename_base),
      output_objects_filename = generate_rdata_filename(output_filename_base)
    )
  )
  
  render_logistic_regression_analyis(
    rmarkdown_file_path = MAIN_RMARKDOWN_RUN_FILE,
    data = cases$`1`$data,
    outcomes = cases$`1`$outcomes,
    output_dir_path = RENDERED_OUTPUT_DIR,
    output_filename = cases$`1`$output_filename,
    output_objects_filename = cases$`1`$output_objects_filename,
    predictors = cases$`1`$predictors,
    covariates = NULL,
    confounders = NULL,
    plot_subtitle = "No control for confounders",
    page_title = "Simple Logistic Regression Tests"
  )

} else if (length(args) > 0 && (args[1] == "--m" || args[1] == "--multiple")) {
  
  # Handle loading data. 
  if (length(args) > 1 && args[2] == "--skip-data-processing") {
    # Load from previously saved objects.
    load(concat_paths(DIR_PATH.TMP_DATA_OUTPUT, "merged_metadata_generic.RData"))
  } else {
    # Fresh run of data processing.
    import.data_loading_manager("merge_metadata.R")
  }
  
  ## Case 2: multiple linear regression, no control for confounding. <disease/phenotype> ~ <single covariate>
  data <- data_processor.select_metadata_by_data_processing_type(.MERGED_METADATA.3m, use_default_options = TRUE)
  outcomes <- data_processor.get_disease_and_phenotype_col_names(data)
  data <- data_processor.set_reference_level_for_response(data, outcomes)
  output_filename_base <- "multiple_logistic_regression"
  
  cases <- list(
    "1" = list(
      data = data,
      outcomes = data_processor.get_disease_and_phenotype_col_names(data),
      predictors = colnames(data[, !(names(data) %in% names(outcomes))]), 
      output_filename = generate_html_filename(output_filename_base),
      output_objects_filename = generate_rdata_filename(output_filename_base)
    )
  )
  
  render_logistic_regression_analyis(
    rmarkdown_file_path = MAIN_RMARKDOWN_RUN_FILE,
    data = cases$`1`$data,
    outcomes = cases$`1`$outcomes,
    output_dir_path = RENDERED_OUTPUT_DIR,
    output_filename = cases$`1`$output_filename,
    output_objects_filename = cases$`1`$output_objects_filename,
    predictors = cases$`1`$predictors[31:35],
    covariates = cases$`1`$predictors[1:3],
    confounders = NULL,
    plot_subtitle = "No control for confounders. With control for covariates that are statistically significant in simple logistic regresssion tests.",
    page_title = "Multiple Logistic Regression Tests"
  )
}
