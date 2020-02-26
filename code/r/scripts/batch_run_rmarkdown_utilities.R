source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_utils()

## UTILITY FUNCTIONS
generate_html_filename <- function(filename_base) {
  return(paste0(filename_base, ".html"))
}
generate_rdata_filename <- function(filename_base) {
  return(paste0(filename_base, "_objects.RData"))
}
generate_output_objects_filepath <- function(output_objects_filename) {
  return(concat_paths(DIR_PATH.TMP_DATA_OUTPUT, output_objects_filename))
}

render_logistic_regression_analyis <- function(
  rmarkdown_file_path, 
  data, outcomes, 
  output_dir_path, output_filename, output_objects_filename = NULL, output_objects_filepath = NULL,
  predictors = NULL, covariates = NULL, confounders = NULL, 
  plot_subtitle = "", page_title = "", 
  idx = NULL) {
  ## Generates RMarkdown files automatically summarizing logistic regression analyses. 

  if (is_string_blank(rmarkdown_file_path)) {
    stop("Rmarkdown file path is blank (NULL or empty).")
  }
  if (is_string_blank(output_dir_path) || is_string_blank(output_filename)) {
    stop("Output directory path or output filename is blank (NULL or empty). Cannot create models without an output file to render to.")
  }
  if (is_dataframe_blank(data)) {
    stop("Data container (containing all the variables) cannot be blank.")
  }
  if (is_vector_blank(outcomes)) {
    stop("Outcomes is blank (NULL or empty). Cannot create models with empty outcomes")
  }
  if (is_vector_blank(predictors)) {
    warning("Predictors vector is blank. Using all variables from the <data> dataframe. If predictors are provided then the models will use only the variables from that.")
  }
  
  output_file_path <- file.path(output_dir_path, output_filename)
  
  rmarkdown::render(rmarkdown_file_path, 
                    params = list(
                      data = data,
                      outcomes = outcomes,
                      predictors = predictors,
                      covariates = covariates,
                      confounders = confounders, 
                      plot_subtitle = plot_subtitle,
                      page_title = page_title,
                      output_objects_filename = output_objects_filename,
                      output_objects_filepath = output_objects_filepath
                    ),
                    output_file = output_file_path)
  
  function_name <- "render_logistic_regression_analyis():"
  output_message <- paste0(function_name, " output files rendered to ", output_dir_path)
  if (!is.null(idx)) {
    output_message <- paste0(function_name, " case #", idx, ", output files rendered to ", output_dir_path)
  }
  message(output_message)
}

render_effect_size_comparator_analysis <- function(
  rmarkdown_file_path, 
  rendered_output_dir_path, rendered_output_filename, 
  diff_tables, diff_plots, cov_summary,
  plot_subtitle = "", page_title = "", 
  idx = NULL) {
  ## Generates RMarkdown files automatically summarizing effect size comparison analysis. 
  
  if (is_string_blank(rmarkdown_file_path)) {
    stop("Rmarkdown file path is blank (NULL or empty).")
  }
  if (is_string_blank(rendered_output_dir_path) || is_string_blank(rendered_output_filename)) {
    stop("Output directory path or output filename is blank (NULL or empty).")
  }
  if (is.null(diff_tables)) {
    stop("Data container (containing all the variables) cannot be blank.")
  }
  
  output_file_path <- file.path(rendered_output_dir_path, rendered_output_filename)
  
  rmarkdown::render(rmarkdown_file_path, 
                    params = list(
                      diff_tables = diff_tables,
                      diff_plots = diff_plots,
                      cov_summary = cov_summary,
                      plot_subtitle = plot_subtitle,
                      page_title = page_title
                    ),
                    output_file = output_file_path)
  
  function_name <- "render_effect_size_comparator_analysis():"
  output_message <- paste0(function_name, " output files rendered to ", rendered_output_dir_path)
  if (!is.null(idx)) {
    output_message <- paste0(function_name, " case #", idx, ", output files rendered to ", rendered_output_dir_path)
  }
  message(output_message)
}