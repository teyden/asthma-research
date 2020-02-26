source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_utils()
import.data_loading_manager("metadata_preprocessed_data_loader.R")

MAIN_RMARKDOWN_RUN_FILE <- concat_paths(DIR_PATH.R_CODE.RMD, "covariate-association", "covariate_logistic_regression.Rmd")
RENDERED_OUTPUT_DIR <- concat_paths(DIR_PATH.DOCS_HTML, "covariate-association", "merged-metadata")

output_filename <- "logistic_regression_models.html"
output_file_path <- concat_paths(RENDERED_OUTPUT_DIR, output_filename)

data <- data_processor.select_metadata_by_data_processing_type(DATA, use_default_options = TRUE)
covariates <- c()
predictors <- colnames(data)
outcomes <- c("diseasestatus_3y_fbinary", "diseasestatus_5y_fbinary")
confounders <- c("StudyCenter_nominal")

rmarkdown::render(rmarkdown_file_path, 
                  params = list(
                    data = data,
                    outcomes = outcomes,
                    predictors = predictors,
                    covariates = covariates,
                    confounders = confounders, 
                    plot_subtitle = "Using full cohort metadata, n = 3264; control for StudyCenter",
                    page_title = "Mass Logistic Regression Models"
                  ),
                  output_file = output_file_path)
