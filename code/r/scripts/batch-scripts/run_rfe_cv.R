source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_utils()

DIR_PATH.DEFAULT_RESULTS_OUTPUT = concat_paths(DIR_PATH.TMP_DATA_OUTPUT, "mirkc-results")

print("#############################################################################")
print("#############################################################################")
print(paste0("RFE operation starting: ", timestamp(prefix="", suffix="")))
print("#############################################################################")
print("#############################################################################")

# Handle arguments
library("optparse")
option_list = list(
  make_option(c("--id"), type="character", default=NULL, help="descriptive ID suffix to append to the output folder", metavar="character"),
  make_option(c("--timepoint"), type="character", default=NULL, help="data timepoint", metavar="character"),
  make_option(c("--outcome"), type="character", default=NULL, help="y outcome variable", metavar="character"),
  make_option(c("--kernel"), type="character", default=NULL, help="kernel function", metavar="character"),
  make_option(c("--output"), type="character", default=DIR_PATH.DEFAULT_RESULTS_OUTPUT, help="output directory path", metavar="character"),
  make_option(c("--selectioncriterion"), type="character", default="pval", help="model selection criterion for feature selection", metavar="character"),
  make_option(c("--usetiny"), type="character", default=NULL, help="whether or not to use a condensed number of features (for testing)", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message("Options provided:")
print(opt)

IN_MOCK_STATE = FALSE
NUM_OPT_DEFAULTS = 3 ## Manual
if (length(opt) == NUM_OPT_DEFAULTS) {
  # Set up a mock version for testing in the IDE.
  IN_MOCK_STATE = TRUE
  opt.id_suffix = "this-is-a-test"
  opt.timepoint = "3m"
  opt.outcome = "diseasestatus_3y_binary"
  opt.kernel_func = "bray"
  
  opt.id_suffix = if (is.null(opt$id)) "" else opt$id
  opt.dir_path.results_output = if (is.null(opt$output)) "" else opt$output
  opt.selection_criterion = opt$selectioncriterion
  opt.use_tiny = "yep"
  
  message("Running RFE variable ranking mock test.")
} else if (!is.null(opt$timepoint) && !is.null(opt$outcome) && !is.null(opt$kernel)) {
  # Required arguments
  opt.timepoint = opt$timepoint
  opt.outcome = opt$outcome
  opt.kernel_func = opt$kernel
  
  # Optional arguments
  opt.id_suffix = if (is.null(opt$id)) "" else opt$id
  opt.dir_path.results_output = if (is.null(opt$output)) "" else opt$output
  opt.selection_criterion = opt$selectioncriterion
  opt.use_tiny = opt$usetiny
} else {
  stop("Required arguments must be supplied [--timepoint, --outcome, --kernel]", call.=TRUE)
}

# Imports
import.project_helper_script("split_phyloseq_train_test.R")
import.data_modeling_utils("mirkat_utilities.R")
import.feature_selection_utils("rfe_variable_ranking.R")
import.mirkat_modeling_dependencies()

# Handle creating and setting the directory path 
dir_path.results_output = opt.dir_path.results_output
if (is_string_blank(dir_path.results_output)) { ## Should never happen?
  dir_path.results_output = DIR_PATH.DEFAULT_RESULTS_OUTPUT
}
if (!dir.exists(dir_path.results_output)) {
  dir.create(dir_path.results_output, recursive = TRUE)
  message("New output directory created:")
  message(dir_path.results_output)
}

## TODO - add "id suffix" handling. 

# Handle creating task IDs
file_path.task_IDs = concat_paths(dir_path.results_output, "task_IDs.csv")
if (!file.exists(file_path.task_IDs)) {
  CURRENT_TASK_ID = "task_ID__1"
  task_IDs = data.frame(ID = c(CURRENT_TASK_ID), created = c(timestamp(prefix="", suffix="")))
  write.csv(task_IDs, file_path.task_IDs, quote = FALSE, row.names = FALSE)
} else {
  task_IDs = read.csv(file_path.task_IDs, stringsAsFactors = FALSE)
  
  # Grab the last ID and increment it.
  last_row = task_IDs[nrow(task_IDs), ]
  new_ID = paste0("task_ID__", as.numeric(strsplit(last_row$ID, "__")[[1]][2]) + 1)
  line = paste(c(new_ID, timestamp(prefix="", suffix="")), collapse = ",")
  write(line, file = file_path.task_IDs, append = TRUE)
  
  CURRENT_TASK_ID = new_ID
}

# Create a directory for this unique task ID
dir_path.results_output.CURRENT_TASK_ID = concat_paths(dir_path.results_output, CURRENT_TASK_ID)
dir.create(dir_path.results_output.CURRENT_TASK_ID, recursive = TRUE)

# Set the data variables
if (opt.timepoint == "3m") {
  pso = pso.rarefied.ra.3m
  metadata = merged_metadata.3m
  timepoint = opt.timepoint
} else if (opt.timepoint == "1y") {
  pso = pso.rarefied.ra.1y
  metadata = merged_metadata.1y
  timepoint = opt.timepoint
} else {
  stop(paste0("Timepoint provided is not valid [", opt.timepoint, "]"))
}
outcome = opt.outcome
if (!(outcome) %in% colnames(metadata)) {
  stop(paste0("Outcome (y) variable provided is not found in the metadata [", outcome, "]"))
}
kernel_function = opt.kernel_func
if (!(kernel_function %in% names(ENUM_KERNEL_FUNCTION()))) {
  stop(paste0("Kernel function provided is not valid [", kernel_function, "]"))
}
selection_criterion = opt.selection_criterion
if (!(selection_criterion %in% names(ENUM_SELECTION_CRITERION()))) {
  stop(paste0("Selection criterion provided is not valid [", selection_criterion, "]"))
}


# Hard code these for now. Very unlikely to change. 
id_variable = "SampleID_continuous"
sample_data(pso)$SampleID_continuous = as.character(sample_data(pso)$SampleID) # Is this necessary?
covariates = NULL
confounders = c("StudyCenter___vancouver_ndnominal",
                "StudyCenter___edmonton_ndnominal",
                "StudyCenter___winnipeg_ndnominal",
                "StudyCenter___toronto_ndnominal")
response_type = "dichotomous"

# Prerprocess the metadata and phyloseq object to match samples
filtered_metadata = filter_empty_samples(metadata, outcome)
matched_data = match_phyloseq_and_metadata_samples(pso, filtered_metadata, id_variable = id_variable)
pso = matched_data$phyloseq_obj
metadata = matched_data$metadata

# Split the data into train and test 
train_test_split_data = create_train_test_objects_by_class(pso = pso, metadata = metadata, 
                                                          sampleID_col_name = id_variable, binary_outcome = outcome)

PSO.TRAIN = train_test_split_data$pso.train
METADATA.TRAIN = train_test_split_data$metadata.train
PSO.TEST = train_test_split_data$pso.test
METADATA.TEST = train_test_split_data$metadata.test

# Create a smaller pso with less features for quicker testing. 
if (!is.null(opt.use_tiny)) {
  condensed_num_features = 3
  random_feature_IDs = sample(rownames(otu_table(PSO.TRAIN)), condensed_num_features)
  tiny_pso = PSO.TRAIN
  otu_table(tiny_pso) = subset(otu_table(tiny_pso), (rownames(otu_table(tiny_pso)) %in% random_feature_IDs))
  otu_table(tiny_pso) <- otu_table(tiny_pso)[1:condensed_num_features] 
  PSO.TRAIN = tiny_pso
}

# Split traiining data into k folds.
k = 10
if (class(PSO.TRAIN) == "phyloseq") {
  FOLDS = split_to_k_folds(t(otu_table(PSO.TRAIN)), k = k, phyloseq_obj = PSO.TRAIN)
} else {
  FOLDS = split_to_k_folds(t(PSO.TRAIN), k = k)
}

#### 1. Train on training data with RFECV to get the optimal number of features.
MiRKC.RFE_results.rfe_cv = MiRKC.RFE.CV_find_optimal(FOLDS = FOLDS,
                                          full_phyloseq_obj = PSO.TRAIN,
                                          y_variable = outcome,
                                          metadata = METADATA.TRAIN,
                                          kernel_function = kernel_function,
                                          covariates = covariates,
                                          confounders = confounders,
                                          selection_criterion = ENUM_SELECTION_CRITERION()$pval,
                                          id_variable = id_variable)

d = dim(otu_table(PSO.TRAIN))[1]
saveRDS(MiRKC.RFE_results.rfe_cv, 
        file = concat_paths(DIR_PATH.TMP_DATA_OUTPUT, MiRKC.RFE.output_filename("rfecv", k, timepoint, outcome, d)))

OPTIMAL_NUMBER_OF_FEATURES = MiRKC.RFE_results.rfe_cv$OPTIMAL_NUMBER_OF_FEATURES

#### 2. Re-train on training data with RFE to get ranked list of features.

MiRKC.RFE_results.train_rfe = MiRKC.RFE.rank_features(.X = PSO.TRAIN,
                                                  y_variable = outcome,
                                                  metadata = METADATA.TRAIN,
                                                  kernel_function = kernel_function,
                                                  covariates = covariates,
                                                  confounders = confounders,
                                                  response_type = response_type,
                                                  id_variable = id_variable,
                                                  selection_criterion = selection_criterion)

OPTIMAL_FEATURES_SELECTED = MiRKC.RFE_results.train_rfe$RANKED_FEATURES[1:OPTIMAL_NUMBER_OF_FEATURES]

#### 3. Test with MiRKC on training data by taking the optimal number of features (step 1) for the 
####    ranked list from step 2.

PSO.TRAIN.OPTIMAL = PSO.TRAIN
ot = otu_table(PSO.TRAIN.OPTIMAL)
otu_table(PSO.TRAIN.OPTIMAL) = ot[rownames(ot) %in% OPTIMAL_FEATURES_SELECTED, ] ## TODO - create helper func.

MiRKC.RFE_results.train_model = MiRKC.RFE.fit_MiRKC_model(PSO.TRAIN.OPTIMAL,
                                                   kernel_function,
                                                   metadata = METADATA.TRAIN,
                                                   y_variable = outcome,
                                                   covariates = covariates,
                                                   confounders = confounders,
                                                   response_type = response_type,
                                                   id_variable = id_variable)
MiRKC.RFE_results.train_model = MiRKC.RFE_results.train_model[[outcome]]

#### 4. Test with MiRKC on testing data by taking the optimal number of features (step 1) for the 
####    ranked list from step 2.

# Subset the test phyloseq object to keep only the optimal features.
PSO.TEST.OPTIMAL = PSO.TEST
ot = otu_table(PSO.TEST.OPTIMAL)
otu_table(PSO.TEST.OPTIMAL) = ot[rownames(ot) %in% OPTIMAL_FEATURES_SELECTED, ] ## TODO - create helper func.

MiRKC.RFE_results.test_model = MiRKC.RFE.fit_MiRKC_model(PSO.TEST.OPTIMAL,
                                        kernel_function,
                                        metadata = METADATA.TEST,
                                        y_variable = outcome,
                                        covariates = covariates,
                                        confounders = confounders,
                                        response_type = response_type,
                                        id_variable = id_variable)
MiRKC.RFE_results.test_model = MiRKC.RFE_results.test_model[[outcome]]

message("#### OPTIMAL NUMBER OF FEATURES IDENTIFIED FROM RFE CV:")
print(OPTIMAL_NUMBER_OF_FEATURES)

message("#### TRAINED TO GENERATE RANKED LIST OF FEATURES (MOST TO LEAST IMPT): ")
print(MiRKC.RFE_results.train_rfe$RANKED_FEATURES)

message("#### OPTIMAL FEATURES ARE:")
print(OPTIMAL_FEATURES_SELECTED)

message("#### TEST USING OPTIMAL FEATURES ON TRAINING SAMPLES.")
print(MiRKC.RFE_results.train_model)

message("#### TEST USING OPTIMAL FEATURES ON TESTING SAMPLES.")
print(MiRKC.RFE_results.test_model)

#### Summarize all the results.
MiRKC.RFE.summarize_results(MiRKC.RFE_results.rfe_cv,
                            MiRKC.RFE_results.train_rfe,
                            MiRKC.RFE_results.train_model,
                            MiRKC.RFE_results.test_model,
                            dir_path.results_output.CURRENT_TASK_ID,
                            CURRENT_TASK_ID,
                            train_test_split_data,
                            dim(otu_table(PSO.TRAIN))[2],
                            dim(otu_table(PSO.TEST))[2])

print("#############################################################################")
print("#############################################################################")
print(paste0("RFE operation complete: ", timestamp(prefix="", suffix="")))
print("#############################################################################")
print("#############################################################################")
