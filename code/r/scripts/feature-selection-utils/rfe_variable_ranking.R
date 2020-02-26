source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_utils()
import.project_helper_script("split_phyloseq_train_test.R")
import.data_modeling_utils("mirkat_utilities.R")

# ENUM Functions for standardizing possible options.
ENUM_KERNEL_FUNCTION <- function() {
  list(
    bray = "bray",
    jaccard = "jaccard",
    unifrac = "UniFrac",
    wunifrac = "wUniFrac",
    rbf = "RBF",
    gaussian = "Gaussian"
  )
}
ENUM_SELECTION_CRITERION <- function() {
  list(
    pval = "pval",
    aic = "aic",
    bic = "bic",
    auc = "auc",
    precision = "precision",
    higher_is_impt = c("pval", "aic"),
    lower_is_impt = c("auc")
  )
}

MiRKC.RFE.rank_features <- function(.X, 
                                    y_variable,
                                    metadata,
                                    kernel_function,
                                    covariates = NULL,
                                    confounders = NULL,
                                    response_type = "dichotomous",
                                    num_optimal_features = NULL, 
                                    features_are_cols = TRUE,
                                    FUNC_DISPLAY_NAME = "MiRKC.RFE.rank_features()",
                                    id_variable = "SubjectNumber_continuous",
                                    selection_criterion = ENUM_SELECTION_CRITERION()$pval) {
  ## Ranks variables for RFE.
  # @param X                      n x d dataframe containing n samples and d features OR a phyloseq object
  # @param y_variable             a string for the name of a numeric response variable
  # @param kernel_function        a string for the kernel function to use, see @ENUM_KERNEL_FUNCTION for options
  # @param response_type          default = "dichotomous", other option is "continuous"
  # @param num_optimal_features   integer for the optimal number of features. Used to subset the ranked listed.
  #                               This optimal number would have more meaning for a forward selection process, as it would act as a stopping 
  #                               criterion.
  # @param features_are_cols      logical for whether or not the features are columns

  USE_PHYLOSEQ = FALSE
  
  # Use either the phyloseq object, or X. 
  if (is.null(.X) || is.null(y_variable) || is_vector_blank(y_variable)) {
    stop("Input data are null.")
  } 
  # Check class of the primary data input. 
  if (class(.X) == "phyloseq") {
    USE_PHYLOSEQ = TRUE
    X = .X
  } else if (class(.X) == "matrix") {
    if (features_are_cols) {
      message("Transposing table so that column features become row features (consistent with phyloseq downstream handling).")
      X = t(.X)
    }
  }
  
  num_samples_total = if(USE_PHYLOSEQ) dim(otu_table(X))[2] else dim(X)[2]
  num_features_total = if (USE_PHYLOSEQ) dim(otu_table(X))[1] else dim(X)[1]
  num_features_eliminated = 0
  RANKED_FEATURES = c()
  ITERATION_SCORES = data.frame(matrix(nrow=0, ncol=3))
  colnames(ITERATION_SCORES) = c("feature_ID", "score", "iteration")
  
  for (i in 1:num_features_total) {
    print_process_message(FUNC_DISPLAY_NAME, paste0("Recursive feature elimination step = [", i, "/", num_features_total, "]"))
    remaining_feature_IDs = if (USE_PHYLOSEQ) rownames(otu_table(X)) else rownames(X) 
    num_features_left = length(remaining_feature_IDs)
    
    # Set up the temporary data frame for evaluating scores subsets.
    SUBSET_SCORES = data.frame(matrix(nrow=num_features_left, ncol=3))
    colnames(SUBSET_SCORES) = c("feature_ID", "score", "iteration")

    # Create individual models for each num_features_left - 1 subset with a feature "temporarily" removed.  
    for (j in 1:num_features_left) {
      
      # Remove the current feature from the data input matrix temporarily
      curr_feature_ID = remaining_feature_IDs[j]
      if (num_features_left > 1) {
        curr_X = MiRKC.RFE.eliminate_feature(X, curr_feature_ID, USE_PHYLOSEQ)
      }
      
      # Create the model. Should probably create a class for this output and standardize it. 
      test_result = MiRKC.RFE.fit_MiRKC_model(curr_X,
                                               kernel_function,
                                               metadata = metadata, 
                                               y_variable = y_variable,
                                               covariates = covariates,
                                               confounders = confounders,
                                               response_type = response_type,
                                               id_variable = id_variable)
      test_result = test_result[[y_variable]]
      
      # Save the test result based on the selection criterion.
      SUBSET_SCORES = MiRKC.RFE.append_test_result(j, SUBSET_SCORES, test_result, curr_feature_ID, selection_criterion)
      print_process_message(FUNC_DISPLAY_NAME, 
                          paste0("j-th feature temporarily removed out of remaining variables [", 
                          j, "/", num_features_left, "]"))   
      print_inner_loop_progress(num_features_eliminated, num_features_left, FUNC_DISPLAY_NAME)  
    }

    # Store in larger data frame for debugging / proof of process. 
    SUBSET_SCORES[, "iteration"] = i
    ITERATION_SCORES = rbind(SUBSET_SCORES, ITERATION_SCORES)
    
    # Find the least important feature based on the selection criterionn
    lowest_ranked_feature = MiRKC.RFE.identify_least_important_feature(SUBSET_SCORES, selection_criterion)
    RANKED_FEATURES = c(lowest_ranked_feature$feature_ID, RANKED_FEATURES)
    
    # Eliminate the variable from the last iteration that is of the lowest rank
    if (num_features_left > 1) {
      X = MiRKC.RFE.eliminate_feature(X, lowest_ranked_feature$feature_ID, USE_PHYLOSEQ)
    }
    num_features_eliminated = num_features_eliminated + 1
  }
  print_process_message(FUNC_DISPLAY_NAME, paste0("Last and top ranked feature: ", remaining_feature_IDs))
  return(list(RANKED_FEATURES = RANKED_FEATURES, 
              ITERATION_SCORES = ITERATION_SCORES, 
              SELECTION_CRITERION = selection_criterion,
              NUM_SAMPLES = num_samples_total,
              NUM_FEATURES = num_features_total,
              MODEL_PARAMETERS = list(
                KERNEL_FUNCTION = kernel_function,
                OUTCOME = y_variable,
                COVARIATES = covariates,
                CONFOUNDERS = confounders)))
} 

MiRKC.RFE.compute_num_tests_total <- function(d) { 
  if (d == 0) { return(d) }
  return(d + MiRKC.RFE.compute_num_tests_total(d-1)) 
}

MiRKC.RFE.eliminate_feature <- function(X, feature_ID, USE_PHYLOSEQ) {
  ## Eliminates a feature from the input data, which is either a matrix or a phyloseq object.
  if (USE_PHYLOSEQ) {
    otu_table(X) = subset(otu_table(X), !(rownames(otu_table(X)) %in% c(feature_ID)))
  } else {
    X = subset(X, !(rownames(X) %in% c(feature_ID))) 
  }
  return(X)
}

MiRKC.RFE.get_score <- function(test_result, selection_criterion) {
  # Retrieve the score based on which selection criterion 
  if (selection_criterion == ENUM_SELECTION_CRITERION()$pval) {
    score = test_result$indivP
    
  } else if (selection_criterion == ENUM_SELECTION_CRITERION()$aic) {
    score = test_result$model$aic
    
  } else if (selection_criterion == ENUM_SELECTION_CRITERION()$auc) {
    # take test data and predict on it...... but for now perhaps just use the training data which was used for
    # estimating the weights for the features
    predictions = predict(object = test_result$model, type = "response")
    score = Metrics::auc(actual = test_result$model$y, predicted = predictions)   
  } else {
    stop("Selection criterion does not match any of the available types. Check @ENUM_SELECTION_CRITERION().")
  }
  return(score)
}

MiRKC.RFE.append_test_result <- function(j, SUBSET_SCORES, test_result, feature_ID, selection_criterion) {
  # Appends the test result (score) to the table.
  SUBSET_SCORES[j, "score"] = MiRKC.RFE.get_score(test_result, selection_criterion)
  SUBSET_SCORES[j, "feature_ID"] = feature_ID
  return(SUBSET_SCORES)
}

MiRKC.RFE.fit_MiRKC_model <- function(X,
                                 kernel_function,
                                 metadata = metadata, 
                                 y_variable,
                                 covariates = NULL,
                                 confounders = NULL,
                                 response_type,
                                 id_variable = id_variable) {
  # Handles fitting the model. Is a wrapper around MiRKC and the kernel construction.
  
  # Construct kernel(s). This function does not yet handle generic matrices, and works on phyloseq objects.
  kernel_construction_output = construct_kernels(X, # A phyloseq object at the moment.
                                                 kernel_functions = c(kernel_function),
                                                 id_variable = id_variable,
                                                 save_kernels_to_file = FALSE,
                                                 FUNC_DISPLAY_NAME = "MiRKC.RFE.fit_MiRKC_model()")
  kernel_matrices = kernel_construction_output$kernels
  kernel_name_order = kernel_construction_output$kernel_order

  # Run regression test for this d - j case.
  MiRKC_output = list()
  MiRKC_output <- subset_and_run_univariate_mirkat(
    metadata_source = metadata, 
    outcomes = c(y_variable),
    kernels = kernel_matrices,
    output_obj = MiRKC_output, 
    covariates = covariates,
    confounders = confounders,
    id_variable = id_variable,
    outcome_type = ifelse(response_type == "dichotomous", "D", "C"),
    test_type = "prediction",
    verbose = FALSE
  )
  return(MiRKC_output)
}

print_inner_loop_progress <- function(num_features_eliminated, num_features_left, FUNC_DISPLAY_NAME) {
  print_process_message(FUNC_DISPLAY_NAME, 
                        paste0("Number of features eliminated so far: ", num_features_eliminated))
  print_process_message(FUNC_DISPLAY_NAME, 
                        paste0("Number of features in current kernel [", num_features_left - 1, "]"))
  print("\n")
}

is_higher_score_impt <- function(selection_criterion) {
  return(selection_criterion %in% ENUM_SELECTION_CRITERION()$higher_is_impt)
}

MiRKC.RFE.identify_least_important_feature <- function(SUBSET_SCORES, selection_criterion) {
  print(SUBSET_SCORES)
  ## Eliminates the lowest ranked feature based on a given criterion.
  if (is_higher_score_impt(selection_criterion)) {
    SUBSET_SCORES = SUBSET_SCORES[order(SUBSET_SCORES$score), ]
    # Smaller score (e.g. pval) for the subset (in absence of a given OTU) means less important feature. 
    # Larger score (e.g. pval) for the subset (in absence of a given OTU) means more important feature.
  } else if (!is_higher_score_impt(selection_criterion)) {
    SUBSET_SCORES = SUBSET_SCORES[-order(SUBSET_SCORES$score), ]
    # Larger score (e.g. accuracy) for the subset (in absence of a given OTU) means less important feature. 
    # Smaller score (e.g. accuracy) for the subset (in absence of a given OTU) means more important feature.
  }
  feature = SUBSET_SCORES[1, ]
  return(feature)
}

are_samples_present <- function(phyloseq_obj, metadata, id_variable = "SubjectNumber_continuous") {
  if (is.null(metadata[[id_variable]])) {
    stop("Provided id_variable for checking is not present in the metadata.")
  }
  phyloseq_IDs = sample_data(phyloseq_obj)[[id_variable]]
  if (is.null(phyloseq_IDs)) {
    stop("Provided id_variable for checking is not present in the metadata.")
  }
  if (!all(phyloseq_IDs %in% metadata[[id_variable]])) {
    stop(paste0("Not all IDs [", id_variable, "] in the phyloseq object are present in the provided metadata."))
  }
  return(TRUE)
}

MiRKC.RFE.CV_find_optimal <- function(FOLDS, 
                                     full_phyloseq_obj,
                                     y_variable, 
                                     metadata, 
                                     kernel_function,
                                     covariates = NULL,
                                     confounders = NULL,
                                     response_type = "dichotomous",
                                     features_are_cols = TRUE,
                                     id_variable = "SubjectNumber_continuous",
                                     selection_criterion = ENUM_SELECTION_CRITERION()$pval,
                                     k = 3) {
  ## Finds the optimal number of features
  # metadata - a dataframe containing the metadata. It could be larger than the .X (phyloseq object). All samples in the primary
  #            data source must be present in the metadata.
  
  ## TODO - create an object that takes in typical parameters for modeling.
  # Pass that around instead of the params directly. 
  # Params: 1. X, 2. y variable, 3. covariates, 4. confounders. 
  
  # Check that all samples in the input data are present in the metadata.
  if (!are_samples_present(full_phyloseq_obj, metadata, id_variable)) {
    stop("Subset the phyloseq object to contain IDs present in the provided metadata.")
  } else {
    message("Confirmed that phyloseq object contains IDs present in the provided metadata.")
  }
  
  d = dim(otu_table(full_phyloseq_obj))[1]
  
  ## Finds the optimal number of features for RFE using CV.
  # for each fold, use the training samples and run RFE on it
  # then do an iterative test using increasing number of variables from the ranked list to test on the test samples
  # save it in a table, the performance
  K_FOLD_RESULTS = list()
  K_FOLD_RESULTS[["FOLDS"]] = list()
  for (i in 1:k) {
    print("*************************************************************************************")
    print(paste("Handling fold = ", i))
    fold_data = FOLDS[[i]]
    
    train_phyloseq_obj = fold_data$train_phyloseq_obj
    test_phyloseq_obj = fold_data$test_phyloseq_obj
    
    print(paste0("Train samples (n = ", dim(otu_table(train_phyloseq_obj))[2], ", d = ", dim(otu_table(train_phyloseq_obj))[1], "): "))
    print(otu_table(train_phyloseq_obj)[, 1:2])
    print(paste0("Test samples (n = ", dim(otu_table(test_phyloseq_obj))[2], ", d = ", dim(otu_table(test_phyloseq_obj))[1], "): "))
    print(otu_table(test_phyloseq_obj)[, 1:2])
    
    # Train the model to identify ranking of variables.
    FEATURE_RANKING_RESULTS = MiRKC.RFE.rank_features(
      .X = train_phyloseq_obj,
      y_variable = y_variable,
      metadata = metadata,
      kernel_function = kernel_function,
      covariates = covariates,
      confounders = confounders,
      response_type = response_type,
      features_are_cols = features_are_cols,
      id_variable = id_variable,
      selection_criterion = selection_criterion)
    
    K_FOLD_RESULTS[["FOLDS"]][[i]] = list()
    K_FOLD_RESULTS[["FOLDS"]][[i]][["RANKED_FEATURES"]] = FEATURE_RANKING_RESULTS$RANKED_FEATURES
    K_FOLD_RESULTS[["FOLDS"]][[i]][["ITERATION_SCORES"]] = FEATURE_RANKING_RESULTS$ITERATION_SCORES
    K_FOLD_RESULTS[["FOLDS"]][[i]][["SELECTION_CRITERION"]] = FEATURE_RANKING_RESULTS$SELECTION_CRITERION
    
    message(paste0("MiRKC.RFE.find_optimal_number_features(): Ranking complete [k = ", i, "], total of ", d, " features ranked (most to least impt)."))
    print(FEATURE_RANKING_RESULTS$RANKED_FEATURES)
    
    scores_for_fold = c()
    # Move this outside the loop.
    for (p in 1:d) {
      test_result = MiRKC.RFE.test_rfe_model(ranked_features = FEATURE_RANKING_RESULTS$RANKED_FEATURES,
                               p = p,
                               test_samples = test_phyloseq_obj,
                               kernel_function = kernel_function,
                               metadata = metadata,
                               y_variable = y_variable,
                               covariates = covariates,
                               confounders = confounders,
                               response_type = response_type,
                               features_are_cols = features_are_cols,
                               id_variable = id_variable)
      
      score = MiRKC.RFE.get_score(test_result, selection_criterion)
      scores_for_fold = c(scores_for_fold, score)
    }
    
    K_FOLD_RESULTS[["FOLDS"]][[i]][["SCORES_FOR_FOLD"]] = scores_for_fold
  }
  K_FOLD_RESULTS[["DESCRIPTION"]] = "The scores for the fold from left to right represent scores for increasing number of top variables tested. E.g. the score in the 1-th index is for 1 top variable and the score in the 2-th index is for the 2 top variables tested."
  K_FOLD_RESULTS[["ORDER_OF_NUM_VARIABLES_TESTED"]] = seq(1, d, by = 1)
  K_FOLD_RESULTS[["k"]] = k
  K_FOLD_RESULTS[["d"]] = d
  
  # Summarize the CV performance and identify the optimal number of features.
  CV_averages = compute_CV_averages_per_num_var_tested(K_FOLD_RESULTS)
  optimal_num_features = find_optimal_num_features(CV_averages, selection_criterion)

  library(ggplot2)
  p = ggplot(data.frame(x = K_FOLD_RESULTS[["ORDER_OF_NUM_VARIABLES_TESTED"]], y = CV_averages), aes(x = x, y = y)) 
  p = p + geom_point() + geom_line(color="blue") + geom_vline(xintercept = optimal_num_features)
  p = p + xlab("Number of Features") + ylab("Average CV score (k-fold)")
  p = p + labs(title = paste0("Identifying the Optimal Number of Features with ", k, "-fold CV (", selection_criterion, ")"))
  
  K_FOLD_RESULTS[["PLOT"]] = p
  K_FOLD_RESULTS[["CV_AVERAGES"]] = CV_averages
  K_FOLD_RESULTS[["OPTIMAL_NUMBER_OF_FEATURES"]] = optimal_num_features

  return(K_FOLD_RESULTS)
}

MiRKC.RFE.test_rfe_model <- function(ranked_features, p, test_samples, fitted_model = NULL,
                           kernel_function, metadata, y_variable, covariates, confounders,
                           id_variable = "SubjectNumber_continuous", 
                           response_type = "dichotomous", 
                           features_are_cols = TRUE) {
  ## Runs the test
  # Use fitted_model when the model with the trained coefficients are given.
  
  # Select the p top features from the ranked list of features
  p_top_features = ranked_features[1:p]
  message(paste0("MiRKC.RFE.test_rfe_model(): Testing p = ", p, " top features"))
  print(p_top_features)
  test_samples.top_p = test_samples
  ot = otu_table(test_samples.top_p)
  otu_table(test_samples.top_p) = ot[rownames(ot) %in% p_top_features, ]

  test_result = MiRKC.RFE.fit_MiRKC_model(test_samples.top_p,
                                             kernel_function,
                                             metadata = metadata, 
                                             y_variable = y_variable,
                                             covariates = covariates,
                                             confounders = confounders,
                                             response_type = response_type,
                                             id_variable = id_variable)
  test_result = test_result[[y_variable]]
  
  return(test_result)
} 

split_to_k_folds <- function(data_samples, k = 3, phyloseq_obj = NULL) {
  # Randomly shuffle the samples. The indices for train and test correspond to order of this. 
  data_samples <- data_samples[sample(nrow(data_samples)), ]
  # Create k equal size folds out of the samples
  folds <- cut(seq(1, nrow(data_samples)), breaks = k, labels = FALSE)
  
  cv_data = list()
  for (i in 1:k){
    # Segment your data by fold using the which() function 
    test_indices <- which(folds == i, arr.ind = TRUE)
    training_indices <- which(folds != i, arr.ind = TRUE)
    test_samples <- data_samples[test_indices, ]
    training_samples <- data_samples[training_indices, ]
    
    cv_data[[i]] = list()
    cv_data[[i]][["test_samples"]] = test_samples
    cv_data[[i]][["test_indices"]] = test_indices
    cv_data[[i]][["training_samples"]] = training_samples
    cv_data[[i]][["training_indices"]] = training_indices
    cv_data[[i]][["fold"]] = i
    cv_data[[i]][["performance"]] = NULL
    cv_data[[i]][["fitted_model"]] = NULL
    
    # Test these.
    if (!is.null(phyloseq_obj)) {
      test_phyloseq_obj = phyloseq_obj
      otu_table(test_phyloseq_obj) = otu_table(test_phyloseq_obj)[, (colnames(otu_table(test_phyloseq_obj)) %in% rownames(test_samples))]
      cv_data[[i]][["test_phyloseq_obj"]] = test_phyloseq_obj
      
      train_phyloseq_obj = phyloseq_obj
      otu_table(train_phyloseq_obj) = otu_table(train_phyloseq_obj)[, (colnames(otu_table(train_phyloseq_obj)) %in% rownames(training_samples))]
      cv_data[[i]][["train_phyloseq_obj"]] = train_phyloseq_obj
    }
  }
  cv_data[["shuffled_data_fixed_indices"]] = data_samples
  return(cv_data)
}

find_optimal_num_features <- function(CV_averages, selection_criterion) {
  if (selection_criterion != ENUM_SELECTION_CRITERION()$pval) {
    best = max(CV_averages)
  } else {
    best = min(CV_averages)
  }
  indices = which(CV_averages %in% best)
  
  # Return the first one that matches. Very unlikely there'll be more than one.
  return(indices[1])
}

compute_CV_averages_per_num_var_tested <- function(k_fold_results) {
  # Determine number of folds
  k = k_fold_results$k
  # Determine number of features
  d = k_fold_results$d
  
  averages = rep(0, d)
  for (i in 1:k) {
    fold_data = k_fold_results[["FOLDS"]][[i]]
    scores_for_fold = fold_data[["SCORES_FOR_FOLD"]]
    averages = averages + scores_for_fold
  }
  averages = averages / d 
  return(averages)
}

filter_empty_samples <- function(metadata, outcome) {
  message(paste("Initial dimensions of metadata:", dim(metadata))[1], " ",  dim(metadata)[2])
  .metadata <- metadata[!is_vector_el_blank(metadata[[outcome]]), ]
  message(paste("Dimensions of metadata after filtering out samples with a blank outcome [", outcome, "]:", dim(.metadata))[1], " ", dim(.metadata)[2])
  return(.metadata)
}

match_phyloseq_and_metadata_samples <- function(phyloseq_obj, metadata, outcome, id_variable = "SubjectNumber_continuous") {
  ## Matches the samples between the phyloseq object and the metadata. 
  
  # The primary objective is that before running rfe or rfecv, the phyloseq object and metadata should be matched
  # beforehand. When CV splits the data into folds, the metadata will have to be split and passed around the same way.
  
  .metadata = metadata
  
  .phyloseq_obj = phyloseq_obj
  .otu_table = otu_table(.phyloseq_obj)
  .otu_table = .otu_table[, (colnames(.otu_table) %in% as.character(.metadata[[id_variable]]))]
  otu_table(.phyloseq_obj) = .otu_table
  
  phyloseq_IDs = colnames(otu_table(.phyloseq_obj))
  
  # Check that the dimensions (num samples) match between the metadata and phyloseq object. 
  if (dim(.metadata)[1] == length(phyloseq_IDs)) {
    message("Metadata and phyloseq object expected number of samples match.")
  }
  # Check that the IDs match between the two.
  metadata_IDs = as.character(.metadata[[id_variable]])
  if ((metadata_IDs %in% phyloseq_IDs) && (phyloseq_IDs %in% metadata_IDs)) {
    message("Metadata and phyloseq object IDs match.")
  }
  return(list(phyloseq_obj = .phyloseq_obj, metadata = .metadata))
}

MiRKC.RFE.output_filename <- function(file_base = "rfecv", k, timepoint, outcome, d) {
  return(paste0(file_base, "_", k, "fold_tp", timepoint, "_", outcome, "_", d, "features", ".RData"))
}

MiRKC.RFE.summarize_results <- function(MiRKC.RFE_results.rfe_cv, 
                                        MiRKC.RFE_results.train_rfe,
                                        MiRKC.RFE_results.train_model,
                                        MiRKC.RFE_results.test_model,
                                        dir_path.output,
                                        task_ID,
                                        MiRKC.RFE_results.train_test_split_data,
                                        num_training_samples,
                                        num_testing_samples) {
  # Save the objects
  file_path.output.objects = concat_paths(dir_path.output, paste0(task_ID, "_objects.RData"))
  save(MiRKC.RFE_results.train_test_split_data,
       MiRKC.RFE_results.rfe_cv, 
       MiRKC.RFE_results.train_rfe,
       MiRKC.RFE_results.train_model,
       MiRKC.RFE_results.test_model,
       list = c("MiRKC.RFE_results.train_test_split_data",
                "MiRKC.RFE_results.rfe_cv",
                "MiRKC.RFE_results.train_rfe",
                "MiRKC.RFE_results.train_model",
                "MiRKC.RFE_results.test_model"), 
       file = file_path.output.objects)
  
  # Save the plot
  file_path.output.plot = concat_paths(dir_path.output, paste0(task_ID, ".png"))
  ggsave(file_path.output.plot, plot = MiRKC.RFE_results.rfe_cv$PLOT)

  # Add a text summary
  OPTIMAL_NUMBER_OF_FEATURES = MiRKC.RFE_results.rfe_cv$OPTIMAL_NUMBER_OF_FEATURES
  MODEL_PARAMS = MiRKC.RFE_results.train_rfe$MODEL_PARAMETERS
  summary = c(
    paste0("Number of training vs. testing samples:\t", paste(num_training_samples, num_testing_samples, collapse = "/")),
    "\n", 
    paste0("Number of features:\t", MiRKC.RFE_results.rfe_cv$d),
    "\n", 
    paste0("Cross validation folds (k-fold):\t", MiRKC.RFE_results.rfe_cv$k),
    "\n", 
    paste0("Model selection criterion:\t", MiRKC.RFE_results.rfe_cv$SELECTION_CRITERION),
    "\n", 
    paste0("Kernel function:\t", MODEL_PARAMS$KERNEL_FUNCTION),
    "\n", 
    paste0("Outcome variable:\t", MODEL_PARAMS$OUTCOME),
    "\n", 
    paste0("Covariates:\t", paste(MODEL_PARAMS$COVARIATES, collapse = ",")),
    "\n", 
    paste0("Confounders:\t", paste(MODEL_PARAMS$CONFOUNDERS, collapse = ",")),
    "\n",
    paste0("Optimal number of features determined from RFE CV on training samples:\t", OPTIMAL_NUMBER_OF_FEATURES),
    "\n",
    paste0("Ranked list of features (most to least impt) from training with training samples:\t", paste(MiRKC.RFE_results.train_rfe$RANKED_FEATURES, collapse = ",")),
    "\n",
    paste0("Optimal features:\t", paste(MiRKC.RFE_results.train_rfe$RANKED_FEATURES[1:OPTIMAL_NUMBER_OF_FEATURES], collapse = ",")),
    "\n",
    paste0("Training performance with optimal features:\t", MiRKC.RFE.get_score(MiRKC.RFE_results.train_model, MiRKC.RFE_results.train_rfe$SELECTION_CRITERION)),
    "\n",
    paste0("Testing performance with optimal features:\t", MiRKC.RFE.get_score(MiRKC.RFE_results.test_model, MiRKC.RFE_results.train_rfe$SELECTION_CRITERION))
  )
  file_path.output.summary_text = concat_paths(dir_path.output, paste0(task_ID, "_summary_info.txt"))
  file_connection = file(file_path.output.summary_text)
  writeLines(summary, file_connection)
  close(file_connection)
  
  print("\n")
  print(summary)
  print("TASK OUTPUT DIRECTORY:")
  print(dir_path.output)
}
