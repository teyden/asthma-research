source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_utils()

install_package_if_missing(package_names = c("MiRKAT", "MiRKC", "rdetools", "KRLS", "kernlab", "rlist", "phyloseq"))
library(rlist)
library(kernlab)
library(KRLS)
library(rdetools)
library(MiRKAT)
library(phyloseq)
library(MiRKC)

subset_and_run_univariate_mirkat <- function(
  metadata_source, 
  outcomes, 
  kernels,
  output_obj,
  confounders = NULL, 
  covariates = NULL, 
  verbose = FALSE, 
  id_variable = "SubjectNumber_continuous",
  outcome_type = "C",
  test_type = "association") {
  # Microbiome kernel-based regression testing (association and classification) wrapper. Optimized implementation.
  
  # Based on data availability of the confounders and response variables in question, the original kernels are subsetted.
  # Subsetting kernels allows optimization by preventing re-computation of kernels when samples are removed. 
  # Samples are removed when the variable being modeled is missing for that sample. 
  # E.g. If 10 subjects are removed because they don't have 
  # any data for the variable "Birth_mode", then those 10 subjects are removed from the kernel (associated row/column).
  
  confounders <- if (is.null(confounders)) c() else confounders
  covariates <- if (is.null(covariates)) c() else covariates
  
  if (sum(!outcomes %in% names(metadata_source))) stop("Variables in outcomes must exist in the metadata dataframe.")
  
  N_outc <- length(outcomes)
  i <- 1
  
  for (outcome in outcomes) {
    if (verbose == TRUE) {
      message(paste0("Handling outcome... ", i, "/", N_outc))
      print(outcome)
    }
    
    i <- i + 1 
    
    # Skip if the outcome is one of the confounders
    if (outcome %in% confounders) next
    
    if (!is.null(covariates) || length(covariates) > 0) {
      output_obj[[outcome]][["model_by_covariates"]] <- list()
      
      N_cov <- length(covariates)
      k <- 1
      for (covariate in covariates) {
        if (verbose == TRUE) {
          message(paste0("Handling covariate... ", k, "/", N_cov))
          print(outcome)
        }
        
        curr_output <- handle_running_mirkat(metadata_source, outcome, confounders, 
                                             kernels, id_variable, outcome_type = outcome_type, covariate = covariate, k = k, 
                                             verbose = verbose, test_type = test_type)
        if (is.null(curr_output)) {
          next
        }
        output_obj[[outcome]][["model_by_covariates"]][[covariate]] <- curr_output
        k <- k + 1
      }
    } else {
      curr_output <- handle_running_mirkat(metadata_source, outcome, confounders, 
                                           kernels, id_variable, outcome_type = outcome_type, covariate = NULL, k = 1, 
                                           verbose = verbose, test_type = test_type)
      if (is.null(curr_output)) {
        next
      }
      output_obj[[outcome]] <- curr_output
    }
  }
  return(output_obj)
}

handle_running_mirkat <- function(metadata_source, outcome, confounders, kernels, 
                                  id_variable, covariate = NULL, k = NULL, 
                                  outcome_type = "D",
                                  verbose = FALSE,
                                  test_type = "association") {
  
  # Kernels must have entries that exist in the metadata. The metadata could have more samples than the kernel.
  # Check at least one of the kernels.
  if (!all(colnames(kernels[[1]]) %in% metadata_source[[id_variable]])) {
    stop("handle_running_mirkat(): Kernels must contain samples that exist in the metadata source.")
  }
  
  # Filter out samples in the metadata that aren't present in the kernels.
  metadata_source.subset = metadata_source[metadata_source[[id_variable]] %in% colnames(kernels[[1]]), ]
  
  # Take out empty rows. MiRKAT does not handle NA's. After this operation, the metadata could be smaller than the kernel.
  metadata_source.subset <- na.exclude(metadata_source.subset[, c(outcome, confounders, covariate, id_variable)])

  # Skip if after excluding, the metadata subset is empty
  if (is.null(metadata_source.subset) || dim(metadata_source.subset)[1] == 0) {
    return(NULL)
  }
  
  y = metadata_source.subset[, c(outcome, id_variable)]
  # Skip if the outcome vector is constant, null or empty
  if (is_constant(y[[outcome]]) || is.null(y[[outcome]]) || length(y[[outcome]])[1] == 0) {
    return(NULL)
  }
  
  X = NULL
  if (!is.null(covariate)) {
    if (length(confounders)[1] == 0) {
      X = metadata_source.subset[, c(covariate, id_variable)]
    } else {
      X = metadata_source.subset[, c(confounders, covariate, id_variable)]
    }  
  }
  
  # Subset the kernel's row and columns (it's nxn) to remove those subjects that have empty values for the
  # variables of interest (outcome, covariates, confounders)  
  subsetted_kernels <- list()
  for (i in 1:length(kernels)) {
    kernel <- kernels[[i]]
    
    # Match the IDs between the kernel and the outcome, covariates, etc. Assumes the metadata IDs are all present in the kernels.
    metadata_IDs = metadata_source.subset[[id_variable]]
    subsetted_kernel <- kernel[metadata_IDs, metadata_IDs]
    subsetted_kernels[[i]] <- subsetted_kernel
    if (verbose == TRUE && i == 1 && k == 1) {
      print(dim(subsetted_kernel))
    }
    
    # Check that y, X and kernel are consistent (dimensions match and that ID orders match).
    if (!(dim(y)[1] == dim(subsetted_kernel)[1])) {
      stop("Number of samples in y and kernel do not match.")
    }
    if (!all(y[[id_variable]] == colnames(subsetted_kernel))) {
      stop("Order of sample IDs between y and kernel do not match.")
    }
    if (!is.null(X)) {
      if (!(dim(X)[1] == dim(y)[1])) {
        stop("Number of samples in X do not match y and kernel.")
      }
      if (!all(X[[id_variable]] == y[[id_variable]])) {
        stop("Order of sample IDs between X and y do not match.")
      }
      if (!all(y[[id_variable]] == colnames(subsetted_kernel))) {
        stop("Order of sample IDs between X and kernel do not match.")
      }
    }
  }
  if (test_type == "association") {
    return(MiRKAT(y = as.matrix(y[[outcome]]), X = X, K = subsetted_kernels, out_type = outcome_type))
  } else if (test_type == "prediction") {
    return(MiRKC(y = as.matrix(y[[outcome]]), X = X, K = subsetted_kernels, out_type = outcome_type))
  }
}

adjust_p_for_kernels <- function(kernel_names, data_table, p.adjust_method = "BH") {
  ## Performs p-value adjustment for different kernels used by MiRKAT
  for (kernel_name in kernel_names) {
    adj_name <- paste0("adj_", kernel_name)
    data_table[[adj_name]] <- p.adjust(data_table[[kernel_name]], method = p.adjust_method)
  }
  data_table <- round_numeric_columns(data_table, 6)
  col_name_order <- c()
  for (kernel_function in kernel_names) {
    col_name_order <- c(col_name_order, kernel_function, paste0("adj_", kernel_function))
  }
  data_table <- data_table[, col_name_order]
  return(data_table)
}

construct_kernels <- function(phyloseq_obj, 
                              kernel_functions = c("bray", "jaccard", "wUniFrac", "UniFrac", "RBF", "Gaussian"),
                              abs_output_dir = concat_paths(DIR_PATH.TMP_DATA_OUTPUT, "mirkat/mirkat-kernels"),
                              output_filename_base = "kernels_taxa_subset_d",
                              id_variable = "SubjectNumber",
                              save_kernels_to_file = TRUE,
                              FUNC_DISPLAY_NAME = CURR_FILE) {
  ## Constructs the kernels needed by MiRKAT, and outputs it to an RData object for use by scripts but also outputs it in a list.
  ## The phyloseq object should be rarefied for proper beta composition analysis.
  ## Perform any subsetting of the phyloseq object prior to providing to the function (if used for RFE).
  
  otu_table <- otu_table(phyloseq_obj)
  num_taxa <- dim(otu_table)[1]
  print_process_message(FUNC_DISPLAY_NAME, paste("Constructing kernel(s), for a total number of taxa = [", num_taxa, "] with functions applied = [", paste(kernel_functions, collapse = ", "), "]"))
  
  subject_ids <- as.character(sample_data(phyloseq_obj)[[id_variable]])
  kernels <- list()
  
  # Compute distances and kernels
  for (kernel_function in kernel_functions) {
    if (kernel_function == "bray" || kernel_function == "jaccard") {
      ## Add a constant temporarily to address rows with 0 counts during selection.
      otu_table(phyloseq_obj) <- otu_table(phyloseq_obj) + 1
    }
    
    if (kernel_function == "bray") {
      dist.BC <- as.matrix(phyloseq::distance(phyloseq_obj, method = "bray")) 
      kernel.BC <- D2K(dist.BC)
      dimnames(kernel.BC) <- list(subject_ids, subject_ids)
      kernels <- list.append(kernels, kernel.BC)
    }
    if (kernel_function == "jaccard") {
      dist.J <- as.matrix(phyloseq::distance(phyloseq_obj, method = "jaccard")) 
      kernel.J <- D2K(dist.J)
      dimnames(kernel.J) <- list(subject_ids, subject_ids)
      kernels <- list.append(kernels, kernel.J)
    }
    if (kernel_function == "wUniFrac") {
      dist.wUniFrac <- as.matrix(phyloseq::distance(phyloseq_obj, method = "wUniFrac", set.seed(100))) 
      kernel.wUniFrac <- D2K(dist.wUniFrac)
      dimnames(kernel.wUniFrac) <- list(subject_ids, subject_ids)
      kernels <- list.append(kernels, kernel.wUniFrac)
    }
    if (kernel_function == "UniFrac") {
      dist.UniFrac <- as.matrix(phyloseq::distance(phyloseq_obj, method = "UniFrac", set.seed(100))) 
      kernel.UniFrac <- D2K(dist.UniFrac)
      dimnames(kernel.UniFrac) <- list(subject_ids, subject_ids)
      kernels <- list.append(kernels, kernel.UniFrac)
    }
    if (kernel_function == "RBF") {
      kernel.RBF = rbfkernel(t(otu_table), sigma = 1, Y = NULL)
      dimnames(kernel.RBF) <- list(subject_ids, subject_ids)
      kernels <- list.append(kernels, kernel.RBF)
    }
    if (kernel_function == "Gaussian") {
      kernel.Gaussian = gausskernel(t(otu_table), sigma = 1)
      dimnames(kernel.Gaussian) <- list(subject_ids, subject_ids)
      kernels <- list.append(kernels, kernel.Gaussian)
    }
  }
  
  # Output to RData objects
  kernel_order <- kernel_functions
  print_process_message(FUNC_DISPLAY_NAME, paste("Kernel matrices created. Functions applied = [", paste(kernel_functions, collapse = ", "), "]"))
  
  file_path.debugger_objs <- concat_paths(abs_output_dir, paste0("last_mirkat_kernels__d", num_taxa, ".RData"))
  save(kernels, kernel_order, list = c("kernels", "kernel_order"), file = file_path.debugger_objs)
  
  if (save_kernels_to_file) {
    file_path.kernel_output <- concat_paths(abs_output_dir, paste0(output_filename_base, num_taxa, ".RData"))
    save(kernels, kernel_order, list = c("kernels", "kernel_order"), file = file_path.kernel_output)
    print_process_message(FUNC_DISPLAY_NAME, paste("Kernel objects saved to RData objects in path:", file_path.kernel_output))
  }

  return(list(kernels = kernels, kernel_order = kernel_order))
}