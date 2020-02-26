source("~/Projects/asthma/asthma-research/code/r/PROJECT_GLOBAL_CONFIGS.R")
import.project_utils()

create_train_test_objects <- function(pso, metadata, train_percentage=0.8, 
                                      dirpath.output=concat_paths(DIR_PATH.TMP_DATA_OUTPUT, "phyloseq-train-test-objects"),
                                      filename.base="this_is_a_test", sampleID_col_name="SampleID_continuous") {
  ## Randomly subsamples train/test sets from the pso object.
  ## Saves the training object and testing object inside of a folder, with corresponding names.
  
  if (!dir.exists(dirpath.output)) {
    dir.create(dirpath.output)
  }
  filepath.output = concat_paths(dirpath.output, paste0(filename.base, "__train_test_pso.RData"))

  # Samples are columns.
  otu_table = otu_table(pso)
  
  .metadata <- metadata[metadata[[sampleID_col_name]] %in% as.character(colnames(otu_table)), ]
  .metadata_class1 <- metadata 
  
  n_samples_train = round(train_percentage * ncol(otu_table))
  n_samples_test = dim(otu_table)[2] - n_samples_train

  sampleIDs_train = sample(colnames(otu_table), size = n_samples_train, replace = FALSE)
  sampleIDs_test = setdiff(colnames(otu_table), sampleIDs_train)
  
  metadata.train <- metadata[metadata[[sampleID_col_name]] %in% as.character(sampleIDs_train), ]
  metadata.test <- metadata[metadata[[sampleID_col_name]] %in% as.character(sampleIDs_test), ]
  
  pso.train = pso
  pso.test = pso
  otu_table(pso.train) <- otu_table(pso)[, sampleIDs_train]
  otu_table(pso.test) <- otu_table(pso)[, sampleIDs_test]

  phyloseq_objs.train_test <- list()
  phyloseq_objs.train_test[["pso.train"]] <- pso.train
  phyloseq_objs.train_test[["pso.test"]] <- pso.test
  phyloseq_objs.train_test[["metadata.train"]] <- metadata.train
  phyloseq_objs.train_test[["metadata.test"]] <- metadata.test
  
  save(phyloseq_objs.train_test, list = c(phyloseq_objs.train_test), file = filepath.output)
  print_process_message("create_train_test_objects()", paste0("Train/test objects output to ", filepath.output))
  
  return(phyloseq_objs.train_test)
}

create_train_test_objects_by_class <- function(pso, metadata, train_percentage=0.80, 
                                      dirpath.output=concat_paths(DIR_PATH.TMP_DATA_OUTPUT, "phyloseq-train-test-objects"),
                                      filename.base="this_is_a_test", sampleID_col_name="SampleID_continuous",
                                      binary_outcome=NULL) {
  ## Randomly subsamples train/test sets from the pso object, handling class imbalance.
  ## Saves the training object and testing object inside of a folder, with corresponding names.
  
  if (!dir.exists(dirpath.output)) {
    dir.create(dirpath.output)
  }
  filepath.output = concat_paths(dirpath.output, paste0(filename.base, "__train_test_pso.RData"))
  
  # Samples are columns.
  otu_table = otu_table(pso)
  
  # Subset the metadata to retain only samples that are in the microbiome data.
  .metadata <- metadata[metadata[[sampleID_col_name]] %in% as.numeric(colnames(otu_table)), ]
  
  # Need to remove the the blank cases. Causes issues. Sample IDs end up taking on an NA value somehow as well.
  # It's ok so they'll actually be left out during modeling anyway. So makes sure there is a more accurate percentage split between classes.
  .metadata <- .metadata[!is_vector_el_blank(.metadata[[binary_outcome]]), ]
  .metadata_class1 <- .metadata[.metadata[[binary_outcome]] == 1, ]
  .metadata_class2 <- .metadata[.metadata[[binary_outcome]] == 0, ]

  # Determine the number of samples for each class, so that subsampling is done in proportion to the class imbalances.
  n_samples_train_class1 = round(train_percentage * nrow(.metadata_class1))
  n_samples_train_class2 = round(train_percentage * nrow(.metadata_class2))
  
  # Sample the training samples for each class.
  sampleIDs_train_class1 = sample(.metadata_class1[[sampleID_col_name]], size = n_samples_train_class1, replace = FALSE)
  sampleIDs_train_class2 = sample(.metadata_class2[[sampleID_col_name]], size = n_samples_train_class2, replace = FALSE)
  
  # The test samples are the complement of the train samples.
  sampleIDs_test_class1 = setdiff(.metadata_class1[[sampleID_col_name]], sampleIDs_train_class1)
  sampleIDs_test_class2 = setdiff(.metadata_class2[[sampleID_col_name]], sampleIDs_train_class2)
  
  # Combine the two classes.
  sampleIDs_train <- as.character(c(sampleIDs_train_class1, sampleIDs_train_class2))
  sampleIDs_test <- as.character(c(sampleIDs_test_class1, sampleIDs_test_class2))
  
  metadata.train <- metadata[metadata[[sampleID_col_name]] %in% sampleIDs_train, ]
  metadata.test <- metadata[metadata[[sampleID_col_name]] %in% sampleIDs_test, ]
  
  pso.train = pso
  pso.test = pso
  otu_table(pso.train) <- otu_table(pso)[, sampleIDs_train]
  otu_table(pso.test) <- otu_table(pso)[, sampleIDs_test]
  
  phyloseq_objs.train_test <- list()
  phyloseq_objs.train_test[["pso.train"]] <- pso.train
  phyloseq_objs.train_test[["pso.test"]] <- pso.test
  phyloseq_objs.train_test[["metadata.train"]] <- metadata.train
  phyloseq_objs.train_test[["metadata.test"]] <- metadata.test
  
  save(phyloseq_objs.train_test, list = c("phyloseq_objs.train_test"), file = filepath.output)
  print_process_message("create_train_test_objects()", paste0("Train/test objects output to ", filepath.output))
  
  return(phyloseq_objs.train_test)
}