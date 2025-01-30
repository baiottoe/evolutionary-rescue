
# Define the function to run a single iteration of the slim spinup model
run_slim_normal <- function(parameters) {
  tryCatch({
    ENV_SEED <- parameters$envSeed
    SLiM_SEED <- parameters$slimSeed
    HL_SEED <- parameters$hlSeed
    hlAcLen <- parameters$acLen
    delta_e <- parameters$delta_e
    sigma_p <- parameters$sigma_p
    lambda_o <- parameters$lambda_o
    sigma_f <- parameters$sigma_f
    sigma_QTL <- parameters$sigma_QTL
    K <- parameters$K
    output_suffix <- parameters$output_suffix
    spinup_suffix <- parameters$spinup_suffix
    l_e <- parameters$l_e
    p_QTL <- parameters$p_QTL
    p_hl <- parameters$p_hl
    
    # Model directory/name
    normal_model_path <- file.path(getwd(), 'script/slim/continuous_n.slim')
    
    # Function to create the command with specified HLDim
    run_slim <- paste0(
      paste0("slim -d ", "'", 'envSeed="', ENV_SEED, '"', "' "),
      paste0("-d slimSeed=", SLiM_SEED, " "),
      paste0("-d hlSeed=", HL_SEED, " "),
      paste0("-d hlAclen=", hlAcLen, " "),
      paste0("-d ", "'", 'CURRENT_DIRECTORY="', getwd(), '"', "' "),
      paste0("-d ", "'", 'SPINUP_SUFFIX="', spinup_suffix, '"', "' "),
      paste0("-d ", "'", 'OUTPUT_SUFFIX="', output_suffix, '"', "' "),
      paste0("-d delta_env=", delta_e, " "),
      paste0("-d sigma_p=", sigma_p, " "),
      paste0("-d lambda_o=", lambda_o, " "),
      paste0("-d sigma_fitness=", sigma_f, " "),
      paste0("-d sigma_qtl=", sigma_QTL, " "),
      paste0("-d carrying_cap=", K, " "),
      paste0("-d envAcLen=", l_e, " "),
      paste0("-d pQTL=", p_QTL, " "),
      paste0("-d pHL=", p_hl, " "),
      paste0("-d HLDim=", 256, " "),
      normal_model_path, " 2>&1"
    )
    
    # First try with HLDim=256
    system(run_slim, intern = TRUE)
    
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error in run for envSeed=%s, hlSeed=%s, slimSeed=%s: %s", 
                    ENV_SEED, HL_SEED, SLiM_SEED, e$message))
    return(FALSE)
  })
}

run_slim_spinup <- function(parameters) {
  tryCatch({
    ENV_SEED <- parameters$envSeed
    SLiM_SEED <- parameters$slimSeed
    sigma_p <- parameters$sigma_p
    lambda_o <- parameters$lambda_o
    sigma_f <- parameters$sigma_f
    sigma_QTL <- parameters$sigma_QTL
    K <- parameters$K
    spinup_suffix <- parameters$spinup_suffix
    l_e <- parameters$l_e
    p_QTL <- parameters$p_QTL
    
    # Define the value of the variable
    spinup_model_path <- file.path(getwd(), 'script/slim/continuous_s.slim')
    
    # Run slim spinup
    run_slim_spinup <- paste0(
      paste0("nohup slim -d ", "'", 'envSeed="', ENV_SEED, '"', "' "),
      paste0("-d slimSeed=", SLiM_SEED, " "),
      paste0("-d ", "'", 'CURRENT_DIRECTORY="', getwd(), '"', "' "),
      paste0("-d ", "'", 'SPINUP_SUFFIX="', spinup_suffix, '"', "' "),
      paste0("-d sigma_p=", sigma_p, " "),
      paste0("-d lambda_o=", lambda_o, " "),
      paste0("-d sigma_fitness=", sigma_f, " "),
      paste0("-d sigma_qtl=", sigma_QTL, " "),
      paste0("-d carrying_cap=", K, " "),
      paste0("-d envAcLen=", l_e, " "),
      paste0("-d pQTL=", p_QTL, " "),
      spinup_model_path, " > /dev/null"
    )
    
    # Run spinup
    system(run_slim_spinup, intern = TRUE)
    
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error in run for envSeed=%s, slimSeed=%s: %s", 
                    ENV_SEED, SLiM_SEED, e$message))
    return(FALSE)
  })
}

# Function to chunk the data and process in batches
process_in_batches <- function(data, batch_size, slim_function, num_processes) {
  total_batches <- ceiling(nrow(data) / batch_size)
  
  for (batch_num in 1:total_batches) {
    message(sprintf("Processing batch %d of %d", batch_num, total_batches))
    
    # Get the current batch
    start_idx <- ((batch_num - 1) * batch_size) + 1
    end_idx <- min(batch_num * batch_size, nrow(data))
    current_batch <- data[start_idx:end_idx, ]
    
    # Set up parallel processing for this batch
    plan(multisession, workers = num_processes)
    
    # Process the batch
    results <- future_lapply(1:nrow(current_batch), function(i) {
      slim_function(current_batch[i, ])
    }, future.seed = TRUE)
    
    # Clean up parallel workers
    plan(sequential)
    gc()
    
    # Save progress
    successful <- sum(unlist(results))
    failed <- length(results) - successful
    message(sprintf("Batch complete: %d successful, %d failed", successful, failed))
    
    # Optional: add a small delay between batches
    Sys.sleep(1)
  }
}