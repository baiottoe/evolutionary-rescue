
get_emergent_properties <- function(PA_df, column_name, dispersal_power_constant, Y_CLINE_WIDTH) {
  for (i in 1:nrow(PA_df)) {
    #######################################################################
    # Setup #
    #######################################################################
    # get this iteration of PA design as vector
    PA_design <- PA_df[i, column_name]
    side_length <- sqrt(nchar(PA_design)) # since we are dealing with square metapopulations
    
    # Only for taxa-specific
    if (dispersal_power_constant){
      dispersal_power <- 1
    } else {
      dispersal_power <- PA_df[i, "dispersal_power"]
    }
    
    # Convert the PA vector to a matrix
    PA_matrix <- matrix(as.numeric(strsplit(PA_design, "")[[1]]), ncol = side_length, byrow = TRUE)
    
    #######################################################################
    # Emergent properties 1+2: mean and variance of opt phenotypes of PAs #
    #######################################################################
    # Create matrix of optimum phenotype offsets from mean
    opt_offsets_as_vector <- rep(seq(from=0-(Y_CLINE_WIDTH/2), to=0+(Y_CLINE_WIDTH/2), length.out=dim(PA_matrix)[1]), dim(PA_matrix)[1])
    
    opt_offsets_as_matrix <- matrix(opt_offsets_as_vector, ncol = side_length, byrow = FALSE)
    
    opt_offsets_as_matrix[PA_matrix == 0] <- NA
    
    # Set value for this PA design in df
    PA_df[i, "mean_opt"] <- mean(opt_offsets_as_matrix, na.rm=TRUE)
    PA_df[i, "var_opt"] <- var(as.vector(opt_offsets_as_matrix), na.rm=TRUE)
    
    #######################################################################
    ######### Emergent property 3:Number of 4-connectivity patches ########
    #######################################################################
    
    # Placeholder for final value of # of patches
    num_patches <- 0
    
    # Matrix of patch IDs
    patch_IDs_4con <- matrix(0, nrow = side_length, ncol = side_length)
    
    # Create list to store all cell indexes
    to_search <- c()
    
    # Fill for all 
    for (col in 1:dim(PA_matrix)[1]) {
      for (row in 1:dim(PA_matrix)[2]) {
        to_search[[paste(col, row, sep = "_")]] <- c(col, row)
      }
    }
    
    # Algorithm to find total # patches
    while(length(to_search) > 0) {
      # Randomly sample patch index from list of all indeces
      patch_list_index <- sample(length(to_search), 1)
      
      # Get patch location and remove from full list
      this_patch <- to_search[patch_list_index]
      to_search <- to_search[-patch_list_index]
      
      # Skip if this patch is not protected
      if (PA_matrix[this_patch[[1]][1], this_patch[[1]][2]] == 0) {
        this_patch <- c()
      }
      
      # If this is part of a patch, add 1 to the number of patches
      if (length(this_patch) > 0){
        num_patches <-num_patches + 1
      }
      
      # find all features of this patch
      while(length(this_patch) > 0) {
        col_index <- this_patch[[1]][2]
        row_index <- this_patch[[1]][1]
        
        # Add index of patch to the cell we are currently on
        patch_IDs_4con[row_index, col_index] <- num_patches
        
        # Create matrix of surrounding 4 cell's indeces
        surrounding <- matrix(c(row_index, col_index-1, row_index, col_index+1, row_index-1, col_index, row_index+1, col_index), ncol = 2, byrow = TRUE)
        
        # Have to replace 10s with 1 and 0s with 9 to represent torus
        surrounding[surrounding == 10] <- 1
        surrounding[surrounding == 0] <- 9
        
        # For each of the surrounding cells:
        for (row_num in 1:nrow(surrounding)) {
          # Get indeces
          row_index_search <- surrounding[row_num, ][1]
          col_index_search <- surrounding[row_num, ][2]
          
          # Add to this patch for further searching if is a PA, remove from further consideration if not
          if ((paste(row_index_search, col_index_search, sep = "_") %in% names(to_search)) && PA_matrix[row_index_search, col_index_search] == 1) {
            # Add to this patch to futher search around edges
            this_patch <- c(this_patch, to_search[paste(row_index_search, col_index_search, sep = "_")])
          }
          # Remove this patch from to_search
          to_search <- to_search[!names(to_search) %in% paste(row_index_search, col_index_search, sep = "_")]
        }
        
        # Also remove the patch we just searched around from this patch, since we have already checked around it
        this_patch <- this_patch[!names(this_patch) %in% paste(row_index, col_index, sep = "_")]
      }
    }
    
    # Add calculated number of patches to df
    PA_df[i, "number_patches_4con"] <- num_patches
    
    #######################################################################
    ######### Emergent property 4:Number of 8-connectivity patches ########
    #######################################################################
    
    # Placeholder for final value of # of patches
    num_patches <- 0
    
    # Matrix of patch IDs
    patch_IDs_8con <- matrix(0, nrow = side_length, ncol = side_length)
    
    # Create list to store all cell indexes
    to_search <- c()
    
    # Fill for all 
    for (col in 1:dim(PA_matrix)[1]) {
      for (row in 1:dim(PA_matrix)[2]) {
        to_search[[paste(col, row, sep = "_")]] <- c(col, row)
      }
    }
    
    # Algorithm to find total # patches
    while(length(to_search) > 0) {
      # Randomly sample patch index from list of all indeces
      patch_list_index <- sample(length(to_search), 1)
      
      # Get patch location and remove from full list
      this_patch <- to_search[patch_list_index]
      to_search <- to_search[-patch_list_index]
      
      # Skip if this patch is not protected
      if (PA_matrix[this_patch[[1]][1], this_patch[[1]][2]] == 0) {
        this_patch <- c()
      }
      
      # If this is part of a patch, add 1 to the number of patches
      if (length(this_patch) > 0){
        num_patches <-num_patches + 1
      }
      
      # find all features of this patch
      while(length(this_patch) > 0) {
        col_index <- this_patch[[1]][2]
        row_index <- this_patch[[1]][1]
        
        # Add index of patch to the cell we are currently on
        patch_IDs_8con[row_index, col_index] <- num_patches
        
        # Create matrix of surrounding 4 cell's indeces
        surrounding <- matrix(c(row_index, col_index-1, row_index, col_index+1, row_index-1, col_index, row_index+1, col_index, row_index-1, col_index-1, row_index-1, col_index+1, row_index+1, col_index-1, row_index+1, col_index+1), ncol = 2, byrow = TRUE)
        
        # Have to replace 10s with 1 and 0s with 9 to represent torus
        surrounding[surrounding == 10] <- 1
        surrounding[surrounding == 0] <- 9
        
        # For each of the surrounding cells:
        for (row_num in 1:nrow(surrounding)) {
          # Get indeces
          row_index_search <- surrounding[row_num, ][1]
          col_index_search <- surrounding[row_num, ][2]
          
          # Add to this patch for further searching if is a PA, remove from further consideration if not
          if ((paste(row_index_search, col_index_search, sep = "_") %in% names(to_search)) && PA_matrix[row_index_search, col_index_search] == 1) {
            # Add to this patch to futher search around edges
            this_patch <- c(this_patch, to_search[paste(row_index_search, col_index_search, sep = "_")])
          }
          # Remove this patch from to_search
          to_search <- to_search[!names(to_search) %in% paste(row_index_search, col_index_search, sep = "_")]
        }
        
        # Also remove the patch we just searched around from this patch, since we have already checked around it
        this_patch <- this_patch[!names(this_patch) %in% paste(row_index, col_index, sep = "_")]
      }
    }
    
    # Add calculated number of patches to df
    PA_df[i, "number_patches_8con"] <- num_patches
    
    #######################################################################
    # Emergent property 5:connectivity,measured as mean global centrality #
    #######################################################################
    ############################### AND ###################################
    #######################################################################
    ###### Emergent property 6:proportion of between patch dispersal ######
    #######################################################################
    # Create adjacency matrix, will have values weighted by distance
    weighted_adjacency_matrix <- matrix(NaN, nrow = side_length**2, ncol = side_length**2)
    
    # Also create a matrix that holds the proportion of dispersing individuals from that subpop that will disperse to subpops outside of their "patch"
    between_patch_dispersal_prop_4con <- matrix(NaN, nrow = side_length, ncol = side_length)
    between_patch_dispersal_prop_8con <- matrix(NaN, nrow = side_length, ncol = side_length)
    
    
    PA_list <- unlist(strsplit(PA_design, split=""))
    for (subpop in 1:side_length**2) {
      # get row and column indexes for matrix for this subpop
      row_center <- ceiling(subpop/side_length)
      col_center <- ((subpop-1) %% side_length) + 1
      
      # if subpop is protected, get distance vector to other protected subpops
      if (PA_list[subpop] == "1") {
        # get distance matrix to all subpops, then multiply by PA_matrix to exclude areas not protected
        distances <- calculate_distances_torus(matrix(0, nrow = side_length, ncol = side_length), row_center, col_center) * PA_matrix
        
        # Get matrix of dispersal probabilities
        p_dispersal <- 1/(distances**(dispersal_power))
        p_dispersal[is.infinite(p_dispersal)] <- 0 #replace infinite values with 0, as they represent patches that individuals cannot migrate to
        p_dispersal <- p_dispersal/sum(p_dispersal) # scale so sum of dispersal probabilities = 1
        
        # Get sum of probabilities for dispersing to other patches (equivalent to proportion)
        this_patch_ID_4con <- patch_IDs_4con[row_center, col_center]
        patch_proportion_4con <- sum(p_dispersal[patch_IDs_4con != this_patch_ID_4con])
        this_patch_ID_8con <- patch_IDs_8con[row_center, col_center]
        patch_proportion_8con <- sum(p_dispersal[patch_IDs_8con != this_patch_ID_8con])
        
        # Store calculated proportion of between patch dispersers in previously created matrix
        between_patch_dispersal_prop_4con[row_center, col_center] <- patch_proportion_4con
        between_patch_dispersal_prop_8con[row_center, col_center] <- patch_proportion_8con
        
        # assign this distance vector to the subpop's row in the adjacency matrix
        weighted_adjacency_matrix[subpop, ] <- as.vector(t(distances))
      }
    }
    
    
    # calculate mean global centrality
    # https://rpubs.com/odenipinedo/network-analysis-in-R
    # https://www.datacamp.com/tutorial/centrality-network-analysis-R
    mean_global_centrality <- mean(rowMeans(weighted_adjacency_matrix), na.rm=TRUE)
    
    PA_df[i, "mean_global_centrality"] <- mean_global_centrality
    
    
    # calculate proportion of between patch dispersal
    prop_between_4con <- mean(between_patch_dispersal_prop_4con, na.rm=TRUE)
    prop_between_8con <- mean(between_patch_dispersal_prop_8con, na.rm=TRUE)
    
    PA_df[i, "prop_between_4con"] <- prop_between_4con
    PA_df[i, "prop_between_8con"] <- prop_between_8con
    
  }
  
  return(PA_df)
}
