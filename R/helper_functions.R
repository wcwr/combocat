#Helper functions that may be used across multiple combocat functions:

#--------------------------------------------------------------
#Function to drop outliers within controls based on z-score
drop_outliers <- function(vector, threshold){
  mean_val <- mean(vector)
  std_val  <- sd(vector)
  z_scores <- (vector - mean_val) / std_val
  outliers <- vector[abs(z_scores) > threshold]
  
  cleaned_vector <- subset(vector, !(vector %in% outliers))
  return(cleaned_vector)
}
#--------------------------------------------------------------





#--------------------------------------------------------------
#Function to calculate Z'
zprime_func <- function(pos_ctrl, neg_ctrl){
  1-(3*(sd(pos_ctrl, na.rm = TRUE)) + 3*(sd(neg_ctrl, na.rm = TRUE))) / 
    abs(mean(pos_ctrl, na.rm=TRUE) - mean(neg_ctrl, na.rm=TRUE))
}
#--------------------------------------------------------------





#--------------------------------------------------------------
#Function to calculate robust Z'
robust_zprime_func <- function(pos_ctrl, neg_ctrl){
  1-3*(mad(pos_ctrl, na.rm = TRUE)+mad(neg_ctrl, na.rm = TRUE))/abs(median(pos_ctrl, na.rm = TRUE)-median(neg_ctrl, na.rm = TRUE))
}
#--------------------------------------------------------------





#--------------------------------------------------------------
#Function to perform normalization
norm_func <- function(raw_vals, neg_ctrl_vals, pos_ctrl_vals){
  100*((raw_vals-mean(neg_ctrl_vals, na.rm=TRUE))/(mean(pos_ctrl_vals, na.rm=TRUE) - mean(neg_ctrl_vals, na.rm=TRUE)))
}
#--------------------------------------------------------------





#--------------------------------------------------------------
#Function to fit a dose-response model to a dataframe
dr_func <- function(input_dataframe, #Dataframe with a column for dose and a column for response
                    dose_column      = "drug_conc",       #Name of dose  column
                    response_column  = "perc_cell_death", #Name of resp. column (should be the mean of replicates)
                    drug_name_column = "drug_name"        #Name of the drug column
){ 
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Initialize variables that will be returned
  dr_model <- NULL
  pred_df <- NULL
  model_success <- FALSE
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Make formula for the following drm step, using the dose and response column names
  my_formula <- as.formula(paste(response_column, "~", dose_column))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Dose-response modeling section
  
  #Use a tryCatch statement to handle what happens if model cannot be fit
  tryCatch({
    
    #================================
    #If model CAN be successfully fit
    #================================
    
    #Generate the dose-response model
    dr_model <- drm(my_formula,
                    data = input_dataframe,
                    fct = LL.4())
    
    #Expand the concentrations
    #Importantly, we expand on a log-scale and then transform back to linear scale
    #This will fit points evenly between the concentrations once plotted on a log scale
    concs_expanded <- 
      seq(log10(min(input_dataframe[[dose_column]])),
          log10(max(input_dataframe[[dose_column]])),
          length.out = 1000) %>%
      10^.
    
    #Predict the $Cell Death response using the model
    pred_df <-
      data.frame("dose" = concs_expanded,
                 "response"  = predict(dr_model, 
                                       newdata = data.frame(dose_column = concs_expanded)))
    
    #Report if model was successfully fit
    model_success <- TRUE
    
    
    
  }, error = function(e){
    
    #================================
    #If model CAN NOT be successfully fit
    #================================
    
    #Print a message detailing which drug failed
    message("Model fitting failed with message: ", 
            e$message, 
            " (drug is ", 
            input_dataframe[[drug_name_column]] %>% unique(),
            ")")
    
    #Set dr_model to NULL
    dr_model <- NULL
    
    #Expand the concentrations
    #Importantly, we expand on a log-scale and then transform back to linear scale
    #This will fit points evenly between the concentrations once plotted on a log scale
    concs_expanded <- 
      seq(log10(min(input_dataframe[[dose_column]])),
          log10(max(input_dataframe[[dose_column]])),
          length.out = 1000) %>%
      10^.
    
    #Just list 0 as the response for pred_df if model cannot be fit
    pred_df <-
      data.frame("dose" = concs_expanded,
                 "response"  = 0)
    
    #Report if model was successfully fit
    model_success <- FALSE
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  return(list(dr_model = dr_model, 
              pred_df = pred_df, 
              model_success = model_success))
  
  
}#End dr_func
#--------------------------------------------------------------





#--------------------------------------------------------------
#Function to calculate position_id from plate_row and plate_col (for sparse mode plates)
calculate_position_id <- function(plate_row, plate_col) {
  # Calculate row group based on plate_row
  row_group <- (plate_row - 2) %% 3 + 1
  
  # Calculate column group based on plate_col
  col_group <- (plate_col - 2) * 3
  
  # Calculate the actual position_id
  position_id <- row_group + col_group
  return(position_id)
}
#--------------------------------------------------------------





#--------------------------------------------------------------
#Function to calculate missing concentration for a specified drug (used in cc_makeMeta)
#When used in a loop over any missing concentrations, this function will fill in the missing concentration values
#NOTE: The dilution factor is provided by user, default is 3 (1:3 dilution)
#This should technically work even if only 1 concentration were provided
#--------------------------------------------------------------
calculate_concentration <- function(nth, df, conc_column, dilution_factor) {
  # Check if the missing concentration is the highest concentration (nth == 1)
  if (nth == 1) {
    # Start looking for the next available concentration
    next_nth <- nth + 1
    # Loop to find the next available concentration if the immediate next one is missing
    while (is.na(df[[conc_column]][df$nth_conc == next_nth]) && next_nth <= 10) {
      next_nth <- next_nth + 1
    }
    if (next_nth <= 10) {
      # Calculate the missing concentration by multiplying the next available concentration
      # by dilution_factor^(next_nth - nth)
      calculated_conc <- df[[conc_column]][df$nth_conc == next_nth] * 
        (dilution_factor ^ (next_nth - nth))
      
      # ROUND the calculated concentration to 6 decimal places
      return(round(calculated_conc, digits = 6))
    } else {
      return(NA)
    }
  } else {
    # If the missing concentration is not the highest, start looking for the previous available concentration
    prev_nth <- nth - 1
    while (is.na(df[[conc_column]][df$nth_conc == prev_nth]) && prev_nth >= 1) {
      prev_nth <- prev_nth - 1
    }
    if (prev_nth >= 1) {
      # Calculate the missing concentration by dividing the previous available concentration
      # by dilution_factor^(nth - prev_nth)
      calculated_conc <- df[[conc_column]][df$nth_conc == prev_nth] / 
        (dilution_factor ^ (nth - prev_nth))
      
      # ROUND the calculated concentration to 6 decimal places
      return(round(calculated_conc, digits = 6))
    } else {
      return(NA)
    }
  }
}
#--------------------------------------------------------------






#--------------------------------------------------------------
#Function to calculate total # of plates for a sparse mode experiment
#(Assumes all-by-all screen)
calc_sparse_plates <- function(num_cell_lines,
                                          num_drugs,
                                          num_sa_reps){
  
  #Requirements:
  #num_cell_lines............The number of total cell lines to be used
  #num_drugs.................The total number of unique drugs (up to 135)
  #num_sa_reps...............The number of replicates of each SINGLE-AGENT plate (combo plates are not repeated)
  
  # #Make sure the num_drugs does not exceed 135
  if(num_drugs>135){print("number of drugs cannot exceed 135")}
  stopifnot(all(num_drugs <= 135))

  first_term  <- num_sa_reps*ceiling((num_drugs/135))
  second_term <- ceiling((choose(num_drugs,2))/135)
  
  result <- num_cell_lines*(first_term + second_term)
  
  return(result)
}
#--------------------------------------------------------------





#--------------------------------------------------------------
#Function to determine the monotonicity of a single-agent 
calc_monotonicity <- function(input_df){
  
  #Background: 
  #Data points not being exactly monotonically increasing aren't always a problem for DR curves 
  #But these points can be problematic for synergy scoring.
  #For a given data point (%CD), we determine how far up or down it is compared to its preceding point. Ideally, it always goes up (when evaluated low-to-high)
  #Or at least, it doesn't go down by too much (some wiggle room can be set with the cutoff_mono parameter in cc_getQC)
  #This function tests for monotonicity in both directions (low-to-high, and high-to-low)
  #When tested low-to-high, the datas hould be monotonically increasing, and when tested high-to-low, the data should be monotonically decreasing
  
  #Requirements:
  #input_df..................A dataframe with columns "drug_conc" &  "perc_cell_death" (should be 10 rows)
  
  
  #If input_df has > 10 rows, stop and render a warning:
  if(nrow(input_df) > 10){
    stop("Input dataframe must have exactly 10 rows")
  }
  
  
  #When evaluating concentrations low-to-high: For a given data point get the preceding data point and calculate the difference
  res_low_to_high <- 
    input_df %>%
    arrange(drug_conc) %>% #Arrange the data as low-to-high concentrations
    mutate("prev_val_low_to_high"=dplyr::lag(perc_cell_death)) %>% #Grab the previous %CD value at any given concentration
    mutate("delta_val_low_to_high" = perc_cell_death - prev_val_low_to_high)  #Calculate the difference between the current and previous %CD values
  
  #When evaluating concentrations high-to-low: For a given data point get the preceding data point and calculate the difference
  res_high_to_low <- 
    input_df %>%
    arrange(-drug_conc) %>% #Arrange the data as high-to-low concentrations 
    mutate("prev_val_high_to_low"=dplyr::lag(perc_cell_death)) %>% #Grab the previous %CD value at any given concentration
    mutate("delta_val_high_to_low" = perc_cell_death - prev_val_high_to_low) #Calculate the difference between the current and previous %CD values
  
  
  res_final <- 
    merge(res_low_to_high,
          res_high_to_low,
          by=c("drug_conc", "perc_cell_death"))
  
  return(res_final)
  
}
#--------------------------------------------------------------





#--------------------------------------------------------------
#Function to calculate Moran's I for spatial autocorrelation of a synergy matrix
calc_moran <- function(input_syn_mat){
  
  #Background:
  #Moran's I is a measure of spatial autocorrelation in a matrix. 
  #This can be useful in helping rank synergy matrices, as seemingly spatially autocorrelated ones seem to be quite believable/validated
  #Moran's I ranges from -1 (totally dispersed values) to 1 (totally clustered values), with 0 indicating random distribution
  #NOTE: This function calculates the non-standardized AND standardized Moran's I (details below)
  
  #Requirements:
  #input_syn_mat..............A square matrix of synergy values where the column names are concs of drug1, and row names are concs of drug2
  #...........................For the purposes of Combocat, this should be 10x10
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #The formula for global Moran's I is:
  #I = (n / W) * (sum(i=1 to n) sum(j=1 to n) w_ij * (x_i - x_bar) * (x_j - x_bar)) / (sum(i=1 to n) (x_i - x_bar)^2)
  
  #Where:
  #n = number of observations
  #w_ij = spatial weight between observations i and j
  #W = The sum of all spatial weights w_ij
  #x_i = value of observation i
  #x_bar = mean of all observations (values in the synergy matrix)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Terms easy to calculate from the input matrix
  n           <- length(input_syn_mat)        #Number of total observations/values in the matrix
  xbar        <- mean(input_syn_mat)          #Mean of all observations
  denom_value <- sum((input_syn_mat-xbar)^2)  #Denominator value for the Moran's I formula
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Prepare the spatial weights matrix w_ij
  
  #Initialize the matrix (all zeros to start)
  w_ij <- matrix(0, nrow = n, ncol = n)
  
  #Row/column labels for the weights matrix (for interpretability)
  row_col_labels <- apply(expand.grid(1:sqrt(n), 1:sqrt(n), KEEP.OUT.ATTRS = FALSE)[, 2:1], 1, function(x) paste0("[", x[1], ",", x[2], "]"))
  
  #Apply these labels to the w_ij matrix
  rownames(w_ij) <- row_col_labels
  colnames(w_ij) <- row_col_labels
  
  #Define the matrix size (in this case: the number of rows or columns of the input matrix)
  matrix_size <- sqrt(n)
  
  #Define the neighbor offsets (Queen's contiguity)
  #This will tell us how far to move in the matrix to find a neighbor
  #Note: Queen's contiguitiy means that cells are neighbors if they share an edge or a corner
  neighbor_offsets <- expand.grid(row_offset = -1:1, col_offset = -1:1) #All possible positions surriounding a given matrix cell
  neighbor_offsets <- neighbor_offsets[!(neighbor_offsets$row_offset == 0 & neighbor_offsets$col_offset == 0), ] #Filters out [0,0] which is the cell itself
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Loop over each cell of the original matrix to find its neighbors and update w_ij
  for (i in 1:n) {
    # Calculate the current cell's row and column in the original matrix
    # The %/% operator returns the integer division (i.e., how many times matrix_size fits in (i-1))
    current_row <- (i - 1) %/% matrix_size + 1
    # The %% operator returns the remainder of division (i.e., which column in the row)
    current_col <- (i - 1) %% matrix_size + 1
    
    # Loop over each possible neighbor offset (Queen's contiguity)
    for (offset in 1:nrow(neighbor_offsets)) {
      # Calculate the row and column indices of the potential neighbor
      new_row <- current_row + neighbor_offsets$row_offset[offset]
      new_col <- current_col + neighbor_offsets$col_offset[offset]
      
      # Check if the new row and column are within the bounds of the matrix
      if (new_row >= 1 && new_row <= matrix_size && new_col >= 1 && new_col <= matrix_size) {
        # Calculate the linear index of the neighbor cell
        neighbor_index <- (new_row - 1) * matrix_size + new_col
        # Set the weight in w_ij to 1, indicating that these two cells are neighbors
        w_ij[i, neighbor_index] <- 1
      }
    }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Generate a dataframe that tells us the contributions of each cell to the Moran's I numerator term
  contributions_df <- data.frame(
    contributing_cell = character(),
    neighbor_cell = character(),
    contribution_value = numeric(),
    total_index_contribution = numeric()
  )
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Loop through each row of w_ij to calculate its numerator contribution
  for (i in rownames(w_ij)) {
    
    #Return the current row values that are greater than zero
    non_zero_columns <- w_ij[i, w_ij[i, ] > 0]
    
    #Extract the names of these columns
    non_zero_column_names <- colnames(w_ij)[w_ij[i, ] > 0]
    
    #For each of the non-zero columns, calculate the contribution to the numerator
    for (k in non_zero_column_names) {
      
      # Extract numeric indices for `i` and `k` by parsing their names
      i_indices <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", i), ",")))
      k_indices <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", k), ",")))
      
      # Access the values in the `input_syn_mat`
      x_i <- input_syn_mat[i_indices[1], i_indices[2]]
      x_k <- input_syn_mat[k_indices[1], k_indices[2]]
      
      #Multiply the weight times (x_i - x_bar) * (x_k - x_bar)
      numerator_contribution <- w_ij[i, k] * (x_i - xbar) * (x_k - xbar)
      
      #Append the contribution to the data frame
      contributions_df <- rbind(
        contributions_df, 
        data.frame(
          contributing_cell = i,
          neighbor_cell = k,
          contribution_value = numerator_contribution,
          stringsAsFactors = FALSE
        ))
    }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Calculate the total index contribution for each cell
  #This is not strictly necessary for the Moran's I calculation, but it can be useful for interpretation
  contributions_df <- 
    contributions_df %>%
    group_by(contributing_cell) %>%
    mutate(total_index_contribution = sum(contribution_value))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Get W, the sum of w_ij
  W <- sum(w_ij)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Calculate the final numerator value
  numerator_value <- sum(contributions_df$contribution_value)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Plug in values to get Moran's I (NON-standardized version)
  moran_i <- (n / (W * denom_value)) * numerator_value
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #====================================
  #======STANDARDIZED MORAN'S I========
  #====================================
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Standardizing the w_ij matrix aims to control for the influence of cells with more/fewer neighbors than other cells
  #e.g., a cell in the center of the matrix will have more neighbors than a cell on the edge/corner
  #This is done by dividing each row of w_ij by the sum of the row (so that all rows sum to 1)
  
  #We start from the w_ij matrix calculated above, then standardize, then repeat the remaining steps:
  w_ij_STND <- w_ij #Make a copy of w_ij which gets standardized
  w_ij_STND <- sweep(w_ij_STND, 1, rowSums(w_ij_STND), "/") #This is the standardizing part. All rows should now sum to 1
  
  W_STND <- sum(w_ij_STND) #Get the sum of the standardized w_ij matrix
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  contributions_df_STND <- data.frame(
    contributing_cell = character(),
    neighbor_cell = character(),
    contribution_value = numeric(),
    total_index_contribution = numeric()
  )
  
  for (i in rownames(w_ij_STND)) {
    non_zero_columns <- w_ij_STND[i, w_ij_STND[i, ] > 0]
    non_zero_column_names <- colnames(w_ij_STND)[w_ij_STND[i, ] > 0]
    for (k in non_zero_column_names) {
      i_indices <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", i), ",")))
      k_indices <- as.numeric(unlist(strsplit(gsub("\\[|\\]", "", k), ",")))
      x_i <- input_syn_mat[i_indices[1], i_indices[2]]
      x_k <- input_syn_mat[k_indices[1], k_indices[2]]
      numerator_contribution <- w_ij_STND[i, k] * (x_i - xbar) * (x_k - xbar)
      
      contributions_df_STND <- rbind(
        contributions_df_STND, 
        data.frame(
          contributing_cell = i,
          neighbor_cell = k,
          contribution_value = numerator_contribution,
          stringsAsFactors = FALSE
        ))
    }
  }
  
  numerator_value_STND <- sum(contributions_df_STND$contribution_value)
  
  moran_i_STND <- (n / (W_STND * denom_value)) * numerator_value_STND
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  return(list(
    "moran_i"=moran_i,
    "moran_i_STND"=moran_i_STND
  ))
  
}
#--------------------------------------------------------------