#--------------------------------------------------------------
#Function to predict 90 cell death values of the sparse matrix
cc_predict <- function(models_list,
                       norm_data
                       
  ){
  
  
  #Input requirements:
  #models_list---------------------The 'models_list.RDS' file output from cc_buildModel
  #norm_data-----------------------Output from cc_norm
  #.
  #.
  #.
  #Note: This function will return norm_data, just with 'mean_perc_cell_death_mat' and 'ref_df' updated to reflect predicted values
  

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Initialize empty list to store the data in preparation for the prediction steps
  prepped_list <- list()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Prepare data for prediction steps
  #Loop along the (sparse) normalized data
  for(i in seq_along(norm_data)){
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Extract info we will use to prepare the data
    d1_name <- norm_data[[i]]$d1_name
    d2_name <- norm_data[[i]]$d2_name
    sample_name <- norm_data[[i]]$sample_name
    current_entry_name <- names(norm_data[i])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Generate the prepared data in a format suitable for the ML prediction steps
    #This should be 30 columns (corresponding to the 30 measured/input squares), and 1 row (corresponding to the entry_name)
    prepped_df <- 
      bind_rows(
        
        #First extract the diagonal squares
        #There should only be 10 rows in this data frame
        norm_data[[i]]$ref_df %>%
          filter(plate_type=="combination") %>% #This filters to just combination data
          filter(replicate==1)  %>%#There should only be 1 replicate of the diagonal squares
          arrange(nth_conc) %>%  #arrange concentrations high-to-low, top-to-bottom
          mutate("label"="diag") %>%
          #Add the respective row and col indices of the corresponding matrix
          #Don't forget that concs are sorted high-to-low top-to-bottom
          mutate("row"=c(1,2, 3,4,5,6,7,8,9,10),
                 "col"=c(11,10,9,8,7,6,5,4,3,2)) ,
        
        #Next extract the drug1only squares
        #There should only be 10 rows in this data frame
        norm_data[[i]]$ref_df %>%
          filter(drug1_name==d1_name & plate_type=="single_agent") %>%
          filter(replicate==1) %>% #This is fine since we're selecting 'perc_cell_death' which is the mean across all replicates
          arrange(nth_conc) %>%  #arrange concentrations high-to-low, top-to-bottom
          mutate("label"="drug1only") %>%
          #Add the respective row and col indices of the corresponding 10x10 matrix
          #(Don't forget concs are sorted high-to-low, top-to-bottom)
          mutate("row"=c(11,11,11,11,11,11,11,11,11,11),
                 "col"=c(11,10,9,8,7,6,5,4,3,2)) ,
        
        #Next extract the drug2only squares
        #There should only be 10 rows in this data frame
        norm_data[[i]]$ref_df %>%
          filter(drug1_name==d2_name & plate_type=="single_agent") %>% #Note the column names here! Drug2 as a single-agent is in the 'drug1' column
          filter(replicate==1) %>% #This is fine since we're selecting 'perc_cell_death' which is the mean across all replicates
          arrange(nth_conc) %>%  #arrange concentrations high-to-low, top-to-bottom
          mutate("label"="drug2only") %>%
          #Add the respective row and col indices of the corresponding 10x10 matrix
          #(Don't forget concs are sorted high-to-low, top-to-bottom)
          mutate("row"=c(1,2,3,4,5,6,7,8,9,10),
                 "col"=c(1,1,1,1,1,1,1,1,1,1))
        
        
      ) %>%
      mutate("ind"=paste0("ind_",row,".",col,"_",label)) %>%
      select(ind, perc_cell_death) %>%
      t() %>%
      as.data.frame() %>%
      #Use the next 2 lines to make the 'ind' column the row names
      {setNames(., unlist(.[1, ]))} %>%
      .[-1, ] %>% #Remove the first row
      #Set rownames to sample_name:
      'rownames<-'(current_entry_name) %>%
      mutate(across(.cols = everything(),
                    .fns = as.numeric))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Deposit result into prepped_list
    #After the loop completes, prepped_list will have 90 elements - one for each square we're trying to predict, aka outcome variable
    prepped_list[[current_entry_name]] <- prepped_df
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
  }#End loop over norm_data
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Merge the prepped data into a single data frame
  #This will contain 30 columns (corresponding to the 30 measured/input squares), and n rows (corresponding to the number of samples/combinations)
  #The strategy is to do the predictions on the combined data once, then split them up into their respective matrices
  merged_prepped_df <- 
    do.call(rbind, prepped_list)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Initialize an empty list for the predicted results to go into
  predicted_list <- list()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Predict 90 values of %Cell death (using the appropriate 90 ML models)
  
  #Loop along the (bundled, serialized) models
  for(i in seq_along(models_list)){
    
    #Define the current index label (e.g. 'ind_10.11_predictor')
    current_index <- names(models_list[i])
    
    #Do the prediction and deposit result into appropriate element of predicted_list
    predicted_list[[current_index]] <- 
      models_list[[i]] %>% 
      unbundle() %>%  #De-serialize the model
      predict(., merged_prepped_df) %>%
      as.data.frame()
    
    #Adjust the col/row names of the predicted data frame
    colnames(predicted_list[[i]]) <- current_index
    rownames(predicted_list[[i]]) <- rownames(merged_prepped_df)
    
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Print the progress (iterations over the norm_data loop)
    percent_complete <- (i/length(models_list) * 100)
    cat(sprintf("\rProgress: %.2f%% complete", percent_complete))
    flush.console()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Combine the observed data (30 columns) and predicted data (90 columns)
  #There will be 120 columns and each row is still a sample/combination
  predicted_data <- 
    cbind(
      merged_prepped_df,
      do.call(cbind,predicted_list)
    )
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Loop along the normalized data and add the newly predicted values to each respective element
  for(i in seq_along(norm_data)){
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Get the current entry name of norm_data
    current_entry_name <- names(norm_data[i])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Subset the predicted data to include rows that match names(norm_data[i])
    predicted_data_subset <- 
      predicted_data %>%
      filter(rownames(.) %in% names(norm_data[i]))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Extract out the mean %Cell death matrix
    #Note: this contains 30 observed values as an 11x11 matrix
    #The values are means of any collected replicates
    mean_perc_cell_death_mat <- norm_data[[i]]$mean_perc_cell_death_mat
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Reformat the data to have the observed and predicted values together, the drug names/concs, and the row/col indices
    predicted_data_subset_clean <- 
      predicted_data_subset %>%
      as.matrix() %>%
      t() %>%
      as.data.frame() %>%
      mutate(ind=rownames(.)) %>%
      mutate(row = str_split(ind, "[.]") %>%
               do.call(rbind, .) %>%
               as.data.frame() %>%
               select(V1) %>%
               mutate(row= as.numeric(gsub(".*_", "", V1))) %>%
               select(-V1)) %>%
      mutate(col = str_split(ind, "[.]") %>%
               do.call(rbind, .) %>%
               as.data.frame() %>%
               select(V2) %>%
               mutate(row= as.numeric(sub("_[^_]+$", "", V2))) %>%
               select(-V2)) %>%
      mutate(across(.cols = everything(),
                    .fns = unlist)) %>%
      rename("perc_cell_death"=1) %>% #Rename first column as 'perc_cell_death'
      mutate("drug1_conc"=as.numeric(colnames(mean_perc_cell_death_mat)[col]),
             "drug2_conc"=as.numeric(rownames(mean_perc_cell_death_mat)[row]),
             "drug1_name"=norm_data[[i]]$d1_name,
             "drug2_name"=norm_data[[i]]$d2_name,
             "label"=str_extract(ind, "(?<=_)[^_]+$"))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Make the new mean %Cell Death matrix populated with the predicted values
    
    #Initialize an empty matrix
    new_mat <- 
      matrix(NA, nrow = 11, ncol = 11) %>%
      as.data.frame()
    
    #Populate the matrix with %Cell death values (this will fill in ALL values, both observe and predicted to form a complete matrix)
    for(j in 1:nrow(predicted_data_subset_clean)){
      
      row_val <- predicted_data_subset_clean[j,c("row")]
      col_val <- predicted_data_subset_clean[j,c("col")]
      pcd_val <- predicted_data_subset_clean[j,c("perc_cell_death")]
      
      new_mat[row_val, col_val] <- pcd_val
      
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Update the 0x0 index to be a value of 0
    new_mat[11,1] <- 0
    
    #Update row/col names to match original mean_perc_cell_death_mat
    rownames(new_mat) <- rownames(mean_perc_cell_death_mat)
    colnames(new_mat) <- colnames(mean_perc_cell_death_mat)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    

    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Extract necessary elements for next section (merging predicted data back to ref_df)
    ref_df <- norm_data[[i]]$ref_df
    
    d1_name <- norm_data[[i]]$d1_name
    d2_name <- norm_data[[i]]$d2_name
    sample_name <- norm_data[[i]]$sample_name
    units <- norm_data[[i]]$units
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Joining predicted results back to ref_df
    #The goal is to join ONLY the predicted data back to the observed data that already exists in ref_df
    #The reason is because the observed data in ref_df has a lot more info associated with it (plate barcode, zprime, etc), we don't want to lose
    
    
    #From the subsetted clean predicted data, extract only the values we predicted (no observed values)
    mydf_tmp <- 
    predicted_data_subset_clean  %>%
      filter(label=="predictor") %>% #Only keep the predicted values
      select(drug1_name, 
             drug1_conc,
             drug2_name,
             drug2_conc,
             perc_cell_death) %>%
      mutate("units"=units,
             "sample_name"=sample_name,
             "plate_type"="combination",
             "combo_id"=paste0(drug1_name, "_", drug2_name),
             "replicate"=1) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Before joining back to ref_df, we need to account for the columns that exist in ref_df but don't exist in mydf_tmp
    #Figure out the missing columns and add them to the mydf_tmp (as all NA's)
    
    new_columns <- setdiff(colnames(ref_df),
                           colnames(mydf_tmp))
    
    new_columns_df <- as.data.frame(matrix(ncol = length(new_columns), nrow = nrow(mydf_tmp)))
    
    colnames(new_columns_df) <- new_columns
    
    #Add the new columns back to mydf_tmp
    #There should still be 90 rows, corresponding to the 90 predicted values
    mydf_tmp <- 
      cbind(mydf_tmp,
            new_columns_df)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Bind the rows of mydf_tmp (predicted values) back to the observed values in ref_df
    #This should have the number of rows of the original (sparse) ref_df + 90
    ref_df_new <- 
      rbind(ref_df, mydf_tmp)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Now replace the mean_perc_cell_death_mat and ref_df with the updated version that contain the predicted values:
    norm_data[[i]]$mean_perc_cell_death_mat <- new_mat
    norm_data[[i]]$ref_df <- ref_df_new
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
   
    

  }#End loop over norm_data
  
  return(norm_data)
  
}#End function
