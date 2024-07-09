#--------------------------------------------------------------
#Function to normalize mapped data
cc_norm <- function(mapped_data){ 
  
  
  #Input requirements:
  #mapped_data------------------Output from cc_map
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Determine whether data is dense or sparse mode
  if(mapped_data$combo_wise_list[[1]]$data_mode == "sparse"){
    data_mode <- "sparse"
    print("Detected data type is SPARSE MODE. Proceeding with sparse mode normalization...")
  } else {
    data_mode <- "dense"
    print("Detected data type is DENSE MODE. Proceeding with dense mode normalization...")}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

  
  
  if(data_mode=="dense"){
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle dense mode data------------------------
    #
    #=================================================================
    #=================================================================
    
    
    #Initialize empty list all results will be stored in (and returned)
    mylist <- list()
    
    
    for(i in seq_along(mapped_data$combo_wise_list)){
      
      
      #Get the current plate name
      current_plate_name <- names(mapped_data$combo_wise_list)[i]
      
  
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #mapped reference dataframe
      ref_df <- mapped_data$combo_wise_list[[i]]$ref_df
      
      #Initialize empty column that normalized data will go into in ref_df
      ref_df$replicate_perc_cell_death <- NA
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Extract the elements that need to be normalized:
      
      #Drug1/2 single-agent replicates:
      rep_drug1_list_raw <- mapped_data$combo_wise_list[[i]]$rep_drug1_list_raw
      rep_drug2_list_raw <- mapped_data$combo_wise_list[[i]]$rep_drug2_list_raw
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
 
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Extract the controls that will be used for normalization
      pos_ctrl <- mapped_data$combo_wise_list[[i]]$raw_controls_list_pos[[current_plate_name]]
      neg_ctrl <- mapped_data$combo_wise_list[[i]]$raw_controls_list_neg[[current_plate_name]]
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
     
  
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Normalize the 'raw_val' column of ref_df, into new column perc_cell_death
      for(k in 1:nrow(ref_df)){
        ref_df$replicate_perc_cell_death[k] <- norm_func(raw_vals = ref_df$raw_val[k],
                                               pos_ctrl_vals = pos_ctrl,
                                               neg_ctrl_vals = neg_ctrl)
      }
      
      #When drug1_conc and drug2_conc = 0, and thus raw_val=0, set perc_cell_death to 0
      ref_df[ref_df$raw_val==0,c("replicate_perc_cell_death")] <- 0
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Get the mean %Cell Death across replicates
      ref_df <- 
        ref_df %>%
        group_by(drug1_name, drug2_name, drug1_conc, drug2_conc)  %>%
        mutate("perc_cell_death"=mean(replicate_perc_cell_death)) %>%
        ungroup()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Make the single-agent normalized dataframes subsetting ref_df
      #These dataframes will have the normalized replicate-level and averaged %Cell Death values
      drug1_df_norm <- 
        ref_df %>%
        filter(drug2_conc==0) %>%
        select(drug1_name, drug1_conc, replicate, replicate_perc_cell_death, perc_cell_death) %>%
        rename("drug_name"=drug1_name,
               "drug_conc"=drug1_conc)
      
      drug2_df_norm <- 
        ref_df %>%
        filter(drug1_conc==0) %>%
        select(drug2_name, drug2_conc, replicate, replicate_perc_cell_death, perc_cell_death) %>%
        rename("drug_name"=drug2_name,
               "drug_conc"=drug2_conc)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Deposit the (normalized) mean %Cell Death values into an 11x11 matrix
      mean_perc_cell_death_mat <- 
        ref_df %>%
        filter(replicate==1) %>% #We're only using (mean) perc_cell_death so it doesn't matter which replicate
        select(drug1_conc, drug2_conc, perc_cell_death) %>%
        reshape2::dcast(drug2_conc ~ drug1_conc) %>%
        .[order(nrow(.):1),] %>% #reorder matrix top-to-bottom (row-wise) so that 0x0 is bottom left
        {rownames(.) <- .$drug2_conc; .} %>% #set rownames to drug2_conc values
        select(-drug2_conc)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Add normalized data to list, and retain appropriate elements from mapped_data
      mylist[[current_plate_name]] <- list(
        "data_mode"=mapped_data$combo_wise_list[[i]]$data_mode,
        "ref_df"=ref_df,
        "drug1_df_norm"=drug1_df_norm,
        "drug2_df_norm"=drug2_df_norm,
        "mean_perc_cell_death_mat"=mean_perc_cell_death_mat,
        "d1_name"=mapped_data$combo_wise_list[[i]]$d1_name,
        "d2_name"=mapped_data$combo_wise_list[[i]]$d2_name,
        "units"=mapped_data$combo_wise_list[[i]]$units,
        "sample_name"=mapped_data$combo_wise_list[[i]]$sample_name
      )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      

    }#End loop over combo_wise_list
    
    
    return(mylist)
    
    
  } else {
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle sparse mode data-----------------------
    #
    #=================================================================
    #=================================================================
    
    #Initialize empty list all results will be stored in (and returned)
    mylist <- list()
    
    
    for(i in seq_along(mapped_data$combo_wise_list)){
      
      
      #Get the current combination name
      current_combo <- names(mapped_data$combo_wise_list[i])
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #mapped reference dataframe
      ref_df <- mapped_data$combo_wise_list[[current_combo]]$ref_df
      
      #Initialize empty column that normalized data will go into in ref_df
      ref_df$replicate_perc_cell_death <- NA
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Normalize ref_df
      
      #Loop over rows of ref_df
      #NOTE: for sparse mode data, raw values for a given drug combination come from different plates
      #It is critical to normalize raw values using the controls from their respective plates
      for(k in 1:nrow(ref_df)){
        
        
        #Extract the filename for the current row
        current_plate_name <- ref_df$filename[k]
        
        #Retrieve the corresponding control values
        pos_ctrl <- mapped_data$combo_wise_list[[current_combo]]$raw_controls_list_pos[[current_plate_name]]
        neg_ctrl <- mapped_data$combo_wise_list[[current_combo]]$raw_controls_list_neg[[current_plate_name]]
        
        #Apply the normalization function
        ref_df$replicate_perc_cell_death[k] <- norm_func(raw_vals = ref_df$raw_val[k],
                                               pos_ctrl_vals = pos_ctrl,
                                               neg_ctrl_vals = neg_ctrl)
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Get the mean %Cell Death across replicates
      ref_df <- 
        ref_df %>%
        group_by(drug1_name, drug2_name, drug1_conc, drug2_conc)  %>%
        mutate("perc_cell_death"=mean(replicate_perc_cell_death)) %>%
        ungroup()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #(Still within each combo)
      
      #Get d1_name (left side of the underscore of current_combo)
      d1_name <- strsplit(current_combo, "_")[[1]][1]
      
      #Get d2_name (right side of the underscore of current_combo)
      d2_name <- strsplit(current_combo, "_")[[1]][2]
      
      #Get single-agent concentrations
      d1_concs <- 
        ref_df %>%
        filter(plate_type == "single_agent") %>%
        filter(drug1_name == d1_name) %>%
        arrange(-drug1_conc) %>% #Descending concs (high-to-low)
        pull(drug1_conc) %>%
        unique()
      
      d2_concs <- 
        ref_df %>%
        filter(plate_type == "single_agent") %>%
        filter(drug1_name == d2_name) %>%         #Note the column is still drug1_name here
        arrange(-drug1_conc) %>%
        pull(drug1_conc) %>%
        unique()
      
      #Make the single-agent normalized dataframes subsetting ref_df
      #These dataframes will have the normalized replicate-level and averaged %Cell Death values
      drug1_df_norm <- 
        ref_df %>%
        filter(plate_type =="single_agent") %>%
        filter(drug1_name == d1_name) %>%
        select(drug1_name, drug1_conc, replicate, replicate_perc_cell_death, perc_cell_death) %>%
        rename("drug_name"=drug1_name,
               "drug_conc"=drug1_conc)
      
      drug2_df_norm <- 
        ref_df %>%
        filter(plate_type =="single_agent") %>%
        filter(drug1_name == d2_name) %>% #Note the column is still drug1_name here
        select(drug1_name, drug1_conc, replicate, replicate_perc_cell_death, perc_cell_death) %>%
        rename("drug_name"=drug1_name,
               "drug_conc"=drug1_conc)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Deposit the (normalized) mean %Cell Death values into a (sparse) 11x11 matrix:
      mytmpdf <- 
        rbind(
          ref_df %>%
            filter(plate_type == "single_agent") %>%
            filter(replicate  == 1), #Replicate number is arbitrary here since we're using the mean %Cell Death
          
          ref_df %>%
            filter(plate_type=="combination")
        )
      
      
      #Obtain the drug1 and drug2 (normalized, mean) single-agent values: 
      #Note: Values correspond to high-to-low concentrations
      d1_norm_vals <- 
        mytmpdf %>%
        filter(plate_type == "single_agent") %>%
        filter(drug1_name == d1_name) %>%
        arrange(nth_conc) %>%
        pull(perc_cell_death)
      
      d2_norm_vals <- 
        mytmpdf %>%
        filter(plate_type == "single_agent") %>%
        filter(drug1_name == d2_name) %>%
        arrange(nth_conc) %>%
        pull(perc_cell_death)
      
      #Initialize an empty matrix that the normalized data for the current replicate will populate
      mytmpmat <- 
        matrix(NA, nrow=11, ncol=11) %>%
        as.data.frame() 
      
      #Set the col/row names of the sparse matrix as the corresponding drug concentrations
      colnames(mytmpmat) <- c(0, rev(d1_concs))
      rownames(mytmpmat) <- c(d2_concs, 0)
      
      #Add the single-agent (normalized) values to the matrix
      mytmpmat[11, 1:11] <- c(0, rev(d1_norm_vals)) #Drug1 single-agent values
      mytmpmat[1:11, 1]  <- c(d2_norm_vals, 0) #Drug2 single-agent values
      
      #Populate the sparse matrix with the combination data
      for(m in 1:nrow(mytmpdf)) {
        # Extract the drug concentrations and raw value for this row
        drug1_conc <- mytmpdf$drug1_conc[m]
        drug2_conc <- mytmpdf$drug2_conc[m]
        norm_value  <- mytmpdf$perc_cell_death[m]
        
        # Determine the appropriate indices in the matrix for these concentrations
        col_index <- which(c(0, rev(d1_concs)) == drug1_conc)
        row_index <- which(c(d2_concs, 0) == drug2_conc)
        
        # Place the raw value in the matrix
        if (length(col_index) > 0 && length(row_index) > 0) {
          mytmpmat[row_index, col_index] <- norm_value
        }
      }
      
      mean_perc_cell_death_mat <- mytmpmat
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Add normalized data to list, and retain appropriate elements from mapped_data
      mylist[[current_combo]] <- list(
        "data_mode"=mapped_data$combo_wise_list[[i]]$data_mode,
        "ref_df"=ref_df,
        "drug1_df_norm"=drug1_df_norm,
        "drug2_df_norm"=drug2_df_norm,
        "mean_perc_cell_death_mat"=mean_perc_cell_death_mat,
        "d1_name"=mapped_data$combo_wise_list[[i]]$d1_name,
        "d2_name"=mapped_data$combo_wise_list[[i]]$d2_name,
        "units"=mapped_data$combo_wise_list[[i]]$units,
        "sample_name"=mapped_data$combo_wise_list[[i]]$sample_name
      )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Print the progress (iterations over the combo_wise_list loop)
      percent_complete <- (i/length(mapped_data$combo_wise_list) * 100)
      cat(sprintf("\rProgress: %.2f%% complete", percent_complete))
      flush.console()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
 
    }#End loop over combo_wise_list
    
    
    return(mylist)
    
    
  } #End handle sparse mode
  
  
}#End function
#--------------------------------------------------------------