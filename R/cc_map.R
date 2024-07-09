#--------------------------------------------------------------
#Function to map plates using metadata
cc_map <- function(metadf,
                   files_dir,
                   control_outlier_threshold = 2.5,
                   save_raw_plate_heatmaps = FALSE,
                   save_transfer_fail_heatmaps = FALSE){
  
  
  
  #Input requirements:
  #metadf-----------------------A simplified metadata file (made by user) or complete metadata file output from cc_makeMeta. (.CSV)
  #files_dir--------------------A directory that contains all raw data files listed in metadf. (.CSV)
  #control_outlier_threshold----Z-score based threshold used for removing outliers.
  #.............................Operates on raw cell death and vehicle control values
  #.............................Higher number is more relaxed (aka fewer values will be dropped)
  #.............................Set this value to 1000 to disable outlier filtering (as nothing will get dropped)
  #.............................Note that Z-prime values are calculated using outlier-filtered raw control values
  #save_raw_plate_heatmaps------Option to save raw plate heatmaps
  #save_transfer_fail_heatmaps--Option to save transfer fail heatmaps (applies only for complete-type metadata)
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Determine whether the metadata type is 'simplified' or 'complete'
  if("failed_transfer" %in% colnames(metadf)){ #Only 'complete'-type metadata has this column
    metadata_type <- "complete"
    print("Complete-type metadata detected")
  } else {
    metadata_type <- "simplified"
    print("simpified-type metadata detected")}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Determine whether data is dense or sparse mode
  if("position_id" %in% colnames(metadf)){ #Only sparse mode has this column
    data_mode <- "sparse"
    print("Detected data type is SPARSE MODE. Proceeding with sparse mode mapping...")
  } else {
    data_mode <- "dense"
    print("Detected data type is DENSE MODE. Proceeding with dense mode mapping...")}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  if(data_mode=="dense"){
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle dense mode data------------------------
    #
    #=================================================================
    #=================================================================
    
    
    setwd(files_dir)

    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize combination-wise list
    #Each element returned will relate to a specific combination (ex: raw mapped data)
    
    combo_wise_list <- list() #Initialize empty list that is returned at the end of the function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize plate-wise lists (or dataframes)
    #Each element will contain metrics for all plates (ex: list of all z-prime values, or plate heatmaps)
    
    raw_plate_heatmap_list <- list()
    
    transfer_fails_heatmap_list <- list()
    
    zprime_df <- data.frame("platename"=character(),
                            "zprime"  =numeric(),
                            "robust_zprime"=numeric())
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Read the raw files (as named in metadata) into a list
    rawfiles_list <- list()
    for(i in unique(metadf$filename)){
      rawfiles_list[[i]] <- read_csv(i, col_names = FALSE)
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(i in seq_along(rawfiles_list)){
      
      mydf <- as.data.frame(rawfiles_list[[i]][1:16,c(1:24)]) #16rows x 24cols (384-wells)
      colnames(mydf) <- as.character(seq(1:24))
      rownames(mydf) <- LETTERS[1:16]
      
      #Extract replicates (without the last 2 transposed rows)
      rep1 <- mydf[1:10,c(1,4,7,10,13,16,19,22)]
      rep2 <- mydf[1:10,c(2,5,8,11,14,17,20,23)]
      rep3 <- mydf[1:10,c(3,6,9,12,15,18,21,24)]
      
      #Extract transposed sections (and back-transpose them so they're now columns)
      rep1_transposed <- t(mydf[c(11,14),c(1:10)])
      rep2_transposed <- t(mydf[c(12,15),c(1:10)])
      rep3_transposed <- t(mydf[c(13,16),c(1:10)])
      
      #Bind the last 2 columns to each replicate
      rep1 <- cbind(rep1, rep1_transposed)
      rep2 <- cbind(rep2, rep2_transposed)
      rep3 <- cbind(rep3, rep3_transposed)
      
      #Extract the single-agent values (aka "drug-only" values)
      #Transpose them so they are in column format
      rep1_drug1_only <- t(mydf[c(11),c(11:20)])
      rep1_drug2_only <- t(mydf[c(14),c(11:20)])
      
      rep2_drug1_only <- t(mydf[c(12),c(11:20)])
      rep2_drug2_only <- t(mydf[c(15),c(11:20)])
      
      rep3_drug1_only <- t(mydf[c(13),c(11:20)])
      rep3_drug2_only <- t(mydf[c(16),c(11:20)])
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Split off the raw control values 
      #Each file (each combination, in dense mode) will have 12 values each for pos and neg controls
      
      #Initialize empty lists for the raw pos/neg controls to go into
      raw_controls_list_pos <- list()
      raw_controls_list_neg <- list()
      
      #Get the current plate name (file name)
      platename <- names(rawfiles_list[i])
      
      #Negative control values (DMSO / vehicle control)
      neg_ctrl <- 
        mydf[c(11:16),c(21:22)] %>%
        stack() %>%
        pull(values) %>%
        drop_outliers(., threshold = control_outlier_threshold) #Drop any outliers
      
      #Positive control values (Staurosporine / cell death control)
      pos_ctrl <- 
        mydf[c(11:16),c(23:24)] %>%
        stack() %>%
        pull(values) %>%
        drop_outliers(., threshold = control_outlier_threshold) #Drop any outliers
      
      
      #Deposit results into list
      raw_controls_list_pos[[platename]] <- pos_ctrl
      raw_controls_list_neg[[platename]] <- neg_ctrl
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate the z-prime and robust-zprime, and deposit into the zprime_df
      zprime_value  <- zprime_func(pos_ctrl = pos_ctrl,
                                   neg_ctrl = neg_ctrl)
      
      robust_zprime_value <- robust_zprime_func(pos_ctrl = pos_ctrl,
                                                neg_ctrl = neg_ctrl)
      
      #Deposit zprime and robust_zprime as new row in the zprime_df
      zprime_df <- 
        rbind(zprime_df,
              data.frame("platename"=platename,
                         "zprime"=zprime_value,
                         "robust_zprime"=robust_zprime_value)
        )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Generate (and optionally, save) heatmaps of the raw data for each plate
      
      #Make a title that will become the plot title
      tmp_title <- 
        paste0(platename, " ",
               "(Z' = ", round(zprime_df %>%
                                 filter(platename == names(rawfiles_list[i])) %>%
                                 pull(zprime),
                               digits = 2),
               ")")
      
      raw_plate_heatmap <- 
        mydf %>%
        as.matrix() %>%
        reshape2::melt() %>%
        as.data.frame() %>%
        rename("Row_Names" = "Var1",
               "Column_Names" = "Var2") %>%
        mutate("Column_Names" = factor(Column_Names, levels = 1:24), 
               "Row_Names" = factor(Row_Names, levels = rev(LETTERS[1:16]))) %>%
        ggplot(aes(x = Column_Names, y = Row_Names, fill = value)) +
        geom_tile(color="gray40") +
        xlab("") +
        ylab("") +
        scale_fill_gradient("", low=c("darkblue", "white") ,high=c("pink", "red3")) +
        ggtitle(tmp_title) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_discrete(expand = c(0, 0)) + # Remove x axis expansion
        scale_y_discrete(expand = c(0, 0))   # Remove y axis expansion
      
      #Deposit plot into list
      raw_plate_heatmap_list[[platename]] <- raw_plate_heatmap
      
      #Save raw plate heatmaps if option is selected
      if(save_raw_plate_heatmaps==TRUE){
        
        #Create the raw_heatmaps/ directory if it doesn't already exist
        if(!dir.exists(file.path(files_dir, "raw_heatmaps"))){
          dir.create(file.path(files_dir, "raw_heatmaps"))
        }
        
        print(paste0("Saving heatmap to file for plate ", platename)) #Print message
        ggsave(paste0("raw_heatmaps/", "heatmap_", platename, ".svg"),
               raw_plate_heatmap,
               width=9,
               height=5)
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Add the drug concentrations to reps
      #Drug1 decreases horizontally (col-wise). Left-to-right is high-to-low (this is done for ease, but gets flipped downstream)
      #Drug2 decreases vertically (row-wise). Top-to-bottom is high-to-low
      
      d1_concs <- 
        metadf %>%
        filter(filename == names(rawfiles_list[i])) %>%
        filter(!drug1_conc==0)    %>%   #Filter out 0 here (which is important for 'complete'-type metadata)
        arrange(desc(drug1_conc)) %>%
        pull(drug1_conc) %>%
        unique()
      
      d2_concs <- 
        metadf %>%
        filter(filename == names(rawfiles_list[i])) %>%
        filter(!drug2_conc==0)    %>%   #Filter out 0 here (which is important for 'complete'-type metadata)
        arrange(desc(drug2_conc)) %>%
        pull(drug2_conc) %>%
        unique()
      
      
      colnames(rep1) <- d1_concs
      rownames(rep1) <- d2_concs
      
      colnames(rep2) <- d1_concs
      rownames(rep2) <- d2_concs
      
      colnames(rep3) <- d1_concs
      rownames(rep3) <- d2_concs
      
      #Add the drug concentrations to the single-agent values
      #Note the single-agents are both columns so they each get new 'rownames'
      rownames(rep1_drug1_only) <- d1_concs
      rownames(rep1_drug2_only) <- d2_concs
      
      rownames(rep2_drug1_only) <- d1_concs
      rownames(rep2_drug2_only) <- d2_concs
      
      rownames(rep3_drug1_only) <- d1_concs
      rownames(rep3_drug2_only) <- d2_concs
      

      #Add the single-agent values (forms an 11x11 matrix)
      rep1 <- rbind(rep1, t(rep1_drug1_only))
      rep1 <- cbind(rep1, rbind(rep1_drug2_only,0))
      
      rep2 <- rbind(rep2, t(rep2_drug1_only))
      rep2 <- cbind(rep2, rbind(rep2_drug2_only,0))
      
      rep3 <- rbind(rep3, t(rep3_drug1_only))
      rep3 <- cbind(rep3, rbind(rep3_drug2_only,0))
      
      #Rename the last column and last row to be '0'
      names(rep1)[length(names(rep1))] <- c("0")
      rownames(rep1)[length(rownames(rep1))] <- c("0")
      
      names(rep2)[length(names(rep2))] <- c("0")
      rownames(rep2)[length(rownames(rep2))] <- c("0")
      
      names(rep3)[length(names(rep3))] <- c("0")
      rownames(rep3)[length(rownames(rep3))] <- c("0")
      
      #Flip the dataframe left to right 
      #Drug1 should be columns where concentration increases left-to-right
      #Drug2 should be rows where concentration increases bottom-to-top
      #No-drug (0x0) should be in the lower left
      rep1 <- rev(rep1)
      rep2 <- rev(rep2)
      rep3 <- rev(rep3)
      
      
      #Make the drug-only data dataframes with columns "drug_conc" and "raw_value"
      rep1_drug1_only <- rep1_drug1_only %>% as.data.frame() %>% rownames_to_column(var = "drug_conc") %>% 'colnames<-'(c("drug_conc", "raw_value"))
      rep1_drug2_only <- rep1_drug2_only %>% as.data.frame() %>% rownames_to_column(var = "drug_conc") %>% 'colnames<-'(c("drug_conc", "raw_value"))
      
      rep2_drug1_only <- rep2_drug1_only %>% as.data.frame() %>% rownames_to_column(var = "drug_conc") %>% 'colnames<-'(c("drug_conc", "raw_value"))
      rep2_drug2_only <- rep2_drug2_only %>% as.data.frame() %>% rownames_to_column(var = "drug_conc") %>% 'colnames<-'(c("drug_conc", "raw_value"))
      
      rep3_drug1_only <- rep3_drug1_only %>% as.data.frame() %>% rownames_to_column(var = "drug_conc") %>% 'colnames<-'(c("drug_conc", "raw_value"))
      rep3_drug2_only <- rep3_drug2_only %>% as.data.frame() %>% rownames_to_column(var = "drug_conc") %>% 'colnames<-'(c("drug_conc", "raw_value"))
      
      
      #Grab the units of the drug doses and the drug names
      units <- 
        metadf %>%
        filter(filename == names(rawfiles_list[i])) %>%
        pull(units) %>%
        unique()
      
      sample_name <- 
        metadf %>%
        filter(filename == names(rawfiles_list[i])) %>%
        pull(sample_name) %>%
        unique()
        
      d1_name <- 
        metadf %>%
        filter(filename == names(rawfiles_list[i])) %>%
        filter(!is.na(drug1_name)) %>% #Important for 'complete'-type metadata (where valid NAs in this column exist, like when it's drug2-only)
        pull(drug1_name) %>%
        unique()
      
      d2_name <- 
        metadf %>%
        filter(filename == names(rawfiles_list[i])) %>%
        filter(!is.na(drug2_name)) %>% #Important for 'complete'-type metadata (where valid NAs in this column exist, like when it's drug1-only)
        pull(drug2_name) %>%
        unique()
      
      
      
      #Make 'combo_id' as a combination of the drug names in a single string
      combo_id <- paste0(d1_name, "_", d2_name)
      
      
      #Deposit the replicates into respective lists for matrices and single-agents
      rep_mats_list_raw   <- list(rep1, rep2, rep3)
      
      rep_drug1_list_raw <- list(rep1_drug1_only,
                                 rep2_drug1_only,
                                 rep3_drug1_only)
      rep_drug2_list_raw <- list(rep1_drug2_only,
                                 rep2_drug2_only,
                                 rep3_drug2_only)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Construct a reference dataframe; an extension of the metadata that has mapped values in a dataframe
      #This is handled differently depending on metadata type
      
      #NOTE: There should be 363 rows in the resulting ref_df
      
      if(metadata_type=="simplified"){
        
        ref_df <- 
          rbind(
            
            rep1 %>%
              rownames_to_column(var = "drug2_conc") %>%
              pivot_longer(cols = -drug2_conc, names_to = "drug1_conc", values_to = "raw_val") %>%
              mutate("replicate" = 1), 
            
            rep2 %>%
              rownames_to_column(var = "drug2_conc") %>%
              pivot_longer(cols = -drug2_conc, names_to = "drug1_conc", values_to = "raw_val") %>%
              mutate("replicate" = 2), 
            
            rep3 %>%
              rownames_to_column(var = "drug2_conc") %>%
              pivot_longer(cols = -drug2_conc, names_to = "drug1_conc", values_to = "raw_val") %>%
              mutate("replicate" = 3)
          ) %>%
          mutate("drug1_name"=d1_name,
                 "drug2_name"=d2_name,
                 "units"=units,
                 "sample_name"=sample_name,
                 "filename"=platename,
                 "combo_id"=combo_id) %>%
          
          #Ensure the drug concs are numerics
          mutate("drug1_conc"=as.numeric(drug1_conc), 
                 "drug2_conc"=as.numeric(drug2_conc)) %>%
          
          #Re-arrange column order
          select(drug1_name, 
                 drug1_conc, 
                 drug2_name,
                 drug2_conc, 
                 combo_id,
                 sample_name,
                 units,
                 filename,
                 raw_val, 
                 replicate)
        
      } else {
        
        
        #Testing just how ref_df would look with complete-type metadata
        ref_df <- 
          rbind(
            
            rep1 %>%
              rownames_to_column(var = "drug2_conc") %>%
              pivot_longer(cols = -drug2_conc, names_to = "drug1_conc", values_to = "raw_val") %>%
              mutate("replicate" = 1), 
            
            rep2 %>%
              rownames_to_column(var = "drug2_conc") %>%
              pivot_longer(cols = -drug2_conc, names_to = "drug1_conc", values_to = "raw_val") %>%
              mutate("replicate" = 2), 
            
            rep3 %>%
              rownames_to_column(var = "drug2_conc") %>%
              pivot_longer(cols = -drug2_conc, names_to = "drug1_conc", values_to = "raw_val") %>%
              mutate("replicate" = 3)) %>%
          
          #Ensure the drug concs are numerics
          mutate(drug1_conc = as.numeric(drug1_conc), 
                 drug2_conc = as.numeric(drug2_conc)) %>%
          
          #Add the drug1/2_names
          mutate("drug1_name"=d1_name,
                 "drug2_name"=d2_name) %>%
          
          full_join(., 
                    metadf %>%
                      filter(filename == names(rawfiles_list[i])),
                      by=c("drug1_name",
                           "drug1_conc",
                           "drug2_name",
                           "drug2_conc",
                           "replicate"))
          
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Merge the zprime_df results back to ref_df
      ref_df <- 
        ref_df %>%
        merge(zprime_df, 
              by.x="filename",
              by.y="platename",
              all.x=TRUE)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Generate (and optionally, save) heatmaps of any plates with failed transfers (complete metadata only)
      
      if(metadata_type=="complete"){
        
        #Notes: 
        #1....There can be 2 failed transfers per well (one for each drug), but that won't be seen in the heatmap
        #2....The control wells are not tracked for transfer failures, so they are blank on the heatmap
        
        #Get the number of failed transfers within the plate
        num_failed_transfers <- sum(ref_df$failed_transfer, na.rm = TRUE)
        
        #ONLY Make the heatmaps for plates with transfer fails
        if(num_failed_transfers > 0){
          
          #Make a title that will become the plot title
          tmp_title_transfer_fails <- 
            paste0(platename, " ",
                   "(", num_failed_transfers, " failed transfers", ")")
          
          
          transfer_fails_heatmap <- 
            ref_df %>%
            filter(!(drug1_conc == 0 & drug2_conc == 0)) %>%
            ggplot(aes(x = factor(plate_col), 
                       y = factor(plate_row_letter, levels=rev(LETTERS[1:16])), 
                       fill = as.factor(failed_transfer))) +
            geom_tile(color="gray40") +
            xlab("") +
            ylab("") +
            scale_fill_manual("failed", values=c("1"="salmon", "0"="gray")) +
            ggtitle(tmp_title_transfer_fails) +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_x_discrete(expand = c(0, 0)) + # Remove x axis expansion
            scale_y_discrete(expand = c(0, 0))   # Remove y axis expansion
          
          transfer_fails_heatmap_list[[platename]] <- transfer_fails_heatmap
          
          #Save the transfer fail heatmaps if option is selected
          if(save_transfer_fail_heatmaps==TRUE){
            
            #Create transfer_fail_heatmaps/ directory if it doesn't already exist
            if(!dir.exists(file.path(files_dir, "transfer_fail_heatmaps"))){
              dir.create(file.path(files_dir, "transfer_fail_heatmaps"))
            }
            
            print(paste0("Saving transfer fail heatmap for plate ", platename)) #Print message
            ggsave(paste0("transfer_fail_heatmaps/", "heatmap_transfer_fail_", platename, ".svg"),
                   transfer_fails_heatmap,
                   width=9,
                   height=5)
          }
          
        }
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Deposit combination-wise data into combo_wise_list
      
      combo_wise_list[[i]] <- list(
        "data_mode"=data_mode,
        "ref_df"=ref_df,
        "rep_mats_list_raw"=rep_mats_list_raw,
        "rep_drug1_list_raw"=rep_drug1_list_raw,
        "rep_drug2_list_raw"=rep_drug2_list_raw,
        "raw_controls_list_pos"=raw_controls_list_pos,
        "raw_controls_list_neg"=raw_controls_list_neg,
        "d1_name"=d1_name,
        "d2_name"=d2_name,
        "units"=units,
        "sample_name"=sample_name
      )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Deposit plate-wise data into plate_metrics_list
      
      plate_metrics_list <- list(
        "raw_plate_heatmap_list"=raw_plate_heatmap_list,
        "transfer_fails_heatmap_list"=transfer_fails_heatmap_list,
        "zprime_df"=zprime_df
      )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Place the combo_wise_list and plate_metrics_list together in a list that gets returned
      
      final_list <- list(
        "combo_wise_list"=combo_wise_list,
        "plate_metrics_list"=plate_metrics_list
      )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
    }#End loop over files
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Set names of resulting list to names of input list
    names(final_list$combo_wise_list) <- names(rawfiles_list)
    
    
    #Return the final list
    return(final_list)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
  } else {
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle sparse mode data-----------------------
    #
    #=================================================================
    #=================================================================
    
    
    
    #Set working directory to files_dir (which contains raw data)
    setwd(files_dir)
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize combination-wise list
    #Each element returned will relate to a specific combination (ex: raw mapped data)
    
    combo_wise_list <- list() #Initialize empty list that is returned at the end of the function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize plate-wise lists (or dataframes)
    #Each element will contain metrics for all plates (ex: list of all z-prime values, or plate heatmaps)
    
    raw_plate_heatmap_list <- list()
    
    zprime_df <- data.frame("platename"=character(),
                            "zprime"  =numeric(),
                            "robust_zprime"=numeric())
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #Read the raw files (as named in metadata) into a list
    rawfiles_list <- list()
    for(i in unique(metadf$filename)){
      rawfiles_list[[i]] <- read_csv(i, col_names = FALSE)
    }
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Extract the raw control values for each plate/file, and use them to get z-prime values
    #Each file will have 14 pos and 14 neg controls
    #Each file will also have a z-prime value and a robust z-prime value
    
    #Initialize empty lists for raw results to go into
    raw_controls_list_pos <- list()
    raw_controls_list_neg <- list()
    
    
    
    for(i in seq_along(rawfiles_list)){
      
      #Get plate name
      platename <- names(rawfiles_list)[i]
      
      #Positive control values (Staurosporine / cell death control)
      pos_ctrl <- 
        rawfiles_list[[i]][2:15,c(47)]  %>% 
        pull() %>%
        drop_outliers(., threshold = control_outlier_threshold) #Drop any outliers
      
      #Negative control values (DMSO / vehicle control)
      neg_ctrl <- 
        rawfiles_list[[i]][18:31,c(47)]  %>%
        pull() %>%
        drop_outliers(., threshold = control_outlier_threshold) #Drop any outliers
      
      #Deposit results into list
      raw_controls_list_pos[[platename]] <- pos_ctrl
      raw_controls_list_neg[[platename]] <- neg_ctrl
      
      
      #Calculate the z-prime and robust-zprime, and deposit into the zprime_df
      zprime_value  <- zprime_func(pos_ctrl = pos_ctrl,
                                   neg_ctrl = neg_ctrl)
      
      robust_zprime_value <- robust_zprime_func(pos_ctrl = pos_ctrl,
                                                neg_ctrl = neg_ctrl)
      
      #Deposit zprime and robust_zprime as new row in the zprime_df
      zprime_df <- 
        rbind(zprime_df,
              data.frame("platename"=platename,
                         "zprime"=zprime_value,
                         "robust_zprime"=robust_zprime_value)
        )
      
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Generate (and optionally, save) heatmaps of the raw data for each plate
    for(i in seq_along(rawfiles_list)){
      
      #Get plate name
      platename <- names(rawfiles_list)[i]
      
      #Figure out if the current plate is a single-agent or combination plate
      #If it's a single-agent plate, denote that on the raw plate heatmap plot with asterisk (*)
      if(
        #This will be TRUE in the case of single-agent plates (where the drug2_name column is NA)
        metadf %>%
         filter(filename==platename) %>%
         count(is.na(drug2_name)) %>%
         .[1] == TRUE){
        platename_modified <- paste0(platename,"*")
      } else {
        platename_modified <- platename #Otherwise, it's a combination plate so keep original name
      }
      
      
        
      
      #Make a title that will become the plot title
      tmp_title <- 
        paste0(platename_modified, " ",
               "(Z' = ", round(zprime_df %>%
                                 filter(platename == names(rawfiles_list[i])) %>%
                                 pull(zprime),
                               digits = 2), 
               ")")
      
      raw_plate_heatmap <- 
        rawfiles_list[[i]] %>%
        as.matrix() %>%
        reshape2::melt() %>%
        as.data.frame() %>%
        rename("Row_Names" = "Var1",
               "Column_Names" = "Var2") %>%
        mutate("Column_Names" = gsub("X", "", Column_Names)) %>%
        mutate("Column_Names" = factor(Column_Names, levels = 1:48),
               "Row_Names" = factor(Row_Names, levels = rev(1:32))) %>%
        ggplot(aes(x = Column_Names, y = Row_Names, fill = value)) +
        geom_tile(color="gray40") +
        xlab("") +
        ylab("") +
        scale_fill_gradient("", low=c("darkblue", "white") ,high=c("pink", "red3")) +
        ggtitle(tmp_title) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_discrete(expand = c(0, 0)) + # Remove x axis expansion
        scale_y_discrete(expand = c(0, 0))   # Remove y axis expansion
      
      #Deposit plot into list
      raw_plate_heatmap_list[[platename]] <- raw_plate_heatmap
      
      #Save raw plate heatmaps if option is selected
      if(save_raw_plate_heatmaps==TRUE){
        
        #Create the raw_heatmaps/ directory if it doesn't already exist
        if(!dir.exists(file.path(files_dir, "raw_heatmaps"))){
          dir.create(file.path(files_dir, "raw_heatmaps"))
        }
        
        print(paste0("Saving heatmap to file for plate ", platename)) #Print message
        ggsave(paste0("raw_heatmaps/", "heatmap_", platename, ".svg"),
               raw_plate_heatmap,
               width=9,
               height=5)
      }
      
    } #End loop to generate raw plate heatmaps
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    #Make a combo_id column (which is not present in simplified metadata)
    metadf <- 
      metadf %>%
      mutate(combo_id = paste0(drug1_name, "_", drug2_name))
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Deposit the raw values into appropriate cells of the sparse matrices
    
    #Get all combinations from data:
    mycombos <- 
      metadf %>%
      filter(plate_type=="combination") %>%
      pull(combo_id) %>%
      unique()
    
    
    #Begin by isolating the current combination
    for(i in seq_along(mycombos)){
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Initialize empty list for the sparse matrices to go into
      rep_mats_list_raw  <- list() 
      
      #Initialize empty lists for single-agent raw values 
      rep_drug1_list_raw <- list()
      rep_drug2_list_raw <- list()
      
      #Obtain the sample name and units for current combo_id
      sample_name <- metadf %>% filter(combo_id == mycombos[i]) %>% pull(sample_name) %>% unique()
      units_name  <- metadf %>% filter(combo_id == mycombos[i]) %>% pull(units) %>% unique()
      
      #Obtain the drug1 and drug2 names from combo_id (format is drug1_drug2)
      d1_name <- strsplit(mycombos[i], "_")[[1]][1]
      d2_name <- strsplit(mycombos[i], "_")[[1]][2]
      
      #Obtain the drug1 and drug2 concentrations: (high-to-low)
      d1_concs <- 
        metadf %>%
        filter(plate_type == "single_agent") %>%
        filter(drug1_name == d1_name) %>%
        arrange(-drug1_conc) %>% #Descending concs (high-to-low)
        pull(drug1_conc) %>%
        unique()
      
      d2_concs <- 
        metadf %>%
        filter(plate_type == "single_agent") %>%
        filter(drug1_name == d2_name) %>%         #Note the column is still drug1_name here
        arrange(-drug1_conc) %>%
        pull(drug1_conc) %>%
        unique()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Make a reference data frame for mapping
      #ref_df will be an extension of the metadata that contains mapped values
      
      
      
      #Create 3 row groups (1,4,5,10...133), (2,5,8,11...134), and (3,6,9,12...135). To help with figuring out which row each position_id belongs to
      rowgroup1 <- seq(from=1, to=135, by=3)
      rowgroup2 <- seq(from=2, to=135, by=3)
      rowgroup3 <- seq(from=3, to=135, by=3)
      
      #NOTE: We explicitly (manually) put in the position_id here, for extra clarity
      #However, it is possible to calculate the position_id using the row and column numbers
      #This can be done with a modulo operation (see 'calculate_position_id' helper function)
      
      
      ref_df <- 
        metadf %>%
        filter(combo_id==mycombos[i] | #Filter to current combo_id and single-agents
                 combo_id==paste0(d1_name,"_NA") | 
                 combo_id==paste0(d2_name,"_NA")) %>%
        group_by(filename, position_id) %>% #Group data filename and position_id
        arrange(-drug1_conc) %>% #sort Descending
        #Add a column called 'nth_conc' which displays which concentration (1:10, high-to-low) a given row corresponds to
        mutate(nth_conc=seq(1:10)) %>%
        ungroup() %>%
        #Figure out the corresponding plate row based on rowgroup and nth_conc
        mutate(plate_row = case_when(
          position_id %in% rowgroup1 & nth_conc == 1 ~ 2,
          position_id %in% rowgroup2 & nth_conc == 1 ~ 3,
          position_id %in% rowgroup3 & nth_conc == 1 ~ 4,
          position_id %in% rowgroup1 & nth_conc == 2 ~ 5,
          position_id %in% rowgroup2 & nth_conc == 2 ~ 6,
          position_id %in% rowgroup3 & nth_conc == 2 ~ 7,
          position_id %in% rowgroup1 & nth_conc == 3 ~ 8,
          position_id %in% rowgroup2 & nth_conc == 3 ~ 9,
          position_id %in% rowgroup3 & nth_conc == 3 ~ 10,
          position_id %in% rowgroup1 & nth_conc == 4 ~ 11,
          position_id %in% rowgroup2 & nth_conc == 4 ~ 12,
          position_id %in% rowgroup3 & nth_conc == 4 ~ 13,
          position_id %in% rowgroup1 & nth_conc == 5 ~ 14,
          position_id %in% rowgroup2 & nth_conc == 5 ~ 15,
          position_id %in% rowgroup3 & nth_conc == 5 ~ 16,
          position_id %in% rowgroup1 & nth_conc == 6 ~ 17,
          position_id %in% rowgroup2 & nth_conc == 6 ~ 18,
          position_id %in% rowgroup3 & nth_conc == 6 ~ 19,
          position_id %in% rowgroup1 & nth_conc == 7 ~ 20,
          position_id %in% rowgroup2 & nth_conc == 7 ~ 21,
          position_id %in% rowgroup3 & nth_conc == 7 ~ 22,
          position_id %in% rowgroup1 & nth_conc == 8 ~ 23,
          position_id %in% rowgroup2 & nth_conc == 8 ~ 24,
          position_id %in% rowgroup3 & nth_conc == 8 ~ 25,
          position_id %in% rowgroup1 & nth_conc == 9 ~ 26,
          position_id %in% rowgroup2 & nth_conc == 9 ~ 27,
          position_id %in% rowgroup3 & nth_conc == 9 ~ 28,
          position_id %in% rowgroup1 & nth_conc == 10 ~ 29,
          position_id %in% rowgroup2 & nth_conc == 10 ~ 30,
          position_id %in% rowgroup3 & nth_conc == 10 ~ 31,
          TRUE ~ NA)) %>%
        #Figure out column that each position_id corresponds to
        #If position_id=1:3, the column should be 2. If position_id=4:6, the column should be 3. etc...
        mutate(plate_col = case_when(
          position_id %in% 1:3 ~ 2,
          position_id %in% 4:6 ~ 3,
          position_id %in% 7:9 ~ 4,
          position_id %in% 10:12 ~ 5,
          position_id %in% 13:15 ~ 6,
          position_id %in% 16:18 ~ 7,
          position_id %in% 19:21 ~ 8,
          position_id %in% 22:24 ~ 9,
          position_id %in% 25:27 ~ 10,
          position_id %in% 28:30 ~ 11,
          position_id %in% 31:33 ~ 12,
          position_id %in% 34:36 ~ 13,
          position_id %in% 37:39 ~ 14,
          position_id %in% 40:42 ~ 15,
          position_id %in% 43:45 ~ 16,
          position_id %in% 46:48 ~ 17,
          position_id %in% 49:51 ~ 18,
          position_id %in% 52:54 ~ 19,
          position_id %in% 55:57 ~ 20,
          position_id %in% 58:60 ~ 21,
          position_id %in% 61:63 ~ 22,
          position_id %in% 64:66 ~ 23,
          position_id %in% 67:69 ~ 24,
          position_id %in% 70:72 ~ 25,
          position_id %in% 73:75 ~ 26,
          position_id %in% 76:78 ~ 27,
          position_id %in% 79:81 ~ 28,
          position_id %in% 82:84 ~ 29,
          position_id %in% 85:87 ~ 30,
          position_id %in% 88:90 ~ 31,
          position_id %in% 91:93 ~ 32,
          position_id %in% 94:96 ~ 33,
          position_id %in% 97:99 ~ 34,
          position_id %in% 100:102 ~ 35,
          position_id %in% 103:105 ~ 36,
          position_id %in% 106:108 ~ 37,
          position_id %in% 109:111 ~ 38,
          position_id %in% 112:114 ~ 39,
          position_id %in% 115:117 ~ 40,
          position_id %in% 118:120 ~ 41,
          position_id %in% 121:123 ~ 42,
          position_id %in% 124:126 ~ 43,
          position_id %in% 127:129 ~ 44,
          position_id %in% 130:132 ~ 45,
          position_id %in% 133:135 ~ 46
        ))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Populate raw_val column with data from a corresponding raw data file
      for(k in 1:nrow(ref_df)){
        
        filename  <- as.character(ref_df[k, "filename"])
        plate_row <- as.numeric(ref_df[k, "plate_row"])
        plate_col <- as.numeric(ref_df[k, "plate_col"])
        
        raw_data <- rawfiles_list[[which(names(rawfiles_list) == filename)]] 
        
        #Get the raw value
        raw_val <- raw_data[plate_row, plate_col] %>% as.numeric()
        
        ref_df[k, "raw_val"] <- raw_val
        
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Isolate current replicate within current combination
      for(k in unique(ref_df$replicate)){
        
        mytmpdf <- 
          rbind(
            
            #Single-agent data for drug1, replicate 'k'
            ref_df %>%
              filter(plate_type == "single_agent") %>%
              filter(replicate  == k),
            
            #Combination data for current combo
            #NOTE: Combinations don't have replicates:
            ref_df %>%
              filter(plate_type == "combination")
            
          )
        
        #Obtain the drug1 and drug2 (raw) single-agent values: 
        #Note: Values correspond to high-to-low concentrations, but the values will likely be low-to-high since high conc produces low viability values
        d1_raw_vals <- 
          mytmpdf %>%
          filter(plate_type == "single_agent") %>%
          filter(drug1_name == d1_name) %>%
          arrange(nth_conc) %>%
          pull(raw_val)
        
        d2_raw_vals <- 
          mytmpdf %>%
          filter(plate_type == "single_agent") %>%
          filter(drug1_name == d2_name) %>%
          arrange(nth_conc) %>%
          pull(raw_val)
        
        #Get the drug1- and drug2-only (raw) data frames for the current replicate
        #Deposit them into the appropriate drug-only lists
        rep_drug1_list_raw[[k]] <- 
          data.frame("drug_conc"=d1_concs,
                     "raw_value"=d1_raw_vals)
        
        rep_drug2_list_raw[[k]] <-
          data.frame("drug_conc"=d2_concs,
                     "raw_value"=d2_raw_vals)
        
        
        #Create a sparse matrix for the current replicate
        mytmpmat <- 
          matrix(NA, nrow = 11, ncol=11) %>%
          as.data.frame() 
        
        #Set the col/row names of the sparse matrix as the corresponding drug concentrations
        colnames(mytmpmat) <- c(0, rev(d1_concs))
        rownames(mytmpmat) <- c(d2_concs, 0)
        
        #Add the single-agent (raw) values to the matrix
        mytmpmat[11, 1:11] <- c(0, rev(d1_raw_vals)) #Drug1 single-agent values
        mytmpmat[1:11, 1]  <- c(d2_raw_vals, 0) #Drug2 single-agent values
        
        
        #Populate the sparse matrix with the combination data
        
        for(m in 1:nrow(mytmpdf)) {
          # Extract the drug concentrations and raw value for this row
          drug1_conc <- mytmpdf$drug1_conc[m]
          drug2_conc <- mytmpdf$drug2_conc[m]
          raw_value  <- mytmpdf$raw_val[m]
          
          # Determine the appropriate indices in the matrix for these concentrations
          col_index <- which(c(0, rev(d1_concs)) == drug1_conc)
          row_index <- which(c(d2_concs, 0) == drug2_conc)
          
          # Place the raw value in the matrix
          if (length(col_index) > 0 && length(row_index) > 0) {
            mytmpmat[row_index, col_index] <- raw_value
          }
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        #Add the sparse matrix to the list
        rep_mats_list_raw[[k]] <- mytmpmat
        
      } #End loop over replicates
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Merge the zprime_df results back to ref_df
      ref_df <- 
        ref_df %>%
        merge(zprime_df, 
              by.x="filename",
              by.y="platename",
              all.x=TRUE)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      
      combo_wise_list[[i]] <- list("data_mode"=data_mode,
                                   "ref_df"=ref_df,
                                   "rep_mats_list_raw"=rep_mats_list_raw,
                                   "rep_drug1_list_raw"=rep_drug1_list_raw,
                                   "rep_drug2_list_raw"=rep_drug2_list_raw,
                                   "raw_controls_list_pos"=raw_controls_list_pos,
                                   "raw_controls_list_neg"=raw_controls_list_neg,
                                   "d1_name"=d1_name,
                                   "d2_name"=d2_name,
                                   "units"=units_name,
                                   "sample_name"=sample_name)
      
      names(combo_wise_list)[[i]] <- mycombos[i]
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Deposit plate-wise data into plate_metrics_list
      
      plate_metrics_list <- list(
        "raw_plate_heatmap_list"=raw_plate_heatmap_list,
        "zprime_df"=zprime_df
      )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Place the combo_wise_list and plate_metrics_list together in a list that gets returned
      
      final_list <- list(
        "combo_wise_list"=combo_wise_list,
        "plate_metrics_list"=plate_metrics_list
      )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Print the progress (iterations over the mycombos loop)
      percent_complete <- (i/length(mycombos) * 100)
      cat(sprintf("\rProgress: %.2f%% complete", percent_complete))
      flush.console()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
    } #End loop over mycombos
    
    
    
    return(final_list)
    
    
  }#End handling of sparse mode data
  
  
}#End of function
#--------------------------------------------------------------