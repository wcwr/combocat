#--------------------------------------------------------------
#' Calculate QC Metrics and Add Summary Statistics for Normalized Data
#'
#' This function calculates Quality Control (QC) metrics, including standard deviation (SD), monotonicity, and residuals, for both dense and sparse mode data. It flags concentrations that exceed user-defined thresholds and adds summary statistics to the reference data frame.
#'
#' NOTE: It is important to run this function AFTER `cc_getSyn`
#'
#'
#' @param norm_data List. Output from `cc_norm`
#' @param dr_data List. Output from `cc_getDR`
#' @param cutoff_sd Numeric. The standard deviation cutoff for flagging concentrations. Default is 29. Values exceeding this are flagged.
#' @param cutoff_mono Numeric. Monotonicity threshold. Should be a negative value because it will be used in the low-to-high evaluation, default is -16.
#'   Should be a negative value because it will be used in the low-to-high evaluation...
#'   Monotonicity is evaluated bi-directionally. When low-to-high, values less than mono_cutoff are flagged
#'   when evaluated high-to-low, values greater than the absolute value of this number will be flagged.
#'   Ex if cutoff_mono = -10:  lowest conc (nth_conc==10) is 50% cell death, then 2nd lowest conc (nth_conc==9) is 5% cell death, nth_conc==9 gets flagged (evaluated low-to-high)
#'   NOTE: when evaluating this high-to-low, nth_conc==10 will also get flagged. A bi-directional evaluation means flagging always occurs in pairs.
#' @param cutoff_resid Numeric. Residual threshold. Concs are flagged if their corresponding response is greater than (the ABSOLUTE VALUE) of this number. Default is 15.
#'
#' @return A list, similar to `norm_data`, with the following added elements:
#' \describe{
#'   \item{`ref_df`}{Data frame. The reference data frame with added flags and summary statistics.}
#'   \item{`flagged_sd`}{Numeric. Flag for concentrations with SD above the `cutoff_sd`.}
#'   \item{`flagged_mono`}{Numeric. Flag for concentrations violating the monotonicity condition.}
#'   \item{`flagged_resid`}{Numeric. Flag for concentrations with residuals greater than the `cutoff_resid`.}
#'   \item{`mean_syn`}{Numeric. Mean synergy for all combinations, excluding single-agent concentrations.}
#'   \item{`median_syn`}{Numeric. Median synergy for all combinations, excluding single-agent concentrations.}
#'   \item{`max_syn`}{Maximum synergy score across all combination wells, excluding single-agent wells.}
#'   \item{`mean_syn_adj`}{Numeric. Mean synergy for all combinations, excluding flagged concentrations.}
#'   \item{`median_syn_adj`}{Numeric. Median synergy for all combinations, excluding flagged concentrations.}
#'   \item{`max_syn_adj`}{Maximum synergy score after excluding flagged concentrations.}
#'   \item{`mean_perc_cell_death`}{Numeric. Mean percent cell death across all concentrations.}
#'   \item{`num_flagged_diag`}{Numeric. The number of flagged diagonal (measured) values.}
#'   \item{`moran_i`}{Numeric. Moran's I value for spatial autocorrelation of synergy scores.}
#'   \item{`moran_i_STND`}{Numeric. Standardized Moran's I value.}
#' }
#'
#' @export
#--------------------------------------------------------------




#--------------------------------------------------------------
#Function to calculate QC metrics and add summary statistics for normalized data
cc_getQC <- function(norm_data,
                     dr_data,
                     cutoff_sd    =  29,
                     cutoff_mono  = -16,
                     cutoff_resid =  15){
  
  
  #Note: It is important to run this function after the synergy (and other) functions
  #Any discrepancies in concentrations caused by floating point precision errors need to be handled beforehand
  #These discrepancies are noted in cc_makeMeta (when metadata is 'complete' type) and are handled in cc_getSyn
  

  #Input requirements:
  #Note: these values are heuristics determined from distributions of all collected dense mode data
  #norm_data-----------------------Output from cc_norm
  #dr_data-------------------------Output from cc_getDR
  #cutoff_sd-----------------------SD of %CD among replicates for all concs of each drug. Above this value concs are flagged
  #cutoff_mono---------------------Monotonicity threshold. Should be a negative value because it will be used in the low-to-high evaluation...
  #................................monotonicity is evaluated bi-directionally. When low-to-high, values less than mono_cutoff are flagged....
  #................................when evaluated high-to-low, values greater than the absolute value of this number will be flagged.
  #................................Ex if cutoff_mono = -10:  lowest conc (nth_conc==10) is 50% cell death, then 2nd lowest conc (nth_conc==9) is 5% cell death, nth_conc==9 gets flagged (evaluated low-to-high)
  #................................NOTE: when evaluating this high-to-low, nth_conc==10 will also get flagged. A bi-directional evaluation means flagging always occurs in pairs.
  #cutoff_resid--------------------Residual threshold. Concs are flagged if their corresponding response is greater than (the ABSOLUTE VALUE) of this number
  
  
  #Determine whether data is dense or sparse mode
  if(norm_data[[1]]$data_mode == "sparse"){
    data_mode <- "sparse"
    print("Detected data type is SPARSE MODE. Proceeding with sparse mode QC...")
  } else {
    data_mode <- "dense"
    print("Detected data type is DENSE MODE. Proceeding with dense mode QC...")}
  
  
  #Print message about running this function after the other functions:
  print("NOTE: This function is meant to be run AFTER synergy calculations")
  
  
  if(data_mode=="dense"){
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle dense mode data------------------------
    #
    #=================================================================
    #=================================================================
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Loop over each element of norm_data
    for(i in seq_along(norm_data)){
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Extract necessary data:
      ref_df <- norm_data[[i]]$ref_df
      
      d1_name <- norm_data[[i]]$d1_name
      d2_name <- norm_data[[i]]$d2_name
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #---------Begin SD flagging----------
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      #Make dataframes of the SD of%CD of each single-agent across replicates
      #Each should contain only 10 rows (one for each concentration).
      
      sd_df_drug1 <- 
        ref_df %>%
        filter(drug2_conc==0) %>% #Filter to the drug1 single-agent
        filter(!(drug1_conc==0 & drug2_conc==0)) %>% #Remove the 0x0 row
        group_by(drug1_conc) %>%
        mutate("single_agent_sd"=sd(replicate_perc_cell_death)) %>%
        ungroup() %>%
        select(drug1_name,
               combo_id,
               drug1_conc,
               single_agent_sd) %>%
        unique()  %>% #Keep unique values (which are the same across replicates)
        mutate("drug2_name"=d2_name,
               "drug2_conc"=0)
       
      
      sd_df_drug2 <- 
        ref_df %>%
        filter(drug1_conc==0) %>% #Filter to the drug2 single-agent
        filter(!(drug1_conc==0 & drug2_conc==0)) %>% #Remove the 0x0 row
        group_by(drug2_conc) %>%
        mutate("single_agent_sd"=sd(replicate_perc_cell_death)) %>%
        ungroup() %>%
        select(drug2_name,
               combo_id,
               drug2_conc,
               single_agent_sd) %>%
        unique() %>% #Keep unique values (which are the same across replicates)
        mutate("drug1_name"=d1_name,
               "drug1_conc"=0)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Stack the sd_df's on top of each other into 1 dataframe
      single_agent_sd_df <- 
        rbind(sd_df_drug1, sd_df_drug2)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Join the SD dataframes back to the ref_df:
      #The single-agents will have values in this column
      #The combinations will have NA
      ref_df <- 
        ref_df %>%
        left_join(.,
                  single_agent_sd_df,
                  by=c("drug1_name", 
                       "combo_id",
                       "drug1_conc",
                       "drug2_name",
                       "drug2_conc"))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Identify the concentrations with SD > cutoff_sd, for each drug
      flagged_sd_concs_d1 <- 
        ref_df %>%
        filter(drug2_conc==0, #Filters to drug1-only
               single_agent_sd > cutoff_sd) %>%
        pull(drug1_conc) %>%
        unique()
      
      flagged_sd_concs_d2 <- 
        ref_df %>%
        filter(drug1_conc==0, #Filters to drug2-only
               single_agent_sd > cutoff_sd) %>%
        pull(drug2_conc) %>%
        unique()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Flag the concentrations that exceed the cutoff_sd
      #This is done for single-agents and combinations
      #The idea is that if a single-agent concentration is flagged, then all combinations containing that concentration should also be flagged
      ref_df <- 
        ref_df %>%
        mutate("flagged_sd" = case_when(
          
          drug1_conc %in% flagged_sd_concs_d1 ~ 1,
          drug2_conc %in% flagged_sd_concs_d2 ~ 1,
          TRUE ~ 0
        ))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #-----------End SD flagging----------
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #----Begin monotonicity flagging-----
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #The DR curve should be monotonically increasing left-to-right (low-to-high concentrations)
      #We will flag any concentrations that do not follow this pattern
      
      #Because the lowest and highest concentrations don't have points before/after them, we check monotonicity in both directions
      #One consequence of this is that any given flagged concentration will also have a neighboring concentration that gets flagged
      
      #The 'calc_monotonicity' helper function will evaluate the difference in %CD between each concentration and its preceding concentration (bidirectionally)
      #The required input is a dataframe with columns 'drug_conc' and 'perc_cell_death'
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Grab monotonicity for each single-agent
      mono_df_drug1 <- 
        ref_df %>%
        filter(drug2_conc   == 0,                  #Filters to drug1-only
               !(drug1_conc == 0 & drug2_conc==0), #Removes the 0x0 row
               replicate    == 1) %>%              #Only uses 1 replicate (fine since we're using averaged %CD as perc_cell_death)
        select(drug1_conc, perc_cell_death) %>%
        rename("drug_conc" = "drug1_conc")  %>%    #Because calc_monotonicity function needs columns 'drug_conc' and 'perc_cell_death'
        calc_monotonicity(.) %>%
        rename("drug1_conc"= "drug_conc")   %>%    #Set conc column back to original name of drug1_conc
        mutate("drug2_conc"= 0,                    #Helps with merging back to ref_df
               "drug1_name" = d1_name,
               "drug2_name" = d2_name,
               "combo_id"= paste0(d1_name, "_", d2_name))                    
        
      
      mono_df_drug2 <- 
        ref_df %>%
        filter(drug1_conc == 0,                    #Filters to drug2-only
               !(drug1_conc == 0 & drug2_conc==0), #Removes the 0x0 row
               replicate    == 1) %>%              #Only uses 1 replicate (fine since we're using averaged %CD as perc_cell_death)
        select(drug2_conc, perc_cell_death) %>%
        rename("drug_conc" = "drug2_conc")  %>%    #Because calc_monotonicity function needs columns 'drug_conc' and 'perc_cell_death'
        calc_monotonicity(.) %>%
        rename("drug2_conc"= "drug_conc")   %>%    #Set conc column back to original name of drug2_conc
        mutate("drug1_conc"= 0,                    #Helps with merging back to ref_df
               "drug1_name" = d1_name,
               "drug2_name" = d2_name,
               "combo_id"= paste0(d1_name, "_", d2_name))    
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Stack the mono_df's on top of each other into 1 dataframe:
      single_agent_mono_df <- 
        rbind(mono_df_drug1, mono_df_drug2) %>%
        select(-perc_cell_death)  #Drop this column so it doesn't get duplicated in the next section's join step
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Join the single_agent_mono_df back to the ref_df
      #This contains the monotonicity values for the single-agents
      ref_df <- 
        ref_df %>%
        left_join(.,
                  single_agent_mono_df,
                  by=c("drug1_name", 
                       "combo_id",
                       "drug1_conc",
                       "drug2_name",
                       "drug2_conc"))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Identify the concentrations violating the cutoff_mono
      #Remember, we will evaluate low-to-high and high-to-low
      flagged_mono_concs_d1 <- 
        mono_df_drug1 %>%
        filter((delta_val_low_to_high < cutoff_mono) | (delta_val_high_to_low > abs(cutoff_mono))) %>%
        pull(drug1_conc)
      
      flagged_mono_concs_d2 <- 
        mono_df_drug2 %>%
        filter((delta_val_low_to_high < cutoff_mono) | (delta_val_high_to_low > abs(cutoff_mono))) %>%
        pull(drug2_conc)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Flag the concentrations that violate the cutoff_mono
      ref_df <- 
        ref_df %>%
        mutate("flagged_mono" = case_when(
          
          drug1_conc %in% flagged_mono_concs_d1 ~ 1,
          drug2_conc %in% flagged_mono_concs_d2 ~ 1,
          TRUE ~ 0
        ))
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #------End monotonicity flagging-----
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #------Begin residual flagging-------
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #The point is to obtain the residuals of the dose-response model fit for each concentration point, for each drug
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #First extract the observed data for the single-agents (we will add residual info to this)
      resid_df_drug1 <- 
        ref_df %>%
        filter(drug2_conc   == 0,                  #Filters to drug1-only
               !(drug1_conc == 0 & drug2_conc==0), #Removes the 0x0 row
               replicate    == 1) %>%              #Only uses 1 replicate (fine since we're using averaged %CD as perc_cell_death)
        select(drug1_name, drug1_conc, perc_cell_death) 
      
      resid_df_drug2 <- 
        ref_df %>%
        filter(drug1_conc == 0,                    #Filters to drug2-only
               !(drug1_conc == 0 & drug2_conc==0), #Removes the 0x0 row
               replicate    == 1)  %>%             #Only uses 1 replicate (fine since we're using averaged %CD as perc_cell_death)
        select(drug2_name, drug2_conc, perc_cell_death)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #We need to construct a name that can correctly index into the predicted_results_list
      #For dense mode, this is the drug name followed by the name of the index in norm_data:
      tmp_name_d1 <- paste0(d1_name, "_", names(norm_data[i]))
      tmp_name_d2 <- paste0(d2_name, "_", names(norm_data[i]))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #If the DR model for Drug1 WAS successfully fit:
      if(dr_data$predicted_results_list[[tmp_name_d1]]$model_output$model_success==TRUE){
        
        #------------------------------------
        #Extract predicted DR data for drug1
        predicted_data_drug1 <- 
          dr_data$predicted_results_list[[tmp_name_d1]]$model_output$pred_df
        
        #Find the nearest predicted concentration for each observed concentration (drug1) :
        resid_df_drug1 <- 
          resid_df_drug1 %>%
          rowwise() %>%
          mutate("nearest_dose" =
                   predicted_data_drug1 %>%
                   filter(abs(dose - drug1_conc) == min(abs(dose - drug1_conc))) %>%
                   pull(dose)) %>%
          mutate("nearest_pred" = 
                   predicted_data_drug1 %>%
                   filter(abs(dose - drug1_conc) == min(abs(dose - drug1_conc))) %>%
                   pull(response)) %>%
          mutate("abs_res" = abs(perc_cell_death - nearest_pred)) %>% #Calculate the absolute value of the residual
          mutate("drug2_name" = d2_name,
                 "drug2_conc" = 0)
      } else {
        #------------------------------------
        #Else, If the DR model for drug1 WAS NOT successfuly fit:
        resid_df_drug1 <- 
          resid_df_drug1 %>%
          mutate("nearest_dose" = NA,
                 "nearest_pred" = NA,
                 "abs_res"      = NA,     #Set the residual to NA
                 #The NA should tell us of curve fitting failure
                 #Note, this means it will NOT get flagged as exceeding 
                 #Justification is that drug1 may be a totally flat line, zero response,
                 #...but could be synergistic when combined with a certain other drug
                 "drug2_name"   = d2_name,
                 "drug2_conc"   = 0)
        
      }
      #------------------------------------
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #If the DR model for Drug2 WAS successfully fit:
      if(dr_data$predicted_results_list[[tmp_name_d2]]$model_output$model_success==TRUE){
        
        #------------------------------------
        #Extract predicted DR data for Drug2
        predicted_data_drug2 <-
          dr_data$predicted_results_list[[tmp_name_d2]]$model_output$pred_df
        
        #Find the nearest predicted concentration for each observed concentration (Drug2) :
        resid_df_drug2 <- 
          resid_df_drug2 %>%
          rowwise() %>%
          mutate("nearest_dose" =
                   predicted_data_drug2 %>%
                   filter(abs(dose - drug2_conc) == min(abs(dose - drug2_conc))) %>%
                   pull(dose)) %>%
          mutate("nearest_pred" = 
                   predicted_data_drug2 %>%
                   filter(abs(dose - drug2_conc) == min(abs(dose - drug2_conc))) %>%
                   pull(response)) %>%
          mutate("abs_res" = abs(perc_cell_death - nearest_pred)) %>% #Calculate the absolute value of the residual
          mutate("drug1_name" = d1_name,
                 "drug1_conc" = 0)
      } else {
        #------------------------------------
        #Else, If the DR model for Drug2 WAS NOT successfully fit:
        resid_df_drug2 <- 
          resid_df_drug2 %>%
          mutate("nearest_dose" = NA,
                 "nearest_pred" = NA,
                 "abs_res"      = NA,     #Set the residual to NA
                 #The NA should tell us of curve fitting failure
                 #Note, this means it will NOT get flagged as exceeding 
                 #Justification is that drug1 may be a totally flat line, zero response,
                 #...but could be synergistic when combined with a certain other drug
                 "drug1_name"   = d1_name,
                 "drug1_conc"   = 0)
        #------------------------------------
      }
        
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Stack the resid_dfs to combine into one dataframe
      single_agent_resid_df <- 
        bind_rows(resid_df_drug1, resid_df_drug2) %>%
        select(-perc_cell_death) #Remove the %CD column (otherwise it woudl get duplicated with the next join step)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Join the single_agent_resid_df back to the ref_df
      ref_df <- 
        ref_df %>%
        left_join(., 
                  single_agent_resid_df,
                  by=c("drug1_name", "drug1_conc", "drug2_name", "drug2_conc"))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Identify the concentrations which exceed the cutoff_resid
      flagged_resid_concs_d1 <- 
        ref_df %>%
        filter(drug2_conc   == 0,
               !(drug1_conc == 0 & drug2_conc == 0),
               replicate    == 1) %>%
        filter(abs_res > cutoff_resid) %>%
        pull(drug1_conc)
      
      flagged_resid_concs_d2 <- 
        ref_df %>%
        filter(drug1_conc == 0,
               !(drug1_conc == 0 & drug2_conc == 0),
               replicate    == 1) %>%
        filter(abs_res > cutoff_resid) %>%
        pull(drug2_conc)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      #Flag the concentrations in ref_df that exceed the cutoff_resid
      ref_df <- 
        ref_df %>%
        mutate("flagged_resid" = case_when(
          drug1_conc %in% flagged_resid_concs_d1 ~ 1,
          drug2_conc %in% flagged_resid_concs_d2 ~ 1,
          TRUE ~ 0
        ))

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #--------End residual flagging-------
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Final flagging
      #Make a final_flag column which flags rows where any of the other flagged columns have a value of 1:
      ref_df <- 
        ref_df %>%
        mutate("flagged_final" = pmax(flagged_sd,
                                      flagged_mono,
                                      flagged_resid))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Add some summary statistics to summarize the synergy with and without QC metrics
      ref_df <- 
        ref_df %>%
        mutate(
          
          #Mean/Median/Max of all synergy values:
          #We filter out single-agents (which are 0 by definition) and filter to replicate==1 to avoid counting the same values multiple times
          "mean_syn"       = mean  (bliss_synergy[(!drug1_conc==0 & !drug2_conc==0 & replicate==1)]), #This excludes the single-agents and filters to replicate==1
          "median_syn"     = median(bliss_synergy[(!drug1_conc==0 & !drug2_conc==0 & replicate==1)]),
          "max_syn"        = max   (bliss_synergy[(!drug1_conc==0 & !drug2_conc==0 & replicate==1)]),
          "second_max_syn" = sort  (bliss_synergy[(!drug1_conc==0 & !drug2_conc==0 & replicate==1)], decreasing=TRUE)[2],

          
          
          #Mean/Median/Max of all synergy values with flagged_final==0:
          "mean_syn_adj"       = mean  (bliss_synergy[(!drug1_conc==0 & !drug2_conc==0  & replicate==1 & flagged_final==0)]),
          "median_syn_adj"     = median(bliss_synergy[(!drug1_conc==0 & !drug2_conc==0  & replicate==1 & flagged_final==0)]),
          "max_syn_adj"        = max   (bliss_synergy[(!drug1_conc==0 & !drug2_conc==0  & replicate==1 & flagged_final==0)]),
          "second_max_syn_adj" = sort  (bliss_synergy[(!drug1_conc==0 & !drug2_conc==0  & replicate==1 & flagged_final==0)], decreasing=TRUE)[2],
        )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate the mean_perc_cell_death
      #Note: this is the mean of ALL 100 squares, so it includes the single-agents
      ref_df <- 
        ref_df %>%
        mutate("mean_perc_cell_death" = mean(perc_cell_death)) #Reminder that 'perc_cell_death' is the mean across replicates for a given concentration pair
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate Moran's I to check spatial autocorrelation among the synergy matrix:
      
      #The input for the calc_moran helper function is a 10x10 synergy matrix (drug1 concs as cols, drug2 concs as rows)
      #So first we need to get the input 10x10 synergy matrix: (There should only be 100 values)
      mydf <- 
        ref_df %>%
        #Filter out single-agents so we're left just with combination rows
        filter( !(drug1_conc == 0 & drug2_conc == 0), 
                !drug1_conc==0,
                !drug2_conc==0,
                replicate==1) %>%
        select(drug1_conc, drug2_conc, bliss_synergy)  #Select only the columns we need
      
      
      synergy_matrix <- 
        reshape2::dcast(mydf, 
                        drug2_conc ~ drug1_conc,           #NOTE: use 'drug2_conc ~ drug1_conc' so that drug2 is rows, drug1 is cols
                        value.var = "bliss_synergy") %>%
        #Make first column the row names:
        column_to_rownames(var="drug2_conc") %>%
        #flip rows top-to-bottom:
        .[order(nrow(.):1),]
      
      
      # Convert row names to numeric if they are not
      rownames(synergy_matrix) <- as.numeric(rownames(synergy_matrix))
      
      # Ensure that all columns are numeric
      synergy_matrix[] <- lapply(synergy_matrix, function(x) as.numeric(as.character(x)))
      
      # Convert synergy_matrix to a numeric matrix
      synergy_matrix <- as.matrix(synergy_matrix)
      
      #Calculate moran's I (and the standardized version of it)
      moran_res <- calc_moran(synergy_matrix)
      
      #Add moran results to the data frame
      ref_df <- 
        ref_df %>%
        mutate("moran_i" = moran_res$moran_i,
               "moran_i_STND"= moran_res$moran_i_STND)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Update ref_df
      norm_data[[i]]$ref_df <- ref_df
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Print progress
      percent_complete <- (i / length(norm_data)) * 100
      cat(sprintf("\rProgress: %.2f%% complete", percent_complete))
      flush.console()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
    }#End loop over norm_data
      
    
    
  } else {
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle sparse mode data-----------------------
    #
    #=================================================================
    #=================================================================
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Loop over each element of norm_data
    for(i in seq_along(norm_data)){
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Extract necessary data:
      ref_df <- norm_data[[i]]$ref_df
      
      d1_name <- norm_data[[i]]$d1_name
      d2_name <- norm_data[[i]]$d2_name
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #---------Begin SD flagging----------
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Make dataframes of the SD of%CD of each single-agent across replicates
      #Each should contain only 10 rows (one for each concentration). 
      
      sd_df_drug1 <- 
        ref_df %>%
        filter(plate_type=="single_agent",
               drug1_name==d1_name) %>%
        group_by(drug1_conc) %>%
        mutate("single_agent_sd"=sd(replicate_perc_cell_death)) %>%
        select(drug1_name, 
               combo_id, 
               drug1_conc,
               single_agent_sd) %>%
        unique() #Keep unique values (which are the same across replicates)
    
      sd_df_drug2 <- 
        ref_df %>%
        filter(plate_type=="single_agent",
               drug1_name==d2_name) %>% #NOTE the use of drug1_name==d2_name here, due to structure of sparse mode single-agents
        group_by(drug1_conc) %>%
        mutate("single_agent_sd"=sd(replicate_perc_cell_death)) %>%
        select(drug1_name,
               combo_id,
               drug1_conc,
               single_agent_sd) %>%
        unique() #Keep unique values (which are the same across replicates)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Stack the sd_df's on top of each other into 1 dataframe
      single_agent_sd_df <- 
        rbind(sd_df_drug1, sd_df_drug2)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Join the SD dataframes back to the ref_df:
      ref_df <- 
        ref_df %>%
        left_join(.,
                  single_agent_sd_df,
                  by=c("drug1_name", "combo_id", "drug1_conc")) 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Identify the concentrations with SD > cutoff_sd, for each drug
      flagged_sd_concs_d1 <- 
        ref_df %>%
        filter(plate_type=="single_agent",
               drug1_name==d1_name,
               single_agent_sd > cutoff_sd) %>%
        pull(drug1_conc) %>%
        unique()
      
      flagged_sd_concs_d2 <- 
        ref_df %>%
        filter(plate_type=="single_agent",
               drug1_name==d2_name, #Note use of 'd2_name' here
               single_agent_sd > cutoff_sd) %>%
        pull(drug1_conc) %>%
        unique()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Flag the concentrations that exceed the cutoff_sd
      #This is done for single-agents and combinations
      #The idea is that if a single-agent concentration is flagged, then all combinations containing that concentration should also be flagged

      ref_df <- 
        ref_df %>%
        mutate(
          "flagged_sd" = case_when(
            
            #First flag the single-agent concentrations that exceed the cutoff_sd
            drug1_name==d1_name & plate_type=="single_agent" & drug1_conc %in% flagged_sd_concs_d1 ~1,
            drug1_name==d2_name & plate_type=="single_agent" & drug1_conc %in% flagged_sd_concs_d2 ~1, #Note use of 'drug1_name/conc'
            
            #Now flag the combinations that contain the flagged single-agent concentrations
            drug1_name==d1_name & plate_type=="combination" & drug1_conc %in% flagged_sd_concs_d1 ~1,
            drug2_name==d2_name & plate_type=="combination" & drug2_conc %in% flagged_sd_concs_d2 ~1,
            TRUE ~ 0))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #-----------End SD flagging----------
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #----Begin monotonicity flagging-----
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #The DR curve should be monotonically increasing left-to-right (low-to-high concentrations)
      #We will flag any concentrations that do not follow this pattern
      
      #Because the lowest and highest concentrations don't have points before/after them, we check monotonicity in both directions
      #One consequence of this is that any given flagged concentration will also have a neighboring concentration that gets flagged
      
      #The 'calc_monotonicity' helper function will evaluate the difference in %CD between each concentration and its preceding concentration (bidirectionally)
      #The required input is a dataframe with columns 'drug_conc' and 'perc_cell_death'
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Grab monotonicity for each single-agent and the diagonal
      mono_df_drug1 <- 
        ref_df %>%
        filter(plate_type == "single_agent", #Filter to just the single-agents
               drug1_name == d1_name,        #Filter to drug1-only
               replicate  == 1) %>%          #Filter to just replicate 1 (fine since we are using averaged perc_cell_death)
        select(drug1_conc, perc_cell_death) %>%
        rename("drug_conc"="drug1_conc") %>% #Rename the drug1_conc column to 'drug_conc' (required for calc_monotonicity function)
        calc_monotonicity(.) %>%
        mutate("drug1_name"=d1_name,
               "combo_id"=paste0(d1_name, "_", NA))
      
      
      mono_df_drug2 <- 
        ref_df %>%
        filter(plate_type == "single_agent", #Filter to just the single-agents
               drug1_name == d2_name,        #Filter to drug2-only , note use of d2_name!
               replicate  == 1) %>%          #Filter to just replicate 1 (fine since we are using averaged perc_cell_death)
        select(drug1_conc, perc_cell_death) %>%
        rename("drug_conc"="drug1_conc") %>% #Rename the drug1_conc column to 'drug_conc' (required for calc_monotonicity function)
        calc_monotonicity(.) %>%
        mutate("drug1_name"=d2_name, #Note use of d2_name here
               "combo_id" = paste0(d2_name, "_", NA)) 
      
      
      mono_df_diag <- 
        ref_df %>%
        filter(plate_type == "combination") %>%
        filter(!is.na(raw_val)) %>% #This is a way to filter to just the measured combinations (predicted, non-measured rows are NA here)
        select(drug1_conc, perc_cell_death) %>%  #NOTE: Since the diagonal/measured combo values are always screened 1:1 relative concentration...
                                                  #......it does not matter which drug's concentrations we use here
                                                  #......even if the drugs have differing concentrations, the relative/'nth_conc' is always the same
                                                  #......Still, let it be convention to use drug1_conc here, 
                                                  #......especially since we are merging these together with the single-agent DF's, which use 'drug1_name' and 'drug1_conc'
        rename("drug_conc"="drug1_conc") %>%      #because calc_monotonicity needs the column name to just be 'drug_conc'
        unique() %>%
        calc_monotonicity(.) %>%
        mutate("drug1_name"=d1_name,
               "combo_id"=paste0(d1_name,"_",d2_name))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Stack the mono_df's on top of each other into 1 dataframe:
      single_agent_mono_df <-
        rbind(mono_df_drug1, mono_df_drug2, mono_df_diag) %>%
        rename("drug1_conc"="drug_conc")  %>% #Set conc column back to original name of drug1_conc
        select(-perc_cell_death) %>% #Drop this column so it doesn't get duplicated in the next section's join step
        
        #We need to make a temporary column just for joining this data to the ref_df
        #The reason is that the single-agent df's would join just fine, but
        #The mono_df_diag would have a problem. It would join in the correct places, 
        #But it would also join in all of the places where drug1_conc, drug1_name, and combo_id match BUT are not from diagonal measurements
        #This temporary column will be used to join only to measured values
        mutate("actually_measured"="yes")
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Join the single_agent_mono_df back to the ref_df
      #This contains the monotonicity values for the single-agents
      ref_df <- 
        ref_df %>%
        #Add the temporary 'actually_measured' column as described above, for joining
        mutate("actually_measured"=case_when(
          !is.na(raw_val) ~ "yes",
          TRUE ~ "no"
        )) %>%
        left_join(.,
                  single_agent_mono_df,
                  by=c("drug1_name", "drug1_conc", "combo_id", "actually_measured")) %>%
        
        select(-actually_measured) #Drop this column now that we are done with it
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Identify the concentrations violating the cutoff_mono
      #Remember, we will evaluate low-to-high and high-to-low
      flagged_mono_concs_d1 <- 
        ref_df %>%
        filter(plate_type == "single_agent",
               drug1_name == d1_name,
               replicate  == 1) %>%
        filter((delta_val_low_to_high < cutoff_mono) | (delta_val_high_to_low > abs(cutoff_mono))) %>%
        pull(drug1_conc)
      
      flagged_mono_concs_d2 <- 
        ref_df %>%
        filter(plate_type == "single_agent",
               drug1_name == d2_name, #Note d2_name
               replicate  == 1) %>%
        filter((delta_val_low_to_high < cutoff_mono) | (delta_val_high_to_low > abs(cutoff_mono))) %>%
        pull(drug1_conc)
      
      flagged_mono_concs_diag <-
        ref_df %>%
        filter(plate_type == "combination",
               !is.na(raw_val)) %>%
        filter((delta_val_low_to_high < cutoff_mono) | (delta_val_high_to_low > abs(cutoff_mono))) %>%
        pull(drug1_conc)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Flag the concentrations that violate the cutoff_mono
      ref_df <- 
        ref_df %>%
        mutate("flagged_mono" = case_when(
          
          #First flag the single-agent concentrations that violate the cutoff_mono
          drug1_name==d1_name & plate_type=="single_agent" & drug1_conc %in% flagged_mono_concs_d1 ~1,
          drug1_name==d2_name & plate_type=="single_agent" & drug1_conc %in% flagged_mono_concs_d2 ~1, #Note use of d2_name here
          
          #Now flag the combinations that contain the flagged single-agent concentrations which violated monnotonicity:
          drug1_name==d1_name & plate_type=="combination" & drug1_conc %in% flagged_mono_concs_d1  ~1,
          drug2_name==d2_name & plate_type=="combination" & drug2_conc %in% flagged_mono_concs_d2  ~1,
          
          #Now flag the combinations on the diagonal which contain the flagged diagonal concs that violated monotonicity:
          #We use '!is.na(raw_val)' to subset to just the diagonal values aka the measured data
          #Then we keep track of the flagged diagonal concs using 'drug1_conc' which we selected previously in making mono_df_diag
          #Note that the diagonal values are always 1:1 relative concentration, so keeping track by drug1_conc is fine as long as drug1_conc was used previously to make mono_df_diag
          (plate_type=="combination" & !is.na(raw_val) & drug1_conc %in% flagged_mono_concs_diag) ~1,
          
          TRUE ~ 0))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #------End monotonicity flagging-----
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #------Begin residual flagging-------
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #The point is to obtain the residuals of the dose-response model fit for each concentration point, for each drug

      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #First extract the observed data for the single-agents (we will add residual info to this)
      resid_df_drug1 <- 
        ref_df %>%
        filter(plate_type == "single_agent", 
               drug1_name == d1_name,
               replicate  == 1) %>%
        select(drug1_name, drug1_conc, perc_cell_death)
      
      resid_df_drug2 <- 
        ref_df %>%
        filter(plate_type == "single_agent", 
               drug1_name == d2_name, #Note d2_name
               replicate  == 1) %>%
        select(drug1_name, drug1_conc, perc_cell_death)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      

      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #If the DR model for Drug1 WAS successfully fit:
      if(dr_data$predicted_results_list[[d1_name]]$model_success==TRUE){
        
        
        #------------------------------------
        #Extract predicted DR data for drug1
        predicted_data_drug1 <- 
          dr_data$predicted_results_list[[d1_name]]$pred_df
        
        #Find the nearest predicted concentration for each observed concentration (drug1) :
        resid_df_drug1 <- 
          resid_df_drug1 %>%
          rowwise() %>%
          mutate("nearest_dose" = 
                   predicted_data_drug1 %>%
                   filter(abs(dose - drug1_conc) == min(abs(dose - drug1_conc))) %>%
                   pull(dose))  %>%
          mutate("nearest_pred" = 
                   predicted_data_drug1 %>%
                   filter(abs(dose - drug1_conc) == min(abs(dose - drug1_conc))) %>%
                   pull(response)) %>%
          mutate("abs_res" = abs(perc_cell_death - nearest_pred)) %>% #Calculate the absolute value of the residual
          mutate("combo_id" = paste0(d1_name, "_", NA)) #This will help with correctly joining back to ref_df
      } else {
        #------------------------------------
        #Else, If the DR model for drug1 WAS NOT successfuly fit:
        resid_df_drug1 <- 
          resid_df_drug1 %>%
          mutate("nearest_dose" = NA,
                 "nearest_pred" = NA,
                 "abs_res"      = NA,
                 "combo_id" = paste0(d1_name, "_", NA)) #This will help with correctly joining back to ref_df
      }
        #------------------------------------
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #If the DR model for Drug2 WAS successfully fit:
      if(dr_data$predicted_results_list[[d2_name]]$model_success==TRUE){
        
        
        #------------------------------------
        #Extract predicted DR data for Drug2
        predicted_data_drug2 <-
          dr_data$predicted_results_list[[d2_name]]$pred_df
        
        #Find the nearest predicted concentration for each observed concentration (Drug2) :
        resid_df_drug2 <- 
          resid_df_drug2 %>%
          rowwise() %>%
          mutate("nearest_dose" = 
                   predicted_data_drug2 %>%
                   filter(abs(dose - drug1_conc) == min(abs(dose - drug1_conc))) %>%
                   pull(dose))  %>%
          mutate("nearest_pred" = 
                   predicted_data_drug2 %>%
                   filter(abs(dose - drug1_conc) == min(abs(dose - drug1_conc))) %>%
                   pull(response)) %>%
          mutate("abs_res" = abs(perc_cell_death - nearest_pred)) %>%  #Calculate the absolute value of the residual
          mutate("combo_id" = paste0(d2_name, "_", NA)) #This will help with correctly joining back to ref_df
      } else {
        #------------------------------------
        #Else, If the DR model for Drug2 WAS NOT successfully fit:
        resid_df_drug2 <-
          resid_df_drug2 %>%
          mutate("nearest_dose" = NA,
                 "nearest_pred" = NA,
                 "abs_res"      = NA,
                 "combo_id" = paste0(d2_name, "_", NA)) #This will help with correctly joining back to ref_df
        #------------------------------------
      }
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Stack the resid_dfs to combine into one dataframe
      single_agent_resid_df <- 
        bind_rows(resid_df_drug1, resid_df_drug2) %>%
        select(-perc_cell_death) #Remove the %CD column (otherwise it woudl get duplicated with the next join step)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Join the single_agent_resid_df back to the ref_df
      ref_df <- 
        ref_df %>%
        left_join(.,
                  single_agent_resid_df,
                  by = c("drug1_name", "drug1_conc", "combo_id"))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Identify the concentrations which exceed the cutoff_resid
      flagged_resid_concs_d1 <- 
        ref_df %>%
        filter(plate_type == "single_agent",
               drug1_name == d1_name,
               replicate  == 1,
               abs_res > cutoff_resid) %>%
        pull(drug1_conc)
      
      flagged_resid_concs_d2 <- 
        ref_df %>%
        filter(plate_type == "single_agent",
               drug1_name == d2_name,        #Note d2_name
               replicate  == 1,
               abs_res > cutoff_resid) %>%
        pull(drug1_conc)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Flag the concentrations in ref_df that exceed the cutoff_resid
      ref_df <- 
        ref_df %>%
        mutate("flagged_resid" = case_when(
          
          #First flag the single-agent concentrations that exceed the cutoff_resid
          plate_type == "single_agent" & drug1_name == d1_name & drug1_conc %in% flagged_resid_concs_d1 ~1,
          plate_type == "single_agent" & drug1_name == d2_name & drug1_conc %in% flagged_resid_concs_d2 ~1,
          
          #Now flag the combinations that contain the flagged single_agent concentrations
          plate_type == "combination" & drug1_name == d1_name & drug1_conc %in% flagged_resid_concs_d1 ~1,
          plate_type == "combination" & drug2_name == d2_name & drug2_conc %in% flagged_resid_concs_d2 ~1, #Note use drug2_name AND drug2_conc for combinations!
          TRUE ~ 0
        ))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #--------End residual flagging-------
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Final flagging
      #Make a final_flag column which flags rows where any of the other flagged columns have a value of 1:
      ref_df <- 
        ref_df %>%
        mutate("flagged_final" = pmax(flagged_sd,
                                      flagged_mono,
                                      flagged_resid))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Add some summary statistics to summarize the synergy with and without QC metrics
      ref_df <- 
        ref_df %>%
        mutate(
          
          #Mean/Median/Max of all synergy values:
          "mean_syn"       = mean  (bliss_synergy[plate_type=="combination"]), #We use plate_type==combination to exclude all single-agent synergy values (which are always 0)
          "median_syn"     = median(bliss_synergy[plate_type=="combination"]),
          "max_syn"        = max   (bliss_synergy[plate_type=="combination"]),
          "second_max_syn" = sort  (bliss_synergy[plate_type=="combination"], decreasing=TRUE)[2],
          
          
          #Mean/Media/Max of all synergy values with flagged_final==0:
          "mean_syn_adj"       = mean  (bliss_synergy[plate_type=="combination" & flagged_final==0]),
          "median_syn_adj"     = median(bliss_synergy[plate_type=="combination" & flagged_final==0]),
          "max_syn_adj"        = max   (bliss_synergy[plate_type=="combination" & flagged_final==0]),
          "second_max_syn_adj" = sort  (bliss_synergy[plate_type=="combination" & flagged_final==0], decreasing=TRUE)[2],
          
          #Mean/Median/Max of diagonal/observed values' synergy:
          "mean_syn_diag"       = mean  (bliss_synergy[plate_type=="combination" & !is.na(nth_conc)]), #This filters to just the 10 measured values
          "median_syn_diag"     = median(bliss_synergy[plate_type=="combination" & !is.na(nth_conc)]),
          "max_syn_diag"        = max   (bliss_synergy[plate_type=="combination" & !is.na(nth_conc)]),
          "second_max_syn_diag" = sort  (bliss_synergy[plate_type=="combination" & !is.na(nth_conc)], decreasing=TRUE)[2],
          
          #Mean/Median/Max of diagonal/observed values' synergy with flagged_final==0:
          "mean_syn_diag_adj"       = mean  (bliss_synergy[plate_type=="combination" & !is.na(nth_conc) & flagged_final==0]),
          "median_syn_diag_adj"     = median(bliss_synergy[plate_type=="combination" & !is.na(nth_conc) & flagged_final==0]),
          "max_syn_diag_adj"        = max   (bliss_synergy[plate_type=="combination" & !is.na(nth_conc) & flagged_final==0]),
          "second_max_syn_diag_adj" = sort  (bliss_synergy[plate_type=="combination" & !is.na(nth_conc) & flagged_final==0], decreasing=TRUE)[2]
        )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #In some cases ALL TEN observed/diagonal values can be flagged
      #For mean/median this results in NA, but for max this results in -Inf, which is bad for downstream analysis

      #Manually change any -Inf values to NA here:
      ref_df <- 
        ref_df %>%
        mutate(
                "max_syn_diag"     = ifelse(is.infinite(max_syn_diag), NA, max_syn_diag),
                "max_syn_diag_adj" = ifelse(is.infinite(max_syn_diag_adj), NA, max_syn_diag_adj))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate the mean_perc_cell_death
      #Note: this is the mean of ALL 100 squares, so it includes the single-agents
      ref_df <- 
        ref_df %>%
        mutate("mean_perc_cell_death" = mean(perc_cell_death)) #Reminder that 'perc_cell_death' is the mean across replicates for a given concentration pair
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Add a column that tells the number of flagged diagonal squares
      #This will be later when filtering (ex: only keeping combos that have >= 3 surviving/non-flagged diagonal squares)
      num_flagged_diag <- 
        ref_df %>%
        filter(plate_type=="combination" & !is.na(raw_val)) %>%
        pull(flagged_final) %>%
        sum()
      
      ref_df <- 
        ref_df %>%
        mutate("num_flagged_diag" = num_flagged_diag)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate Moran's I to check spatial autocorrelation among the synergy matrix:
      
      #The input for the calc_moran helper function is a 10x10 synergy matrix (drug1 concs as cols, drug2 concs as rows)
      #So first we need to get the input 10x10 synergy matrix:
      mydf <- 
        ref_df %>%
        filter(plate_type=="combination", replicate==1) %>% #Subset to just replicate 1
        select(drug1_conc, drug2_conc, bliss_synergy)       #Select only the columns we need
      
      synergy_matrix <- 
        reshape2::dcast(mydf, 
                        drug2_conc ~ drug1_conc,           #NOTE: use 'drug2_conc ~ drug1_conc' so that drug2 is rows, drug1 is cols
                        value.var = "bliss_synergy") %>%
        #Make first column the row names:
        column_to_rownames(var="drug2_conc") %>%
        #flip rows top-to-bottom:
        .[order(nrow(.):1),] %>%
        as.matrix()
      
      #Calculate moran's I (and the standardized version of it)
      moran_res <- calc_moran(synergy_matrix)
      
      #Add moran results to the data frame
      ref_df <- 
        ref_df %>%
        mutate("moran_i" = moran_res$moran_i,
               "moran_i_STND"= moran_res$moran_i_STND)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
        

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Update ref_df
      norm_data[[i]]$ref_df <- ref_df
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Print progress
      percent_complete <- (i / length(norm_data)) * 100
      cat(sprintf("\rProgress: %.2f%% complete", percent_complete))
      flush.console()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
    }#End loop over each element of norm_data
    
    
  }#End handle sparse mode data
  
  
  return(norm_data)
  
  
}#End function
#--------------------------------------------------------------
