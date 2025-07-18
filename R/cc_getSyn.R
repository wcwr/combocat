#--------------------------------------------------------------
#' Calculate Synergy Scores
#'
#' This function calculates Bliss and Loewe synergy scores for both dense and sparse mode data, based on normalized input data.
#' The results are added to the reference data frame.
#'
#' @param norm_data List. Output from `cc_norm`, containing normalized data and metadata.
#' @param dr_data List. Output from `cc_getDR`, containing dose-response and IC50 data.
#' @param conc_rounding_factor Numeric. Number of decimal places to round concentrations to. Default is `6`. 
#'   Rounding is important for merging drugs with slightly different concentrations due to floating-point errors
#'   This is particularly important due to floating point errors that may arise from `cc_makeMeta`
#'
#' @export
#--------------------------------------------------------------




#--------------------------------------------------------------
#Function to calculate synergy
cc_getSyn <- function(norm_data,
                      dr_data,
                      conc_rounding_factor=6){
  
  
  #Input requirements:
  #norm_data-----------------------Output from cc_norm
  #dr_data-------------------------Output from cc_getDR
  #conc_rounding_factor------------Number of decimal places to round concentrations to. Default is 6.
  #................................Note: this is important for merging drugs which may have slightly different concentrations due 
  #................................to floating point decimal errors from calculating missing concentrations in cc_makeMeta
  #.
  #Note: This function returns the norm_data but adds synergy results to the reference dataframe (ref_df)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Determine whether data is dense or sparse mode
  if(norm_data[[1]]$data_mode == "sparse"){
    data_mode <- "sparse"
    print("Detected data type is SPARSE MODE. Calculating synergy for sparse mode data...")
  } else {
    data_mode <- "dense"
    print("Detected data type is DENSE MODE. Calculating synergy for dense mode data...")}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
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
      #Extract out the ref_df
      ref_df <- norm_data[[i]]$ref_df
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Round the concentrations
      #NOTE: This is critical for merging concentrations which might be very slightly different due to floating point errors...
      #For example: 
      #DrugA has a concentration of 0.111111112
      #DrugB was missing at this concentration, but through divisions (in cc_makeMeta) we calculate it to be 0.111111113
      #They now cannot be joined due to not being exactly the same, so we round them to make them the same
      ref_df <- 
        ref_df %>%
        mutate("drug1_conc"=round(drug1_conc, conc_rounding_factor),
               "drug2_conc"=round(drug2_conc, conc_rounding_factor))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Add single-agent effects to the data frame for corresponding concentrations
      ref_df <-
        ref_df %>%
        mutate("single_effect_drug1" = ifelse(drug2_conc == 0, perc_cell_death, NA),
               "single_effect_drug2" = ifelse(drug1_conc == 0, perc_cell_death, NA))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate single-agent effects separately
      single_effect_drug1_df <- 
        ref_df %>%
        filter(drug2_conc == 0) %>%
        select(drug1_conc, single_effect_drug1 = perc_cell_death) %>%
        distinct()
      
      single_effect_drug2_df <- 
        ref_df %>%
        filter(drug1_conc == 0) %>%
        select(drug2_conc, single_effect_drug2 = perc_cell_death) %>%
        distinct()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Merge single-agent effects back into the main dataframe and handle column names
      ref_df <- 
        ref_df %>%
        left_join(single_effect_drug1_df, by = "drug1_conc", suffix = c("", ".single1")) %>%
        left_join(single_effect_drug2_df, by = "drug2_conc", suffix = c("", ".single2"))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Explicitly select and rename the columns to ensure correct column names
      ref_df <- 
        ref_df %>%
        mutate(single_effect_drug1 = coalesce(single_effect_drug1, single_effect_drug1.single1),
               single_effect_drug2 = coalesce(single_effect_drug2, single_effect_drug2.single2)) %>%
        select(-ends_with(".single1"), -ends_with(".single2"))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate Bliss expectation with the correct column names
      ref_df <- 
        ref_df %>%
        mutate(bliss_exp = 100 - ((100 - single_effect_drug1) * (100 - single_effect_drug2) / 100))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate synergy (as Excess over Bliss)
      ref_df <- 
        ref_df %>%
        mutate(bliss_synergy = perc_cell_death - bliss_exp)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Adjust bliss_synergy to be zero when drug1_conc or drug2_conc is zero
      ref_df <- 
        ref_df %>%
        mutate(bliss_synergy = ifelse(drug1_conc == 0 | drug2_conc == 0, 0, bliss_synergy))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate Loewe synergy
      
      #Extract the dose-response data for the single-agents
      #Get the correct indices of dr_data
      #This will ensure the correct format is used to get the right dr data. The format for a given index is <drug_name>_<filename>
      drug1_dr_index <- paste0(unique(ref_df$drug1_name), "_", unique(ref_df$filename))
      drug2_dr_index <- paste0(unique(ref_df$drug2_name), "_", unique(ref_df$filename))
      
      #Now we can extract the single-agent dose-response dataframes for these drugs:
      drug1_dr_df <- dr_data$predicted_results_list[[drug1_dr_index]]$model_output$pred_df
      drug2_dr_df <- dr_data$predicted_results_list[[drug2_dr_index]]$model_output$pred_df
      
      #Note: We can't calculate Loewe in cases where one or both of the DR curves failed - so set to NA
      if(is.null(drug1_dr_df) || is.null(drug2_dr_df)){
        ref_df$loewe_synergy <- NA_real_
        message(paste0(
          "Warning: Loewe synergy cannot be calculated for ", 
          unique(ref_df$drug1_name), " and ", 
          unique(ref_df$drug2_name), 
          " (file: ", unique(ref_df$filename), 
          ") due to missing dose-response data."
        ))
        
      } else {
        
        #Otherwise if both DR curves are valid:
        #Calculate Loewe synergy using the `calc_loewe` helper function
        ref_df <- 
          calc_loewe(
            reference_df = ref_df,
            drug1_conc_column = "drug1_conc",
            drug2_conc_column = "drug2_conc",
            response_column   = "perc_cell_death",
            drug1_dr = drug1_dr_df,
            drug2_dr = drug2_dr_df
          )
      }
      
      #Adjust loewe_synergy to be zero when drug1_conc or drug2_conc is zero
      ref_df <- 
        ref_df %>%
        mutate(loewe_synergy = ifelse(drug1_conc == 0 | drug2_conc == 0, NA_real_, loewe_synergy))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Replace the original ref_df with the updated version containing synergy data
      norm_data[[i]]$ref_df <- ref_df
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
      #Extract out the ref_df and drug names
      ref_df <- norm_data[[i]]$ref_df
      d1_name <- norm_data[[i]]$d1_name
      d2_name <- norm_data[[i]]$d2_name
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Round the concentrations
      #NOTE: This is critical for merging concentrations which might be very slightly different due to floating point errors...
      #For example: 
      #DrugA has a concentration of 0.111111112
      #DrugB was missing at this concentration, but through divisions (in cc_makeMeta) we calculate it to be 0.111111113
      #They now cannot be joined due to not being exactly the same, so we round them to make them the same
      ref_df <- 
        ref_df %>%
        mutate("drug1_conc"=round(drug1_conc, conc_rounding_factor),
               "drug2_conc"=round(drug2_conc, conc_rounding_factor))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Split the ref_df into single_agent and combination plates
      #This is the first step to adding the single-agent effects of each drug (which we later use to calculate synergy)
      ref_df_single_agent_plates <- 
        ref_df %>%
        filter(plate_type == "single_agent")
      
      ref_df_combination_plates <- 
        ref_df %>%
        filter(plate_type == "combination")
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Add single-agent effects to the single_agents data frame for corresponding concentrations
      #'plate_type' is used to determine single-agent effects in sparse mode
      ref_df_single_agent_plates <-
        ref_df_single_agent_plates %>%
        mutate("single_effect_drug1" = ifelse(drug1_name == d1_name, perc_cell_death, NA),
               "single_effect_drug2" = ifelse(drug1_name == d2_name, perc_cell_death, NA)) #Note use of 'drug1_name' column for d2_name
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Set the single_effect_drug1/2 values to 0 for single_agent plates in cases where the value is current NA:
      #This completes the single_effect_drug1/2 columns for single-agents
      ref_df_single_agent_plates <- 
        ref_df_single_agent_plates %>% 
        mutate("single_effect_drug1" = ifelse(is.na(single_effect_drug1), 0, single_effect_drug1),
               "single_effect_drug2" = ifelse(is.na(single_effect_drug2), 0, single_effect_drug2))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Begin figuring out single_effects for combination plates
      #First, obtain single-agent effects separately
      single_effects_d1 <- 
        ref_df_single_agent_plates %>%
        filter(drug1_name == d1_name) %>%
        select(drug1_conc, perc_cell_death) %>%
        rename("drug_conc"=drug1_conc,
               "single_effect_drug1"=perc_cell_death) %>%
        distinct()
      
      single_effects_d2 <- 
        ref_df_single_agent_plates %>%
        filter(drug1_name == d2_name) %>%
        select(drug1_conc, perc_cell_death) %>%
        rename("drug_conc"=drug1_conc,
               "single_effect_drug2"=perc_cell_death) %>%
        distinct()
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Add the single_effects_d1/d2 to the combination plates
      ref_df_combination_plates <-
        ref_df_combination_plates%>%
        left_join(.,
                  single_effects_d1,
                  by = c("drug1_conc" = "drug_conc")) %>% #Merge in single_effect_drug1 values
        left_join(.,
                  single_effects_d2,
                  by = c("drug2_conc" = "drug_conc"))     #Merge in single_effect_drug2 values
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Merge the single_agent and combination plate types together, now that all single-agent effects have been added
      ref_df <- 
        bind_rows(ref_df_single_agent_plates, ref_df_combination_plates)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate Bliss expectation with the correct column names
      ref_df <- 
        ref_df %>%
        mutate(bliss_exp = 100 - ((100 - single_effect_drug1) * (100 - single_effect_drug2) / 100))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate synergy (as Excess over Bliss)
      ref_df <- 
        ref_df %>%
        mutate(bliss_synergy = perc_cell_death - bliss_exp)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Adjust bliss_synergy to be zero for any single-agents
      ref_df <- 
        ref_df %>%
        mutate(bliss_synergy = ifelse(plate_type=="single_agent", 0, bliss_synergy))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Calculate Loewe synergy
      
      #Extract the dose-response data for the single-agents
      #Note this is indexed in a slightly different way than dense mode
      drug1_dr_df <- dr_data$predicted_results_list[[d1_name]]$pred_df
      drug2_dr_df <- dr_data$predicted_results_list[[d2_name]]$pred_df
      
      #Note: We can't calculate Loewe in cases where one or both of the DR curves failed - so set to NA
      if(is.null(drug1_dr_df) || is.null(drug2_dr_df)){
        ref_df$loewe_synergy <- NA_real_
        message(paste0(
          "Warning: Loewe synergy cannot be calculated for combination ", 
          d1_name, " and ", d2_name, 
          " due to missing dose-response data."
        ))
        
      } else {
        
        #Otherwise if both DR curves are valid:
        #Calculate Loewe synergy using the `calc_loewe` helper function
        ref_df <- 
          calc_loewe(
            reference_df = ref_df,
            drug1_conc_column = "drug1_conc",
            drug2_conc_column = "drug2_conc",
            response_column   = "perc_cell_death",
            drug1_dr = drug1_dr_df,
            drug2_dr = drug2_dr_df
          )
      }
      
      #Adjust loewe_synergy to be zero when drug1_conc or drug2_conc is zero
      ref_df <- 
        ref_df %>%
        mutate(loewe_synergy = ifelse(drug1_conc == 0 | drug2_conc == 0, NA_real_, loewe_synergy))
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Replace the original ref_df with the updated version containing synergy data
      norm_data[[i]]$ref_df <- ref_df
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
    }#End loop over norm_data
    
  }#End handle sparse mode data
  
  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Return results
  return(norm_data)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
}#End function
#--------------------------------------------------------------