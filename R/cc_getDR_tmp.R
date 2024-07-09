#--------------------------------------------------------------
#Function to output dose-response data
cc_getDR <- function(norm_data,
                     name_char_limit = 16,
                     rounding_value_IC50 = 2){
  
  #Input requirements:
  #norm_data-----------------------Output from cc_norm
  #name_char_limit-----------------Character limit for drug names in plots
  #rounding_value_IC50-------------Rounding the drug concentrations for the IC50 in the plot title
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Determine whether data is dense or sparse mode
  if(norm_data[[1]]$data_mode == "sparse"){
    data_mode <- "sparse"
    print("Detected data type is SPARSE MODE. Proceeding with sparse mode dose-response modeling...")
  } else {
    data_mode <- "dense"
    print("Detected data type is DENSE MODE. Proceeding with dense mode dose-response modeling...")}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  if(data_mode=="dense"){
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle dense mode data------------------------
    #
    #=================================================================
    #=================================================================
  
  
    
    #In dense mode, each plate is 1 combination.
    #If plate1 has combo A_B and plate2 has combo A_C, then single-agent A has been tested twice
    #(And it will get two dose-response curves)
    
    #Initialize empty lists for drug1 and drug2
    #These will get combined downstream
    d1_list <- list()
    d2_list <- list()
    
    
    for(i in seq_along(norm_data)){
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Extract the current combination data
      drug1_df_norm <- norm_data[[i]]$drug1_df_norm
      drug2_df_norm <- norm_data[[i]]$drug2_df_norm
      
      d1_name <- drug1_df_norm$drug_name %>% unique()
      d2_name <- drug2_df_norm$drug_name %>% unique()
      
      units <- norm_data[[i]]$units
      
      plate_name <- names(norm_data[i])
      
      sample_name <- norm_data[[i]]$sample_name
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Make unique list entry names for each single-agent
      #This is necessary because two elements can't have the same name ('DrugA' twice, for example)
      d1_entry_name <- paste0(d1_name, "_", plate_name)
      d2_entry_name <- paste0(d2_name, "_", plate_name)
      
      #d1 model
      d1_model_output <- 
        drug1_df_norm %>%
        filter(!drug_conc==0) %>% #Remove where drug_conc==0
        dr_func(.,
                dose_column     = "drug_conc",
                response_column = "perc_cell_death") #note 'perc_cell_death' is the mean of replicates
      
      #Deposit the results into the d1_list
      #NOTE: We modify some names so they can be accessed downstream in a generalizable way 
      #drug1_df_norm becomes drug_df_norm, d1_name becomes drug_name, etc...
      d1_list[[d1_entry_name]] <- list("drug_df_norm"=drug1_df_norm,   
                                       "drug_name"=d1_name,
                                       "units"=units,
                                       "plate_name"=plate_name,
                                       "sample_name"=sample_name,
                                       "model_output"=d1_model_output)
      
      
      #d2 model
      d2_model_output <-
        drug2_df_norm %>%
        filter(!drug_conc==0) %>% #Remove where drug_conc==0
        dr_func(.,
                dose_column     = "drug_conc",
                response_column = "perc_cell_death") #note 'perc_cell_death' is the mean of replicates
      
      #Deposit the results into the d1_list
      #NOTE: We modify some names so they can be accessed downstream in a generalizable way 
      #drug1_df_norm becomes drug_df_norm, d1_name becomes drug_name, etc...
      d2_list[[d2_entry_name]] <- list("drug_df_norm"=drug2_df_norm,    #Note: we name it 'drug_df_norm' not 'drug2_df_norm' (To be accessed in a generalizable way downstream)
                                       "drug_name"=d2_name,
                                       "units"=units,
                                       "plate_name"=plate_name,
                                       "sample_name"=sample_name,
                                       "model_output"=d2_model_output)  #Note: we name it 'model_output' not 'drug2_model_output' (To be accessed in a generalizable way downstream)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    }
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Combine d1_list and d2_list into a single list:
    predicted_results_list <- c(d1_list, d2_list)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize empty list for DR plots, and a df for IC50s
    plots_list <- list()
    IC50_df <- data.frame("drug_name"=character(),
                          "IC50"=numeric(),
                          "units"=character(),
                          "sample_name"=character(),
                          "plate_name"=character(),
                          "comment"=character())
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
 
    for(i in seq_along(predicted_results_list)){
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Set up a 'comment' that will be included in the IC50_df, and will be filled in if IC50 cannot be calculated
      my_comment <- NA_real_
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Extract necessary data for plotting the DR curve
      #This is handled differently depending on whether the model was successful
      if(predicted_results_list[[i]]$model_output$model_success==TRUE){
        
        
        #If model was successful....
        
        #Subset the original data
        og_data <- predicted_results_list[[i]]$drug_df_norm
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Get the predicted data
        pred_df   <- predicted_results_list[[i]]$model_output$pred_df
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        drug_name <- predicted_results_list[[i]]$drug_name
        units <- predicted_results_list[[i]]$units
        plate_name <- predicted_results_list[[i]]$plate_name
        sample_name <- predicted_results_list[[i]]$sample_name
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Extract the IC50 value
        #True IC50 is the closest concentration to y=50, NOT simply half of the max response (aka half max)
        drug_IC50 <- pred_df[which.min(abs(pred_df$response-50)), "dose"]
        
        #If the maximum response is <50% cell death, set IC50 to NA
        #This will clearly indicate when the maximum concentration of a drug does not kill 50% of the cells
        if(max(pred_df$response) < 50){
          drug_IC50 <- NA
          my_comment <- "Max response < 50% cell death"
          }
        
        #Add the IC50 to the IC50_df
        IC50_df <- 
          rbind(IC50_df,
                data.frame("drug_name"=drug_name,
                           "IC50"=drug_IC50,
                           "units"=units,
                           "plate_name"=plate_name,
                           "sample_name"=sample_name,
                           "comment"=my_comment))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make labels for plot
        
        #Cap the drug  name at chosen character limit
        drug_name_capped <- drug_name %>% substr(., 1, name_char_limit)
        
        #Make plot title
        my_plot_title <- paste0(drug_name_capped, " IC50 = ", round(drug_IC50, rounding_value_IC50), units)
        
        #Make x-axis label
        my_xlab <- paste0("Log10(", drug_name, ", ", units, ")")
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make the y-axis limits for the plot based on the data
        
        #Figure out which data produces the lower min value, either the predicted data or the original response data
        min_of_pred_df <- min(pred_df$response)
        min_of_og_data <- min(og_data$perc_cell_death)
        smaller_min    <- min(min_of_pred_df, min_of_og_data)
        
        #Take the lower min value between the two to figure out what the ymin of the plot is. If greater than 0, set to 0.
        my_ymin <- ifelse(smaller_min>0, 0, smaller_min)
        
        
        #Figure out which data produces the higher max value, either the predicted data or the original response data
        max_of_pred_df <- max(pred_df$response)
        max_of_og_data <- max(og_data$perc_cell_death)
        bigger_max     <- max(max_of_pred_df, max_of_og_data)
        
        #Take the bigger max value between the two to figure out what the ymax of the plot is. If less than 100, set to 100. 
        my_ymax <- ifelse(bigger_max<100, 100, bigger_max)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
 
        
        
      } else {
        
        #Else handle if model was NOT successful...
        
        
        #Subset the original data
        og_data <- predicted_results_list[[i]]$drug_df_norm
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        drug_name <- predicted_results_list[[i]]$drug_name
        units <- predicted_results_list[[i]]$units
        plate_name <- predicted_results_list[[i]]$plate_name
        sample_name <- predicted_results_list[[i]]$sample_name
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        drug_name <- predicted_results_list[[i]]$drug_name
        
        #Set pred_df to just zero 
        pred_df <- 
          data.frame("dose"=1,
                     "response"=1)
        
        #Set drug_IC50 to NA
        drug_IC50 <- NA
        
        #Set my_comment to the model failure message
        my_comment <- "model could not be fit to data"
        
        #Add the (NA) IC50 to IC50_df
        IC50_df <- 
          rbind(IC50_df,
                data.frame("drug_name"=drug_name,
                           "IC50"=drug_IC50,
                           "units"=units,
                           "plate_name"=plate_name,
                           "sample_name"=sample_name,
                           "comment"=my_comment))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make labels for plot
        #In case the model failed, this should just plot the data with no curve (and IC50 will be NA)
        
        #Cap the drug  name at chosen character limit
        drug_name_capped <- drug_name %>% substr(., 1, name_char_limit)
        
        #Make plot title
        my_plot_title <- paste0(drug_name_capped, " IC50 = ", round(drug_IC50, rounding_value_IC50), units)
        
        #Make x-axis label
        my_xlab <- paste0("Log10(", drug_name, ", ", units, ")")
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make the y-axis limits for the plot based on the data
        
        #Figure out which data produces the lower min value, either the predicted data or the original response data
        min_of_pred_df <- min(pred_df$response)
        min_of_og_data <- min(og_data$perc_cell_death)
        smaller_min    <- min(min_of_pred_df, min_of_og_data)
        
        #Take the lower min value between the two to figure out what the ymin of the plot is. If greater than 0, set to 0.
        my_ymin <- ifelse(smaller_min>0, 0, smaller_min)
        
        
        #Figure out which data produces the higher max value, either the predicted data or the original response data
        max_of_pred_df <- max(pred_df$response)
        max_of_og_data <- max(og_data$perc_cell_death)
        bigger_max     <- max(max_of_pred_df, max_of_og_data)
        
        #Take the bigger max value between the two to figure out what the ymax of the plot is. If less than 100, set to 100. 
        my_ymax <- ifelse(bigger_max<100, 100, bigger_max)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
    }#End if model was not successful
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Generate the plot
      dr_plot <- 
        
        #Add the replicate data points
        og_data %>%
        filter(!drug_conc==0) %>% #Remove where drug_conc==0
        ggplot(aes(x=drug_conc, y=replicate_perc_cell_death)) +
        geom_point(color="gray80") +
        
        #Add the dose-response curve data
        geom_line(data=pred_df,
                  aes(x=dose, y=response, color=response), linewidth=1.5) +
        scale_color_gradient2(low= "#6991a7ff", 
                              mid= "#ffc4b0ff",
                              high="#f25651ff",
                              midpoint=50) +
        
        #Add the mean data points
        geom_point(data=og_data %>% filter(!drug_conc==0),
                   aes(x=drug_conc, y=perc_cell_death, fill=perc_cell_death),
                   size=3, shape=21, color="black", stroke=1.5)+
        scale_fill_gradient2(low= "#6991a7ff", 
                             mid= "#ffc4b0ff",
                             high="#f25651ff",
                             midpoint=50) +
        
        scale_x_log10() +
        geom_hline(yintercept = 50, linetype="dashed") +
        #Change the y-axis limits: if lowest value is >0, set to 0; if highest value is <100, set to 100  
        ylim(my_ymin, my_ymax) +
        xlab(my_xlab) +  
        ylab("%Cell Death") +
        ggtitle(my_plot_title) +
        theme_minimal() +
        theme(axis.ticks = element_line(color="black"),
              legend.position = "none",
              plot.title=element_text(hjust=0.5)) +
        geom_vline(xintercept=drug_IC50, linetype="dashed")
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #If my_comment is anything other than NA, then no IC50 was calcualted
      #So make the plot_title just the drug name without an 'IC50=NA' label
      if(!is.na(my_comment)){
        dr_plot <- 
          dr_plot + ggtitle(drug_name_capped)
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      

      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Deposit into list
      plots_list[[i]] <- dr_plot
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      

    }#End loop along predicted results
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Name elements of plots_list the same as predicted_results_list
    names(plots_list) <- names(predicted_results_list)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    #Return data
    return(list(
      "plots_list" = plots_list, 
      "IC50_df" = IC50_df,
      "predicted_results_list" = predicted_results_list))
    
    

    
  } else {
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle sparse mode data-----------------------
    #
    #=================================================================
    #=================================================================
    
    #In sparse mode, each single-agent is recycled each time it is used in a new combination
    #We do not want to redundantly make DR curves for the same single-agent
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Get the unique drugs from norm_data
    
    #Create a list where each element contains the drug names from each combination in norm_data
    drug_names_list <- lapply(norm_data, function(combo) {
      c(combo$d1_name, combo$d2_name)
    })
    
    #Convert the list of drug names into a single vector
    all_drug_names <- unlist(drug_names_list)
    
    #Get the unique drug names from the vector
    unique_drugs <- unique(all_drug_names)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Create a list to hold drug data from norm_data
    dr_list <- vector("list", length(unique_drugs))
    names(dr_list) <- unique_drugs
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize an empty dataframe for each element in the single_agent_list
    for(i in unique_drugs){
      dr_list[[i]] <- 
        data.frame("drug_conc"=numeric(),
                   "perc_cell_death"=numeric(),
                   "drug_name"=character())
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Populate dr_list with combined data for each unique drug
    for(i in seq_along(norm_data)){
      
      #Extract the current combination data
      tmplist <- norm_data[i][[1]]
      
      #Extract drug names and units
      d1_name <- tmplist$d1_name
      d2_name <- tmplist$d2_name
      units   <- tmplist$units
      sample_name <- tmplist$sample_name
      
      #Extract drug1 data (and add units)
      drug1_data <- 
        tmplist$drug1_df_norm %>%
        mutate("units"=units,
               "sample_name"=sample_name)
      
      #Extract drug2 data (and add units)
      drug2_data <- 
        tmplist$drug2_df_norm %>%
        mutate("units"=units,
               "sample_name"=sample_name)
      
      #Add the drug data to the dr_list
      dr_list[[d1_name]] <- 
        rbind(dr_list[[d1_name]],
              drug1_data)
      
      dr_list[[d2_name]] <- 
        rbind(dr_list[[d2_name]],
              drug2_data)
      
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Clean the dr_list elements so only unique rows remain
    for(i in unique_drugs){
      dr_list[[i]] <- 
        dr_list[[i]] %>% 
        unique()
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize empty list for predicted results to go into
    predicted_results_list <- list()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Generate the dose-response curve outputs
    for(i in seq_along(dr_list)){
      
      model_output <- 
        dr_list[[i]] %>%
        dr_func(., 
                dose_column     = "drug_conc",
                response_column = "perc_cell_death") #note 'perc_cell_death' is the mean of replicates
      
      predicted_results_list[[i]] <- model_output
    }
    
    names(predicted_results_list) <- names(dr_list)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize empty list for DR plots, and a df for IC50s
    plots_list <- list()
    IC50_df <- data.frame("drug_name"=character(),
                          "IC50"=numeric(),
                          "units"=character(),
                          "sample_name"=character(),
                          "comment"=character())
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #Loop along the predicted results (each is a single-agent)
    for(i in seq_along(predicted_results_list)){
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Set up a 'comment' that will be included in the IC50_df, and will be filled in if IC50 cannot be calculated
      my_comment <- NA_real_
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Extract data that is independent of the model success
      og_data <- dr_list[[i]] %>% unique()    #Original collected data
      drug_name <- unique(og_data$drug_name)  #Drug name
      units <- unique(og_data$units)  
      sample_name <- unique(og_data$sample_name)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #Extract necessary data for plotting the DR curve
      #This is handled differently depending on whether the model was successful
      if(predicted_results_list[[i]]$model_success==TRUE){
        
        
        #If model was successful....
        
        
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       #Get the predicted data
        pred_df   <- predicted_results_list[[i]]$pred_df
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Extract the IC50 value
        #True IC50 is the closest concentration to y=50, NOT simply half of the max response (aka half max)
        drug_IC50 <- pred_df[which.min(abs(pred_df$response-50)), "dose"]
        
        #If the maximum response is <50% cell death, set IC50 to NA
        #This will clearly indicate when the maximum concentration of a drug does not kill 50% of the cells
        if(max(pred_df$response) < 50){
          drug_IC50 <- NA
          my_comment <- "Max response < 50% cell death"
        }
        
        #Add the IC50 to the IC50_df
        IC50_df <- 
          rbind(IC50_df,
                data.frame("drug_name"=drug_name,
                           "IC50"=drug_IC50,
                           "units"=units,
                           "sample_name"=sample_name, 
                           "comment"=my_comment))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make labels for plot
        
        #Cap the drug  name at chosen character limit
        drug_name_capped <- drug_name %>% substr(., 1, name_char_limit)
        
        #Make plot title
        my_plot_title <- paste0(drug_name_capped, " IC50 = ", round(drug_IC50, rounding_value_IC50), units)
        
        #Make x-axis label
        my_xlab <- paste0("Log10(", drug_name, ", ", units, ")")
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make the y-axis limits for the plot based on the data
        
        #Figure out which data produces the lower min value, either the predicted data or the original response data
        min_of_pred_df <- min(pred_df$response)
        min_of_og_data <- min(og_data$perc_cell_death)
        smaller_min    <- min(min_of_pred_df, min_of_og_data)
        
        #Take the lower min value between the two to figure out what the ymin of the plot is. If greater than 0, set to 0.
        my_ymin <- ifelse(smaller_min>0, 0, smaller_min)
        
        
        #Figure out which data produces the higher max value, either the predicted data or the original response data
        max_of_pred_df <- max(pred_df$response)
        max_of_og_data <- max(og_data$perc_cell_death)
        bigger_max     <- max(max_of_pred_df, max_of_og_data)
        
        #Take the bigger max value between the two to figure out what the ymax of the plot is. If less than 100, set to 100. 
        my_ymax <- ifelse(bigger_max<100, 100, bigger_max)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
      } else {
        
        #Else handle if model was NOT successful
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Set pred_df to just zero 
        pred_df <- 
          data.frame("dose"=1,
                     "response"=1)
        
        #Set drug_IC50 to NA
        drug_IC50 <- NA
        
        #Set my_comment to the model failure message
        my_comment <- "model could not be fit to data"
        
        #Add the (NA) IC50 to IC50_df
        IC50_df <- 
          rbind(IC50_df,
                data.frame("drug_name"=drug_name,
                           "IC50"=drug_IC50,
                           "units"=units,
                           "sample_name"=sample_name,
                           "comment"=my_comment))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make labels for plot
        #In case the model failed, this should just plot the data with no curve (and IC50 will be NA)
        
        #Cap the drug  name at chosen character limit
        drug_name_capped <- drug_name %>% substr(., 1, name_char_limit)
        
        #Make plot title
        my_plot_title <- paste0(drug_name_capped, " IC50=", round(drug_IC50, rounding_value_IC50), units)
        
        #Make x-axis label
        my_xlab <- paste0("Log10(", drug_name, ", ", units, ")")
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make the y-axis limits for the plot based on the data
        
        #Figure out which data produces the lower min value, either the predicted data or the original response data
        min_of_pred_df <- min(pred_df$response)
        min_of_og_data <- min(og_data$perc_cell_death)
        smaller_min    <- min(min_of_pred_df, min_of_og_data)
        
        #Take the lower min value between the two to figure out what the ymin of the plot is. If greater than 0, set to 0.
        my_ymin <- ifelse(smaller_min>0, 0, smaller_min)
        
        
        #Figure out which data produces the higher max value, either the predicted data or the original response data
        max_of_pred_df <- max(pred_df$response)
        max_of_og_data <- max(og_data$perc_cell_death)
        bigger_max     <- max(max_of_pred_df, max_of_og_data)
        
        #Take the bigger max value between the two to figure out what the ymax of the plot is. If less than 100, set to 100. 
        my_ymax <- ifelse(bigger_max<100, 100, bigger_max)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
        
      }#End if model was not successful
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Generate the plot
        dr_plot <- 
          
          #Add the replicate data points
          og_data %>%
          ggplot(aes(x=drug_conc, y=replicate_perc_cell_death)) +
          geom_point(color="gray80") +
          
          #Add the dose-response curve data
          geom_line(data=pred_df,
                    aes(x=dose, y=response, color=response), linewidth=1.5) +
          scale_color_gradient2(low= "#6991a7ff", 
                                mid= "#ffc4b0ff",
                                high="#f25651ff",
                                midpoint=50) +
          
          #Add the mean data points
          geom_point(data=og_data,
                     aes(x=drug_conc, y=perc_cell_death, fill=perc_cell_death),
                     size=3, shape=21, color="black", stroke=1.5)+
          scale_fill_gradient2(low= "#6991a7ff", 
                               mid= "#ffc4b0ff",
                               high="#f25651ff",
                               midpoint=50) +
          
          scale_x_log10() +
          geom_hline(yintercept = 50, linetype="dashed") +
          #Change the y-axis limits: if lowest value is >0, set to 0; if highest value is <100, set to 100  
          ylim(my_ymin, my_ymax) +
          xlab(my_xlab) +  
          ylab("%Cell Death") +
          ggtitle(my_plot_title) +
          theme_minimal() +
          theme(axis.ticks = element_line(color="black"),
                legend.position = "none",
                plot.title=element_text(hjust=0.5)) +
          geom_vline(xintercept=drug_IC50, linetype="dashed")
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #If my_comment is anything other than NA, then no IC50 was calcualted
        #So make the plot_title just the drug name without an 'IC50=NA' label
        if(!is.na(my_comment)){
          dr_plot <- 
            dr_plot + ggtitle(drug_name_capped)
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Deposit into list
        plots_list[[i]] <- dr_plot
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
    }#End loop along predicted results
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Use drug names as the element names of this list
    #NOTE: In dense mode, 'DrugA' may be used twice (DrugA_DrugB and DrugA_DrugC, for example)
    #......so in dense mode, the elements of the list are named as [drug_name]_[plate_name]
    #......in sparse mode however, it is sufficient to just use the drug names as the names of this list
    #......because single-agents are not re-measured like they are in dense moe
    names(plots_list) <- names(dr_list)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return(list(
      "plots_list" = plots_list, 
      "IC50_df" = IC50_df,
      "predicted_results_list" = predicted_results_list))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  }#End handling sparse mode data
  
}#End function
#--------------------------------------------------------------
