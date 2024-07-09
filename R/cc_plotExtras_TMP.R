#--------------------------------------------------------------
#Function to make various supplemental plots
cc_plotExtras <- function(syn_data,
                          color_midpoint = 20,
                          rounding_value = 4,
                          name_char_limit = 16
                          ){
  
  
  #Input requirements:
  #syn_data---------------------Output from cc_getSyn
  #color_midpoint---------------Midpoint for the color scale
  #rounding_value---------------Number of decimal places to round the drug concentrations for plotting
  #name_char_limit--------------Maximum number of characters for the drug names for plotting
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Initialize empty lists:
  syn_barplots_list <- list()
  cd_scatterplots_list <- list()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Loop over each element in syn_data
  for(i in seq_along(syn_data)){
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Subset the ref_df
    ref_df <- syn_data[[i]]$ref_df
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Get the maximum synergy value
    max_syn_value <- 
      max(ref_df$bliss_synergy[!ref_df$bliss_synergy==0]) %>% 
      unique() #Unique in case of ties
    
    #Get the expected %Cell Death corresponding to the max synergy
    bliss_exp_value <- 
      ref_df %>%
      filter(bliss_synergy == max_syn_value) %>%
      pull(bliss_exp) %>%
      unique()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Extract necessary values for plotting
    
    #Get the units
    units <- 
      ref_df %>%
      pull(units) %>%
      unique() %>%
      .[1]
    
    #Get the drug names
    d1_name <- 
      syn_data[[i]]$d1_name
    
    d2_name <- 
      syn_data[[i]]$d2_name
    
    #Get the drug concentration of each drug that corresponds to the max synergy
    d1_best_conc <-
      ref_df %>%
      filter(bliss_synergy == max_syn_value) %>%
      pull(drug1_conc) %>%
      unique() %>%
      round(., rounding_value) %>%
      .[1] #In case of ties or value of 0, take the first entry
    
    d2_best_conc <- 
      ref_df %>%
      filter(bliss_synergy == max_syn_value) %>%
      pull(drug2_conc) %>%
      unique() %>%
      round(., rounding_value) %>%
      .[1] #In case of ties or value of 0, take the first entry
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Make labels for the plot -- e.g. "DrugA (0.1 uM)"
    drug1_label <- 
      paste0(d1_name, " (", d1_best_conc, " ", units, ")")
    
    drug2_label <-
      paste0(d2_name, " (", d2_best_conc, " ", units, ")")
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Make a long format dataframe for the barplot
    ref_df_long <- 
      ref_df %>%
      filter(bliss_synergy == max_syn_value) %>%
      filter(replicate == 1) %>% #Filter to first replicate because we're dealing with perc_cell_death which is the average of replicates
      select(single_effect_drug1, single_effect_drug2, perc_cell_death) %>%
      pivot_longer(cols = everything(), names_to = "treatment", values_to = "perc_cell_death")
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Adjust the factor levels in 'treatment' column (to keep consistent ordering of the bars in the plot)
    ref_df_long$treatment <- 
      factor(ref_df_long$treatment, 
             levels = c("single_effect_drug1", "single_effect_drug2", "perc_cell_death"),
             labels = c(drug1_label, drug2_label, "Combined"))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Make the barplot
    syn_barplot <- 
      ref_df_long %>%
      ggplot(aes(x=treatment, y=perc_cell_death, fill=treatment)) +
      geom_hline(yintercept = bliss_exp_value, linetype="dashed", color="red") +
      geom_bar(stat="identity", width=0.5, color="black") +
      theme_minimal() +
      scale_fill_manual(values = c("gray70", "gray70", "forestgreen")) +
      xlab("") +
      ylab("%Cell Death") +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Deposit into list
    syn_barplots_list[[i]] <- syn_barplot
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Make the observed vs. expected cell death plots (cd_scatterplot)
    cd_scatterplot <- 
      ref_df %>%
      ggplot(aes(x=bliss_exp, y=perc_cell_death, color=bliss_synergy)) +
      geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed") +  
      geom_point() +
      scale_color_gradient2("Bliss Synergy",
                            low=  c("#6991a7ff", "#a1b5beff", "#bbc7c9ff", "#e8e6ddff"),
                            mid=  c("#ffc4b0ff", "#ffa491ff"),
                            high= c("#ff8f77ff", "#f5795eff", "#f25651ff", "#ee3d3cff"),
                            midpoint=color_midpoint) +
      theme_minimal() +
      xlab("%Cell Death [expected]") +
      ylab("%Cell Death [observed]")
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
   
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Deposit cd_scatterplot into list
    cd_scatterplots_list[[i]] <- cd_scatterplot
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
  }#End loop over syn_data
  
  names(syn_barplots_list) <- names(syn_data)
  names(cd_scatterplots_list) <- names(syn_data)
  
  return(list(
    "syn_barplots_list" = syn_barplots_list,
    "cd_scatterplots_list" = cd_scatterplots_list
  ))
  
 

}#End function
#--------------------------------------------------------------