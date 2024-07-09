#--------------------------------------------------------------
#Function to plot a matrix of cell death or synergy
cc_plotMat <- function(input_data,                            
                       plotting_variable = "perc_cell_death",
                       color_midpoint = 50, 
                       rounding_value = 6,  
                       name_char_limit = 16){
                       
  
  
  #Input requirements:
  #input_data-----------------------Output from cc_norm (for plotting cell death) or cc_getSyn (for plotting cell death or synergy)
  #plotting_variable----------------Variable to plot ('perc_cell_death' or 'bliss_synergy')
  #color_midpoint-------------------Midpoint of the color scale (50 for cell death)
  #rounding_value-------------------Number of digits to round concentration values in plots
  #name_char_limit------------------Character limit for drug names in plots
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Determine whether data is dense or sparse mode
  if(input_data[[1]]$data_mode == "sparse"){
    data_mode <- "sparse"
    print("Detected data type is SPARSE MODE. Generating plots for sparse mode results...")
  } else {
    data_mode <- "dense"
    print("Detected data type is DENSE MODE. Generating plots for dense mode results...")}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Initialize empty list 
  mats_list <- list()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Loop over the input_data
  for(i in seq_along(input_data)){
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Extract the reference dataframe from input_data
    #Subset to just replicate 1 (because we're using the 'perc_cell_death' or 'bliss_synergy' values, which are independent of replicates)
    ref_df <- 
      input_data[[i]]$ref_df %>%
      filter(replicate==1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Cap the drug names at the character limit
    d1_name_capped <- 
      input_data[[i]]$d1_name %>%
      substr(., 1, name_char_limit)
    
    d2_name_capped <- 
      input_data[[i]]$d2_name %>%
      substr(., 1, name_char_limit)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Round concentrations (dense mode format)
    ref_df <- 
      ref_df %>% 
      mutate(drug1_conc = round(drug1_conc, rounding_value),
             drug2_conc = round(drug2_conc, rounding_value))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #If the data is sparse mode, we need to reformat the matrix to resemble dense mode data
    #This keeps the plotting generalizable across either data mode
    if(data_mode=="sparse"){
      
      #First create a 0x0 row (drug1_conc==0 and drug2_conc==0)
      #This will get added to our modified ref_df, and will just have NA in all non-relevant columns
      template_row <- 
        ref_df[NA,] %>%
        slice(1) %>%
        mutate("drug1_name"=input_data[[i]]$d1_name,
               "drug2_name"=input_data[[i]]$d2_name,
               "drug1_conc"=as.numeric(0),
               "drug2_conc"=as.numeric(0),
               "perc_cell_death"=as.numeric(0))
      
      #Add the bliss_exp and bliss_synergy columns IF input_data already has synergy data
      if("bliss_exp" %in% colnames(ref_df)){
        template_row <- 
          template_row %>%
          mutate("bliss_exp"=as.numeric(0),
                 "bliss_synergy"=as.numeric(0))
      }
      
      ref_df <- 
        rbind(
          
          #Drug1-only data (10 rows)
          ref_df %>%
            filter(plate_type=="single_agent") %>%
            filter(drug1_name==input_data[[i]]$d1_name) %>%
            
            mutate("drug2_name"=input_data[[i]]$d2_name,
                   "drug2_conc"=0),
          
          #Drug2-only data (10 rows)
          ref_df %>%
            filter(plate_type=="single_agent") %>%
            filter(drug1_name==input_data[[i]]$d2_name) %>%
            mutate("drug2_name" = drug1_name,
                   "drug2_conc" = drug1_conc,
                   "drug1_name" = input_data[[i]]$d1_name,
                   "drug1_conc" = 0),
          
          #Combination data (100 rows)
          ref_df %>%
            filter(plate_type=="combination"),
          
          #0x0 row
          template_row)
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Plot the matrix
    p <- 
      ref_df %>%
      ggplot(aes(x=factor(as.numeric(drug1_conc)), y=factor(as.numeric(drug2_conc)), fill=.data[[plotting_variable]])) +
      geom_tile(color="black") +
      scale_fill_gradient2(
        low=  c("#6991a7ff", "#a1b5beff", "#bbc7c9ff", "#e8e6ddff"),
        mid=  c("#ffc4b0ff", "#ffa491ff"),
        high= c("#ff8f77ff", "#f5795eff", "#f25651ff", "#ee3d3cff"),
        midpoint = color_midpoint) +
      theme(axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(hjust=0.5)) +
      xlab(paste0(d1_name_capped, " (", unique(ref_df$units), ")")) +
      ylab(paste0(d2_name_capped, " (", unique(ref_df$units), ")")) +
      ggtitle(unique(ref_df$sample_name))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Change legend title based on plotting_variable
    if(plotting_variable == "perc_cell_death") {
      p <- p + labs(fill="%Cell Death")
    } else if(plotting_variable == "bliss_synergy") {
      p <- p + labs(fill="Bliss Synergy")
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Deposit plot into mats_list
    mats_list[[i]] <- p
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
  }#End loop over input_data elements
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Keep elements of the list same as input_data
  names(mats_list) <- names(input_data)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Return the list of plots:
  return(mats_list)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

}#End function
 