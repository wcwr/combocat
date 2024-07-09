#--------------------------------------------------------------
#Function to report a summary of all results
cc_report <- function(syn_data,
                      dr_data,
                      cd_plots,
                      syn_plots,
                      extra_plots,
                      save_summary_file = TRUE,
                      save_summary_plots = TRUE
){
  
  
  #Input requirements:
  #syn_data---------------------Output from cc_getSyn
  #dr_data----------------------Output from cc_getDr
  #cd_plots---------------------Output from cc_plotMat (with cell death as input)
  #syn_plots--------------------Output from cc_plotMat (with synergy as input)
  #extra_plots------------------Output from cc_plotExtras
  #save_summary_file------------Option to save the combined ref_df as a .csv
  #save_summary_plots-----------Option to save the summary plots as .svg files (NOTE! can take a long time)
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Initialize lists
  summary_plots_list <- list()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
   
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Combine all the ref_df's together
  combined_ref_df <- 
    syn_data %>%
    lapply(., function(x) x$ref_df) %>%
    bind_rows()
 
  
  #Save summary file if option is selected
  if(save_summary_file==TRUE){
    
    #Create the summary_file/ directory if it doesn't already exist
    if(!dir.exists(file.path( "summary_file"))){
      dir.create(file.path("summary_file"))
    }
    
    write.csv(combined_ref_df, file.path("summary_file", "combined_ref_df.csv"))
    
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Generate the combined observed vs. expected cell death scatter plot
  combined_cd_scatterplot <- 
    combined_ref_df  %>%
    filter(replicate==1) %>% #Only use the first replicate since perc_cell_death is the average of replicates
    ggplot(aes(x=bliss_exp, y=perc_cell_death, color=bliss_synergy)) +
    geom_point() +
    scale_color_gradient2("Bliss Synergy",
                          low=  c("#6991a7ff", "#a1b5beff", "#bbc7c9ff", "#e8e6ddff"),
                          mid=  c("#ffc4b0ff", "#ffa491ff"),
                          high= c("#ff8f77ff", "#f5795eff", "#f25651ff", "#ee3d3cff"),
                          midpoint=20) +
    theme_minimal() +
    xlab("%Cell Death [expected]") +
    ylab("%Cell Death [observed]") +
    geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed")
  
  ggsave("combined_cd_scatterplot.svg",
         combined_cd_scatterplot,
         width=5,
         height=3.5)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  

  if(syn_data[[1]]$data_mode=="dense"){
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle dense mode data------------------------
    #
    #=================================================================
    #=================================================================
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Get the unique filenames
    my_filenames <- 
      unique(combined_ref_df$filename) %>%
      .[!is.na(.)] #filter out any NA values
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Loop along the unique filenames
    for(i in my_filenames){
      
      
      #====================================
      #Subset dose-response data to any elements whose names contain i
      #The reason we do this is: 
      #e.g. if 'Plate1' contains DrugA and DrugB, and 'Plate2' contains DrugA and DrugC, then DrugA will have 2 DR curves
      #However, the DrugA in 'Plate1' is unique to that plate, so we use the plate/filename to subset the DR data
      subset_dr_list <- 
        dr_data$plots_list[grepl(i, names(dr_data$plots_list))]
      #====================================
      
      
      
      
      #====================================
      #Combine the DR curves into one plot:
      combined_curves <- subset_dr_list[[1]] / subset_dr_list[[2]]
      #====================================
      
      
      
      
      #====================================
      #Make the summary plot by combining several plots:
      summary_plot <- 
        wrap_plots(
          combined_curves, 
          cd_plots[[i]],
          syn_plots[[i]],
          extra_plots$syn_barplots_list[[i]],
          extra_plots$cd_scatterplots_list[[i]],
          nrow=1
        )
      #====================================
      
      
      
      
      #====================================
      #Deposit summary plot into list
      summary_plots_list[[i]] <- summary_plot
      #====================================
      
      
      
      
      #====================================
      #Save raw plate heatmaps if option is selected
      if(save_summary_plots==TRUE){
        
        #Create the raw_heatmaps/ directory if it doesn't already exist
        if(!dir.exists(file.path( "summary_plots"))){
          dir.create(file.path("summary_plots"))
        }
        
        print(paste0("Saving summary plot for plate ", i)) #Print message
        ggsave(paste0("summary_plots/", "summary_plots", i, ".svg"),
               summary_plot,
               width=21,
               height=4.6)
      }
      #====================================
      
      
    }#End loop along unique filenames

    
  } else {
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle sparse mode data-----------------------
    #
    #=================================================================
    #=================================================================
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Define a list that holds all unique combo_id's
    combo_id_list <- list()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Populate the list with unique combo_id's
    for(i in seq_along(syn_data)){
      combo_id_list[[i]] <- 
        syn_data[[i]]$ref_df %>%
        filter(plate_type=="combination") %>%
        pull(combo_id) %>%
        unique()
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Turn combo_id_list into a data frame and isolate drug1/2 names
    combo_id_df <- 
      do.call(rbind, combo_id_list) %>%
      as.data.frame() %>%
      rename("combo_id"=1) %>%
      #Set 'drug1_name' and 'drug2_name' as whatever comes before and after the underscore, respectively, keeping "combo_id" column":
      separate(combo_id, c("drug1_name", "drug2_name"), sep="_", remove=FALSE)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Loop along rows of combi_id_df
    for(i in 1:nrow(combo_id_df)){
      
      
      #====================================
      #Get the current drug names and combo_id
      d1_name <- combo_id_df$drug1_name[i]
      d2_name <- combo_id_df$drug2_name[i]
      
      current_combo_id <- combo_id_df$combo_id[i]
      #====================================
      
      
      
      
      #====================================
      #Combine the DR curves into one plot:
      combined_curves <- 
        dr_data$plots_list[[d1_name]] / dr_data$plots_list[[d2_name]]
      #====================================
      
      
      
      
      #====================================
      #Make the summary plot by combining several plots:
      summary_plot <- 
        wrap_plots(
          combined_curves, 
          cd_plots[[current_combo_id]],
          syn_plots[[current_combo_id]],
          extra_plots$syn_barplots_list[[current_combo_id]],
          extra_plots$cd_scatterplots_list[[current_combo_id]],
          nrow=1
        )
      #====================================
      
      
      
      
      #====================================
      #Deposit summary plot into list
      summary_plots_list[[i]] <- summary_plot
      #====================================
      
      
      
      
      #====================================
      #Save raw plate heatmaps if option is selected
      if(save_summary_plots==TRUE){
        
        #Create the raw_heatmaps/ directory if it doesn't already exist
        if(!dir.exists(file.path( "summary_plots"))){
          dir.create(file.path("summary_plots"))
        }
        
        print(paste0("Saving summary plot for plate ", i, ": ", current_combo_id)) #Print message
        ggsave(paste0("summary_plots/", "summary_plots", current_combo_id, ".svg"),
               summary_plot,
               width=21,
               height=4.6)
      }
      #====================================
      
    }#End loop along rows of combo_id_df
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
  }#End handling sparse mode data
  
  
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Return the combined_cd_scatterplot and summary_plots_list
  return(list(
    "combined_cd_scatterplot"=combined_cd_scatterplot,
    "summary_plots_list" = summary_plots_list
  ))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
}#End function
#--------------------------------------------------------------