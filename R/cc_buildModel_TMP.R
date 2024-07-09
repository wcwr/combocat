#--------------------------------------------------------------
#Function to build, tune, analyze, and save ML models using dense mode data
cc_buildModel <- function(summary_file,
                          min_mean_perc_cell_death=10,
                          training_prop=0.75,
                          tune_models=TRUE,
                          seed = NA,
                          save_models_file=TRUE,
                          color_midpoint=0.2,
                          save_cor_plots=FALSE,
                          save_importance_plots=FALSE,
                          animate_importance_plots=FALSE){
  
  
  #Input requirements:
  #summary_file--------------------Merged file containing ALL 'combined ref_df' files (outputs from cc_report)
  #min_mean_perc_cell_death--------Optional filtering of all data at the beginning, to only keep samples with mean %CD above this value
  #training_prop-------------------Proportion of data to use for training (default=0.75 aka 75%)
  #tune_models---------------------Option to use hyperparameter tuning for each of the 90 models
  #................................Note: This can increase the time to build the models by 10-20x!
  #seed----------------------------Integer for reproducibility when testing outputs (Don't use this if not doing tests!)
  #................................If using a seed, replace NA with an integer
  #save_models_file----------------Option to save the list of (bundled/serialized) 90 models as an .RDS file
  #color_midpoint------------------The midpoint for the color scale in the importance plots (default=0.2)
  #save_cor_plots------------------Option to save %cd correlation plots for each model (note: produces 90 .svg files!)
  #save_importance_plots-----------Option to save feature importance barplots and heatmaps for each model (produces 180 total .svg files!)
  #animate_importance_plots--------Option to produce & save animation of the importance barplots/heatmaps across all models (produces 3 .gif files)
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Set the seed if doing tests (like directly comparing results from different hyperparameters)
  if(!is.na(seed)){
    set.seed(seed)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Print message if hyperparameter tuning is enabled
  if(tune_models==TRUE){
    print("Hyperparameter tuning is ENABLED")
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create/set the working directory as a folder named with the current date and time:
  current_date_time <- format(Sys.time(), "%d_%b_%Y_%I.%M%p")
  
  base_dir_name <- paste0("Model_outputs_", current_date_time)
  
  dir.create(base_dir_name)
  setwd(base_dir_name)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create directories if certain plot saving options are selected
  
  #For saving %cd cor plots
  if(save_cor_plots==TRUE){
    
    #Make a directory to save correlation plots
    cor_plots_dir <- paste0(getwd(), "/cor_plots")
    
    dir.create(cor_plots_dir, showWarnings = FALSE)
  }
  
  
  #For saving importance_plots
  if(save_importance_plots==TRUE){
    
    #Make directories to save importance plots
    importance_barplots_dir <- paste0(getwd(), "/importance_barplots")
    importance_heatmaps_dir <- paste0(getwd(), "/importance_heatmaps")
    
    dir.create(importance_barplots_dir, showWarnings = FALSE)
    dir.create(importance_heatmaps_dir, showWarnings = FALSE)
  }
  
  
  #For saving summary plots:
  summary_plots_dir <- paste0(getwd(), "/summary_plots")
  dir.create(summary_plots_dir, showWarnings = FALSE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Add a mean %CD column, and filter out samples with mean %CD below min_mean_perc_cell_death
  summary_file <- 
    summary_file %>%
    group_by(filename, combo_id) %>%
    mutate("mean_perc_cell_death" = mean(perc_cell_death)) %>%
    ungroup()
  
  #Get any samples which need to be dropped
  dropped_samples <- 
    summary_file %>%
    select(filename, combo_id, sample_name, mean_perc_cell_death) %>%
    filter(mean_perc_cell_death <= min_mean_perc_cell_death) %>%
    unique()
  
  #Filter out dropped samples
  summary_file <- 
    summary_file %>%
    filter(mean_perc_cell_death >= min_mean_perc_cell_death)
  
  #Print message and write file alerting user to dropped samples
  if(!is.na(min_mean_perc_cell_death)){
    print(paste0("Dropping ", nrow(dropped_samples), " samples with mean %CD below ", min_mean_perc_cell_death, "%"))
    write.csv(dropped_samples, "dropped_samples.csv", row.names=FALSE)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Explicitly define all necessary indices of a total 11x11 matrix (a 10x10 matrix + the single-agent wells)
  
  #Define the indices corresponding to the 10x10 matrix diagonal (n=10) 
  #NOTE: This is lower-left to upper-right diagonal - different from "diag" function in R which is opposite
  #NOTE: we DO NOT include the 0x0 square as an input feature. A model will not be created for the 0x0 square. 
  diag_squares <- data.frame("row"=c(10,9,8,7,6,5,4,3,2,1),
                             "col"=c(2,3,4,5,6,7,8,9,10,11),
                             "label"=c("diag"))
  
  #Define the indices corresponding to drug1only squares (column-wise) (n=10)
  drug1only_squares <- data.frame("row"=c(11,11,11,11,11,11,11,11,11,11),
                                  "col"=c(2,3,4,5,6,7,8,9,10,11),
                                  "label"=c("drug1only"))
  
  #Define the indices corresponding to drug2only squares (row-wise) (n=10)
  drug2only_squares <- data.frame("row"=c(1,2,3,4,5,6,7,8,9,10),
                                  "col"=c(1,1,1,1,1,1,1,1,1,1),
                                  "label"=c("drug2only"))
  
  
  #Define the indices corresponding to the squares we're trying to predict (n=90)
  #Note: We try to predict these one at a time (producing 90 total models, one for each of these squares)???
  predictor_squares <- data.frame("row"=c(1,2,3,4,5,6,7,8,9, 
                                          1,2,3,4,5,6,7,8,10,
                                          1,2,3,4,5,6,7,9,10,
                                          1,2,3,4,5,6,8,9,10,
                                          1,2,3,4,5,7,8,9,10,
                                          1,2,3,4,6,7,8,9,10,
                                          1,2,3,5,6,7,8,9,10,
                                          1,2,4,5,6,7,8,9,10,
                                          1,3,4,5,6,7,8,9,10,
                                          2,3,4,5,6,7,8,9,10),
                                  
                                  "col"=c(2,2,2,2,2,2,2,2,2,
                                          3,3,3,3,3,3,3,3,3,
                                          4,4,4,4,4,4,4,4,4,
                                          5,5,5,5,5,5,5,5,5,
                                          6,6,6,6,6,6,6,6,6,
                                          7,7,7,7,7,7,7,7,7,
                                          8,8,8,8,8,8,8,8,8,
                                          9,9,9,9,9,9,9,9,9,
                                          10,10,10,10,10,10,10,10,10,
                                          11,11,11,11,11,11,11,11,11),
                                  
                                  "label"=c("predictor"))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Initialize lists that contain outputs producing during the following loop stage
  importance_list         <- list() #For feature importance
  model_metrics_list      <- list() #For the metrics of the final model (hyperparameters)() #For final (bundled) models
  prediction_list         <- list() #For prediction results (obs vs. pred %Cell death)
  models_list             <- list() #For final (bundled) models
  cor_plots_list          <- list() #For correlation plots
  importance_barplot_list <- list() #For importance barplots
  importance_heatmap_list <- list() #For importance heatmaps
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Begin prediction square-wise loop
  for(i in 1:nrow(predictor_squares)){
  
    
    #====================================
    #Print the progress message
    percent_complete <- round((i / nrow(predictor_squares)) * 100, 2)
    
    # Print the progress message
    cat(sprintf("Building Model %d of %d.....%s %% complete\n", 
                i, nrow(predictor_squares), percent_complete))
    #====================================
    
    
    
    #====================================
    #Subset the current square for which we're building a model
    current_square <- predictor_squares[i,]
    
    
    #Combine the data frames, into a single df 
    #This will contain the 30 features + the single current prediction square (total n=31 rows)
    input_squares <- rbind(current_square,
                           diag_squares,
                           drug1only_squares,
                           drug2only_squares)
    
    #Add a column to the input that gets populated by all the corresponding values
    input_squares$value <- NA
    #====================================
    
    
    
    
    #====================================
    #Populate the 'value' column of input_squares
    #This will get the 31 input squares the model needs for each 10x10 matrix
    
    #Initialize list for all the input data to go into
    input_list <- list() 
    
    #Make a dataframe that contains unique filename-combo_id-sample_name combinations
    #We will use this to make sure we're looping over unique samples
    #Each of these must be unique because they'll all be row names in a resulting data frame
    #NOTE! While each dense mode experiment itself requires unique filenames...
    #......the filenames when combining ALL dense mode summary files may not be unique. E.g., DrugA_DrugB.csv (in CHP134) vs. DrugA_DrugB.csv (BE2C)
    #......that's why we must consider the filename and sample_name (and combo_id)
    samples_df <- 
      summary_file %>%
      select(filename, combo_id, sample_name) %>%
      mutate("concat"=paste0(filename,"_",combo_id,"_",sample_name)) %>% 
      unique() %>%
      filter(!is.na(filename))
    
    
    #Loop over each unique filename-combo_id-sample_name combination
    for(j in unique(samples_df$concat)){
      

      #____________________________________
      #Get the current filename, combo_id, and sample_name
      current_filename <- 
        samples_df %>%
        filter(concat==j) %>%
        pull(filename)
      
      current_combo_id <- 
        samples_df %>%
        filter(concat==j) %>%
        pull(combo_id)
      
      current_sample_name <-
        samples_df %>%
        filter(concat==j) %>%
        pull(sample_name)
      #____________________________________
      
      
      
      
      #____________________________________
      #Extract the current matrix
      #We 'dcast' the current tabular-form data into the 10x10 (11x11) matrix
      #Note: we use 'perc_cell_death' which mean the input is the MEAN matrix of the 3 replicates (per combo)
      current_mat <- 
        summary_file %>%
        filter(filename   == current_filename,
               combo_id   == current_combo_id,
               sample_name== current_sample_name) %>%
        select(drug1_conc, drug2_conc, perc_cell_death) %>%
        unique() %>%
        #'dcast' into matrix form
        #CRITICAL: We do 'drug2_conc ~ drug1_conc' so that drug1 is the x-axis (columns) and drug2 is the y-axis (rows)
        dcast(drug2_conc ~ drug1_conc, value.var="perc_cell_death") %>%
        .[rev(seq_len(nrow(.))), ] %>% #Flip matrix top-to-bottom so that 0x0 is at the bottom-left
        #Set the first column as the row names
        'rownames<-'(., .[,1]) %>% .[,-1]
      #____________________________________
      
      
      
      
      #____________________________________
      #Define a variable 'current_input' that will be the input_squares before populating the 'value' column
      current_input <- input_squares
      #____________________________________
      
      
      
      
      #____________________________________
      #Loop along the input squares to populate the current_input dataframe
      for(k in 1:nrow(input_squares)){
        
        #Populate the 'value' column of current_input with the corresponding values
        current_input[k,]$value <- current_mat[current_input[k,]$row, current_input[k,]$col]
        
      }
      #____________________________________
      
      
      
      
      #____________________________________
      #Combine the row,col and label data into a column called 'index'
      #The string "ind_1.4_diag" corresponds to square [1,4] which is in the diagonal group
      current_input$index <- paste0("ind_",current_input$row,".",current_input$col,"_",current_input$label)
      #____________________________________
      
      
      
      
      #____________________________________
      #Get the current_input into a format compatible with the ML modeling:
      #There should be 31 columns (corresponding to the 30 features + 1 outcome variable/predictor square)
      #There should be 1 row, which should be named as the 'concat' name (filename+combo_id+sample_name)
      #The 'predictor' column should be the first one - it's what we're trying to predict
      current_input <- 
        current_input %>%
        select(index, value) %>%
        t() %>%
        as.data.frame() %>%
        'colnames<-'(.[1,]) %>% #Set colnames to the first row
        filter(rownames(.)=="value") %>%
        'rownames<-'(., j)      #Set rownames to the current 'concat' name
      #____________________________________
      
      
      
      
      #____________________________________
      #Deposit the dataframe into input_list
      input_list[[j]] <- current_input
      #____________________________________
      
      
    }#End loop over unique samples
    
    
    
    
    #====================================
    #Combine all dataframes of input_list into a single dataframe
    #There should still be just 31 columns, corresponding to the 31 input squares
    #The first column should be the 'predictor' square, aka what we're trying to predict
    #Each row is a different sample
    input_df <- do.call(rbind, input_list)
    
    #Also need to make sure all columns are numeric:
    input_df <- input_df %>% mutate_all(as.numeric)
    #====================================
    
    
    
    
    #====================================
    #Define the square name
    #Only one column at this point has the word "predictor" in it, which should be the 1st column. That's the current square name
    current_square_name <- colnames(input_df)[grepl("predictor", colnames(input_df))]
    
    
    #Also define the outcome variable index more formally in the format of "[1,2]"  [row,col]
    out_var <- paste0("[", predictor_squares[i,]$row, ",", predictor_squares[i,]$col, "]")
    #====================================
    
    
    
    
    #====================================
    #Randomly shuffle the data rows before splitting into train and test
    input_df <- 
      input_df %>%
      sample_frac(1)
    #====================================
    
    
    
    
    #====================================
    #Split the data into train and test
    input_split <- initial_split(input_df, prop = training_prop)
    
    training_set <- training(input_split)
    test_set     <- testing(input_split)
    #====================================
    
    
    
    
    #====================================
    if(tune_models==TRUE){
      
      #Hyperparameter tuning
      
      #____________________________________
      #Set up the model specification
      #The hyperparameters will be tuned
      xgb_spec <- 
        boost_tree(
        trees = 1000,
        tree_depth = tune(),
        min_n = tune(),
        loss_reduction = tune(),                    
        sample_size = tune(),
        mtry = tune(),   
        learn_rate = tune()) %>%
        set_engine("xgboost") %>%
        set_mode("regression")
      #____________________________________
      
      
      
      
      #____________________________________
      #Set up a space-filling grid design to cover the hyperparameter space as well as possible
      xgb_grid <-
        grid_latin_hypercube(
        tree_depth(),
        min_n(),
        loss_reduction(),
        sample_size = sample_prop(),
        finalize(mtry(), training_set), #gets treated differently b/c it depends on actual # of predictors in data
        learn_rate(),
        size = 40)
      #____________________________________
      
      
      
      
      #____________________________________
      #Put the model specification into a workflow:
      
      #Try to predict current_square_name, using all other columns whose name is not current_square_name
      my_char <- paste(current_square_name, "~ .")
      my_formula <- as.formula(my_char)
      
      xgb_wf <- 
        workflow() %>%
        add_formula(my_formula) %>% #Equivalent to add_formula(ind_10.11_predictor ~.)
        add_model(xgb_spec)
      #____________________________________
      
      
      
      
      #____________________________________
      #Create cross-validation resamples
      input_folds <- vfold_cv(training_set, strata = current_square_name)
      #____________________________________
      
      
      
      
      #____________________________________
      #Prepare for parallelization before tuning:
      num_cores <- parallelly::availableCores()      #Detect # of cores
      cl <- parallel::makeCluster(num_cores-2)       #Use the total # of available cores MINUS 2 (to be free for other processes)
                                                     #...If a crash occurs, the last log file will correspond to the iterand that crashed
      doParallel::registerDoParallel(cl)
      #____________________________________
      
      
      
      
      #____________________________________
      #Use the tuneable workflow to tune
      #This step can be time consuming!
      xgb_res <- 
        tune_grid(
        xgb_wf,
        resamples = input_folds,
        grid = xgb_grid,
        control = control_grid(save_pred = TRUE))
      #____________________________________
      
      
      
      
      #____________________________________
      #Stop parallelization
      stopCluster(cl)
      #____________________________________
      
      
      #____________________________________
      #Collect model metrics
      xgb_res_metrics <- 
        xgb_res %>%
        collect_metrics() %>%
        mutate("current_square_name"=current_square_name,
               "out_var"=out_var)
        
      #Deposit into model_metrics_list
      model_metrics_list[[current_square_name]] <- xgb_res_metrics
      #____________________________________
      
      
      
      
      #____________________________________
      #Select the best parameters based on RMSE
      #NOTE: The best parameters can also be selected based on R^2 
      #Fitting models based on best parameters chosen by R^2 or RMSE produce slightly different results
      best_rmse <- select_best(xgb_res, "rmse")
      #____________________________________
      
      
      
      
      #____________________________________
      #Finalize the tuneable workflow using the best parameters
      final_xgb <- 
        finalize_workflow(
          xgb_wf,
          best_rmse
        )
      #____________________________________

      
     
      
    } else {
      #Else if NOT tuning hyperparameters:
      
      
      
      #____________________________________
      #Set up the model specification
      xgb_spec <- 
        boost_tree(trees = 1000) %>%
        set_mode("regression") %>%
        set_engine("xgboost")
      #____________________________________
      
      
      
      
      #____________________________________
      #Put the model specification into a workflow:
      
      #Try to predict current_square_name, using all other columns whose name is not current_square_name
      my_char <- paste(current_square_name, "~ .")
      my_formula <- as.formula(my_char)
      
      xgb_wf <- 
        workflow() %>%
        add_formula(my_formula) %>% #Equivalent to add_formula(ind_10.11_predictor ~.)
        add_model(xgb_spec)
      #____________________________________
      
      
      
      
      #____________________________________
      #Create cross-validation resamples
      input_folds <- vfold_cv(training_set, strata = current_square_name)
      #____________________________________
      
      
      
      
      #____________________________________
      #Fit and evaluate the model on each fold
      #Note we do it directly on the folds without tuning since we are not tuning hyperparameters here
      cv_results <- 
        fit_resamples(
          xgb_wf,
          resamples = input_folds
        )
      #____________________________________
      
      
      

      #____________________________________
      #Collect model metrics
      #This is the same as the tuning case, but we are collecting metrics from the cv_results
      xgb_res_metrics <- 
        cv_results %>%
        collect_metrics() %>%
        mutate("current_square_name"=current_square_name,
               "out_var"=out_var)
      
      #Deposit into model_metrics_list
      model_metrics_list[[current_square_name]] <- xgb_res_metrics
      #____________________________________
      
      
      
      
      #____________________________________
      #Set the final_xgb to simply the xgb_wf
      #Note that in the case of tuning, this final_xgb is the final tuned model
      final_xgb <- xgb_wf
      #____________________________________

      
    }#End of model generation WITHOUT tuning
    
    
    
    #====================================
    #Extract the variable importance of all 30 features for the current model
    importance_df <- 
      final_xgb %>%
      fit(data = training_set) %>%
      extract_fit_parsnip() %>%
      vip::vi() %>% 
      
      #Transform e.g. 'ind_11.10_drug1only' into columns 'row' (11), 'col' (10), and 'label' (drug1only)
      mutate(
        # Extract the row number by finding digits before the period
        "row" = as.numeric(str_extract(Variable, "(?<=ind_)\\d+")),
        # Extract the col number by finding digits after the period and before the underscore
        "col" = as.numeric(str_extract(Variable, "(?<=\\.)\\d+")),
        # Extract the label by finding everything after the last underscore
        "label" = str_extract(Variable, "(?<=_)[^_]+$")
      ) %>%
      
      mutate("index"=paste0("[",row,",",col,"]")) %>% #add the index column where each row should be formatted as [row,col]
      mutate("out_var"=out_var) %>% #Add the current outcome variable
      mutate("current_square_name"=current_square_name)  #Add the current_square_name
    #====================================
    
    
    
    
    #====================================
    #Deposit the importance_df into the importance_list
    importance_list[[current_square_name]] <- importance_df
    #====================================
    
    
    
    
    #====================================
    #Fit the final best model to training set and evaluate the test set
    final_res <- last_fit(final_xgb, input_split)
    #====================================
    
    
    
    
    #====================================
    #Save the final model (bundled/serialized)
    res_bundle <- 
      final_res %>%
      extract_workflow() %>%
      butcher() %>%
      bundle()
    
    #Deposit into list
    models_list[[current_square_name]] <- res_bundle
    #====================================
    
    
    
    
    #====================================
    #Get the model-predicted values of the test set
    pred_df <- 
      final_res %>%
      collect_predictions() %>%
      as.data.frame() %>%
      mutate("sample_name"=rownames(test_set), .before=1) %>%
      rename("obs_cd" =current_square_name,
             "pred_cd"=".pred") %>%
      mutate("current_square_name"=current_square_name,
             "out_var"=out_var)
    
 
      #Add pred_df to the prediction_list
      prediction_list[[current_square_name]] <- pred_df
    #====================================
    
    
    
    
    #====================================
    #Create a scatterplot of the observed vs predicted values
    
    #Calculate the correlation between the observed and predicted values
    mycor <- cor.test(pred_df$obs_cd, pred_df$pred_cd, method="pearson")
    
    #Create the correlation scatterplot of obs vs. pred %Cell Death
    cd_cor_plot <- 
      pred_df %>%
      ggplot(aes(x=obs_cd, y=pred_cd)) +
      geom_point() +
      geom_smooth(method="lm", color="forestgreen", se=FALSE)  +
      xlab("%Cell Death (observed)") +
      ylab("%Cell Death (predicted)") +
      theme_minimal() +
      ggtitle(paste0("Model ", 
                     out_var, 
                     " r=",round(mycor$estimate,2),
                     " (n=", nrow(pred_df), ")")) +
      theme(plot.title = element_text(hjust=0.5))
    
    
    #Deposit into the cor_plot_list
    cor_plots_list[[current_square_name]] <- cd_cor_plot
    #====================================
    
    
    
    
    #====================================
    #If save_cor_plots is TRUE, save the %cd correlation plots:
    if(save_cor_plots==TRUE){
      
      ggsave(
        filename = paste0(cor_plots_dir, "/", out_var, ".svg"),
        plot = cd_cor_plot,
        width = 4.5,
        height = 3)
    }
    #====================================
    
    
    
    
    #====================================
    #Create a barplot of the variable importance
    importance_barplot <- 
      importance_df %>%
      ggplot(aes(x=Importance, y=reorder(index, Importance), fill=Importance)) +
      geom_col(color="black") +
      scale_fill_gradient2(
        low=  c("#6991a7ff", "#a1b5beff", "#bbc7c9ff", "#e8e6ddff"),
        mid=  c("#ffc4b0ff", "#ffa491ff"),
        high= c("#ff8f77ff", "#f5795eff", "#f25651ff", "#ee3d3cff"),
        midpoint = 0.2) +
      theme_minimal() +
      xlab("Importance") +
      ylab("Feature Index [row,column]") +
      ggtitle(paste0("Model ", out_var)) +
      theme(plot.title = element_text(hjust=0.5))
    #====================================
    
    
    

    #====================================
    #Create a matrix heatmap of the variable importance

    #Extract the row and col from current_square_name:
    current_row <- as.numeric(str_extract(current_square_name, "(?<=ind_)\\d+"))
    current_col <- as.numeric(str_extract(current_square_name, "(?<=\\.)\\d+"))
    
    #First, create a reference dataframe with all possible combinations of row and column
    ref_df <- expand_grid("row"=seq(1:11), 
                          "col"=seq(1:11))
    
    #Then, merge the reference dataframe with the importance_df
    ref_df <- 
      ref_df %>%
      left_join(importance_df, by=c("row"="row", "col"="col")) %>%
      mutate(Importance=ifelse(is.na(Importance), NA, Importance))
    
    #Finally, create the heatmap
    importance_heatmap <- 
      ref_df %>%
      ggplot(aes(x=factor(col), y=factor(-row), fill=Importance)) +
      geom_tile(color="black") +
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0), labels=rev(seq(1:11))) +
      scale_fill_gradient(
        low=  c("#6991a7ff", "#a1b5beff", "#bbc7c9ff"),
        high= c("#ff8f77ff", "#f5795eff", "#f25651ff", "#ee3d3cff"),
        na.value = "white") +
      geom_point(data=ref_df %>% filter(row==current_row, col==current_col),
                 fill="#ebe45b",
                 color="black",
                 shape=21,
                 size=6) +
      xlab("Feature Column") +
      ylab("Feature Row") +
      ggtitle(paste0("VI for predicting ", out_var)) +
      theme(plot.title = element_text(hjust=0.5))
    #====================================
    
    
    
    
    #====================================
    #If save_importance_plots is TRUE, save the barplots/heatmaps in respective folders:
    if(save_importance_plots==TRUE){
      
      ggsave(
        filename = paste0(importance_barplots_dir, "/", out_var, ".svg"),
        plot = importance_barplot,
        width = 4.5,
        height = 5)
      
      ggsave(
        filename = paste0(importance_heatmaps_dir, "/", out_var, ".svg"),
        plot = importance_heatmap,
        width = 6,
        height = 5)
    }
    #====================================
    
    
  }#End loop over predictor squares
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Save the (bundled/serialized) models as a .RDS file (if option is selected)
  if(save_models_file==TRUE){
    
    saveRDS(models_list, 
            paste0("m", #Start it with a letter to avoid opening the file later in RStudio which doesn't like files that start with numbers
                   current_date_time,
                   "_Train",
                   nrow(training_set),
                   "_Test",
                   nrow(test_set),
                   "_Tune",
                   if(tune_models==TRUE){"T"}else{"F"},
            ".RDS")
    )
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Combine various lists into single data frames for summarizing, and save them
  combined_prediction_df    <- do.call(rbind, prediction_list)
  combined_model_metrics_df <- do.call(rbind, model_metrics_list)
  combined_importance_df    <- do.call(rbind, importance_list) %>% as.data.frame()
  
  
  write.csv(combined_prediction_df, "combined_prediction_df.csv", row.names=FALSE)
  write.csv(combined_model_metrics_df, "combined_model_metrics_df.csv", row.names=FALSE)
  write.csv(combined_importance_df, "combined_importance_df.csv", row.names=FALSE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Modify current_model_metrics_df to have 'row', 'col', 'label', for plotting:
  #Split the 'current_square_name' column into columns row/col/label
  combined_model_metrics_df <- 
    combined_model_metrics_df %>%
    mutate(
      # Extract the row number by finding digits before the period
      "row" = as.numeric(str_extract(current_square_name, "(?<=ind_)\\d+")),
      # Extract the col number by finding digits after the period and before the underscore
      "col" = as.numeric(str_extract(current_square_name, "(?<=\\.)\\d+")),
      # Extract the label by finding everything after the last underscore
      "label" = str_extract(current_square_name, "(?<=_)[^_]+$")
    )
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create a summary scatterplot of observed vs. predicted %cd values
  
  #Calculate the correlation between observed and predicted %cd values
  mycor_combined <- cor.test(combined_prediction_df$obs_cd, combined_prediction_df$pred_cd)
  
  combined_cd_cor_plot <- 
    combined_prediction_df %>%
    ggplot(aes(x=obs_cd, y=pred_cd)) +
    geom_point() +
    geom_smooth(method="lm", color="forestgreen", se=FALSE) +
    xlab("%Cell Death (observed)") +
    ylab("%Cell Death (predicted)") +
    theme_minimal() +
    ggtitle(paste0(" r=", round(mycor_combined$estimate,2),
                   " (",
                   
                   #Include the number of samples and the number of models
                   #number of models is always 90
                   #number of samples = number of rows in combined_prediction_df/90
                   #The number of points in the plot is 90*number of samples
                   #Ex: 15 samples in the test set will produce 90*15=1350 data points
                   
                   nrow(combined_prediction_df)/90, " samples, 90 models)")) +
    theme(plot.title = element_text(hjust=0.5))
  
  
  ggsave(filename=paste0(summary_plots_dir, "/", "summary_cd_cor_plot.svg"),
         plot=combined_cd_cor_plot,
         width=4.5,
         height=3)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create a summary heatmap plot of the mean variable importance
  #This can show which features are most important across all models (aka most important squares overall)
  mean_importance_heatmap <- 
    combined_importance_df %>%
    group_by(Variable) %>%
    mutate("mean_importance"=mean(Importance)) %>%
    ungroup() %>%
    merge(.,
          expand_grid("row"=seq(1:11),
                      "col"=seq(1:11)),
          by=c("row","col"),
          all.y=TRUE) %>%
    ggplot(aes(x=factor(col), y=factor(-row), fill=mean_importance)) +
    geom_tile(color="black") +
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0), labels=rev(seq(1:11))) +
    scale_fill_gradient(
      "Mean importance",
      low=  c("#6991a7ff", "#a1b5beff", "#bbc7c9ff"),
      high= c("#ff8f77ff", "#f5795eff", "#f25651ff", "#ee3d3cff"),
      na.value = "white") +
    xlab("Feature Column") +
    ylab("Feature Row") +
    ggtitle("Mean VI across all models") +
    theme(plot.title = element_text(hjust=0.5))
  
  
  ggsave(filename=paste0(summary_plots_dir, "/", "mean_importance_heatmap.svg"),
         plot=mean_importance_heatmap,
         width=6,
         height=5)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create a summary heatmap of the mean RMSE
  #This is the mean RMSE across all folds, for each model
  mean_rmse_heatmap <- 
    combined_model_metrics_df %>%
    filter(.metric=="rmse") %>%
    filter(!is.na(.metric)) %>%
    merge(.,
          expand_grid("row"=seq(1:11),
                      "col"=seq(1:11)),
          by=c("row","col"),
          all.y=TRUE) %>%
    ggplot(aes(x=factor(col), y=factor(-row), fill=mean)) +
    geom_tile(color="black") +
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0), labels=rev(seq(1:11))) +
    scale_fill_gradientn(
      "RMSE",
      colors=c("#ee3d3cff", "#f25651ff", "#f5795eff", "#ff8f77ff", "#e8e6ddff", "#bbc7c9ff", "#a1b5beff", "#6991a7ff"),
      na.value="white"
    ) +
    xlab("Feature Column") +
    ylab("Feature Row") +
    ggtitle("Model-wise RMSE (mean across CV folds)") +
    theme(plot.title = element_text(hjust=0.5))
  
   
  
  ggsave(filename=paste0(summary_plots_dir, "/", "mean_rmse_heatmap.svg"),
         plot=mean_rmse_heatmap,
         width=6,
         height=5)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Create a summary heatmap of the mean r^2
  #This is the mean r^2 across all folds, for each model
  mean_rsq_heatmap <- 
    combined_model_metrics_df %>%
    filter(.metric=="rsq") %>%
    filter(!is.na(.metric)) %>%
    merge(.,
          expand_grid("row"=seq(1:11),
                      "col"=seq(1:11)),
          by=c("row","col"),
          all.y=TRUE) %>%
    ggplot(aes(x=factor(col), y=factor(-row), fill=mean)) +
    geom_tile(color="black") +
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0), labels=rev(seq(1:11))) +
    scale_fill_gradientn(
      "r^2",
      colors=c("#ee3d3cff", "#f25651ff", "#f5795eff", "#ff8f77ff", "#e8e6ddff", "#bbc7c9ff", "#a1b5beff", "#6991a7ff"),
      na.value="white"
    ) +
    xlab("Feature Column") +
    ylab("Feature Row") +
    ggtitle("Model-wise r^2 (mean across CV folds)") +
    theme(plot.title = element_text(hjust=0.5))
  
  
  ggsave(filename=paste0(summary_plots_dir, "/", "mean_rsq_heatmap.svg"),
         plot=mean_rsq_heatmap,
         width=6,
         height=5)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Animate the importance plots if selected
  if(animate_importance_plots==TRUE){
    
    print("Animating variable importance")
    
    
    #====================================
    #Parse out the row and col correspond to out_var
    combined_importance_df <- 
      combined_importance_df %>%
      mutate("out_var_split" = str_extract_all(current_square_name, "\\d+"), # Extract all numbers
             "out_var_row" = as.numeric(map_chr(out_var_split, 1)), # First number is row
             "out_var_col" = as.numeric(map_chr(out_var_split, 2)))  # Second number is col

    
    #Arrange by row and then by column to get the desired order
    #This ensures the animation proceeds in a more logical order
    combined_importance_df <-
      combined_importance_df %>%
      arrange(out_var_row, out_var_col) %>%
      mutate("out_var_key"=paste0("[", out_var_row, ",", out_var_col, "]"))
    
    
    #Create a unique key based on the row and column, which will be used for animation ordering
    combined_importance_df$out_var_key <- 
      factor(combined_importance_df$out_var_key, levels = unique(combined_importance_df$out_var_key))
    #====================================
    
    
    
    
    #====================================
    #Create the desired order which will be the y-axis of the bar plot
    #We do this to make the order of indices fairly logical (starting with [1,1] and ending with [11,11])
    #Note, this only has to do with the y-axis ordering, nothing to do with ordering of the animation
    desired_order <- 
      expand.grid(row = 1:11, col = 1:11) %>%
      mutate("index" = paste0("[", row, ",", col, "]")) %>%
      pull(index)
    
    #Animate the bar plot:
    animated_importance_barplot <- 
      combined_importance_df %>%
      mutate("index"=factor(index, levels = rev(desired_order))) %>%
      ggplot(aes(x=Importance,  y=index, fill=Importance)) +
      geom_col(color="black") +
      scale_fill_gradient(
        low=  c("#6991a7ff", "#a1b5beff", "#bbc7c9ff"),
        high= c("#ff8f77ff", "#f5795eff", "#f25651ff", "#ee3d3cff"),
        na.value = "white") +
      theme_minimal() +
      xlab("Importance") +
      ylab("Feature Index [row,column]") +
      theme(plot.title = element_text(hjust=0.5)) +
      ggtitle(paste0("VI for predicting ", "{closest_state}")) +
      transition_states(out_var_key) +
      ease_aes('linear')
    #====================================
    
    
    
    
    #====================================
    #Begin animating importance heatmap
    
    
    
    #Create the reference dataframe with all possible combinations of rows and columns
    #This is needed to ensure the heatmap has the empty squares when a feature is not present
    ref_df <- expand_grid("row"=seq(1:11), 
                          "col"=seq(1:11))
    
    #Initialize an empty dataframe that will be used for plotting
    current_df_final <- data.frame()
    #====================================
    
    
    
    
    #Loop through each square to be predict (aka outcome variable) and modify dataframe as needed:
    for(i in unique(combined_importance_df$current_square_name)){
      
      #____________________________________
      #Extract the current square (this is what is being predictor)
      the_current_square <- i
      
      #Extract the row and column that the current square belongs to 
      the_current_row    <- str_extract(i, "(?<=ind_)\\d+") %>% as.numeric()
      the_current_col    <- str_extract(i, "(?<=\\.)\\d+")  %>% as.numeric()
      #____________________________________
      
      
      
      
      #____________________________________
      #Subset the importance dataframe to only include the current square
      current_df <- 
        combined_importance_df %>%
        filter(current_square_name == the_current_square)
      #____________________________________
      
      
      
      
      #____________________________________
      #Extract out the current outcome variable as well 
      current_out_var <- 
        current_df %>%
        filter(!is.na(out_var)) %>%
        pull(out_var) %>%
        unique()
      #____________________________________
        
      
      
      
      #____________________________________
      #Merge the current dataframe with the reference dataframe
      current_df <- 
        ref_df %>%
        merge(., 
              current_df,
              by=c("row", "col"), 
              all.x=TRUE)
      #____________________________________
      
      
      
      
      #____________________________________
      #Now fix the NA gaps introduced by the merge in 'current_square_name' and 'out_var' columns
      current_df$current_square_name <- i
      current_df$out_var <- current_out_var
      
      #Add a column denoting whether or not the feature square indices ('row', 'col') correspond to the outcome variable square indices:
      current_df$is_current <- current_df$row == the_current_row & current_df$col == the_current_col
      
      #Bind the current dataframe to the final dataframe
      current_df_final <- rbind(current_df_final, current_df)
      #____________________________________
      
      }#End loop over 'current_square_name' / outcome variables
    
    
    #====================================
    #Parse out the row and col correspond to out_var
    current_df_final <- 
      current_df_final %>%
      mutate(
        out_var_split = str_extract_all(out_var, "\\d+"), # Extract all numbers
        out_var_row = as.numeric(map_chr(out_var_split, 1)), # First number is row
        out_var_col = as.numeric(map_chr(out_var_split, 2)))  # Second number is col
    
    #Arrange by row and then by column to get the desired order
    #This ensures the animation proceeds in a more logical order
    current_df_final <- current_df_final %>%
      arrange(out_var_row, out_var_col)
    
    #Create a unique key based on the row and column, which will be used for animation ordering
    current_df_final <- current_df_final %>%
      mutate("out_var_key" = paste0("[", out_var_row, ",", out_var_col, "]"))
    
    # Make sure that 'out_var_key' is a factor with the levels in the order you want
    current_df_final$out_var_key <- factor(current_df_final$out_var_key, levels = unique(current_df_final$out_var_key))
    #====================================
    
    
    
    
    #====================================
    #Animate the heatmap
    animated_importance_heatmap <- 
      current_df_final %>%
      ggplot(aes(x=factor(col), y=factor(-row), fill=Importance)) +
      geom_tile(color="black") +
      scale_x_discrete(expand=c(0,0))+
      scale_y_discrete(expand=c(0,0), labels=rev(seq(1:11))) +
      scale_fill_gradient(
        low=  c("#6991a7ff", "#a1b5beff", "#bbc7c9ff"),
        high= c("#ff8f77ff", "#f5795eff", "#f25651ff", "#ee3d3cff"),
        na.value = "white") +
      theme_minimal() +
      theme(plot.title = element_text(hjust=0.5)) +
      xlab("Feature Column") +
      ylab("Feature Row") +
      geom_point(data = 
                  current_df_final %>%
                  filter(is_current==TRUE),
                fill="#ebe45b",
                color="black",
                shape=21,
                size=6) +
      transition_states(out_var_key) +
      ggtitle(paste0("VI for predicting ", "{closest_state}")) +
      ease_aes("linear")
    #====================================
    
    
    
    
    #====================================
    #Use 'magick' package to combine animated gifs
    a_plot <- animate(animated_importance_heatmap, height=4, width=5, units="in", res=250, nframes=100)
    b_plot <- animate(animated_importance_barplot, height=4, width=5, units="in", res=250, nframes=100)
    
    a_plot_file <- anim_save("importance_heatmap.gif", a_plot)
    b_plot_file <- anim_save("importance_barplot.gif", b_plot)
    
    a_mgif <- image_read("importance_heatmap.gif")
    b_mgif <- image_read("importance_barplot.gif")
    
    #Combine the animated heatmap and barplot into 1 animation
    new_gif <- image_append(c(a_mgif[1], b_mgif[1]), stack = FALSE)
    for(i in 2:100){ #Because there are 100 frames. This would change if each gif has a different number of frames
      combined <- image_append(c(a_mgif[i], b_mgif[i]), stack=FALSE)
      new_gif <- c(new_gif, combined)
    }
    
    #Save the combined animation
    anim_save("feature_importance.gif", 
              animation = new_gif)
    #====================================
   
    
  }#End animate importance plots
  
  

  #====================================
  #Return results
  return(list(
    "combined_prediction_df"   =   combined_prediction_df,
    "combined_model_metrics_df"=combined_model_metrics_df,
    "combined_importance_df"   =   combined_importance_df,
    "models_list"              =              models_list,
    "prediction_list"          =          prediction_list,
    "cor_plots_list"           =           cor_plots_list,
    "importance_barplot_list"  =  importance_barplot_list,
    "importance_heatmap_list"  =  importance_heatmap_list,
    "model_metrics_list"       =       model_metrics_list,
    "importance_list"          =          importance_list)
  )
  #====================================
  
  
  
}#End function
#--------------------------------------------------------------