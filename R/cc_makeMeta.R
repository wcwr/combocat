#--------------------------------------------------------------
#' Generate 'Complete' Metadata for Combocat Experiments
#'
#' This function generates "complete" type metadata for a Combocat experiment (either sparse or dense mode). 
#' It requires a very specific registration file containing plate registration and Echo transfer logs. The raw data files must be named as `[Plate Barcode].csv` where `[Plate Barcode]` is the corresponding column in the registration file.
#'
#' @param reg_file Data frame. A file containing plate registration and Echo transfer logs (`.csv` format).
#' @param data_mode Character. Specifies the data mode: `"dense"` or `"sparse"`. 
#' @param dilution_factor Numeric. The dilution factor used for all drugs. Default is `3` (i.e., 1:3 dilution). Multiple dilution factors are not supported for this function.
#' @param assay_ready_volume Character or Numeric. The total expected volume of the assay-ready plate in nanoliters (nL). Defaults to `"auto_detect"`, which uses the most frequent volume in the `Amount` column. In sparse mode, the standard volume is 20 nL. In dense mode, the standard volume is 200nL.
#' @param final_volume Character or Numeric. The final volume of the assay in NANOliters. Default is `"standard"` which is 40,000nL (40uL) in dense mode or 4,000nL (4uL) sparse mode. 
#' @param sample_name Character. The name of the sample. Defaults to `"not_provided"`.
#'
#' @export
#--------------------------------------------------------------




#--------------------------------------------------------------
#A function to generate 'complete' metadata for a Combocat experiment (either sparse or dense)
#NOTE: This function is dependent on a very specific registration file (containing plate registration + Echo transfer logs)
#NOTE: The associated raw data files MUST be named [`Plate Barcode`].csv where [`Plate Barcode`] is the column of the same name
#Contact Charlie Wright (charlie.wright@stjude.org) for an example input file. 
cc_makeMeta <- function(reg_file,
                        data_mode="dense",
                        dilution_factor=3,
                        assay_ready_volume="auto_detect",
                        final_volume="standard",
                        sample_name="not_provided"){
  
  
  #Input requirements:
  #reg_file------------A file containing plate registration and Echo transfer logs (.csv)
  #data_mode-----------Either "dense" or "sparse"
  #dilution_factor-----The dilution factor used for all drugs. Default is 3 (aka 1:3 dilution). If multiple dilution factors, user is responsible to adjust the output manually
  #assay_ready_volume--The expected total volume of the assay-ready plate in NANOliters. This is used to flag transfer failures. 'auto_detect' is default, which takes the mode, aka most frequent volume. (20nL in standard sparse mode)
  #final_volume--------The final volume of the assay in NANOliters. Default is "standard" which is 40,000nL (40uL) in dense mode or 4000nL (4uL) sparse mode. 
  #sample_name---------The name of the sample.
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Determine whether data is dense or sparse mode
  if(data_mode == "sparse"){
    print("Detected data type is SPARSE MODE. Proceeding with sparse mode metadata generation")
  } else {
    print("Detected data type is DENSE MODE. Proceeding with dense mode metadata generation")}
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
    #If assay_ready_volume is 'auto_detect' (default), set it to the mode : most frequent volume in 'Amount' column
    if(assay_ready_volume=="auto_detect"){
      assay_ready_volume <- names(sort(-table(reg_file$Amount)))[1] %>% as.numeric()
    }
    print(paste0("Assay-ready plate volume is set to ", assay_ready_volume, " nL"))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #If final_volume is 'standard' (default), set it to 40000nL (40uL) for dense mode
    if(final_volume=="standard"){
      final_volume <- 40000
    }
    print(paste0("Final volume is set to ", final_volume, " nL"))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize an empty data frame
    #This is what the final results will go into
    final_dataframe <- data.frame()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Assign 'sample_name' column
    reg_file$sample_name <- sample_name
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Assign 'filename' column (`Plate Barcode` + ".csv")
    #NOTE: Be sure the csv files are named in this way
    reg_file$filename <- paste0(reg_file$`Plate Barcode`, ".csv")
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Assign 'units' column
    #The calculation is set up to provide unit concentrations in uM
    reg_file$units <- "uM"  #Use 'u' instead of 'mu'/micro symbol to make sure nothing crashes
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Calculate the 'final' concentrations of each drug
    
    #NOTE: We have to make some modifications to get the correct values before we do the calculation:
    #1.... We need to double the 'Substance1/2 Concentration' (they are halved in the report)
    #2.... We multiply the 'Substance1/2 Concentration' by 1000 to convert from mM to uM
    #3.... We need to halve the 'Amount' column (since the report is the TOTAL transferred volume, not actual volume of the given drug)
    
    #Then we use C1V1=C2V2 to calculate the final concentration
    
    #Note, this will produce NAs in drug2 wells where there is no drug (single-agents)
    reg_file <- 
      reg_file %>%
      mutate("drug1_conc"= (2*(1000*`Substance1 Concentration`)) * (Amount/2) / final_volume) %>% 
      mutate("drug2_conc"= (2*(1000*`Substance2 Concentration`)) * (Amount/2) / final_volume)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Change the drug name columns to the easier-to-read "drug1_name/drug2_name"
    reg_file <- 
      reg_file %>%
      rename("drug1_name" = `Substance1 Compound Multipart 1`,
             "drug2_name" = `Substance2 Compound Multipart 1`)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Concatenate the drug names to form a unique 'combo_id'
    reg_file <- 
      reg_file %>%
      mutate("combo_id"=paste0(`drug1_name`, "_", `drug2_name`))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Convert the column 'Well Alpha' from current format (e.g. 'AA10') to columns 'plate_row' and 'plate_col'
    #Also keep the letters just for sanity check, as a column "plate_row_letter"
    
    
    #Make a dataframe that will be used to convert the row letters to numbers
    row_conversion_df <- data.frame(
      "row_letter" = c(LETTERS[1:26], "AA", "AB", "AC", "AD", "AE", "AF"),
      "row_number" = seq(1:32)
    )
    
    #convert 'Well Alpha' to 'plate_row_letter' and 'plate_col'
    reg_file <- 
      reg_file %>%
      mutate("plate_row_letter"=stringr::str_extract(`Well Alpha`, "^[A-Za-z]+"),
             "plate_col" = as.numeric(str_extract(`Well Alpha`, "[0-9]+$")))
    
    # Join the reg_file with row_conversion_df
    reg_file <- 
      reg_file %>%
      left_join(row_conversion_df, by = c("plate_row_letter" = "row_letter")) %>%
      rename("plate_row" = row_number)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Make a reference data frame of all the wells that are expected to have data
    #For dense mode, each plate will be full and we expect data from 384 wells
    reference_df <- 
      expand.grid(
        plate_row = 1:16,
        plate_col = 1:24
      ) %>%
      mutate("drug1_nth_conc"=NA,
             "drug2_nth_conc"=NA,
             "transfer_type"=NA)
    
    #NOTE: Remove control wells
    #We assume controls transferred correctly, or if not they will be caught by outlier removal
    #Therefore, we expect 360 wells/values of data
    reference_df <- 
      reference_df %>%
      filter(!(plate_row %in% 11:16 & plate_col %in% 21:24))
    
    
    #Explicitly define where each drug1_nth_conc and drug2_nth_conc should be: 
    reference_df <-
      reference_df %>%
      mutate("drug1_nth_conc" = case_when(
        
        #Handle the drug1_nth_conc values (for combination wells)
        plate_row %in% 1:10  & plate_col %in% 1:3   ~1,
        plate_row %in% 1:10  & plate_col %in% 4:6   ~2,
        plate_row %in% 1:10  & plate_col %in% 7:9   ~3,
        plate_row %in% 1:10  & plate_col %in% 10:12 ~4,
        plate_row %in% 1:10  & plate_col %in% 13:15 ~5,
        plate_row %in% 1:10  & plate_col %in% 16:18 ~6,
        plate_row %in% 1:10  & plate_col %in% 19:21 ~7,
        plate_row %in% 1:10  & plate_col %in% 22:24 ~8,
        plate_row %in% 11:13 & plate_col %in% 1:10  ~9,   #Transposed section (drug1)
        plate_row %in% 14:16 & plate_col %in% 1:10  ~10,  #Transposed section (drug1)
        
        #Handle the drug1_nth_conc values (for single-agent wells)
        plate_row %in% 11:13 & plate_col==11 ~ 1, 
        plate_row %in% 11:13 & plate_col==12 ~ 2, 
        plate_row %in% 11:13 & plate_col==13 ~ 3, 
        plate_row %in% 11:13 & plate_col==14 ~ 4, 
        plate_row %in% 11:13 & plate_col==15 ~ 5, 
        plate_row %in% 11:13 & plate_col==16 ~ 6, 
        plate_row %in% 11:13 & plate_col==17 ~ 7, 
        plate_row %in% 11:13 & plate_col==18 ~ 8, 
        plate_row %in% 11:13 & plate_col==19 ~ 9, 
        plate_row %in% 11:13 & plate_col==20 ~ 10)) %>%
      
      mutate("drug2_nth_conc" = case_when(
        
        #Handle the drug2_nth_conc values (for combination wells)
        plate_row==1  & plate_col %in% 1:24 ~ 1,
        plate_row==2  & plate_col %in% 1:24 ~ 2,
        plate_row==3  & plate_col %in% 1:24 ~ 3,
        plate_row==4  & plate_col %in% 1:24 ~ 4,
        plate_row==5  & plate_col %in% 1:24 ~ 5,
        plate_row==6  & plate_col %in% 1:24 ~ 6,
        plate_row==7  & plate_col %in% 1:24 ~ 7,
        plate_row==8  & plate_col %in% 1:24 ~ 8,
        plate_row==9  & plate_col %in% 1:24 ~ 9,
        plate_row==10 & plate_col %in% 1:24 ~ 10,
        plate_row %in% 11:16 & plate_col==1 ~ 1,  #Transposed section (drug2)
        plate_row %in% 11:16 & plate_col==2 ~ 2,  #Transposed section (drug2)
        plate_row %in% 11:16 & plate_col==3 ~ 3,  #Transposed section (drug2)
        plate_row %in% 11:16 & plate_col==4 ~ 4,  #Transposed section (drug2)
        plate_row %in% 11:16 & plate_col==5 ~ 5,  #Transposed section (drug2)
        plate_row %in% 11:16 & plate_col==6 ~ 6,  #Transposed section (drug2)
        plate_row %in% 11:16 & plate_col==7 ~ 7,  #Transposed section (drug2)
        plate_row %in% 11:16 & plate_col==8 ~ 8,  #Transposed section (drug2)
        plate_row %in% 11:16 & plate_col==9 ~ 9,  #Transposed section (drug2)
        plate_row %in% 11:16 & plate_col==10 ~ 10, #Transposed section (drug2)
        
        #Handle the drug2_nth_conc values (for single-agent wells)
        plate_row %in% 14:16 & plate_col==11 ~1,
        plate_row %in% 14:16 & plate_col==12 ~2,
        plate_row %in% 14:16 & plate_col==13 ~3,
        plate_row %in% 14:16 & plate_col==14 ~4,
        plate_row %in% 14:16 & plate_col==15 ~5,
        plate_row %in% 14:16 & plate_col==16 ~6,
        plate_row %in% 14:16 & plate_col==17 ~7,
        plate_row %in% 14:16 & plate_col==18 ~8,
        plate_row %in% 14:16 & plate_col==19 ~9,
        plate_row %in% 14:16 & plate_col==20 ~10)) %>%
      
      mutate("replicate" = case_when(
        #Handle the replicate values:
        plate_row %in% 1:10 & plate_col %in% c(1,4,7,10,13,16,19,22) ~ 1,
        plate_row %in% 1:10 & plate_col %in% c(2,5,8,11,14,17,20,23) ~ 2,
        plate_row %in% 1:10 & plate_col %in% c(3,6,9,12,15,18,21,24) ~ 3,
        
        plate_row %in% c(11,14) & plate_col %in% 1:10 ~ 1,  #Transposed section
        plate_row %in% c(12,15) & plate_col %in% 1:10 ~ 2,  #Transposed section
        plate_row %in% c(13,16) & plate_col %in% 1:10 ~ 3,  #Transposed section
        
        #Handle replicates for single-agents:
        plate_row %in% c(11,14) & plate_col %in% 11:20 ~1,      #Single-agents rep1
        plate_row %in% c(12,15) & plate_col %in% 11:20 ~2,      #Single-agents rep2
        plate_row %in% c(13,16) & plate_col %in% 11:20 ~3)) %>% #Single-agents rep3
      
      mutate("transfer_type"= case_when(
        #Handle the transfer_type (single-agent or combination)
        plate_row %in% 11:13 & plate_col %in% 11:20 ~ "single_agent_drug1", 
        plate_row %in% 14:16 & plate_col %in% 11:20 ~ "single_agent_drug2",
        plate_row %in% 1:10  & plate_col %in% 1:24  ~ "combination",
        plate_row %in% 11:16 & plate_col %in% 1:10  ~ "combination"   #Transposed section (all combination wells)
      
      ))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Loop over each plate
    for(i in unique(reg_file$`Plate Barcode`)){
      
      
      
      #====================================
      #Subset the reg_file to the current plate
      reg_file_subset <- 
        reg_file %>%
        filter(`Plate Barcode`==i)
      #====================================
      
      
      
      
      #====================================
      #Merge the reference_df with the reg_file_subset
      #This will keep all the reference wells, even if they don't have a transfer
      reg_file_subset2 <- 
        merge(reference_df,
              reg_file_subset,
              by=c("plate_row", "plate_col"),
              all.x=TRUE) %>% #Keep all reference wells, even if they don't have a transfer
        select(plate_row,
               plate_col,
               plate_row_letter,
               transfer_type,
               replicate,
               drug1_name,
               drug1_conc,
               drug1_nth_conc,
               drug2_name,
               drug2_conc,
               drug2_nth_conc,
               units,
               combo_id,
               `Plate Barcode`,
               sample_name,
               Amount,
               filename)
      #====================================
      
      
      
      
      #====================================
      #Figure out what the expected drug1/2 names are
      #We use the most common names in the drug1-only or drug2-only wells to figure this out
      
      #Figure out what d1_name should be: 
      d1_name <- 
        reg_file_subset2 %>%
        filter(transfer_type=="single_agent_drug1") %>%
        add_count(drug1_name) %>% #Count number of times each drug1_name appears
        arrange(desc(n)) %>%      #Sort by the count (descending)
        slice(1) %>%              #Select the first row (which has the highest count of drug1_name)
        pull(drug1_name)
      
      #Figure out what d2_name should be:
      d2_name <- 
        reg_file_subset2 %>%
        filter(transfer_type=="single_agent_drug2") %>%
        add_count(drug1_name) %>% #Count number of times each drug1_name appears #NOTE: IT'S DRUG1 NOT DRUG2! (because Echo places it into drug1 column when only 1 drug is transferred)
        arrange(desc(n)) %>%      #Sort by the count (descending)
        slice(1) %>%              #Select the first row (which has the highest count of drug2_name)
        pull(drug1_name)
      #====================================
      
      
      
      
      #====================================
      #Assign the drug names and concentrations in appropriate places
      #By default, some of these are in the wrong place even when transfers are successful
      #This is due to how the Echo reports the transfers
      
      #These are the four conditions where drug names and concentrations are correct or need to be swapped:
      
      #1. If transfer_type is "single_agent_drug1", then drug1_name, drug1_conc, drug2_name, & drug2_conc are CORRECT
      #2. If transfer_type is "single_agent_drug2", then switch drug1_name, drug1_conc with drug2_name, drug2_conc
      #3. If drug1_name==d1_name, then drug1_name and drug1_conc are CORRECT
      #4. If drug1_name==d2_name, then switch drug1_name, drug1_conc with drug2_name, drug2_conc
      
      #(This is done row-wise and can take some time to complete)
      reg_file_subset3 <- 
        reg_file_subset2 %>%
        rowwise() %>%
        mutate(
          new_drug1_name = case_when(
            transfer_type == "single_agent_drug1" ~ drug1_name,
            transfer_type == "single_agent_drug2" ~ drug2_name,
            drug1_name == d1_name ~ drug1_name,
            drug1_name == d2_name ~ drug2_name,
            TRUE ~ drug1_name
          ),
          new_drug1_conc = case_when(
            transfer_type == "single_agent_drug1" ~ drug1_conc,
            transfer_type == "single_agent_drug2" ~ drug2_conc,
            drug1_name == d1_name ~ drug1_conc,
            drug1_name == d2_name ~ drug2_conc,
            TRUE ~ drug1_conc
          ),
          new_drug2_name = case_when(
            transfer_type == "single_agent_drug1" ~ drug2_name,
            transfer_type == "single_agent_drug2" ~ drug1_name,
            drug1_name == d1_name ~ drug2_name,
            drug1_name == d2_name ~ drug1_name,
            TRUE ~ drug2_name
          ),
          new_drug2_conc = case_when(
            transfer_type == "single_agent_drug1" ~ drug2_conc,
            transfer_type == "single_agent_drug2" ~ drug1_conc,
            drug1_name == d1_name ~ drug2_conc,
            drug1_name == d2_name ~ drug1_conc,
            TRUE ~ drug2_conc
          )
        ) %>%
        ungroup() %>%
        # Replace the old columns with the new ones
        select(-drug1_name, -drug1_conc, -drug2_name, -drug2_conc) %>%
        rename(
          drug1_name = new_drug1_name, 
          drug1_conc = new_drug1_conc,
          drug2_name = new_drug2_name, 
          drug2_conc = new_drug2_conc
        ) %>%
        #Put columns in preferred order
        select(plate_row,
               plate_col,
               plate_row_letter,
               drug1_name,
               drug1_conc,
               drug2_name,
               drug2_conc,
               everything())
      #====================================
        
      
  
      
      #====================================
      #Flag the wells with a failed transfer
      #Any well with an 'Amount' different than assay_ready_volume is a failed transfer
      #For single-agent transfers, a sign of failed transfer is when 'Amount' is NA
      #For combination transfers, a sign of failed transfer is when 'Amount' is not assay_ready_volume (or NA).
      #Just to note: In combination transfers, if 1 drug fails to transfer the 'Amount' will be 1/2 of assay_ready_volume
      
      reg_file_subset3 <- 
        reg_file_subset3 %>%
        mutate("failed_transfer" = ifelse(is.na(Amount) | !(Amount == assay_ready_volume), 1, 0))
      #====================================
      
      
      
      
      #====================================
      #Fix the various missing/incorrect values caused by transfer failures
      
      #The drug1_name becomes what we previously defined as d1_name
      #This is done ONLY if the transfer failed in 'combination' or 'single_agent_drug1' transfers 
      #............(because drug1_name should be NA in cases of transfer_type==single_agent_drug2)
      reg_file_subset3 <-
        reg_file_subset3 %>%
        mutate("drug1_name"= case_when(
          !transfer_type=="single_agent_drug2" & failed_transfer == 1 ~ d1_name,
          TRUE ~ drug1_name
        )) 
      
      #Same approach for drug2_name
      reg_file_subset3 <-
        reg_file_subset3 %>%
        mutate("drug2_name"= case_when(
          !transfer_type=="single_agent_drug1" & failed_transfer == 1 ~ d2_name,
          TRUE ~ drug2_name
        ))

  
      #Fix any missing entries in 'filename'
      #There should only be 1 filename per plate in dense mode
      
      #Get the known filename
      known_filename <- 
        reg_file_subset3 %>%
        filter(!is.na(filename)) %>%
        pull(filename) %>%
        unique()
      
      reg_file_subset3 <-
        reg_file_subset3 %>%
        mutate("filename"= case_when(
          failed_transfer == 1 ~ known_filename,
          TRUE ~ filename
        ))
      
      
      #Same method to fix any missing sample_name
      known_sample_name <-
        reg_file_subset3 %>%
        filter(!is.na(sample_name)) %>%
        pull(sample_name) %>%
        unique()
      
      reg_file_subset3 <-
        reg_file_subset3 %>%
        mutate("sample_name"= case_when(
          failed_transfer == 1 ~ known_sample_name,
          TRUE ~ sample_name
        ))
      
      
      #Same method to fix any missing `Plate Barcode`
      known_plate_barcode <-
        reg_file_subset3 %>%
        filter(!is.na(`Plate Barcode`)) %>%
        pull(`Plate Barcode`) %>%
        unique()
      
      reg_file_subset3 <-
        reg_file_subset3 %>%
        mutate("Plate Barcode"= case_when(
          failed_transfer == 1 ~ known_plate_barcode,
          TRUE ~ `Plate Barcode`
        ))
      
      
      #Same method to fix any missing units
      known_units <- 
        reg_file_subset3 %>%
        filter(!is.na(units)) %>%
        pull(units) %>%
        unique()
      
      reg_file_subset3 <-
        reg_file_subset3 %>%
        mutate("units"= case_when(
          failed_transfer == 1 ~ known_units,
          TRUE ~ units
        ))
      
      
      #Fix missing plate_row_letter values
      reg_file_subset3 <- 
        reg_file_subset3 %>%
        mutate("plate_row_letter" = ifelse(is.na(plate_row_letter),
                                           LETTERS[plate_row],
                                           plate_row_letter))
      
      
      #====================================
      
      
      
      
      #====================================
      #Fill in the missing concentrations
      
      #Figure out which nth_conc columns are associated with missing drug1/2_conc values
      mising_nth_conc_drug1 <- 
        reg_file_subset3 %>% 
        filter(!transfer_type=="single_agent_drug2")  %>% 
        filter(is.na(drug1_conc)) %>%
        pull(drug1_nth_conc) %>% 
        unique()
      
      mising_nth_conc_drug2 <- 
        reg_file_subset3 %>% 
        filter(!transfer_type=="single_agent_drug1")  %>% 
        filter(is.na(drug2_conc)) %>%
        pull(drug2_nth_conc) %>% 
        unique()
      
      #Make a temporary df that is in the compatible format for the calculate_concentration function
      #The dataframe should have the columns, 'drug1/2_name', 'drug1/2_conc', and 'nth_conc'
      tmp_df_drug1 <- 
        reg_file_subset3 %>%
        filter(!is.na(drug1_nth_conc)) %>% #Removes single_agent_drug1 (where drug1_nth_conc is supposed to be NA)
        select(drug1_name, drug1_conc, drug1_nth_conc) %>%
        unique() %>%
        rename("nth_conc" = drug1_nth_conc)
      
      
      tmp_df_drug2 <- 
        reg_file_subset3 %>%
        filter(!is.na(drug2_nth_conc)) %>% #Removes single_agent_drug2 (where drug2_nth_conc is supposed to be NA)
        select(drug2_name, drug2_conc, drug2_nth_conc) %>%
        unique() %>%
        rename("nth_conc" = drug2_nth_conc)
      
      
      #Fill in each missing concentration
      #(Nothing happens to the dataframes if there are no missing concentrations)
      for(nth in mising_nth_conc_drug1) {
        calculated_conc <- calculate_concentration(nth, tmp_df_drug1, "drug1_conc", dilution_factor=dilution_factor)
        tmp_df_drug1$drug1_conc[tmp_df_drug1$nth_conc == nth] <- calculated_conc
      }
      
      for(nth in mising_nth_conc_drug2) {
        calculated_conc <- calculate_concentration(nth, tmp_df_drug2, "drug2_conc", dilution_factor=dilution_factor)
        tmp_df_drug2$drug2_conc[tmp_df_drug2$nth_conc == nth] <- calculated_conc
      }
      
      
      #Change the 'nth_conc' column names back to 'drug1/2_nth_conc'
      tmp_df_drug1 <- 
        tmp_df_drug1 %>%
        rename("drug1_nth_conc" = nth_conc)
      
      tmp_df_drug2 <- 
        tmp_df_drug2 %>%
        rename("drug2_nth_conc" = nth_conc)
      
      
      # Update the drug1_conc and drug2_conc in reg_file_subset3 using the calculated values
      #This will fill in the missing drug1/2_conc values
      #Note that due to floating point precision errors, the calculated values may be e.g. 24.9999998 rather than 25
      reg_file_subset3 <- 
        reg_file_subset3 %>%
        mutate(
          drug1_conc = ifelse(!is.na(drug1_nth_conc) & drug1_nth_conc %in% mising_nth_conc_drug1, tmp_df_drug1$drug1_conc[match(drug1_nth_conc, tmp_df_drug1$drug1_nth_conc)], drug1_conc),
          drug2_conc = ifelse(!is.na(drug2_nth_conc) & drug2_nth_conc %in% mising_nth_conc_drug2, tmp_df_drug2$drug2_conc[match(drug2_nth_conc, tmp_df_drug2$drug2_nth_conc)], drug2_conc)
        )
      #====================================
      
      
      
      
      #====================================
      #Some adjustments to make this version of metadata more similar to 'simplified' metadata:
      
      #Fill in all 'NA' values in drug1/2_name columns
      #We have already done all corrections at this point, so it is safe to assign all rows as d1/d2_name
      reg_file_subset3 <- 
        reg_file_subset3 %>%
        mutate(
          drug1_name = ifelse(is.na(drug1_name), d1_name, drug1_name),
          drug2_name = ifelse(is.na(drug2_name), d2_name, drug2_name)
        )
      
      #Fill in the 'NA' values in the drug1/2_conc columns
      #We have already corrected any 'true' missing concentrations
      #Now any 'NA' concentrations will correspond to single-agent transfers, so we set them to 0
      reg_file_subset3 <-
        reg_file_subset3 %>%
        mutate(
          drug1_conc = case_when(
            transfer_type=="single_agent_drug2" ~ 0,
            TRUE ~ drug1_conc
          )) %>%
        mutate(
          drug2_conc = case_when(
            transfer_type=="single_agent_drug1" ~ 0,
            TRUE ~ drug2_conc
          ))
      
      #Make the 'combo_id' column as [d1_name]_[d2_name], even for single-agents, as opposed to [d1_name]_NA
      #This is to make the column compatible with 'simplified' type metadata
      reg_file_subset3 <- 
        reg_file_subset3 %>%
        mutate("combo_id"=paste0(d1_name, "_", d2_name))
      #====================================
      
        
      
      
      #====================================
      #Deposit the corrected results into the final dataframe
      final_dataframe <- bind_rows(final_dataframe, reg_file_subset3)
      #====================================
      

      
    }#End loop over each plate
    
    
    
  } else {
    
    #=================================================================
    #=================================================================
    #
    #-------------------Handle sparse mode data-----------------------
    #
    #=================================================================
    #=================================================================
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #If assay_ready_volume is 'auto_detect' (default), set it to the mode : most frequent volume in 'Amount' column
    if(assay_ready_volume=="auto_detect"){
      assay_ready_volume <- names(sort(-table(reg_file$Amount)))[1] %>% as.numeric()
    }
    print(paste0("Assay-ready plate volume is set to ", assay_ready_volume, " nL"))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #If final_volume is 'standard' (default), set it to 4000nL (4uL) for sparse mode
    if(final_volume=="standard"){
      final_volume <- 4000
    }
    print(paste0("Final volume is set to ", final_volume, " nL"))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Initialize an empty data frame
    #This is what the final results will go into
    final_dataframe <- data.frame()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Assign 'sample_name' column
    reg_file$sample_name <- sample_name
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Assign 'filename' column (`Plate Barcode` + ".csv")
    #NOTE: Be sure the csv files are named in this way
    reg_file$filename <- paste0(reg_file$`Plate Barcode`, ".csv")
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Assign 'units' column
    #The calculation is set up to provide unit concentrations in uM
    reg_file$units <- "uM"  #Use 'u' instead of 'mu'/micro symbol to make sure nothing crashes
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Calculate the 'final' concentrations of each drug
    
    #NOTE: We have to make some modifications to get the correct values before we do the calculation:
    #1.... We need to double the 'Substance1/2 Concentration' (they are halved in the report)
    #2.... We multiply the 'Substance1/2 Concentration' by 1000 to convert from mM to uM
    #3.... We need to halve the 'Amount' column (since the report is the TOTAL transferred volume, not actual volume of the given drug)
    
    #Then we use C1V1=C2V2 to calculate the final concentration
    
    #Note, this will produce NAs in drug2 wells where there is no drug (single-agents)
    reg_file <- 
      reg_file %>%
      mutate("drug1_conc"= (2*(1000*`Substance1 Concentration`)) * (Amount/2) / final_volume) %>% 
      mutate("drug2_conc"= (2*(1000*`Substance2 Concentration`)) * (Amount/2) / final_volume)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Change the drug name columns to the easier-to-read "drug1_name/drug2_name"
    reg_file <- 
      reg_file %>%
      rename("drug1_name" = `Substance1 Compound Multipart 1`,
             "drug2_name" = `Substance2 Compound Multipart 1`)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Concatenate the drug names to form a unique 'combo_id'
    reg_file <- 
      reg_file %>%
      mutate("combo_id"=paste0(`drug1_name`, "_", `drug2_name`))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Determine whether plate_type is single_agent or combination
    #Grouping by Plate Barcode, only if ALL drug2_conc are NA, then the plate is single_agent
    reg_file <- 
      reg_file %>%
      group_by(`Plate Barcode`) %>%
      mutate("plate_type" = ifelse(all(is.na(drug2_conc)), "single_agent", "combination")) %>%
      ungroup() 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Convert the column 'Well Alpha' from current format (e.g. 'AA10') to columns 'plate_row' and 'plate_col'
    #Also keep the letters just for sanity check, as a column "plate_row_letter"
    
    
    #Make a dataframe that will be used to convert the row letters to numbers
    row_conversion_df <- data.frame(
      "row_letter" = c(LETTERS[1:26], "AA", "AB", "AC", "AD", "AE", "AF"),
      "row_number" = seq(1:32)
    )
    
    #convert 'Well Alpha' to 'plate_row_letter' and 'plate_col'
    reg_file <- 
      reg_file %>%
      mutate("plate_row_letter"=stringr::str_extract(`Well Alpha`, "^[A-Za-z]+"),
             "plate_col" = as.numeric(str_extract(`Well Alpha`, "[0-9]+$")))
    
    # Join the reg_file with row_conversion_df
    reg_file <- 
      reg_file %>%
      left_join(row_conversion_df, by = c("plate_row_letter" = "row_letter")) %>%
      rename("plate_row" = row_number)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Make a reference data frame of all the wells that are expected to have data
    #For full plates, there should be 1350 total wells (we still ignore controls and outer wells)
    reference_df <- 
      expand.grid(
        plate_row = 2:31,
        plate_col = 2:46
      ) %>%
      mutate("position_id" = calculate_position_id(plate_row, plate_col)) #Also get the position_id
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Loop over each plate
    for(i in unique(reg_file$`Plate Barcode`)){
      
      
      
      #====================================
      #Subset the data to get rid of wells we don't care about: Outer wells and control wells:
      #This should leave 1350 rows per plate (or fewer if some transfers failed)
      reg_file_subset <- 
        reg_file %>%
        filter(`Plate Barcode`==i)        %>% #Subset to the current plate
        filter(!(plate_row %in% c(1,32))) %>% #Remove outer rows
        filter(!(plate_col %in% c(1,48))) %>% #Remove outer columns
        filter(!plate_col==47)                #Removes control wells column
      #====================================
      
      
      
      
      #====================================
      #Use the calculate_position_id function to get the position_id column
      #Note: for a full plate, the smallest position_id should be 1 and the largest 135
      reg_file_subset <- 
        reg_file_subset %>%
        mutate("position_id" = calculate_position_id(plate_row, plate_col))
      #====================================
      
      
      
      
      #====================================
      #Merge the reference_df with the reg_file_subset
      #This will keep all the reference wells, even if they don't have a transfer
      reg_file_subset2 <- 
        merge(reference_df,
              reg_file_subset %>%
                select(plate_row,
                       plate_col,
                       plate_row_letter, 
                       `drug1_name`, 
                       drug1_conc,
                       `drug2_name`,
                       drug2_conc,
                       units,
                       combo_id,
                       `Plate Barcode`,
                       sample_name,
                       Amount,
                       plate_type,
                       position_id,
                       filename),
              by=c("plate_row", "plate_col", "position_id"),
              all.x = TRUE) #Keep all reference wells, even if they don't have a transfer
      #====================================
      
      
      
      
      #====================================
      #Add the nth_conc column which tells which relative concentration is expected (1 being highest, 10 being lowest)
      #We go by row (if plate_row is 1:3, nth_conc should be 1 .... if plate_row is 29:31, nth_conc should be 10)
      reg_file_subset2 <- 
        reg_file_subset2 %>%
        mutate(nth_conc = case_when(
          plate_row %in% 2:4   ~ 1,
          plate_row %in% 5:7   ~ 2,
          plate_row %in% 8:10  ~ 3,
          plate_row %in% 11:13 ~ 4,
          plate_row %in% 14:16 ~ 5,
          plate_row %in% 17:19 ~ 6,
          plate_row %in% 20:22 ~ 7,
          plate_row %in% 23:25 ~ 8,
          plate_row %in% 26:28 ~ 9,
          plate_row %in% 29:31 ~ 10
        ))
      #====================================
      
      
      
      
      #====================================
      #Flag the wells with a failed transfer
      
      #Any well with an 'Amount' different than assay_ready_volume is a failed transfer
      #For single-agent plates, a sign of failed transfer is when 'Amount' is NA
      #For combination plates,  a sign of failed transfer is when 'Amount' is not assay_ready_volume (or NA).
      #Just to note: In combination plates, if 1 drug fails to transfer the 'Amount' will be 1/2 of assay_ready_volume
      
      reg_file_subset2 <- 
        reg_file_subset2 %>%
        mutate("failed_transfer" = ifelse(is.na(Amount) | !(Amount == assay_ready_volume), 1, 0))
      #====================================
      
      
      
      
      #---(Still within plate-wise loop)---
      #Loop over the position_ids
      for(k in unique(reg_file_subset2$position_id)){
        
        
        #====================================
        #Subset the current plate to the current position_id
        reg_file_subset3 <- 
          reg_file_subset2 %>%
          filter(position_id==k)
        #====================================
        
        
        
        
        #====================================
        #Fix the various missing/incorrect values caused by transfer failures
        
        #The drug1_name becomes whatever is the most common in that column
        #Ex: if names are (a,a,a,a,b,a,a,a,a,a), then the drug1_name becomes "a"
        reg_file_subset3 <- 
          reg_file_subset3 %>%
          add_count(drug1_name) %>%
          arrange(desc(n)) %>%
          mutate("drug1_name" = first(drug1_name)) %>%
          select(-n)
        
        #Same approach for drug2_name
        #(Note that in single_agent plates, drug2_name should end up being NA)
        reg_file_subset3 <- 
          reg_file_subset3 %>%
          add_count(drug2_name) %>%
          arrange(desc(n)) %>%
          mutate("drug2_name" = first(drug2_name)) %>%
          select(-n)
        
        #Same approach for combo_id
        reg_file_subset3 <- 
          reg_file_subset3 %>%
          add_count(combo_id) %>%
          arrange(desc(n)) %>%
          mutate("combo_id" = first(combo_id)) %>%
          select(-n)
        
        #Same approach for filename
        reg_file_subset3 <- 
          reg_file_subset3 %>%
          add_count(filename) %>%
          arrange(desc(n)) %>%
          mutate("filename" = first(filename)) %>%
          select(-n)
        
        #Same approach for sample_name
        reg_file_subset3 <- 
          reg_file_subset3 %>%
          add_count(sample_name) %>%
          arrange(desc(n)) %>%
          mutate("sample_name" = first(sample_name)) %>%
          select(-n)
        
        #Same approach for units
        reg_file_subset3 <- 
          reg_file_subset3 %>%
          add_count(units) %>%
          arrange(desc(n)) %>%
          mutate("units" = first(units)) %>%
          select(-n)
        #====================================
        
        
        
        
        #====================================
        #Also fix any missing information in Plate Barcode or plate_type columns
        #(This should really just affect the single_agent plates, since those plates have cases where entire rows are not output for failed transfers)
        
        #Same approach for Plate Barcode
        reg_file_subset3 <- 
          reg_file_subset3 %>%
          add_count(`Plate Barcode`) %>%
          arrange(desc(n)) %>%
          mutate(`Plate Barcode` = first(`Plate Barcode`)) %>%
          select(-n)
        
        #Same approach for plate_type
        reg_file_subset3 <- 
          reg_file_subset3 %>%
          add_count(plate_type) %>%
          arrange(desc(n)) %>%
          mutate(plate_type = first(plate_type)) %>%
          select(-n)
        #====================================
        
        
        
        
        #====================================
        #Fill in the missing concentrations
        #We do this for drug1_conc in either plate_type, but only 'combination' plates apply this for drug2_conc
        
        #CRITICAL NOTE: 
        #Some concentrations that get filled in might be very slightly different than what they are intended to be
        #This is due to floating point precision issues, and occurs in wells with failed transfers
        #For example, 0.111111112 might be the intended concentration, but after division we might get a result of 0.111111113
        #This shouldn't be a problem at this stage, but will be a problem downstream when concentrations need to be matched exactly
        #For example: 
        #.............Combining DrugA at 0.111111112uM + DrugB at 0.111111112uM (same concs) could give 80% cell death, 
        #.............So to calculate synergy, we would first grab the single-agent effects,
        #.............First look for DrugA alone at 0.111111112uM then DrugB alone at 0.111111112uM
        #.............But if there is no DrugA concentration of 0.111111112uM, because it was calculated to be 0.111111113, 
        #.............This will cause an error
        #.............We handle this by rounding to the 6th decimal place in the calculate_concentration() function
        #.............Rounding to 6 decimal places should be sufficient to handle floating point precision issues (can be changed)
        
        
        
        
        #First define which nth_conc's are missing (for drug1_conc)
        #CRITICAL: For missing nth_concs for drug1, we actually consider drug1_conc AND drug2_conc
        #Reason : 
        #If drug1_conc is na, then it is definitely missing. But ALSO...
        #If drug2_conc is NA, that could be from one of two cases: 
        #.....................First case is that drug2 simply did not transfer
        #.....................Second case is that drug1 did not transfer, and the information of drug2 got shifted to drug1
        #.....................Therefore, any drug2_conc that is NA (only in the case of combination plates), that nth_conc also gets marked as a missing nth_conc for drug1
        missing_nth_conc_drug1 <- 
          c(
            #nth_conc where drug1_conc is NA
            reg_file_subset3 %>%
              filter(is.na(drug1_conc)) %>%
              pull(nth_conc),
            
            #nth_conc where drug2_conc is NA (only in the case of combination plates)
            reg_file_subset3 %>%
              filter(plate_type=="combination") %>%
              filter(is.na(drug2_conc)) %>%
              pull(nth_conc)
          )
        
        
        #Now get the missing nth_conc drug2
        #This will be just like above, grab any nth_concs where drug2_conc is NA (for combinatin plates, only)
        missing_nth_conc_drug2 <- 
          reg_file_subset3 %>%
          filter(plate_type=="combination") %>%
          filter(is.na(drug2_conc)) %>%
          pull(nth_conc)
        
        
        
        
        # Loop over each missing nth_conc value to calculate the missing concentration
        for(nth in missing_nth_conc_drug1) {
          calculated_conc <- calculate_concentration(nth, reg_file_subset3, "drug1_conc", dilution_factor=dilution_factor)
          # Now, update the concentration in reg_file_subset3 for the corresponding nth_conc
          reg_file_subset3$drug1_conc[reg_file_subset3$nth_conc == nth] <- calculated_conc
        }
        
        
        #Repeat this specific for combination plate types only for drug2_conc
        if(all(reg_file_subset3$plate_type == "combination", na.rm = TRUE)){
          
          
          for(nth in missing_nth_conc_drug2) {
            calculated_conc <- calculate_concentration(nth, reg_file_subset3, "drug2_conc", dilution_factor=dilution_factor)
            reg_file_subset3$drug2_conc[reg_file_subset3$nth_conc == nth] <- calculated_conc
          }
        }
        #====================================
        
        
        
        
        #====================================
        #Deposit the corrected results into the final dataframe
        final_dataframe <- bind_rows(final_dataframe, reg_file_subset3)
        #====================================
        
        
      }#End loop over position_ids
      

    }#End loop over plate barcodes
    
    
    #====================================
    #Assign a 'replicate' column to the final_dataframe
    #only plate_type=="single_agent" would potentially have multiple replicates. "combination" plates will have replicate==1
    
    #Get the unique filenames for single-agent plates and assign them a replicate number
    single_agent_plates_df <- 
      data.frame("filename" = 
                   final_dataframe %>%
                   filter(plate_type=="single_agent") %>%
                   pull(filename) %>%
                   unique()) %>%
      mutate("replicate"=rownames(.))
    
    #Join the replicates (which are just for single-agent plates) to the final_dataframe
    final_dataframe <- 
      final_dataframe %>%
      left_join(single_agent_plates_df, by="filename") 
    
    #If plate_type=="combination" make replicate a value of 1:
    final_dataframe$replicate[final_dataframe$plate_type=="combination"] <- 1
    #====================================
    
    
  }#End handling of sparse mode data
  
  
  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Print transfer fail/success rates message
  num_failed_transfers <- sum(final_dataframe$failed_transfer)
  num_total_transfers  <- nrow(final_dataframe)
  
  print(paste0(num_failed_transfers, " of ", num_total_transfers, " transfers failed."))
  print(paste0("Failure rate: ", round(num_failed_transfers/num_total_transfers*100, 2), "% --- ",
  "Success rate: ", round(100 - num_failed_transfers/num_total_transfers*100, 2), "%"))
  print("(ONLY considers drug transfers, no controls or backfills)")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Return final_dataframe 
  print("Generation of metadata complete. It is recommended to spot check result for accuracy!")
  return(final_dataframe)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
}#End function
#--------------------------------------------------------------