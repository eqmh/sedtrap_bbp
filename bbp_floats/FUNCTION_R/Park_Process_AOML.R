# ************************************************************************
# PURPOSE:
# This function processes raw message files for a given AOML APEX float, 
# calculates useful values from the raw signals observed at park depth
# (CHl, Bbp, etc.) and merges these results for a given profile.
#
# INPUTS:
#   float_id = Float ID as a string
#   calibration info:
#     BetabDC = bbp dark count
#     BetabScale = bbp scale factor
#     ChlDC = chlorophyll dark count
#     ChlScale = chlorophyll scale factor
#
# OUTPUTS:
#   PARK = list containing park depth data and some metadata
# ************************************************************************

Park_Process_AOML<-function(floatid,
                            path_data,
                            path_save){
  #floatid <- "ua19097"
  
  # Enter calibration info for floats of interest
  if (floatid == "19073") {
    BetabDC <- 46
    BetabScale <- 1.608e-06
    ChlDC <- 47
    ChlScale <- 0.0073
  } else if (floatid == "19097") {
    BetabDC <- 50
    BetabScale <- 1.809e-06
    ChlDC <- 47
    ChlScale <- 0.0066
  } else if (floatid == "19821") {
    BetabDC <- 45
    BetabScale <- 1.706e-06
    ChlDC <- 47
    ChlScale <- 0.0073
  } else if (floatid == "1469") {
    BetabDC <- 48
    BetabScale <- 3.078e-07
    ChlDC <- 51
    ChlScale <- 0.001869
  }
  
  # Set directories
  floatdir <- paste0(path_data,floatid,"/")
  #msgdir <- paste0(floatdir,"msg/")
  msgdir <-floatdir
  
  # List msg files in directory
  msg_list <- list.files(path = msgdir, pattern = "*.msg", full.names = TRUE)
  
  # Placeholder for PARK data storage
  PARK <- list()
  
  # Process each msg file
  for (msg_file in msg_list) {
    print(msg_file)
    
    # Extract profile number from the filename
    cast_num <- sub(".*/[0-9]+\\.([0-9]+)\\.msg$", "\\1", msg_file)
    
    # Parse the message file (Placeholder for actual function)
    d <- parse_APEXmsg4ARGO(msg_file)
    
    # Fill values
    fv <- list(bio = 99999, QC = 99)
    
    # Valid ranges
    RCR <- list(P = c(0, 12000), 
                T = c(-2.5, 40), # from argo parameter list
                CHL = c(-0.1, 150), # argoBGC QC manual 09July2016
                BB700 = c(-0.000025, 0.1)) # argoBGC QC manual 09July2016
    RC <- list(CHL = c(-0.1, 50), 
               BB700 = c(-0.000025, 0.1)
               )
    
    # Assuming d.pk_d contains data (otherwise, set empty logic)
    if (!is.null(d$pk_d)) {
      
      # Saving data in PARK
      PARK$cast <- d$cast
      PARK$gps <- d$gps
      PARK$SDN <- d$sdn
      
      # Identify headers of interest
      iPT <- which(d$pk_hdr == "t")
      iPChl <- which(d$pk_hdr == "Fsig")
      iPkPres <- which(d$pk_hdr == "p")
      iPBb <- which(d$pk_hdr == "Bbsig")
      iPsdn <- which(d$pk_hdr == "Date")
      
      # Set up general vars
      PARK$PRES <- d$pk_d[, iPkPres]
      PARK$PRES_ADJUSTED <- PARK$PRES
      PARK$SDN <- d$pk_d[, iPsdn]
      
      pkfill0 <- rep(0, length(PARK$PRES))
      
      PARK$PRES_QC <- pkfill0 + fv$QC
      PARK$PRES_ADJUSTED_QC <- pkfill0 + fv$QC
      PARK$TEMP <- d$pk_d[, iPT]
      PARK$TEMP_ADJUSTED <- PARK$TEMP
      PARK$TEMP_QC <- pkfill0 + fv$QC
      PARK$TEMP_ADJUSTED_QC <- pkfill0 + fv$QC
      
      # Check for bad pressure values
      PKQF_P <- (PARK$PRES < RCR$P[1]) | (PARK$PRES > RCR$P[2])
      PARK$PRES_QC[PKQF_P] <- 4 #BAD
      PARK$PRES_QC[!PKQF_P & PARK$PRES != fv$bio] <- 1 # GOOD
      
      # Check for bad temp values
      PKQF_T <- (PARK$TEMP < RCR$T[1]) | (PARK$TEMP > RCR$T[2])
      PARK$TEMP_QC[PKQF_T] <- 4 #BAD
      PARK$TEMP_QC[!PKQF_T & PARK$TEMP != fv$bio] <- 1 #GOOD
      
      # Calculate CHL Concentration
      if (!is.null(iPChl)) {
        t_nan <- is.na(d$pk_d[, iPChl])
        PARK$FLUORESCENCE_CHLA <- pkfill0 + fv$bio
        PARK$FLUORESCENCE_CHLA_QC <- pkfill0 + fv$QC
        PARK$CHLA <- pkfill0 + fv$bio
        PARK$CHLA_QC <- pkfill0 + fv$QC
        PARK$CHLA_ADJUSTED <- pkfill0 + fv$bio
        PARK$CHLA_ADJUSTED_QC <- pkfill0 + fv$QC
        
        PARK$FLUORESCENCE_CHLA[!t_nan] <- d$pk_d[!t_nan, iPChl]
        PARK$FLUORESCENCE_CHLA_QC[!t_nan] <- fv$QC
        
        PARK$CHLA[!t_nan] <- (d$pk_d[!t_nan, iPChl] - ChlDC) * ChlScale
        PARK$CHLA_QC[!t_nan] <- 3 # QC flag 3 do not use w/o adjusting
        
        PARK$CHLA_ADJUSTED[!t_nan] <- (d$pk_d[!t_nan, iPChl] - ChlDC) * ChlScale/2
        PARK$CHLA_ADJUSTED_QC[!t_nan] <- 1 # process chl data
        
        PARK$CHLA_ADJUSTED_ERROR[!t_nan]<-abs(PARK$CHLA_ADJUSTED[!t_nan]*2)
        
        # Range check
        t_bio <- PARK$CHLA != fv$bio
        t_chk <- t_bio & (PARK$CHLA < RCR$CHL[1] | PARK$CHLA > RCR$CHL[2])
        PARK$CHLA_QC[t_chk] <- 4
        PARK$FLUORESCENCE_CHLA_QC[t_chk] <- 4
        
        t_bio <- PARK$CHLA_ADJUSTED != fv$bio
        t_chk <- t_bio & (PARK$CHLA_ADJUSTED < RCR$CHL[1] | PARK$CHLA_ADJUSTED > RCR$CHL[2])
        PARK$CHLA_ADJUSTED_QC[t_chk] <- 4
        PARK$FLUORESCENCE_CHLA_ADJUSTED_QC[t_chk] <- 4
      }
      
      # Calculate Particle Backscatter Coefficient
      if (!is.null(iPBb)) {
        t_nan <- is.na(d$pk_d[, iPBb])
        pVSF <- pkfill0 + fv$bio
        BETA_SW <- pkfill0 + fv$bio
        PARK$BETA_BACKSCATTERING700 <- pkfill0 + fv$bio
        PARK$BBP700 <- pkfill0 + fv$bio
        PARK$BBP700_ADJUSTED <- pkfill0 + fv$bio
        PARK$BETA_BACKSCATTERING700_QC <- pkfill0 + fv$QC
        PARK$BBP700_QC <- pkfill0 + fv$bQC
        PARK$BBP700_ADJUSTED_QC <- pkfill0 + fv$QC
        
        # Volume Scattering Function (VSF)
        pVSF[!t_nan] <- (d$pk_d[!t_nan, iPBb] - BetabDC) * BetabScale
        
        # Constants
        X <- 1.097 * 2 * pi           # FLBB, BBP processing doc, July 2016
        LAMBDA <- 700
        DELTA <- 0.039                # Depolarization ratio
        THETA <- 142                  # FLBBAP2, DAC manual says 120
        SALest <- 33.5                # Estimated salinity for BSW calculation
        
        message("CAUTION: Using estimated sal value of 33.5 in BSW calculation")
        
        # Identify seawater data (particles only)
        BETA_SW_ind <- which(t_nan == 0)
        
        if (length(BETA_SW_ind) > 0) {
          for (ct in seq_along(BETA_SW_ind)) {
            # Call the equivalent function to betasw_ZHH2009
            results <- betasw_ZHH2009(LAMBDA, d$pk_d[BETA_SW_ind[ct], iPT], THETA, SALest, DELTA)
            BETA_SW[BETA_SW_ind[ct]] <- results$betasw
            b90sw <- results$beta90sw
            bsw <- results$bsw
          }
        }
        
        # Populate PARK structure with calculated values
        PARK$BETA_BACKSCATTERING700[!t_nan] <- d$pk_d[!t_nan, iPBb]              # counts
        PARK$BETA_BACKSCATTERING700_QC[!t_nan] <- fv$QC                        # counts
        PARK$BBP700[!t_nan] <- (pVSF[!t_nan] - BETA_SW[!t_nan]) * X           # b_bp m^-1
        PARK$BBP700_QC[!t_nan] <- 2 # 3 do not use w/o adjusting ... 6/10/21 modify qcraw flag from 3 to 2.                                          # Update QC flag
        
        PARK$BBP700_ADJUSTED[!t_nan] <- PARK$BBP700[!t_nan]
        PARK$BBP700_ADJUSTED_QC[!t_nan] <- 1
        PARK$BBP700_ADJUSTED_ERROR[!t_nan] <- fv$bio                          # Placeholder
        
        # Final range check on values
        t_bio <- PARK$BBP700 != fv$bio                                        # Check if BBP data is present
        t_chk <- t_bio & (PARK$BBP700 < RCR$BB700[1] | PARK$BBP700 > RCR$BB700[2])
        PARK$BBP700_QC[t_chk] <- 4
        PARK$BETA_BACKSCATTERING700_QC[t_chk] <- 4
        
        t_bio <- PARK$BBP700_ADJUSTED != fv$bio                               # Check for adjusted BBP data
        t_chk <- t_bio & (PARK$BBP700_ADJUSTED < RC$BB700[1] | PARK$BBP700_ADJUSTED > RC$BB700[2])
        PARK$BBP700_ADJUSTED_QC[t_chk] <- 4
        
        # Clear temporary variables
        rm(BETA_SW, X, pVSF, ct, b90sw, bsw)
      }
    } else{print(c("No park data in message file for",msg_file))}
    # Save data (Placeholder)
    save(PARK, file = file.path(path_save, paste0("PARK_",floatid,"_",cast_num, ".RData")))
  }
}

