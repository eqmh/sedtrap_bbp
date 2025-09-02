# PURPOSE:
#   This function parses an UW/MBARI APEX biochemical float *.msg file and
#   returns a structure of raw data. The # of columns are determined from
#   the low resolution header line
#
# USAGE:
#	data = parse_APEXmsg(file_name)
#
# INPUTS:
#	msg_name  = string of file name or  path\file name
#
# OUTPUTS:
#       data =  a structure of the data and headers. Size varies depending
#               on data flag.
#           data.FwRev  = Firmware revision number
#           data.CpActivationP =  Activation pressure (dbar) for cp?
#           data.FlbbMode = Flag in newer firmware 1=Yes,0=NO,[] no exist
#           data.lr_hdr = low res header cell array
#           data.lr_d   = low res data matrix
#           data.hr_hdr = high res header cell array
#           data.hr_d   = high res data matrix
#           data.pk_hdr = park sample header cell array               
#           data.pk_d   = park sample data matrix                            
#           data.cast   = profile number
#           data.sdn    = profile termination date [matlab sdn]
#           data.ice_flag = UnderIceEvasion y/n
#           data.ice_hex = UnderIceEvasion status for the previous 8 cycles
#           data.gps    = gps location fix [lon, lat]
#			data.irid   = iridium location fix (newer floats)
#           data.air    = Air oxygen measurements [Temp phase (Rphase)]
#           data.EOT    = flag count of EOT occurances
# EXAMPLES:
#   parse_APEXmsg4ARGO('c:\temp\7601.003.msg')
#   parse_APEXmsg4ARGO('7601.003.msg')
#
# CHANGES
#   01/18/2017 - cp header line search more generic & extracted CTD
#               type and serial number ultimately for ODV meta info
#   01/23/2017 - if no PTS remove line
#   04/06/2017 - add code to catch bad cp data lines, "00000000000000[2]",
#       in middle of cp section of msg file ~line 230
#   05/16/2017 - add code to look for <EOT> and record # of instances
#
#   09/28/2020 - Code was skipping over first line of optode air-cal
#   sequence.  Added a fix at the end of cp-mode hex extraction. (may be
#   revised at a later date).
#
#   10/26/20 - TM, Modified gps fix extraction to carry over datetime. Was
#              not included in original code.
#   10/28/20 - JP, Modified Surface Obs parseing code to include PAL sensor
#              float formats
#   04/28/21 - JP fixed 12712 issue processing cyles 61 & 111 with no
#              profile terminted time line in msg file. Get time from
#              SBE cp header line
#   05/13/21 - TM added check for hi-res pressure level of 0; exists in msg
#             file 12892.064.msg but bad data, should be removed (AOML screening for this)
#   06/03/21 - TM added code in support of SBE83 optode on APEX test float (air-cal format spec)
#
#   02/07/22 - EC added parser section for ParkPt and ParkPtFLBB
#   02/19/22 - EC added IceEvasionRecord in 8-digit and current status flag
#   04/27/22 - TM; incorporated the iridium fix positioning into parser (newer floats).
#   02/09/25 - MC: Fixed the issue of the GPS lon/lat position and time wrongly extracted

#----------------
# FOR TESTING
# msg_file = '\\atlas\chemwebdata\floats\alternate\f9095\9095.011.msg';

# Full PkPt
# msg_file = 'C:\Users\eclark\Documents\Matlab\ARGO_PROCESSING\MFILES\EKC commented tests\12878.054.msg'

# Full ParkPtFLBB
# wmo 5906471
# msg_file = 'C:\Users\eclark\Documents\Matlab\ARGO_PROCESSING\MFILES\EKC commented tests\20532.003.msg'
#msg_file = 'c:\temp\9031.003.msg'; # TESTING

# 10/28/2020 JP testing
# msg_file = 'c:\temp\12652.036.msg';
#----------------




# ************************************************************************
# FORMATS & VARIABLES
# ************************************************************************
# HIGH RESOLUTION SAMPLES FOR UW / MBARI APEX FLOATS
# The first 4-bytes of the encoded sample represents the pressure in
# centibars.  The second 4-bytes represents the temperature in
# millidegrees.  The third 4-bytes represent the salinity in parts per
# million.  The final 2-bytes represent the number of samples collected in
# the 2dbar pressure bin. 2's compliment

parse_APEXmsg4ARGO<-function(msg_file){             # Define the main parser function taking a path to a .msg file
  # PREDIMENSION OUTPUT STRUCTURE
  data <- list(                                      # Initialize the output list with all expected fields
    pk_hdr = NULL,  # park data header cell array     # Placeholder for ParkPt header (APEX-style)
    pk_d = NULL,    # park data matrix                # Placeholder for ParkPt data rows
    lr_hdr = NULL,  # low res header cell array       # Placeholder for low-resolution header tokens
    lr_d = NULL,    # low res data matrix             # Placeholder for low-resolution data matrix
    hr_hdr = NULL,  # high res header                 # Placeholder for high-resolution header tokens
    hr_d = NULL,    # high res data                   # Placeholder for high-resolution (CP-mode) data
    cast = NA,      # profile #                       # Profile (cast) number inferred from filename
    sdn = NA,       # profile termination time        # Profile end time (POSIXct)
    sdni = NA,      # iridium timestamp               # Iridium fix time (POSIXct) if present
    gps = NULL,     # gps location fix                # GPS fixes: time, lon, lat
    irid = NULL,    # iridium location fix (only some floats) # Iridium fixes: time, lon, lat
    air = NULL,     # Air oxygen measurements         # Surface obs (temp/phase etc.)
    aircal = NULL,  # air calibration measurements    # Optode air-calibration entries
    FwRev = NULL,                                   # Firmware revision code
    CpActivationP = NULL,                           # Constant profiling activation pressure (if reported)
    FlbbMode = NA,                                  # FLBB mode flag if present
    CTDtype = '',                                   # Parsed CTD type from SBE header
    CTDsn = '',                                     # Parsed CTD serial number from SBE header
    EOT = 0,        # complete message file flag      # Count of <EOT> markers (one per comms cycle)
    ice_flag = NULL,                                 # Under-ice evasion active now (0/1)
    ice_hex = NULL                                   # 8-bit history of ice evasion status
  )
  
  f_aircal <- '%*s%s%s%s%s%*f%f%f%f%f%f'            # Legacy (unused now) scanf-like format for OptodeAirCal
  f83_aircal <- '%*s%s%s%s%s%*f%f%f%f%f'             # Legacy format string for SBE83 air cal lines
  
  # Open file
  fid <- file(msg_file, "r")                         # Open the message file for reading
  
  tline <- ""                                         # Initialize current line buffer
  # The header line should start with 'p t s'
  # Find out if the file is complete / whether or not to proceed
  while (!grepl("p\\s+t\\s+s\\s+", tline)) {      # Loop until low-res header tokens (p t s ...) are found
    #print(tline)
    tline <- readLines(fid, n = 1)                    # Read next line
    if (length(tline)==0) {                           # If EOF before header -> incomplete file
      message("Incomplete message file! No profile data for ", msg_file)
      close(fid)                                      # Close handle
      return(NULL)                                    # Abort parse
    }
  }
  close(fid)                                          # Found header in this line; close handle for fresh pass
  
  # Find the header line starting with '$ p t s'
  ind1 <- regexpr("\\$\\s+p\\s+t\\s+s", tline)     # Locate the '$ p t s' marker within the header line
  
  # Extract headers into a character vector
  msg_hdr <- strsplit(substring(tline, ind1 + 1), "\\s+")[[1]]  # Tokenize the header fields after '$'
  float_vars <- msg_hdr                                 # Copy header tokens to float_vars
  float_vars <- float_vars[float_vars != ""]           # Remove empty tokens
  
  # BUILD LOW RES DATA FORMAT STRING - VARIES W/ FLOAT
  #low_res_format <- paste(rep("%f", length(float_vars)), collapse = " ") # (Not used: reading via scan)
  
  # Clean up variables
  rm(ind1, tline, fid, msg_hdr)                        # Free temporary variables
  
  # ************************************************************************
  # ************************************************************************
  #                         PARSE MESSAGE FILE
  # ************************************************************************
  # ORDER = header, parkPt, termination time, cp data, gps fix, surf obs
  # ************************************************************************
  fid <- file(msg_file, "r")                           # Reopen file to iterate from the start
  tline <- readLines(fid, n = 1)                        # Prime the read buffer
  
  # FILE EXISTS BUT NO DATA OR NON-STANDARD FILE FORMAT & FILE ENDED
  if (length(tline) == 0) {                             # Empty file guard
    message("File exists but empty inside - moving to next message file.")
    close(fid)
    data <- list()                                      # Return empty list to signal no data
    return(data)
  }
  
  # CAST # FROM MESSAGE FILE NAME
  str <- regmatches(msg_file, regexpr("\\d{3}(?=\\.msg)", msg_file, perl = TRUE)) # Extract 3-digit cast number
  data <- list(cast = as.numeric(str))                  # Initialize return object overriding with cast only
  
  # Clean up
  rm(str)
  
  # ************************************************************************
  # ************************************************************************
  # PARSING MSG FILE
  # ************************************************************************
  # ************************************************************************
  
  # Create checkpoint variables
  data_chk <- 0                                         # State of parsing (0: before low-res data, 1: got term time, 2: in cp)
  alt_sdn_chk <- 0                                      # Flag to indicate we must infer sdn from SBE header
  msg_task <- "profile time"                             # Current section we are looking for
  
  # Predimension variables
  low_res <- NULL                                       # Accumulator for low-resolution discrete samples
  high_res <- NULL                                      # Accumulator for CP-mode hex samples
  CpActP <- NULL                                        # Parsed activation pressure
  FwRev <- NULL                                         # Parsed firmware revision
  FlbbMode <- NA                                        # Parsed FLBB mode flag
  pk_data <- NULL                                       # Accumulator for ParkPt data
  ice_flag <- 0                                         # Under-ice now flag (set later via ice_hex LSB)
  ice_hex <- NULL                                       # 8-bit history of ice evasion
  
  # Step through lines to parse data
  while (length(tline) > 0) {                           # Iterate until EOF
    #print(tline)
    tline <- trimws(tline)                              # Trim whitespace for robust matching
    
    # GET SOME INFO FOR ANNIE
    # Use regexpr to find the position of "FwRev"
    FwRev_ind <- regexpr("FwRev", tline)               # Look for firmware revision token
    # Check if "FwRev" was found
    if (FwRev_ind != -1) {                              # If present on this line
      # Extract the firmware version using sub based on the position of "FwRev"
      FwRev <- as.numeric(sub(".*FwRev\\s+(\\d+):.*", "\\1", tline))  # Grab integer after 'FwRev'
      
      # Check if FwRev is valid
      if (is.na(FwRev)) {
        #cat("No firmware version found.\n")           # Silent if parse failed
      } else {
        #cat("Firmware version extracted:", FwRev, "\n") # Optional debug
      }
    }
    
    
    CpAct_ind <- regexpr("\\$ CpAct", tline)           # Look for CpAct (activation pressure) line
    if (CpAct_ind != -1) { # constant profiling activation P
      CpActP <- as.numeric(regmatches(tline, regexpr("\\d+", tline, perl = TRUE)))  # First integer on line
    }
    
    FlbbMode_ind <- regexpr("\\$ FlbbMode", tline)      # Look for FLBB mode flag line
    if (FlbbMode_ind != -1) { # flag may exist stating whether flbb on
      FlbbMode <- as.numeric(regmatches(tline, regexpr("\\d+", tline, perl = TRUE))) # Parse integer
    }
    
    
    # --------------------------------------------------------------------
    
    # GET PARK DATA
    # PARKPT
    if (grepl("^ParkPt:", tline)) {                    # APEX ParkPt line (pressure, temp at park)
      
      data$pk_d <- strsplit(tline, "\\s+")[[1]]       # Split the line by whitespace into tokens
      
      
      if (length(data$pk_d) == 9) {                     # Expected token count for ParkPt format
        
        data$pk_hdr <- c("Date", "p", "t")           # Define ParkPt header columns
        pkdat <- data$pk_d                               # Alias tokens
        
        # Collect date chunks
        d_str <- paste(pkdat[2], pkdat[3], pkdat[4], pkdat[5])  # Assemble "Mon dd YYYY HH:MM:SS"
        
        # Attempt to convert date
        sdn <- tryCatch(                                 # Parse date-time to POSIXct (UTC)
          as.POSIXct(d_str, format = "%b %d %Y %H:%M:%S", tz = "UTC"),
          error = function(e) NA
        )
        
        # Convert strings to numeric values
        Pk_p <- as.numeric(pkdat[8])                     # Park pressure
        Pk_t <- as.numeric(pkdat[9])                     # Park temperature
        pk_keep <- c(Pk_p, Pk_t)                         # Vector of keepers
        
        # Append to pk_data
        pk_data <- rbind(pk_data, c(sdn, pk_keep))       # Accumulate ParkPt row
        
      } else {
        message("Uncharacteristic number of ParkPt columns, check for incomplete msg file")
      }
      
      # PARKPT FLBB
    } else if (grepl("^ParkPtFlbb:", tline)) {          # ParkPt line including FLBB signals
      
      data$pk_d <- strsplit(tline, "\\s+")[[1]]       # Tokenize
      
      if (length(data$pk_d) == 12) {                     # Expected token count for ParkPtFlbb format
        
        data$pk_hdr <- c("Date", "p", "t", "Fsig", "Bbsig", "Tsig") # Define header incl. optical signals
        pkdat <- data$pk_d                               # Alias
        
        # Collect date chunks
        d_str <- paste(pkdat[2], pkdat[3], pkdat[4], pkdat[5])  # Date-time string
        sdn <- as.POSIXct(d_str, format = "%b %d %Y %H:%M:%S", tz = "UTC") # Parse time
        
        # Convert strings to numeric values
        Pk_p <- as.numeric(pkdat[8])                     # Pressure
        Pk_t <- as.numeric(pkdat[9])                     # Temperature
        Pk_Fs <- as.numeric(pkdat[10])                   # FLBB fluorescence signal
        Pk_Bbs <- as.numeric(pkdat[11])                  # FLBB backscatter signal
        Pk_Ts <- as.numeric(pkdat[12])                   # FLBB temperature signal
        
        pk_keep <- c(Pk_p, Pk_t, Pk_Fs, Pk_Bbs, Pk_Ts)   # Assemble output fields
        
        # Append to pk_data
        pk_data <- rbind(pk_data, c(sdn, pk_keep))       # Accumulate row
        
      } else {
        message("Uncharacteristic number of ParkPt columns, check for incomplete msg file")
      }
    }
    
    if (data_chk == 0 && grepl("^\\$ Profile", tline)) {   # Look for profile termination time line
      
      # Define the format for extracting date and time
      msg_time_format <- "%*s %*s %*s %*s %*s %s %s %s %s" # (Unused; kept for readability)
      
      # Extract date components
      d_str <- strsplit(tline, "\\s+")[[1]][6:9]       # Grab tokens with day, month, time, year
      
      # Arrange in a date-time string format similar to MATLAB
      s1 <- paste(d_str[2], "-", d_str[1], "-", d_str[4], " ", d_str[3], sep = "") # dd-Mon-YYYY HH:MM:SS
      
      # Convert to POSIXct to store as profile end time
      data$sdn <- as.POSIXct(s1, format = "%d-%b-%Y %H:%M:%S", tz = "UTC")  # Save termination time
      
      data_chk <- 1                                        # Transition state: termination time found
      rm(d_str, s1)
      
    } else if (data_chk == 0 && grepl("^\\$ Discrete samples", tline)) { # No termination time present
      
      message("WARNING: Profile terminated line not found for ", msg_file)
      message("Attempt to estimate termination time from SBE CP header.")
      data_chk <- 1                                        # Proceed as if termination time found
      alt_sdn_chk <- 1                                     # Mark that we must infer later
      
    } else if (data_chk == 1) {                            # After termination time, look for first data row
      
      ind1 <- regexpr("^(\\d+\\.\\d+)|^(-\\d+\\.\\d+)", tline)  # Numeric line indicates low-res sample row
      if (ind1 != -1) {
        msg_task <- "profile data"                          # Switch to low-res data extraction
        break                                               # Exit header scanning loop to start reading data
      }
    }
    
    tline <- readLines(fid, n = 1)                         # Read next line
  }
  
  # EXTRACT PROFILE DATA ----
  
  while (length(tline) > 0) {                               # Iterate over the remainder of the file
    #print(tline)
    tline <- trimws(tline)                                  # Trim for robust matching
    
    
    # ICE INFO        
    IceEvas_ind <- regexpr("IceEvasionRecord", tline)      # Look for hex-encoded ice evasion history
    if (IceEvas_ind != -1) {
      
      # Old, did not worl: ihex <- regmatches(tline, regexpr("(?<=IceEvasionRecord.+)\\w+", tline))
      ihex <- str_match(tline, "IceEvasionRecord.*?(\\w+)$")[, 2]  # Extract trailing hex string
      
      if (length(ihex) > 0) {
        ice_hex <- as.integer(intToBits(strtoi(ihex, base = 16))) # Convert hex to bits (LSB-first vector)
        lastnum <- ice_hex[8]  # Extract last digit of hex (8th position)
        
        if (lastnum == 1) {
          ice_flag <- 1  # Change flag to 1 if under ice now
        }
      }
      
      rm(IceEvas_ind)
    }
    
    # GET SOME INFO FOR ANNIE - Older float, found in footer
    FwRev_ind <- regexpr("FwRev", tline)                   # Firmware in footer for some floats
    if (FwRev_ind != -1) { # Firmware version
      FwRev <- as.numeric(regmatches(tline, regexpr("(?<=FwRev=)\\d+", tline, perl = TRUE))) # Parse digits after 'FwRev='
    }
    
    switch(msg_task,
           
           'profile data'= { # EXTRACT DISCRETE SAMPLE DATA
             if (regexpr("^#", tline) == -1) { # # means end of low res
               ind1 <- regexpr("\\(Park Sample\\)$", tline)      # Skip Park Sample banner lines inside data block
               if (ind1 == -1 && nchar(tline) > 0) {
                 tline <- gsub("nan", "NA", tline)                 # Normalize 'nan' to R's NA
                 tmp <- scan(text = tline, what = numeric(), quiet = TRUE) # Parse numeric row (p t s ...)
                 low_res <- rbind(low_res, tmp)                      # Append to low-res matrix
               }
             } else if (regexpr("^#.+sbe", tline, ignore.case = TRUE) != -1) { # cp hdr
               if (alt_sdn_chk == 1) { # 04/28/2021 JP no terminated time, get from SBE hdr
                 str <- regmatches(tline, regexpr("(?<=#\\s+)[\\w\\s:]+(?=\\s+Sbe)", tline)) # Extract date-time before 'Sbe'
                 if (length(str) > 0) {
                   data$sdn <- as.POSIXct(str, format = "%b %d %Y %H:%M:%S", tz = "UTC") # Use SBE header time as sdn
                 }
               }
               # MC: modification to read lines
               # Install stringr if not already installed
               if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
               
               # Load the package
               library(stringr)
               
               # Revised code snippet using str_extract from stringr
               CTDtype <- str_extract(tline, "(?i)sbe\\w+(?=serno)")      # Parse CTD model (e.g., sbe41/61/83)
               CTDsn <- str_extract(tline, "(?i)(?<=serno\\[)\\d+(?=\\])") # Parse serial number within brackets
               # data$CTDtype <- regmatches(tline, regexpr("sbe\\w+(?=serno)", tline, ignore.case = TRUE))
               # data$CTDsn <- regmatches(tline, regexpr("(?<=serno\\s)\\d+", tline))
               msg_task <- 'cp_data'  # finished low res, start high res
             } else { # No high resolution data found - look for a fix
               msg_task <- 'GPS data'
             }
           },
           
           'cp_data'= { # EXTRACT CP SAMPLE DATA
             # find cp data line & make sure hex characters only
             if (nchar(tline) == 14) {  # Check if the length of tline is 14 (right number of chars)
               
               high_res_format <- "%04x%04x%04x%02x"  # Format for hexadecimal input (reference only)
               ind1 <- regexpr("[A-Z0-9]{14}", tline)  # Check for 14 hex characters
               
               if (ind1 != -1) {  # If ind1 is not -1, hex characters were found
                 # Convert the hex string to numeric using strtoi
                 tmp <- c(strtoi(substr(tline, 1, 4), base = 16),
                          strtoi(substr(tline, 5, 8), base = 16),
                          strtoi(substr(tline, 9, 12), base = 16),
                          strtoi(substr(tline, 13, 14), base = 16))
                 
                 high_res <- rbind(high_res, tmp)  # Append tmp to high_res
                 
                 if (data_chk == 1) {
                   data_chk <- 2  # Indicate that we're now in the cp data lines
                 }
               }
             } else if (data_chk == 2 && regexpr("0{14}\\[\\d+\\]", tline) != -1) {
               # Comment out this else if/display when running through fleet to save time/power
               # message(paste('No data hex line in middle of cp data:', tline, 'for', msg_file))
             } else if (data_chk == 2 && nchar(tline) != 14) { # end of cp
               # TM 9/28/20 Code was skipping first line of aircal data.
               # Add OptodeAirCal block here after end of cpmode extraction.
               if (regexpr("^OptodeAirCal", tline) != -1) {
                 # MC: 30/10/2024 : modify to extract data
                 split_tline <- unlist(strsplit(tline, "\\s+"))           # Tokenize air-cal line
                 
                 # Extract relevant components based on expected positions
                 month <- split_tline[2]
                 day <- split_tline[3]
                 year <- split_tline[4]
                 time <- split_tline[5]  # This should be treated as a character
                 
                 # Convert remaining elements to numeric, keeping the first element as timestamp
                 ac_tmp <- as.numeric(split_tline[c(6, 7, 8, 9, 10, 11)])  # Numeric portion
                 
                 # Combine date parts to form the datetime string
                 d_str <- paste(month, day, year, time)
                 
                 # Check if the second element from ac_tmp is valid before appending to data
                 if (!is.na(ac_tmp[1]) && nchar(as.character(ac_tmp[1])) > 0) {
                   # Convert d_str to POSIXct datetime
                   d_str<-as.POSIXct(d_str, format = "%b %d %Y %H:%M:%S", tz = "UTC")
                   
                   # Use ac_tmp[1] as the timestamp (1692487525) instead of the numeric value (113)
                   data$aircal <- rbind(data$aircal, 
                                        c(d_str, ac_tmp[2]))  # Store datetime and next numeric value
                 }
               }
               if (regexpr("^Sbe83AirCal", tline) != -1) {
                 ac_tmp <- as.data.frame(matrix(scan(text = tline, what = character(), nmax = length(f83_aircal), sep = ""), ncol = length(f83_aircal), byrow = TRUE))
                 d_str <- paste(ac_tmp[1, 1], ac_tmp[1, 2], ac_tmp[1, 3], ac_tmp[1, 4])
                 if (!is.na(ac_tmp[1, 2]) && nchar(ac_tmp[1, 2]) > 0) {
                   data$aircal <- rbind(data$aircal, c(as.POSIXct(d_str, format = "%b %d %Y %H:%M:%S", tz = "UTC"), as.numeric(ac_tmp[1, 2])))
                 }
               }
               msg_task <- 'GPS data';                            # Transition to GPS/surface obs block
             }
           },
           
           'GPS data'= { # GET GPS POSITION & AIR OXYGEN VALUES
             if (grepl("^Fix:", tline)) {                        # GPS fix line
               gps <- as.numeric(scan(text = tline, what = character(), quiet = TRUE)[2:3])   # Read lon, lat
               if (!is.na(gps[1]) && !is.na(gps[2])) {
                 tokens <- scan(text = tline, what = character(), quiet = TRUE)
                 gps_sdn <- paste(tokens[4], tokens[5]) # [mm/dd/yyyy; hhmmss]
                 if (nchar(gps_sdn) == 17) {
                   Gsdn <- as.POSIXct(gps_sdn, format = "%m/%d/%Y %H%M%S", tz = "UTC")   # Parse time
                   data$gps <- rbind(data$gps, c(Gsdn, gps))                                  # Append row
                 } else {
                   data$gps <- rbind(data$gps, c(data$sdn, NA, NA))                           # Fallback to termination time
                 }
               } else {
                 data$gps <- rbind(data$gps, c(data$sdn, NA, NA))                             # Missing values -> NA
               }
             }
             
             if (grepl("^IridiumFix:", tline)) {                 # Iridium fix line (newer floats)
               irid <- as.numeric(scan(text = tline, what = character(), quiet = TRUE)[3:4])   # Read lon, lat
               if (!is.na(irid[1]) && !is.na(irid[2])) {
                 tokensi <- scan(text = tline, what = character(), quiet = TRUE)
                 gps_sdni <- paste(tokensi[7], tokensi[8]) # [mm/dd/yyyy; hh:mm:ss]
                 if (nchar(gps_sdni) == 19) {
                   Isdn <- as.POSIXct(gps_sdni, format = "%m/%d/%Y %H:%M:%S", tz = "UTC")  # Parse time
                   data$irid <- rbind(data$irid, c(Isdn, irid))                                # Append row
                 } else {
                   data$irid <- rbind(data$irid, c(data$sdni, NA, NA))                         # Fallback
                 }
               } else {
                 data$irid <- rbind(data$irid, c(data$sdni, NA, NA))                           # Missing -> NA
               }
             }
             
             # complete surf line is bounded by curly braces, count to check
             if (grepl("^(SurfaceObs)", tline) && length(gregexpr("[{}]", tline)[[1]]) == 2) {
               # use slash as delimiter to break up string
               stmp <- unlist(strsplit(tline, "/"))             # Split into segments
               if (length(stmp) == 4) { # count substrings
                 O2_str <- trimws(stmp[3])                        # Usual position of O2 values
               } else if (length(stmp) == 2) { # PAL SENSOR FLOATS!
                 O2_str <- trimws(stmp[2])
                 O2_str <- gsub("}", "", O2_str)                 # Remove closing brace
               } else {
                 message('SurfaceObs line format not recognized. Could not parse line')
                 O2_str <- ''
               }
               
               # space character divides numbers, number count = space count + 1
               num_ct <- sum(strsplit(O2_str, '')[[1]] == ' ') + 1  # Count numbers present
               tf_endC <- grepl("C$", O2_str)                    # Temperature in Celsius at end?
               if (num_ct == 3 && tf_endC) { # 3 #'s & "C" at end of string
                 tmp <- as.numeric(unlist(strsplit(gsub("C$", "", O2_str), " ")))  # Parse numbers
                 tmp <- tmp[c(3, 1, 2)] # T TPh RPh                             # Reorder to T, TempPhase, RPhase
                 
               } else if (num_ct == 3 && grepl("\\s[\\d\\.]+C\\s", O2_str)) {
                 tmp <- as.numeric(unlist(strsplit(gsub("C$", "", O2_str), " ")))
               } else if (num_ct == 3) { # 3 #'s & "C" must be at beginning
                 tmp <- as.numeric(unlist(strsplit(gsub("C$", "", O2_str), " ")))
               } else if (num_ct == 2 && tf_endC) { # 2 #'s & "C" at end of string
                 tmp <- as.numeric(unlist(strsplit(gsub("C$", "", O2_str), " ")))
                 tmp <- tmp[c(2, 1)] # temp, phase                            # Reorder to Temp, Phase
               } else if (num_ct == 2) {
                 tmp <- as.numeric(unlist(strsplit(gsub("C$", "", O2_str), " ")))
               } else {
                 message('Could not parse surf obs')
               }
               data$air <- rbind(data$air, tmp)                 # Accumulate surface obs values
             }
             
             if (grepl("^OptodeAirCal", tline)) {              # OptodeAirCal line outside CP block
               # MC: 30/10/2024 : modify to extract data
               split_tline <- unlist(strsplit(tline, "\\s+"))  # Tokenize
               
               # Extract relevant components based on expected positions
               month <- split_tline[2]
               day <- split_tline[3]
               year <- split_tline[4]
               time <- split_tline[5]  # This should be treated as a character
               
               # Convert remaining elements to numeric, keeping the first element as timestamp
               ac_tmp <- as.numeric(split_tline[c(6, 7, 8, 9, 10, 11)])  # Numeric portion
               
               # Combine date parts to form the datetime string
               d_str <- paste(month, day, year, time)
               
               # Check if the second element from ac_tmp is valid before appending to data
               if (!is.na(ac_tmp[1]) && nchar(as.character(ac_tmp[1])) > 0) {
                 # Convert d_str to POSIXct datetime
                 d_str<-as.POSIXct(d_str, format = "%b %d %Y %H:%M:%S", tz = "UTC")
                 
                 # Use ac_tmp[1] as the timestamp (1692487525) instead of the numeric value (113)
                 data$aircal <- rbind(data$aircal, 
                                      c(d_str, ac_tmp[2]))        # Store time + value
               }
             }
             
             if (grepl("^Sbe83AirCal", tline)) {               # SBE83-specific air-cal line
               ac_tmp <- as.data.frame(matrix(scan(text = tline, what = character(), nmax = length(f83_aircal), sep = ""), ncol = length(f83_aircal), byrow = TRUE))
               d_str <- paste(ac_tmp[1, 1], ac_tmp[1, 2], ac_tmp[1, 3], ac_tmp[1, 4])
               if (!is.na(ac_tmp[1, 2])) {
                 data$aircal <- rbind(data$aircal, c(as.POSIXct(d_str, format = "%b %d %Y %H:%M:%S", tz = "UTC"), as.numeric(ac_tmp[1, 2])))
               }
             }
             
             # check for file termination, one per satellite comms/GPS cycle
             if (grepl("^<EOT>", tline)) {                      # End-of-transmission marker
               data$EOT <- data$EOT + 1;                         # Increment EOT count
             }
           }
    )
    tline <- readLines(fid, n = 1)                        # Read next line for this loop
  }
  close(fid)                                              # Close file after main parse
  
  # Check if GPS data is empty and initialize
  if (is.null(data$gps) || nrow(data$gps) == 0) {         # Ensure gps field exists even if missing
    data$gps <- data.frame(sdn = data$sdn, lon = NaN, lat = NaN)
  }
  
  # ************************************************************************
  # CONVERT CP TO USEFUL NUMBERS
  # ************************************************************************
  if (!is.null(high_res) && nrow(high_res) > 0) {         # If CP hex samples were collected
    # Conversion constants
    hex_conv <- matrix(c(32768, 10,                        # Pressure: 2's complement cutoff & scale
                         61440, 1000,                      # Temperature: cutoff & scale
                         61440, 1000), ncol = 2, byrow = TRUE) # Salinity: cutoff & scale
    
    tmp <- matrix(NaN, nrow = nrow(high_res), ncol = 3)   # Working matrix for sign flags
    t0 <- matrix(1, nrow = nrow(tmp), ncol = 1)           # Column of ones for broadcasting
    
    # Build test matrix: > CUTOFF, NEG #'s, = 1; < CUTOFF, POS #'s, = 0; = CUTOFF, BAD VALUE, = NaN
    t_hi <- (high_res[, 1:3] - t0 %*% t(hex_conv[, 1]) > 0)  # Identify negative 2's-complement values
    tmp[t_hi] <- 1                                        # Mark negatives
    t_low <- (high_res[, 1:3] - t0 %*% t(hex_conv[, 1]) < 0) # Identify positives
    tmp[t_low] <- 0                                       # Mark positives
    
    # Convert hex values
    high_res[, 1:3] <- high_res[, 1:3] - (t0 %*% matrix(c(65536, 65536, 65536), ncol = 3)) * tmp # Apply 2's-comp sign
    high_res[, 1:3] <- high_res[, 1:3] / (t0 %*% t(hex_conv[, 2]))                                # Scale to engineering units
  }
  
  # Clear temporary variables (R automatically cleans up unreferenced variables)
  # ************************************************************************
  # FILL OUTPUT STRUCTURE
  # ************************************************************************
  if (is.null(low_res) && is.null(high_res)) {            # No data detected in either block
    message(sprintf("No data found in %s", msg_file))
  } else {
    data$FwRev         <- FwRev                           # Store firmware rev
    data$CpActivationP <- CpActP                          # Store activation pressure
    data$FlbbMode      <- FlbbMode                        # Store FLBB flag
    data$ice_flag      <- ice_flag                        # Under-ice flag now
    data$ice_hex       <- ice_hex                         # 8-bit history
    
    if (!is.null(low_res) && nrow(low_res) > 0) {         # If low-res data collected
      tnan <- rowSums(is.na(low_res[, 1:3])) == 3         # Identify rows missing P, T, S entirely
      low_res <- low_res[!tnan, ]                         # Drop empty PTS rows
      data$lr_hdr <- float_vars                           # Save header tokens
      data$lr_d <- low_res                                # Save low-res matrix
    } else {
      message(sprintf("No low resolution data found for %s", msg_file))
    }
    
    if (!is.null(high_res) && nrow(high_res) > 0) {       # If CP-mode data present
      # Check for any zero pressure level - remove
      high_res <- high_res[high_res[, 1] != 0, ]          # Remove obviously bad P==0 rows
      data$hr_d  <- high_res[, 1:4]                       # Keep p, t, s, nbin
      data$hr_hdr <- c(float_vars[1:3], 'nbin ctd')       # Header for high-res block
    } else {
      message(sprintf("No high resolution data found for %s", msg_file))
    }
    
    # Uncomment for testing ParkPt data status
    # if (!is.null(pk_data) && nrow(pk_data) > 0) {
    data$pk_d <- pk_data                                  # Save ParkPt accumulator (may be NULL)
    # } else {
    #     message(sprintf("%s park data is empty", msg_file))
    # }
  }
  
  # Clear workspace of unnecessary variables (optional, R generally handles this automatically)
  # List all objects in the environment
  all_objects <- ls()                                     # Enumerate objects in current scope
  
  # Identify objects to remove, excluding 'data'
  objects_to_remove <- setdiff(all_objects, "data")      # Everything except the return payload
  
  # Remove selected objects
  rm(list = objects_to_remove)                            # Clean up memory within function scope
  
  return(data)                                            # Return parsed structure
  
}                                                         # End of parse_APEXmsg4ARGO()
