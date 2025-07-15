##########################################################################
# Script: PARK Data Processing and BBP Drift Plotting for BGC-Argo Floats
# Author: M. Cornec
# Date: 14 July 2025
#
# Description:
# This script processes BGC-Argo PARK data for specified WMO floats. 
# It extracts BBP700_ADJUSTED values from PARK messages, applies quality control, 
# and generates drift plots of BBP values over time, color-coded by depth.
#
# Notes:
# - Paths must be adapted by the user before running the script.
# - Required libraries and functions are automatically checked/loaded.
##########################################################################

# Clean workspace and clear console
rm(list = ls())
cat("\014")

# Define file paths (to be updated by user)
path <- "D:/WIP/AOML/GOM_PRO/For Github/bbp_floats/"
dir_function <- paste0(path, "FUNCTION_R/")
path_data <- paste0(path, "Tables/")
path_data_msg <- paste0(path_data, "msg/")
path_fig <- paste0(path, "Figs/")

# Required functions to source ----
funX <- c(
  "check_and_install.R",
  "betasw_ZHH2009.R",
  "Park_Process_AOML.R",
  "parse_APEXmsg4ARGO.R"
)

for (i in funX) {
  source(paste0(dir_function, i))
}

# Required libraries ----
pack <- c(
  "lubridate",  # Easy date-time handling
  "geosphere",
  "ggplot2",
  "gridExtra",
  "viridis",
  "cowplot",
  "scales"
)

for (i in pack) {
  check_and_install(i)
}

# Load float and corresponding message IDs ----
wmoz <- c(
  4903622,
  4903624,
  4903625
)

msg_wmoz <- c(
  19073,
  19097,
  19821
)

# Extract the PARK data for each float ----
original_locale <- Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "C")  # Ensure consistent month abbreviations

for (woko in wmoz) {
  print(woko)
  corr_msg <- msg_wmoz[which(wmoz == woko)]
  
  Park_Process_AOML(
    corr_msg,
    path_data_msg,
    paste0(path_data, "PARK/")
  )
}

Sys.setlocale("LC_TIME", original_locale)  # Restore original locale

# Process and plot BBP drift per float ----
for (woko in wmoz) {
  print(woko)
  corr_msg <- msg_wmoz[which(wmoz == woko)]
  
  dat_list <- list.files(paste0(path_data, "PARK/"))
  dat_list <- dat_list[grep(corr_msg, dat_list)]
  
  # Extract dataframe with QC = 1 only
  data_park <- NULL
  for (i in dat_list) {
    load(paste0(path_data, "PARK/", i))
    
    time <- as.POSIXct(PARK$SDN)
    depth <- PARK$PRES_ADJUSTED
    Pbbp <- PARK$BBP700_ADJUSTED
    Pbbp_qc <- PARK$BBP700_ADJUSTED_QC
    
    remove(PARK)
    
    data_park <- rbind(data_park, data.frame(cbind(time, depth, Pbbp, Pbbp_qc)))
  }
  
  data_park <- data_park[which(data_park$Pbbp_qc == 1), ]
  data_park$time <- as.Date(as.POSIXct(data_park$time))
  
  # Plot BBP drift if valid data exist ----
  if (all(is.na(data_park$Pbbp)) == FALSE) {
    fig_bbp <- ggplot() +
      geom_point(data = data_park, aes(x = time, y = Pbbp, color = depth), size = .7) +
      scale_y_continuous(
        limits = c(0, 0.01),
        name = expression(b[bp] ~ (m^{-1}))
      ) +
      scale_color_viridis_c(
        option = "turbo",
        direction = -1,
        name = expression(atop(depth, (m)))
      ) +
      scale_x_date(date_labels = "%b %Y", name = "Time") +
      theme_bw() +
      theme(
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text.y = element_text(size = 15)
      )
    
    ggsave(paste0(path_fig, woko, "_bbp_drift.png"), fig_bbp, width = 10, height = 5, dpi = 300)
  }
}
