#!/usr/bin/env Rscript

# Script Name: install_pavianCoreTools_packages.R
# Purpose: To install required packages for pavianCoreTools.R

# Define the required packages for pavianCoreTools
standard_packages <- c("tools", "dplyr", "writexl", "readr", "purrr", "argparse", 
                       "utils", "plyr", "stats", "htmlwidgets", "webshot")

special_packages <- list(sankeyD3 = "fbreitwieser/sankeyD3")

# Function to install missing standard packages from CRAN
install_standard_packages <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(missing_packages) > 0) {
    message("Installing missing standard packages: ", paste(missing_packages, collapse=", "))
    install.packages(missing_packages, repos="https://cloud.r-project.org/")
  } else {
    message("All standard packages are already installed.")
  }
}

# Function to install special packages from sources like GitHub
install_special_packages <- function(packages) {
  for (pkg in names(packages)) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message("Installing special package: ", pkg, " from source: ", packages[[pkg]])
      if(pkg == "sankeyD3" && !require("devtools", character.only = TRUE, quietly = TRUE)) {
        install.packages("devtools", repos="https://cloud.r-project.org/")
      }
      devtools::install_github(packages[[pkg]])
    }
  }
}

# Setting up dependencies for webshot
setup_webshot_dependencies <- function() {
  if (require("webshot", character.only = TRUE, quietly = TRUE)) {
    message("Setting up additional dependencies for webshot.")
    phantom_path <- webshot:::find_phantom()
    if (is.null(phantom_path) || !file.exists(phantom_path)) {
      webshot::install_phantomjs()
    }
  }
}

# Check if pandoc is installed
check_pandoc <- function() {
  pandoc_version <- system("pandoc --version", intern = TRUE)
  
  # Check OS type
  os_type <- .Platform$OS.type

  # Conditionally check for pandoc based on OS type
  if (os_type == "windows" && Sys.info()["sysname"] == "Windows") {
    message("You're on a Windows system. RStudio typically comes with pandoc. Make sure you're using RStudio if you encounter pandoc-related issues.")
  } else if (os_type == "unix" && Sys.info()["sysname"] == "Linux") {
    if (length(pandoc_version) == 0 || !grepl("pandoc", pandoc_version[1])) {
      message("pandoc is not detected on your Linux system. It's essential for saving self-contained HTML widgets.")
      message("Please download and install it from the Pandoc official website: https://pandoc.org/installing.html")
      message("After installation, make sure it's added to your system PATH so R can find and use it.")
    } else {
      message("pandoc is detected: ", pandoc_version[1])
    }
  }
}


# Call the functions
install_standard_packages(standard_packages)
install_special_packages(special_packages)
setup_webshot_dependencies()
check_pandoc()
