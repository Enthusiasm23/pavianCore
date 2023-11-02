#!/usr/bin/env Rscript

################################################################################

suppressPackageStartupMessages(library(tools))       # 用于文件名操作
suppressPackageStartupMessages(library(dplyr))       # 用于数据处理
suppressPackageStartupMessages(library(writexl))     # 用于excel文件写入
suppressPackageStartupMessages(library(tidyverse))   # 用于tsv文件写入
suppressPackageStartupMessages(library(purrr))       # 用于处理列表和函数
suppressPackageStartupMessages(library(argparse))    # 用于处理命令行参数

################################################################################

print_full_description <- function() {
  ascii_art <- paste(
    "                    _              ____              _____           _     ",
    "  _ __   __ ___   _(_) __ _ _ __  / ___|___  _ __ __|_   _|__   ___ | |___ ",
    " | '_ \\ / _` \\ \\ / / |/ _` | '_ \\| |   / _ \\| '__/ _ \\| |/ _ \\ / _ \\| / __|",
    " | |_) | (_| |\\ V /| | (_| | | | | |__| (_) | | |  __/| | (_) | (_) | \\__ \\",
    " | .__/ \\__,_| \\_/ |_|\\__,_|_| |_|\\____\\___/|_|  \\___||_|\\___/ \\___/|_|___/",
    " |_|                                                                        ", 
    sep="\n"
  )
  
  cat(
    "----------------------------------------------------------------------------\n",
    ascii_art, "\n",
    "----------------------------------------------------------------------------\n",
    "Description:\n",
    "  This script is designed for analyzing, organizing, and summarizing\n",
    "  Kraken report data. It facilitates the visualization of the analysis\n",
    "  results. This script is a modification of the Pavian tool, an interactive\n",
    "  browser application for analyzing metagenomics classification results.\n",
    "  Original Pavian can be found at: https://github.com/fbreitwieser/pavian\n",
    "  This modified command-line version is available at:\n",
    "  https://github.com/Enthusiasm23/pavianCore\n",
    "\n",
    "Author: Enthusiasm23\n",
    "Completion Date: 2023-10-31\n",
    "\n",
    "Dependencies:\n",
    "  tools, dplyr, writexl, tidyverse, purrr, argparse.\n",
    "  For a full list of dependencies, see the pavianCore repository at:\n",
    "  https://github.com/Enthusiasm23/pavianCore\n",
    "\n",
    "Usage:\n",
    "  This script is meant to be run from the command line.\n",
    "  It accepts parameters to specify the input data and options for the\n",
    "  analysis and output formatting.\n",
    "\n",
    "License:\n",
    "  Please refer to the original Pavian license and the modifications at\n",
    "  the pavianCore repository.\n",
    "----------------------------------------------------------------------------\n",
    sep = ""
  )
}


dmessage <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), " ", ...)
}


read_report <- function(myfile, has_header=NULL, check_file = FALSE) {
  
  first.line <- tryCatch( readLines(myfile,n=1, warn=FALSE),
                          error = function(e) { warning("Error reading ",myfile); return() })
  isASCII <-  function(txt) {
    if (length(txt) == 0)
      return(FALSE)
    raw <- charToRaw(txt)
    all(raw <= as.raw(127) & (raw >= as.raw(32) | raw == as.raw(9)))
  }
  if (length(first.line) == 0) { 
    dmessage("Could not read ", myfile, ".")
    return(NULL)
  }
  tryCatch({
    if (nchar(first.line) == 0) {
      dmessage("First line of ", myfile, " is empty")
      return(NULL)
    }
  }, error = function(e) { 
    dmessage(e)
    return(NULL)
  })
  
  if (!isTRUE(isASCII(first.line))) {
    dmessage(myfile," is not a ASCII file")
    return(NULL)
  }
  
  if (is.null(has_header)) {
    has_header <- grepl("^[a-zA-Z#%\"]",first.line)
  }
  is_metaphlan3_fmt <- grepl("^#mpa_v", first.line)
  is_metaphlan2_fmt <- grepl("Metaphlan2_Analysis$", first.line)
  is_krakenu_fmt <- grepl("^.?%\treads\ttaxReads\tkmers", first.line)
  is_kaiju_fmt <- grepl("^  *%\t  *reads", first.line)
  
  ntabs <- lengths(regmatches(first.line, gregexpr("\t", first.line)))
  
  nrows <- ifelse(isTRUE(check_file), 5, -1)
  if (!is_krakenu_fmt && is_kaiju_fmt) {
    cont <- readLines(myfile)
    cont <- cont[!grepl("^-", cont)]
    cont <- sub(".*\t  *","", cont)
    cont <- sub("; ?$","", cont)
    report <- utils::read.delim(textConnection(cont), stringsAsFactors = FALSE)
    colnames(report) <- c("taxonReads", "taxLineage")
    report$cladeReads <- report$taxonReads
    
    report$taxLineage <- gsub("^","-_",report$taxLineage)
    report$taxLineage <- gsub("; ","|-_",report$taxLineage)
    report$taxLineage <- gsub("-_Viruses", "d_Viruses", report$taxLineage, fixed=T)
    report$taxLineage <- gsub("-_cellular organisms|-_Bacteria", "-_cellular organisms|d_Bacteria", report$taxLineage, fixed=T)
    report$taxLineage <- gsub("-_cellular organisms|-_Eukaryota", "-_cellular organisms|d_Eukaryota", report$taxLineage, fixed=T)
    report$taxLineage <- gsub("-_cellular organisms|-_Archaea", "-_cellular organisms|d_Archaea", report$taxLineage, fixed=T)
    report$taxLineage[1:(length(report$taxLineage)-1)] <- paste0("-_root|", report$taxLineage[1:(length(report$taxLineage)-1)])
    
    report$taxLineage[report$taxLineage=="-_unclassified"] <- "u_unclassified"
    
    new_counts <- integer(length = 0)
    for (j in seq_len(nrow(report))) {
      count <- report$cladeReads[j]
      tl <- report$taxLineage[j]
      tl2 <- sub("\\|[^|]*$","", tl)
      while (tl2 != tl) {
        if (tl2 %in% names(new_counts)) {
          new_counts[tl2] <- new_counts[tl2] + count
        } else {
          new_counts[tl2] <- count
        }
        tl <- tl2
        tl2 <- sub("\\|[^|]*$","", tl)
      }
    }
    report <- rbind(report,
                    data.frame(taxonReads=0,taxLineage=names(new_counts),cladeReads=as.integer(new_counts)))
    tl_order <- order(report$taxLineage)
    tl_order <- c(tl_order[length(tl_order)],tl_order[-length(tl_order)])
    report <- report[tl_order, c("taxLineage", "taxonReads", "cladeReads")]
  } else if (is_metaphlan3_fmt) {
    report <- tryCatch({
      utils::read.table(myfile,sep="\t",header = F,
                        quote = "",stringsAsFactors=FALSE,
                        comment.char = "#", nrows = nrows,
                        col.names = c("taxLineage", "taxID", "cladeReads", "additional_species"),
                        check.names=FALSE)
    }, error = function(x) NULL, warning = function(x) NULL)
    if (is.null(report)) { return(NULL); }
    
    report$taxID <- sub(".*\\|", "", report$taxID)
    
  } else if (has_header) {
    report <- tryCatch({
      utils::read.table(myfile,sep="\t",header = T,
                        quote = "",stringsAsFactors=FALSE,
                        comment.char = ifelse(is_metaphlan2_fmt, "", "#"), nrows = nrows, 
                        check.names=FALSE)
    }, error = function(x) NULL, warning = function(x) NULL)
    if (is.null(report)) { return(NULL); }
    
    # harmonize column names. TODO: Harmonize them in the scripts!
    colnames(report)[colnames(report) %in% c("#%","%","clade_perc","perc","percReadsClade")] <- "percentage"
    colnames(report)[colnames(report) %in% c("reads","numReadsClade","n_reads_clade","n.clade","n-clade")] <- "cladeReads"
    colnames(report)[colnames(report) %in% c("taxReads","numReadsTaxon","n_reads_taxo","n.stay","n-stay")] <- "taxonReads"
    colnames(report)[colnames(report) %in% c("rank","tax_taxRank","level")] <- "taxRank"
    colnames(report)[colnames(report) %in% c("tax","taxonid")] <- "taxID"
    colnames(report)[colnames(report) %in% c("indentedName","taxName")] <- "name"
    colnames(report)[colnames(report) %in% c("dup")] <- "kmerDuplicity"
    colnames(report)[colnames(report) %in% c("cov")] <- "kmerCoverage"
  } else {
    report <- NULL
    if (ntabs == 5) {
      col_names <-  c("percentage","cladeReads","taxonReads","taxRank","taxID","name")
    } else if (ntabs == 7) {
      col_names <- c("percentage","cladeReads","taxonReads", "n_unique_kmers","n_kmers", "taxRank","taxID","name")
    } 
    report <- tryCatch({
      utils::read.table(myfile,sep="\t",header = F,
                        col.names = col_names,
                        quote = "",stringsAsFactors=FALSE,
                        nrows = nrows)
    }, error=function(x) NULL, warning=function(x) NULL)
    
    if (is.null(report)) {
      dmessage(paste("Warning: File",myfile,"does not have the required format"))
      return(NULL); 
    }
    
  }
  
  if (ncol(report) < 2) {
    dmessage(paste("Warning: File",myfile,"does not have the required format"))
    return(NULL) 
  }
  if (is_metaphlan2_fmt || is_metaphlan3_fmt) {
    # Metaphlan report
    colnames(report)[1] <- "taxLineage"
    colnames(report)[colnames(report) == "Metaphlan2_Analysis"] <- "cladeReads"
    report <- report[order(report$taxLineage), ]
    report$taxLineage <- gsub("_"," ",report$taxLineage)
    report$taxLineage <- gsub("  ","_",report$taxLineage)
    report$taxLineage <- paste0("-_root|", report$taxLineage)
    root_lines <- data.frame(taxLineage=c("u_unclassified","-_root"),"cladeReads"=c(0,100), stringsAsFactors = F)
    
    if ("taxID" %in% colnames(report)) {
      root_lines <- cbind(root_lines, "taxID" = c(0, 1))
      report <- report[, c("taxLineage", "cladeReads", "taxID")]
    }
    
    report <- rbind(root_lines, report)
  }
  
  if (all(c("name","taxRank") %in% colnames(report)) && !"taxLineage" %in% colnames(report)) {
    # Kraken report
    report$depth <- nchar(gsub("\\S.*","",report$name))/2
    if (!all(report$depth == floor(report$depth))) {
      warning("Depth doesn't work out!")
      return(NULL)
    }
    report$name <- gsub("^ *","",report$name)
    
    # 'fix' taxRank
    table(report$taxRank)
    allowed_taxRanks <- c("U", "S", "G", "F", "C", "D", "O", "K", "P")
    report$taxRank[report$taxRank=="class"] <- "C"
    report$taxRank[report$taxRank=="family"] <- "F"
    report$taxRank[report$taxRank=="genus"] <- "G"
    report$taxRank[report$taxRank=="superkingdom"] <- "D"
    report$taxRank[report$taxRank=="kingdom"] <- "K"
    report$taxRank[report$taxRank=="order"] <- "O"
    report$taxRank[report$taxRank=="phylum"] <- "P"
    report$taxRank[report$taxRank=="species"] <- "S"
    report$taxRank[report$name=="unclassified"] <- "U"
    report$taxRank[!report$taxRank %in% allowed_taxRanks] <- "-"
    
    report$name <- paste(tolower(report$taxRank),report$name,sep="_")
    
    rownames(report) <- NULL
    
    # make taxLineage path
    report$taxLineage <- report$name
    n <- nrow(report)
    depths <- report$depth
    taxLineages <- report$name
    taxLineages_p <- as.list(seq_along(report$name))
    
    depth_row_tmp <- c(1:25)
    
    for (current_row in seq(from=1, to=nrow(report))) {
      dcr <- depths[current_row]
      depth_row_tmp[dcr+1] <- current_row
      if (dcr >= 1) {
        prev_pos <- depth_row_tmp[[dcr]]
        taxLineages_p[[current_row]] <- c(taxLineages_p[[prev_pos]], current_row)
      }
    }
    report$taxLineage <- sapply(taxLineages_p, function(x) paste0(taxLineages[x], collapse="|"))
    
  } else if ("taxLineage" %in% colnames(report)) {
    taxLineages <- strsplit(report$taxLineage, "|", fixed=TRUE)
    
    if (!"name" %in% colnames(report))
      report$name <- sapply(taxLineages, function(x) x[length(x)])
    
    if (!"depth" %in% colnames(report)) {
      report$depth <- sapply(taxLineages, length) - 1
    }
    if (!"taxRank" %in% colnames(report))
      report$taxRank <- toupper(substr(report$name, 0, 1))
  }
  
  if (!all(c("name","taxRank") %in% colnames(report)) ||
      nrow(report) < 2 ||
      !((report[1,"name"] == "u_unclassified" && report[2,"name"] == "-_root") || report[1,"name"] == "-_root")) {
    dmessage(paste("Warning: File",myfile,"does not have the required format"))
    str(report)
    return(NULL)
  }
  
  if (!"taxonReads" %in% colnames(report)) {
    parent <- sub("^\\(.*\\)\\|.*$", "\\1", report$taxLineage)
    taxLineages <- strsplit(report$taxLineage, "|", fixed=TRUE)
    # fix taxonReads
    report$taxonReads <- report$cladeReads - sapply(report$name, function(x) sum(report$cladeReads[parent == x]))
    report$taxonReads[report$taxonReads <= 0.00001] <- 0  # fix for rounding in percentages by MetaPhlAn
  }
  
  report$percentage <- signif(report$cladeReads/sum(report$taxonReads),6) * 100
  if ('n_unique_kmers'  %in% colnames(report))
    report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers,na.rm=T),6) * 100
  
  if ("taxID" %in% colnames(report)) {
    std_colnames <- c("percentage","cladeReads","taxonReads","taxRank", "taxID","name")
  } else {
    std_colnames <- c("percentage","cladeReads","taxonReads","taxRank","name")
  }
  stopifnot(all(std_colnames %in% colnames(report)))
  report[, c(std_colnames, setdiff(colnames(report), std_colnames))]
}


read_reports <- function(report_files, report_names = NULL) {
  if (length(report_files) == 0) {
    return()
  }
  
  if (length(report_files) == 1 && isTRUE(file.info(report_files)$isdir)) {
    report_files <- read_sample_data(report_files, ext = NULL)
  }
  
  if (is.data.frame(report_files) && all(c("ReportFilePath", "Name") %in% colnames(report_files))) {
    report_names <- report_files$Name
    report_files <- report_files$ReportFilePath
  }
  
  if (is.null(report_names)) {
    report_names = basename(report_files)
  }
  
  if (any(duplicated(report_names))) {
    report_names = report_files
  }
  
  dmessage("Reading ", length(report_files), " reports ...")
  n_reports <- length(report_files)
  
  my_reports <- lapply(seq_along(report_files), function(i) {
    report <- tryCatch({
      read_report(report_files[i])
    }, error = function(e) {
      stop(paste("Error reading file", report_files[i], ": ", e$message))
    })
    
    if (is.null(report) || nrow(report) == 0) {
      stop(paste("Error reading file", report_files[i]))
    }
    report
  })
  
  names(my_reports) <- report_names
  
  my_reports[sapply(my_reports, length) > 0]
}


summarize_report <- function(my_report) {
  my_report <- my_report[!duplicated(my_report$name),]
  row.names(my_report) <- my_report[["name"]]
  unidentified_reads <- my_report["u_unclassified","cladeReads"]
  identified_reads <- my_report["-_root","cladeReads"]
  artificial_reads <- zero_if_na1(my_report["s_synthetic construct","cladeReads"])
  human_reads <- zero_if_na1(my_report["s_Homo sapiens","cladeReads"])
  chordate_reads <- zero_if_na1(my_report["p_Chordata","cladeReads"])
  root_reads <- zero_if_na1(my_report["-_root","taxonReads"])
  
  data.frame(
    number_of_raw_reads=unidentified_reads+identified_reads,
    classified_reads=identified_reads,
    chordate_reads=chordate_reads,
    artificial_reads=artificial_reads,
    unclassified_reads=unidentified_reads,
    microbial_reads=identified_reads-chordate_reads-artificial_reads-root_reads,
    bacterial_reads=zero_if_na1(my_report["d_Bacteria","cladeReads"]) + zero_if_na1(my_report["k_Bacteria","cladeReads"]), # MetaPhLan reports bacteria as kingdom; Kraken as domain. Sum them
    viral_reads=zero_if_na1(my_report["d_Viruses","cladeReads"]) + zero_if_na1(my_report["k_Viruses","cladeReads"]), # same as for Bacteria
    fungal_reads=zero_if_na1(my_report["k_Fungi","cladeReads"]),
    protozoan_reads=sum(zero_if_na1(my_report[names(protist_taxids),"cladeReads"]))
  )
}


zero_if_na1 <- function(x) {
  x[is.na(x)] <- 0
  x
}


summarize_reports <- function(reports) {
  do.call(rbind, lapply(reports, summarize_report))
}


beautify_string <- function(x) {
  x <- gsub("[\\._]"," ",x)
  x <- sub("^([[:alpha:]])", "\\U\\1", x, perl=TRUE)
  x
}


filter_taxon <- function(report, filter_taxon, rm_clade = TRUE, do_message=FALSE) {
  taxon_depth <- NULL
  taxonReads <- 0
  
  pos.taxons <- which(sub("._","",report$name) %in% filter_taxon)
  if (length(pos.taxons) == 0) {
    return(report)
  }
  
  row_seq <- seq_len(nrow(report))
  rows_to_delete <- rep(FALSE,nrow(report))
  
  taxon_depths <- report[pos.taxons,"depth"]
  if (isTRUE(rm_clade)) {
    taxonReads <- report[pos.taxons,"cladeReads"]
  } else {
    taxonReads <- report[pos.taxons,"taxonReads"]
    report[pos.taxons,"taxonReads"] <- 0
  }
  
  
  for (i in seq_along(pos.taxons)) {
    pos.taxon <- pos.taxons[i]
    if (pos.taxon == 1) {
      rows_to_delete[1] <- TRUE
      next
    }
    taxon_depth <- taxon_depths[i]
    taxonReads <- taxonReads[i]
    
    if (rm_clade) {
      tosum_below <-  row_seq >= pos.taxon & report$depth <= taxon_depth
      taxons_below <- cumsum(tosum_below) == 1
      rows_to_delete[taxons_below] <- TRUE
    }
    rows_to_update <- c(pos.taxon)
    
    taxons_above <- seq_len(nrow(report)) < pos.taxon & report$depth == taxon_depth
    
    any_stays <- FALSE
    prev_taxon_depth <- taxon_depth
    taxons_above <- c()
    for (i in seq(from=(pos.taxon-1),to=1)) {
      curr_taxon_depth <- report[i,"depth"]
      if (curr_taxon_depth < prev_taxon_depth) {
        if (!any_stays) {
          if (report[i,"cladeReads"] == taxonReads) {
            rows_to_delete[i] <- TRUE
            if (do_message)
              dmessage("Deleting ",report[i,"name"])
          } else {
            any_stays <- TRUE
          }
        }
        if (!rows_to_delete[i]) {
          rows_to_update <- c(rows_to_update, i)
          if (do_message)
            dmessage("Updating ",report[i,"name"])
        }
        prev_taxon_depth <- curr_taxon_depth
      } else {
        any_stays <- TRUE
      }
    }
    report[rows_to_update, "cladeReads"] <- report[rows_to_update, "cladeReads"] - taxonReads
  }
  
  if (rm_clade)
    report[!rows_to_delete,]
  else
    report
}


build_sankey_network <- function(my_report, taxRanks =  c("D","K","P","C","O","F","G","S"), maxn=10,
                         zoom = F, title = NULL,
                         ...) {
  stopifnot("taxRank" %in% colnames(my_report))
  if (!any(taxRanks %in% my_report$taxRank)) {
    warning("report does not contain any of the taxRanks - skipping it")
    return()
  }
  my_report <- subset(my_report, taxRank %in% taxRanks)
  my_report <- plyr::ddply(my_report, "taxRank", function(x) x[utils::tail(order(x$cladeReads,-x$depth), n=maxn), , drop = FALSE])
  
  my_report <- my_report[, c("name","taxLineage","taxonReads", "cladeReads","depth", "taxRank")]
  
  my_report <- my_report[!my_report$name %in% c('-_root'), ]
  
  splits <- strsplit(my_report$taxLineage, "\\|")
  
  root_nodes <- sapply(splits[sapply(splits, length) ==2], function(x) x[2])
  
  sel <- sapply(splits, length) >= 3
  splits <- splits[sel]
  
  links <- data.frame(do.call(rbind,
                              lapply(splits, function(x) utils::tail(x[x %in% my_report$name], n=2))), stringsAsFactors = FALSE)
  colnames(links) <- c("source","target")
  links$value <- my_report[sel,"cladeReads"]
  
  my_taxRanks <- taxRanks[taxRanks %in% my_report$taxRank]
  taxRank_to_depth <- stats::setNames(seq_along(my_taxRanks)-1, my_taxRanks)
  
  
  nodes <- data.frame(name=my_report$name,
                      depth=taxRank_to_depth[my_report$taxRank],
                      value=my_report$cladeReads,
                      stringsAsFactors=FALSE)
  
  for (node_name in root_nodes) {
    diff_sum_vs_all <- my_report[my_report$name == node_name, "cladeReads"] - sum(links$value[links$source == node_name])
    if (diff_sum_vs_all > 0) {
      nname <- paste("other", sub("^._","",node_name))
    }
  }
  
  names_id = stats::setNames(seq_len(nrow(nodes)) - 1, nodes[,1])
  links$source <- names_id[links$source]
  links$target <- names_id[links$target]
  links <- links[links$source != links$target, ]
  
  nodes$name <- sub("^._","", nodes$name)
  links$source_name <- nodes$name[links$source + 1]
  
  if (!is.null(links))
    sankeyD3::sankeyNetwork(
      Links = links,
      Nodes = nodes,
      doubleclickTogglesChildren = TRUE,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      NodeGroup = "name",
      NodePosX = "depth",
      NodeValue = "value",
      dragY = TRUE,
      xAxisDomain = my_taxRanks,
      numberFormat = "pavian",
      title = title,
      nodeWidth = 15,
      linkGradient = TRUE,
      nodeShadow = TRUE,
      nodeCornerRadius = 5,
      units = "cladeReads",
      fontSize = 12,
      iterations = maxn * 100,
      align = "none",
      highlightChildLinks = TRUE,
      orderByPath = TRUE,
      scaleNodeBreadthsByString = TRUE,
      zoom = zoom,
    )
}


save_sankey_plot <- function(sankey_plot, 
                             format = "html", 
                             filename = "sankey_plot", 
                             width = 1920, 
                             height = 1080) {
  
  # 确保sankey_plot和format有有效的值
  stopifnot(!missing(sankey_plot))
  stopifnot(format %in% c("html", "png", "jpg", "pdf"))
  
  # 组合文件名和格式
  filename <- paste0(filename, ".sankey", ".", format)
  
  if (format == "html") {
    htmlwidgets::saveWidget(sankey_plot, filename, selfcontained = TRUE)
  } else {
    temp_file <- tempfile(fileext = ".html")
    htmlwidgets::saveWidget(sankey_plot, temp_file, selfcontained = TRUE)
    
    if (format %in% c("png", "jpg")) {
      webshot::webshot(url = temp_file, 
                       file = filename, 
                       vwidth = width, 
                       vheight = height, 
                       zoom = 4)
    } else if (format == "pdf") {
      webshot::webshot(url = temp_file, 
                       file = filename, 
                       vwidth = width, 
                       vheight = height)
    }
  }
  
  dmessage(paste0("Sankey plot (based on filtered Kraken report) saved as: ", filename))
}


merge_reports <- function(my_reports, col_names = NULL, fix_taxnames = TRUE, update_progress = FALSE,
                           id_cols = c("name", "taxRank", "taxID", "taxLineage"),
                           numeric_cols = c("cladeReads","taxonReads")) {
  common_colnames <- Reduce(intersect, lapply(my_reports, colnames))
  
  if (is.null(my_reports) || length(my_reports) == 0)
    return(NULL)
  
  if (!"taxID" %in% common_colnames)
    id_cols <- c("name", "taxRank", "taxLineage")
  
  if (!all(c(id_cols,numeric_cols) %in% common_colnames)) {
    stop("Not all required columns are present Required: ",
         paste0(sort(c(id_cols,numeric_cols)), collapse=", "),". Present: ",
         paste0(sort(common_colnames), collapse=", "))
  }
  
  if (fix_taxnames && length(my_reports) > 1) {
    if (!"taxID" %in% common_colnames) {
      dmessage("Can't fix taxnames without taxID!")
    } else {
      c_id_to_name1 <- lapply(my_reports, function(r) {
        r[,c("taxID","name")]
      }) %>% do.call(rbind, .)
      rownames(c_id_to_name1) <- NULL
      c_id_to_name <- unique(c_id_to_name1[order(c_id_to_name1$taxID),])
      rownames(c_id_to_name) <- NULL
      sel_diff_taxid <- c_id_to_name$name %in% c_id_to_name$name[duplicated(c_id_to_name$name)]
      sel_diff_name <- c_id_to_name$taxID %in% c_id_to_name$taxID[duplicated(c_id_to_name$taxID)]
      if (any(sel_diff_name)) {
        dmessage("The following taxons have the same taxIDs but differing names:")
        print(c_id_to_name[sel_diff_name,])
      }
    }
  }
  
  my_reports <- lapply(seq_along(my_reports), function(i) {
    mm <- my_reports[[i]][,c(id_cols, numeric_cols)]
    colnames(mm)[colnames(mm) %in% numeric_cols] <- sprintf("%s.%s", colnames(mm)[colnames(mm) %in% numeric_cols], i)
    mm
  })
  
  if (length(my_reports) > 1) {
    if (isTRUE(update_progress)) {
      dmessage("Merging Kraken reports progress...")
      pb <- txtProgressBar(min = 0, max = length(my_reports) - 1, style = 3)
    }
    
    merged_reports <- Reduce(
      function(merged_rep, rep_index) {
        if (isTRUE(update_progress)) {
          setTxtProgressBar(pb, rep_index - 1)
          flush.console()
        }
        dplyr::full_join(merged_rep, my_reports[[rep_index]], by = id_cols)
      }, 
      seq(from=2, to=length(my_reports)),
      init=my_reports[[1]]
    )
    
    if (isTRUE(update_progress)) {
      close(pb)
    }
  } else {
    merged_reports <- my_reports[[1]]
  }
  
  tax_data <- merged_reports[, id_cols, drop = FALSE]
  tax_data[, 1] <- sub("^[a-z-]_", "", tax_data[, 1])
  
  idx_cladeReads <- seq(from = length(id_cols) + 1, to = ncol(merged_reports), by = 2)
  idx_taxonReads <- seq(from = length(id_cols) + 2, to = ncol(merged_reports), by = 2)
  
  cladeReads <- as.matrix(merged_reports[, idx_cladeReads, drop = FALSE])
  taxonReads <- as.matrix(merged_reports[, idx_taxonReads, drop = FALSE])
  
  if (!is.null(col_names)) {
    stopifnot(length(col_names) == ncol(cladeReads))
    colnames(cladeReads) <- col_names
    colnames(taxonReads) <- col_names
  }
  
  list(tax_data = tax_data, cladeReads = cladeReads, taxonReads = taxonReads)
}

################################################################################

# 记录开始时间
start_time <- Sys.time()

# 默认过滤物种列表
default_filt <- c("Chordata", "other sequences")

# 创建简写和完整名称的映射
abbreviation_map <- list(
  a = "artificial sequences",
  c = "Chordata",
  ca = "Cutibacterium acnes",
  e = "Enterobacteria phage phiX174 sensu lato",
  en = "Enterobacteriales",
  ec = "Escherichia coli",
  h = "Homo sapiens",
  m = "Mus musculus",
  r = "Ralstonia pickettii",
  s = "Saccharomyces cerevisiae",
  u = "unclassified",
  o = "other sequences"
)

# 将 default_filt 转换为一个由逗号分隔的字符串
default_filt_str <- paste(default_filt, collapse = ",")

# 创建一个字符串，该字符串包含缩写和全称的映射
abbreviation_help <- paste(
  unlist(lapply(names(abbreviation_map), function(abbr) {
    sprintf("%s: %s", abbr, abbreviation_map[[abbr]])
  })),
  collapse = ", "
)

# 创建命令行参数解析器
parser <- ArgumentParser(description='This script processes kraken2 pathogen analysis reports. Use --info to view more detailed information (strongly recommended).')

# 添加命令行选项
parser$add_argument('--info', action='store_true', help='Print full description information and exit.')
parser$add_argument("-i", "--input", type="character", required=TRUE,
                    dest="input",
                    help="Path to the directory containing kraken2 pathogen analysis reports or a specific .k2report file. This script is designed to process .k2report files.",
                    metavar="PATH")
parser$add_argument("-o", "--output", type="character", default=getwd(),
                    dest="output",
                    help="Path to the output directory. If not provided, defaults to the current working directory.",
                    metavar="PATH")
parser$add_argument("-f", "--filter", type="character", default=default_filt_str,
                    dest="filter",
                    help=paste("Comma-separated list of abbreviations to filter (e.g., 'c,h,m'). Default is '", default_filt_str, "'. Valid abbreviations: ", abbreviation_help),
                    metavar="ABBREVS")
parser$add_argument("--sankey-format", type="character", default="html",
                    dest="sankey_format",
                    help="Sankey plot format: html, pdf, jpg, or png. If not provided, default is html.",
                    metavar="FORMAT")
parser$add_argument("--limit-rows", action="store_true", default=FALSE, 
                    dest="limit_rows",
                    help="Limit the number of rows in the biological classification data output. If selected data exceeds 100 rows, only the top 100 rows related to biological classifications are shown. Default is FALSE.")
parser$add_argument('-v', '--version', action='version', version='%(prog)s v2.0')

# 先获取命令行参数
cmd_args <- commandArgs(trailingOnly = TRUE)

# 检查是否有 --info 参数
if ("--info" %in% cmd_args) {
  print_full_description()
  q(status = 0)  # 正常退出
}

# 如果没有提供任何参数，显示帮助信息并退出
if (length(cmd_args) == 0) {
  parser$print_help()
  q(status = 1)  # 退出状态 1 表示错误
}

# 解析命令行参数
args <- parser$parse_args()

dmessage("Analysis started ...")

# 检查输入路径是目录还是文件
if (file.info(args$input)$isdir) {
  # 输入路径是目录
  report_dir <- args$input
  all_files <- list.files(path = report_dir, full.names = TRUE)
  report_files <- all_files[grepl("\\.k2report$", all_files)]
  dmessage(length(report_files), " Kraken report file(s) provided in directory: ", paste(basename(report_files), collapse = ", "))
} else if (grepl("\\.k2report$", args$input)) {
  # 输入路径是.k2report文件
  report_files <- args$input
  dmessage("Single Kraken report file provided: ", basename(args$input))
} else {
  stop("Invalid input. Please provide a directory or a .k2report file.")
}

# 检查是否找到.k2report的报告文件
if (length(report_files) == 0) {
  stop("No report files with the extension .k2report were found.")
}

# 检测用户是否提供了过滤条件并且这些条件都是已知的
if (!is.null(args$filter) && args$filter != default_filt_str) {
  # 用户提供了自定义过滤条件
  dmessage("Custom filter provided: ", args$filter)
  # 分割用户提供的过滤条件缩写
  abbrevs <- unlist(strsplit(args$filter, ","))
  # 将缩写转换为完整形式，如果不是已知缩写，保持原样
  filt <- sapply(abbrevs, function(abbr) {
    if (abbr %in% names(abbreviation_map)) {
      return(abbreviation_map[[abbr]])
    } else {
      return(abbr)  # 保留未知缩写
    }
  })
} else {
  # 用户没有提供过滤条件，使用默认值
  dmessage("Using default filter: ", default_filt_str)
  # 默认值已经是缩写形式，无需转换
  filt <- unlist(strsplit(default_filt_str, ","))
}

# 检查是否所有的过滤条件都是已知的
unknown_filt <- filt[!filt %in% abbreviation_map]
if(length(unknown_filt) > 0) {
  # 创建并打印一个格式化的缩写映射表
  formatted_map <- paste(sapply(names(abbreviation_map), function(k) {
    paste(k, "->", abbreviation_map[[k]])
  }), collapse = "\n")
  dmessage("The available abbreviation mappings are:\n", formatted_map)
  stop("Unknown filter abbreviations: ", paste(unknown_filt, collapse = ", "))
}

# 打印过滤条件
dmessage("Filter conditions applied: ", paste(filt, collapse = ", "))

# 验证 --sankey-format 选项的值
valid_formats <- c("html", "pdf", "jpg", "png")
if (!args$sankey_format %in% valid_formats) {
  stop(sprintf("Invalid sankey format: %s. Valid formats are: %s.",
               args$sankey_format, paste(valid_formats, collapse=", ")))
}

# 将 args$sankey-format 的值赋给 sankey_format
sankey_format <- args$sankey_format

# 如果提供了 -o/--output 选项
if (!is.null(args$output) && args$output != getwd()) {
  # 将相对路径转换为绝对路径
  absolute_output_path <- normalizePath(args$output, mustWork = FALSE)
  
  # 检查路径是否存在
  if (file.exists(absolute_output_path)) {
    # 检查它是否是一个目录
    if (!isTRUE(file.info(absolute_output_path)$isdir)) {
      stop(sprintf("The specified output path is not a directory: %s", absolute_output_path))
    } else {
      dmessage(sprintf("Using existing output directory: %s", absolute_output_path))
    }
  } else {
    # 如果路径不存在，则尝试创建它
    dir_create_error <- try(dir.create(absolute_output_path, recursive = TRUE), silent = TRUE)
    
    # 如果创建路径失败，停止脚本并显示错误消息
    if (inherits(dir_create_error, "try-error")) {
      stop(sprintf("Failed to create output directory: %s", absolute_output_path))
    } else {
      dmessage(sprintf("Created output directory: %s", absolute_output_path))
    }
  }
  # 设置 output_dir 变量的值
  output_dir <- absolute_output_path
} else {
  # 使用当前工作目录作为输出路径
  output_dir <- getwd()
  dmessage(sprintf("Using current working directory as output path: %s", output_dir))
}

# 样本分类分类数据输出展示限制
limit_rows <- args$limit_rows

# 打印limit_rows的值
if (limit_rows) {
  dmessage("Row limiting is enabled. Only the top 100 rows will be shown.")
} else {
  dmessage("Row limiting is not enabled. All rows will be shown.")
}

################################################################################
# 1. Sample set summary 

# 从文件路径中提取文件名（不带扩展名）作为报告名称
report_names <- sapply(report_files, function(file) {
  file_base <- basename(file)  # 获取文件的基本名称
  file_name <- file_path_sans_ext(file_base)  # 移除扩展名
  return(file_name)
})

# 读取报告
reports <- read_reports(report_files = report_files, report_names = report_names)

protist_taxids <- c("-_Diplomonadida"=5738,
                    "-_Amoebozoa"=554915,
                    "-_Alveolata"=33630)

# 调用函数summarize_reports来汇总多个报告的数据
samples_summary <- summarize_reports(reports)

# 将汇总数据的行名（通常是自动生成的）添加为一个新列“Name”
samples_summary$Name <- rownames(samples_summary)

# 定义需要特别处理的列名
extra_cols <- c("Name")

# 重新组织数据框列的顺序，确保“Name”列在前，其余列按原顺序排列
samples_summary <- samples_summary[,c(extra_cols, setdiff(colnames(samples_summary),extra_cols))]

# 使用beautify_string函数美化数据框的列名，使其更易读
colnames(samples_summary) <- beautify_string(colnames(samples_summary))

# 创建一个新的数据框用于存储百分比形式的汇总数据
samples_summary_percent <- samples_summary

# 修改samples_summary_percent的列名，从第三列开始添加" (%)"
colnames(samples_summary_percent)[3:ncol(samples_summary_percent)] <- 
  paste0(colnames(samples_summary_percent)[3:ncol(samples_summary_percent)], " (%)")

# 计算每个样本的百分比数据，这里是将原始数据转换为百分比形式
# 'sweep' 函数用于对数据应用一些函数操作，这里是从每行的数值中减去第二列的数值，然后每个数值除以第二列的对应数值
# 'signif' 函数用于四舍五入处理计算结果，保留4位有效数字
# 结果乘以100，转换为百分比
samples_summary_percent[, 3:ncol(samples_summary)] <- 
  100 * signif(
    sweep(
      samples_summary[, 3:ncol(samples_summary)], 
      1, 
      samples_summary[, 2], 
      `/`
    ), 
    4
  )

# 列表的名称将成为Excel工作表的名称
sheets <- list("Classification summary" = samples_summary_percent,
               "Raw read numbers" = samples_summary)

# 指定要保存的文件名
summary_file <- file.path(output_dir, "sample_summary.xlsx")

# 检查文件是否存在
if (file.exists(summary_file)) {
  warning(sprintf("File '%s' already exists. Overwriting...", summary_file))
}

# 使用write_xlsx函数将数据写入文件
# 使用try函数来捕获错误
result <- try(write_xlsx(sheets, path = summary_file), silent = TRUE)

# 检查是否有错误，如果有，停止脚本执行
if (inherits(result, "try-error")) {
  stop("Failed to write data to file: ", summary_file)
}

# 如果没有错误，打印成功消息
dmessage(paste0("Kraken report summary successfully written to: ", summary_file))

################################################################################
# 2. Classification results

# 2.1 桑基图
all_names <- sub("^._","",sort(unique(unlist(sapply(reports,function(x) x$name[x$taxRank != "-"])))))
colourScale <- sankeyD3::JS(sprintf("d3.scaleOrdinal().range(d3.schemeCategory20b).domain([%s])", 
                                    paste0('"',c(all_names,"other"),'"',collapse=",")))

# 使用purrr包的walk函数，它是lapply的一个变体，但不返回任何值。
walk(names(reports), function(n) {
    my_report <- reports[[n]]
    
    for (f in filt)
      my_report <- filter_taxon(my_report, f)
    
    # 保存 my_report 为 tsv 文件
    report_file_path = file.path(output_dir, paste0(n, ".report", ".tsv"))
    write_tsv(my_report, report_file_path)
    dmessage(paste0("Kraken report source file saved to: ", report_file_path))
    
    sankey_plot <- build_sankey_network(my_report, nodePadding=13, xScalingFactor=.9, nodeStrokeWidth=0, zoom = FALSE, colourScale = colourScale, LinkGroup = "source_name")
    
    save_sankey_plot(sankey_plot, format=sankey_format, filename=file.path(output_dir,n))
})

# 2.2 样本分类
merged_reports <- merge_reports(reports, col_names = samples_summary_percent$Name, update_progress = TRUE)
taxonReads <- merged_reports$taxonReads
cladeReads <- merged_reports$cladeReads
tax_data <- merged_reports[["tax_data"]]

# 创建逻辑向量，向量定义了哪些行（即数据集中的条目）属于特定的生物类别
# 通过检查“taxLineage”列中的分类信息以及“taxRank”列中的分类等级来实现的

# 每个条件的解释：
# sel_bacteria:grepl("d_Bacteria", tax_data[,"taxLineage"]) 检查分类谱系中是否存在“细菌”域（domain）。
# tax_data[,"taxRank"] == 'S' 确保只选择种（species）级别的条目。结果是一个逻辑向量，指示哪些行属于细菌种。
# sel_viruses:类似地，这个条件寻找谱系中的“病毒”域，并且是种级别的条目。
# sel_fungi:这里，代码寻找属于“真菌”界（kingdom）的条目，这是一个比域更具体的分类级别。
# sel_euk:这个条件选择所有属于“真核生物”域的条目，不论它们属于哪个更低的分类级别。
# sel_protists:这个条件稍微复杂一些，因为它首先选择所有“真核生物”域的条目，但排除了真菌（由 !sel_fungi 实现）和脊索动物（通过 !grepl("p_Chordata", tax_data[,"taxLineage"]) 排除）。

sel_bacteria = grepl("d_Bacteria",tax_data[,"taxLineage"]) & tax_data[,"taxRank"] == 'S'
sel_viruses = grepl("d_Viruses",tax_data[,"taxLineage"])  & tax_data[,"taxRank"] == 'S'
sel_fungi = grepl("k_Fungi",tax_data[,"taxLineage"])  & tax_data[,"taxRank"] == 'S'
sel_euk = grepl("d_Eukaryota",tax_data[,"taxLineage"]) & tax_data[,"taxRank"] == 'S'
sel_protists = grepl("d_Eukaryota",tax_data[,"taxLineage"]) & !sel_fungi &!grepl("p_Chordata",tax_data[,"taxLineage"]) & tax_data[,"taxRank"] == 'S'

my_df <- data.frame(Name=tax_data$name,Max=apply(cladeReads,1,max,na.rm=TRUE),cladeReads,Lineage=pavian:::beautify_taxLineage(tax_data$taxLineage), stringsAsFactors = FALSE)

process_data <- function(df, sel, category, limit_rows = TRUE) {
  # 筛选数据
  selected_data <- df[sel, ]
  
  # 根据limit_rows参数决定是否限制行数
  if (limit_rows) {
    # 如果选择的数据超过100行，我们只取前100行。
    if (nrow(selected_data) > 100) {
      selected_data <- selected_data[order(selected_data$Max, decreasing = TRUE), ]
      selected_data <- head(selected_data, 100)
      dmessage(sprintf("Showing top 100 of %s species for category: %s.", sum(sel), category))
    } else {
      dmessage(sprintf("Showing all %s species for category: %s.", nrow(selected_data), category))
    }
  } else {
    # 如果不限制行数，显示所有选定的数据
    dmessage(sprintf("Showing all %s species for category: %s.", sum(sel), category))
  }
  
  # 返回筛选后的数据
  return(selected_data)
}

# 对不同的类别调用函数，并将结果存储在不同的变量中
data_bacteria <- process_data(my_df, sel_bacteria, "Bacteria", limit_rows)
data_viruses <- process_data(my_df, sel_viruses, "Viruses", limit_rows)
data_eukaryotes <- process_data(my_df, sel_euk, "Eukaryotes", limit_rows)
data_fungi <- process_data(my_df, sel_fungi, "Eukaryotes/Fungi", limit_rows)
data_protists <- process_data(my_df, sel_protists, "Eukaryotes/Protists", limit_rows)

# 创建一个命名列表来保存你的数据表
data_list <- list(
  Bacteria = data_bacteria,
  Viruses = data_viruses,
  Eukaryotes = data_eukaryotes,
  `Eukaryotes Fungi` = data_fungi,
  `Eukaryotes Protists` = data_protists
)

# 指定要保存的文件路径和文件名
classification_file <- file.path(output_dir, "classification_results.xlsx")

# 使用try函数来捕获错误
result <- try(write_xlsx(data_list, path = classification_file), silent = TRUE)

# 检查是否有错误，如果有，停止脚本执行
if (inherits(result, "try-error")) {
  stop("Failed to write data to file: ", classification_file)
}

# 如果没有错误，打印成功消息
dmessage(paste0("Classification results successfully written to: ", classification_file))

# 记录结束时间
end_time <- Sys.time()

# 结束分析
duration <- end_time - start_time
dmessage(paste0("Analysis completed successfully in ", round(duration, 2), " seconds."))

