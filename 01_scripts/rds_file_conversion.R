# =========================
# Copy 02_results -> 02_resultsv2
# Then find all .rds files inside 02_resultsv2
# and convert each one to a .csv in the same folder.
# If a .rds cannot be converted, delete the copied .rds.
# =========================

# Define source and destination folders 
src_dir <- file.path(project_root, "02_results")
dst_dir <- file.path(project_root, "02_resultsv2")

# Basic checks
if (!dir.exists(src_dir)) {
  stop("Source folder does not exist: ", src_dir)
}

# Copy 02_results to 02_resultsv2 
if (dir.exists(dst_dir)) {
  unlink(dst_dir, recursive = TRUE, force = TRUE)
}

dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)

copy_ok <- file.copy(
  from = list.files(src_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE),
  to = dst_dir,
  recursive = TRUE,
  copy.mode = TRUE,
  copy.date = TRUE
)

if (!all(copy_ok)) {
  warning("Some files/folders may not have copied successfully.")
} else {
  message("Copy completed: ", src_dir, " -> ", dst_dir)
}

# Helper: convert different R object types into a data frame
object_to_df <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  
  if (is.matrix(x)) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  
  if (is.vector(x) || is.factor(x)) {
    return(data.frame(value = x, stringsAsFactors = FALSE))
  }
  
  if (is.list(x)) {
    # Case 1: list of equal-length atomic vectors -> data frame
    lens <- lengths(x)
    atomic_like <- vapply(x, function(z) is.atomic(z) || is.factor(z), logical(1))
    
    if (length(x) > 0 && all(atomic_like) && length(unique(lens)) == 1) {
      return(as.data.frame(x, stringsAsFactors = FALSE))
    }
    
    # Case 2: named list of uneven lengths -> pad to longest length
    if (length(x) > 0 && all(atomic_like)) {
      max_len <- max(lens)
      padded <- lapply(x, function(z) {
        length(z) <- max_len
        z
      })
      return(as.data.frame(padded, stringsAsFactors = FALSE))
    }
  }
  
  # Unsupported structure
  return(NULL)
}

# Find the top-level folders inside 02_resultsv2
top_level_dirs <- list.dirs(dst_dir, full.names = TRUE, recursive = FALSE)

if (length(top_level_dirs) == 0) {
  warning("No top-level folders found inside: ", dst_dir)
} else {
  message("Top-level folders found:")
  print(basename(top_level_dirs))
}

# Find all .rds files under top-level folders
rds_files <- unlist(
  lapply(top_level_dirs, function(folder) {
    list.files(
      path = folder,
      pattern = "\\.rds$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
  }),
  use.names = FALSE
)

rds_files <- unique(rds_files)

if (length(rds_files) == 0) {
  message("No .rds files found anywhere under ", dst_dir)
} else {
  message("Found ", length(rds_files), " .rds file(s). Beginning conversion to .csv ...")
}

# Convert each .rds to .csv in the same folder
conversion_log <- data.frame(
  rds_file = character(),
  csv_file = character(),
  status = character(),
  stringsAsFactors = FALSE
)

for (rds_path in rds_files) {
  csv_path <- sub("\\.rds$", ".csv", rds_path, ignore.case = TRUE)
  
  result <- tryCatch({
    obj <- readRDS(rds_path)
    df <- object_to_df(obj)
    
    if (is.null(df)) {

      unlink(rds_path, force = TRUE)
      
      data.frame(
        rds_file = rds_path,
        csv_file = NA_character_,
        status = "RDS DELETED: unsupported R object for CSV conversion",
        stringsAsFactors = FALSE
      )
    } else {
      write.csv(df, csv_path, row.names = FALSE, na = "")
      
      unlink(rds_path, force = TRUE)
      
      data.frame(
        rds_file = rds_path,
        csv_file = csv_path,
        status = "OK: converted to CSV and deleted RDS",
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {

    unlink(rds_path, force = TRUE)
    
    data.frame(
      rds_file = rds_path,
      csv_file = NA_character_,
      status = paste("RDS DELETED AFTER ERROR:", conditionMessage(e)),
      stringsAsFactors = FALSE
    )
  })
  
  conversion_log <- rbind(conversion_log, result)
}

# Print summary
message("Conversion finished.")
print(conversion_log)

# Save the log
log_path <- file.path(dst_dir, "rds_to_csv_conversion_log.csv")
write.csv(conversion_log, log_path, row.names = FALSE)
message("Conversion log saved to: ", log_path)