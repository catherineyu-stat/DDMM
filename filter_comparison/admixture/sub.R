

# --- 1. Define Output Directories ---
# It's good practice to define these at the top for easy configuration.

# Base directory for all output
output_base_dir <- "slurm_output"

# Subdirectories for different file types
scripts_dir <- file.path(output_base_dir, "submission_scripts")
log_dir     <- file.path(output_base_dir, "log_files")
err_dir     <- file.path(output_base_dir, "error_files")
rout_dir    <- file.path(output_base_dir, "rout_files")

# --- 2. Create Directories ---
# This ensures the folders exist before we try to write files to them.
# The 'recursive = TRUE' argument allows it to create nested directories (e.g., slurm_output/log_files).
# The 'showWarnings = FALSE' argument prevents messages if the directory already exists.
dir.create(scripts_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(err_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(rout_dir,    showWarnings = FALSE, recursive = TRUE)


# Define your tasks (replace with your actual task list)
tasks <- 1:50

# --- 3. Loop and Generate Scripts with Organized Paths ---
for(state in tasks)
{
  task_name <- paste('c', state, sep='')
  jobname   <- task_name
  
  # --- Construct full paths for each file ---
  # Use file.path() to correctly join the directory and filename.
  qsub_file <- file.path(scripts_dir, paste0(task_name, '.txt'))
  logname   <- file.path(log_dir,     paste0(task_name, 'state.log'))
  errfile   <- file.path(err_dir,     paste0(task_name, 'state.err'))
  
  # Define the full path for the .Rout file
  # Note: The R CMD BATCH command needs the full path for its output file.
  rout_file <- file.path(rout_dir, paste0(task_name, '-state.Rout'))
  
  # --- Write the Slurm submission script ---
  sink(qsub_file)
  
  cat('#!/bin/bash\n')
  cat('\n')
  cat(sprintf("#SBATCH -J  %s\n", jobname))
  cat('#SBATCH -N 1 -c 1\n')
  cat('#SBATCH --time=72:59:59\n')
  cat('\n')
  
  # Use the full, absolute paths for log and error files.
  # We use normalizePath to ensure Slurm gets an unambiguous absolute path.
  cat(sprintf("#SBATCH -o %s\n", normalizePath(logname)))
  cat(sprintf("#SBATCH -e %s\n", normalizePath(errfile)))
  cat('\n')
  
  # The 'cd' command is still useful if your R script needs to read other files
  # from the directory where you launched this generator script.
  dir <- getwd()
  cat(paste('cd', dir, '\n'))
  cat('\n')
  
  # The R command now writes its .Rout file to the specified directory.
  # We also use normalizePath() here for robustness.
  cat(sprintf("R CMD BATCH --vanilla '--args %d' compare_new_extend.R %s", state, normalizePath(rout_file)))
  
  sink()
  
  # --- Submit the job to Slurm ---
  com <- paste('sbatch', qsub_file)
  system(com)
  
  cat(paste("Submitted job for state:", state, "with script:", qsub_file, "\n"))
}

  
  
