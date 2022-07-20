# Setup Script


create_dir <- function(dir) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

data <- "data"
results <- "results"
figures <- "figures"

dirs <- c(data, results, figures)
sapply(dirs, create_dir)
