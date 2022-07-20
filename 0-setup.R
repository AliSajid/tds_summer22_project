# Setup Script


create_dir <- function(dir) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

data <- "data"
results <- "results"
dge_dir <- file.path(results, "dge")
regions <- file.path(dge_dir, c(nAcc = "nAcc", DLPFC = "DLPFC", AnCg = "AnCg"))
figures <- "figures"

dirs <- c(data, regions, figures)
sapply(dirs, create_dir)
