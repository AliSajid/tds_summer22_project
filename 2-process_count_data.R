# Process the Gene Count Data

library(tidyverse)

metadata <- read_csv("data/dataset_metadata.csv")

count_data <- read_csv("raw/GSE80655-Count-Data.csv.csv") |>
    select(gene_id, any_of(metadata$sample_id)) |>
    write_csv("data/dataset_count_data.csv")
