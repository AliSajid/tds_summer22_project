# Process the Gene Count Data

library(tidyverse)
library(org.Hs.eg.db)

metadata <- read_csv("data/dataset_metadata.csv")

raw_count_data <- read_csv("raw/GSE80655-Count-Data.csv.csv") |>
    dplyr::select(gene_id, any_of(metadata$sample_id))

gene_metadata <-
    AnnotationDbi::select(
        org.Hs.eg.db,
        keys = raw_count_data$gene_id,
        keytype = "ENSEMBL",
        columns = c("ENSEMBL", "GENENAME", "GENETYPE", "SYMBOL", "ENTREZID")
    ) |>
    dplyr::rename(
        gene_id = ENSEMBL,
        gene_name = GENENAME,
        gene_type = GENETYPE,
        hgnc_symbol = SYMBOL,
        entrez_id = ENTREZID
    )



count_data <- gene_metadata |>
    inner_join(raw_count_data, by = "gene_id") |>
    write_csv("data/dataset_count_data.csv")
