# Differential Gene Expression Analysis

library(tidyverse)
library(edgeR)
library(org.Hs.eg.db)

get_dge <- function(region) {
    metadata <- read_csv("data/dataset_metadata.csv") |>
        filter(brain_region == region)
    counts <- read_csv("data/dataset_count_data.csv") |>
        dplyr::select(gene_id:entrez_id, any_of(metadata$sample_id))

    dge <-
        DGEList(
            counts = counts[, 6:ncol(counts)],
            genes = counts[, 1:5],
            samples = metadata,
            group = metadata$clinical_diagnosis
        )

    keep <- filterByExpr(dge)

    dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

    diag <- metadata$clinical_diagnosis

    design <- model.matrix(~ 0 + diag)

    contrasts <- makeContrasts(
        contrasts = c(BPD = "diagBPD - diagCTL",
                      MDD = "diagMDD - diagCTL",
                      SCZ = "diagSCZ - diagCTL"),
        levels = colnames(design)
    )

    dge_filtered <- calcNormFactors(dge_filtered)

    dge_filtered <- estimateDisp(dge_filtered, design)

    fit <- glmQLFit(dge_filtered, design)

    qlf_bpd <- glmQLFTest(fit, contrast = contrasts[, 1])

    qlf_mdd <- glmQLFTest(fit, contrast = contrasts[, 2])

    qlf_scz <- glmQLFTest(fit, contrast = contrasts[, 3])

    qlf_all <- glmQLFTest(fit, contrast = contrasts)

    out <- list(
        BPD = topTags(qlf_bpd, n = Inf)$table,
        MDD = topTags(qlf_mdd, n = Inf)$table,
        SCZ = topTags(qlf_scz, n = Inf)$table,
        ALL = topTags(qlf_all, n = Inf)$table
    )
}


write_results <- function(dataset, name) {

    split_name <- str_split(name, "\\.", simplify = TRUE)
    region <- split_name[1]
    diff_name <- split_name[2]

    outfile <- file.path(
        "results",
        "dge",
        region,
        str_glue("{region}-{diff_name}-dge.csv")
    )

    write_csv(dataset, outfile)

}


regions <- c(nAcc = "nAcc", DLPFC = "DLPFC", AnCg = "AnCg")

diff_genes <- regions |>
    map(get_dge) |>
    unlist(FALSE) %>%
    map2(names(.), write_results)
