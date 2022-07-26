# Get drugfindR Results

library(tidyverse)
library(drugfindR)
library(furrr)

investigate_signature_custom <-
    function(expr,
             prop = 0.95,
             similarity_threshold = 0.2,
             discordant = FALSE,
             gene_column = "hgnc_symbol",
             logfc_column = "logFC",
             pval_column = "PValue",
             source_name = "Input") {
        expr_signature <-
            expr %>% prepare_signature(
                gene_column = gene_column,
                logfc_column = logfc_column,
                pval_column = pval_column
            )

        signature_id <- unique(expr_signature$signatureID)

        up_threshold <-
            quantile(expr_signature$Value_LogDiffExp, prop)

        down_threshold <-
            abs(quantile(expr_signature$Value_LogDiffExp, {1 - prop}))

        filtered_up <-
            expr_signature %>% filter_signature(direction = "up",
                                                threshold = up_threshold)
        filtered_down <-
            expr_signature %>% filter_signature(direction = "down",
                                                threshold = down_threshold)
        concordant_up <-
            filtered_up %>% get_concordants(library = "CP")

        concordant_down <-
            filtered_down %>% get_concordants(library = "CP")

        consensus_targets <- consensus_concordants(
            concordant_up,
            concordant_down,
            paired = TRUE,
            cell_line = NULL,
            discordant = discordant,
            cutoff = similarity_threshold
        )

        augmented <-
            consensus_targets %>% dplyr::mutate(
                SourceSignature = signature_id,
                Source = source_name,
                SourceCellLine = "NA",
                SourceTime = "NA",
                UpThreshold = up_threshold,
                UpThresholdGeneCount = nrow(filtered_up),
                DownThreshold = down_threshold,
                DownThresholdGeneCount = nrow(filtered_down),
                Quantile = prop
            ) %>% dplyr::select(
                .data$Source,
                .data$Target,
                .data$Similarity,
                .data$SourceSignature,
                .data$SourceCellLine,
                dplyr::any_of(c("SourceConcentration")),
                .data$SourceTime,
                .data$TargetSignature,
                .data$TargetCellLine,
                dplyr::any_of(c("TargetConcentration")),
                .data$TargetTime,
                .data$UpThreshold,
                .data$UpThresholdGeneCount,
                .data$DownThreshold,
                .data$DownThresholdGeneCount,
                .data$Quantile
            )

        augmented
    }

process_signature <- function(signature, data_source = NULL) {
    quantile <- c(0.90, 0.95, 0.97)
    discordant <- c(TRUE, FALSE)


    parameters <-
        expand_grid(quantile,
                    discordant)

    results <- parameters |>
        future_pmap_dfr(
            ~ investigate_signature_custom(
                expr = signature,
                prop = ..1,
                discordant = ..2,
                source_name = data_source
            )
        )

    results

}


filenames <- list.files("results/dge/", ".csv", recursive = TRUE) |>
    str_remove(".*/") |>
    str_remove("-dge.csv") |>
    keep( ~ str_detect(.x, "ALL", negate = TRUE))

filepaths <-
    list.files("results/dge/",
               ".csv",
               recursive = TRUE,
               full.names = TRUE) |>
    keep( ~ str_detect(.x, "ALL", negate = TRUE))

plan(multisession)

if (!file.exists("results/full_dataset.RDS")) {
    full_dataset <- filepaths |>
        set_names(filenames) |>
        map(~ read_csv(.x, show_col_types = FALSE)) |>
        future_imap(~ process_signature(.x, .y))

    saveRDS(full_dataset, file = "results/full_dataset.RDS")
} else {
    full_dataset <- readRDS("results/full_dataset.RDS")
}

