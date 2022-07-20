# Process Sample Metadata
#

library(tidyverse)

metadata <- read_csv("raw/GSE80655-Sample-Metadata.csv.csv") |>
    group_by(`Sample-ID`) |>
    filter(Bases == max(Bases)) |>
    ungroup() |>
    rename_with(~ str_replace_all(.x, "\\s", "_")) |>
    select(`Sample-ID`, Sample_Name, Run, clinical_diagnosis,
           brain_region, gender, age_at_death,
           Brain_pH, `post-mortem_interval`, Ethnicity) |>
    rename_with(~ str_replace_all(.x, "-", "_")) |>
    rename_with(~ str_to_lower(.x)) |>
    mutate(clinical_diagnosis = case_when(
        clinical_diagnosis == "Control" ~ "CTL",
        clinical_diagnosis == "Schizophrenia" ~ "SCZ",
        clinical_diagnosis == "Bipolar Disorder" ~ "BPD",
        clinical_diagnosis == "Major Depression" ~ "MDD"
    )) |>
    write_csv("data/dataset_metadata.csv")

