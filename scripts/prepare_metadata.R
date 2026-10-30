#!/usr/bin/env Rscript

library(readxl)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)

# Hent kommandolinjeparametre
args <- commandArgs(trailingOnly = TRUE)
excel_input <- args[1]
tsv_input <- args[2]
tsv_output <- args[3]

# Les inn filer
requests <- read_excel(excel_input)
hav <- read_tsv(tsv_input,
                col_types = cols(.default = "c"))

# Sørg for at manglende verdier blir tom streng
requests <- requests %>%
  mutate(across(everything(), ~replace_na(as.character(.), "")))

hav <- hav %>%
  mutate(across(everything(), ~replace_na(as.character(.), "")))

# Data fra Excel-filen
metadata_requests <- requests %>%
  transmute(
    id = SAMPLE_NUMBER,
    genotype = str_extract(Resultat, "(?<=HAV GENOTYPE\\s)\\S+"),
    variant = Outbreak_variant,
    sample_date = Year
  ) %>%
  mutate(
    genotype = replace_na(genotype, ""),
    sample_date = replace_na(sample_date, "")
  )

# Data fra TSV-filen
# Sjekk om SAMPLED_DATE kolonne finnes
has_sampled_date <- "SAMPLED_DATE" %in% names(hav)

metadata_hav <- hav %>%
  transmute(
    id = paste(SAMPLE_NUMBER, SekvensID, sep = "_"),
    genotype = str_extract(Resultat, "(?<=HAV GENOTYPE\\s)\\S+"),
    variant = Outbreak_variant,
    sample_date = if (has_sampled_date) {
      if_else(
        !is.na(SAMPLED_DATE) & SAMPLED_DATE != "",
        SAMPLED_DATE,
        ""
      )
    } else {
      ""
    }
  ) %>%
  mutate(
    genotype = replace_na(genotype, ""),
    sample_date = replace_na(sample_date, "")
  )

# Slå sammen
metadata <- bind_rows(
  metadata_requests,
  metadata_hav
) %>%
  mutate(across(everything(), ~replace_na(., "")))

# Skriv ut TSV
write_tsv(metadata, tsv_output, na = "")