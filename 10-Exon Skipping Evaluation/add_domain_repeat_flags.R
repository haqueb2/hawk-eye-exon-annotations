#!/usr/bin/env Rscript

# -------------------------------------------------------------
# HAWK-EYE: Add protein-level domain/repeat counts and flags
# -------------------------------------------------------------
#
# Purpose:
#   This script adds protein-level UniProt domain/repeat counts to an
#   exon-level HAWK-EYE annotation file using a supplied lookup table.
#   It does NOT query the UniProt API.
#
# Required input files:
#   1) input_file:
#        Exon-level HAWK-EYE CSV containing at least:
#          - Uniprot_ID
#          - Exon_Specific_Count
#          - Exon_Specific_Count_Repeat
#
#   2) lookup_file:
#        Protein-level lookup CSV containing one row per UniProt ID:
#          - Uniprot_ID
#          - Total_Protein_Domain_Count
#          - Total_Protein_Repeat_Count
#
# Output:
#   Exon-level HAWK-EYE CSV with lookup counts and flags:
#          - Total_All_Protein_Domains
#          - Total_All_Protein_Repeat_Domains
#          - OnlyDomains
#          - OnlyRepeats
#          - MultipleDomains
#          - MultipleRepeats
#          - MultipleDomainsANDRepeats
#
# Example:
#   Rscript OnlyDomainandMultipleFilter_lookup_table.R \
#     input_file="HAWKEYE_Database_Domain_Repeat_Counts.csv" \
#     lookup_file="lookup/domain_repeat_lookup.csv" \
#     output_file="HAWKEYE_Database_Domain_Repeat_Counts_Flagged.csv"
#
# -------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# -------------------------------------------------------------
# 1. User inputs with command-line overrides
# -------------------------------------------------------------

input_file  = "HAWKEYE_Database_Domain_Repeat_Counts.csv"
lookup_file = "lookup/domain_repeat_lookup.csv"
output_file = "HAWKEYE_Database_Domain_Repeat_Counts_Flagged.csv"

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (a in args) {
    if (grepl("=", a)) {
      parts = strsplit(a, "=", fixed = TRUE)[[1]]
      key = parts[1]
      value = paste(parts[-1], collapse = "=")
      value = gsub('^"|"$', "", value)
      value = gsub("^'|'$", "", value)
      if (key %in% c("input_file", "lookup_file", "output_file")) {
        assign(key, value)
      } else {
        warning("Ignoring unknown argument: ", key)
      }
    }
  }
}

# -------------------------------------------------------------
# 2. Helper functions
# -------------------------------------------------------------

check_required_columns = function(df, required_cols, file_label) {
  missing_cols = setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      file_label, " is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

sum_semicolon_numbers = function(x) {
  if (is.na(x) || x == "" || x == "NA") {
    return(NA_real_)
  }

  parts = unlist(strsplit(as.character(x), ";", fixed = TRUE))
  parts = trimws(parts)
  nums = suppressWarnings(as.numeric(parts))

  if (all(is.na(nums))) {
    return(NA_real_)
  }

  sum(nums, na.rm = TRUE)
}

# -------------------------------------------------------------
# 3. Load input files
# -------------------------------------------------------------

hawkeye = read_csv(input_file, show_col_types = FALSE)
lookup = read_csv(lookup_file, show_col_types = FALSE)

check_required_columns(
  hawkeye,
  c("Uniprot_ID", "Exon_Specific_Count", "Exon_Specific_Count_Repeat"),
  "input_file"
)

check_required_columns(
  lookup,
  c("Uniprot_ID", "Total_Protein_Domain_Count", "Total_Protein_Repeat_Count"),
  "lookup_file"
)

# Keep exactly one lookup row per UniProt ID.
# If duplicates are present, keep the first and warn the user.
duplicate_ids = lookup$Uniprot_ID[duplicated(lookup$Uniprot_ID)]
if (length(duplicate_ids) > 0) {
  warning(
    "Duplicate UniProt IDs found in lookup table. Keeping the first row for each ID. Example duplicate(s): ",
    paste(head(unique(duplicate_ids), 10), collapse = ", ")
  )
}

lookup_clean = lookup %>%
  mutate(
    Uniprot_ID = as.character(Uniprot_ID),
    Total_Protein_Domain_Count = suppressWarnings(as.numeric(Total_Protein_Domain_Count)),
    Total_Protein_Repeat_Count = suppressWarnings(as.numeric(Total_Protein_Repeat_Count))
  ) %>%
  distinct(Uniprot_ID, .keep_all = TRUE) %>%
  select(
    Uniprot_ID,
    Total_Protein_Domain_Count,
    Total_Protein_Repeat_Count
  )

# -------------------------------------------------------------
# 4. Join lookup counts and calculate flags
# -------------------------------------------------------------

hawkeye_out = hawkeye %>%
  mutate(Uniprot_ID = as.character(Uniprot_ID)) %>%
  select(
    -any_of(c(
      "Total_Protein_Domain_Count",
      "Total_Protein_Repeat_Count",
      "Total_All_Protein_Domains",
      "Total_All_Protein_Repeat_Domains",
      "OnlyDomains",
      "OnlyRepeats",
      "MultipleDomains",
      "MultipleRepeats",
      "MultipleDomainsANDRepeats"
    ))
  ) %>%
  left_join(lookup_clean, by = "Uniprot_ID") %>%
  mutate(
    # Retain old output column names for compatibility with previous versions
    Total_All_Protein_Domains = Total_Protein_Domain_Count,
    Total_All_Protein_Repeat_Domains = Total_Protein_Repeat_Count,

    # Exon-specific sums from semicolon-separated count columns
    Exon_Specific_Domain_Sum = vapply(
      Exon_Specific_Count,
      sum_semicolon_numbers,
      numeric(1)
    ),
    Exon_Specific_Repeat_Sum = vapply(
      Exon_Specific_Count_Repeat,
      sum_semicolon_numbers,
      numeric(1)
    ),

    # Y when this exon accounts for all protein-level UniProt domains
    OnlyDomains = case_when(
      is.na(Total_All_Protein_Domains) | is.na(Exon_Specific_Domain_Sum) ~ NA_character_,
      Total_All_Protein_Domains == 0 & Exon_Specific_Domain_Sum == 0 ~ NA_character_,
      Total_All_Protein_Domains == Exon_Specific_Domain_Sum ~ "Y",
      TRUE ~ NA_character_
    ),

    # Y when this exon accounts for all protein-level UniProt repeats
    OnlyRepeats = case_when(
      is.na(Total_All_Protein_Repeat_Domains) | is.na(Exon_Specific_Repeat_Sum) ~ NA_character_,
      Total_All_Protein_Repeat_Domains == 0 & Exon_Specific_Repeat_Sum == 0 ~ NA_character_,
      Total_All_Protein_Repeat_Domains == Exon_Specific_Repeat_Sum ~ "Y",
      TRUE ~ NA_character_
    ),

    # Y when this exon overlaps more than one UniProt domain
    MultipleDomains = case_when(
      is.na(Exon_Specific_Domain_Sum) ~ NA_character_,
      Exon_Specific_Domain_Sum > 1 ~ "Y",
      TRUE ~ NA_character_
    ),

    # Y when this exon overlaps more than one UniProt repeat
    MultipleRepeats = case_when(
      is.na(Exon_Specific_Repeat_Sum) ~ NA_character_,
      Exon_Specific_Repeat_Sum > 1 ~ "Y",
      TRUE ~ NA_character_
    ),

    # Y when this exon overlaps at least one domain and at least one repeat
    MultipleDomainsANDRepeats = case_when(
      !is.na(Exon_Specific_Domain_Sum) &
        !is.na(Exon_Specific_Repeat_Sum) &
        Exon_Specific_Domain_Sum >= 1 &
        Exon_Specific_Repeat_Sum >= 1 ~ "Y",
      TRUE ~ "N"
    )
  ) %>%
  select(
    -Exon_Specific_Domain_Sum,
    -Exon_Specific_Repeat_Sum
  )

# -------------------------------------------------------------
# 5. Report unmatched UniProt IDs
# -------------------------------------------------------------

unmatched_ids = hawkeye_out %>%
  filter(!is.na(Uniprot_ID), is.na(Total_Protein_Domain_Count), is.na(Total_Protein_Repeat_Count)) %>%
  distinct(Uniprot_ID) %>%
  pull(Uniprot_ID)

if (length(unmatched_ids) > 0) {
  warning(
    length(unmatched_ids),
    " UniProt ID(s) in the HAWK-EYE input were not found in the lookup table. Example unmatched ID(s): ",
    paste(head(unmatched_ids, 10), collapse = ", ")
  )
}

# -------------------------------------------------------------
# 6. Write output
# -------------------------------------------------------------

# Replace blank values in newly-created flag columns with "NA"
flag_cols = c(
  "OnlyDomains",
  "OnlyRepeats",
  "MultipleDomains",
  "MultipleRepeats",
  "MultipleDomainsANDRepeats"
)

hawkeye_out = hawkeye_out %>%
  mutate(
    across(
      all_of(flag_cols),
      ~ ifelse(is.na(.) | trimws(as.character(.)) == "", "NA", as.character(.))
    )
  )

write_csv(hawkeye_out, output_file, na = "")

message("Done. Wrote: ", output_file)
