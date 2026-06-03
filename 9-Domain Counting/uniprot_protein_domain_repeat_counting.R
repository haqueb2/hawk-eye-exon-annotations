library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# ---- User inputs ----

input_file = "HAWKEYE_Database_PSI.csv"
output_file = "HAWKEYE_Database_Domain_Repeat_Counts.csv"

# ---- Override with command-line key=value args (e.g. input_file="path") ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (a in args) {
    if (grepl("=", a)) {
      try(eval(parse(text = a)), silent = TRUE)
    }
  }
}

# ---- Load data ----

hawkeye = read_csv(input_file, show_col_types = FALSE)

# ---- Helper functions ----

clean_na = function(x) {
  x = as.character(x)
  x[x %in% c("", "NA", "Na", "na")] = NA_character_
  x
}

collapse_values = function(x) {
  x = x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = "; ")
}

# ---- Add Exon_ID column ----

hawkeye = hawkeye %>%
  mutate(
    Exon_ID = paste0(Gene, " exon ", Exon)
  ) %>%
  relocate(Exon_ID, .after = Gene)

# ---- Domain counts ----

domain_list = hawkeye %>%
  select(Transcript_ID, Exon_ID, Functional_Domains) %>%
  mutate(Functional_Domains = clean_na(Functional_Domains)) %>%
  separate_rows(Functional_Domains, sep = "\\|") %>%
  mutate(
    Functional_Domains = str_remove(Functional_Domains, ";.*"),
    Functional_Domains = trimws(Functional_Domains),
    Join_Key = tolower(str_remove(Functional_Domains, "\\s+\\d+$")),
    Domain_Number_Str = str_extract(Functional_Domains, "\\s+\\d+$"),
    Domain_Index = if_else(
      is.na(Domain_Number_Str),
      1,
      as.numeric(trimws(Domain_Number_Str))
    )
  ) %>%
  filter(!is.na(Join_Key), Join_Key != "na", Join_Key != "")

domain_totals = domain_list %>%
  group_by(Transcript_ID, Join_Key) %>%
  summarise(
    Total_In_Protein = max(Domain_Index, na.rm = TRUE),
    .groups = "drop"
  )

domain_summary = domain_list %>%
  left_join(domain_totals, by = c("Transcript_ID", "Join_Key")) %>%
  group_by(Transcript_ID, Exon_ID, Join_Key) %>%
  summarise(
    Local_Exon_Count = max(Domain_Index, na.rm = TRUE) - min(Domain_Index, na.rm = TRUE) + 1,
    Total_In_Protein = dplyr::first(Total_In_Protein),
    Exon_Contribution = round(Local_Exon_Count / Total_In_Protein, 2),
    .groups = "drop"
  ) %>%
  group_by(Transcript_ID, Exon_ID) %>%
  summarise(
    Encoded_Domains = collapse_values(Join_Key),
    Exon_Specific_Count = collapse_values(as.character(Local_Exon_Count)),
    Total_In_Protein = collapse_values(as.character(Total_In_Protein)),
    Exon_Contribution = collapse_values(as.character(Exon_Contribution)),
    .groups = "drop"
  )

# ---- Repeat counts ----

repeat_list = hawkeye %>%
  select(Transcript_ID, Exon_ID, Repeat) %>%
  mutate(Repeat = clean_na(Repeat)) %>%
  separate_rows(Repeat, sep = "\\|") %>%
  mutate(
    Repeat = str_remove(Repeat, ";.*"),
    Repeat = if_else(str_detect(Repeat, "^\\d+$"), paste("Repeat", Repeat), Repeat),
    Repeat = trimws(Repeat),
    Join_Key = tolower(str_remove(Repeat, "\\s+\\d+$|(?<=\\d-)\\d+$")),
    Repeat_Number_Str = str_extract(Repeat, "\\s+\\d+$|(?<=\\d-)\\d+$"),
    Repeat_Index = if_else(
      is.na(Repeat_Number_Str),
      1,
      as.numeric(trimws(Repeat_Number_Str))
    )
  ) %>%
  filter(!is.na(Join_Key), Join_Key != "na", Join_Key != "")

repeat_totals = repeat_list %>%
  group_by(Transcript_ID, Join_Key) %>%
  summarise(
    Total_Repeats_In_Protein = max(Repeat_Index, na.rm = TRUE),
    .groups = "drop"
  )

repeat_summary = repeat_list %>%
  left_join(repeat_totals, by = c("Transcript_ID", "Join_Key")) %>%
  group_by(Transcript_ID, Exon_ID, Join_Key) %>%
  summarise(
    Local_Exon_Count_Repeat = max(Repeat_Index, na.rm = TRUE) - min(Repeat_Index, na.rm = TRUE) + 1,
    Total_Repeats_In_Protein = dplyr::first(Total_Repeats_In_Protein),
    Exon_Contribution_Repeat = round(Local_Exon_Count_Repeat / Total_Repeats_In_Protein, 2),
    .groups = "drop"
  ) %>%
  group_by(Transcript_ID, Exon_ID) %>%
  summarise(
    Encoded_Repeats = collapse_values(Join_Key),
    Exon_Specific_Count_Repeat = collapse_values(as.character(Local_Exon_Count_Repeat)),
    Total_Repeats_In_Protein = collapse_values(as.character(Total_Repeats_In_Protein)),
    Exon_Contribution_Repeat = collapse_values(as.character(Exon_Contribution_Repeat)),
    .groups = "drop"
  )

# ---- Merge annotations ----

hawkeye = hawkeye %>%
  left_join(domain_summary, by = c("Transcript_ID", "Exon_ID")) %>%
  left_join(repeat_summary, by = c("Transcript_ID", "Exon_ID"))

# ---- Export ----

write_csv(hawkeye, output_file)