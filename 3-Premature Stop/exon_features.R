library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(Biostrings)

# ---- User inputs ----

exon_file = "HAWKEYE_Database_ClinVar_Counts.csv"

ccds_sequence_file = "premature_stop/CCDS_nucleotide.current.fna"
ccds_to_sequence_file = "premature_stop/CCDS2Sequence.current.txt"
ccds_summary_file = "premature_stop/CCDS.current.txt"

output_file = "HAWKEYE_Database_Exon_Features.csv"

# ---- Override with command-line key=value args (e.g. exon_file="path") ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (a in args) {
    if (grepl("=", a)) {
      try(eval(parse(text = a)), silent = TRUE)
    }
  }
}

# ---- Load exon data ----

exonData = read_csv(exon_file, show_col_types = FALSE) %>%
  mutate(
    Exon_Start = as.integer(Exon_Start),
    Exon_End = as.integer(Exon_End),
    Exon_Length = as.numeric(Exon_Length),
    Exon_Frame = as.integer(Exon_Frame)
  ) %>%
  arrange(Gene, Transcript_ID, Exon)

# ---- Determine first coding exon and last exon ----

exonData = exonData %>%
  group_by(Gene, Transcript_ID) %>%
  mutate(
    First_Coding_Exon_Number = suppressWarnings(min(Exon[Exon_Frame %in% 0:2], na.rm = TRUE)),
    First_Coding_Exon_Number = if_else(is.infinite(First_Coding_Exon_Number), NA_real_, First_Coding_Exon_Number),
    
    First_Coding_Exon = if_else(
      !is.na(First_Coding_Exon_Number) & Exon == First_Coding_Exon_Number,
      "Y",
      "N"
    ),
    
    Last_Exon = if_else(
      Exon == max(Exon, na.rm = TRUE),
      "Y",
      "N"
    )
  ) %>%
  ungroup()

# ---- Load CCDS transcript mapping ----

ccds_info = read_tsv(ccds_to_sequence_file, show_col_types = FALSE) %>%
  mutate(
    Transcript_ID = str_remove(nucleotide_ID, "\\..*")
  ) %>%
  filter(Transcript_ID %in% unique(exonData$Transcript_ID)) %>%
  dplyr::select(Transcript_ID, CCDS_ID = ccds)

exonData = exonData %>%
  left_join(ccds_info, by = "Transcript_ID")

# ---- Load CCDS summary ----

ccds_summary = read_tsv(ccds_summary_file, show_col_types = FALSE) %>%
  filter(ccds_id %in% unique(exonData$CCDS_ID)) %>%
  dplyr::select(
    CCDS_ID = ccds_id,
    CCDS_Gene = gene,
    cds_strand,
    cds_from,
    cds_to
  )

exonData = exonData %>%
  left_join(ccds_summary, by = "CCDS_ID")

# ---- Identify coding start position ----

exonData = exonData %>%
  mutate(
    Start_Pos = case_when(
      First_Coding_Exon == "Y" & cds_strand == "+" ~ as.numeric(cds_from) + 1,
      First_Coding_Exon == "Y" & cds_strand == "-" ~ as.numeric(cds_to) + 1,
      TRUE ~ NA_real_
    )
  )

# ---- Adjust first coding exon length before sequence extraction ----

exonData = exonData %>%
  mutate(
    Exon_Length = case_when(
      
      # Positive strand:
      # CDS starts at Start_Pos and continues toward Exon_End
      First_Coding_Exon == "Y" &
        Strand == "+" &
        !is.na(Start_Pos) ~
        abs(Exon_End - Start_Pos + 1),
      
      # Negative strand:
      # CDS starts at Start_Pos and continues toward Exon_Start
      First_Coding_Exon == "Y" &
        Strand == "-" &
        !is.na(Start_Pos) ~
        abs(Start_Pos - Exon_Start + 1),
      
      TRUE ~ Exon_Length
    )
  )

# ---- Load CCDS coding sequences ----

seqCCDS = readDNAStringSet(ccds_sequence_file, format = "fasta")

seqCCDS_df = tibble(
  SequenceName = names(seqCCDS),
  Sequence = as.character(seqCCDS)
) %>%
  mutate(
    CCDS_ID = sapply(strsplit(SequenceName, "\\|"), function(x) x[1])
  ) %>%
  filter(CCDS_ID %in% unique(exonData$CCDS_ID)) %>%
  dplyr::select(CCDS_ID, Sequence)

exonData = exonData %>%
  left_join(seqCCDS_df, by = "CCDS_ID")

# ---- Extract coding exon sequence from CCDS sequence ----

exonData = exonData %>%
  group_by(Gene, Transcript_ID) %>%
  mutate(
    Coding_Exon_Length = if_else(Exon_Frame %in% 0:2, Exon_Length, NA_real_),
    
    CCDS_Start = if_else(
      Exon_Frame %in% 0:2,
      cumsum(replace_na(Coding_Exon_Length, 0)) - Coding_Exon_Length + 1,
      NA_real_
    ),
    
    CCDS_End = if_else(
      Exon_Frame %in% 0:2,
      CCDS_Start + Coding_Exon_Length - 1,
      NA_real_
    ),
    
    Exon_Sequence = if_else(
      Exon_Frame %in% 0:2 & !is.na(Sequence),
      substr(Sequence, CCDS_Start, CCDS_End),
      NA_character_
    )
  ) %>%
  ungroup()

# ---- Normalize sequences and correct exon length/CDS annotations ----

exonData = exonData %>%
  mutate(
    Exon_Sequence = trimws(Exon_Sequence),
    Exon_Sequence = na_if(Exon_Sequence, ""),
    Exon_Sequence = na_if(Exon_Sequence, "NA"),
    Exon_Sequence = na_if(Exon_Sequence, "Na"),
    Exon_Sequence = na_if(Exon_Sequence, "na")
  ) %>%
  group_by(Gene, Transcript_ID) %>%
  mutate(
    Last_Coding_Exon_Number = suppressWarnings(max(Exon[!is.na(Exon_Sequence)], na.rm = TRUE)),
    Last_Coding_Exon_Number = if_else(is.infinite(Last_Coding_Exon_Number), NA_real_, Last_Coding_Exon_Number),
    
    Last_Coding_Exon = if_else(
      !is.na(Exon_Sequence) & Exon == Last_Coding_Exon_Number,
      "Y",
      "N"
    ),
    
    Exon_Length = case_when(
      is.na(Exon_Sequence) ~ NA_real_,
      First_Coding_Exon == "Y" ~ as.numeric(nchar(Exon_Sequence)),
      Last_Coding_Exon == "Y" ~ as.numeric(nchar(Exon_Sequence)),
      TRUE ~ Exon_Length
    ),
    
    CDS_Total = sum(Exon_Length[!is.na(Exon_Sequence)], na.rm = TRUE),
    
    CDS = if_else(
      CDS_Total > 0 & !is.na(Exon_Length),
      Exon_Length / CDS_Total,
      NA_real_
    ),
    
    Coding_Exon_Length = if_else(
      !is.na(Exon_Sequence),
      Exon_Length,
      NA_real_
    )
  ) %>%
  ungroup()

# ---- Determine in-frame exons ----

exonData = exonData %>%
  mutate(
    In_Frame_Exon = if_else(
      !is.na(Exon_Length) & Exon_Length %% 3 == 0,
      "Y",
      "N"
    )
  )

# ---- Predict codon generated after exon skipping ----

exonData = exonData %>%
  group_by(Gene, Transcript_ID) %>%
  mutate(
    Coding_Length_Before_Exon = lag(cumsum(replace_na(Coding_Exon_Length, 0)), default = 0),
    Previous_Exon_Sequence = lag(Exon_Sequence),
    Next_Exon_Sequence = lead(Exon_Sequence),
    
    Skipped_Codon = case_when(
      In_Frame_Exon == "N" ~ NA_character_,
      Exon_Frame == -1 ~ NA_character_,
      First_Coding_Exon == "Y" ~ NA_character_,
      Last_Coding_Exon == "Y" ~ NA_character_,
      is.na(Previous_Exon_Sequence) ~ NA_character_,
      is.na(Next_Exon_Sequence) ~ NA_character_,
      
      Coding_Length_Before_Exon %% 3 == 0 ~ "0-0",
      
      Coding_Length_Before_Exon %% 3 == 1 ~ paste0(
        substr(
          Previous_Exon_Sequence,
          nchar(Previous_Exon_Sequence),
          nchar(Previous_Exon_Sequence)
        ),
        substr(Next_Exon_Sequence, 1, 2)
      ),
      
      Coding_Length_Before_Exon %% 3 == 2 ~ paste0(
        substr(
          Previous_Exon_Sequence,
          nchar(Previous_Exon_Sequence) - 1,
          nchar(Previous_Exon_Sequence)
        ),
        substr(Next_Exon_Sequence, 1, 1)
      ),
      
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# ---- Assess premature stop codons ----

exonData = exonData %>%
  mutate(
    Premature_Stop = case_when(
      is.na(Skipped_Codon) ~ NA_character_,
      Skipped_Codon %in% c("TAA", "TAG", "TGA") ~ "Y",
      TRUE ~ "N"
    )
  )

# ---- Remove intermediate columns ----

exonData = exonData %>%
  dplyr::select(
    -First_Coding_Exon_Number,
    -Last_Coding_Exon_Number,
    -CCDS_Gene,
    -cds_strand,
    -cds_from,
    -cds_to,
    -Sequence,
    -Coding_Exon_Length,
    -CCDS_Start,
    -CCDS_End,
    -Previous_Exon_Sequence,
    -Next_Exon_Sequence,
    -Coding_Length_Before_Exon,
    -CDS_Total
  )

# ---- Export ----

write_csv(exonData, output_file)