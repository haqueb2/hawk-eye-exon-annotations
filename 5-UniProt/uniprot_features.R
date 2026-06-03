library(dplyr)
library(readr)
library(stringr)
library(httr)
library(jsonlite)
library(tidyr)

# ---- User inputs ----

input_file = "HAWKEYE_Database_PhyloP_Exon_Conservation.csv"
idmapping_file = "uniprot_files/HUMAN_9606_idmapping.dat"
output_file = "HAWKEYE_Database_UniProt_Features.csv"

# ---- Override with command-line key=value args (e.g. input_file="path") ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (a in args) {
    if (grepl("=", a)) {
      try(eval(parse(text = a)), silent = TRUE)
    }
  }
}

# ---- Load exon data ----

exonData = read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    Exon_Start = as.integer(Exon_Start),
    Exon_End = as.integer(Exon_End),
    Exon_Length = as.integer(Exon_Length),
    Exon_Frame = as.integer(Exon_Frame)
  ) %>%
  arrange(Gene, Transcript_ID, Exon)

# ---- Map gene symbols to first UniProt ID ----

idMap = read_tsv(
  idmapping_file,
  col_names = c("Uniprot_ID", "ID_Type", "Gene"),
  show_col_types = FALSE
)

geneIdMap = idMap %>%
  filter(ID_Type == "Gene_Name") %>%
  group_by(Gene) %>%
  summarise(
    Uniprot_ID = dplyr::first(Uniprot_ID),
    .groups = "drop"
  )

exonData = exonData %>%
  left_join(geneIdMap, by = "Gene")

# ---- Calculate amino acid positions per exon ----

exonData = exonData %>%
  group_by(Gene, Transcript_ID) %>%
  mutate(
    Codon_Count = case_when(
      is.na(Exon_Sequence) ~ NA_integer_,
      Exon_Frame == 1 ~ floor((nchar(Exon_Sequence) - 2) / 3) + 1,
      Exon_Frame == 2 ~ floor((nchar(Exon_Sequence) - 1) / 3) + 1,
      Exon_Frame == 0 ~ floor(nchar(Exon_Sequence) / 3),
      TRUE ~ NA_integer_
    ),
    AA_Start = if_else(
      !is.na(Codon_Count),
      cumsum(replace_na(Codon_Count, 0)) - Codon_Count + 1,
      NA_integer_
    ),
    AA_End = if_else(
      !is.na(Codon_Count),
      AA_Start + Codon_Count - 1,
      NA_integer_
    )
  ) %>%
  mutate(
    AA_End = if_else(
      lead(Exon_Frame) %in% c(1, 2) & !is.na(AA_End),
      AA_End + 1,
      AA_End
    ),
    AA_Pos = if_else(
      !is.na(AA_Start) & !is.na(AA_End),
      paste0("[", AA_Start, "-", AA_End, "]"),
      NA_character_
    )
  ) %>%
  ungroup()

# ---- Helper functions ----

collapse_unique = function(x) {
  x = x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(unique(x), collapse = ";")
}

get_uniprot_features = function(uniprot_id) {
  request_url = paste0(
    "https://rest.uniprot.org/uniprotkb/",
    uniprot_id,
    ".json?fields=ft_region%2Cft_zn_fing%2Cft_compbias%2Cft_motif"
  )
  
  response = GET(request_url, accept("application/json"))
  
  if (status_code(response) != 200) {
    warning(paste("UniProt request failed for", uniprot_id))
    return(NULL)
  }
  
  json_content = fromJSON(rawToChar(response$content), flatten = TRUE)
  
  if (is.null(json_content$features) || length(json_content$features) == 0) {
    return(NULL)
  }
  
  json_content$features %>%
    as_tibble() %>%
    filter(type %in% c("Region", "Compositional bias", "Zinc finger", "Motif")) %>%
    transmute(
      Uniprot_ID = uniprot_id,
      Feature_Type = type,
      Feature_Description = description,
      Feature_Start = location.start.value,
      Feature_End = location.end.value
    )
}

# ---- Download UniProt features ----

uniprot_features = bind_rows(
  lapply(
    unique(na.omit(exonData$Uniprot_ID)),
    get_uniprot_features
  )
)

# ---- Annotate exon-level UniProt feature overlaps ----

exonData = exonData %>%
  mutate(
    Region = NA_character_,
    Compositional_Bias = NA_character_,
    Motif = NA_character_,
    Zinc_Finger = NA_character_
  )

for (i in seq_len(nrow(exonData))) {
  if (
    is.na(exonData$Uniprot_ID[i]) ||
    is.na(exonData$AA_Start[i]) ||
    is.na(exonData$AA_End[i]) ||
    nrow(uniprot_features) == 0
  ) {
    next
  }
  
  overlaps = uniprot_features %>%
    filter(
      Uniprot_ID == exonData$Uniprot_ID[i],
      Feature_Start <= exonData$AA_End[i],
      Feature_End >= exonData$AA_Start[i]
    )
  
  exonData$Region[i] = collapse_unique(
    overlaps$Feature_Description[overlaps$Feature_Type == "Region"]
  )
  
  exonData$Compositional_Bias[i] = collapse_unique(
    overlaps$Feature_Description[overlaps$Feature_Type == "Compositional bias"]
  )
  
  exonData$Motif[i] = collapse_unique(
    overlaps$Feature_Description[overlaps$Feature_Type == "Motif"]
  )
  
  exonData$Zinc_Finger[i] = collapse_unique(
    overlaps$Feature_Description[overlaps$Feature_Type == "Zinc finger"]
  )
}

# ---- Remove intermediate columns if desired ----

exonData = exonData %>%
  select(
    -Codon_Count,
    -AA_Start,
    -AA_End
  )

# ---- Export ----

write_csv(exonData, output_file)