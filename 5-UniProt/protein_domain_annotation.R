library(dplyr)
library(readr)
library(stringr)
library(httr)
library(jsonlite)

# ---- User inputs ----

input_file = "HAWKEYE_Database_UniProt_Features.csv"
output_file = "HAWKEYE_Database_UniProt_All_Protein_Domain_Features.csv"

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

exonData = read_csv(input_file, show_col_types = FALSE)

# ---- Helper functions ----

collapse_unique = function(x, sep = "|") {
  x = x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  paste(unique(x), collapse = sep)
}

parse_aa_position = function(aa_pos) {
  aa_pos = str_remove_all(as.character(aa_pos), "\\[|\\]")
  
  if (is.na(aa_pos) || aa_pos == "NA" || aa_pos == "") {
    return(c(NA_real_, NA_real_))
  }
  
  values = str_split(aa_pos, "-", simplify = TRUE)
  
  if (ncol(values) != 2) {
    return(c(NA_real_, NA_real_))
  }
  
  as.numeric(values[1, ])
}

get_uniprot_domain_repeat_features = function(uniprot_id) {
  request_url = paste0(
    "https://rest.uniprot.org/uniprotkb/",
    uniprot_id,
    ".json?fields=ft_domain%2Cft_repeat"
  )
  
  response = GET(request_url, accept("application/json"))
  
  if (status_code(response) != 200) {
    warning(paste("UniProt domain/repeat request failed for", uniprot_id))
    return(NULL)
  }
  
  json_content = fromJSON(rawToChar(response$content), flatten = TRUE)
  
  if (is.null(json_content$features) || length(json_content$features) == 0) {
    return(NULL)
  }
  
  json_content$features %>%
    as_tibble() %>%
    filter(type %in% c("Domain", "Repeat")) %>%
    transmute(
      Uniprot_ID = uniprot_id,
      Feature_Type = type,
      Feature_Description = description,
      Feature_Start = as.numeric(location.start.value),
      Feature_End = as.numeric(location.end.value),
      Feature_Position = paste0(
        "amino acids ",
        location.start.value,
        "-",
        location.end.value,
        " on protein ",
        uniprot_id
      )
    )
}

# ---- Parse exon amino acid positions ----

aa_positions = t(sapply(exonData$AA_Pos, parse_aa_position))

exonData = exonData %>%
  mutate(
    AA_Start = aa_positions[, 1],
    AA_End = aa_positions[, 2],
    Functional_Domains = NA_character_,
    Functional_Domains_Pos = NA_character_,
    Repeat = NA_character_,
    Repeat_Pos = NA_character_
  )

# ---- Download UniProt domain and repeat features ----

uniprot_features = bind_rows(
  lapply(
    unique(na.omit(exonData$Uniprot_ID)),
    get_uniprot_domain_repeat_features
  )
)

# ---- Annotate exon-level overlaps ----

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
  
  exonData$Functional_Domains[i] = collapse_unique(
    overlaps$Feature_Description[overlaps$Feature_Type == "Domain"]
  )
  
  exonData$Functional_Domains_Pos[i] = collapse_unique(
    overlaps$Feature_Position[overlaps$Feature_Type == "Domain"]
  )
  
  exonData$Repeat[i] = collapse_unique(
    overlaps$Feature_Description[overlaps$Feature_Type == "Repeat"]
  )
  
  exonData$Repeat_Pos[i] = collapse_unique(
    overlaps$Feature_Position[overlaps$Feature_Type == "Repeat"]
  )
}

# ---- Remove intermediate columns ----

exonData = exonData %>%
  dplyr::select(-AA_Start, -AA_End)

# ---- Export ----

write_csv(exonData, output_file)
