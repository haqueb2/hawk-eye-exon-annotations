library(dplyr)
library(readr)
library(stringr)
library(httr)
library(jsonlite)

# ---- User inputs ----

input_file = "HAWKEYE_Database_UniProt_All_Protein_Domain_Features.csv"
output_file = "HAWKEYE_Database_UniProt_PTM_Other_Features.csv"

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
    Exon_Length = as.integer(Exon_Length)
  )

# ---- Helper functions ----

collapse_unique = function(x) {
  x = x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  paste(unique(x), collapse = ";")
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

append_unique_value = function(existing_value, new_value) {
  if (is.na(new_value) || new_value == "") return(existing_value)
  
  if (is.na(existing_value) || existing_value == "" || existing_value == "NA") {
    return(new_value)
  }
  
  collapse_unique(c(str_split(existing_value, ";", simplify = FALSE)[[1]], new_value))
}

# ---- Parse exon amino acid positions ----

aa_positions = t(sapply(exonData$AA_Pos, parse_aa_position))

exonData = exonData %>%
  mutate(
    AA_Start = aa_positions[, 1],
    AA_End = aa_positions[, 2]
  )

# ---- UniProt feature columns ----

feature_to_col = c(
  "Modified residue" = "Modifs",
  "Signal" = "Signal_Peptide",
  "Transit peptide" = "Transit_Peptide",
  "Cross-link" = "Crosslink",
  "Disulfide bond" = "Disulfide_Bond",
  "Glycosylation" = "Glycosylation_Site",
  "Initiator methionine" = "Initiator_Methionine",
  "Lipidation" = "Lipidation_Site",
  "Peptide" = "Peptide",
  "Propeptide" = "Propeptide",
  "Active site" = "Active_Site",
  "Binding site" = "Binding_Site",
  "DNA binding" = "DNA_Binding_Region",
  "Site" = "Other_Site",
  "Topological domain" = "Topological_Domain",
  "Transmembrane" = "Transmembrane",
  "Coiled coil" = "Coiled_coil"
)

for (col in unique(feature_to_col)) {
  exonData[[col]] = NA_character_
}

# ---- Fetch UniProt protein features ----

get_uniprot_features = function(uniprot_id) {
  request_url = paste0(
    "https://rest.uniprot.org/uniprotkb/",
    uniprot_id,
    ".json?fields=",
    paste(
      c(
        "ft_mod_res",
        "ft_signal",
        "ft_transit",
        "ft_crosslnk",
        "ft_disulfid",
        "ft_carbohyd",
        "ft_init_met",
        "ft_lipid",
        "ft_peptide",
        "ft_propep",
        "ft_act_site",
        "ft_binding",
        "ft_dna_bind",
        "ft_site",
        "ft_topo_dom",
        "ft_transmem",
        "ft_coiled"
      ),
      collapse = "%2C"
    )
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
  
  features = as_tibble(json_content$features)
  
  features %>%
    filter(type %in% names(feature_to_col)) %>%
    transmute(
      Uniprot_ID = uniprot_id,
      Feature_Type = type,
      Feature_Column = feature_to_col[type],
      Feature_Description = description,
      Feature_Start = as.numeric(location.start.value),
      Feature_End = as.numeric(if_else(
        is.na(location.end.value),
        location.start.value,
        location.end.value
      )),
      Ligand_Name = if ("ligand.name" %in% names(features)) ligand.name else NA_character_
    )
}

uniprot_features = bind_rows(
  lapply(
    unique(na.omit(exonData$Uniprot_ID)),
    get_uniprot_features
  )
)

# ---- Format feature labels ----

if (nrow(uniprot_features) > 0) {
  uniprot_features = uniprot_features %>%
    mutate(
      Feature_Label = case_when(
        Feature_Type == "Modified residue" ~ paste0(Feature_Description, "(", Feature_Start, ")"),
        Feature_Type == "Binding site" & !is.na(Ligand_Name) & Ligand_Name != "" ~ Ligand_Name,
        Feature_Type == "Coiled coil" ~ "Coiled coil",
        TRUE ~ Feature_Description
      )
    )
}

# ---- Annotate exon-level UniProt feature overlaps ----

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
  
  if (nrow(overlaps) == 0) next
  
  for (feature_col in unique(overlaps$Feature_Column)) {
    labels = overlaps %>%
      filter(Feature_Column == feature_col) %>%
      pull(Feature_Label)
    
    exonData[[feature_col]][i] = collapse_unique(labels)
  }
}

# ---- Clean feature columns ----

for (col in unique(feature_to_col)) {
  exonData[[col]] = sapply(
    exonData[[col]],
    function(x) {
      if (is.na(x) || x == "" || x == "NA") return(NA_character_)
      collapse_unique(str_split(as.character(x), ";", simplify = FALSE)[[1]])
    }
  )
}

# ---- Remove intermediate columns ----

exonData = exonData %>%
  dplyr::select(-AA_Start, -AA_End)

# ---- Export ----

write_csv(exonData, output_file)
