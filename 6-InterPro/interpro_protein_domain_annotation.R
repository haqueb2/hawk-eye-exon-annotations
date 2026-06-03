library(dplyr)
library(readr)
library(stringr)
library(httr)
library(jsonlite)
library(purrr)
library(tidyr)

# ---- User inputs ----

input_file = "HAWKEYE_Database_UniProt_PTM_Other_Features.csv"
output_file = "HAWKEYE_Database_InterPro_Features.csv"

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

collapse_unique = function(x) {
  x = x[!is.na(x) & x != "" & x != "NA"]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(unique(x), collapse = ";")
}

parse_aa_position = function(aa_pos) {
  aa_pos = str_remove_all(as.character(aa_pos), "\\[|\\]")
  
  if (is.na(aa_pos) || aa_pos == "" || aa_pos == "NA") {
    return(c(NA_real_, NA_real_))
  }
  
  values = str_split(aa_pos, "-", simplify = TRUE)
  
  if (ncol(values) != 2) {
    return(c(NA_real_, NA_real_))
  }
  
  as.numeric(values[1, ])
}

safe_value = function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(NA_character_)
  }
  
  as.character(x[1])
}

fetch_interpro_features = function(uniprot_id) {
  
  request_url = paste0(
    "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/",
    uniprot_id,
    "/?page_size=200"
  )
  
  all_features = list()
  
  repeat {
    
    response = GET(request_url)
    
    if (status_code(response) != 200) {
      warning(paste("InterPro request failed for", uniprot_id))
      return(NULL)
    }
    
    response_text = content(response, "text", encoding = "UTF-8")
    
    if (is.na(response_text) || response_text == "") {
      return(NULL)
    }
    
    json_data = fromJSON(
      response_text,
      simplifyVector = FALSE
    )
    
    if (!is.null(json_data$results) && length(json_data$results) > 0) {
      
      for (result in json_data$results) {
        
        if (is.null(result$metadata)) next
        
        entry_name = result$metadata$name
        entry_type = result$metadata$type
        
        if (is.null(entry_name) || length(entry_name) == 0) {
          entry_name = "InterPro feature"
        }
        
        if (is.null(entry_type) || length(entry_type) == 0) {
          next
        }
        
        if (tolower(entry_type) %in% c("family", "homologous_superfamily")) {
          next
        }
        
        if (is.null(result$proteins) || length(result$proteins) == 0) {
          next
        }
        
        for (protein in result$proteins) {
          
          locations = protein$entry_protein_locations
          
          if (is.null(locations) || length(locations) == 0) {
            next
          }
          
          for (location in locations) {
            
            fragments = location$fragments
            
            if (is.null(fragments) || length(fragments) == 0) {
              next
            }
            
            for (fragment in fragments) {
              
              if (is.null(fragment$start) || is.null(fragment$end)) {
                next
              }
              
              all_features[[length(all_features) + 1]] = tibble(
                Uniprot_ID = uniprot_id,
                InterPro = paste0(entry_name, " (", entry_type, ")"),
                InterPro_Start = as.numeric(fragment$start),
                InterPro_End = as.numeric(fragment$end)
              )
            }
          }
        }
      }
    }
    
    if (is.null(json_data$`next`) || is.na(json_data$`next`) || json_data$`next` == "") {
      break
    }
    
    request_url = json_data$`next`
  }
  
  if (length(all_features) == 0) {
    return(NULL)
  }
  
  bind_rows(all_features)
}

# ---- Parse exon amino acid positions ----

aa_positions = t(sapply(exonData$AA_Pos, parse_aa_position))

exonData = exonData %>%
  mutate(
    AA_Start = as.numeric(aa_positions[, 1]),
    AA_End = as.numeric(aa_positions[, 2]),
    InterPro = NA_character_
  )

# ---- Download InterPro features ----

uniprot_ids = unique(na.omit(exonData$Uniprot_ID))
uniprot_ids = uniprot_ids[uniprot_ids != "" & uniprot_ids != "NA"]

message("Number of UniProt IDs to query: ", length(uniprot_ids))

interpro_features = bind_rows(
  lapply(seq_along(uniprot_ids), function(i) {
    
    if (i %% 100 == 0) {
      message("Queried ", i, " of ", length(uniprot_ids), " UniProt IDs")
    }
    
    fetch_interpro_features(uniprot_ids[i])
  })
)

message("Number of InterPro features fetched: ", nrow(interpro_features))

if (nrow(interpro_features) == 0) {
  stop("No InterPro features were fetched. Check UniProt IDs, API response, or parsing.")
}

# ---- Annotate exon-level InterPro overlaps ----

for (i in seq_len(nrow(exonData))) {
  
  if (
    is.na(exonData$Uniprot_ID[i]) ||
    exonData$Uniprot_ID[i] == "" ||
    exonData$Uniprot_ID[i] == "NA" ||
    is.na(exonData$AA_Start[i]) ||
    is.na(exonData$AA_End[i])
  ) {
    next
  }
  
  overlaps = interpro_features %>%
    filter(
      Uniprot_ID == exonData$Uniprot_ID[i],
      InterPro_Start <= exonData$AA_End[i],
      InterPro_End >= exonData$AA_Start[i]
    )
  
  exonData$InterPro[i] = collapse_unique(overlaps$InterPro)
}

# ---- Remove intermediate columns ----

exonData = exonData %>%
  dplyr::select(-AA_Start, -AA_End)

# ---- Export ----

write_csv(exonData, output_file)