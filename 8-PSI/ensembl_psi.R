library(dplyr)
library(readr)
library(stringr)
library(biomaRt)

# ---- User inputs ----

input_file = "HAWKEYE_Database_GTEx_Expression.csv"
output_file = "HAWKEYE_Database_PSI.csv"

ensembl_version = 115

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
    Transcript_ID = str_remove(Transcript_ID, "\\.\\d+$")
  )

# ---- Standardize Ensembl gene ID column ----

if (!"ensembl_gene_id_version" %in% names(exonData)) {
  
  if ("gencodeId" %in% names(exonData)) {
    
    names(exonData)[names(exonData) == "gencodeId"] = "ensembl_gene_id_version"
    
  } else {
    
    exonData$ensembl_gene_id_version = NA_character_
  }
}

# ---- Connect to Ensembl (try mirrors with fallback) ----

try_use_ensembl <- function(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = NULL) {
  mirrors <- c("www", "useast", "uswest", "asia")
  for (m in mirrors) {
    res <- tryCatch({
      if (m == "www") {
        useEnsembl(biomart = biomart, dataset = dataset, version = version)
      } else {
        useEnsembl(biomart = biomart, dataset = dataset, version = version, mirror = m)
      }
    }, error = function(e) {
      message(sprintf("useEnsembl failed for mirror=%s: %s", m, conditionMessage(e)))
      NULL
    })

    if (!is.null(res)) {
      message(sprintf("Connected to Ensembl using mirror=%s", m))
      return(res)
    }
  }

  stop("Failed to connect to Ensembl via available mirrors (www, useast, uswest, asia).")
}

ensembl <- try_use_ensembl(version = ensembl_version)

# ---- Fill missing Ensembl exon IDs where possible ----

missing_exons = exonData %>%
  filter(is.na(ensembl_exon_id) | ensembl_exon_id == "NA")

if (nrow(missing_exons) > 0) {
  
  transcript_conversion = getBM(
    attributes = c(
      "ensembl_gene_id_version",
      "ensembl_transcript_id_version",
      "transcript_mane_select",
      "refseq_mrna"
    ),
    filters = "refseq_mrna",
    values = unique(na.omit(missing_exons$Transcript_ID)),
    mart = ensembl
  ) %>%
    mutate(
      transcript_mane_select = str_remove(transcript_mane_select, "\\.\\d+$")
    ) %>%
    distinct()
  
  missing_exons_updated = missing_exons %>%
    left_join(
      transcript_conversion,
      by = c("Transcript_ID" = "transcript_mane_select"),
      relationship = "many-to-many",
      suffix = c("", "_new")
    )
  
  exon_conversion = getBM(
    attributes = c(
      "ensembl_gene_id_version",
      "ensembl_transcript_id_version",
      "ensembl_exon_id",
      "exon_chrom_start",
      "exon_chrom_end",
      "rank"
    ),
    filters = "ensembl_transcript_id_version",
    values = unique(na.omit(missing_exons_updated$ensembl_transcript_id_version_new)),
    mart = ensembl
  ) %>%
    distinct()
  
  missing_exons_updated = missing_exons_updated %>%
    left_join(
      exon_conversion,
      by = c(
        "ensembl_gene_id_version_new" = "ensembl_gene_id_version",
        "ensembl_transcript_id_version_new" = "ensembl_transcript_id_version",
        "Exon" = "rank"
      ),
      suffix = c("", "_filled")
    )
  
  exonData = exonData %>%
    left_join(
      missing_exons_updated %>%
        select(
          Gene,
          Transcript_ID,
          Exon,
          ensembl_gene_id_version_new,
          ensembl_transcript_id_version_new,
          ensembl_exon_id_filled
        ),
      by = c("Gene", "Transcript_ID", "Exon")
    ) %>%
    mutate(
      ensembl_gene_id_version = if_else(
        is.na(ensembl_gene_id_version) | ensembl_gene_id_version == "NA",
        ensembl_gene_id_version_new,
        ensembl_gene_id_version
      ),
      ensembl_transcript_id_version = if_else(
        is.na(ensembl_transcript_id_version) | ensembl_transcript_id_version == "NA",
        ensembl_transcript_id_version_new,
        ensembl_transcript_id_version
      ),
      ensembl_exon_id = if_else(
        is.na(ensembl_exon_id) | ensembl_exon_id == "NA",
        ensembl_exon_id_filled,
        ensembl_exon_id
      )
    ) %>%
    select(
      -ensembl_gene_id_version_new,
      -ensembl_transcript_id_version_new,
      -ensembl_exon_id_filled
    )
}

# ---- Retrieve Ensembl exon-transcript relationships ----

exon_transcript_map = getBM(
  attributes = c(
    "ensembl_gene_id_version",
    "ensembl_transcript_id_version",
    "transcript_mane_select",
    "ensembl_exon_id",
    "external_gene_name"
  ),
  filters = "external_gene_name",
  values = unique(exonData$Gene),
  mart = ensembl
) %>%
  distinct()

# ---- Calculate exon inclusion across transcripts ----

psi_table = exon_transcript_map %>%
  group_by(ensembl_gene_id_version, ensembl_exon_id) %>%
  summarise(
    n_transcripts_including = n_distinct(ensembl_transcript_id_version),
    .groups = "drop"
  ) %>%
  left_join(
    exon_transcript_map %>%
      group_by(ensembl_gene_id_version) %>%
      summarise(
        total_transcripts = n_distinct(ensembl_transcript_id_version),
        .groups = "drop"
      ),
    by = "ensembl_gene_id_version"
  ) %>%
  mutate(
    percent_spliced_in = 100 * n_transcripts_including / total_transcripts
  )

# ---- Annotate HAWK-EYE exons with PSI ----

exonData = exonData %>%
  left_join(
    psi_table,
    by = c("ensembl_gene_id_version", "ensembl_exon_id")
  ) %>%
  mutate(
    exon_class = case_when(
      percent_spliced_in == 100 ~ "Constitutive",
      percent_spliced_in < 100 ~ "Alternative",
      is.na(percent_spliced_in) ~ "Unknown",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(Gene, Transcript_ID, Exon)

# ---- Export ----

write_csv(exonData, output_file)