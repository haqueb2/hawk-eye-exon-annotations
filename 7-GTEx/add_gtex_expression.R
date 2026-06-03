library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(httr)
library(jsonlite)
library(biomaRt)
library(purrr)

# ---- User inputs ----

input_file = "HAWKEYE_Database_InterPro_Features.csv"
output_file = "HAWKEYE_Database_GTEx_Expression.csv"

gtex_dataset = "gtex_v10"
gencode_version = "v39"
ensembl_version = 105

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
    Transcript_ID = str_remove(Transcript_ID, "\\.\\d+$"),
    Exon_Start = as.numeric(Exon_Start),
    Exon_End = as.numeric(Exon_End)
  )

# ---- Helper functions ----

clean_column_names = function(x) {
  str_replace_all(x, "[^A-Za-z0-9_]", "_")
}

get_gtex_gencode_id = function(gene) {
  response = GET(
    "https://gtexportal.org/api/v2/reference/gene",
    query = list(
      geneId = gene,
      gencodeVersion = gencode_version
    )
  )
  
  if (status_code(response) != 200) {
    warning(paste("GTEx gene lookup failed for", gene))
    return(NA_character_)
  }
  
  data = fromJSON(content(response, "text", encoding = "UTF-8"))
  
  if (is.null(data$data) || !"gencodeId" %in% names(data$data)) {
    return(NA_character_)
  }
  
  data$data$gencodeId[1]
}

fetch_gtex_paginated = function(endpoint, gencode_id, items_per_page = 1000) {
  page = 0
  
  response = GET(
    endpoint,
    query = list(
      gencodeId = gencode_id,
      datasetId = gtex_dataset,
      page = page,
      itemsPerPage = items_per_page
    )
  )
  
  if (status_code(response) != 200) {
    warning(paste("GTEx request failed for", gencode_id))
    return(NULL)
  }
  
  data = fromJSON(content(response, "text", encoding = "UTF-8"))
  total_pages = data$paging_info$numberOfPages
  
  if (is.null(total_pages) || total_pages == 0) {
    return(NULL)
  }
  
  results = list()
  
  for (page in 0:(total_pages - 1)) {
    response = GET(
      endpoint,
      query = list(
        gencodeId = gencode_id,
        datasetId = gtex_dataset,
        page = page,
        itemsPerPage = items_per_page
      )
    )
    
    if (status_code(response) != 200) {
      warning(paste("GTEx page request failed for", gencode_id, "page", page))
      next
    }
    
    data = fromJSON(content(response, "text", encoding = "UTF-8"))
    
    if (!is.null(data$data) && length(data$data) > 0) {
      results[[length(results) + 1]] = as_tibble(data$data) %>%
        mutate(gencodeId = gencode_id)
    }
  }
  
  bind_rows(results)
}

# ---- Get GENCODE gene IDs from GTEx ----

gene_data = tibble(
  Gene = unique(exonData$Gene)
) %>%
  mutate(
    gencodeId = sapply(Gene, get_gtex_gencode_id)
  )

exonData = exonData %>%
  left_join(gene_data, by = "Gene")

# ---- Get MANE RefSeq to Ensembl transcript mapping ----

ensembl = useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = ensembl_version
)

mane_conversion = getBM(
  attributes = c(
    "ensembl_gene_id_version",
    "refseq_mrna",
    "ensembl_transcript_id_version",
    "transcript_mane_select"
  ),
  filters = "ensembl_gene_id_version",
  values = unique(na.omit(exonData$gencodeId)),
  mart = ensembl
) %>%
  distinct()

exonData = exonData %>%
  left_join(
    mane_conversion,
    by = c(
      "Transcript_ID" = "refseq_mrna",
      "gencodeId" = "ensembl_gene_id_version"
    )
  )

# ---- Fetch median GTEx transcript expression ----

transcript_expression_results = bind_rows(
  lapply(
    unique(na.omit(exonData$gencodeId)),
    function(gencode_id) {
      fetch_gtex_paginated(
        endpoint = "https://gtexportal.org/api/v2/expression/medianTranscriptExpression",
        gencode_id = gencode_id,
        items_per_page = 250
      )
    }
  )
)

# ---- Keep tissues of interest ----

tissues_of_interest = c(
  "Liver",
  "Cells_Cultured_fibroblasts",
  "Cells_EBV-transformed_lymphocytes"
)

brain_tissues = unique(transcript_expression_results$tissueSiteDetailId[
  grepl("Brain", transcript_expression_results$tissueSiteDetailId, ignore.case = TRUE)
])

tissues_of_interest = c(tissues_of_interest, brain_tissues)

transcript_expression_wide = transcript_expression_results %>%
  filter(tissueSiteDetailId %in% tissues_of_interest) %>%
  dplyr::select(transcriptId, tissueSiteDetailId, median) %>%
  pivot_wider(
    names_from = tissueSiteDetailId,
    values_from = median,
    names_prefix = "GTEx_Transcript_TPM_"
  )

names(transcript_expression_wide) = clean_column_names(names(transcript_expression_wide))

exonData = exonData %>%
  left_join(
    transcript_expression_wide,
    by = c("ensembl_transcript_id_version" = "transcriptId")
  )

# ---- Retrieve Ensembl exon IDs and gene-model ranks ----

exon_conversion = getBM(
  attributes = c(
    "ensembl_transcript_id_version",
    "ensembl_gene_id_version",
    "ensembl_exon_id",
    "strand",
    "exon_chrom_start",
    "exon_chrom_end",
    "rank"
  ),
  filters = "ensembl_gene_id_version",
  values = unique(na.omit(exonData$gencodeId)),
  mart = ensembl
) %>%
  mutate(
    exon_chrom_start = as.numeric(exon_chrom_start),
    exon_chrom_end = as.numeric(exon_chrom_end),
    strand_norm = case_when(
      strand %in% c("+", "1", 1) ~ 1L,
      strand %in% c("-", "-1", -1) ~ -1L,
      TRUE ~ NA_integer_
    )
  )

gene_models = exon_conversion %>%
  group_by(ensembl_gene_id_version, ensembl_transcript_id_version) %>%
  summarise(max_rank = max(rank), .groups = "drop") %>%
  group_by(ensembl_gene_id_version) %>%
  slice_max(max_rank, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(
    ensembl_gene_id_version,
    gene_model_transcript = ensembl_transcript_id_version
  )

gene_model_exons = exon_conversion %>%
  inner_join(gene_models, by = "ensembl_gene_id_version") %>%
  filter(ensembl_transcript_id_version == gene_model_transcript) %>%
  group_by(ensembl_gene_id_version, strand_norm) %>%
  arrange(
    if_else(strand_norm == 1L, exon_chrom_start, -exon_chrom_start),
    .by_group = TRUE
  ) %>%
  mutate(gene_model_rank = row_number()) %>%
  ungroup() %>%
  dplyr::select(
    ensembl_gene_id_version,
    ensembl_exon_id,
    exon_chrom_start,
    exon_chrom_end,
    gene_model_rank
  )

exon_conversion = exon_conversion %>%
  rowwise() %>%
  mutate(
    gene_model_rank = {
      gm = gene_model_exons %>%
        filter(ensembl_gene_id_version == .data$ensembl_gene_id_version)
      
      exact_match = gm %>%
        filter(ensembl_exon_id == .data$ensembl_exon_id)
      
      if (nrow(exact_match) > 0) {
        exact_match$gene_model_rank[1]
      } else {
        overlap_match = gm %>%
          filter(
            exon_chrom_start <= .data$exon_chrom_end,
            exon_chrom_end >= .data$exon_chrom_start
          )
        
        if (nrow(overlap_match) > 0) {
          overlap_match$gene_model_rank[1]
        } else {
          NA_integer_
        }
      }
    }
  ) %>%
  ungroup()

exonData = exonData %>%
  left_join(
    exon_conversion,
    by = c(
      "gencodeId" = "ensembl_gene_id_version",
      "ensembl_transcript_id_version" = "ensembl_transcript_id_version",
      "Exon_Start" = "exon_chrom_start",
      "Exon_End" = "exon_chrom_end"
    )
  ) %>%
  mutate(
    exonId = if_else(
      !is.na(gencodeId) & !is.na(gene_model_rank),
      paste0(gencodeId, "_", gene_model_rank),
      NA_character_
    )
  )

# ---- Fetch median GTEx exon expression ----

exon_expression_results = bind_rows(
  lapply(
    unique(na.omit(exonData$gencodeId)),
    function(gencode_id) {
      fetch_gtex_paginated(
        endpoint = "https://gtexportal.org/api/v2/expression/medianExonExpression",
        gencode_id = gencode_id,
        items_per_page = 1000
      )
    }
  )
)

exon_expression_wide = exon_expression_results %>%
  filter(tissueSiteDetailId %in% tissues_of_interest) %>%
  dplyr::select(exonId, tissueSiteDetailId, median) %>%
  pivot_wider(
    names_from = tissueSiteDetailId,
    values_from = median,
    names_prefix = "GTEx_Exon_read_count_"
  )

names(exon_expression_wide) = clean_column_names(names(exon_expression_wide))

exonData = exonData %>%
  left_join(exon_expression_wide, by = "exonId") %>%
  arrange(Gene, Transcript_ID, Exon)

# ---- Export ----

write_csv(exonData, output_file)