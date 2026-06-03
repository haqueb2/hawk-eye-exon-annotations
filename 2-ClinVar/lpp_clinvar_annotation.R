library(data.table)
library(readr)
library(stringr)

# ---- User inputs ----

exon_file = "HAWKEYE_Database_Parse_ExonCalculator.csv"
clinvar_summary_file = "clinvar_db/variant_summary.txt"
clinvar_inframe_file = "clinvar_db/CLINVAR_20251106_INS_DEL_DUP_INV.csv"

output_all_lpp = "clinvar_LPP_variants.csv"
output_lpp_exonic = "clinvar_LPP_variants_in_HAWKEYE_exons.csv"
output_exon_counts = "HAWKEYE_Database_ClinVar_Counts.csv"

# ---- Override with command-line key=value args (e.g. exon_file="path") ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (a in args) {
    if (grepl("=", a)) {
      try(eval(parse(text = a)), silent = TRUE)
    }
  }
}

# ---- Helper functions ----

add_chr_prefix = function(x) {
  x = as.character(x)
  fifelse(str_detect(x, "^chr"), x, paste0("chr", x))
}

collapse_variant_ids = function(x) {
  x = x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  paste(unique(x), collapse = ";")
}

count_overlaps_by_exon = function(exon_dt, variant_dt, count_col) {
  
  if (nrow(variant_dt) == 0) {
    exon_dt[, (count_col) := 0L]
    return(exon_dt)
  }
  
  vars = copy(variant_dt)
  vars[, `:=`(
    var_start = as.integer(Start),
    var_end = as.integer(Stop)
  )]
  
  setkey(exon_dt, Gene, Chromosome, Exon_Start, Exon_End)
  setkey(vars, GeneSymbol, Chromosome, var_start, var_end)
  
  overlaps = foverlaps(
    vars,
    exon_dt,
    by.x = c("GeneSymbol", "Chromosome", "var_start", "var_end"),
    by.y = c("Gene", "Chromosome", "Exon_Start", "Exon_End"),
    type = "any",
    nomatch = 0
  )
  
  counts = overlaps[, .N, by = Original_Exon_Row]
  setnames(counts, "N", count_col)
  
  exon_dt[counts, on = "Original_Exon_Row", (count_col) := get(paste0("i.", count_col))]
  exon_dt[is.na(get(count_col)), (count_col) := 0L]
  
  return(exon_dt)
}

# ---- Load exon data ----

exonData = fread(exon_file)

exonData[, `:=`(
  Exon_Start = as.integer(Exon_Start),
  Exon_End = as.integer(Exon_End),
  Exon_Length = as.integer(Exon_Length),
  Chromosome = add_chr_prefix(Chromosome),
  Original_Exon_Row = .I
)]

# ---- Load and filter ClinVar summary ----

clinvar_summary = fread(clinvar_summary_file)

genes_to_keep = unique(exonData$Gene)

clinLPP = clinvar_summary[
  str_detect(ClinicalSignificance, regex("pathogenic", ignore_case = TRUE)) &
    !str_detect(ClinicalSignificance, regex("conflicting", ignore_case = TRUE)) &
    str_detect(OriginSimple, regex("germline", ignore_case = TRUE)) &
    !str_detect(OriginSimple, regex("somatic", ignore_case = TRUE)) &
    Assembly == "GRCh38" &
    GeneSymbol %in% genes_to_keep
]

clinLPP[, Chromosome := add_chr_prefix(Chromosome)]
clinLPP[, `:=`(
  Start = as.integer(Start),
  Stop = as.integer(Stop)
)]

# ---- Identify ClinVar LPP variants overlapping HAWK-EYE exons ----

exon_intervals = exonData[, .(
  Gene,
  Chromosome,
  Exon_Start,
  Exon_End,
  Original_Exon_Row
)]

clinLPP_intervals = copy(clinLPP)
clinLPP_intervals[, `:=`(
  var_start = Start,
  var_end = Stop,
  Original_ClinVar_Row = .I
)]

setkey(exon_intervals, Gene, Chromosome, Exon_Start, Exon_End)
setkey(clinLPP_intervals, GeneSymbol, Chromosome, var_start, var_end)

overlap_hits = foverlaps(
  clinLPP_intervals,
  exon_intervals,
  by.x = c("GeneSymbol", "Chromosome", "var_start", "var_end"),
  by.y = c("Gene", "Chromosome", "Exon_Start", "Exon_End"),
  type = "any",
  nomatch = 0
)

overlap_row_ids = unique(overlap_hits$Original_ClinVar_Row)

overlapLPP = unique(clinLPP[overlap_row_ids])

# ---- Classify variant consequence categories ----

missenseLPP = overlapLPP[
  str_detect(Type, regex("single nucleotide variant", ignore_case = TRUE)) &
    str_detect(Name, "\\(p\\.[A-Za-z]+\\d+[A-Za-z]+\\)") &
    !str_detect(Name, regex("Ter", ignore_case = TRUE))
]

nonsenseLPP = overlapLPP[
  str_detect(Type, regex("single nucleotide variant", ignore_case = TRUE)) &
    str_detect(Name, "\\(p\\.[A-Za-z]+\\d+Ter\\)")
]

snvLPP = overlapLPP[
  str_detect(Type, regex("single nucleotide variant", ignore_case = TRUE))
]

# ---- Add exon-level ClinVar SNV counts ----

exonData = count_overlaps_by_exon(exonData, overlapLPP, "LPP_ClinVar")
exonData = count_overlaps_by_exon(exonData, missenseLPP, "LPP_Missense_ClinVar")
exonData = count_overlaps_by_exon(exonData, nonsenseLPP, "LPP_Nonsense_ClinVar")
exonData = count_overlaps_by_exon(exonData, snvLPP, "LPP_SNV_ClinVar")

# ---- Load and filter ClinVar in-frame deletions ----

clinvar_inframe = fread(clinvar_inframe_file)

clinvar_inframe = clinvar_inframe[
  FRAME == "inframe" &
    Type == "Deletion" &
    str_detect(ClinicalSignificance, regex("pathogenic", ignore_case = TRUE)) &
    !str_detect(ClinicalSignificance, regex("conflict", ignore_case = TRUE))
]

clinvar_inframe[, `:=`(
  Chromosome = add_chr_prefix(Chromosome),
  POS = as.integer(POS),
  VariationID = as.character(VariationID),
  pos_start = as.integer(POS),
  pos_end = as.integer(POS)
)]

# ---- Add exon-level ClinVar in-frame deletion counts and IDs ----

if (nrow(clinvar_inframe) > 0) {
  
  setkey(exon_intervals, Gene, Chromosome, Exon_Start, Exon_End)
  setkey(clinvar_inframe, GeneSymbol, Chromosome, pos_start, pos_end)
  
  inframe_hits = foverlaps(
    clinvar_inframe,
    exon_intervals,
    by.x = c("GeneSymbol", "Chromosome", "pos_start", "pos_end"),
    by.y = c("Gene", "Chromosome", "Exon_Start", "Exon_End"),
    type = "within",
    nomatch = 0
  )
  
  inframe_summary = inframe_hits[, .(
    LPP_ClinVar_Inframe = .N,
    LPP_ClinVar_Inframe_VarID = collapse_variant_ids(VariationID)
  ), by = Original_Exon_Row]
  
  exonData[inframe_summary, on = "Original_Exon_Row", `:=`(
    LPP_ClinVar_Inframe = i.LPP_ClinVar_Inframe,
    LPP_ClinVar_Inframe_VarID = i.LPP_ClinVar_Inframe_VarID
  )]
}

exonData[is.na(LPP_ClinVar_Inframe), LPP_ClinVar_Inframe := 0L]
exonData[is.na(LPP_ClinVar_Inframe_VarID), LPP_ClinVar_Inframe_VarID := NA_character_]

# ---- Remove helper row index before export ----

exonData[, Original_Exon_Row := NULL]

# ---- Export outputs ----

write_csv(as.data.frame(clinLPP), output_all_lpp)
write_csv(as.data.frame(overlapLPP), output_lpp_exonic)
write_csv(as.data.frame(exonData), output_exon_counts)