## ============================================================
##  HAWKEYE FILTERING PIPELINE
## ============================================================

library(readxl)
library(openxlsx)

args = commandArgs(trailingOnly = TRUE)

input_file = ifelse(length(args) >= 1, args[1], "HAWKEYE_Database_Domain_Repeat_Counts_Flagged.csv")
output_file = ifelse(length(args) >= 2, args[2], "HAWK-EYE_Database_Final_Filtered.xlsx")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

## ============================================================
##  0. HELPERS
## ============================================================

to_num = function(x) suppressWarnings(as.numeric(x))

is_na_str = function(x) {
  x = as.character(x)
  is.na(x) | trimws(x) == "" | trimws(x) == "NA"
}

na_false = function(x) { x[is.na(x)] = FALSE; x }

get_max_fraction = function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return(NA_real_)
    val = trimws(val)
    if (val == "" || val == "NA") return(NA_real_)
    parts = trimws(unlist(strsplit(val, ";", fixed = TRUE)))
    nums = suppressWarnings(as.numeric(parts))
    nums = nums[!is.na(nums)]
    if (length(nums) == 0) return(NA_real_)
    max(nums)
  })
}

yn_signature = function(flag_df, order_cols = colnames(flag_df), sep = "; ") {
  flag_df = flag_df[, order_cols, drop = FALSE]
  apply(flag_df, 1, function(r) {
    paste0(names(r), " = ", ifelse(as.logical(r), "Y", "N"), collapse = sep)
  })
}

## ============================================================
##  1. LOAD
## ============================================================

HAWKEYE = read.csv(input_file, stringsAsFactors = FALSE)
HAWKEYE = as.data.frame(HAWKEYE)

required_cols = c(
  "CDS", "Exon_Frame", "LPP_ClinVar_Inframe",
  "First_Coding_Exon", "Last_Coding_Exon", "In_Frame_Exon",
  "Premature_Stop", "OnlyDomains", "OnlyRepeats",
  "MultipleDomains", "MultipleRepeats", "MultipleDomainsANDRepeats",
  "LPP_Missense_ClinVar"
)

missing_cols = setdiff(required_cols, colnames(HAWKEYE))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in input file: ",
    paste(missing_cols, collapse = ", ")
  )
}

HAWKEYE$CDS_raw = as.character(HAWKEYE$CDS)
HAWKEYE$CDS = to_num(HAWKEYE$CDS_raw)
HAWKEYE$Exon_Frame = to_num(HAWKEYE$Exon_Frame)
HAWKEYE$LPP_ClinVar_Inframe = to_num(HAWKEYE$LPP_ClinVar_Inframe)

## ============================================================
##  PASS 1 — INDISPENSABLE / FIRST PASS
## ============================================================

cond_indispensable = na_false(
  is.na(HAWKEYE$CDS) |
    HAWKEYE$First_Coding_Exon == "Y" |
    HAWKEYE$Last_Coding_Exon == "Y" |
    HAWKEYE$Exon_Frame < 0 |
    HAWKEYE$In_Frame_Exon == "N" |
    HAWKEYE$Premature_Stop == "Y" |
    HAWKEYE$OnlyDomains == "Y" |
    HAWKEYE$OnlyRepeats == "Y" |
    HAWKEYE$LPP_ClinVar_Inframe >= 1 |
    (
      HAWKEYE$CDS >= 0.1 &
        (HAWKEYE$MultipleDomains == "Y" |
           HAWKEYE$MultipleRepeats == "Y" |
           HAWKEYE$MultipleDomainsANDRepeats == "Y")
    )
)

idx_indispensable = which(cond_indispensable)

Indispensable_HAWKEYE = HAWKEYE[idx_indispensable, ]
FirstPass_HAWKEYE     = HAWKEYE[-idx_indispensable, ]

## ============================================================
##  PASS 2 — DISPENSABLE
## ============================================================

cols_all_NA = c("Region","Compositional_Bias","Motif","Zinc_Finger","Modifs",
                 "Functional_Domains","Functional_Domains_Pos","Repeat","Repeat_Pos",
                 "Signal_Peptide","Transit_Peptide","Crosslink","Disulfide_Bond",
                 "Glycosylation_Site","Initiator_Methionine","Lipidation_Site",
                 "Peptide","Propeptide","Active_Site","Binding_Site",
                 "DNA_Binding_Region","Topological_Domain","Transmembrane",
                 "Coiled_coil","Other_Site","InterPro")

missing_cols = setdiff(cols_all_NA, colnames(FirstPass_HAWKEYE))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns needed for PASS 2: ",
    paste(missing_cols, collapse = ", ")
  )
}

rows_all_NA = apply(
  sapply(cols_all_NA, function(col) {
    is_na_str(FirstPass_HAWKEYE[[col]])
  }),
  1,
  all
)

lpp_zero = to_num(FirstPass_HAWKEYE$LPP_Missense_ClinVar) == 0
cds_small = to_num(FirstPass_HAWKEYE$CDS) < 0.10

Dispensable_idx = rows_all_NA & !is.na(lpp_zero) & lpp_zero & cds_small

Dispensable_HAWKEYE = FirstPass_HAWKEYE[Dispensable_idx, ]

## ============================================================
##  PASS 3 — INDETERMINATE - UNLIKELY DISPENSABLE
## ============================================================

FP = FirstPass_HAWKEYE

FP$Max_Domain_Fraction = get_max_fraction(FP$Exon_Contribution)
FP$Max_Repeat_Fraction = get_max_fraction(FP$Exon_Contribution_Repeat)

cds_vals = to_num(FP$CDS)
cds_big = !is.na(cds_vals) & cds_vals >= 0.10
cds_small = !is.na(cds_vals) & cds_vals < 0.10

dom_big = na_false(FP$Max_Domain_Fraction > (2/3))
rep_big = na_false(FP$Max_Repeat_Fraction > (2/3))

Indeterminate_Unlikely_Dispensable_idx =
  cds_big | dom_big | rep_big

Indeterminate_Unlikely_Dispensable_HAWKEYE =
  FP[Indeterminate_Unlikely_Dispensable_idx, ]

## ============================================================
##  PASS 4 — INDETERMINATE - POSSIBLY DISPENSABLE
## ============================================================

class_cols = c("Region","Compositional_Bias","Motif","Zinc_Finger","Modifs",
               "Functional_Domains","Functional_Domains_Pos","Repeat",
               "Repeat_Pos","Topological_Domain","Transmembrane","Coiled_coil")

feat_matrix = sapply(class_cols, function(col) {
  !is_na_str(FP[[col]])
})

has_feature = apply(feat_matrix, 1, any)

dom_ok = is.na(FP$Max_Domain_Fraction) | FP$Max_Domain_Fraction <= (2/3)
rep_ok = is.na(FP$Max_Repeat_Fraction) | FP$Max_Repeat_Fraction <= (2/3)

Indeterminate_Possibly_Dispensable_idx =
  cds_small & has_feature & dom_ok & rep_ok & !Dispensable_idx

Indeterminate_Possibly_Dispensable_HAWKEYE =
  FP[Indeterminate_Possibly_Dispensable_idx, ]

## ============================================================
##  PASS 5 — INDETERMINATE - PROBABLY DISPENSABLE
## ============================================================

no_feature = apply(!feat_matrix, 1, all)

Indeterminate_Probably_Dispensable_idx =
  cds_small & no_feature & dom_ok & rep_ok & !Dispensable_idx

Indeterminate_Probably_Dispensable_HAWKEYE =
  FP[Indeterminate_Probably_Dispensable_idx, ]

## ============================================================
##  EXPLANATIONS (UNCHANGED)
## ============================================================

## ------------------------------------------------------------
##  Explanation for Indispensable
## ------------------------------------------------------------

Indispensable_HAWKEYE$Explanation = yn_signature(data.frame(
  "First Coding Exon" = na_false(Indispensable_HAWKEYE$First_Coding_Exon == "Y"),
  "Last Coding Exon"  = na_false(Indispensable_HAWKEYE$Last_Coding_Exon == "Y"),
  "UTR" = na_false(to_num(Indispensable_HAWKEYE$Exon_Frame) < 0),
  "Out-of-Frame" = na_false(Indispensable_HAWKEYE$In_Frame_Exon == "N"),
  "Creates a PTC when Skipped" = na_false(Indispensable_HAWKEYE$Premature_Stop == "Y"),
  "Codes for the Only Functional Domain" = na_false(Indispensable_HAWKEYE$OnlyDomains == "Y"),
  "Codes for the Only Repeat Domain" = na_false(Indispensable_HAWKEYE$OnlyRepeats == "Y"),
  "Contains a LP/P In-Frame Deletion" = na_false(to_num(Indispensable_HAWKEYE$LPP_ClinVar_Inframe) >= 1),
  "≥ 10% of CDS and Codes for Multiple Domains/Repeats" =
    na_false(
      to_num(Indispensable_HAWKEYE$CDS) >= 0.10 &
        (Indispensable_HAWKEYE$MultipleDomains == "Y" |
           Indispensable_HAWKEYE$MultipleRepeats == "Y" |
           Indispensable_HAWKEYE$MultipleDomainsANDRepeats == "Y")
    ),
  check.names = FALSE
))

## ------------------------------------------------------------
##  Explanation for Indeterminate - Unlikely Dispensable
## ------------------------------------------------------------

Indeterminate_Unlikely_Dispensable_HAWKEYE$Explanation = yn_signature(data.frame(
  "≥ 10% of CDS" = na_false(to_num(Indeterminate_Unlikely_Dispensable_HAWKEYE$CDS) >= 0.10),
  "Encodes a Non-Redundant Functional Domain" =
    na_false(get_max_fraction(Indeterminate_Unlikely_Dispensable_HAWKEYE$Exon_Contribution) > (2/3)),
  "Encodes a Non-Redundant Repeat Domain" =
    na_false(get_max_fraction(Indeterminate_Unlikely_Dispensable_HAWKEYE$Exon_Contribution_Repeat) > (2/3)),
  check.names = FALSE
))

## ------------------------------------------------------------
##  Explanation for Indeterminate - Possibly Dispensable
## ------------------------------------------------------------

IP2_flags = data.frame(
  
  "Encodes a Functional Element (Excluding Domains/Repeats)" =
    na_false(apply(
      sapply(c(
        "Region","Compositional_Bias","Motif","Zinc_Finger","Modifs",
        "Topological_Domain","Transmembrane","Coiled_coil"
      ), function(col) {
        !is_na_str(Indeterminate_Possibly_Dispensable_HAWKEYE[[col]])
      }),
      1, any
    )),
  
  "Encodes a Redundant Functional Domain" =
    na_false(
      !is_na_str(Indeterminate_Possibly_Dispensable_HAWKEYE$Functional_Domains) &
        get_max_fraction(Indeterminate_Possibly_Dispensable_HAWKEYE$Exon_Contribution) <= (2/3)
    ),
  
  "Encodes a Redundant Repeat Domain" =
    na_false(
      !is_na_str(Indeterminate_Possibly_Dispensable_HAWKEYE$Repeat) &
        get_max_fraction(Indeterminate_Possibly_Dispensable_HAWKEYE$Exon_Contribution_Repeat) <= (2/3)
    ),
  
  check.names = FALSE
)

Indeterminate_Possibly_Dispensable_HAWKEYE$Explanation =
  yn_signature(
    IP2_flags,
    order_cols = c(
      "Encodes a Functional Element (Excluding Domains/Repeats)",
      "Encodes a Redundant Functional Domain",
      "Encodes a Redundant Repeat Domain"
    )
  )

## ------------------------------------------------------------
##  Explanation for Indeterminate - Probably Dispensable (FIXED)
## ------------------------------------------------------------

IP1_flags = data.frame(
  
  "Encodes PTM Sites and Other General Protein Features" =
    na_false(
      apply(
        sapply(c(
          "Signal_Peptide",
          "Transit_Peptide",
          "Crosslink",
          "Disulfide_Bond",
          "Glycosylation_Site",
          "Initiator_Methionine",
          "Lipidation_Site",
          "Peptide",
          "Propeptide",
          "Active_Site",
          "Binding_Site",
          "DNA_Binding_Region",
          "Other_Site"
        ), function(col) {
          !is_na_str(Indeterminate_Probably_Dispensable_HAWKEYE[[col]])
        }),
        1, any
      )
    ),
  
  "Encodes a Functional Element Predicted by Interpro" =
    na_false(
      !is_na_str(Indeterminate_Probably_Dispensable_HAWKEYE$InterPro)
    ),
  
  "Contains Pathogenic or Likely Pathogenic Missense Variants" =
    na_false(
      to_num(Indeterminate_Probably_Dispensable_HAWKEYE$LPP_Missense_ClinVar) > 0
    ),
  
  check.names = FALSE
)

Indeterminate_Probably_Dispensable_HAWKEYE$Explanation =
  yn_signature(
    IP1_flags,
    order_cols = c(
      "Encodes PTM Sites and Other General Protein Features",
      "Encodes a Functional Element Predicted by Interpro",
      "Contains Pathogenic or Likely Pathogenic Missense Variants"
    )
  )

## ------------------------------------------------------------
##  Explanation for Dispensable
## ------------------------------------------------------------

Dispensable_HAWKEYE$Explanation =
  "Encodes No Functional Elements = Y; Encodes No PTM Sites or Other General Protein Features = Y; Contains No Pathogenic or Likely Pathogenic Missense Variants = Y"

## ============================================================
##  ASSESSMENT
## ============================================================

Indispensable_HAWKEYE$Assessment = "Indispensable"
Dispensable_HAWKEYE$Assessment = "Dispensable"
Indeterminate_Unlikely_Dispensable_HAWKEYE$Assessment = "Indeterminate - Unlikely Dispensable"
Indeterminate_Possibly_Dispensable_HAWKEYE$Assessment = "Indeterminate - Possibly Dispensable"
Indeterminate_Probably_Dispensable_HAWKEYE$Assessment = "Indeterminate - Probably Dispensable"

## ============================================================
##  EXPLANATION SUMMARY
## ============================================================

summarize_explanations = function(df, category_name, explanation_col = "Explanation") {
  if (is.null(df) || nrow(df) == 0 || !(explanation_col %in% colnames(df))) {
    return(data.frame(
      Category = character(0),
      Explanation = character(0),
      Exon_Count = integer(0),
      Percent_of_Category = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  expl = as.character(df[[explanation_col]])
  expl[is.na(expl) | trimws(expl) == "" | trimws(expl) == "NA"] = NA_character_
  
  row_id = seq_len(nrow(df))
  split_list = strsplit(expl, ";", fixed = TRUE)
  
  out = do.call(rbind, lapply(row_id, function(i) {
    if (is.na(expl[i])) return(NULL)
    reasons = trimws(split_list[[i]])
    reasons = reasons[!(is.na(reasons) | reasons == "")]
    if (length(reasons) == 0) return(NULL)
    data.frame(row_id = i, Explanation = reasons, stringsAsFactors = FALSE)
  }))
  
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame(
      Category = character(0),
      Explanation = character(0),
      Exon_Count = integer(0),
      Percent_of_Category = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  out = unique(out)
  
  counts = aggregate(row_id ~ Explanation, data = out,
                      FUN = function(x) length(unique(x)))
  
  colnames(counts)[colnames(counts) == "row_id"] = "Exon_Count"
  
  counts$Category = category_name
  counts$Percent_of_Category = (counts$Exon_Count / nrow(df)) * 100
  
  counts = counts[order(-counts$Exon_Count, counts$Explanation), ]
  counts = counts[, c("Category", "Explanation", "Exon_Count", "Percent_of_Category")]
  
  counts
}

Explanation_Summary = rbind(
  summarize_explanations(Indispensable_HAWKEYE, "Indispensable"),
  summarize_explanations(Dispensable_HAWKEYE, "Dispensable"),
  summarize_explanations(Indeterminate_Unlikely_Dispensable_HAWKEYE, "Indeterminate - Unlikely Dispensable"),
  summarize_explanations(Indeterminate_Possibly_Dispensable_HAWKEYE, "Indeterminate - Possibly Dispensable"),
  summarize_explanations(Indeterminate_Probably_Dispensable_HAWKEYE, "Indeterminate - Probably Dispensable")
)

Explanation_Summary$Percent_of_Category =
  round(Explanation_Summary$Percent_of_Category, 2)

## ============================================================
##  OVERLAP CHECKS BETWEEN ALL CLASSES
## ============================================================

overlap_count = function(df1, df2, name1, name2) {
  n = length(intersect(rownames(df1), rownames(df2)))
  cat(sprintf("%s vs %s: %d\n", name1, name2, n))
}

cat("\n========== OVERLAP SUMMARY ==========\n")

# Indispensable vs others
overlap_count(Indispensable_HAWKEYE, Dispensable_HAWKEYE,
              "Indispensable", "Dispensable")

overlap_count(Indispensable_HAWKEYE, Indeterminate_Unlikely_Dispensable_HAWKEYE,
              "Indispensable", "Indeterminate_Unlikely")

overlap_count(Indispensable_HAWKEYE, Indeterminate_Possibly_Dispensable_HAWKEYE,
              "Indispensable", "Indeterminate_Possibly")

overlap_count(Indispensable_HAWKEYE, Indeterminate_Probably_Dispensable_HAWKEYE,
              "Indispensable", "Indeterminate_Probably")

# Dispensable vs indeterminate groups
overlap_count(Dispensable_HAWKEYE, Indeterminate_Unlikely_Dispensable_HAWKEYE,
              "Dispensable", "Indeterminate_Unlikely")

overlap_count(Dispensable_HAWKEYE, Indeterminate_Possibly_Dispensable_HAWKEYE,
              "Dispensable", "Indeterminate_Possibly")

overlap_count(Dispensable_HAWKEYE, Indeterminate_Probably_Dispensable_HAWKEYE,
              "Dispensable", "Indeterminate_Probably")

# Indeterminate vs indeterminate
overlap_count(Indeterminate_Unlikely_Dispensable_HAWKEYE,
              Indeterminate_Possibly_Dispensable_HAWKEYE,
              "Indeterminate_Unlikely", "Indeterminate_Possibly")

overlap_count(Indeterminate_Unlikely_Dispensable_HAWKEYE,
              Indeterminate_Probably_Dispensable_HAWKEYE,
              "Indeterminate_Unlikely", "Indeterminate_Probably")

overlap_count(Indeterminate_Possibly_Dispensable_HAWKEYE,
              Indeterminate_Probably_Dispensable_HAWKEYE,
              "Indeterminate_Possibly", "Indeterminate_Probably")

## ============================================================
##  EXPORT
## ============================================================

fix_cds_for_export = function(df) {
  df$CDS = ifelse(is.na(df$CDS), "NA", as.character(df$CDS))
  df = df[, !names(df) %in% c("CDS_raw"), drop = FALSE]
  df
}


write.xlsx(
  list(
    Unfiltered_HAWKEYE = fix_cds_for_export(HAWKEYE),
    Dispensable = fix_cds_for_export(Dispensable_HAWKEYE),
    Indeterminate_Probably = fix_cds_for_export(Indeterminate_Probably_Dispensable_HAWKEYE),
    Indeterminate_Possibly = fix_cds_for_export(Indeterminate_Possibly_Dispensable_HAWKEYE),
    Indeterminate_Unlikely = fix_cds_for_export(Indeterminate_Unlikely_Dispensable_HAWKEYE),
    Indispensable = fix_cds_for_export(Indispensable_HAWKEYE),
    Explanation_Summary = Explanation_Summary
  ),
  file = output_file,
  rowNames = FALSE
)

