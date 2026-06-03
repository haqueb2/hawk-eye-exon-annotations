library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# ---- User inputs ----

input_file = "ExonLength_output.csv"
cds_file = "CDS_output.csv"

output_file = "HAWKEYE_Database_Parse_ExonCalculator.csv"

# ---- Override with command-line key=value args (e.g. input_file="path") ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (a in args) {
    if (grepl("=", a)) {
      try(eval(parse(text = a)), silent = TRUE)
    }
  }
}

# ---- Load exon calculator output ----

exonCalc_out = read_csv(input_file, show_col_types = FALSE) %>%
  distinct()

# ---- Load CDS lengths ----

cdsData = read_csv(cds_file, show_col_types = FALSE) %>%
  mutate(
    Transcript_ID = str_remove(transcript_id, "\\..*"),
    CDS_length = as.numeric(CDS_length)
  ) %>%
  dplyr::select(Transcript_ID, CDS_length) %>%
  distinct()

# ---- Clean semicolon-delimited exon fields ----

clean_exon_field = function(x) {
  x = str_replace_all(x, ";;+", ";")
  x = str_remove(x, "^;")
  x = str_remove(x, ";$")
  return(x)
}

exonCalc_out = exonCalc_out %>%
  mutate(
    exonStarts = clean_exon_field(exonStarts),
    exonEnds = clean_exon_field(exonEnds),
    exonLengths = clean_exon_field(str_replace_all(exonLengths, ";0;", ";")),
    exonFrames = clean_exon_field(exonFrames),
    transcript_id = str_remove(transcript_id, "\\..*")
  )

# ---- Helper functions ----

split_exon_values = function(x) {
  values = unlist(strsplit(as.character(x), ";"))
  values = values[values != ""]
  return(values)
}

trim_to_exon_count = function(values, exon_count) {
  if (length(values) == exon_count) {
    return(values)
  }
  
  if (length(values) > exon_count) {
    first_item_length = nchar(values[1])
    last_item_length = nchar(values[length(values)])
    
    if (first_item_length >= last_item_length) {
      return(values[seq_len(exon_count)])
    } else {
      return(tail(values, exon_count))
    }
  }
  
  return(c(values, rep(NA, exon_count - length(values))))
}

make_exon_table = function(row) {
  exon_count = row$exonCount
  
  exon_starts = trim_to_exon_count(split_exon_values(row$exonStarts), exon_count)
  exon_ends = trim_to_exon_count(split_exon_values(row$exonEnds), exon_count)
  exon_lengths = trim_to_exon_count(split_exon_values(row$exonLengths), exon_count)
  exon_frames = trim_to_exon_count(split_exon_values(row$exonFrames), exon_count)
  
  if (row$strand == "-") {
    exon_starts = rev(exon_starts)
    exon_ends = rev(exon_ends)
    exon_lengths = rev(exon_lengths)
    exon_frames = rev(exon_frames)
  }
  
  data.frame(
    Gene = row$gene,
    Transcript_ID = row$transcript_id,
    Chromosome = row$chrom,
    Strand = row$strand,
    Exon = seq_len(exon_count),
    Exon_Start = exon_starts,
    Exon_End = exon_ends,
    Exon_Length = exon_lengths,
    Exon_Frame = exon_frames
  )
}

# ---- Build exon-level dataframe safely ----

exonData = bind_rows(
  lapply(seq_len(nrow(exonCalc_out)), function(i) {
    make_exon_table(exonCalc_out[i, ])
  })
) %>%
  mutate(
    Exon_Start = as.numeric(Exon_Start) + 1,
    Exon_End = as.numeric(Exon_End),
    Exon_Length = as.numeric(Exon_Length),
    Exon_Frame = as.numeric(Exon_Frame)
  )

# ---- Final cleanup ----

exonData = exonData %>%
  arrange(Transcript_ID, Exon)

# ---- Export ----

write_csv(exonData, output_file)