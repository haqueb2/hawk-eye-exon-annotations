#!/usr/bin/env python3

import logging
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import pyBigWig

# ---- Logging ----

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s"
)

# ---- User inputs (can be overridden via CLI) ----

parser = argparse.ArgumentParser(description="Compute mean/median/max/min phyloP scores per exon")
parser.add_argument("--input", "-i", default="HAWKEYE_Database_Exon_Features.csv", help="Exon features CSV input")
parser.add_argument("--bigwig", "-b", default="phylop/hg38.phyloP100way.bw", help="phyloP bigwig file")
parser.add_argument("--output", "-o", default="HAWKEYE_Database_PhyloP_Exon_Conservation.csv", help="Output CSV file")

args = parser.parse_args()

input_file = Path(args.input)
bigwig_file = Path(args.bigwig)
output_file = Path(args.output)

# ---- Load inputs ----

df = pd.read_csv(input_file)

bw = pyBigWig.open(str(bigwig_file))

# ---- Helper functions ----

def format_chromosome(chromosome):
    chromosome = str(chromosome)
    if chromosome.startswith("chr"):
        return chromosome
    return f"chr{chromosome}"


def get_exon_interval(row):
    """
    Returns 1-based inclusive exon interval.

    For the first coding exon, restrict the interval to the coding portion only.
    """
    if row["First_Coding_Exon"] == "Y" and row["Strand"] == "+":
        start = row["Start_Pos"]
        end = row["Exon_End"]
    elif row["First_Coding_Exon"] == "Y" and row["Strand"] == "-":
        start = row["Exon_Start"]
        end = row["Start_Pos"]
    else:
        start = row["Exon_Start"]
        end = row["Exon_End"]

    return int(start), int(end)


def get_phylop_scores(row):
    chromosome = format_chromosome(row["Chromosome"])

    try:
        start_1based, end_1based = get_exon_interval(row)

        # pyBigWig expects 0-based, half-open intervals.
        start_0based = start_1based - 1
        end_0based = end_1based

        if start_0based < 0 or end_0based <= start_0based:
            logging.warning(
                f"Invalid interval skipped: {chromosome}:{start_1based}-{end_1based}"
            )
            return pd.Series([np.nan, np.nan, np.nan, np.nan])

        mean_score = bw.stats(
            chromosome, start_0based, end_0based, type="mean", exact=True
        )[0]

        max_score = bw.stats(
            chromosome, start_0based, end_0based, type="max", exact=True
        )[0]

        min_score = bw.stats(
            chromosome, start_0based, end_0based, type="min", exact=True
        )[0]

        values = bw.values(chromosome, start_0based, end_0based)
        values = [value for value in values if value is not None and not np.isnan(value)]

        median_score = np.median(values) if len(values) > 0 else np.nan

        return pd.Series([mean_score, median_score, max_score, min_score])

    except Exception as error:
        logging.warning(
            f"Could not retrieve phyloP scores for "
            f"{row.get('Gene', 'NA')} {row.get('Transcript_ID', 'NA')} exon {row.get('Exon', 'NA')}: {error}"
        )
        return pd.Series([np.nan, np.nan, np.nan, np.nan])


# ---- Calculate exon-level phyloP conservation scores ----

df[
    [
        "Mean_Conservation",
        "Median_Conservation",
        "Max_Conservation",
        "Min_Conservation"
    ]
] = df.apply(get_phylop_scores, axis=1)

# ---- Export ----

df.to_csv(output_file, index=False)

bw.close()

logging.info(f"Conservation features written to: {output_file}")