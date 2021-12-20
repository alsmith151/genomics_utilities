import os
import sys
import numpy as np
import pandas as pd
from pybedtools import BedTool
import subprocess
from typing import Union, Iterable


def create_design_matrix(counts, **kwargs):
    """Generates a template dataframe from counts table.
       Keyword arguments are passed to pandas assign"""
    df = counts.columns.to_frame()
    df.columns = ["sample"]
    return df.assign(**kwargs)


def format_feature_counts_dataframe(fn):
    """Reads and formats featureCounts dataframe"""

    df = (
        pd.read_csv(fn, skiprows=1, sep="\t")
        .rename(columns=lambda col: col.lower())
        .set_index("geneid")
    )

    rf_columns = []
    for col in df.columns:
        if "/" in col:
            col = col.split("/")[-1]
        if ".bam" in col:
            col = col.split(".")[0]

        rf_columns.append(col)

    df.columns = rf_columns

    return df


def calculate_tpm(counts, lengths):
    rpk = counts.div(lengths, axis=0)
    scale_factor = rpk.sum() / 1e6
    tpm = rpk / scale_factor
    return tpm


def convert_bed_to_saf(bed, output="out.saf", input_format="narrowPeak", save=False):

    if isinstance(bed, BedTool):
        df = bed.to_dataframe()
    elif isinstance(bed, str):
        df = BedTool(bed).to_dataframe()
    elif isinstance(bed, pd.DataFrame):
        df = bed.copy()

    if "name" in df.columns:
        if "strand" in df.columns:
            saf = df[["name", "chrom", "start", "end", "strand"]]
        else:
            saf = df[["name", "chrom", "start", "end"]].assign(strand=".")

    else:
        df["name"] = df["chrom"].str.cat(df[["start", "end"]].astype(str), sep="|")

        saf = df[["name", "chrom", "start", "end"]].assign(strand=".")

    if input_format == "bed":
        saf["start"] = saf["start"].astype(int) + 1

    if save:
        saf.to_csv(output, sep="\t", header=None, index=None)

    return saf


def run_feature_counts(
    regions: Union[BedTool, pd.DataFrame],
    bam_files: Iterable,
    output: os.PathLike = "featureCounts_out.tsv",
    threads: int = 8,
    paired_end: bool = True,
    allow_multioverlap: bool = True,
    return_dataframe: bool = False,
    **kwargs,
):

    
    print("Converting SAF format")
    saf = convert_bed_to_saf(regions)
    saf.to_csv("regions.saf", header=None, index=None, sep="\t")
    regions = "regions.saf"

    cmd = [
        "featureCounts",
        "-T",
        f"{threads}",
        "-a",
        f"{regions}",
        "-o",
        f"{output}",
        "-F",
        os.path.basename(regions).split(".")[-1].upper(),
    ]

    if paired_end:
        cmd.append("-p")
    if allow_multioverlap:
        cmd.append("-O")

    # Appends arbitrary keywords
    for k, v in kwargs.items():
        cmd.append(k)
        cmd.append(v)

    # Append files
    cmd.append(" ".join(bam_files))

    print(f'Running featureCounts with command: {" ".join(cmd)}')
    subprocess.run(" ".join(cmd), shell=True)

    if os.path.exists("regions.saf"):
        os.remove("regions.saf")

    if return_dataframe:
        return format_feature_counts_dataframe(output)


def identify_expressed_genes(
    counts, min_count=0, min_n_samples=1, replicates=None, min_n_replicates=2
):

    counts_over_thresh = (counts >= min_count).astype(int)

    if not replicates:
        return counts_over_thresh.sum(axis=1) > min_n_samples
    else:
        replicates_over_thresh = counts_over_thresh.groupby(replicates, axis=1).sum()

        genes_over_thresh = (replicates_over_thresh >= min_n_replicates).sum(axis=1)
        return genes_over_thresh > 0
