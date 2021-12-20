import os
import sys
import subprocess as sub
from pybedtools import BedTool
from typing import Union
import numpy as np
import pandas as pd


def compute_matrix_deeptools(
    bed: str,
    bigwigs: list,
    bin_size: int = 100,
    bp_before: int = 1000,
    bp_after: int = 1000,
    reference: str = "center",
    n_processes: int = 12,
    output: str = "matrix.gz",
    **compute_matrix_kwargs,
):

    cmd = [
        "computeMatrix",
        "reference-point",
        "--referencePoint",
        reference,
        "-R",
        bed,
        "-a",
        str(bp_after),
        "-b",
        str(bp_before),
        "-bs",
        str(bin_size),
        "-S",
        " ".join(bigwigs),
        "-o",
        output,
        "-p",
        str(n_processes),
    ]

    for k, v in compute_matrix_kwargs.items():
        cmd.append(str(k))
        cmd.append(str(v))

    cmd = " ".join(cmd)
    print(cmd)

    sub.run(cmd, shell=True)
    return output


def plot_heatmap_deeptools(matrix: str, output="heatmap.png", **plot_heatmap_kwargs):

    cmd = [
        "plotHeatmap",
        "-m",
        matrix,
        "--outFileName",
        output,
    ]

    for k, v in plot_heatmap_kwargs.items():
        cmd.append(str(k))
        cmd.append(str(v))

    sub.run(" ".join(cmd), shell=True)
    return output


def multibigwig_summary_deeptools(
    bigwigs: list,
    bed: Union[str, BedTool],
    output: str = "out.npz",
    n_processors=8,
    use_env=None,
):

    if isinstance(bed, BedTool):
        bed = bed.fn

    cmd = [
        "multiBigwigSummary",
        "BED-file",
        "-b",
        " ".join(bigwigs),
        "-o",
        output,
        "--BED",
        bed,
        "-p",
        str(n_processors),
    ]

    if use_env:
        cmd = run_with_conda_env(" ".join(cmd), use_env)
    else:
        cmd = " ".join(cmd)

    sub.run(cmd, shell=True)
    return output


def multibam_summary_deeptools(
    bams: list,
    bed: Union[str, BedTool],
    output: str = "out.npz",
    n_processors=8,
    scaling_factors: str = "",
    raw_counts: str = "",
):

    if isinstance(bed, BedTool):
        bed = bed.fn

    cmd = [
        "multiBamSummary",
        "BED-file",
        "-b",
        " ".join(bams),
        "-o",
        output,
        "--BED",
        bed,
        f"--scalingFactors {scaling_factors}" if scaling_factors else "",
        "-e",
        f"--outRawCounts {raw_counts}" if raw_counts else "",
        "-p",
        str(n_processors),
    ]

    sub.run(" ".join(cmd), shell=True)
    return output


def bigwig_compare_deeptools(bigwigs: list, output: str = "out.npz", n_processes=8):

    cmd = [
        "bigwigCompare",
        "BED-file",
        "-b1",
        bigwigs[0],
        "-b2",
        bigwigs[1],
        "-o",
        output,
        "-p",
        str(n_processors),
    ]

    sub.run(" ".join(cmd), shell=True)
    return output


def extract_dataframe_from_deeptools_matrix(mat_fn, index=None):

    counts = np.load(mat_fn)
    df = pd.DataFrame(data=counts["matrix"], columns=counts["labels"], index=index)
    return df


def run_with_conda_env(cmd, env):

    activate = f"bash -c 'source /t1-data/user/asmith/Software/miniconda3/bin/activate /t1-data/user/asmith/Software/miniconda3/envs/{env}"
    deactivate = "source deactivate'"
    return " && ".join([activate, cmd, deactivate])

