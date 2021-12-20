import itertools
import functools
import glob
#import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool
import os
import subprocess as sub
import wget


def calculate_rpkm(counts, lengths=None):
    a = counts.mul((1 * 10 ** 9))
    b = lengths.mul(counts.sum())
    return a / b


def categorise_differential_peaks(
    row,
    fc_col="log2FoldChange",
    pval_col="padj",
    up_label="SEM",
    down_label="RS411",
    lfc_thresh=0,
    sig_thresh=0.05,
):
    if row[pval_col] < sig_thresh:
        if row[fc_col] > lfc_thresh:
            return up_label
        elif row[fc_col] < -lfc_thresh:
            return down_label
        else:
            return "Unchanged"
    else:
        return "Unchanged"


def generate_ma_plot(
    df,
    kind="deseq2",
    color_col="",
    colors=["red", "blue", "grey"],
    label=False,
    min_x=0,
    max_x=10,
    font_size=10,
):

    import altair as alt
    parameters = {
        "deseq2": {"x": "log2BaseMean", "y": "log2FoldChange"},
        "edger": {"x": "logCPM", "y": "FDR"},
    }
    params_selected = parameters.get(kind)

    if color_col:
        ma = (
            alt.Chart(df, height=500, width=500)
            .mark_circle(size=20)
            .encode(
                x=alt.X(
                    params_selected["x"],
                    title="log2(Mean readcount)",
                    scale=alt.Scale(domain=(min_x, max_x), clamp=True),
                ),
                y=alt.Y(params_selected["y"], title="log2(Fold Change)"),
                color=alt.Color(color_col, scale=alt.Scale(range=colors)),
            )
        )

        if label:
            df_text = (
                df.groupby(f"{color_col}")
                .size()
                .reset_index()
                .rename(columns={0: "count"})
            )

            df_text["x"] = np.repeat(max_x - 1, len(df_text))
            df_text["y"] = df.groupby(color_col)[params_selected["y"]].mean().to_list()

            ma += (
                alt.Chart(df_text)
                .mark_text(fontSize=font_size)
                .encode(x="x", y="y", text="count")
            )

    else:
        ma = (
            alt.Chart(df, height=500, width=500)
            .mark_circle(size=20)
            .encode(
                x=alt.X(
                    params_selected["x"],
                    title="log2(Mean readcount)",
                    scale=alt.Scale(domain=(min_x, max_x), clamp=True),
                ),
                y=alt.Y(params_selected["y"], title="log2(Fold Change)"),
            )
        )

    ln = (
        alt.Chart(df)
        .mark_rule(strokeDash=[3, 1])
        .encode(y="a:Q")
        .transform_calculate(a="0")
    )

    return ma + ln


def convert_bed_to_saf(bed, output="out.saf", input_format="narrowPeak"):

    if isinstance(bed, BedTool):
        df = bed.to_dataframe()
    elif isinstance(bed, str):
        df = BedTool(bed).to_dataframe()

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

    saf.to_csv(output, sep="\t", header=None, index=None)

    return saf


def run_feature_counts(
    regions,
    bam_files,
    output="featureCounts_out.tsv",
    threads=8,
    paired_end=True,
    allow_multioverlap=True,
):

    cmd = [
        "featureCounts",
        "-T",
        f"{threads}",
        "-a",
        f"{regions}",
        "-o",
        f"{output}",
        " ".join(bam_files),
        "-F",
        os.path.basename(regions).split(".")[-1].upper(),
    ]

    if paired_end:
        cmd.append("-p")
    if allow_multioverlap:
        cmd.append("-O")

    sub.run(" ".join(cmd), shell=True)


def get_overlap_sets(a, b):
    return [(a - b).count(), (b - a).count(), (a + b).count()]


# def calculate_tpm(counts, lengths):
#     rpk = (counts / lengths * 1e3)
#     per_million_scaling_factor = (rpk.sum() / 1e6)
#     return rpk / per_million_scaling_factor


def standardise_peak_width(peaks, output="bt", length=1000):

    """Uses the midpoint of each interval to standardise each interval to a
       predefined length. Output is returned as a dataframe or a Bedtool object"""

    length = int(length)
    print(f"Standardising peaks to {length/1000}kb")

    if isinstance(peaks, BedTool):
        peaks = peaks.to_dataframe()
    elif isinstance(peaks, str):
        peaks = BedTool(peaks).to_dataframe()

    # Generate a new dataframe and determine start and end coordinates based on the midpoints
    df_coords = pd.DataFrame()
    df_coords["chrom"] = peaks["chrom"]
    df_coords["midpoint"] = (
        peaks["start"] + (peaks["end"] - peaks["start"]) / 2
    ).astype(int)
    df_coords["start"] = (df_coords["midpoint"] - length / 2).astype(int)
    df_coords["end"] = (df_coords["midpoint"] + length / 2).astype(int)
    df_coords["name"] = peaks["name"]

    # Drop unwanted columns
    df_coords = df_coords.drop(columns=["midpoint"])

    if output == "df":
        return df_coords
    else:
        return BedTool.from_dataframe(df_coords)


def get_bins(start=1, stop=10, step=1):
    bins = list(
        itertools.chain.from_iterable(
            [(-np.inf, 0), range(start, stop, step), (np.inf,)]
        )
    )
    labels = pd.Series([f"{bins[i]}-{bins[i+1]}" for i in range(len(bins) - 1)])

    labels = (
        labels.str.replace("-inf-", "<")
        .str.replace("0-1", "0")
        .str.replace("-inf", "+")
    )

    return (bins, labels)


def get_groupby_counts(df, cat1, cat2):
    return (
        df.groupby([cat1, cat2])
        .size()
        .to_frame()
        .rename(columns={0: "count"})
        .reset_index()
        .pivot(index=cat1, columns=cat2, values="count")
        .fillna(0)
    )


def bin_series(ser, **kwargs):
    """Uses the get_bins function to bin a pandas series.
       provide: start, stop and/or step to specify bins."""
    bins, labels = get_bins(**kwargs)
    return pd.cut(ser, bins=bins, labels=labels, right=False)


def set_up_dirs(base_path: os.PathLike, nb_name: str) -> dict:

    dirs = {}
    dirs["processed"] = f"processed_data/{nb_name}"
    dirs["results"] = f"results/{nb_name}"
    dirs["figures"] = f"figures/{nb_name}"

    try:
        for d in dirs.values():
            os.makedirs(d, exist_ok=True)
    except Exception as e:
        print(e)

    return dirs


def download(url, out=""):
    if os.path.isfile(out):
        os.remove(out)

    fn = wget.download(url, out=out)
    return fn


def prepare_peaks(bed, merge_distance=2000):
    bt = BedTool(bed).sort().merge(d=merge_distance)
    df = bt.to_dataframe()
    if "name" in df.columns.to_list():
        return bt
    else:
        return df.assign(
            name=os.path.basename(bed).split(".")[0] + "_" + df.index.astype(str)
        ).pipe(BedTool.from_dataframe)


def split_gtf_attributes(attributes):
    """Generates a dataframe from field=value pairs of attributes
       """
    kv_pairs_re = re.compile(r"(.+?)=(.+?)(;|$)")

    res = Parallel(n_jobs=8)(
        delayed(get_attribute_dict)(kv_pairs_re, att) for att in attributes
    )
    return pd.DataFrame(res)


def get_attribute_dict(exp, att):
    return {e[0]: e[1] for e in exp.findall(att)}


def gtf_to_dataframe(gtf: str):
    """Converts gtf to dataframe using pybedtools.
       
       Remember that bed files are 0 based
    """
    return (
        BedTool(gtf)
        .to_dataframe()
        .dropna()
        .assign(
            start=lambda df: df["start"].astype(int),
            end=lambda df: df["end"].astype(int),
        )
    )


class RefSeqGtf:
    def __init__(
        self, gtf_url=None, conv_url=None, download_location="data/misc/", genome="hg19"
    ):

        self.gtf_url, self.conv_url = gtf_url, conv_url
        default_urls = self.get_default_urls()

        if not gtf_url:
            self.gtf_url = default_urls[genome]["refseq_gtf"]
            self.conv_url = default_urls[genome]["conversion_table"]

        self.gtf_path = os.path.join(download_location, os.path.basename(self.gtf_url))
        self.conv_path = os.path.join(
            download_location, os.path.basename(self.conv_url)
        )

    @property
    def gtf_df_converted(self):
        df_conv = self._load_conversion_file()
        df_gtf = self.gtf_df
        return df_gtf.merge(
            df_conv[["RefSeq-Accn", "UCSC-style-name"]],
            left_on="seqname",
            right_on="RefSeq-Accn",
        ).rename(columns={"UCSC-style-name": "chrom"})

    @property
    def gtf_file(self):
        return RefSeqGtf._download(self.gtf_url, self.gtf_path)

    @property
    def chrom_conv_file(self):
        return RefSeqGtf._download(self.conv_url, self.conv_path)

    @property
    def gtf_df(self):
        return RefSeqGtf.gtf_to_dataframe(self.gtf_file)

    @property
    def chrom_conv_df(self):
        return self._load_conversion_file()

    @staticmethod
    def _download(url, out):
        return wget.download(url, out)

    @staticmethod
    def gtf_to_dataframe(gtf: str):
        return (
            BedTool(gtf)
            .to_dataframe()
            .dropna()
            .assign(
                start=lambda df: df["start"].astype(int),
                end=lambda df: df["end"].astype(int),
            )
        )

    def get_default_urls(self):
        return {
            "hg19": {
                "refseq_gtf": "".join(
                    [
                        "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/",
                        "annotation/GRCh37_latest/refseq_identifiers/",
                        "GRCh37_latest_genomic.gff.gz",
                    ]
                ),
                "conversion_table": "".join(
                    [
                        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/",
                        "000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13",
                        "_GRCh37_assembly_report.txt",
                    ]
                ),
            }
        }

    def _load_conversion_file(self):
        return pd.read_csv(
            self.chrom_conv_file,
            comment="#",
            sep="\t",
            header=None,
            names=[
                "Sequence-Name",
                "Sequence-Role",
                "Assigned-Molecule",
                "Assigned-Molecule-Location/Type",
                "GenBank-Accn",
                "Relationship",
                "RefSeq-Accn",
                "Assembly-Unit",
                "Sequence-Length",
                "UCSC-style-name",
            ],
        )


def find_tss_from_gtf(df_gtf, features=["transcript", "gene"]):
    return df_gtf.query("feature in @features").assign(
        tss_start=lambda df: np.where(df["strand"] == "+", df["start"], df["end"]),
        tss_end=lambda df: df["tss_start"] + 1,
    )


def train_test_validation_split_genome(df, validation_chroms, test_chroms):
    df.rename(columns={"chr": "chrom"})
    return {
        "train": df.loc[
            lambda df: ~(df["chrom"].isin([*validation_chroms, *test_chroms]))
        ],
        "validation": df.loc[lambda df: df["chrom"].isin(validation_chroms)],
        "test": df.loc[lambda df: df["chrom"].isin(test_chroms)],
    }


def merge_bams(bams: list, output: str, n_threads=4) -> str:
    cmd = f'samtools merge {output} {" ".join(bams)} -@ {n_threads}'
    sub.run(cmd, shell=True)
    return output


def macs2_call_peaks(
    bam: str, outdir=".", name="called_peaks", bam_format="BAMPE", genome="hs"
):
    cmd = f"macs2 callpeak -t {bam} --outdir {outdir} -n {name} -f {bam_format} -g {genome}"
    sub.run(cmd, shell=True)
    return os.path.join(outdir, f"{name}_peaks.narrowPeak")
