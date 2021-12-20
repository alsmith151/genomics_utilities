import rpy2
from rpy2.robjects import pandas2ri, Formula, DataFrame, r
from rpy2.robjects.packages import importr
import subprocess as sub
import pandas as pd
import numpy as np

deseq = importr("DESeq2")
summarised_experiment = importr("SummarizedExperiment")
pandas2ri.activate()

class py_DESeq2:
    """
    DESeq2 object through rpy2
    
    input:
    count_matrix: should be a pandas dataframe with each column as count, and index == geneid
        example:
            sampleA1    sampleA2
        id    
        geneA    5    1
        geneB    4    5
        geneC    1    2
   
   design_matrix: an design matrix in the form of pandas dataframe, see DESeq2 manual, samplenames as rownames
                treatment
    sampleA1        A
    sampleA2        A
    sampleB1        B
    sampleB2        B
    
    design_formula: see DESeq2 manual, example: "~ treatment""
    """

    def __init__(
        self,
        count_matrix: pd.DataFrame,
        design_matrix: pd.DataFrame,
        design_formula: str,
    ):

        self.deseq_result = None
        self.comparison = None
        self.normalized_count_matrix = None
        self.count_matrix = count_matrix
        self.design_matrix = design_matrix
        self.design_formula = Formula(design_formula)
        self.dds = self.get_dds()

    def get_dds(self, **kwargs):
        self.dds = deseq.DESeqDataSetFromMatrix(
            countData=self.count_matrix.values,
            colData=self.design_matrix,
            design=self.design_formula,
        )
        return self.dds

    def run_deseq(self, **kwargs):
        self.dds = deseq.DESeq(self.dds, **kwargs)
        self.normalized_count_matrix = deseq.counts_DESeqDataSet(
            self.dds, normalized=True
        )
        return self

    def get_deseq_results(self, contrast=None, **kwargs):

        """Returns result of running DESeq2 with the specified contrast.
           Contrast in the format ["design_col_name", "Con_A", "Con_B"]
           Alternatively pass keyword arguments to perform other DESeq2 functions"""

        if contrast and not isinstance(contrast, np.ndarray):
            contrast = np.array(contrast)
            self.comparison = deseq.resultsNames(self.dds)
            self.deseq_result = deseq.results(self.dds, contrast, **kwargs)
        else:
            self.deseq_result = deseq.results(self.dds)

        self.deseq_result = py_DESeq2.to_dataframe(self.deseq_result)
        return self.deseq_result

    def get_deseq_lfcshrink_result(self, coef=2, model_type="apeglm"):
        """Returns deseq2 results with lfcShrink applied.
           Note: R is one based so adjust coef appropriately"""
        return py_DESeq2.to_dataframe(
            deseq.lfcShrink(self.dds, coef=coef, type=model_type)
        )

    @property
    def vst_matrix(self, **kwargs):
        matrix = summarised_experiment.assay(deseq.vst(self.dds))
        return pd.DataFrame(
            matrix, columns=self.design_matrix.index, index=self.count_matrix.index
        )

    @property
    def rlog_matrix(self, **kwargs):
        matrix = summarised_experiment.assay(deseq.rlogTransformation(self.dds))
        return pd.DataFrame(
            matrix, columns=self.design_matrix.index, index=self.count_matrix.index
        )

    @staticmethod
    def to_dataframe(x):
        func = r("function(x) data.frame(x)")
        return func(x)


