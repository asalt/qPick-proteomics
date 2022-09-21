from dataclasses import dataclass
from collections import defaultdict
import pandas as pd
from utils import get_logger

logger = get_logger(__name__)


peptide_gene_mapper = dict()


@dataclass()
class PeptideObserver:
    _peptidecounts_dict: dict = None
    _genecounts_dict: dict = None
    _peptidegene_mapping: dict = None

    @property
    def peptidecounts_dict(self):
        if self._peptidecounts_dict is None:
            self._peptidecounts_dict = {}
        return self._peptidecounts_dict

    @property
    def genecounts_dict(self):
        if self._genecounts_dict is None:
            self._genecounts_dict = {}
        return self._genecounts_dict

    @property
    def peptidegene_mapping(self):
        if self._peptidegene_mapping is None:
            self._peptidegene_mapping = defaultdict(set)
        return self._peptidegene_mapping

    @staticmethod
    def update_counter_dict(d: dict, counts: pd.Series):

        for k, v in counts.items():
            if k in d:
                d[k] += v
                # logger.info(f"Adding {v} to {k}")
            else:
                d[k] = v
                # logger.info(f"Adding {k} then adding {v}")
        return d

    @staticmethod
    def update_one_to_many_dict(d: dict, one_to_many: pd.Series):

        for k, v in one_to_many.items():
            if k in d:
                d[k] |= set(v)
                # logger.info(f"Adding {v} to {k}")
            else:
                d[k] = set(
                    v,
                )
                # logger.info(f"Adding {k} then adding {v}")
        return d

    def count_from_e2g(self, df):

        logger.info(f"dataframe of length: {len(df)}")
        df["PeptideSet"] = df["PeptidePrint"].str.split("_")
        dfl = df.explode("PeptideSet").drop_duplicates(subset=["PeptideSet"])

        local_pept_gene_mapping = dfl.groupby("PeptideSet").GeneID.unique().to_dict()
        self.update_one_to_many_dict(self.peptidegene_mapping, local_pept_gene_mapping)
        # self.peptidegene_mapping.update(pept_gene_mapping)

        counts = dfl.PeptideSet.value_counts()
        counts.index = counts.index.str.upper()
        # counts.index.name = "PeptidePrint"
        logger.info(f"counts dataframe of length: {len(counts)}")
        self.update_counter_dict(self.peptidecounts_dict, counts)

        self.update_counter_dict(
            self.genecounts_dict, df["GeneID"].drop_duplicates().value_counts()
        )

        return counts

    def make_summary(self):

        peptide_counts_df = pd.Series(self.peptidecounts_dict).to_frame(
            name="PeptideCount"
        )
        gene_counts_df = pd.Series(self.genecounts_dict).to_frame(name="GeneCount")

        peptide_gene_mapping = pd.Series(self.peptidegene_mapping)
        peptide_gene_mapping.index = peptide_gene_mapping.index.str.upper()
        peptide_gene_mapping.index.name = "Peptide"
        peptide_gene_mapping = (
            peptide_gene_mapping.to_frame(name="GeneID").explode("GeneID").reset_index()
        )
        peptide_gene_mapping = peptide_gene_mapping.merge(
            gene_counts_df, left_on="GeneID", right_index=True, how="left"
        )

        out = peptide_counts_df.merge(
            peptide_gene_mapping, left_index=True, right_on="Peptide", how="left"
        )

        out["RelPeptideCount"] = out["PeptideCount"] / out["GeneCount"]
        # out.index.name = "Peptide"
        return out
