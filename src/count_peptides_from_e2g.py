from dataclasses import dataclass
import ipdb
import re

# count
from pathlib import Path
import pandas as pd
import typer
from collections import Counter
from io import StringIO
from datetime import datetime

date_ = datetime.now().strftime("%Y%m%d%H%M%S%p")
import logging

logging.basicConfig(level=logging.INFO)
# app = typer.Typer()
run_app = typer.Typer(chain=True)
# app.add_typer(run_app, name="run")

from containers import PeptideObserver

def _get_logger():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)
    try:
        fh = logging.FileHandler(f"{__file__}.log")
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    except PermissionError:
        pass

    # fh.setLevel(logger.DEBUG)
    # create console handler with a higher log level
    return logger


logger = _get_logger()

# shelf
from datetime import datetime

date = datetime.now().strftime("%Y%m%d%H%M%S")


import pytest
from dataclasses import dataclass

from containers import PeptideObserver

#     geneid: int
#     symbol: str
#     peptides = None
#     def add_peptide(self, peptide):
#         if self.peptides is None:
#             self.peptides = []
#         self.peptides.append(peptide)
#
# class Peptide:
#     pass

def keep_e2g_file(p: Path, require_regex_match=True):
    """
    keeps e2g qual only?
    """
    logger.debug(f"maybe keep file {p.name}")
    if "QUANT" in p.name:
        return False
    if "TMT12" in p.name:
        return False
    if "TMT13" in p.name:
        return False
    if "iSPEC_import" in p.name or "ispec_import" in p.name:
        return False
    _match = re.match("\\d{5,}_\\d{1,}_.*", p.name)
    if (
        _match is None and require_regex_match == True
    ):  # TODO: skip this if running test
        return False
    if "QUAL" not in p.name:
        return False
    logger.debug(f"kept")
    logger.info(f"kept {p}")
    return True


def count_from_e2g(df):
    res = df["PeptidePrint"].str.split("_").explode()
    import ipdb; ipdb.set_trace()
    counts = res.value_counts()
    counts.index = counts.index.str.upper()
    counts.index.name = "PeptidePrint"
    logger.info(f"counts dataframe of length: {len(counts)}")
    return counts


def update(d: dict, counts: pd.Series):
    for k, v in counts.items():
        if k in d:
            d[k] += v
            # logger.info(f"Adding {v} to {k}")
        else:
            d[k] = v
            # logger.info(f"Adding {k} then adding {v}")
    return d


def read(path: Path):
    logger.info(f"Reading {path}")
    df = pd.read_csv(path, sep="\t")
    if "PeptidePrint" not in df.columns:
        logger.warn(f"PeptidePrint not in {path}")
        return None
    return df


def read_and_count(p: Path):
    df = pd.read_csv(p, sep="\t")
    return count(df)


def get_count_dict():
    d = dict()
    return d


def shape_counts(counts_dict: dict):
    out = pd.Series(counts_dict).to_frame(name="Count")
    out["RelCount"] = out["Count"] / out["Count"].sum()
    out.index.name = "Peptide"
    return out


def clean(df: pd.DataFrame):
    front = [
        "geneid",
        "symbol",
        "peptide",
        "miscuts_trypsin",
        "unique_human",
        "Count",
        "RelCount",
        "before",
        "after",
    ]
    assert all(x in df for x in front)
    if "Peptide" in df.columns:
        df = df.drop(columns=["Peptide"], axis=1)  # the column is called peptide
    if "protein_index" in df.columns:
        df = df.drop(columns=["protein_index"], axis=1)  # the column is called peptide
    if "gi" in df.columns:
        df = df.drop(columns=["gi"], axis=1)  # the column is called peptide
    df = df.drop_duplicates(subset=front, inplace=False)
    df = df.reindex(columns=front + [x for x in df.columns if x not in front])
    df = df.sort_values(by=["geneid", "RelCount"], ascending=False)
    return df


# def test_main()
@run_app.command()
def merge(
    remove_missing: bool = typer.Option(
        True, help="Remove peptides that are not seen empirically"
    ),
    peptidome: Path = typer.Argument(
        ..., help="Path to the peptidome file", exists=True
    ),
    counts: Path = typer.Argument(..., help="Path to the counts file", exists=True),
):
    pept_df = pd.read_table(peptidome)
    counts_df = pd.read_table(counts)

    assert all(x in pept_df for x in ("peptide", "before", "after"))
    import ipdb

    ipdb.set_trace()
    df = pept_df.merge(counts_df, right_on="Peptide", left_on="peptide", how="left")
    if remove_missing:
        df = df.dropna(subset=["Count"])

    df = clean(df)

    outname = f"qpick_peptidome_counts_filtered_{remove_missing}_{date_}.tsv"
    logger.info(f"Writing {outname}")
    df.to_csv(outname, sep="\t", index=False)


@run_app.command()
def compile(p: Path):

    counts_dict = dict()
    po = PeptideObserver()

    g0 = p.glob("**/*e2g*tsv")
    g1 = filter(keep_e2g_file, g0)
    g2 = map(read, g1)


    # count = count_from_e2g
    # po.count_from_e2g()
    g3 = map(po.count_from_e2g, g2)
    # g4 = map(lambda x: update(d=counts_dict, counts=x), g3)

    #g4 = map(lambda x: update(d=counts_dict, counts=x), g3)
    g5 = list(g3)  #  walk

    peptides_out = po.make_summary()
    #peptides_out = shape_counts(po.peptidecounts_dict)

    if peptides_out.empty:
        logger.warning("No peptides found")
        return
    basepath = Path('out')
    if not basepath.exists():
        basepath.mkdir()
    name = basepath / f"peptide_counts_{date}.tsv"
    peptides_out.to_csv(name, sep="\t")
    logger.info(f"Saved {name}")
    # print(out)


if __name__ == "__main__":
    run_app()
