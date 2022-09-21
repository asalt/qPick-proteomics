# 2022 05 26 this is the current version of working peptidome code
# missing some features
from pyteomics import parser
from pyteomics.mass import mass
import csv
import pandas as pd
from collections import defaultdict
import numpy as np
import time
import RefProtDB
from RefProtDB import utils
import subprocess
import re

# date_ = time.strftime("%m/%d/%Y %H:%M:%S %p")
date_ = time.strftime("%Y%m%d%H%M%S%p")
# subprocess.call('RefProtDB download 9606 10090')

mass_dict = dict(mass.std_aa_mass)
mass_dict["X"] = 0  # np.mean([x for x in mass_dict.values()])
mass_dict["B"] = 0  # np.mean([mass_dict["D"], mass_dict["N"]])
mass_dict["J"] = 0  # np.mean([mass_dict["I"], mass_dict["L"]])
mass_dict["Z"] = 0  # np.mean([mass_dict["E"], mass_dict["Q"]])
mass_dict["U"] = 167.057  # SELENOCYSTEINE
mass_dict["O"] = 255.318  # pyrrolysine


def digest_protein(protein=None):
    digest = parser.cleave(protein, "[KR]", missed_cleavages=3, min_length=6)
    return digest


f = "/mnt/e/reference_databases/gpGrouper_Homo_sapiens_2020_03_24_refseq_GCF_000001405.39_GRCh38.p13_protein.fa"
d = RefProtDB.utils.fasta_dict_from_file(f, "specific")
ref_total = pd.DataFrame(d)
d = RefProtDB.utils.fasta_dict_from_file(f)
ref = pd.DataFrame(d)
ref["sequence"] = ref.sequence.str.upper()


f = "/mnt/e/reference_databases/gpGrouper_10090_2020_03_24_refseq_GCF_000001635.26_GRCm38.p6_protein.fa"
d = RefProtDB.utils.fasta_dict_from_file(f)
ref_mouse = pd.DataFrame(d)


targets = [
    "SFXN1",
    "TOM1L1",
    "DYNC1LI2",
    "NAA35",
    "CHRNB4",
    "ISOC1",
    "XRCC1",
    "HECTD4",
    "EMC10",
    "ITGA1",
    "TPD52L2",
    "LRRC47",
    "GLUD1",
    "PPP1R12A",
    "VPS35",
    "NUCB1",
    "ACTN1",
    "DHX36",
    "DIDO1",
    "USO1",
    "RDX",
    "BRD3",
    "F13B",
    "RPL36",
    "CISD2",
    "CUL1",
    "SLITRK6",
    "CDV3",
    "NCBP1",
    "CLINT1",
    "ZNF207",
    "SNAP29",
]

targets = [
    "CDK4",
    "POLD2",
    "SRPK1",
    "TRIM28",
]


digest_func = lambda x: parser.cleave(x, "[KR]", missed_cleavages=3, min_length=6)

subref = ref.query("symbol in @targets")

digest_out = subref.sequence.apply(digest_func).apply(tuple)
digest_out = digest_out.reset_index(name="peptide").rename(
    columns={"index": "protein_index"}
)
digest_out = digest_out.explode("peptide")
_nosequence = [x for x in subref if x != "sequence"]
peptidome = digest_out.merge(
    subref[_nosequence], left_on="protein_index", right_index=True
).reset_index(drop=True)


def get_surrounding(peptide, protein):
    start = protein.find(peptide)
    end = len(peptide)
    full_len = len(protein)

    before = list()
    after = list()

    for x in range(start - 4, start):
        if x < 0:
            before.append("-")
            continue
        before.append(protein[x])

    for x in range(end + 1, end + 5):
        if x > full_len:
            after.append("-")
            continue
        after.append(protein[x])

    return "".join(before), "".join(after)


def calculate_miscuts(seq):  # From Alex's code
    miscuts = seq.count("K") + seq.count("R")
    if not any(seq[-1] == x for x in "KR"):  # then at the C terminal
        return miscuts
    return miscuts - 1


def calculate_miscuts_withP(seq):  # From Alex's code
    miscuts = seq.count("KP") + seq.count("RP")
    return miscuts


peptidome["miscuts_trypsinP"] = peptidome.peptide.apply(calculate_miscuts)
peptidome["miscuts_trypsin"] = peptidome.miscuts_trypsinP - peptidome.peptide.apply(
    calculate_miscuts_withP
)
from tqdm import tqdm

before_all, after_all = list(), list()
for ix, row in tqdm(peptidome.iterrows(), total=len(peptidome)):
    prot_seq = subref.loc[row.protein_index].sequence
    before, after = get_surrounding(row.peptide, prot_seq)
    before_all.append(before)
    after_all.append(after)
peptidome["before"] = before_all
peptidome["after"] = after_all

unique_human = list()
for ix, row in tqdm(peptidome.iterrows(), total=len(peptidome)):
    gid = row.geneid
    peptide = row.peptide
    others = ref_total.query("geneid != @gid").sequence.unique()
    is_shared = any(peptide in y for y in others)
    unique_human.append(not is_shared)


unique_mouse = list()
for ix, row in tqdm(peptidome.iterrows(), total=len(peptidome)):
    gid = row.geneid
    peptide = row.peptide
    others = ref_mouse.query("geneid != @gid").sequence.unique()
    is_shared = any(peptide in y for y in others)
    unique_mouse.append(not is_shared)


peptidome["unique_human"] = unique_human
peptidome["unique_mouse"] = unique_mouse

peptidome.to_csv(f"peptidome_genesel_insilico_{date_}.tsv", sep="\t", index=False)
