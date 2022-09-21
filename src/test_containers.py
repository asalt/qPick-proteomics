from containers import PeptideObserver
from dataclasses import dataclass
import pandas as pd
import pytest
from io import StringIO


@pytest.fixture
def counts():
    return pd.Series([1, 2, 3], index=["a", "b", "c"])


@pytest.fixture()
def testdata_e2g():
    s = "GeneID,PeptidePrint\n54617,gnnvpgnpk_ledsstqr_rqqeetnr_saidenqlsr_satsslr_sgtgfgesyslanpsir_vlspfapdyiqr\n51393,lfnavsr\n8722,nswgtdwgek_sdvpfwaik\n3073,agavaer_edipvnymk_elelvtk_evieyar_gfgedfk_gletfsqlvwk_gllldtsr_gsynpvthiytaqdvk_iqpdtiiqvwr_isygpdwk_ltsdltfayer_snpeiqdfmr_teiedfpr\n166,fttsdscdr_qvtapelnsiir"
    handle = StringIO(s)
    handle.seek(0)
    out = pd.read_csv(handle)
    out["GeneID"] = out["GeneID"].astype(str)
    return out


def test_create():
    po = PeptideObserver()
    assert po is not None


def test_update_dict(counts):
    po = PeptideObserver()
    po.update_counter_dict(po.peptidecounts_dict, counts)
    assert po.peptidecounts_dict["a"] == 1

    po.update_counter_dict(po.peptidecounts_dict, {"a": 3})
    assert po.peptidecounts_dict["a"] == 4


def test_e2g_counting(testdata_e2g):
    assert len(testdata_e2g) == 5
    po = PeptideObserver()

    po.count_from_e2g(testdata_e2g)
    assert len(po.peptidecounts_dict) > 1
    assert len(po.genecounts_dict) > 1
    assert len(po.genecounts_dict) == 5
    # counts = count(testdata_e2g)
    # assert len(counts) > 1
    # assert counts.index.name == "PeptidePrint"


def test_e2g_counting_multiple(testdata_e2g):
    po = PeptideObserver()
    po.count_from_e2g(testdata_e2g)

    gid0 = testdata_e2g.iloc[0]["GeneID"]
    assert len(po.genecounts_dict) == 5
    assert po.genecounts_dict.get(gid0) == len(
        testdata_e2g[testdata_e2g.GeneID == gid0]
    )

    po.count_from_e2g(testdata_e2g)
    assert len(po.genecounts_dict) == 5
    assert (
        po.genecounts_dict.get(gid0)
        == len(testdata_e2g[testdata_e2g.GeneID == gid0]) * 2
    )
