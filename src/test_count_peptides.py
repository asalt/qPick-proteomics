#
import pytest
from pathlib import Path
import pandas as pd
from io import StringIO

from count_peptides_from_e2g import keep_e2g_file

# from qpi
@pytest.fixture(scope="session")
def testpath(tmp_path_factory):
    # tmp_path_factory.mktemp("test")
    testdir = tmp_path_factory.getbasetemp() / "test"
    testdir.mkdir()
    (testdir / "11111_1_1_QUANT.tsv").touch()
    (testdir / "11111_1_1_QUAL.tsv").touch()
    return testdir


def test_testpath_fixture(testpath):
    p = Path(testpath)
    res = list(p.glob("*"))
    assert len(res) == 2


@pytest.fixture()
def testdata():
    s = "GeneID,PeptidePrint\n54617,gnnvpgnpk_ledsstqr_rqqeetnr_saidenqlsr_satsslr_sgtgfgesyslanpsir_vlspfapdyiqr\n51393,lfnavsr\n8722,nswgtdwgek_sdvpfwaik\n3073,agavaer_edipvnymk_elelvtk_evieyar_gfgedfk_gletfsqlvwk_gllldtsr_gsynpvthiytaqdvk_iqpdtiiqvwr_isygpdwk_ltsdltfayer_snpeiqdfmr_teiedfpr\n166,fttsdscdr_qvtapelnsiir"
    handle = StringIO(s)
    handle.seek(0)
    out = pd.read_csv(handle)
    return out


def test_keep_e2g_file():

    p = Path("QUANT.tsv")
    assert keep_e2g_file(p) == False
    p = Path("QUAL.tsv")
    assert keep_e2g_file(p) == False
    assert keep_e2g_file(p, require_regex_match=False) == True
    p = Path("11111_1_QUAL.tsv")
    assert keep_e2g_file(p) == True


def test_filter_files(testpath):
    p = Path(testpath)
    res = list(p.glob("*"))
    assert len(res) == 2
    res = list(filter(keep_e2g_file, res))
    assert len(res) == 1


# def test_count(testdata):
#     assert len(testdata) == 5
#     counts = count(testdata)
#     assert len(counts) > 1
#     assert counts.index.name == "PeptidePrint"
