"""This script test the diff_exp_annotations.py with various conditions.

Author: Jia Yi Terri Shen
Date: November 2019
"""

import pytest
from diff_exp_annotations import (transcript_protein_dict, gene_go_dict,
                                  go_name_dict, wrtie_annotations)


def test_transcript_protein_dict_blast():
    """Test if the filename provided is .outfmt6 file."""

    transcript_to_protein = transcript_protein_dict("blastp.outfmt6")
    assert transcript_to_protein


def test_transcript_protein_dict_notblast():
    """Test if the filename provided is not .outfmt6 file."""

    transcript_to_protein = transcript_protein_dict("gene_association_subset.gaf")
    assert not transcript_to_protein


def test_transcript_protein_dict_none():
    """Test if the file is not provided."""

    with pytest.raises(TypeError):
        transcript_to_protein = transcript_protein_dict()


def test_gene_go_dict_gaf():
    """Test if the filename provided is .gaf file."""

    gene_to_go = gene_go_dict("gene_association_subset.gaf")
    assert gene_to_go


def test_gene_go_dict_notgaf():
    """Test if the filename provided is not .gaf file."""

    gene_to_go = gene_go_dict("go-basic.obo")
    assert not gene_to_go


def test_gene_go_dict_none():
    """Test if the file is not provided."""

    with pytest.raises(TypeError):
        gene_to_go = gene_go_dict()


def test_go_name_dict_obo():
    """Test if the filename provided is .obo file."""
    go_to_desc = go_name_dict("go-basic.obo")
    assert go_to_desc


def test_go_name_dict_notobo():
    """Test if the filename provided is not .obo file."""

    go_to_desc = go_name_dict("blastp.outfmt6")
    assert not go_to_desc


def test_go_name_dict_none():
    """Test if the file is not provided."""

    with pytest.raises(TypeError):
        go_to_desc = go_name_dict()


def test_wrtie_annotations_nodicts():
    """Test if one of the file or dictionary is not provided"""

    with pytest.raises(TypeError):
        report_file = wrtie_annotations("diffExpr.P1e-3_C2.matrix",
                                        "gene_association_subset.gaf",
                                        "diffExpr.P1e-3_C2.matrix",
                                        "go-basic.obo")
