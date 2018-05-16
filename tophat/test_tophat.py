import pytest
from .tophat import *

def test_Topology_molecules():
    strlist = [
        "  [ system ]   ",
        "; this is for gmx btw",
        "Test system - microbilayer!",
        "[ molecules ]  ; Molecule section time!!!",
        ";name ;number",
        "protein   1; the protein itself",
        "protein_with_a_very_long_name   1",
        "lipid_a 2",
        " lipid_b  5",
        ";lipid_c     4\n",
        "lipid_a      2",
        "lipid_b\t7",
        "",
        "water    12843;solvent!?!"
    ]
    top = Topology()
    top.read(strlist)

    compstr = '\n'.join([
        "[ molecules ]",
        ";name                          count",
        "protein                            1",
        "protein_with_a_very_long_name      1",
        "lipid_a                            2",
        "lipid_b                            5",
        "lipid_a                            2",
        "lipid_b                            7",
        "water                          12843"
    ])
    assert str(top.molecules) == compstr
    assert top.unparsed == [
        "[ system ]",
        "; this is for gmx btw",
        "Test system - microbilayer!",
    ]
    compstr = '\n'.join([
        "[ system ]",
        "; this is for gmx btw",
        "Test system - microbilayer!",
        "",
        "[ molecules ]",
        ";name                          count",
        "protein                            1",
        "protein_with_a_very_long_name      1",
        "lipid_a                            2",
        "lipid_b                            5",
        "lipid_a                            2",
        "lipid_b                            7",
        "water                          12843"
    ])
    assert str(top) == compstr

    top.molecules *= 10
    compstr = '\n'.join([
        "[ system ]",
        "; this is for gmx btw",
        "Test system - microbilayer!",
        "",
        "[ molecules ]",
        ";name                           count",
        "protein                            10",
        "protein_with_a_very_long_name      10",
        "lipid_a                            20",
        "lipid_b                            50",
        "lipid_a                            20",
        "lipid_b                            70",
        "water                          128430"
    ])
    assert str(top) == compstr

