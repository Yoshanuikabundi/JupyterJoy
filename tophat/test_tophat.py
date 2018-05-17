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
    compstr = '\n'.join([
        "[ system ]",
        "; name",
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
        "water                          12843",
        ""
    ])
    assert str(top) == compstr

    top.molecules *= 10
    compstr = '\n'.join([
        "[ system ]",
        "; name",
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
        "water                          128430",
        ""
    ])
    assert str(top) == compstr


def test_Topology_hashcommands():
    strlist = [
        '#include "./martini_v2.2refP.itp" ; comment',
        "#define RUBBERBANDS",
        "#include './martini_protein.itp'",
        "",
        "#include martini_v2.0_ions.itp",
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

    assert top.includes == [
        '#include "./martini_v2.2refP.itp" ; comment',
        "#include './martini_protein.itp'",
        "#include martini_v2.0_ions.itp"
    ]
    assert top.defines == ["#define RUBBERBANDS"]
    top.includes.append('#include martini_v2.0_sugars.itp')

    compstr = '\n'.join([
        '#include "./martini_v2.2refP.itp" ; comment',
        "#define RUBBERBANDS",
        "#include './martini_protein.itp'",
        "",
        "#include martini_v2.0_ions.itp",
        "#include martini_v2.0_sugars.itp",
        "",
        "[ system ]",
        "; name",
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
        "water                          12843",
        ""
    ])
    assert str(top) == compstr


def test_Topology_system():
    strlist = [
        '#include "./martini_v2.2refP.itp"; comment',
        "#define RUBBERBANDS",
        "#include './martini_protein.itp'",
        "",
        "#include martini_v2.0_ions.itp",
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

    with pytest.raises(ValueError):
        str(top)

    top.name = "Test system - microbilayer!"

    compstr = '\n'.join([
        '#include "./martini_v2.2refP.itp" ; comment',
        "#define RUBBERBANDS",
        "#include './martini_protein.itp'",
        "",
        "#include martini_v2.0_ions.itp",
        "",
        "[ system ]",
        "; name",
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
        "water                          12843",
        ""
    ])
    assert str(top) == compstr

    strlist = [
        '#include "./martini_v2.2refP.itp" ; comment',
        "#define RUBBERBANDS",
        "#include './martini_protein.itp'",
        "",
        "#include martini_v2.0_ions.itp",
        "[ system ]",
        "; name",
        "Test system - microbilayer! - but inconsistent name",
        "",
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
    with pytest.raises(ValueError):
        top.read(strlist)

