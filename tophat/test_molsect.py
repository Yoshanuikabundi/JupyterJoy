import pytest
from .molsect import *


def test_MoleculesSection_sequence():
    protein = "protein"
    lipid_a = "lipid_a"
    lipid_b = "lipid_b"
    water = "water"

    lst = [
        (protein,   1),
        (lipid_a,   4),
        (lipid_b,  12),
        (water, 12843)
    ]

    mol_sect1 = MoleculesSection(*lst)

    mol_sect2 = MoleculesSection(
        (protein,   1),
        (lipid_a,   4),
        (lipid_b,  12),
        (water, 12843)
    )

    # Orders must be preserved!
    for (ak, av), (bk, bv), (lk, lv) in zip(mol_sect1, mol_sect2, lst):
        assert ak == bk == lk
        assert av == bv == lv

    # Can have duplicate entries
    lst3 = [
        (protein,   1),
        (lipid_a,   2),
        (lipid_b,   5),
        (lipid_a,   2),
        (lipid_b,   7),
        (water, 12843)
    ]
    mol_sect3 = MoleculesSection(*lst3)

    for (ak, av), (lk, lv) in zip(mol_sect3, lst3):
        assert ak == lk
        assert av == lv

    # Indexing by integers returns the nth item
    assert mol_sect3[0] == (protein, 1)
    assert mol_sect3[0].idx == 0
    assert mol_sect3[3] == (lipid_a, 2)
    assert mol_sect3[3].idx == 3
    with pytest.raises(IndexError):
        mol_sect3[6]
    assert mol_sect3[-1] == (water, 12843)
    assert mol_sect3[-1].idx == 5
    assert list(mol_sect3[::-1]) == [
        (water, 12843),
        (lipid_b,   7),
        (lipid_a,   2),
        (lipid_b,   5),
        (lipid_a,   2),
        (protein,   1)
    ]
    assert [e.idx for e in mol_sect3[::-1]] == list(range(5, -1, -1))

    # Indexing by string or MoleculeType returns those items
    assert list(mol_sect3[protein]) == [(protein, 1)]
    assert list(e.idx for e in mol_sect3[protein]) == [0]
    assert list(mol_sect3['protein']) == [(protein, 1)]
    assert list(e.idx for e in mol_sect3['protein']) == [0]
    assert list(mol_sect3[lipid_b]) == [(lipid_b, 5), (lipid_b, 7)]
    assert list(e.idx for e in mol_sect3[lipid_b]) == [2, 4]
    assert list(mol_sect3['lipid_b']) == [(lipid_b, 5), (lipid_b, 7)]
    assert list(e.idx for e in mol_sect3['lipid_b']) == [2, 4]
    with pytest.raises(KeyError):
        mol_sect3['Na']
    assert list(mol_sect3[lipid_a, lipid_b]) == [
        (lipid_a, 2),
        (lipid_a, 2),
        (lipid_b, 5),
        (lipid_b, 7)
    ]
    assert [e.idx for e in mol_sect3[lipid_a, lipid_b]] == [1, 3, 2, 4]

    mol_sect3[protein] = 2
    mol_sect3[lipid_a] = 2, 3
    mol_sect3[lipid_b] = 6
    assert list(mol_sect3) == [
        (protein,   2),
        (lipid_a,   2),
        (lipid_b,   6),
        (lipid_a,   3),
        (lipid_b,   6),
        (water, 12843)
    ]

    # Can insert new lines
    na = 'NA'
    cl = 'CL'
    mol_sect3.append((cl, 32))
    mol_sect3.insert(-1, (na, 32))
    assert list(mol_sect3) == [
        (protein,   2),
        (lipid_a,   2),
        (lipid_b,   6),
        (lipid_a,   3),
        (lipid_b,   6),
        (water, 12843),
        (na,       32),
        (cl,       32)
    ]
    assert [e.idx for e in mol_sect3] == list(range(0, 8))


def test_MoleculesSection_numeric():
    protein = "protein"
    lipid_a = "lipid_a"
    lipid_b = "lipid_b"
    water = "water"
    lst3 = [
        (protein,   1),
        (lipid_a,   2),
        (lipid_b,   5),
        (lipid_a,   2),
        (lipid_b,   7),
        (water, 12843)
    ]
    mol_sect3 = MoleculesSection(*lst3)

    # Can add, subtract and multiply, even from indexed stuff
    assert mol_sect3[0] == (protein, 1)
    assert mol_sect3[4] == (lipid_b, 7)
    assert mol_sect3[2] == (lipid_b, 5)
    mol_sect3[0] += 4
    assert mol_sect3[0] == (protein, 5)
    mol_sect3[lipid_b] -= 1
    assert mol_sect3[4] == (lipid_b, 6)
    assert mol_sect3[2] == (lipid_b, 4)
    mol_sect3 *= 10
    assert mol_sect3[:] == [
        (protein,   50),
        (lipid_a,   20),
        (lipid_b,   40),
        (lipid_a,   20),
        (lipid_b,   60),
        (water, 128430)
    ]


def test_MoleculesSection_str():
    protein = "protein"
    lipid_a = "lipid_a"
    lipid_b = "lipid_b"
    water = "water"
    lst3 = [
        (protein,   1),
        (lipid_a,   2),
        (lipid_b,   5),
        (lipid_a,   2),
        (lipid_b,   7),
        (water, 12843)
    ]
    mol_sect3 = MoleculesSection(*lst3)


    compstr = '\n'.join([
        "[ molecules ]",
        ";name    count",
        "protein      1",
        "lipid_a      2",
        "lipid_b      5",
        "lipid_a      2",
        "lipid_b      7",
        "water    12843"
    ])
    assert str(mol_sect3) == compstr
