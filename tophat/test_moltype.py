import pytest
from .moltype import *
from traitlets import TraitError

def test_MoleculeName():
    molname = MoleculeName()
    molname.validate(None, "protein_a")
    with pytest.raises(TraitError):
        molname.validate(None, "protein a")
    with pytest.raises(TraitError):
        molname.validate(None, "protein\ta")

def test_MoleculeType():
    moltype = MoleculeType('protein_a')
    assert moltype.name == 'protein_a'
    assert isinstance(moltype.name, str)

    with pytest.raises(TraitError):
        moltype = MoleculeType('protein a')
    with pytest.raises(TraitError):
        moltype = MoleculeType('protein\ta')
    assert moltype.name == 'protein_a'
    assert isinstance(moltype.name, str)

    moltype2 = MoleculeType('protein_a')
    moltype3 = MoleculeType('protein_b')

    assert moltype2 is not moltype
    assert moltype2 == moltype
    assert moltype == moltype2
    assert moltype2 != moltype3
    assert moltype3 != moltype2
    assert moltype3 != moltype