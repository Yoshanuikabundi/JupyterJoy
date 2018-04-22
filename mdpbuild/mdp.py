from collections import OrderedDict
from .rstparser import parse_rst_file

class MDPBase():
    """Base class for creating MDP files"""
    @classmethod
    def create_subclass_from_doc(cls, name, docpath):
        """Create an MDP class from a mdp-options.rst file distributed with GROMACS"""
        mdp_entries, doc = parse_rst_file(docpath)
        attr = {}
        attr['options'] = mdp_entries
        attr['doc'] = doc
        return type(name, (cls,), attr)