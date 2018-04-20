from collections import OrderedDict
import docutils.nodes
import docutils.parsers.rst
import docutils.utils
from copy import deepcopy

mdp_entries = OrderedDict()
class MDPDirective(docutils.parsers.rst.Directive):
    required_arguments=1
    optional_arguments=4
    def run(self):
        mdp_entries[self.arguments[0]] = []
        if len(self.arguments) != 1:
            print(self.arguments)
        return []
    
class MDPValueDirective(docutils.parsers.rst.Directive):
    required_arguments=1
    optional_arguments=4
    def run(self):
        list(mdp_entries.values())[-1].append(self.arguments[0])
        if len(self.arguments) != 1:
            print(self.arguments)
        return []
    
def not_a_role(*args, **kwargs):
    return ([],[])

docutils.parsers.rst.directives.register_directive('mdp', MDPDirective)
docutils.parsers.rst.directives.register_directive('mdp-value', MDPValueDirective)

docutils.parsers.rst.roles.register_canonical_role('mdp', not_a_role)
docutils.parsers.rst.roles.register_canonical_role('mdp-value', not_a_role)
docutils.parsers.rst.roles.register_canonical_role('ref', not_a_role)

def parse_rst(path):
    parser = docutils.parsers.rst.Parser()
    components = (docutils.parsers.rst.Parser,)
    settings = docutils.frontend.OptionParser(components=components).get_default_values()
    document = docutils.utils.new_document('<rst-doc>', settings=settings)
    with open(path, 'r') as file:
        for line in file:
            parser.parse(line, document)
    return document

class MDPBase():
    """Base class for creating MDP files"""
    @classmethod
    def create_subclass_from_doc(cls, name, docpath):
        """Create an MDP class from a mdp-options.rst file distributed with GROMACS"""
        parse_rst(docpath)
        attr = deepcopy(mdp_entries)
        return type(name, (cls,), attr)

