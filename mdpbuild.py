from collections import OrderedDict
import docutils.nodes
import docutils.parsers.rst
import docutils.utils
from copy import deepcopy

_mdp_entries = OrderedDict()

class MDPDirectiveBase(docutils.parsers.rst.Directive):
    required_arguments=1
    optional_arguments=4
    @property
    def mdp_entries(self):
        """Refers to self.odict if it exists, otherwise the module-wide mdp_entries odict

        Subclasses can specify cls.odict to allow them a particular odict to store mdp entries"""
        try:
            return self.odict
        except AttributeError:
            return _mdp_entries


class MDPDirective(MDPDirectiveBase):
    def run(self):
        """Save the first argument of each parsed mdp directive and create an associated list for its values"""
        self.mdp_entries[self.arguments[0]] = []
        return []
    
class MDPValueDirective(MDPDirectiveBase):
    def run(self):
        """Save the first argument of each parsed mdp-value directive to the last mdp directive parsed"""
        list(self.mdp_entries.values())[-1].append(self.arguments[0])
        return []
    
def not_a_role(*args, **kwargs):
    return ([],[])

class DirectivesRegistered():
    """Context manager for registration of directives"""
    def __init__(self, *directives):
        """Set up the directives to register and unregister

        Arguments should be tuples of name,class pairs for each new directive"""
        self.directives = directives
        self.cached_registry = None
    def __enter__(self):
        """Cache the directive registry and register new directives"""
        self.cached_registry = deepcopy(docutils.parsers.rst.directives._directives)
        for name,klass in self.directives:
            docutils.parsers.rst.directives.register_directive(name, klass)
    def __exit__(self, exception_type, exception_value, traceback):
        """Restore cached directive registry"""
        docutils.parsers.rst.directives._directives = self.cached_registry

# docutils.parsers.rst.directives.register_directive('mdp', MDPDirective)
# docutils.parsers.rst.directives.register_directive('mdp-value', MDPValueDirective)

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
        this_mdp_entries = OrderedDict()
        ThisMDPDirective = type('ThisMDPDirective', (MDPDirective,), {'odict': this_mdp_entries})
        ThisMDPValueDirective = type('ThisMDPValueDirective', (MDPValueDirective,), {'odict': this_mdp_entries})
        with DirectivesRegistered(('mdp', ThisMDPDirective), ('mdp-value', ThisMDPValueDirective)):
            parse_rst(docpath)

        attr = this_mdp_entries
        return type(name, (cls,), attr)

