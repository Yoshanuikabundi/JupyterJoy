from collections import OrderedDict, Mapping
import docutils.nodes
import docutils.parsers.rst
import docutils.utils
from copy import deepcopy
import re

class ReadOnlyDict(Mapping):
    def __init__(self, data):
        self.__data__ = data
    def __getitem__(self, key): 
        return self.__data__[key]
    def __len__(self):
        return len(self.__data__)
    def __iter__(self):
        return iter(self.__data__)
    def __repr__(self):
        return 'ReadOnlyDict({})'.format(repr(self.__data__))

_mdp_entries = OrderedDict()

def content_to_str(content):
    out = ""
    for string in content:
        out += string + "\n"
    return out

def process_name(name):
    return name.replace('_', '').replace('-', '').lower()


class MDPDirectiveBase(docutils.parsers.rst.Directive):
    required_arguments=1
    optional_arguments=4
    has_content = True
    _re_single_newlines = re.compile(r'\n(?!\n)')
    @property
    def mdp_entries(self):
        """Refers to self.odict if it exists, otherwise the module-wide mdp_entries odict

        Subclasses can specify cls.odict to allow them a particular odict to store mdp entries"""
        try:
            return self.odict
        except AttributeError:
            return _mdp_entries


class MDPDirective(MDPDirectiveBase):
    _re_def_and_units_line = re.compile(r"(?:\((?P<default>.+)\) )?(?P<units>\[.+\])?")

    def run(self):
        """Save the first argument of each parsed mdp directive and create an associated list for its values"""
        content = content_to_str(self.content)
        options_dict = OrderedDict()
        ThisMDPValueDirective = type('ThisMDPValueDirective', (MDPValueDirective,), {'odict': options_dict})
        with DirectivesRegistered(('mdp-value', ThisMDPValueDirective)):
            content_doc = parse_rst_string(content)
        docstring, default, units = self.process_content_doc(content_doc)
        name = process_name(self.arguments[0])
        self.mdp_entries[name] = ReadOnlyDict({
            "docname": self.arguments[0],
            "docstring": docstring, 
            "options": options_dict,
            "default": default,
            "units": units
        })
        return []

    def process_content_doc(self, content_doc):
        """Split a content doc tree into the docstring, default and units"""
        content = content_doc.astext()
        match = self._re_def_and_units_line.match(content)
        if match:
            default, units = match.group("default"), match.group("units")
        else:
            default, units = None, None
        docstring = re.sub(self._re_single_newlines , ' ', content)
        return docstring, default, units
    
class MDPValueDirective(MDPDirectiveBase):
    def run(self):
        """Read and process mdp-values into dicts with docstrings"""
        content = content_to_str(self.content)
        content_doc = parse_rst_string(content)
        docstring = re.sub(self._re_single_newlines , ' ', content_doc.astext())
        try:
            options_dict = self.odict
        except AttributeError:
            options_dict = list(self.mdp_entries.values())[-1]["options"]
        options_dict[self.arguments[0]] = docstring
        return []
    
def ignore_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    node = docutils.nodes.Text(text)
    return ([node],[])

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

docutils.parsers.rst.roles.register_local_role('mdp', ignore_role)
docutils.parsers.rst.roles.register_local_role('mdp-value', ignore_role)
docutils.parsers.rst.roles.register_local_role('ref', ignore_role)

def parse_rst_string(string):
    parser = docutils.parsers.rst.Parser()
    components = (docutils.parsers.rst.Parser,)
    settings = docutils.frontend.OptionParser(components=components).get_default_values()
    document = docutils.utils.new_document('<rst-doc>', settings=settings)
    parser.parse(string, document) 
    return document

def parse_rst_file(path):
    mdp_entries = OrderedDict()
    ThisMDPDirective = type('ThisMDPDirective', (MDPDirective,), {'odict': mdp_entries})
    new_directives = (
        ('mdp', ThisMDPDirective), 
        ('MDP', ThisMDPDirective)            
    )
    with DirectivesRegistered(*new_directives), open(path, 'r') as file:
        doc = parse_rst_string(file.read())
    return mdp_entries, doc



