from collections import OrderedDict, Mapping
import docutils.nodes
import docutils.parsers.rst
import docutils.utils
from copy import deepcopy
import re

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
        self.mdp_entries[name] = {
            "docname": self.arguments[0],
            "docstring": docstring,
            "options": options_dict,
            "default": default,
            "units": units
        }
        return []

    @staticmethod
    def get_default_and_units(content):
        if not content:
            return None, None
        firstline = content.splitlines()[0]
        words = firstline.split()
        defaults, units = [], []
        for s in words:
            s = s.strip()
            if s[0] in '(' and s[-1]==')' and not units:
                defaults.append(s[1:-1])
            elif s[0]=='[' and s[-1]==']':
                units.append(s[1:-1])
            elif s in ['/', '*', 'or']:
                units.append(s)
            else:
                break

        if len(defaults) == 0:
            try:
                float(units[0])
            except ValueError:
                default = None
            except IndexError:
                default = None
            else:
                default = units.pop(0)
        elif len(defaults) == 1:
            default = defaults[0]
        else:
            raise ValueError(f'Multiple candidate defaults in {content}')

        units = ''.join(units).replace('[', '').replace(']', '').replace('or', ' or ').strip()
        units = units if units else None

        return default, units

    def process_content_doc(self, content_doc):
        """Split a content doc tree into the docstring, default and units"""
        content = content_doc.astext()
        default, units = self.get_default_and_units(content)
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
