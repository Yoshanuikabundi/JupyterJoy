from .molsect import MoleculesSection
from .moltype import MoleculeType
import re
from ..listviews import REListView

_re_directive = re.compile(r'\[ +([a-zA-Z0-9_]+) +\]')
_re_includes = re.compile(r'\#include .+')
_re_defines = re.compile(r'\#define .+')

class Topology:
    def __init__(self, fname=None):
        self.unparsed = []
        self.name = ''
        self.molecules = MoleculesSection()

        if fname:
            with open(fname) as f:
                self.read(f)

    def __str__(self):
        if not self.name:
            raise ValueError("GROMACS topology must have a name!")

        out = list(self.unparsed)
        out += ['']
        out += ['[ system ]']
        out += ['; name']
        out += [self.name]
        out += ['']
        out += [str(self.molecules)]
        out += [''] # I am very proud of this line
        
        out = [i for n,i in enumerate(out) if i != "" or out[n-1] != ""] 
        return '\n'.join(out)

    @property
    def includes(self):
        return REListView(self.unparsed, _re_includes)

    @property
    def defines(self):
        return REListView(self.unparsed, _re_defines)

    def read(self, f):
        current_directive = None
        for line in f:
            # Strip out comments
            line,*comments = line.split(sep=';', maxsplit=1)
            # Strip trailing and leading whitespace
            line = line.strip()
            comments = [s.strip() for s in comments]

            # We need to set up our directives now
            casedict = {
                    'molecules': self._read_molecules,
                    'system': self._read_system
                }

            # Is this line a directive?
            match = _re_directive.match(line)
            if match and match.group(1).lower() in casedict:
                current_directive = match.group(1).lower()
                continue

            # And now we just run the appropriate directive function on the line
            readerfunc = casedict.get(current_directive, self._read_default)
            readerfunc(line, comments)

    def _read_system(self, line, _):
        if line[:1] == "#":
            raise ValueError("Hashcommand after [ system ] directive not supported")
        if self.name and line and self.name != line:
            raise ValueError(f'System name defined ambiguously: "{self.name}" and "{line}"')
        elif line:
            self.name = line

    def _read_molecules(self, line, _):
        if line[:1] == "#":
            raise ValueError("Hashcommand after [ system ] directive not supported")
        if line:
            name, count = line.split()
            name = MoleculeType(name)
            count = int(count)
            self.molecules.append(name, count)

    def _read_default(self, line, comments):
        if comments:
            self.unparsed.append(' ; '.join([line] + comments))
            return
        else:
            self.unparsed.append(line)
            return

    







