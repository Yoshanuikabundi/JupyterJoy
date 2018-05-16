from .molsect import MoleculesSection
from .moltype import MoleculeType
import re

_re_directive = re.compile(r'\[ +([a-zA-Z0-9_]+) +\]')


class Topology:
    def __init__(self, fname=None):
        self.unparsed = []
        self.molecules = MoleculesSection()
        if fname:
            with open(fname) as f:
                self.read(f)
        else:
            self.molecules = MoleculesSection()

    def __str__(self):
        return '\n'.join(self.unparsed + ['', str(self.molecules)])


    def read(self, f):
        current_directive = None
        for line in f:
            # Strip out comments
            line,*comments = line.split(sep=';', maxsplit=1)
            # Strip trailing and leading whitespace
            line = line.strip()

            # We need to set up our directives now
            casedict = {
                    'molecules': self._read_molecules
                }

            # Is this line a directive?
            match = _re_directive.match(line)
            if match and match.group(1).lower() in casedict:
                current_directive = match.group(1).lower()
                continue

            # And now we just run the appropriate directive function on the line
            readerfunc = casedict.get(current_directive, self._read_default)
            readerfunc(line, comments)

    def _read_molecules(self, line, _):
        if line:
            name, count = line.split()
            name = MoleculeType(name)
            count = int(count)
            self.molecules.append(name, count)

    def _read_default(self, line, comments):
        if comments:
            return self.unparsed.append(';'.join([line] + comments))
        else:
            return self.unparsed.append(line)

    







