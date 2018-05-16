from collections.abc import MutableSequence, Sequence
from .moltype import MoleculeType

def _allow_value_tuple_idx(func):
    """Decorator to allow name and count arguments to be given together as a tuple"""
    def out(self, idx, *args, **kwargs):
        if (isinstance(args[0], tuple) and len(args[0]) == 2) or isinstance(args[0], MolSectEntry):
            (name, count), *args = args
        elif isinstance(args[0], _MolSectEntryList) and self[idx] == args[0]:
            # Probably just signals that idx refers to multiple things 
            # that have just been set by __iadd__ or similar
            # TODO: allow __setitem__ to take iterables of entries
            return
        else:
            name, count, *args = args
        return func(self, idx, name, count, *args, **kwargs)
    return out

def _allow_value_tuple(func):
    """Decorator to allow name and count arguments to be given together as a tuple"""
    def out(self, *args, **kwargs):
        if isinstance(args[0], tuple) and len(args[0]) == 2:
            (name, count), *args = args
        elif isinstance(args[0], MolSectEntry):
            (name, _, count), *args = args
        else:
            name, count, *args = args
        return func(self, name, count, *args, **kwargs)
    return out

class MolSectEntry(Sequence):
    def __init__(self, name, count, parent):
        self.name = name
        self.count = count
        self.parent = parent

    @property
    def idx(self):
        for idx, entry in enumerate(self.parent._entries):
            if entry is self:
                return idx
        raise ValueError(f"MolSectEntry is not in parent")

    def __getitem__(self, idx):
        items = (
            self.name,
            self.count
        )
        return items[idx]

    def __len__(self):
        return 3

    def __repr__(self):
        try:
            idx = self.idx
        except ValueError:
            idx = None
        return f"MolSectEntry(name={self.name}, count={self.count})"

    def __eq__(self, other):
        if isinstance(other, MolSectEntry):
            return self is other
        else:
            return tuple(self) == other

    def __hash__(self):
        return hash([self.name, self.count])

    def __iadd__(self, other):
        self.count += other
        return self
    def __isub__(self, other):
        self.count -= other
        return self
    def __imul__(self, other):
        self.count *= other
        return self
    def __itruediv__(self, other):
        self.count /= other
        return self
    def __ifloordiv__(self, other):
        self.count //= other
        return self

    def __radd__(self, other):
        return other + self.count
    def __rsub__(self, other):
        return other - self.count
    def __rmul__(self, other):
        return other * self.count
    def __rtruediv__(self, other):
        return other / self.count
    def __rfloordiv__(self, other):
        return other // self.count

class _MolSectEntryList(list):
    def __repr__(self):
        return f"_MolSectEntryList({super().__repr__()})"
    def __iadd__(self, other):
        for e in self:
            e += other
        return self
    def __isub__(self, other):
        for e in self:
            e -= other
        return self
    def __imul__(self, other):
        for e in self:
            e *= other
        return self
    def __itruediv__(self, other):
        for e in self:
            e /= other
        return self
    def __ifloordiv__(self, other):
        for e in self:
            e //= other
        return self


class MoleculesSection(MutableSequence):
    """Store the number and order of molecules in a topology"""
    def __init__(self, *args):
        self.__init_done = False

        self._entries = []
        self._name2entries = {}

        for idx,(name,count) in enumerate(args):
            self.append(name, count)
        self.update_name2entries()

        self.__init_done = True

    def __repr__(self):
        return f"MoleculesSection({', '.join(str(tuple(e)) for e in self._entries)})"

    def _key_to_entries(self, key, recursing=False):
        name = None
        if isinstance(key, MoleculeType):
            name = key.name
        if isinstance(key, str):
            name = key

        if name:
            try:
                entries = self._name2entries[name]
            except KeyError:
                raise KeyError(f"MoleculeType {key} not in {self}")
            else:
                if isinstance(entries, MolSectEntry):
                    return entries
                else:
                    return _MolSectEntryList(entries)
        
        try:
            entries = (x for k in key for x in self._key_to_entries(k, recursing=True))
        except TypeError:
            pass
        else:
            if recursing:
                raise TypeError(f"MoleculesSection keys must be an integer, slice, or name, or an iterable of the above, not an iterable of {type(key)}s")
            else:
                return _MolSectEntryList(entries)

        try:
            entries = self._entries[key]
        except TypeError:
            raise TypeError(f"MoleculesSection keys must be an integer, slice, or name, or an iterable of the above, not {type(key)}")
        except IndexError:
            raise IndexError("MoleculesSection index out of range")
        else:
            if isinstance(entries, MolSectEntry):
                return entries
            else:
                return _MolSectEntryList(entries)
        
    def _key_to_idx(self, key):
        entries = self[key]
        try:
            return (entry.idx for entry in entries)
        except TypeError:
            return entries.idx

    def __getitem__(self, key):
        return self._key_to_entries(key)

    @_allow_value_tuple_idx
    def __setitem__(self, key, name, count):
        if not (isinstance(name, MoleculeType) or isinstance(name, str)):
            raise TypeError(f'name must be MoleculeType or str, not {name}')
        if not isinstance(count, int):
            raise TypeError(f'count must be int, not {count}')

        entries = self._key_to_entries(key)

        entry = None
        if isinstance(entries, MolSectEntry):
            entry = entries
        elif len(entries) == 1:
            entry = entries[0]

        if entry:
            entry.name = name
            entry.count = count
            return

        raise ValueError(f'Cannot disambiguate key {key}')

    def __delitem__(self, key):
        entries = self._key_to_entries(key)
        for entry in entries:
            del self._entries[entry.idx]
        self.update_name2idx()


    def __len__(self):
        return len(self._entries)


    @_allow_value_tuple_idx
    def insert(self, idx, name, count):
        if not isinstance(name, MoleculeType):
            raise TypeError(f'count must be MoleculeType, not {name}')
        if not isinstance(count, int):
            raise TypeError(f'count must be int, not {count}')

        entry = MolSectEntry(name, count, self)
        self._entries.insert(idx, entry)

        if self.__init_done:
            self.update_name2entries()

    def update_name2entries(self):
        self._name2entries.clear()
        for entry in self._entries:
            try:
                self._name2entries[str(entry.name)].append(entry)
            except KeyError:
                self._name2entries[str(entry.name)] = [entry]

    @_allow_value_tuple
    def append(self, name, count, *args, **kwargs):
        super().append((name, count), *args, **kwargs)

    @_allow_value_tuple
    def remove(self, name, count, *args, **kwargs):
        super().remove((name, count), *args, **kwargs)

    @_allow_value_tuple
    def __contains__(self, name, count, *args, **kwargs):
        super().__contains__((name, count), *args, **kwargs)

    @_allow_value_tuple
    def index(self, name, count, *args, **kwargs):
        super().index((name, count), *args, **kwargs)

    @_allow_value_tuple
    def count(self, name, count, *args, **kwargs):
        super().count((name, count), *args, **kwargs)

    def __imul__(self, other):
        for e in self:
            e *= other
        return self

    def __str__(self):
        out = ["[ molecules ]"]
        name_comment = ';name'
        count_comment = 'count'
        name_lens = [len(str(e.name)) for e in list(self)] + [len(name_comment)]
        count_lens = [len(str(e.count)) for e in list(self)] + [len(count_comment)]
        fmt = f"{{0: <{max(name_lens)}}}  {{1: >{max(count_lens)}}}"
        out.append(fmt.format(name_comment, count_comment))
        for name, count in self:
            out.append(fmt.format(name, count))
        return '\n'.join(out)


