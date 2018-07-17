from typing import (MutableSequence, Sequence, Callable, Any, Tuple,
                    Union, Dict, Iterable, overload)
from functools import singledispatch


class MolSectEntry(Sequence):
    def __init__(
            self,
            name: str,
            copies: int,
            parent: 'MoleculesSection'
            ) -> None:
        self.name = name
        self.copies = copies
        self.parent = parent

    @property
    def idx(self):
        for idx, entry in enumerate(self.parent._entries):
            if entry is self:
                return idx
        raise ValueError(f"MolSectEntry is not in parent")

    @overload
    def __getitem__(self, idx: int) -> Union[str, int]:
        pass
    @overload
    def __getitem__(self, idx: slice) -> Sequence[Union[str, int]]:
        pass
    def __getitem__(self, idx):
        items = (
            self.name,
            self.copies
        )
        return items[idx]

    def __len__(self) -> int:
        return 2

    def __repr__(self) -> str:
        try:
            idx = self.idx
        except ValueError:
            idx = None
        return f"MolSectEntry(name={self.name}, copies={self.copies})"

    def __eq__(self, other: object) -> bool:
        return tuple(self) == other

    def __hash__(self) -> int:
        return hash([self.name, self.copies])

    def __iadd__(self, other: int) -> 'MolSectEntry':
        self.copies += other
        return self

    def __isub__(self, other: int) -> 'MolSectEntry':
        self.copies -= other
        return self

    def __imul__(self, other: int) -> 'MolSectEntry':
        self.copies *= other
        return self

    def __itruediv__(self, other: int) -> 'MolSectEntry':
        self.copies //= other
        return self

    def __ifloordiv__(self, other: int) -> 'MolSectEntry':
        self.copies //= other
        return self

    def __radd__(self, other: float) -> float:
        return other + self.copies

    def __rsub__(self, other: float) -> float:
        return other - self.copies

    def __rmul__(self, other: float) -> float:
        return other * self.copies

    def __rtruediv__(self, other: float) -> float:
        return other / self.copies

    def __rfloordiv__(self, other: float) -> float:
        return other // self.copies


class _MolSectEntryList(list):
    def __repr__(self) -> str:
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
    def __init__(self, *args: Tuple[str, int]) -> None:
        self.__init_done = False

        self._entries: MutableSequence[MolSectEntry] = []
        self._name2entries: Dict[str, Sequence[MolSectEntry]] = {}

        for idx, (name, copies) in enumerate(args):
            self.append((name, copies))
        self.update_name2entries()

        self.__init_done = True

    def __repr__(self) -> str:
        entries = ', '.join(str(tuple(e)) for e in self._entries)
        return f"MoleculesSection({entries})"

    @overload
    def _key_to_entries(self,
                        key: slice,
                        recursing: bool) -> Sequence[MolSectEntry]:
        pass

    @overload
    def _key_to_entries(self, key: int, recursing: bool) -> MolSectEntry:
        pass

    @overload
    def _key_to_entries(self,
                        key: str,
                        recursing: bool
                        ) -> Union[Sequence[MolSectEntry], MolSectEntry]:
        pass

    def _key_to_entries(self, key, recursing=False):
        name = None
        if isinstance(key, str):
            name = key

        if name:
            try:
                entries = self._name2entries[name]
            except KeyError:
                raise KeyError(f"Molecule type {key} not in {self}")
            else:
                if isinstance(entries, MolSectEntry):
                    return entries
                else:
                    return _MolSectEntryList(entries)

        try:
            entries = (
                x for k in key for x in self._key_to_entries(k, recursing=True)
            )
        except TypeError:
            pass
        else:
            if recursing:
                raise TypeError(
                    f"MoleculesSection keys must be an integer, slice, or "
                    f"name, or an iterable of the above, not an iterable of "
                    f"{type(key)}s"
                )
            else:
                return _MolSectEntryList(entries)

        try:
            entries = self._entries[key]
        except TypeError:
            raise TypeError(
                f"MoleculesSection keys must be an integer, slice, or name, "
                f"or an iterable of the above, not {type(key)}")
        except IndexError:
            raise IndexError("MoleculesSection index out of range")
        else:
            if isinstance(entries, MolSectEntry):
                return entries
            else:
                return _MolSectEntryList(entries)

    @overload
    def __getitem__(self, key: int) -> MolSectEntry:
        pass

    @overload
    def __getitem__(self, key: slice) -> Sequence[MolSectEntry]:
        pass

    def __getitem__(self, key):
        return self._key_to_entries(key)

    @overload
    def __setitem__(self, key: int, tup: Tuple[str, int]) -> None:
        pass

    @overload
    def __setitem__(self, key: slice, tup: Iterable[Tuple[str, int]]) -> None:
        pass

    @overload
    def __setitem__(self, key: str, tup: int) -> None:
        pass

    @overload
    def __setitem__(self, key: str, tup: Iterable[int]) -> None:
        pass

    @overload
    def __setitem__(self, key: str, tup: Tuple[str, int]) -> None:
        pass

    @overload
    def __setitem__(self, key: str, tup: Iterable[Tuple[str, int]]) -> None:
        pass

    def __setitem__(self, key, tup):
        entries = self._key_to_entries(key)

        if isinstance(entries, MolSectEntry):
            entry = entries
            if isinstance(tup, int) and isinstance(key, str):
                name, copies = key, tup
            else:
                name, copies = tup
            if not isinstance(name, str):
                raise TypeError(f'name must be str, not {name}')
            if not isinstance(copies, int):
                raise TypeError(f'copies must be int, not {copies}')
            entry.name = name
            entry.copies = copies
        elif isinstance(tup, int):
            for e in entries:
                e.copies = tup
        elif len(entries) == len(tup):
            for e, t in zip(entries, tup):
                if isinstance(t, int) and isinstance(key, str):
                    name, copies = key, t
                else:
                    name, copies = t
                if not isinstance(name, str):
                    raise TypeError(f'name must be str, not {name}')
                if not isinstance(copies, int):
                    raise TypeError(f'copies must be int, not {copies}')
                e.name = name
                e.copies = copies
        else:
            raise ValueError(f'Cannot disambiguate key {key}')

    def __delitem__(self, key) -> None:
        entries = self._key_to_entries(key)
        if isinstance(entries, MolSectEntry):
            idx = self._entries.index(entries)
            del self._entries[idx]
        else:
            for entry in entries:
                idx = self._entries.index(entry)
                del self._entries[idx]
        self.update_name2entries()

    def __len__(self):
        return len(self._entries)

    def insert(self, idx, tup):
        name, copies = tup
        if not isinstance(name, str):
            raise TypeError(f'name must be str, not {name}')
        if not isinstance(copies, int):
            raise TypeError(f'copies must be int, not {copies}')

        entry = MolSectEntry(name, copies, self)
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

    def __imul__(self, other):
        for e in self:
            e *= other
        return self

    def __str__(self):
        out = ["[ molecules ]"]
        name_comment = ';name'
        copies_comment = 'count'
        name_lens = [len(str(e.name)) for e in list(self)]
        name_lens += [len(name_comment)]
        copies_lens = [len(str(e.copies)) for e in list(self)]
        copies_lens += [len(copies_comment)]
        fmt = f"{{0: <{max(name_lens)}}}  {{1: >{max(copies_lens)}}}"
        out.append(fmt.format(name_comment, copies_comment))
        for name, copies in self:
            out.append(fmt.format(name, copies))
        return '\n'.join(out)
