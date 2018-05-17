from collections.abc import MutableSequence
from abc import abstractmethod
from contextlib import suppress
import re

class ListView(MutableSequence):
    def __init__(self, alist, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.alist = alist

    @abstractmethod
    def map_view(self):
        """Map the view to the desired elements in alist

        map_view must return a sequence of pairs of original alist
        indices and their corresponding values."""
        return enumerate(self.alist)

    @property
    def _filtered(self):
        return list(self.map_view())

    def __getitem__(self, key):
        return self._filtered[key][1]

    def __setitem__(self, key, value):
        idx = self._filtered[key][0]
        self.alist[idx] = value

    def __delitem__(self, key):
        idx = self._filtered[key][0]
        del self.alist[idx]

    def __len__(self):
        return len(self._filtered)

    def __eq__(self, other):
        return list(self) == other

    def insert(self, key, value):
        if key == len(self):
            idx = self._filtered[-1][0] + 1
        else:
            idx = self._filtered[key][0]
        self.alist.insert(idx, value)

class ValidateListView(ListView):
    @abstractmethod
    def validate(self, value):
        """Only values that return true here are included in the view"""
        return True

    def __setitem__(self, key, value):
        if not self.validate(value):
            raise ValueError(f'{value} is not a valid element of {self}')
        super().__setitem__(key, value)
        
    def insert(self, key, value):
        if not self.validate(value):
            raise ValueError(f'{value} is not a valid element of {self}')
        super().insert(key, value)

    def map_view(self):
        return ((i,v) for i,v in super().map_view() if self.validate(v))

class MutateListView(ListView):
    @abstractmethod
    def mutate(self, value):
        """Values in the view are mutated by this method"""
        return value

    @abstractmethod
    def demutate(self, value):
        """Values modified in the view are mutated by this method"""
        return value

    def __setitem__(self, key, value):
        super().__setitem__(key, self.demutate(value))
        
    def insert(self, key, value):
        super().insert(key, self.demutate(value))

    def map_view(self):
        return ((i,self.mutate(v)) for i,v in super().map_view())


class SliceListView(ListView):
    def __init__(self, alist, *args, **kwargs):
        if 'stop' in kwargs:
            start = kwargs.pop('start', 0)
            stop  = kwargs.pop('stop')
            step  = kwargs.pop('step', 1)
        elif len(args) >= 3:
            start, stop, step, *args = args
        elif len(args) == 2:
            start, stop = args
            step, args  = 1, []
        elif len(args) == 1:
            stop,  = args
            start, step, args  = 0, 1, []
        else:
            raise ValueError("Need more arguments!")

        super().__init__(alist, *args, **kwargs)

        self.slice = slice(start, stop, step)

    def map_view(self):
        indices = self.slice.indices(len(self.alist))
        values = self.alist.__getitem__(self.slice)
        return zip(range(*indices), values)

    def __repr__(self):
        return f'SliceListView({self.alist}, {self.slice.start}, {self.slice.stop}, {self.slice.step})'


class REListView(ValidateListView):
    def __init__(self, alist, regex, *args, **kwargs):
        super().__init__(alist, *args, **kwargs)
        if isinstance(regex, str):
            regex = re.compile(regex)
        self.regex = regex

    def validate(self, value):
        return self.regex.match(value) and super().validate(value)

    def __repr__(self):
        return f'REListView({self.alist}, {self.regex})'