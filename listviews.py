from collections.abc import MutableSequence
from abc import abstractmethod, ABCMeta
from contextlib import suppress
import re

class ListViewMeta(ABCMeta):
    def __init__(cls, *args, **kwargs):
        super().__init__(*args, **kwargs)
        cls._SubView = None

    @property
    def SubView(cls):
        if cls._SubView is None:
            cls._SubView = type(
                f"{cls.__name__}.SubView", 
                (cls, SubView),
                {}
            )

        return cls._SubView

class ListView(MutableSequence, metaclass=ListViewMeta):
    def __init__(self, alist, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.alist = alist

    @property
    def SubView(self):
        return self.__class__.SubView

    @property
    @abstractmethod
    def indices(self):
        """Return a Sequence (eg, list, tuple) of indices in the view"""
        return tuple(range(len(self.alist)))

    def create_subview(self, idx, *args, **kwargs):
        """Create a view to return for slicing

        This is a bit tricky. ListView classes automatically create a
        class for their own subviews called cls.SubView. These are
        returned when a slice is requested of the ListView. cls.SubView
        is a subclass of the ListView class itself as well as 
        ExplicitListView so it can be given arbitrary indices. What this
        means is that when a ListView creates a subview , it needs to be 
        able to provide all the arguments it would need to create a new 
        instance of its own class. This method lets it do that.

        When a subclass of ListView changes the argument signature of its
        __init__method, it must also redefine this method.
        For an example, see REListView."""
        return self.SubView(self.alist, indices=idx, *args, **kwargs)

    def __getitem__(self, key):
        idx = self.indices[key]
        try:
            return self.alist[idx]
        except TypeError:
            return self.create_subview(idx)

    def __setitem__(self, key, value):
        # print(f"__setitem__({repr(key)}, {repr(value)})")
        idx = self.indices[key]
        try:
            self.alist[idx] = value
            return
        except TypeError:
            view = self.create_subview(idx)
            value = list(value)

        assert len(idx) == len(view)
        assert isinstance(key, slice)

        start, stop, step = key.start, key.stop, key.step
        if len(value) == len(idx):
            for i,v in enumerate(value):
                view[i] = v
        elif step not in (1, None):
            raise ValueError(f"attempt to assign sequence of size {len(value)} to extended slice of size {len(idx)}")
        else:
            # Standard slice
            if start is None:
                start = 0
            if stop is None:
                stop = len(self)

            # Express indices positively between 0 and len(self)
            if start < -len(self):
                start = 0
            elif start < 0:
                start += len(self)
            if stop < -len(self):
                stop = 0
            if stop < 0:
                stop += len(self)
            start = min(start, len(self))
            stop = min(stop, len(self))

            assert 0 <= start <= len(self)
            assert 0 <= stop <= len(self)

            del self[key]
            for v in reversed(value):
                if len(self) == 0:
                    if start < 0:
                        idx = -1
                    else:
                        idx = 0
                elif start >= len(self):
                    idx = self.indices[-1] + 1
                else:
                    idx = self.indices[start]

                self.alist.insert(idx, v)

    def __delitem__(self, key):
        idx = self.indices[key]
        if isinstance(idx, int):
            del self.alist[idx]
        else:
            for n,i in enumerate(sorted(set(idx))):
                del self.alist[i-n]
        return idx

    def __len__(self):
        return len(self.indices)

    def __eq__(self, other):
        return list(self) == other

    def insert(self, key, value):
        if len(self) == 0:
            if key < 0:
                idx = -1
            else:
                idx = 0
        elif key >= len(self):
            idx = self.indices[-1] + 1
        else:
            idx = self.indices[key]

        self.alist.insert(idx, value)
        return idx

    def __repr__(self):
        clsname = self.__class__.__name__
        values = type(self.indices)(self)
        return f"<{clsname} object bearing {values}>"

class ExplicitListView(ListView):
    def __init__(self, alist, indices, *args, **kwargs):
        super().__init__(alist, *args, **kwargs)
        self._indices = list(indices)

    @property
    def indices(self):
        return self._indices

    @indices.setter
    def indices(self, other):
        self._indices = list(other)

    @property
    def SubView(self):
        return self.__class__

    def __getitem__(self, key):
        if isinstance(key, slice):
            raise ValueError("ExplicitListView does not support slicing (yet?)")
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        if isinstance(key, slice):
            raise ValueError("ExplicitListView does not support slicing (yet?)")
        super().__setitem__(key, value)

    def __delitem__(self, key):
        if isinstance(key, slice):
            raise ValueError("ExplicitListView does not support slicing (yet?)")
        super().__delitem__(key)

    def insert(self, key, value):
        if isinstance(key, slice):
            raise ValueError("ExplicitListView does not support slicing (yet?)")
        super().insert(key, value)



class SubView(ExplicitListView):
    pass

class ValidateListView(ListView):
    @abstractmethod
    def validate(self, value):
        """Only values that return true here are included in the view"""
        return True

    def _check_validation(self, value):
        if not self.validate(value):
            raise ValueError(f'{value} is not a valid element of {self}')

    def __setitem__(self, key, value):
        if isinstance(key, slice):
            for v in value:
                self._check_validation(v)
        else:
            self._check_validation(value)
        super().__setitem__(key, value)
        
    def insert(self, key, value):
        self._check_validation(value)
        super().insert(key, value)

    @property
    def indices(self):
        validate = self.validate
        alist = self.alist
        return tuple(i for i in super().indices if validate(alist[i]))

class MutateListView(ListView):
    @abstractmethod
    def mutate(self, value):
        """Values in the view are mutated by this method"""
        return value

    @abstractmethod
    def demutate(self, value):
        """Values modified in the view are mutated by this method"""
        return value

    def __getitem__(self, key):
        out = super().__getitem__(key)
        if isinstance(out, self.SubView):
            # cls.SubView has the same mutate method, so out is already done
            return out
        else:
            return self.mutate(out)

    def __setitem__(self, key, value):
        if isinstance(key, slice):
            values = (self.demutate(v) for v in value)
            super().__setitem__(key, values)
        else:
            super().__setitem__(key, self.demutate(value))

    def insert(self, key, value):
        super().insert(key, self.demutate(value))

    @property
    def indices(self):
        return super().indices


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

    @property
    def indices(self):
        start, stop, step = self.slice.indices(len(self.alist))
        return range(start, stop, step)


class REListView(ValidateListView):
    def __init__(self, alist, regex, *args, **kwargs):
        super().__init__(alist, *args, **kwargs)
        if isinstance(regex, str):
            regex = re.compile(regex)
        self.regex = regex

    def validate(self, value):
        return self.regex.match(value) and super().validate(value)

    def create_subview(self, *args, **kwargs):
        """Create a view to return for slicing

        I added a new required argument 'regex' to __init__, so I need 
        to show REListView how to create a subview of itself. This is
        the method in which that happens."""
        return super().create_subview(regex=self.regex, *args, **kwargs)