import pytest
from listviews import *
import re

class TrivialListView(ListView):
    @property
    def indices(self):
        return super().indices

def test_ListView_append_insert():
    a = [1,2,3,4,5,6]
    b = [1,2,3,4,5,6]
    aview = TrivialListView(a)

    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a.append(7)
    b.append(7)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview.append(8)
    b.append(8)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview.insert(5, 4.5)
    b.insert(5, 4.5)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a.insert(5, 3.5)
    b.insert(5, 3.5)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview.insert(-1, [25])
    b.insert(-1, [25])
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview.insert(0, 68)
    b.insert(0, 68)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a.insert(-1, '80')
    b.insert(-1, '80')
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a.insert(-8, 38)
    b.insert(-8, 38)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a = []
    b = []
    aview = TrivialListView(a)
    aview.append(3)
    b.append(3)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a = []
    b = []
    aview = TrivialListView(a)
    aview.insert(0,3)
    b.insert(0,3)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a = []
    b = []
    aview = TrivialListView(a)
    aview.insert(1,3)
    b.insert(1,3)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a = []
    b = []
    aview = TrivialListView(a)
    aview.insert(-1,3)
    b.insert(-1,3)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)


def test_ListView_get_set_del():
    a = [0,1,2,3,4,5,6,7,8,9,10]
    b = [0,1,2,3,4,5,6,7,8,9,10]
    aview = TrivialListView(a)

    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    assert a[3] == 3
    assert aview[3] == 3
    assert b[3] == 3
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    assert a[-1] == 10
    assert aview[-1] == 10
    assert b[-1] == 10
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    del a[3]
    del b[3]
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    del aview[4]
    del b[4]
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a[4] = "hi"
    b[4] = "hi"
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview[-1] = "world"
    b[-1] = "world"
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

<<<<<<< HEAD
def slices():
    numbers = [0, 1, 2, 3, 9, 10, 11, 12, 20]
    for start in numbers:
        for stop in numbers:
            for step in numbers:
                combos = [
                    ( start,  stop,  step),
                    ( start,  stop, -step),
                    ( start,  stop,  None),
                    ( start, -stop,  step),
                    ( start, -stop, -step),
                    ( start, -stop,  None),
                    ( start,  None,  step),
                    ( start,  None, -step),
                    ( start,  None,  None),
                    (-start,  stop,  step),
                    (-start,  stop, -step),
                    (-start,  stop,  None),
                    (-start, -stop,  step),
                    (-start, -stop, -step),
                    (-start, -stop,  None),
                    (-start,  None,  step),
                    (-start,  None, -step),
                    (-start,  None,  None),
                    (  None,  stop,  step),
                    (  None,  stop, -step),
                    (  None,  stop,  None),
                    (  None, -stop,  step),
                    (  None, -stop, -step),
                    (  None, -stop,  None),
                    (  None,  None,  step),
                    (  None,  None, -step),
                    (  None,  None,  None)
                ]
                yield from combos

def test_ListView_slicing():
    a = [0,1,2,3,4,5,6,7,8,9,10,11]
    b = [0,1,2,3,4,5,6,7,8,9,10,11]
    aview = TrivialListView(a)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    for start, stop, step in slices():
        a = [0,1,2,3,4,5,6,7,8,9,10,11]
        b = [0,1,2,3,4,5,6,7,8,9,10,11]
        aview = TrivialListView(a)
        try:
            bsl = b[start:stop:step]
        except Exception as e:
            with pytest.raises(type(e)):
                avsl = aview[start:stop:step]
        else:
            try:
                asl = a[start:stop:step]
                avsl = aview[start:stop:step]
                assert asl == bsl == avsl
            except AssertionError as e:
                print(f"seq[{start}:{stop}:{step}]")
                raise e

def test_ListView_delslicing():
    for start, stop, step in slices():
        a = [0,1,2,3,4,5,6,7,8,9,10,11]
        b = [0,1,2,3,4,5,6,7,8,9,10,11]
        aview = TrivialListView(a)
        try:
            del b[start:stop:step]
        except Exception as e:
            with pytest.raises(type(e)):
                del aview[start:stop:step]
        else:
            try:
                del aview[start:stop:step]
                assert b == a
                assert a == b == aview
                assert len(a) == len(b) == len(aview)
            except Exception as e:
                print(f"del seq[{start}:{stop}:{step}]")
                raise e

def test_ListView_setslicing():
    vals = [
        '',
        'a',
        'ab',
        'abc',
        'abcd',
        'abcde',
        'abcdefghijklmnopqrstuvwxyz'
    ]
    for start, stop, step in slices():
        for value in vals:
            a = [0,1,2,3,4,5,6,7,8,9,10,11]
            b = [0,1,2,3,4,5,6,7,8,9,10,11]
            aview = TrivialListView(a)
            try:
                b[start:stop:step] = value
            except Exception as e:
                with pytest.raises(type(e)):
                    aview[start:stop:step] = value
            else:
                try:
                    aview[start:stop:step] = value
                    assert b == a
                    assert a == b == aview
                    assert len(a) == len(b) == len(aview)
                except Exception as e:
                    print(f"seq[{start}:{stop}:{step}] = {repr(value)}")
                    print(f"len(idx): {len(aview.indices[start:stop:step])}")
                    raise e

def test_SliceListView_bigwun():
    a = list(range(1_100_000))
    b = list(range(1_100_000))
=======
def test_SliceListView_bigwun():
    a = list(range(1_000_000))
    b = list(range(1_000_000))
>>>>>>> 501c13d80c57e1b10caf0be3fa9d88f977a0e4a5
    start, stop, step = None, None, 100_000
    aview = SliceListView(a, start, stop, step)

    idx = 5
    del aview[idx]
    del b[idx*step]

    assert a == b 
    assert a[::step] == b[::step] == aview
    assert len(a[::step]) == len(b[::step]) == len(aview)

    idx = 7
    del a[idx*step]
    del b[idx*step]

    assert a == b 
    assert a[::step] == b[::step] == aview
    assert len(a[::step]) == len(b[::step]) == len(aview)


class TrivialValidateListView(ValidateListView):
    def validate(self, value):
        return super().validate(value)

def test_ValidateListView_append_insert():
    a = [1,2,3,4,5,6]
    b = [1,2,3,4,5,6]
    aview = TrivialValidateListView(a)

    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a.append(7)
    b.append(7)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview.append(8)
    b.append(8)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview.insert(5, 4.5)
    b.insert(5, 4.5)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a.insert(5, 3.5)
    b.insert(5, 3.5)
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview.insert(-1, [25])
    b.insert(-1, [25])
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a.insert(-1, '80')
    b.insert(-1, '80')
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

def test_ValidateListView_get_set_del():
    a = [1,2,3,4,5,6,7,8,9,10]
    b = [1,2,3,4,5,6,7,8,9,10]
    aview = TrivialValidateListView(a)

    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    assert a[3] == 4
    assert aview[3] == 4
    assert b[3] == 4
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    assert a[-1] == 10
    assert aview[-1] == 10
    assert b[-1] == 10
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    del a[3]
    del b[3]
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    del aview[4]
    del b[4]
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    a[4] = "hi"
    b[4] = "hi"
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview[-1] = "world"
    b[-1] = "world"
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

class OddValidateListView(ValidateListView):
    def validate(self, value):
        return super().validate(value) and value % 2 

def test_OddValidateListView_slicing():
    a = [0,1,2,3,4,5,6,7,8,9,10,11]
    b = [1,3,5,7,9,11]
    aview = OddValidateListView(a)
    assert b == aview
    assert len(b) == len(aview)

    for start, stop, step in slices():
        a = [0,1,2,3,4,5,6,7,8,9,10,11]
        b = [1,3,5,7,9,11]
        aview = OddValidateListView(a)
        try:
            bsl = b[start:stop:step]
        except Exception as e:
            with pytest.raises(type(e)):
                avsl = aview[start:stop:step]
        else:
            try:
                avsl = aview[start:stop:step]
                assert bsl == avsl
            except AssertionError as e:
                print(f"seq[{start}:{stop}:{step}]")
                raise e

def test_OddValidateListView_delslicing():
    for start, stop, step in slices():
        a = [0,1,2,3,4,5,6,7,8,9,10,11]
        b = [1,3,5,7,9,11]
        aview = OddValidateListView(a)
        try:
            del b[start:stop:step]
        except Exception as e:
            with pytest.raises(type(e)):
                del aview[start:stop:step]
        else:
            try:
                del aview[start:stop:step]
                assert b == aview
                assert len(b) == len(aview)
            except Exception as e:
                print("a:", a)
                print("aview.indices:", aview.indices)
                print(f"del seq[{start}:{stop}:{step}]")
                raise e

def test_OddValidateListView_setslicing():
    vals = [
        [],
        [901],
        [901,903],
        [901,903,905],
        [901,903,905,907],
        [901,903,905,907],
        [901,903,905,907,909,911,913,915,917,919],
    ]
    for start, stop, step in slices():
        for value in vals:
            a = [0,1,2,3,4,5,6,7,8,9,10,11]
            b = [1,3,5,7,9,11]
            aview = OddValidateListView(a)
            print(a, aview)
            try:
                b[start:stop:step] = value
            except Exception as e:
                with pytest.raises(type(e)):
                    aview[start:stop:step] = value
            else:
                try:
                    aview[start:stop:step] = value
                    assert b == aview
                    assert len(b) == len(aview)
                    for i in a:
                        if i % 2:
                            assert b.pop(0) == i
                except Exception as e:
                    print(f"seq[{start}:{stop}:{step}] = {repr(value)}")
                    print(f"len(idx): {len(aview.indices[start:stop:step])}")
                    raise e

def test_REListView():
    l = ["hello","gday","yo","sup","howzit"]
    re_pattern = re.compile(r'[hy]')
    hyview = REListView(l, re_pattern)
    gyview = REListView(l, r'[gy]')

    assert hyview == ["hello","yo","howzit"]
    assert gyview == ["gday","yo"]

    with pytest.raises(ValueError):
        gyview.append("hi")
    with pytest.raises(ValueError):
        gyview[1] = ("hi")
    gyview.append("yodelehihoo")
    assert l == ["hello","gday","yo","yodelehihoo","sup","howzit"]
    assert hyview == ["hello","yo","yodelehihoo","howzit"]
    assert gyview == ["gday","yo","yodelehihoo"]
    assert l != hyview != gyview

    l.insert(-1, "hi")
    assert l == ["hello","gday","yo","yodelehihoo","sup","hi","howzit"]
    assert hyview == ["hello","yo","yodelehihoo","hi","howzit"]
    assert gyview == ["gday","yo","yodelehihoo"]
    assert l != hyview != gyview

    hyview[3] = "heya"
    assert l == ["hello","gday","yo","yodelehihoo","sup","heya","howzit"]
    assert hyview == ["hello","yo","yodelehihoo","heya","howzit"]
    assert gyview == ["gday","yo","yodelehihoo"]
    assert l != hyview != gyview

class DoubleListView(MutateListView):
    def mutate(self, value):
        return super().mutate(value * 2)
    def demutate(self, value):
        return super().mutate(value / 2)

def test_MutateListView_append_insert():
    a = [1,2,3,4,5,6]
    b = [2,4,6,8,10,12]
    aview = DoubleListView(a)

    assert b == aview
    assert len(a) == len(b) == len(aview)

    a.append(7)
    b.append(14)
    assert b == aview
    assert len(a) == len(b) == len(aview)

    aview.append(8)
    b.append(8)
    assert b == aview
    assert len(a) == len(b) == len(aview)

    aview.insert(5, 9)
    b.insert(5, 9)
    assert b == aview
    assert len(a) == len(b) == len(aview)

    a.insert(5, 3.5)
    b.insert(5, 7)
    assert b == aview
    assert len(a) == len(b) == len(aview)

    aview.insert(-1, 50)
    b.insert(-1, 50)
    assert b == aview
    assert len(a) == len(b) == len(aview)

    a.insert(-1, 80)
    b.insert(-1, 160)
    assert b == aview
    assert len(a) == len(b) == len(aview)

def test_MutateListView_get_set_del():
    a = [1,2,3,4,5,6,7,8,9,10]
    b = [2,4,6,8,10,12,14,16,18,20]
    aview = DoubleListView(a)

    assert b == aview
    assert len(a) == len(b) == len(aview)

    assert a[3] == 4
    assert aview[3] == 8
    assert b[3] == 8
    assert b == aview
    assert len(a) == len(b) == len(aview)

    assert a[-1] == 10
    assert aview[-1] == 20
    assert b[-1] == 20
    assert b == aview
    assert len(a) == len(b) == len(aview)

    del a[3]
    del b[3]
    assert b == aview
    assert len(a) == len(b) == len(aview)

    del aview[4]
    del b[4]
    assert b == aview
    assert len(a) == len(b) == len(aview)

    a[4] = 132
    b[4] = 264
    assert b == aview
    assert len(a) == len(b) == len(aview)

    aview[-1] = 196
    b[-1] = 196
    assert b == aview
    assert len(a) == len(b) == len(aview)

class ExcitedREListView(REListView, MutateListView):
    def validate(self, value):
        if not isinstance(value, str):
            raise TypeError("Doesn't validate ExcitedREListView")
        return isinstance(value, str) and super().validate(value)
    def mutate(self, value):
        return super().mutate(value + '!')
    def demutate(self, value):
        if value[-1] != '!':
            raise ValueError("View must take values ending in '!'")
        return super().mutate(value[:-1])

class NotAString:
    def __bool__(self):
        return True
    def __add__(self, other):
        return self
    def __getitem__(self, key):
        return 'g!'

def test_ListView_multiple_inheritance():
    l = ["hello","gday","yo","sup","howzit"]
    re_pattern = re.compile(r'[hy]')
    hyview = ExcitedREListView(l, re_pattern)
    gyview = ExcitedREListView(l, r'[gy]')

    assert hyview == ["hello!","yo!","howzit!"]
    assert gyview == ["gday!","yo!"]

    with pytest.raises(ValueError):
        gyview.append("hi")
    with pytest.raises(ValueError):
        gyview[1] = ("hi")
    with pytest.raises(ValueError):
        gyview.append("yodelehihoo")
    assert NotAString() + '!'
    with pytest.raises(TypeError):
        gyview.append(NotAString())
    gyview.append("yodelehihoo!")
    assert l == ["hello","gday","yo","yodelehihoo","sup","howzit"]
    assert hyview == ["hello!","yo!","yodelehihoo!","howzit!"]
    assert gyview == ["gday!","yo!","yodelehihoo!"]
    assert l != hyview != gyview

    l.insert(-1, "hi")
    assert l == ["hello","gday","yo","yodelehihoo","sup","hi","howzit"]
    assert hyview == ["hello!","yo!","yodelehihoo!","hi!","howzit!"]
    assert gyview == ["gday!","yo!","yodelehihoo!"]
    assert l != hyview != gyview

    hyview[3] = "heya!"
    assert l == ["hello","gday","yo","yodelehihoo","sup","heya","howzit"]
    assert hyview == ["hello!","yo!","yodelehihoo!","heya!","howzit!"]
    assert gyview == ["gday!","yo!","yodelehihoo!"]
    assert l != hyview != gyview

    assert hyview[0:2] == ["hello!","yo!"]
    assert hyview[3:3] == []
    hyview[3:3] = ["yaaaassss!"]
    assert hyview == ["hello!", "yo!", "yodelehihoo!", "yaaaassss!","heya!","howzit!"]
    possibilities = [
        ["hello","gday","yo","yodelehihoo","sup","yaaaassss","heya","howzit"],
        ["hello","gday","yo","yodelehihoo","yaaaassss","sup","heya","howzit"]
    ]
    assert l in possibilities
    with pytest.raises(ValueError):
        hyview[3:7] = ["hi"]
