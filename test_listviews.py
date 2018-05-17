import pytest
from listviews import *
import re

class TrivialListView(ListView):
    def map_view(self):
        return super().map_view()

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

    a.insert(-1, '80')
    b.insert(-1, '80')
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

def test_ListView_get_set_del():
    a = [1,2,3,4,5,6,7,8,9,10]
    b = [1,2,3,4,5,6,7,8,9,10]
    aview = TrivialListView(a)

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

    # assert a[3:4] == b[3:4] == aview[3:4]

    # aview[3:6] = ['a', 'b', 'c']
    # b[3:6] = ['a', 'b', 'c']
    # assert a == b == aview
    # assert len(a) == len(b) == len(aview)

def test_SliceListView_bigwun():
    a = list(range(1_000_000))
    b = list(range(1_000_000))
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

    print(a)
    a[4] = "hi"
    b[4] = "hi"
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview[-1] = "world"
    b[-1] = "world"
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

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

    print(a)
    a[4] = 132
    b[4] = 264
    assert b == aview
    assert len(a) == len(b) == len(aview)

    aview[-1] = 196
    b[-1] = 196
    assert b == aview
    assert len(a) == len(b) == len(aview)

class ExcitedREListView(REListView, MutateListView):
    def mutate(self, value):
        return super().mutate(value + '!')
    def demutate(self, value):
        if value[-1] != '!':
            raise ValueError("View must take values ending in '!'")
        return super().mutate(value[:-1])

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

