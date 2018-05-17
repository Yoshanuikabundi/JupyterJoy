import pytest
from listviews import *
import re

class TrivialListView(ListView):
    def __filtered__(self):
        return enumerate(self.alist)

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

    print(a)
    a[4] = "hi"
    b[4] = "hi"
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

    aview[-1] = "world"
    b[-1] = "world"
    assert a == b == aview
    assert len(a) == len(b) == len(aview)

def test_ListView_bigwun():
    a = list(range(1000))
    b = list(range(1000))
    aview = TrivialListView(a)

    del aview[500]
    del b[500]

    assert a == b == aview


class TrivialValidateListView(ValidateListView):
    def __filter__(self, value):
        return True

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
