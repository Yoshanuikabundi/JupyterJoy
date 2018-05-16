from traitlets import TraitType, HasTraits, TraitError
import six

class MoleculeName(TraitType):
    info_text = 'a non-empty unicode string without whitespace'

    def validate(self, obj, value):
        to_set = None
        if isinstance(value, six.text_type):
            to_set = value
        if isinstance(value, bytes):
            try:
                to_set = value.decode('ascii', 'strict')
            except UnicodeDecodeError:
                msg = "Could not decode {!r} for unicode trait '{}' of {} instance."
                raise TraitError(msg.format(value, self.name, class_of(obj)))
        if to_set and " " not in to_set and "\t" not in to_set:
            return to_set
        self.error(obj, value)

class MoleculeType(HasTraits):
    """A single molecule type"""
    name = MoleculeName()

    def __init__(self, name, **kwargs):
        kwargs['name'] = name
        super().__init__(**kwargs)

    def __repr__(self):
        return f"MoleculeType({self.name})"

    def __str__(self):
        return self.name

    def __format__(self, *args, **kwargs):
        return str.__format__(self.name, *args, **kwargs)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return str(self) == other