import bqplot as bq
from traitlets import Float, Tuple, observe

@bq.register_mark('bqplotools.VertLine')
class VertLine(bq.Lines):
    position = Float(default_value=0.)

    def __init__(self, position, *args, **kwargs):
        kwargs.setdefault('scales', {'x': bq.LinearScale()})
        kwargs['scales'].update({'y': bq.LinearScale()})
        kwds = {
            'y': [-9999, 9999],
            'x': [position, position],
            'position': position
        }
        if set(kwds).intersection(kwargs):
            raise ValueError("Don't define kwarg x or y")
        kwargs.update(kwds)

        super().__init__(*args, **kwargs)

    @observe('position')
    def update_position(self, change):
        self.x = [self.position, self.position]
