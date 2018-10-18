import bqplot as bq
from functools import partial

class EZFigure(bq.Figure):
    def __init__(self, *args, **kwargs):
        """Create scales and axes if they're not defined

        If axes is defined, see if you can identify unique scales for either axis.
        If you can, that's the scale for that axis.
        If axes is undefined, make new axes, using scale_x and scale_y if they're defined.

        If something's defined, that thing goes. We just try to fill in as many gaps as possible."""
        if 'axes' in kwargs:
            axes_c = set()
            scales_x = set()
            scales_y = set()
            for a in kwargs['axes']:
                if isinstance(a, bq.ColorAxis):
                    axes_c.add(a)
                elif a.orientation == 'vertical':
                    scales_y.add(a.scale)
                else:
                    scales_x.add(a.scale)
            scale_x = scales_x.pop() if len(scales_x) == 1 else bq.LinearScale()
            scale_y = scales_y.pop() if len(scales_y) == 1 else bq.LinearScale()

            kwargs.setdefault('scale_x', scale_x)
            kwargs.setdefault('scale_y', scale_y)
        else:
            scale_x = kwargs.setdefault('scale_x', bq.LinearScale())
            scale_y = kwargs.setdefault('scale_y', bq.LinearScale())
            axis_x = bq.Axis(
                scale = scale_x,
                label = kwargs.pop('label_x', '')
            )
            axis_y = bq.Axis(
                scale=scale_y,
                label = kwargs.pop('label_y', ''),
                orientation = 'vertical'
            )
            kwargs['axes'] = [axis_x, axis_y]

        super().__init__(*args, **kwargs)


    def add_mark(self, klass, *args, **kwargs):
        """Add a mark by instantiating klass with self's scales."""
        kwargs.setdefault('scales', {'x': self.scale_x, 'y': self.scale_y})
        mark = klass(*args, **kwargs)
        self.marks = self.marks + [mark]
        return mark

    @property
    def _available_marks(self):
        return {k.lower().split('.')[-1]:v for k, v in bq.Mark.mark_types.items()}


    def __getattribute__(self, name):
        marks = super().__getattribute__('_available_marks')
        try:
            klass = marks[name]
        except KeyError:
            return super().__getattribute__(name)
        else:
            return partial(self.add_mark, klass)

    def __dir__(self):
        dir_dict = set(super().__dir__())
        marks = self._available_marks
        dir_dict.update(marks)
        return dir_dict

    @staticmethod
    def _is_x_axis(axis):
        is_color = isinstance(axis, bq.ColorAxis)
        is_vertical = axis.orientation == 'vertical'
        return not is_color and not is_vertical

    @staticmethod
    def _is_y_axis(axis):
        is_color = isinstance(axis, bq.ColorAxis)
        is_vertical = axis.orientation == 'vertical'
        return not is_color and is_vertical

    @property
    def axes_x(self):
        return [a for a in self.axes if self._is_x_axis(a)]
    @property
    def axis_x(self):
        axes_x = self.axes_x
        if len(axes_x) == 1:
            return axes_x.pop()
        else:
            raise ValueError('No unique x axis')

    @property
    def axes_y(self):
        return [a for a in self.axes if self._is_y_axis(a)]
    @property
    def axis_y(self):
        axes_y = self.axes_y
        if len(axes_y) == 1:
            return axes_y.pop()
        else:
            raise ValueError('No unique y axis')


# for name, klass in bq.Mark.mark_types.items():
#     name = name.lower().split('.')[-1]
#     method = partialmethod(EZFigure.add_mark, klass)
#     setattr(EZFigure, name, method)
