import mdtraj as md
import bqplot as bq
import ipywidgets as w
from traitlets import Integer, Float, link, observe
from traittypes import Array
import nglview as nv
from .ezfigure import EZFigure
import numpy as np


class TrajPlot(w.Box):
    """Class responsible for compositing a traj with a figure"""
    frame = Integer()
    time = Float()
    t_scale = Float(1000)

    def __init__(self, traj, figure, *args, **kwargs):
        super().__init__(*args, **kwargs)

        view = nv.show_mdtraj(traj, **kwargs)

        self.traj = traj
        self.view = view
        self.figure = figure
        self.children = [figure, view]

        link((self, 'frame'), (self.view, 'frame'))

        self.view._remote_call(
            "setSize",
            target='Widget',
            args=['100%', '100%']
        )

        self.view.layout = w.Layout(max_width="1000px")
        self.figure.layout = w.Layout()
        self.layout = w.Layout(
            display='flex',
            flex_flow='row wrap',
            align_items='flex-start'
        )

    @observe('frame')
    def frame2time(self, changes):
        frame = changes['new']
        dt = self.traj.timestep / self.t_scale
        time = frame * dt
        self.time = float(time)

    @observe('time')
    def time2frame(self, changes):
        time = changes['new']
        dt = self.traj.timestep / self.t_scale
        frame = time / dt
        self.frame = int(frame)


class TrajPlotTime(TrajPlot):
    selected = Array(None, allow_none=True)
    stride = Integer(1)

    def __init__(self, traj, data_y, *args, stride=1, **kwargs):
        if not isinstance(traj, md.Trajectory):
            raise ValueError('traj must be an MDTraj Trajectory')

        unit = {
            1e-3: 'fs',
            1e0: 'ps',
            1e3: 'ns',
            1e6: 'us',
            1e9: 'ms',
            1e12: 's'
        }[self.t_scale]
        kwargs.setdefault('label_x', f'Time ({unit})')
        figure = EZFigure(**kwargs)

        data_y = np.asarray(data_y)
        
        if data_y.ndim == 1:
            data_y = np.expand_dims(data_y, 0)

        if not (data_y.ndim == 2 and data_y.shape[1] == len(traj)):
            raise ValueError('traj and data_y should have same lengths')

        for y_line in data_y:
            scatter = figure.scatter(
                x=traj.time[::stride] / self.t_scale,
                y=y_line[::stride]
            )

        line = figure.vertline(self.frame)

        selector = bq.interacts.IndexSelector(
            line_width=0,
            scale=figure.scale_x,
            marks=[scatter]
        )
        figure.interaction = selector

        super().__init__(traj, figure, *args, **kwargs)

        self.stride = stride
        link((self, 'selected'), (scatter, 'selected'))
        link((self, 'time'), (line, 'position'))

    @observe('selected')
    def sele2time(self, change):
        new = change['new']
        if new is None:
            return
        else:
            new = new[0]
        x = self.traj.time[::self.stride][new] / self.t_scale
        self.time = float(x)


# tp = TrajPlotTime(traj, props[prop], stride=1000, label_y=prop)
# tp
