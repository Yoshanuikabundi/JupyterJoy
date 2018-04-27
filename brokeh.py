"""
brokeh.py
Part of the PyMDTools package
Josh Mitchell 2018

Bokeh does not play well with Jupyter widgets. Here are some classes
and functions to help it to be nicer.

It is essential that cds_hack be in the notebook's global namespace. Until
I figure out a better hack, the easiest way to accomplish this is with the following 
import line, which can be used in addition to whatever other imports you want:
from PyMDTools.brokeh import cds_hack
"""
from __future__ import print_function
import ipywidgets as widgets
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, Range1d
from bokeh.plotting import figure, show
from bokeh.plotting.figure import Figure
from bokeh.io import push_notebook, output_notebook
import nglview as nv
import traitlets as tl
import traittypes as tt
from copy import copy, deepcopy
from bokeh.palettes import colorblind 
import numpy as np
from itertools import zip_longest
from mdtraj import Trajectory
from pandas import DataFrame
from collections import Mapping

cds_hack = {'outwidget': widgets.Output()}

def _self_updating_js(ident):
    """Create JS code for a Callback that keeps a 1d CDS selection up to date"""
    return """
        // Define a callback to capture errors on the Python side
        function callback(msg){
            console.log("Python callback returned unexpected message:", msg)
        }
        callbacks = {iopub: {output: callback}};
  
        var selected_str = JSON.stringify(cb_obj.selected['1d'])
        var cmd = "cds_hack['""" + ident + """']['source'].selected['1d'].update(" + selected_str + ")"

        // Execute the command on the Python kernel
        var kernel = IPython.notebook.kernel;
        kernel.execute(cmd, callbacks, {silent : false});
        
        var cmd = "for link in cds_hack['""" + ident + """']['linked']: link.update()"
        kernel.execute(cmd, callbacks, {silent : false});
    """

def nb_column_data_source(*args, **kwargs):
    """Creates a ColumnDataSource who's selected attribute actually updates

    Don't have to use this - Linking a normal CDS to a SeleNGLWidget does all this - 
    but might be useful for future exploits."""
    
    source = ColumnDataSource(*args, **kwargs)
    
    ident = source.ref['id']
    
    # Create/update the global registry of linked datasources
    global cds_hack
    cds_hack[ident] = {'source': source,
                       'linked': []}
    
    # Define the callback
    code = _self_updating_js(ident)
    source.callback = CustomJS(code=code)
    return source

@widgets.register()
class SeleNGLWidget(nv.NGLWidget):
    """Wrapper widget for NGLWidget that can be given selections of frames

    frame_stride: The frame in the selected list corresponds to a frame
    in the trajectory with this factor"""
    selected = tl.Dict()
    frame_stride = tl.Integer()
    
    def __init__(self, *args, orig=None, **kwargs):
        self.original_trajectory = orig
        super().__init__(*args, **kwargs)
        self._current_ind = None
        self._colour = None

        # Make widgets resize sanely:
        self._remote_call("setSize", target='Widget',
                          args=['100%', '100%'])
    
    @tl.observe('selected', type=tl.All)
    def _new_selection(self, change):
        """Update the current frame when a new selection is set"""
        self.update()
        
    def link_to_bokeh_ds(self, datasource):
        """Link this widget to an existing datasource with the appropriate callback"""
        global cds_hack
        ident = datasource.ref['id']
        if ident not in cds_hack:
            cds_hack[ident] = {'source': datasource,
                               'linked': []}
            code = _self_updating_js(ident)
            datasource.callback = CustomJS(code=code)

        cds_hack[ident]['linked'].append(self)
        
        # For some reason all datasources share a selected['1d'] object
        datasource.selected['1d'] = copy(datasource.selected['1d'])
        self.selected = datasource.selected['1d']
        self.update
    
    def update(self):
        """Reset the frame to the first in the selection"""
        self._current_ind = 0
        try:
            self.frame = self.selected['indices'][0] * self.frame_stride
        except IndexError:
            pass
        
    def selected_iter(self):
        """Generator to iterate over selected frames"""
        for frame in self.selected['indices']:
            yield frame * self.frame_stride
            
    def next_selected(self):
        """Display the next frame in the selected list"""
        if self._current_ind is None:
            self._current_ind = 0
        else:
            self._current_ind += 1
            
        try:
            idx = self.selected['indices'][self._current_ind]
        except IndexError:
            if self.selected['indices']:
                self._current_ind = 0
                idx = self.selected['indices'][self._current_ind]
            else:
                idx = self.frame // self.frame_stride
        
        self.frame = idx * self.frame_stride
        return self.frame

def show_mdtraj(mdtraj_trajectory, **kwargs):
    '''Show mdtraj trajectory.

    Examples
    --------
    >>> import nglview as nv
    >>> import mdtraj as md
    >>> t = md.load(nv.datafiles.XTC, top=nv.datafiles.GRO)
    >>> w = nv.show_mdtraj(t)
    >>> w # doctest: +SKIP
    '''
    structure_trajectory = nv.MDTrajTrajectory(mdtraj_trajectory)
    return SeleNGLWidget(structure_trajectory, orig=mdtraj_trajectory, **kwargs)


class FigureTrait(tl.Instance):
    klass = Figure

class BokehWrapper(widgets.Output):
    """A wrapper for a Bokeh figure that lets it act like a widget

    #TODO: It would be super nice if things like selected were traitlets
    in this wrapper, so that they could be easily linked to other widgets.
    Can probably be done very similarly to how we did it for nglview - 
    in which case we probably wouldn't need the special nglview class!"""
    figure = FigureTrait(read_only=True, allow_none=True)
    showing = tl.Bool(default_value=True)

    def __init__(self, figure, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_trait('figure', figure)
        if self.showing:
            self.refresh()

    def show(self):
        """Start showing!"""
        self.showing = True

    @tl.observe('showing')
    def _showing_changed(self, _):
        if self.showing:
            self.refresh()
        else:
            self.clear_output()

    def refresh(self):
        """Refresh where the widget and bokeh need to talk, then redisplay"""
        if not self.showing:
            raise ValueError("Only refresh when widget is showing!")
        self.clear_output()
        self._set_size_from_layout()
        with self:
            _ = show(self.figure, notebook_handle=True)

    def __repr__(self):
        return "BokehWrapper({}, ...)".format(self.figure)

    @tl.observe('layout')
    def _new_layout(self, _):
        """Update the Bokeh height+width when layout changes"""
        print("Layout changed", _)

        self._set_size_from_layout()

        if self.showing:
            self.refresh()

    def _set_size_from_layout(self):
        p = self.figure
        newwidth = int(self.layout.width)
        newheight = int(self.layout.height)

        if p is not None:
            if newwidth and newwidth!=p.plot_width:
                p.plot_width = newwidth
            if newheight and newheight!=p.plot_height:
                p.plot_height = newheight


class NGLViewTrait(tl.Instance):
    klass = SeleNGLWidget

class CDSTrait(tl.Instance):
    klass = ColumnDataSource


class FrameLinkedPlot(widgets.Box):
    ploty = tl.List(tl.Union([tl.Integer(), tl.Unicode()]), default_value=[0])
    plotx = tl.Union([tl.Integer(), tl.Unicode()], default_value=1)
    stride = tl.Integer(default_value=1)
    title = tl.Unicode(default_value='')

    colvars = tl.List(trait=tl.Union([tt.DataFrame(), tt.Array()]), read_only=True)
    bokeh = tl.Instance(klass=BokehWrapper, read_only=True)
    views = tl.List(trait=NGLViewTrait(), read_only=True)
    sources = tl.List(trait=CDSTrait(), read_only=True)

    _reuse_colvars = tl.Bool(read_only=True)

    def_button_layout = widgets.Layout(height='30px')
    def_view_layout   = dict(width='250px', height='300px', flex='0 0 auto')
    def_figure_layout = widgets.Layout(width='400', height='550')
    def_box_layout    = widgets.Layout(display='flex',
                                    flex_flow='column wrap',
                                    height='650px',
                                    width='100%',
                                    align_items='center')

    def __init__(self, 
        colvars, 
        mdtraj_trajectories, 
        plotx_label=None, 
        ploty_label=None, 
        ploty=None, 
        legendlocation='top_right',
        **kwargs):
        super().__init__(**kwargs)

        # colvars and mdtraj_trajectories can be specified either as a list or a single item
        if isinstance(mdtraj_trajectories, Trajectory):
            mdtrajs = [mdtraj_trajectories]
        else:
            mdtrajs = list(mdtraj_trajectories)

        if isinstance(colvars, DataFrame) or isinstance(colvars, Mapping):
            cvs = [colvars]
        else:
            cvs = list(colvars)

        if len(mdtrajs) != 1 and len(cvs) != len(mdtrajs) and len(ploty) != len(mdtrajs):
            print(len(mdtrajs), len(cvs), len(ploty))
            raise ValueError("Should have either 1 mdtraj trajectory, or one for every colvar")

        # ploty wasn't given to super, so it's currently the default per the traitlet
        # So if it was user-specified, we need to make sure its a list and then
        # throw it in
        if isinstance(ploty, str):
            print("ploty is a str")
            self.ploty = [ploty] * len(cvs)
        elif ploty is not None:
            try:
                self.ploty = list(iter(ploty))
            except TypeError:
                self.ploty = [ploty] * len(cvs)
            else:
                if len(cvs) == 1: cvs = cvs * len(self.ploty)
        if len(cvs) != len(self.ploty):
            raise ValueError("Couldn't broadcast colvars to plotys: {}, {}".format(len(cvs), len(ploty)))

        self.set_trait('colvars', cvs)
        
        self.layout = copy(self.def_box_layout)

        plotx = self.plotx
        ploty = self.ploty

        TOOLS="crosshair,pan,zoom_in,zoom_out,box_zoom,box_select,undo,reset,"

        hover = HoverTool(tooltips=[
            ("x", "@x (@plotx_label)"),
            ("y", "@y (@ploty_label)"),
            ("time", "@time ns"),
            ("run", "@run")
        ])

        if plotx_label is None:
            plotx_label = "Collective variable {}".format(plotx+1)
        if plotx_label is None:
            ploty_label = "Collective variable {}".format(ploty+1)

        p = figure(
            title=self.title,
            x_axis_label=plotx_label,
            y_axis_label=ploty_label,
            tools=TOOLS)
        p.add_tools(hover)
        figure_layout = copy(self.def_figure_layout)
        bokeh = BokehWrapper(p, layout=figure_layout, showing=False)
        self.set_trait('bokeh', bokeh)

        sources, views = self._init_plots_views(mdtrajs, 
                                                plotx_label, 
                                                ploty_label)

        p.legend.click_policy = "hide"
        p.legend.location = legendlocation
        if len(cvs) == 1 or len(cvs) == len(mdtrajs):
            p.legend.visible = False

        self.bokeh.show()

        self.set_trait('views', views)
        self.set_trait('sources', sources)

        button = widgets.Button(
            description='Next selected frame',
            disabled=False,
            button_style='', # 'success', 'info', 'warning', 'danger' or ''
            tooltip='Set NGLView frame to the next selected point',
            layout=copy(self.def_button_layout)
        )

        button.on_click(self._on_button_clicked)

        self.children =tuple([self.bokeh, button] + self.views)

    def _on_button_clicked(self, b):
        for view in self.views:
            view.next_selected()

    def _init_plots_views(self, mdtraj_trajectories, plotx_label, ploty_label):
        colvars = self.colvars
        stride = self.stride
        plotx = self.plotx
        plotys = self.ploty
        p = self.figure

        n = len(colvars)
        if n >= 3:
            palette = colorblind['Colorblind'][n]
        else:
            palette = colorblind['Colorblind'][3]
        palette = palette[::-1]

        view_layout = widgets.Layout(**self.def_view_layout)
        if len(mdtraj_trajectories) == 1:
            view_layout.width='500px'
            view_layout.height='600px'

        sources = []
        views = []

        for n,(colvar, traj, ploty) in enumerate(zip_longest(colvars, mdtraj_trajectories, plotys)):
            if traj is not None:
                working_traj = traj

            times = working_traj.time[::stride]
            
            if len(colvar) != len(working_traj):
                raise ValueError("Colvar and trajectory should have same number of frames")
            
            if isinstance(colvar, np.ndarray):
                x = colvar[::stride, plotx]
                y = colvar[::stride, ploty]
            else:
                x = colvar[plotx][::stride]
                y = colvar[ploty][::stride]

            if isinstance(ploty, str):
                this_ploty_label = ploty
            elif ploty_label is not None:
                this_ploty_label = ploty_label
            else:
                this_ploty_label = str(ploty)
            ploty_label_list = [this_ploty_label] * len(y)

            source = ColumnDataSource(data={
                'run': [n]*len(x),
                'plotx_label': [plotx_label]*len(x),
                'ploty_label': ploty_label_list,
                'time': times/1000,
                'x': x,
                'y': y,
                'alphas': [(t)/(times[-1]) for t in times]
            })
            sources.append(source)
            
            colour = palette[n-1]
            
            if traj is not None:
                view = show_mdtraj(traj, gui=False)
                view._colour = colour
                if len(traj.top.select('protein')):
                    view.clear_representations()
                    view.add_cartoon(selection='polymer', color=colour)
                view.frame_stride = stride
                view.layout = view_layout
                view._set_sync_camera()
                views.append(view)

            view.link_to_bokeh_ds(source)
            
            p.scatter(x='x', y='y', 
                      source=source, 
                      color=colour, 
                      fill_alpha='alphas',
                      legend=this_ploty_label)

        return sources, views

    @property
    def figure(self):
        return self.bokeh.figure
    
    @property
    def mdtraj_trajectories(self):
        return (view.original_trajectory for view in self.views)

    @property
    def view(self):
        if len(self.views) == 1:
            return self.views[0]
        else:
            raise ValueError("FrameLinkedPlot has multiple views")

    @property
    def source(self):
        if len(self.views) == 1:
            return self.sources[0]
        else:
            raise ValueError("FrameLinkedPlot has multiple sources")
    
def plot_vs_time(toplot, df, traj, timeaxis='Time', **kwargs):
    """Plot df's columns with headings in toplot against time"""
    try:
        df_time = df[timeaxis]
    except KeyError:
        # Dict or dataframe or similar
        try:
            df[timeaxis] = traj.time
        except AttributeError:
            try:
                times = [t.time for t in traj.values()]
            except AttributeError:
                times = [t.time for t in traj]
            for t in times:
                if not np.isclose(times[0], t).all():
                    raise ValueError("All trajs must have same times!")
            df[timeaxis] = times[0]
    except (IndexError, TypeError):
        # Array or list or similar
        try: 
            data = df[toplot]
        except (IndexError, TypeError):
            data = df
        df = DataFrame({timeaxis: traj.time, toplot: data})
    else:
        if not np.isclose(df_time, traj.time).all():
            raise ValueError("Inconsistent times between df and traj")

    if 'ploty_label' not in kwargs:
        kwargs['ploty_label'] = str(toplot)
    if 'title' not in kwargs:
        kwargs['title'] = kwargs['ploty_label'] + " over time"
    
    if isinstance(toplot, str):
        plotys = toplot
        trajs = traj
    else:
        try:
            plotys = list(toplot)
        except TypeError:
            plotys = toplot
            trajs = traj
        else:
            try:
                trajs = [traj[k] for k in plotys]
            except:
                trajs = traj


    
    linkwidg = FrameLinkedPlot(colvars=df, 
                               mdtraj_trajectories=trajs, 
                               plotx_label='Time (ps)', 
                               plotx=timeaxis,
                               ploty=plotys,
                               **kwargs)

    return linkwidg


output_notebook()