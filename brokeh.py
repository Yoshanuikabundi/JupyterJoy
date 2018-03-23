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
from bokeh.io import push_notebook
import nglview as nv
import traitlets as tl
import traittypes as tt
from copy import copy
from bokeh.palettes import colorblind 
import numpy as np

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





def linked_bokeh_ngl(figure, nglviews, 
    button_layout=None, nglview_layout=None, 
    figure_layout=None, box_layout=None):
    if button_layout is None:
        button_layout  = widgets.Layout(height='30px')
    if nglview_layout is None:
        nglview_layout = widgets.Layout(width='250px', height='300px', flex='0 0 auto')
    if figure_layout is None:
        figure_layout  = widgets.Layout(width='400', height='550')
    if box_layout is None:
        box_layout     = widgets.Layout(display='flex',
                                        flex_flow='column wrap',
                                        height='650px',
                                        width='100%',
                                        align_items='center')

    out = BokehWrapper(figure, layout=figure_layout)

    button = widgets.Button(
        description='Next selected frame',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Set NGLView frame to the next selected point',
        layout=button_layout
    )
    for nglview in nglviews:
        nglview.layout = nglview_layout
        nglview._remote_call("setSize", target='Widget',
                             args=["100%", "100%"])


    def on_button_clicked(b):
        for nglview in nglviews:
            nglview.next_selected()

    button.on_click(on_button_clicked)

    # nvs = widgets.VBox([button] + nglviews)

    # return widgets.HBox([out, nvs])
    return widgets.Box([out, button] + nglviews, layout=box_layout)


class NGLViewTrait(tl.Instance):
    klass = SeleNGLWidget

class CDSTrait(tl.Instance):
    klass = ColumnDataSource


class FrameLinkedPlot(widgets.Box):
    plotx = tl.Integer(default_value=0)
    ploty = tl.Integer(default_value=1)
    stride = tl.Integer(default_value=1)
    title = tl.Unicode(default_value='')

    reduced = tt.Array(read_only=True)
    bokeh = tl.Instance(klass=BokehWrapper, read_only=True)
    views = tl.List(trait=NGLViewTrait(), read_only=True)
    sources = tl.List(trait=CDSTrait(), read_only=True)

    def_button_layout = widgets.Layout(height='30px')
    def_view_layout   = widgets.Layout(width='250px', height='300px', flex='0 0 auto')
    def_figure_layout = widgets.Layout(width='400', height='550')
    def_box_layout    = widgets.Layout(display='flex',
                                    flex_flow='column wrap',
                                    height='650px',
                                    width='100%',
                                    align_items='center')



    def __init__(self, reduced, mdtraj_trajectories, plotx_label=None, ploty_label=None, **kwargs):
        super().__init__(**kwargs)
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
            plotx_label = "Reduced coordinate {}".format(plotx+1)
        if plotx_label is None:
            ploty_label = "Reduced coordinate {}".format(ploty+1)

        p = figure(
            title='Cartesian coordinate PCA',
            x_axis_label=plotx_label,
            y_axis_label=ploty_label,
            tools=TOOLS)
        p.add_tools(hover)
        figure_layout = copy(self.def_figure_layout)
        bokeh = BokehWrapper(p, layout=figure_layout, showing=False)

        self.set_trait('reduced', reduced)
        self.set_trait('bokeh', bokeh)

        sources, views = self._init_plots_views(mdtraj_trajectories, 
                                                plotx_label, 
                                                ploty_label)
        self.set_trait('views', views)
        self.set_trait('sources', sources)
        self.bokeh.show()

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
        mdtraj_trajectories
        reduced = self.reduced
        stride = self.stride
        plotx = self.plotx
        ploty = self.ploty
        p = self.figure
        palette = colorblind['Colorblind'][len(mdtraj_trajectories)]
        view_layout = copy(self.def_view_layout)

        sources = []
        views = []
        for n,(xyz, traj) in enumerate(zip(reduced, mdtraj_trajectories)):
            times = traj.time[::stride]
            
            x = xyz[::stride, plotx]
            y = xyz[::stride, ploty]
            
            source = ColumnDataSource(data={
                'run': [n]*len(x),
                'plotx_label': [plotx_label]*len(x),
                'ploty_label': [ploty_label]*len(y),
                'time': times/1000,
                'x': x,
                'y': y,
                'alphas': [(t)/(times[-1]) for t in times]
            })
            sources.append(source)
            
            colour = palette[n-1]
            
            view = show_mdtraj(traj, gui=False)
            view.link_to_bokeh_ds(source)
            view.add_cartoon(color=colour)
            view.frame_stride = stride
            view.layout = view_layout
            view._set_sync_camera()
            views.append(view)
            
            p.scatter(x='x', y='y', source=source, color=colour, fill_alpha='alphas')

        return sources, views

    @property
    def figure(self):
        return self.bokeh.figure

    @property
    def mdtraj_trajectories(self):
        return [view.original_trajectory for view in self.views]
    