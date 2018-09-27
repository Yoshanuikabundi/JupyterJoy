from copy import copy
import mdtraj as md


def make_martini_ndx(top, f=None):
    """Construct a gromacs ndx file from a martini MDTraj topology"""
    ndx = {}
    for res in top.residues:
        name = res.name.upper()
        indices = [a.index+1 for a in res.atoms]
        try:
            ndx[name] += indices
        except KeyError:
            ndx[name] = indices
    ndx['system'] = [a.index+1 for a in top.atoms]
    ndx['solvent'] = sorted(
        ndx.get('PW', [])
        + ndx.get('W', [])
        + ndx.get('WF', [])
        + ndx.get('NA+', [])
        + ndx.get('CL-', [])
        + ndx.get('NC3+', [])
        + ndx.get('CA+', [])
    )
    ndx['non-solvent'] = [
        a.index+1 for a in top.atoms if a.index+1 not in ndx['solvent']
    ]
    if f:
        for k, v in ndx.items():
            f.write("[ {} ]\n".format(k))
            for n, i in enumerate(v):
                f.write("{: >6} ".format(i))
                if not (n+1) % 5:
                    f.write('\n')
            f.write("\n\n")
    return ndx


class Generator():
    """Easy-to-write generator classes for extensibility

    Just define __gen__ with at least one yield statement in a subclass!"""
    def __init__(self):
        self._gen = self.__gen__()

    def __gen__(self):
        return NotImplemented

    def __next__(self):
        return(next(self._gen))


class BatchNGLViewRenderer(Generator):
    """Helps render a bunch of images from an NGLView view

    When NGLView renders an image, python needs some time (and a new cell
    execution) to catch up to JS. This generator guides you through this
    process and lets you do it by repeatedly executing a single cell.

    Parameters:
    view: an NGLView object we're taking pictures of.
    batch: A dict mapping names to the parameters of `view.render_image()`,
            which is called as `view.render_image(**kwargs)`

    Examples:

    >>> # Render images of frames 3 and 4 of an NGLView 'view'

    >>> # First, create the class instance
    >>> bnr = BatchNGLViewRenderer(view, {
    ...     "pic_of_frame_3": dict(frame=3, factor=1, trim=True),
    ...     "pic_of_frame_42": dict(frame=42, factor=1, trim=True)
    ... })

    >>> # Now, in a different cell so you don't keep making new instances,
    >>> # call the instance. You'll probably want to keep the output too
    >>> images = bnr()

    Now just follow the instructions in the second cell! Output is also
    available as the `out` attribute of the instance, and is returned every
    time the instance is called.
    """

    def __init__(self, view, batch):
        super().__init__()
        self.view = view
        self.batch = dict(batch)
        self.out = {}

    def __gen__(self):
        view = self.view
        batch = self.batch
        out = self.out

        for name, kwargs in batch.items():
            prev_image_data = copy(view._image_data)

            view.frame = kwargs.get("frame", view.frame)
            yield (
                f'Check that the view looks OK for {name}, then execute '
                f'this cell again.'
            )

            while view._image_data == prev_image_data:
                view.render_image(**kwargs)
                yield (
                    'Image (re-)rendered. Execute this cell again in a '
                    'sec. If you see this a lot, try a longer sec!'
                )

            image = view._display_image()
            out[name] = image

    def __call__(self):
        try:
            msg = next(self)
        except StopIteration:
            msg = "All done! You can stop executing now."
        print(msg)
        return self.out


def md_load(*args, **kwargs):
    '''Call mdtraj.load with sane defaults'''
    pdbkwargs = dict(standard_names=False, no_boxchk=True)
    pdbkwargs.update(kwargs)
    try:
        # If it's a pdb, just load it without trying to be clever
        return md.load(*args, **pdbkwargs)
    except TypeError:
        # Not a pdb, so use the defaults
        return md.load(*args, **kwargs)
