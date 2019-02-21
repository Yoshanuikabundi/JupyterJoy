from copy import copy
import mdtraj as md
import subprocess
import os
from urllib.request import urlretrieve
from tempfile import TemporaryDirectory


def get_env_from_sh(scriptfn, shell='sh', inheritenv=False):
    """Get a dict containing the environment after sourcing scriptfn"""
    if not os.path.isfile(scriptfn):
        raise ValueError(f'{scriptfn}: file does not exist')

    if inheritenv:
         command = f'{shell} -c ". {scriptfn} && env"'
    else:
         command = f'env -i {shell} -c ". {scriptfn} && env"'

    env_out = {}
    for line in subprocess.getoutput(command).split("\n"):
        key, value = line.split("=")
        env_out[key]= value
    return env_out

def source(scriptfn, shell='sh'):
    """Modify the environment as if we'd sourced scriptfn"""
    os.environ.update(get_env_from_sh(scriptfn, shell=shell, inheritenv=True))


class WorkingDirectory():
    """Context manager that cds in and out on enter/exit"""
    def __init__(self, target_dir):
        os.makedirs(target_dir, exist_ok=True)
        self.target_dir = target_dir
        self.init_dir = os.getcwd()

    def __enter__(self):
        os.chdir(self.target_dir)
        return os.getcwd()

    def __exit__(self, *args):
        os.chdir(self.init_dir)

def addpath(line):
    os.environ['PATH'] = line + ':' + os.environ['PATH']
    return os.environ['PATH']

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


def md_load(
    filename_or_filenames,
    discard_overlapping_frames=True,
    selection=None,
    chunks=0,
    **kwargs
):
    '''Call mdtraj.(iter)load with more convenience

    If you want speed/efficiency, do it manually - this is designed to
    be flexible and convenient. For instance, if you ask for a selection
    without a topology, the file will be loaded twice!'''
    if chunks:
        loader = md.iterload
        kwargs["chunks"] = chunks
    else:
        loader = md.load

    kwargs["discard_overlapping_frames"] = discard_overlapping_frames
    pdbkwargs = dict(standard_names=False, no_boxchk=True)

    # Load the topology with the same logic
    top = kwargs.pop("top", None)
    if top:
        if isinstance(top, md.Trajectory):
            top = top.top
        elif not isinstance(top, md.Topology):
            top = md_load(top).top
        kwargs["top"] = top

    # Allow only a subset to be loaded by selection
    atom_indices = kwargs.pop("atom_indices", None)
    if selection and atom_indices:
        raise ValueError("Specify selection or indices, not both")
    elif selection and not top:
        atom_indices = md_load(
            filename_or_filenames,
            stride=999999
        ).top.select(selection)
    elif selection:
        atom_indices = top.select(selection)
    kwargs["atom_indices"] = atom_indices

    pdbkwargs.update(kwargs)
    try:
        # If it's a pdb, just load it without trying to be clever
        return loader(filename_or_filenames, **pdbkwargs)
    except TypeError:
        # Not a pdb, so use the defaults
        return loader(filename_or_filenames, **kwargs)


def unwrap(iterator):
    '''Get the sole element of iterator or raise a ValueError'''
    i = iter(iterator)
    try:
        out = next(i)
    except StopIteration:
        raise ValueError("i has no elements")
    try:
        next(i)
    except StopIteration:
        return out
    raise ValueError("i has more than one element")


def get_pdbs(*pdbs, **kwargs):
    with TemporaryDirectory() as td:
        for pdb in pdbs:
            url = f"http://files.rcsb.org/view/{pdb}.pdb"
            path = f"{td}/{pdb}.pdb"
            urlretrieve(url, path)
            traj = md.load_pdb(path, **kwargs)
            yield traj

def pdbget(pdb, **kwargs):
    return unwrap(get_pdbs(pdb))

