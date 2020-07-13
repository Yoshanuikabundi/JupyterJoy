from copy import copy
import mdtraj as md
import subprocess
import os
from urllib.request import urlretrieve
from tempfile import TemporaryDirectory
import numpy as np
from itertools import tee, count
import mdtraj.core.element as elem
from re import sub
import contextlib
import pathlib
from pathlib import Path

aa_tlc = [
    'ala', 'arg', 'asn', 'asp', 'cys', 'glu', 'gln', 'gly', 'his',
    'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp',
    'tyr', 'val'
]

aa_olc = [
    'a', 'r', 'n', 'd', 'c', 'e', 'q', 'g', 'h',
    'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w',
    'y', 'v'
]

aa_olc2tlc = {olc: tlc for olc, tlc in zip(aa_olc, aa_tlc)}
aa_tlc2olc = {tlc: olc for olc, tlc in zip(aa_olc, aa_tlc)}


def rand_rotation_matrices(num=1, deflection=1.0):
    """
    Creates an array of random rotation matrices.

    num: number of rotation matrices to generate
    deflection: the magnitude of the rotation. For 0, no rotation; for 1,
    completely randomrotation. Small deflection => small perturbation.

    Adapted from
    http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
    """
    randnums = np.random.uniform(size=(3, num))

    theta, phi, z = randnums

    theta = theta * 2.0*deflection*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0*deflection  # For magnitude of pole deflection.

    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.

    r = np.sqrt(z)
    Vx, Vy, Vz = V = np.stack([
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
    ])

    st = np.sin(theta)
    ct = np.cos(theta)

    row1 = np.stack([ct, st, np.zeros(num)], axis=-1)
    row2 = np.stack([-st, ct, np.zeros(num)], axis=-1)
    row3 = np.stack([np.zeros(num), np.zeros(num), np.ones(num)], axis=-1)

    R = np.stack([row1, row2, row3], axis=-1)

    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    M = (V.T[..., None] * V.T[:, None, :] - np.eye(3)) @ R
    return M

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
        key, line_is_proper, value = line.partition("=")
        if line_is_proper:
            env_out[key]= value
        else:
            print("WARNING: Line", line, "has no value")
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

def write_dict_to_ndx(ndxdict, f):
    for k, v in ndxdict.items():
        f.write("[ {} ]\n".format(k))
        for n, i in enumerate(v):
            f.write("{: >6} ".format(i))
            if not (n+1) % 5:
                f.write('\n')
        f.write("\n\n")

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
    for chain in top.chains:
        name = "chain" + chr(chain.index + 65)
        indices = [a.index+1 for a in chain.atoms]
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

    protein_idcs = []
    for aa  in aa_tlc:
        protein_idcs += ndx.get(aa.upper(), [])
    ndx['protein'] = sorted(protein_idcs)

    ndx['non-solvent'] = sorted(in_b_but_not_a(ndx['solvent'], ndx['system']))
    ndx['non-protein'] = sorted(in_b_but_not_a(ndx['protein'], ndx['system']))

    if f:
        write_dict_to_ndx(ndx, f)
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

    if isinstance(filename_or_filenames, pathlib.PurePath):
        filename_or_filenames = str(filename_or_filenames)

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

class Unreachable(Exception):
    pass

def in_b_but_not_a(a, b):
    a = iter(sorted(a))
    b = iter(sorted(b))
    try:
        thisa = next(a)
        thisb = next(b)
    except StopIteration:
        yield from b
    else:
        while True:
            step_a = False
            step_b = False
            if thisb < thisa:
                yield thisb
                step_b = True
            elif thisb == thisa:
                step_a = True
                step_b = True
            elif thisb > thisa:
                step_a = True
            else:
                raise Unreachable()

            if step_a:
                try:
                    thisa = next(a)
                except StopIteration:
                    if thisb != thisa:
                        yield thisb
                    yield from b
                    break
            if step_b:
                try:
                    thisb = next(b)
                except StopIteration:
                    break


def window_gen(iterable, n):
    "s, n -> (s0, s1, ...), (s1, s2, ...), ..., (sn, sn+1, ...)"
    for _ in range(n-1):
        out, iterable = tee(iterable)
        next(iterable, None)
        yield out
    yield iterable

def window(iterable, n):
    "s, n -> (s0, s1, ..., sn), (s1, s2, ..., sn+1), (s2, s3, ..., sn+2), ..."
    return zip(*window_gen(iterable, n))

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def load_gro(filename):
    filename = os.path.abspath(filename)
    traj = md.load(filename)
    velxyz = np.zeros_like(traj.xyz)
    coord_indices = None
    top = md.Topology()
    chain = top.add_chain()
    residue = None
    with open(filename) as f:
        lines = iter(f)
        for model in count():
            try:
                title = next(lines)
            except StopIteration:
                assert model == len(traj)
                break
            _, _, time = title.partition('t=')
            n_atoms = int(next(lines))
            assert n_atoms == traj.n_atoms
            for i in range(n_atoms):
                line = next(lines)
                if model == 0:
                    resnum = int(line[0:5])
                    resname = line[5:10].strip()
                    atomname = line[10:15].strip()
                    atomnum = int(line[15:20])

                    if residue is None or resnum != residue.resSeq:
                        residue = top.add_residue(resname, chain, resSeq=resnum)

                    if len(atomname) > 1:
                        elem_symbol = atomname[0] + sub('[A-Z0-9]','',atomname[1:])
                    else:
                        elem_symbol = atomname

                    try:
                        element = elem.get_by_symbol(elem_symbol)
                    except KeyError:
                        element = elem.virtual
                    top.add_atom(atomname, element=element, residue=residue, serial=atomnum)


                if coord_indices is None:
                    decs = (i for i, v in enumerate(line[20:], start=20) if v == '.')
                    decidist = abs(next(decs) - next(decs))
                    coord_indices = list(pairwise(range(20, 20 + decidist * 7, decidist)))

                x, y, z, *vel = (float(line[a:b]) for a, b in coord_indices if b <= len(line))
                assert np.isclose(traj.xyz[model, i, :], (x, y, z)).all()
                if vel:
                    velxyz[model, i, :] = vel
            boxvecs = next(lines)
    traj.top = top
    if np.all(velxyz == 0.0):
        velxyz = None
    return traj, velxyz


def save_gro(filename, traj, velxyz=None, force_overwrite=True, precision=3):
    if velxyz is not None and traj.xyz.shape != velxyz.shape:
        raise ValueError('traj.xyz and velxyz have different shapes')
    filename = os.path.abspath(filename)
    traj.save_gro(filename, force_overwrite=force_overwrite, precision=precision)

    if velxyz is not None:
        fstr = f'{{: >{precision+5}.{precision+1}f}}' * 3
        coord_indices = list(pairwise(range(20, 20 + (precision + 5) * 7, (precision + 5))))
        with open(filename, 'r') as f:
            lines = f.readlines()
        for model, modelvels in enumerate(velxyz):
            for atom, atomvels in enumerate(modelvels):
                velstr = fstr.format(*atomvels)
                linenum = model * (velxyz.shape[1] + 2) + 2 + atom
                lines[linenum] = lines[linenum][:-1] + velstr + lines[linenum][-1]
        with open(filename, 'w') as f:
            for line in lines:
                f.write(line)



@contextlib.contextmanager
def modified_environ(*remove, **update):
    """
    Temporarily updates the ``os.environ`` dictionary in-place.

    The ``os.environ`` dictionary is updated in-place so that the modification
    is sure to work in all situations.

    :param remove: Environment variables to remove.
    :param update: Dictionary of environment variables and values to add/update.
    """
    env = os.environ
    update = update or {}
    remove = remove or []

    # List of environment variables being updated or removed.
    stomped = (set(update.keys()) | set(remove)) & set(env.keys())
    # Environment variables and values to restore on exit.
    update_after = {k: env[k] for k in stomped}
    # Environment variables and values to remove on exit.
    remove_after = frozenset(k for k in update if k not in env)


    try:
        env.update({k: str(v) for k, v in update.items()})
        [env.pop(k, None) for k in remove]
        yield env
    finally:
        env.update(update_after)
        [env.pop(k) for k in remove_after]
