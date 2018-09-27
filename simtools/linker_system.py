import os
from collections import namedtuple
import nglview as nv
import mdtraj as md
from copy import copy,deepcopy
import numpy as np
from scipy.spatial.distance import cdist, pdist
from panedr import edr_to_df
from tempfile import TemporaryDirectory
import tarfile
import scipy
import time

import sys

from JupyterJoy.mdpbuild.mdp import MDP20183
from JupyterJoy.tophat import Topology
from JupyterJoy.simtools.prep_remd import calc_temps
from JupyterJoy.bqploteins import TrajPlotTime
from JupyterJoy.funfuncs import md_load

import JupyterJoy.pbash

import pandas as pd

import PeptideBuilder as pb
import Bio.PDB

import subprocess as sp
import shlex

from math import sqrt

import io

def make_extended_peptide(seq):
    geo = pb.Geometry.geometry
    seq_iter = iter(seq)
    structure = pb.initialize_res(geo(next(seq_iter)))
    for aa in seq_iter:
        structure = pb.add_residue(structure, geo(aa))
    return structure

def unwrap(iterator):
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

def replace_children(entity, new_children):
    for child in entity.get_list():
        child.detach_parent()
    entity.child_list = []
    entity.child_dict = {}
    for child in new_children:
        entity.add(child)

def rename_atom(atom, name):
    atom.id = name
    atom.name = name
    atom.fullname = f'{name: >3}'

def convert_to_ace(res):
    res.resname = 'ACE'
    ch3, c, o = res.child_dict['CA'], res.child_dict['C'], res.child_dict['O']
    rename_atom(ch3, 'CH3')
    replace_children(res, [ch3, c, o])

def convert_to_nme(res):
    res.resname = 'NME'
    n, ch3 = res.child_dict['N'], res.child_dict['CA']
    rename_atom(ch3, 'CH3')
    replace_children(res, [n, ch3])

def make_extended_capped_peptide(seq):
    """Make an extended and capped model of seq, and save it"""
    structure = make_extended_peptide(f'G{seq}G')
    rs = unwrap(structure.get_chains()).child_list
    convert_to_ace(rs[0])
    convert_to_nme(rs[-1])
    return structure

class PDBTrajectoryVirtFile(md.formats.PDBTrajectoryFile):
    def __init__(self, file, standard_names=False):
        self._open = True
        self._file = file
        self._topology = None
        self._positions = None
        self._mode = 'r'
        self._last_topology = None
        self._standard_names = standard_names

        self._read_models()

def load_virt_pdb(
    file,
    stride=None,
    atom_indices=None,
    frame=None,
    no_boxchk=False,
    standard_names=False
):
    atom_indices = md.utils.cast_indices(atom_indices)
    with PDBTrajectoryVirtFile(file, standard_names=standard_names) as f:
        atom_slice = slice(None) if atom_indices is None else atom_indices
        if frame is not None:
            coords = f.positions[[frame], atom_slice, :]
        else:
            coords = f.positions[::stride, atom_slice, :]
        assert coords.ndim == 3, 'internal shape error'
        n_frames = len(coords)

        topology = f.topology
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        if f.unitcell_angles is not None and f.unitcell_lengths is not None:
            unitcell_lengths = np.array([f.unitcell_lengths] * n_frames)
            unitcell_angles = np.array([f.unitcell_angles] * n_frames)
        else:
            unitcell_lengths = None
            unitcell_angles = None

        md.utils.in_units_of(coords, f.distance_unit, md.Trajectory._distance_unit, inplace=True)
        md.utils.in_units_of(unitcell_lengths, f.distance_unit, md.Trajectory._distance_unit, inplace=True)

    time = np.arange(len(coords))
    if frame is not None:
        time *= frame
    elif stride is not None:
        time *= stride

    traj = md.Trajectory(xyz=coords, time=time, topology=topology,
                      unitcell_lengths=unitcell_lengths,
                      unitcell_angles=unitcell_angles)

    if not no_boxchk and traj.unitcell_lengths is not None:
        # Only one CRYST1 record is allowed, so only do this check for the first
        # frame. Some RCSB PDB files do not *really* have a unit cell, but still
        # have a CRYST1 record with a dummy definition. These boxes are usually
        # tiny (e.g., 1 A^3), so check that the particle density in the unit
        # cell is not absurdly high. Standard water density is ~55 M, which
        # yields a particle density ~100 atoms per cubic nm. It should be safe
        # to say that no particle density should exceed 10x that.
        particle_density = traj.top.n_atoms / traj.unitcell_volumes[0]
        if particle_density > 1000:
            warnings.warn('Unlikely unit cell vectors detected in PDB file likely '
                          'resulting from a dummy CRYST1 record. Discarding unit '
                          'cell vectors.')
            traj._unitcell_lengths = traj._unitcell_angles = None

    return traj

def bio_struct_to_mdtraj(structure):
    writer = Bio.PDB.PDBIO()
    writer.set_structure(structure)
    vfile = io.StringIO()
    writer.save(vfile)
    vfile.seek(0)
    traj = load_virt_pdb(vfile)
    vfile.close()
    return traj



def mindist(traj, frame=0):
    xyz = traj.xyz[frame]
    vecs = traj.unitcell_vectors[frame]
    image_mins = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                if (0,0,0) == (i,j,k): continue
                image_offset = np.sum(
                    vecs * np.asarray([(float(i),), (float(j),), (float(k),)]),
                    axis=0
                )
                image_xyz = xyz + image_offset
                dists = cdist(xyz, image_xyz)
                image_mindist = dists.min()
                image_mins.append(image_mindist)

    return min(image_mins)

def maxdist(traj, frame=0):
    '''Maximum internal distance'''
    return pdist(traj.xyz[frame]).max()

class VirtualFiles():
    def __init__(self, path):
        path = path[:-1] if path[-1] == '/' else path
        self._rawfile = io.BytesIO()
        tgz = tarfile.open(
            fileobj=self._rawfile,
            mode='w:gz'
        )
        tgz.add(
            path,
            filter=self.filter
        )
        tgz.close()
        self._rawfile.seek(0)
        self._tgz = tarfile.open(
            fileobj=self._rawfile,
            mode='r'
        )
        self.name = os.path.basename(path)

    def write(self, path='.'):
        self._tgz.extractall(f'{path}')

    def filter(self, tarinfo):
        return tarinfo


class GMXForceField(VirtualFiles):
    def __init__(self, ffdir):
        ffdir = ffdir[:-1] if ffdir[-1] == '/' else ffdir
        if ffdir[-3:] != '.ff':
            raise ValueError('ffdir should end in ".ff"')
        self.name = os.path.basename(ffdir[:-3])
        super().__init__(ffdir)
        self.name = os.path.basename(ffdir[:-3])

    def write(self, path='.'):
        self._tgz.extractall(f'{path}')
        return f'{path}/{self.name}.ff'

    def filter(self, tarinfo):
        path = tarinfo.name
        dirname, basename = os.path.split(path)
        if tarinfo.isdir():
            tarinfo.name = f'{self.name}.ff'
        else:
            tarinfo.name = f'{self.name}.ff/{basename}'
        return tarinfo

class TopolWithItps(Topology):
    def __init__(self, fname, **kwargs):
        dirname = os.path.dirname(fname)
        self._rawfile = io.BytesIO()
        tgz = tarfile.open(
            fileobj=self._rawfile,
            mode='w:gz'
        )
        tgz.add(dirname, filter=self.itpfilter)
        tgz.close()
        self._rawfile.seek(0)
        self._tgz = tarfile.open(
            fileobj=self._rawfile,
            mode='r'
        )

        super().__init__(fname=fname, **kwargs)

    def write(self, f):
        try:
            fname = f.name
        except AttributeError:
            fname = f
        path = os.path.dirname(fname)
        self._tgz.extractall(path)
        super().write(f)

    @staticmethod
    def itpfilter(tarinfo):
        path = tarinfo.name
        if path.endswith('.itp'):
            tarinfo.name = os.path.basename(path)
            return tarinfo
        elif tarinfo.isdir() and not path.endswith('.ff'):
            tarinfo.name = '.'
            return tarinfo
        else:
            return None


def align_to(a, b):
    """Return a rotation matrix R such that Ra is parallel to b

    from https://math.stackexchange.com/a/476311"""
    norm = np.linalg.norm
    # Convert to unit vectors
    a = np.asarray(a) / norm(a)
    b = np.asarray(b) / norm(b)

    # Define the identity matrix
    I = np.array([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]])

    # If the unit vectors are (anti)parallel, return the (negated)
    # identity matrix
    if (a == b).all():
        return I

    if (a == -b).all():
        return -I

    # Otherwise, use the magic from above

    v = np.cross(a,b)
    s = norm(v)
    c = np.dot(a,b)

    v_cross = np.array([[    0, -v[2],  v[1]],
                        [ v[2],     0, -v[0]],
                        [-v[1],  v[0],     0]])

    R = I + v_cross + ((v_cross @ v_cross) * (1/(1+c)))

    return R

WriteInfo = namedtuple('WriteInfo', ['pdb', 'top'])

GMXJobLog = namedtuple('GMXJob', ['files', 'stdin', 'stdout', 'stderr'])
BoxLog = namedtuple('BoxLog', ['boxtype', 'd', 'vectors'])

class GMXLinkerSystem():
    def __init__(self, sequence, name=None, wd=None):
        self.seq = sequence

        self.name = sequence if name is None else name

        structure = make_extended_capped_peptide(seq)
        self._traj = bio_struct_to_mdtraj(structure)

        self.ff = GMXForceField('/store/joshmitchell/linkers/charmm22star_kcx.ff')
        self.watermodel = 'tips3p'
        boxbuffer = 1.2

        self.log = []

        with TemporaryDirectory() as td:
            self.ff.write(td)
            pdbin = 'extended.pdb'
            pdbout = 'extended_gmx.pdb'
            topout = 'topol.top'
            self.traj.save(f'{td}/{pdbin}')
            _, stdout_data, stderr_data = self.call_gmx(
                cmd='pdb2gmx',
                stdin='3\n5\n',
                cwd=td,
                f=pdbin,
                o=pdbout,
                p=topout,
                water=self.watermodel,
                ter=True,
                ignh=True,
                ff=self.ff.name
            )
            print(stderr_data)
            print(stdout_data)
            self.topol = TopolWithItps(f'{td}/{topout}')
            self.traj = md.load_pdb(f'{td}/{pdbout}', no_boxchk=True, standard_names=False)

    def orient_protein(self, target_vec):
        if 'SOL' in self.topol.molecules:
            raise ValueError("Can't orient solvated system")
        pdists = scipy.spatial.distance.squareform(pdist(self.coords))
        idx_a, idx_b = np.unravel_index(np.argmax(pdists), pdists.shape)
        vec_a, vec_b = self.coords[idx_a], self.coords[idx_b]
        long_axis = vec_a - vec_b
        long_axis = long_axis / np.linalg.norm(long_axis)
        rot_mat = align_to(long_axis, target_vec)
        new_coords = [rot_mat @ c for c in self.coords]
        self.traj.xyz = np.array([new_coords])

    def optimise_rdsq_box(self, target_mindist, tolerance=0.01):
        # Orient the protein's long axis with the x-axis - this is the worst-case orientation
        self.orient_protein([1, 0, 0])

        # Optimise the box
        ubound = target_mindist
        lbound = 0.0

        while True:
            self.traj.unitcell_vectors = None
            guess = (ubound + lbound) / 2.0
            d = 2 * guess + self.calc_maxdist()
            self.add_rdsq_box(d)
            this_mindist = self.calc_mindist()
            print(ubound, guess, lbound, this_mindist)
            if abs(this_mindist - target_mindist) <= tolerance:
                break
            elif this_mindist < target_mindist:
                lbound = guess
            elif this_mindist > target_mindist:
                ubound = guess
            else:
                raise ValueError('stuck')
            del self.log[-1]

        print(f'A buffer of {guess} nm gives a box with d={d}. The minimum PI distance is {this_mindist}.')

        # And reorient to the z-axis - best-case orientation
        self.orient_protein([0, 0, 1])

        print(f'After reorientation, the PI distance is {self.calc_mindist()}.')


    def add_rdsq_box(self, d, log=True):
        """Give traj a rhombic dodecahedral box with square in the xy plane and vector length d"""
        d = float(d)
        a = (d, 0.0, 0.0)
        b = (0.0, d, 0.0)
        c = (d/2.0, d/2.0, sqrt(2.0)*d/2.0)
        vecs = np.array(len(self.traj) * [[a, b, c]])
        if self.traj.unitcell_vectors is not None:
            raise ValueError(f'self.traj already has a unit cell: {traj.unitcell_vectors}')
        self.traj.unitcell_vectors = vecs
        self.log.append(BoxLog(
            'rhombic dodecahedron (xy-square)',
            d,
            vecs
        ))
        return self.traj

    def add_rdhex_box(self, d):
        """Give traj a rhombic dodecahedral box with hexagon in the xy plane and vector length d"""
        d = float(d)
        a = (d, 0.0, 0.0)
        b = (d/2.0, sqrt(3.0)*d/2.0, 0.0)
        c = (d/2.0, sqrt(3.0)*d/6.0, sqrt(6.0)*d/3.0)
        vecs = np.array(len(self.traj) * [[a, b, c]])
        if self.traj.unitcell_vectors:
            raise ValueError(f'self.traj already has a unit cell')
        self.traj.unitcell_vectors = vecs
        self.log.append(BoxLog(
            'rhombic dodecahedron (xy-hexagon)',
            d,
            vecs
        ))
        return self.traj

    def calc_maxdist(self):
        return maxdist(self.traj)

    def calc_mindist(self):
        xyz = self.coords
        vecs = self.box
        image_mins = []
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                for k in (-1, 0, 1):
                    if (0,0,0) == (i,j,k): continue
                    image_offset = np.sum(
                        vecs * np.asarray([(float(i),), (float(j),), (float(k),)]),
                        axis=0
                    )
                    image_xyz = xyz + image_offset
                    dists = cdist(xyz, image_xyz)
                    image_mindist = dists.min()
                    image_mins.append(image_mindist)
        return min(image_mins)

    @property
    def coords(self):
        return self.traj.xyz[0]

    @property
    def box(self):
        return self.traj.unitcell_vectors[0]

    @property
    def traj(self):
        return self._traj

    @traj.setter
    def traj(self, traj):
        if not isinstance(traj, md.Trajectory):
            raise ValueError('traj should be an mdtraj trajectory')
        if len(traj) != 1:
            raise ValueError('traj should have only one frame')
        self._traj = traj

    def call_gmx(self, cmdline=None, cmd=None, stdin='', timeout=None, cwd=None, mpiranks=None, **kwargs):
        cmd = cmd.strip()
        mdrun_synonyms = ['mdrun', 'mdrun_mpi']
        if cmdline and cmd:
            raise ValueError('Specify only one of cmd and cmdline')
        if not (cmdline or cmd):
            raise ValueError('Specify either cmdline or cmd')
        if cmdline and kwargs:
            raise ValueError('Specify arguments in cmdline or as kwargs, not both')
        if cmdline and mpiranks:
            raise ValueError('Specify mpi in cmdline or as mpiranks, not both')
        if mpiranks and cmd not in mdrun_synonyms:
            raise ValueError('mpi is only compatible with cmd=mdrun_mpi')

        if cmdline:
            args = ['gmx'] + list(shlex.split(cmdline))
        else:
            args = [('gmx', cmd)] + [self.kwarg_to_arg(k,v) for k,v in kwargs.items()]
            args = [item for sublist in args for item in sublist]

        cmd_is_mdrun = args[1] in mdrun_synonyms

        if mpiranks:
            args = ['mpirun', '-np', str(mpiranks), 'mdrun_mpi'] + args[2:]

        path = cwd if cwd is not None else '.'
        if 'deffnm' in kwargs:
            outname = f'{kwargs["deffnm"]}'
        else:
            outname = f'gmx_{args[1]}'
        outname = f'{path}/{outname}'
        i=0
        while True:
            eout = f'{outname}.e.{i:02}'
            oout = f'{outname}.o.{i:02}'
            if (
                os.path.exists(eout)
                or os.path.exists(oout)
            ):
                i += 1
            else:
                break
        with open(eout, mode='x+') as e, open(oout, mode='x+') as o:
            p = sp.Popen(args, stdin=sp.PIPE, stdout=o, stderr=e, cwd=cwd, encoding='utf-8')
            if stdin:
                p.communicate(stdin, timeout=timeout)
            elif cmd_is_mdrun:
                with open(eout, mode='rb') as f:
                    while p.poll() is None:
                        s = f.read(1000)
                        sys.stdout.write(s)
                        if not s:
                            time.sleep(1)
                    s = f.read()
                    sys.stdout.write(s)

            else:
                p.wait()
            e.seek(0)
            o.seek(0)
            stdout_data = o.read()
            stderr_data = e.read()
        if p.returncode:
            raise ValueError(f'GROMACS returned error code {p.returncode}; stderr was:\n {stderr_data}; contents of cwd are {os.listdir(cwd)}')

        self.log.append(GMXJobLog(
            files=VirtualFiles(path),
            stdin=stdin,
            stdout=stdout_data,
            stderr=stderr_data
        ))

        return (p, stdout_data, stderr_data)

    def write(self, path):
        topfn = self.write_top(path)
        pdbfn = self.write_pdb(path)
        return WriteInfo(pdb=pdbfn, top=topfn)

    def write_top(self, path):
        self.ff.write(path=path)
        topfn = f'{path}/{self.name}.top'
        self.topol.write(topfn)
        return topfn

    def write_pdb(self, path):
        pdbfn = f'{path}/{self.name}.pdb'
        self.traj.save(pdbfn)
        return pdbfn



    @staticmethod
    def kwarg_to_arg(key, value):
        if value is False or value is "no":
            return [f"-no{key}"]
        elif value is True:
            return [f"-{key}"]
        elif isinstance(value, str):
            return [f"-{key}", value]
        else:
            try:
                return [f"-{key}"] + list(value)
            except TypeError:
                return [f"-{key}", str(value)]

    def solvate(self):
        with TemporaryDirectory() as path:
            pdbin, topinout = self.write(path)
            pdbout = 'extended_sol.pdb'
            _, stdout_data, stderr_data = self.call_gmx(
                cmd='solvate',
                stdin='',
                cwd=path,
                cp=pdbin,
                o=pdbout,
                p=topinout,
                cs=f'{self.ff.name}.ff/tips3p884.gro'
            )
            print(stderr_data)
            print(stdout_data)
            self.topol = TopolWithItps(topinout)
            self.load_pdb(f'{path}/{pdbout}')

    def salt(self, conc=0.15):
        with TemporaryDirectory() as path:
            pdbin, topinout = self.write(path)
            pdbout = 'extended_ion.pdb'
            tpr_tmp = 'genion.tpr'

            n_wat = unwrap(self.topol.molecules['SOL']).copies
            mol_water = 55.5
            desired_salt_conc = conc
            n_salt = int(round(n_wat * desired_salt_conc / mol_water))
            mdp = MDP20183(f'{path}/{self.ff.name}.ff/em.mdp')
            mdp.write(f'{path}/genion.mdp')

            _, stdout_data, stderr_data = self.call_gmx(
                cmd='grompp',
                stdin='',
                cwd=path,
                f='genion.mdp',
                c=pdbin,
                p=topinout,
                o=tpr_tmp
            )
            _, stdout_data, stderr_data = self.call_gmx(
                cmd='genion',
                stdin='SOL',
                cwd=path,
                s=tpr_tmp,
                o=pdbout,
                p=topinout,
                np=str(n_salt),
                neutral=True
            )
            print(stderr_data)
            print(stdout_data)
            self.topol = TopolWithItps(topinout)
            self.load_pdb(f'{path}/{pdbout}')



    def em(self):
        with TemporaryDirectory() as path:
            pdbin, topinout = self.write(path)
            pdbout = 'em.pdb'
            tpr_tmp = 'em.tpr'

            mdp = MDP20183(f'{path}/{self.ff.name}.ff/em.mdp')
            mdp.write(f'{path}/em.mdp')

            print(f'Energy minimising with .MDP file:\n{mdp.write()}')

            _, stdout_data, stderr_data = self.call_gmx(
                cmd='grompp',
                stdin='',
                cwd=path,
                f='em.mdp',
                c=pdbin,
                p=topinout,
                o=tpr_tmp
            )
            print(stderr_data)
            print(stdout_data)
            _, stdout_data, stderr_data = self.call_gmx(
                cmd='mdrun',
                stdin='',
                cwd=path,
                deffnm='em',
                c=pdbout,
                v=True
            )
            self.topol = TopolWithItps(topinout)
            self.load_pdb(f'{path}/{pdbout}')

    def trajvis(self, filename, cwd=None):
        pathname, _, extension = filename.rpartition('.')
        trajout_tmp = f'{pathname}.tmp.{extension}'
        trajout = f'{pathname}.vis.{extension}'
        pdbout = f'{pathname}.vis.first.pdb'

        with TemporaryDirectory() as path:
            pdbin, topin = self.write(path)
            mdp_tmp=f'{path}/trjconv.mdp'
            tpr_tmp=f'{path}/trjconv.tpr'
            mdp = MDP20183(f'{path}/{self.ff.name}.ff/em.mdp')
            mdp.write(mdp_tmp)
            self.call_gmx(
                cmd='grompp',
                stdin='',
                cwd=cwd,
                f=mdp_tmp,
                c=pdbin,
                p=topin,
                o=tpr_tmp
            )
            self.call_gmx(
                cmd='trjconv',
                stdin='Protein\nSystem\n',
                cwd=cwd,
                f=filename,
                s=tpr_tmp,
                o=trajout_tmp,
                pbc='mol',
                ur='compact',
                center=True
            )
            if extension == 'xtc':
                self.call_gmx(
                    cmd='trjconv',
                    stdin='System\n',
                    cwd=cwd,
                    f=trajout_tmp,
                    s=tpr_tmp,
                    o=pdbout,
                    dump=0
                )
                self.call_gmx(
                    cmd='trjconv',
                    stdin='C-alpha\nsystem\n',
                    cwd=cwd,
                    f=trajout_tmp,
                    s=tpr_tmp,
                    o=trajout,
                    fit='trans'
                )
                os.remove(trajout_tmp)
            else:
                os.rename(trajout_tmp, trajout)




    def load_pdb(self, filename):
        self.traj = md.load_pdb(filename, no_boxchk=True, standard_names=False)
        with TemporaryDirectory() as path:
            pdbin, topin = self.write(path)
            pdbout = f'{self.name}.vis.pdb'
            mdp = MDP20183(f'{path}/{self.ff.name}.ff/em.mdp')
            mdp.write(f'{path}/trjconv.mdp')
            tpr_tmp='trjconv.tpr'
            _, stdout_data, stderr_data = self.call_gmx(
                cmd='grompp',
                stdin='',
                cwd=path,
                f='trjconv.mdp',
                c=pdbin,
                p=topin,
                o=tpr_tmp
            )
            _, stdout_data, stderr_data = self.call_gmx(
                cmd='trjconv',
                stdin='Protein\nSystem\n',
                cwd=path,
                f=pdbin,
                s=tpr_tmp,
                o=pdbout,
                pbc='mol',
                ur='compact',
                center=True
            )
            self.traj = md.load_pdb(f'{path}/{pdbout}', no_boxchk=True, standard_names=False)

    def run_sim(self, mdp, deffnm, cwd='.'):
        pdbin, topinout = self.write(f'{cwd}')
        tprinout = f'{cwd}/{deffnm}.tpr'
        mdpinout = f'{cwd}/{deffnm}.mdp'

        mdp.write(mdpinout)
        print(f'Running simulation with .MDP file:\n{mdp.write()}')

        _, stdout_data, stderr_data = self.call_gmx(
            cmd='grompp',
            stdin='',
            cwd=cwd,
            f=mdpinout,
            c=pdbin,
            p=topinout,
            o=tprinout
        )
        print(stderr_data)
        print(stdout_data)
        _, stdout_data, stderr_data = self.call_gmx(
            cmd='mdrun',
            stdin='',
            cwd=cwd,
            deffnm=deffnm,
            v=True
        )
        self.trajvis(f'{cwd}/{deffnm}.xtc')

    def load_xtc(self, filename):
        return md.load_xtc(filename, top=self.traj.top)
