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

from JupyterJoy.mdpbuild.mdp import MDP20181
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

    def optimise_rdsq_box(self, target_mindist, tolerance=0.01, init_guess=None):
        init_guess = target_mindist/2 if init_guess is None else init_guess
        ubound = target_mindist
        lbound = 0.0

        while True:
            self.traj.unitcell_vectors = None
            guess = (ubound + lbound) / 2.0
            self.add_rdsq_box(2 * guess + self.calc_maxdist())
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


    def add_rdsq_box(self, d, log=True):
        """Give traj a rhombic dodecahedral box with square in the xy plane and vector length d"""
        d = float(d)
        a = (d, 0.0, 0.0)
        b = (0.0, d, 0.0)
        c = (d/2.0, d/2.0, sqrt(2.0)*d/2.0)
        vecs = np.array(len(traj) * [[a, b, c]])
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
        vecs = np.array(len(traj) * [[a, b, c]])
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

    def call_gmx(self, cmdline=None, cmd=None, stdin='', timeout=None, cwd=None, **kwargs):
        if cmdline and cmd:
            raise ValueError('Specify only one of cmd and cmdline')
        if not (cmdline or cmd):
            raise ValueError('Specify either cmdline or cmd')
        if cmdline and kwargs:
            raise ValueError('Specify arguments in cmdline or as kwargs, not both')

        if cmdline:
            args = ['gmx'] + list(shlex.split(cmdline))
        else:
            args = [('gmx', cmd)] + [self.kwarg_to_arg(k,v) for k,v in kwargs.items()]
            args = [item for sublist in args for item in sublist]

        p = sp.Popen(args, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd)
        stdout_data, stderr_data = p.communicate(bytes(stdin, 'utf-8'), timeout=timeout)
        stdout_data = str(stdout_data, encoding='utf-8')
        stderr_data = str(stderr_data, encoding='utf-8')
        if p.returncode:
            raise ValueError(f'GROMACS returned error code {p.returncode}; stderr was:\n {stderr_data}; contents of cwd are {os.listdir(cwd)}')

        self.log.append(GMXJobLog(
            files=VirtualFiles(cwd),
            stdin=stdin,
            stdout=stdout_data,
            stderr=stderr_data
        ))

        return (p, stdout_data, stderr_data)

    def write(self, path):
        self.ff.write(path=path)
        topfn = f'{path}/{self.name}.top'
        self.topol.write(topfn)
        pdbfn = f'{path}/{self.name}.pdb'
        self.traj.save(pdbfn)
        return WriteInfo(pdb=pdbfn, top=topfn)


    @staticmethod
    def kwarg_to_arg(key, value):
        if value is False or value is "no":
            return (f"-no{key}", )
        elif value is True:
            return (f"-{key}", )
        else:
            return (f"-{key}", value)

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
                p=topinout
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
            mdp = MDP20181(f'{path}/{self.ff.name}.ff/em.mdp')
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

    def load_pdb(self, filename):
        self.traj = md.load_pdb(filename, no_boxchk=True, standard_names=False)
        with TemporaryDirectory() as path:
            pdbin, topin = self.write(path)
            pdbout = f'{self.name}.vis.pdb'
            mdp = MDP20181(f'{path}/{self.ff.name}.ff/em.mdp')
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
