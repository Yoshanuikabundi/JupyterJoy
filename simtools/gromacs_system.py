from collections import namedtuple
import subprocess as sp
import shlex
from tempfile import TemporaryDirectory
import tarfile
from time import sleep
import pickle
from datetime import datetime
import sys
import io
import os
from itertools import cycle
from shutil import copyfile

from pathlib import Path

import mdtraj as md

import JupyterJoy
from JupyterJoy.mdpbuild.mdp import MDP20183
from JupyterJoy.tophat import Topology
from JupyterJoy.simtools.prep_remd import calc_temps
from JupyterJoy.simtools.rest2 import temper_from_file
from JupyterJoy.funfuncs import unwrap, make_martini_ndx, write_dict_to_ndx, load_gro, save_gro

import PeptideBuilder as pb
import Bio.PDB

import numpy as np
from numpy import sqrt
import scipy
from scipy.spatial.distance import cdist, pdist

from panedr import edr_to_df
from urllib.request import urlretrieve

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
        path = Path(path)
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
        self.name = path.name

    def write(self, path='.'):
        self._tgz.extractall(f'{path}')

    def filter(self, tarinfo):
        return tarinfo


class GMXForceField(VirtualFiles):
    def __init__(self, ffdir):
        ffdir = Path(ffdir)
        if ffdir.suffix != '.ff':
            raise ValueError('ffdir should end in ".ff"')
        self.name = ffdir.stem
        super().__init__(ffdir)
        self.name = ffdir.stem

    def write(self, path='.'):
        path = Path(path).absolute()
        self._tgz.extractall(path)
        return path / f'{self.name}.ff'

    def filter(self, tarinfo):
        path = Path(tarinfo.name)
        dirname, basename = os.path.split(path)
        if tarinfo.isdir():
            tarinfo.name = f'{self.name}.ff'
        else:
            tarinfo.name = f'{self.name}.ff/{path.name}'
        return tarinfo

class TopolWithItps(Topology):
    def __init__(self, fname, **kwargs):
        dirname = Path(fname).absolute().parent
        self._rawfile = io.BytesIO()

        tgz = self.open('w')
        tgz.add(
            dirname,
            arcname='.',
            filter=self.itpfilter
        )
        tgz.close()

        super().__init__(fname=fname, **kwargs)

    def open(self, mode='r'):
        self._rawfile.seek(0)
        if mode == 'w':
            mode = 'w:'
        if mode == 'r':
            mode = 'r:'
        return tarfile.open(
            fileobj=self._rawfile,
            mode=mode
        )

    def write(self, f):
        if isinstance(f, io.IOBase):
            fname = f.name
        else:
            fname = f
        path = os.path.dirname(fname)
        self.open().extractall(path)
        super().write(f)

    @staticmethod
    def itpfilter(tarinfo):
        path = tarinfo.name
        if path.lower().endswith('.itp'):
            return tarinfo
        elif tarinfo.isdir() and path == '.':
            return tarinfo
        elif tarinfo.isdir() and path.endswith('.ff'):
            return tarinfo
        else:
            return None

    def add(self, path, arcname=None):
        tgz = self.open('a')
        tgz.add(
            path,
            arcname=arcname,
            filter=self.itpfilter
        )
        tgz.close()

    def include(self, fn):
        if os.path.exists(fn):
            self.add(fn)
        return super().include(fn)




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

WriteInfo = namedtuple('WriteInfo', ['pdb', 'top', 'pickle', 'index'])


class GmxCaller():
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

    @classmethod
    def call_gmx(cls, cmdline=None, cmd=None, stdin='', timeout=None, cwd=None, mpiranks=None, **kwargs):
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
            args = [('gmx', cmd)] + [cls.kwarg_to_arg(k,v) for k,v in kwargs.items()]
            args = [item for sublist in args for item in sublist]

        cmd = args[1]
        cmd_is_mdrun = cmd in mdrun_synonyms

        if mpiranks:
            args = [
                'mpirun',
                '--use-hwthread-cpus',
                '-np', str(mpiranks),
                'mdrun_mpi'
            ] + args[2:]

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
                            sleep(1)
                    s = f.read()
                    sys.stdout.write(s)

            else:
                p.wait()
            e.seek(0)
            o.seek(0)
            stdout_data = o.read()
            stderr_data = e.read()
        if p.returncode:
            raise ValueError(f'GROMACS returned error code {p.returncode}; stderr was:\n{stderr_data}\n\nContents of cwd are:\n{os.listdir(cwd)}')

        return (p, stdout_data, stderr_data)


    @classmethod
    def mdrun(cls, cwd='.', **kwargs):
        cwd = Path(cwd).absolute()

        cmd = 'mdrun_mpi' if 'mpiranks' in kwargs else 'mdrun'

        try:
            kwargs['multidir'] = list(kwargs['multidir'])
        except KeyError:
            pass

        try:
            trajout = kwargs['x']
        except KeyError:
            try:
                trajout = kwargs['o']
            except KeyError:
                try:
                    trajout = kwargs['deffnm'] + '.xtc'
                except KeyError:
                    trajout = 'traj_comp.xtc'

        name, _, ext = trajout.rpartition('.')
        if ext.lower() not in ('xtc', 'trr', 'tng'):
            raise ValueError('Extension for output trajectory not specified')

        try:
            tprin = kwargs['s']
        except KeyError:
            try:
                tprin = kwargs['deffnm'] + '.tpr'
            except KeyError:
                tprin = 'topol.tpr'
        if tprin[-4:] != '.tpr':
            tprin = tprin + '.tpr'

        kwargs.setdefault('pin', 'on')
        kwargs.setdefault('nb', 'gpu')
        kwargs.setdefault('pme', 'gpu')
        kwargs.setdefault('v', 'True')

        cls.call_gmx(
            cmd=cmd,
            cwd=cwd,
            **kwargs
        )

        try:
            multidir = kwargs['multidir']
        except KeyError:
            try:
                trajs = (cwd / f'{name}{i}.{ext}' for i in range(int(kwargs['multi'])))
                tprs = (cwd / f'{name}{i}.tpr' for i in range(int(kwargs['multi'])))
            except KeyError:
                trajs = (cwd / trajout, )
                tprs = (cwd / tprin, )
        else:
            if isinstance(multidir, str):
                multidir = multidir.split()
            trajs = (cwd / d / trajout for d in multidir)
            tprs = (cwd / d / tprin for d in multidir)

        for traj, tpr in zip(trajs, tprs):
            GmxCaller.trajvis(traj, tpr)


    @classmethod
    def trajvis(cls, filename, tprfnm):
        filename = Path(filename).absolute()
        tprfnm = Path(tprfnm).absolute()
        trajout_tmp = f'{filename.stem}.tmp{filename.suffix}'
        trajout = f'{filename.stem}.vis{filename.suffix}'
        pdbout = f'{filename.stem}.vis.first.pdb'

        with TemporaryDirectory() as path:
            cls.call_gmx(
                cmd='trjconv',
                stdin='Protein\nSystem\n',
                f=filename,
                s=tprfnm,
                o=trajout_tmp,
                pbc='mol',
                ur='compact',
                center=True
            )
            if filename.suffix.lower() == '.xtc':
                cls.call_gmx(
                    cmd='trjconv',
                    stdin='System\n',
                    f=trajout_tmp,
                    s=tprfnm,
                    o=pdbout,
                    dump=0
                )
                cls.call_gmx(
                    cmd='trjconv',
                    stdin='C-alpha\nSystem\n',
                    f=trajout_tmp,
                    s=tprfnm,
                    o=trajout,
                    fit='trans'
                )
                os.remove(trajout_tmp)
            else:
                os.rename(trajout_tmp, trajout)


class GmxSystem(GmxCaller):
    def __init__(
        self,
        name,
        traj,
        ffpath,
        mdmdp,
        emmdp=None,
        watermodel='tips3p'
    ):
        super().__init__()

        self.name = name
        self._traj = traj
        self._velxyz = None

        if ffpath is None:
            self.ff_path = None
            self.ff = None
        else:
            self.ff_path = Path(ffpath).absolute()
            self.ff = GMXForceField(self.ff_path)
        self.watermodel = watermodel
        boxbuffer = 1.2

        self._pickle_path = None

        self.mdp = mdmdp.copy()

        if emmdp is None:
            emmdp = mdmdp.copy()
            emmdp.integrator = 'steep'

        self.em_mdp = emmdp.copy()

        self.ndx = None


    def martinize_m22(self, name=None, cwd='.', dssp='mkdssp', ss=None):
        cwd = Path(cwd).absolute()
        if name is None:
            name = self.name

        forcefield_itp = '''
        ;
        ; Martini 2.2 force field
        ;

        #include "martini_v2.2.itp"
        #include "martini_v2.2_aminoacids.itp"
        #include "martini_v2.0_ions.itp"
        '''.strip()

        self.ff_path = cwd / 'martini22.ff'
        os.makedirs(self.ff_path, exist_ok=True)
        urlretrieve(
            "http://cgmartini.nl/images/parameters/ITP/martini_v2.2.itp",
            self.ff_path / "martini_v2.2.itp"
        )
        urlretrieve(
            "http://cgmartini.nl/images/parameters/ITP/martini_v2.2_aminoacids.itp",
            self.ff_path / "martini_v2.2_aminoacids.itp"
        )
        urlretrieve(
            "http://cgmartini.nl/images/parameters/ITP/martini_v2.0_ions.itp",
            self.ff_path / "martini_v2.0_ions.itp"
        )
        with open(self.ff_path / 'forcefield.itp', 'w') as f:
            print(forcefield_itp, file=f)
        self.ff = GMXForceField(self.ff_path)

        pdbin = cwd / f'{name}.pdb'
        pdbout = cwd / f'{name}_cg.pdb'
        topout = cwd / f'{name}.top'
        ndxout = cwd / f'{name}_map.ndx'

        self.save_traj(pdbin)

        scriptname = 'martinize.py'
        copyfile(
            Path(unwrap(JupyterJoy.__path__)) / 'bin' / 'martinize.py',
            scriptname
        )

        args = [
            'python3', scriptname,
            '-f', pdbin,
            '-o', topout,
            '-x', pdbout,
            '-nmap', ndxout,
            '-name', name,
            '-ff', 'martini22',
            '-elastic',
            '-p', 'backbone'
        ]

        if ss is None:
            args += ['-dssp', dssp]
        else:
            args += ['-ss', ss]

        p = sp.run(
            args,
            encoding='utf-8',
            stdout=sp.PIPE,
            stderr=sp.STDOUT,
            cwd=cwd
        )
        print(p.stdout)
        p.check_returncode()

        self.topol = TopolWithItps(topout)
        self.topol.includes[0] = '#include "martini22.ff/forcefield.itp"'
        Topology.write(self.topol, topout)

        self.load_traj(pdbout)

        self.ndx = make_martini_ndx(self.traj.top)


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

    def optimise_rdsq_box(self, target_mindist):
        # Orient the protein's long axis with the x-axis - this is the worst-case orientation
        self.orient_protein([1, 0, 0])

        self.traj.unitcell_vectors = None
        guess = target_mindist / 2.0
        d = 2 * guess + self.calc_maxdist()
        self.add_rdsq_box(d)
        this_mindist = self.calc_mindist()

        print(f'A buffer of {guess} nm gives a box with d={d}. The minimum PI distance is {this_mindist}.')

        # And reorient to the z-axis - best-case orientation
        self.orient_protein([0, 0, 1])

        print(f'After reorientation, the PI distance is {self.calc_mindist()}.')


    def add_rdsq_box(self, d):
        """Give traj a rhombic dodecahedral box with square in the xy plane and vector length d"""
        d = float(d)
        a = (d, 0.0, 0.0)
        b = (0.0, d, 0.0)
        c = (d/2.0, d/2.0, sqrt(2.0)*d/2.0)
        vecs = np.array(len(self.traj) * [[a, b, c]])
        if self.traj.unitcell_vectors is not None:
            raise ValueError(f'self.traj already has a unit cell: {self.traj.unitcell_vectors}')
        self.traj.unitcell_vectors = vecs
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

    @property
    def velxyz(self):
        return self._velxyz

    @velxyz.setter
    def velxyz(self, velxyz):
        if velxyz is not None and velxyz.shape != self.traj.xyz.shape:
            raise ValueError('velxyz should have same shape as self.traj.xyz')
        self._velxyz = velxyz

    def write_all(self, path, suffix=None):
        topfn, pdbfn, *_ = self.write(path, suffix=suffix)
        indexfn = self.write_ndx(path, suffix=suffix)
        picklefn = self.pickle(path)
        return WriteInfo(pdb=pdbfn, top=topfn, pickle=picklefn, index=indexfn)

    def write(self, path, suffix=None):
        topfn = self.write_top(path, suffix=suffix)
        pdbfn = self.write_traj(path, suffix=suffix)
        return WriteInfo(pdb=pdbfn, top=topfn, pickle=None, index=None)

    def write_top(self, path, suffix=None):
        path = Path(path).absolute()
        if suffix is None:
            suffix = ''
        else:
            suffix = '.' + str(suffix)
        self.ff.write(path=path)
        topfn = path / f'{self.name}{suffix}.top'
        self.topol.write(topfn)
        return topfn

    def write_traj(self, path, suffix=None):
        if self.velxyz is None:
            return self.write_pdb(path, suffix=suffix)
        else:
            return self.write_gro(path, suffix=suffix)

    def write_pdb(self, path, suffix=None):
        path = Path(path).absolute()
        if suffix is None:
            suffix = ''
        else:
            suffix = '.' + str(suffix)
        pdbfn = path / f'{self.name}{suffix}.pdb'
        self.save_traj(pdbfn)
        return pdbfn

    def write_gro(self, path, suffix=None):
        path = Path(path).absolute()
        if suffix is None:
            suffix = ''
        else:
            suffix = '.' + str(suffix)
        grofn = path / f'{self.name}{suffix}.gro'
        self.save_traj(grofn)
        return grofn

    def write_ndx(self, path, suffix=None):
        path = Path(path).absolute()
        if suffix is None:
            suffix = ''
        else:
            suffix = '.' + str(suffix)
        ndxfn = path / f'{self.name}{suffix}.ndx'
        with open(ndxfn, 'w') as f:
            if isinstance(self.ndx, dict):
                write_dict_to_ndx(self.ndx, f)
            else:
                print(self.ndx, file=f)
        return ndxfn

    def pickle(self, path):
        path = Path(path).absolute()
        picklefn = path / f'{self.name}.pickle'
        self._pickle_path = str(picklefn)
        with open(picklefn, 'w+b') as f:
            pickle.dump(self, f, protocol=4)
        return picklefn

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['ff']
        del state['topol']
        return state

    def __setstate__(self, state):
        path, _, basename = state['_pickle_path'].rpartition('/')
        name, _, ext = basename.rpartition('.')
        if ext != 'pickle' or name != state['name']:
            raise ValueError(f'{self.__class__} object not stored via its pickle method')
        ff_name = os.path.basename(state['ff_path'])
        state['ff'] = GMXForceField(f'{path}/{ff_name}')
        state['topol'] = TopolWithItps(f'{path}/{name}.top')
        state['_pickle_path'] = None
        state.setdefault('_velxyz', None)
        self.__dict__.update(state)

    @classmethod
    def read(cls, path):
        with open(path, 'r+b') as f:
            self = pickle.load(f)
        return self

    def solvate(self, **kwargs):
        with TemporaryDirectory() as path:
            pdbin, topinout, *_ = self.write(path)
            pdbout = 'extended_sol.pdb'
            _, stdout_data, stderr_data = self.call_gmx(
                cmd='solvate',
                stdin='',
                cwd=path,
                cp=pdbin,
                o=pdbout,
                p=topinout,
                **kwargs
            )
            print(stderr_data)
            print(stdout_data)
            self.topol = TopolWithItps(topinout)
            self.load_traj(f'{path}/{pdbout}')
            self.ndx = make_martini_ndx(self.traj.top)

    def salt_m22(self, conc=0.15, solname='W', pname='NA+', nname='CL-', **kwargs):
        return self.salt(conc*4, solname=solname, pname=pname, nname=nname, **kwargs)

    def salt(self, conc=0.15, solname='SOL', **kwargs):
        with TemporaryDirectory() as path:
            pdbin, topinout, *_ = self.write(path)
            pdbout = 'extended_ion.pdb'
            tpr_tmp = 'genion.tpr'

            n_wat = sum(m.copies for m in self.topol.molecules[solname])
            mol_water = 55.5
            desired_salt_conc = conc
            n_salt = int(round(n_wat * desired_salt_conc / mol_water))
            with open(f'{path}/genion.mdp', 'w') as f:
                print('coulombtype = cutoff', file=f)

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
                stdin=solname,
                cwd=path,
                s=tpr_tmp,
                o=pdbout,
                p=topinout,
                np=str(n_salt),
                neutral=True,
                **kwargs
            )
            print(stderr_data)
            print(stdout_data)
            self.topol = TopolWithItps(topinout)
            self.load_traj(f'{path}/{pdbout}')
            self.ndx = make_martini_ndx(self.traj.top)

    def add_solute(self, proportion=0.1, solvent='W', solute='WF', **kwargs):
        with TemporaryDirectory() as path:
            pdbin, topinout, *_ = self.write(path)
            pdbout = 'extended_ion.pdb'
            tpr_tmp = 'genion.tpr'

            n_solvent = unwrap(self.topol.molecules[solvent]).copies
            n_solute = int(n_solvent * proportion)
            self.em_mdp.write(f'{path}/genion.mdp')

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
                stdin=solvent,
                cwd=path,
                s=tpr_tmp,
                o=pdbout,
                p=topinout,
                np=str(n_solute),
                nn='0',
                pname=solute,
                nname=solvent,
                neutral=False,
                **kwargs
            )
            print(stderr_data)
            print(stdout_data)
            self.topol = TopolWithItps(topinout)
            self.load_traj(f'{path}/{pdbout}')
            self.ndx = make_martini_ndx(self.traj.top)



    def em(self):
        with TemporaryDirectory() as path:
            pdbin, topinout, *_ = self.write(path)
            pdbout = 'em.pdb'
            tpr_tmp = 'em.tpr'

            self.em_mdp.write(f'{path}/em.mdp')

            print(f'Energy minimising with .MDP file:\n{self.em_mdp.write()}')

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
            self.load_traj(f'{path}/{pdbout}')

    def trajvis(self, filename, tprfnm=None):
        filename = Path(filename).absolute()

        with TemporaryDirectory() as path:
            path = Path(path)
            if tprfnm is None:
                pdbin, topin, *_ = self.write(path)
                mdp_tmp = path / 'trjconv.mdp'
                tprfnm = path / 'trjconv.tpr'
                open(mdp_tmp, 'a').close()
                self.call_gmx(
                    cmd='grompp',
                    stdin='',
                    f=mdp_tmp,
                    c=pdbin,
                    p=topin,
                    o=tprfnm
                )
            out = super().trajvis(filename, tprfnm)
        return out


    def load_pdb(self, filename):
        return self.load_traj(self, filename)


    def load_traj(self, filename):
        filename = Path(filename).absolute()
        with TemporaryDirectory() as path:
            path = Path(path)
            topin = self.write_top(path)
            groout = f'{self.name}.vis.gro'
            self.em_mdp.write(path / 'trjconv.mdp')
            tpr_tmp='trjconv.tpr'
            _, stdout_data, stderr_data = self.call_gmx(
                cmd='grompp',
                stdin='',
                cwd=path,
                f='trjconv.mdp',
                c=filename,
                p=topin,
                o=tpr_tmp,
                maxwarn=1
            )
            _, stdout_data, stderr_data = self.call_gmx(
                cmd='trjconv',
                stdin='Protein\nSystem\n',
                cwd=path,
                f=filename,
                s=tpr_tmp,
                o=groout,
                pbc='mol',
                ur='compact',
                center=True
            )
            self.traj, self.velxyz = load_gro(path / groout)

    def save_traj(self, filename):
        filename = Path(filename).absolute()
        if filename.suffix.lower() == '.gro':
            return self.save_gro(filename)
        elif self.velxyz is not None:
            print(f'Warning: {filename} will not include velocities')
        return self.traj.save(str(filename))


    def save_gro(self, filename, *args, **kwargs):
        return save_gro(filename, self.traj, velxyz=self.velxyz, *args, **kwargs)

    def setup_sim(self, mdp, deffnm, cwd='.', usecoordsforposre=False, **kwargs):
        cwd = Path(cwd).absolute()
        pdbin, topinout, *_ = self.write(cwd)
        tprinout = cwd / f'{deffnm}.tpr'
        mdpinout = cwd / f'{deffnm}.mdp'

        if mdp is None:
            mdp = self.mdp

        mdp.write(mdpinout)
        print(f'Running simulation with .MDP file:\n{mdp.write()}')

        defkwargs = dict(
            f=mdpinout,
            c=pdbin,
            p=topinout,
            o=tprinout
        )
        if self.ndx is not None:
            ndxin = self.write_ndx(path=cwd)
            defkwargs['n'] = ndxin

        defkwargs.update(kwargs)
        kwargs = defkwargs

        if usecoordsforposre and 'r' not in kwargs:
            kwargs['r'] = pdbin
        elif usecoordsforposre:
            raise ValueError('Specify r in gromppkwargs or usecoordsforposre, not both')

        _, stdout_data, stderr_data = self.call_gmx(
            cmd='grompp',
            stdin='',
            cwd=cwd,
            **kwargs
        )
        print(stderr_data)
        print(stdout_data)

    def run_sim(self, mdp, deffnm, cwd='.', usecoordsforposre=False, nreps=1, gromppkwargs=None, **mdrunkwargs):
        cwd = Path(cwd).absolute()
        if gromppkwargs is None:
            gromppkwargs = {}

        if nreps <= 0:
            raise ValueError('nreps should be positive')
        elif nreps == 1:
            workingdirs = [cwd]
        else:
            workingdirs = [cwd / f'rep{i}' for i in range(nreps)]
            mdrunkwargs['multidir'] = workingdirs
            mdrunkwargs['mpiranks'] = nreps
        for path in workingdirs:
            os.makedirs(cwd, exist_ok=True)
            self.setup_sim(mdp, deffnm, cwd=path, usecoordsforposre=usecoordsforposre, **gromppkwargs)

        self.mdrun(
            cwd=cwd,
            deffnm=deffnm,
            **mdrunkwargs
        )

    def load_xtc(self, filename, **kwargs):
        return md.load_xtc(filename, top=self.traj.top, **kwargs)


    def get_properties(self, deffnm, cwd='.', stride=1, mdp=None, mindist=True):
        if mdp is None:
            mdp = self.mdp

        if mindist:
            self.call_gmx(
                cmd='mindist',
                stdin='Protein',
                cwd=cwd,
                f=f'{deffnm}.xtc',
                s=f'{deffnm}.tpr',
                od=f'{deffnm}.mindist.xvg',
                pi=True,
                dt=float(mdp.nstenergy) * float(mdp.dt) * stride
            )
            with open(f'{cwd}/{deffnm}.mindist.xvg') as f:
                data = np.array([
                    [float(s) for s in l.split()]
                    for l in f if l[0] not in '#@'
                ]).T
        df = edr_to_df(f'{cwd}/{deffnm}.edr')[::stride]
        traj = self.load_xtc(f'{cwd}/{deffnm}.vis.xtc', stride=stride)

        if mindist:
            minlen = min(data.shape[1], len(df), len(traj))
        else:
            minlen = min(len(df), len(traj))

        df = df.head(minlen)
        if mindist:
            data = data[..., :minlen]
        traj = traj[:minlen]
        if not (
            ((not mindist) or np.array_equal(df['Time'], data[0]))
            and np.array_equal(df['Time'], traj.time)
        ):
            raise ValueError("Could not match times across different inputs")

        calpha_atom_indices = traj.top.select_atom_indices('alpha')
        rmsd = md.rmsd(traj, self.traj, atom_indices=calpha_atom_indices)

        if mindist:
            df['Min. PI dist'] = data[1]
            df['Max. int dist'] = data[2]
        df['RMSD'] = rmsd

        return (traj, df)

class GmxRest2System(GmxSystem):
    def __init__(
        self,
        name,
        traj,
        ffpath,
        mdmdp,
        watermodel='tips3p',
        emmdp=None,
        min_temp=300.0,
        max_temp=300.0,
        num_reps=1,
        exchange_freq=10,
    ):
        self._ladder = []
        self._min_temp = min_temp
        self._max_temp = max_temp
        self._num_reps = num_reps
        self._ladder_method = 'GEOMETRIC'
        self.exchange_freq = exchange_freq
        self._nstlist = None

        super().__init__(
            name=name,
            traj=traj,
            ffpath=ffpath,
            mdmdp=mdmdp,
            emmdp=emmdp,
            watermodel=watermodel
        )
        self.mdp.set_temperature(self.min_temp)


    @property
    def ladder(self):
        if self._ladder:
            return self._ladder

        return calc_temps(
            self.min_temp,
            self.max_temp,
            self.num_reps,
            method=self.ladder_method
        )

    @property
    def min_temp(self):
        if self._ladder:
            return self._ladder[0]

        return self._min_temp

    @property
    def max_temp(self):
        if self._ladder:
            return self._ladder[-1]

        return self._max_temp

    @property
    def num_reps(self):
        if self._ladder:
            return len(self._ladder)

        return self._num_reps

    @property
    def ladder_method(self):
        if self._ladder:
            return 'SPECIFIED'

        return self._ladder_method

    @property
    def nstlist(self):
        ex_freq = self.exchange_freq
        nstlist = self._nstlist

        if nstlist and ex_freq % nstlist == 0:
            return nstlist
        elif nstlist:
            raise ValueError(
                "nstlist {nstlist} doesn't divide exchange_freq {ex_freq}"
            )
        elif ex_freq % 100 == 0:
            return 100
        else:
            raise ValueError('You should manually specify nstlist')

        return self._ladder_method

    @nstlist.setter
    def nstlist(self, value):
        ex_freq = self.exchange_freq
        if value % ex_freq == 0:
            self._nstlist = value
        else:
            raise ValueError(
                "nstlist {value} doesn't divide exchange_freq {ex_freq}"
            )

    @nstlist.deleter
    def nstlist(self, value):
        self._nstlist = None

    def prep_rest2(self, deffnm, path, mdp, startframes):
        rpath = Path(path).absolute()

        if isinstance(startframes, md.Trajectory) and len(startframes) == 1:
            frame_iter = cycle(startframes)
        elif self.num_reps == len(startframes):
            frame_iter = iter(startframes)
        else:
            raise ValueError(f'Give either one frame or {self.num_reps} frames')

        rpath.mkdir(parents=True, exist_ok=True)
        for t in self.ladder:
            if isinstance(t, str):
                tstr = t
            else:
                tstr = f'{t:.2f}'
            tpath = rpath /tstr
            tpath.mkdir(parents=True, exist_ok=True)

            coordin = tpath / f'{deffnm}.start.gro'
            frame = next(frame_iter)
            try:
                frame.save(coordin)
            except AttributeError:
                coordin = Path(frame).absolute()

            mdpin = tpath / f'{deffnm}.mdp'
            if mdp is None:
                mdp = self.mdp
            mdp = mdp.copy()
            mdp.set_temperature(self.min_temp)
            mdp.write(mdpin)

            topin = self.write_top(tpath)
            topinout = tpath / f'{deffnm}.rest2_{tstr}{topin.suffix}'
            tprinout = tpath / f'{deffnm}.tpr'

            with TemporaryDirectory() as cwd:
                cwd = Path(cwd).absolute()
                preprocessed_top = 'preproc.top'
                _, stdout_data, stderr_data = self.call_gmx(
                    cmd='grompp',
                    stdin='',
                    cwd=cwd,
                    f=mdpin,
                    c=coordin,
                    p=topin,
                    pp=cwd / preprocessed_top,
                    maxwarn=1
                )

                temper_from_file(
                    cwd / preprocessed_top,
                    topinout,
                    ['Protein_chain_A'],
                    float(self.min_temp),
                    float(t)
                )

            with open(tpath / 'plumed.dat', 'w') as f:
                f.write('# Necessary for hrex')

            _, stdout_data, stderr_data = self.call_gmx(
                cmd='grompp',
                stdin='',
                cwd=tpath,
                f=mdpin,
                c=coordin,
                p=topinout,
                o=tprinout,
                maxwarn=1
            )
            print(stderr_data)
            print(stdout_data)

    def take_starting_strucs(self, traj, mdp=None, skiptime_ps=100):
        """Take starting structures spaced throughout traj"""
        if mdp is None:
            mdp = self.mdp

        nframes = len(traj)
        skipframes = int(skiptime_ps / mdp.dt) // int(mdp.nstxout_compressed)
        strideframes = (nframes - skipframes) // (self.num_reps - 1)
        temp_ladder_idcs = (strideframes * i + skipframes for i in range(self.num_reps))
        return [traj[i] for i in temp_ladder_idcs]



