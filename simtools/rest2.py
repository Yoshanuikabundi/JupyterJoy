import os
from collections import defaultdict, namedtuple
from itertools import chain, repeat, combinations, zip_longest
import re
from copy import copy, deepcopy

import pandas as pd
import numpy as np


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


class MoleculeType:
    def __init__(self, data, top=None):
        name, nrexcl = data.strip().split()
        self.name = name
        self.nrexcl = nrexcl
        self.top = defaultdict(list) if top is None else top

    def __repr__(self):
        return f'MoleculeType("{self.name} {self.nrexcl}", {self.top})'


class CMapType:
    def __init__(self, *args):
        if len(args) == 1 and isinstance(args[0], str):
            data = args[0].split()
        else:
            data = args
        (
            self.i,
            self.j,
            self.k,
            self.l,
            self.m,
            self.func,
            self.gridx,
            self.gridy,
            *self.values
        ) = data
        assert int(self.gridx) * int(self.gridy) == len(self.values)
        self.values = [float(v) for v in self.values]

    def __repr__(self):
        return (
            f'<CMapType for atoms {self.i} {self.j} {self.k} {self.l} '
            f'{self.m} on a {self.gridx}x{self.gridy} grid at {hex(id(self))}>'
        )


def infer_dtypes(df):
    new_df = pd.DataFrame()
    for key, series in df.items():
        try:
            series_out = pd.to_numeric(series)
        except ValueError:
            series_out = series.astype(str)
        new_df[key] = series_out
    return new_df


def top_dict_to_dfs(top_dict, recursing=False):
    top_df_dict = {}
    for directive, value in top_dict.items():
        if recursing and directive in ['moleculetype', 'system', 'molecules']:
            raise ValueError(f'Should be no "[ {directive} ]" in a moleculetype')
        if directive in ['cmaptypes']:
            top_df_dict[directive] = value
        elif directive in ['defaults', 'system']:
            assert len(value) == 1
            top_df_dict[directive] = value[0]
        elif directive == 'moleculetype':
            top_df_dict[directive] = {mol.name: MoleculeType(f'{mol.name} {mol.nrexcl}', top_dict_to_dfs(mol.top, recursing=True)) for mol in value}
        else:
            top_df_dict[directive] = infer_dtypes(pd.DataFrame(value))

    return top_df_dict


def read_top(fname):
    with open(fname) as f:
        top = f.readlines()

    top_dict = defaultdict(list)

    current_directive = ''
    current_moleculetype = ''
    linecontinue = ''
    dtype = None

    dtype_dct = {}

    for line in top:
        line = linecontinue + line
        if line.endswith('\\\n'):
            linecontinue = f'{line[:-2]} '
            continue
        else:
            linecontinue = ''

        s = line.strip()
        if s.startswith('*'):
            data, comment = '', line
        else:
            data, _, comment = s.partition(';')
        data = data.strip()

        if dtype == "GET_FROM_COMMENT":
            if comment and not data or current_directive in dtype_dct:
                if current_directive in dtype_dct:
                    namedtup = dtype_dct[current_directive]
                else:
                    fieldnames, _, comment = comment.partition(';')
                    fieldnames = fieldnames.strip().split()
                    fields = list(map(lambda s: re.sub(r'\W|^(?=\d)', '_', s), fieldnames))
                    try:
                        namedtup = namedtuple(current_directive, fields)
                    except ValueError as err:
                        print(err)
                        continue
                    dtype_dct[current_directive] = namedtup
                def dtype(data):
                    l = data.split()
                    n_fields = len(namedtup._fields)
                    n_nones = n_fields - len(l)
                    out = namedtup(*chain(l, repeat(None, n_nones)))
                    return out
            else:
                dtype = str.split

        if data.startswith('[') and data.endswith(']'):
            current_directive = data[1:-1].strip()
            if current_directive == 'system':
                dtype = str
            elif current_directive == 'cmaptypes':
                dtype = CMapType
            else:
                dtype = "GET_FROM_COMMENT"
            continue

        if not data:
            continue
        elif not current_directive:
            raise ValueError(f'Line "{line}" appears to have data outside a directive')

        if current_directive == 'moleculetype':
            current_moleculetype = data.split()[0]
            dtype = MoleculeType
        elif current_directive == 'system':
            current_moleculetype = ''

        if current_moleculetype and current_directive != 'moleculetype':
            top_dict['moleculetype'][-1].top[current_directive].append(dtype(data))
        else:
            top_dict[current_directive].append(dtype(data))

    return top_dict_to_dfs(top_dict)


def temper_df(df, name_columns, energy_columns, factor):
    tempered_df = df.copy()

    for col in name_columns:
        tempered_df[col] += '🔥'
        tempered_df[col].replace('X🔥', 'X', inplace=True)

    for col in energy_columns:
        tempered_df[col] *= factor

    return tempered_df


def temper_combinations(df, name_columns, energy_columns, b0, bm, all_combos=False):
    tempered_df_pp = temper_df(df, name_columns, energy_columns, bm/b0)

    dfs = [df, tempered_df_pp]

    if all_combos:
        for n in range(2, len(name_columns) + 1):
            for t in combinations(name_columns, n):
                dfs.append(temper_df(df, t, energy_columns, np.sqrt(bm/b0)))

    return pd.concat(dfs)


def temper_atomtypes(df, t0, tm):
    b0, bm = 1/t0, 1/tm

    tempered_df = df.copy()
    tempered_df.name += '🔥'
    tempered_df.charge *= np.sqrt(bm/b0)
    tempered_df.epsilon *= bm/b0

    return pd.concat([df, tempered_df])


def temper_pairtypes(df, t0, tm, all_combos=False):
    b0, bm = 1/t0, 1/tm

    names = ['i', 'j']
    energies = ['epsilon1_4']

    return temper_combinations(df, names, energies, b0, bm, all_combos=all_combos)


def temper_bondtypes(df, t0, tm, all_combos=False):
    b0, bm = 1/t0, 1/tm

    names = ['i', 'j']
    energies = ['kb']

    return temper_combinations(df, names, energies, b0, bm, all_combos=all_combos)


def temper_angletypes(df, t0, tm, all_combos=False):
    b0, bm = 1/t0, 1/tm

    if (df.func != 5).any():
        raise ValueError("temper_angletypes presently only supports urey-bradley angles")

    names = ['i', 'j', 'k']
    energies = ['cth', 'cub']

    return temper_combinations(df, names, energies, b0, bm, all_combos=all_combos)


def temper_dihedraltypes(df, t0, tm, all_combos=False):
    b0, bm = 1/t0, 1/tm

    if not df.func.isin([2,9]).all():
        raise ValueError("temper_dihedraltypes presently only supports CHARMM-style dihedrals")

    names = ['i', 'j', 'k', 'l']
    energies = ['cp']

    return temper_combinations(df, names, energies, b0, bm, all_combos=all_combos)


def temper_cmaptypes(cmaps, t0, tm, all_combos=False):
    b0, bm = 1/t0, 1/tm

    cmaps_out = list(cmaps)

    for cmap in cmaps:
        cmap_list = [
            cmap.i,
            cmap.j,
            cmap.k,
            cmap.l,
            cmap.m,
            cmap.func,
            cmap.gridx,
            cmap.gridy
        ]

        lit_cmap = copy(cmap_list)
        for i in range(5):
            lit_cmap[i] += '🔥'
        lit_cmap.extend(v * np.sqrt(bm/b0) for v in cmap.values)
        cmaps_out.append(CMapType(*lit_cmap))

        if all_combos:
            for n in range(1, 5):
                for t in combinations(range(5), n):
                    lit_cmap = copy(cmap_list)
                    for i in t:
                        lit_cmap[i] += '🔥'
                    lit_cmap.extend(v * bm/b0 for v in cmap.values)
                    cmaps_out.append(CMapType(*lit_cmap))
    return cmaps_out


def temper_forcefield(top_dict, t_0, t_m, all_combos=False):
    """Copy everything outside moleculetypes and scale energies"""
    assert top_dict['defaults'].nbfunc == '1'
    assert top_dict['defaults'].comb_rule == '2'
    dct = deepcopy(top_dict)
    dct['atomtypes'] = temper_atomtypes(top_dict['atomtypes'], t_0, t_m)
    dct['pairtypes'] = temper_pairtypes(top_dict['pairtypes'], t_0, t_m, all_combos=all_combos)
    dct['bondtypes'] = temper_bondtypes(top_dict['bondtypes'], t_0, t_m, all_combos=all_combos)
    dct['angletypes'] = temper_angletypes(top_dict['angletypes'], t_0, t_m, all_combos=all_combos)
    dct['dihedraltypes'] = temper_dihedraltypes(top_dict['dihedraltypes'], t_0, t_m, all_combos=all_combos)
    dct['cmaptypes'] = temper_cmaptypes(top_dict['cmaptypes'], t_0, t_m, all_combos=all_combos)
    return dct


def temper_moltype(moltype, t0, tm):
    b0, bm = 1/t0, 1/tm

    moltype.top['atoms']['type'] += '🔥'
    moltype.top['atoms']['charge'] *= np.sqrt(bm/b0)

    return moltype


def temper_top(top_dict, temper_mols, t0, tm, all_combos=False):
    tempered = temper_forcefield(top_dict, t0, tm, all_combos=all_combos)

    for mol in temper_mols:
        tempered['moleculetype'][mol] = temper_moltype(tempered['moleculetype'][mol], t0, tm)

    del tempered['implicit_genborn_params']

    return tempered


def format_float(f):
    if np.isnan(f):
        return ''
    else:
        return f'  {f: > .6f}'


def write_df(df, f):
    s = df.to_string(
        index=False,
        justify='left',
        float_format=format_float
    )
    if s[0] != ' ':
        lst = [" " + line for line in s.splitlines()]
        s = '\n'.join(lst)
    s = ';' + s[1:] + '\n'
    f.write(s)


def write_namedtup(nt, f):
    header = ';'
    data = ' '
    for key, value in zip(nt._fields, nt):
        size = max(len(key), len(value)) + 2
        header += f'{{: >{size}}}'.format(key)
        data += f'{{: >{size}}}'.format(value)
    f.write(f'{header}\n{data}\n')


def write_cmaptypes(cmaps, f):
    for cmap in cmaps:
        s = ' '.join(f'{s}' for s in [
            cmap.i,
            cmap.j,
            cmap.k,
            cmap.l,
            cmap.m,
            cmap.func,
            cmap.gridx,
            cmap.gridy
        ])
        f.write(f'{s}\\\n')
        lines = [' '.join(f'{v: > .8f}' for v in t) for t in grouper(cmap.values, int(cmap.gridx))]
        f.write('\\\n'.join(lines))
        f.write('\n\n')


def write_moleculetypes(moltypes, f):
    for i, moltype in enumerate(moltypes.values()):
        if i > 0:
            f.write('[ moleculetype ]\n')
        f.write('; Name            nrexcl\n')
        f.write(f'{moltype.name} {moltype.nrexcl}\n\n')

        for directive, data in moltype.top.items():
            write_directive(directive, data, f)
        f.write('\n')


def write_directive(directive, data, f):
    f.write(f'[ {directive} ]\n')

    write_funcs = dict(
        defaults=write_namedtup,
        cmaptypes=write_cmaptypes,
        moleculetype=write_moleculetypes,
        system=lambda data, f: f.write(data + '\n\n')
    )

    if isinstance(data, pd.DataFrame):
        write_df(data, f)
    elif directive in write_funcs:
        write_funcs[directive](data, f)
    else:
        raise ValueError(f"Don't know how to write {directive}")

    f.write('\n')


def temper_from_file(fname, fname_out, temper_mols, t0, tm, all_combos=False):
    top_dict = read_top(fname)
    tempered = temper_top(top_dict, temper_mols, t0, tm, all_combos=all_combos)

    with open(fname_out, 'w') as f:
        for directive, data in tempered.items():
            write_directive(directive, data, f)
