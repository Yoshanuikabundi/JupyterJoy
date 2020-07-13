from collections import OrderedDict, Mapping
from .rstparser import parse_rst_file, process_name
import os
import re
from numbers import Number
from copy import copy, deepcopy
# from JupyterJoy.pbash import *
import mdtraj as md
import shlex

tail_comment = "__TAIL_COMMENT__"

class Default():
    """Class to signal that a setting hasn't been specified"""
    def __init__(self, object):
        self.obj = object

    def __repr__(self):
        return "Default({})".format(repr(self.obj))

    def __bool__(self):
        return False

    @classmethod
    def unwrap(cls, obj):
        if isinstance(obj, cls):
            return obj.obj
        else:
            return obj

class MDPBase():
    """Base class for creating MDP files"""
    obsoletes = set([])
    @classmethod
    def create_subclass_from_doc(cls, name, docpath):
        """Create an MDP class from a mdp-options.rst file distributed with GROMACS"""
        mdp_entries, doc = parse_rst_file(docpath)
        attr = {
            'options': mdp_entries,
            'doc': doc,
            'obsoletes': copy(cls.obsoletes)
        }
        return type(name, (cls,), attr)

    def __init__(self, data, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__dict__['values'] = {}
        self.__dict__['comments'] = {}
        self.__dict__['meta'] = {}
        if data is None: return

        try:
            with open(data) as f:
                data = f.read()
                self.meta['filename'] = f.name
        except (FileNotFoundError, TypeError, OSError):
            pass

        self.read_mdp(data)



    def __getattr__(self, name):
        """Return the value of the setting, or a Default object"""
        name, mdp_entry = self._get_mdp_entry(name)
        try:
            return self.values[name]
        except KeyError:
            if mdp_entry['default']:
                return Default(mdp_entry['default'])
            elif mdp_entry['options']:
                return Default(list(mdp_entry['options'])[0])
            else:
                return Default(None)

    @staticmethod
    def _check_number(value):
        if isinstance(value, str):
            itervalues = []
            value_to_save = value
            for val in value.split():
                try:
                    itervalues.append(int(val))
                except ValueError:
                    try:
                        itervalues.append(float(val))
                    except ValueError:
                        pass
        else:
            try:
                itervalues = list(value)
            except TypeError:
                itervalues = [value]
                value_to_save = value
            else:
                value_to_save = ' '.join(str(x) for x in value)
        for val in itervalues:
            if not isinstance(val, Number):
                raise AttributeError("{}.{} should be a Number expressed in units of {}, not {}".format(self, name, mdp_entry['units'], repr(value)))
        return value_to_save

    def __setattr__(self, name, value):
        """Set the value of the setting, as long as its appropriate"""
        name, mdp_entry = self._get_mdp_entry(name)

        if mdp_entry['units']:
            # If the mdp_entry has units specified, then its a number,
            # which means any options set are probably ranges. We should
            # make sure value is numeric and then just set it.
            self.values[name] = self._check_number(value)
        elif mdp_entry['options']:
            if process_name(value) in (process_name(v) for v in mdp_entry['options']):
                self.values[name] = value
            else:
                raise AttributeError("Acceptable values for {}.{} are {}, not '{}'".format(self, name, list(mdp_entry['options']), value))
        else:
            self.values[name] = value

    def __delattr__(self, name):
        """Return a setting to its default value"""
        name, mdp_entry = self._get_mdp_entry(name)

        try:
            del self.values[name]
        except KeyError:
            pass

    def __dir__(self):
        out = super().__dir__()
        opts = self.options.items()
        out += [process_name(v['docname']) for k,v in opts if k not in self.obsoletes]
        return out

    @classmethod
    def _get_mdp_entry(cls, name):
        """Return the options dict entry for name"""
        try:
            options_dict = cls.options[process_name(name)]
        except KeyError:
            raise AttributeError("{} object has no attribute '{}'".format(cls, name))
        else:
            return process_name(options_dict['docname']), options_dict

    @classmethod
    def help(cls, option, value=None):
        """Return the docstring for the option, or its value"""
        option, mdp_entry = cls._get_mdp_entry(option)
        if value is None:
            out = mdp_entry['docstring']
            if mdp_entry['options']:
                out += "\nPossible values are:"
                for k,v in mdp_entry['options'].items():
                    out += "\n    {}:".format(k)
                    out += "\n        {}".format(v)
        else:
            out = mdp_entry['options'][value]
        return out

    def write(self, tofile=None):
        """Return all the changed settings and comments in mdp format"""
        out = []
        keypad = 0
        for k in self.options:
            if process_name(k) in self.obsoletes:
                continue
            key, mdp_entry = self._get_mdp_entry(k)
            try:
                out.append((key, self.values[key], mdp_entry['docname']))
            except KeyError:
                pass
            else:
                keypad = max(keypad, len(mdp_entry['docname']))

        outstr = ""
        fmt = "{{: <{}}} = {{}}".format(keypad)
        for key,value,docname in out:
                comments = self.get_comments(key)
                if 'before' in comments:
                    outstr += '; ' + str(comments['before']).replace('\n', '\n; ') + "\n"
                outstr += fmt.format(docname, value)
                if 'on' in comments:
                    outstr += ' ; ' + str(comments['on']).replace('\n', '\n; ')
                outstr += "\n"
                if 'after' in comments:
                    outstr += '; ' + str(comments['after']).replace('\n', '\n; ') + "\n"

        comments = self.get_comments(tail_comment)
        if comments:
            outstr += '; ' + str(comments['after']).replace('\n', '\n; ')

        outstr = outstr.replace('; \n', '\n').strip() + '\n'
        if tofile:
            with open(tofile, 'w') as f:
                print(outstr, file=f, flush=True)
        else:
            return outstr

    def comment(self, option, comment, placement="on"):
        """Add a comment to the option"""
        if option == tail_comment:
            placement = 'after'
        else:
            option, _ = self._get_mdp_entry(option)

        try:
            comment_dict = self.comments[option]
        except KeyError:
            comment_dict = {}
            self.comments[option] = comment_dict

        if placement not in ['before', 'on', 'after']:
            raise ValueError("placement should be 'before', 'on' or 'after'")

        comment_dict[placement] = str(comment)

    def get_comments(self, option):
        if option != tail_comment:
            option, _ = self._get_mdp_entry(option)

        try:
            return self.comments[option]
        except KeyError:
            return {}

    _re_mdp_line_match = re.compile(r'^\s*(?P<key>\S*?)\s*=\s*(?P<value>.*?)\s*$')

    def read_mdp(self, mdp):
        try:
            values = deepcopy(mdp.values)
            comments = deepcopy(mdp.comments)
        except AttributeError:
            hangover_comment = ''
            values = {}
            comments = {}
            for line in mdp.splitlines():
                try:
                    command, comment = line.split(';', maxsplit=1)
                except ValueError:
                    command, comment = line, ''
                else:
                    if comment.startswith(' '):
                        comment = comment[1:]
                    if comment == '':
                        comment = ' '
                match = self._re_mdp_line_match.match(command)
                if match:
                    key = match.group('key')
                    value = match.group('value')
                    key, _ = self._get_mdp_entry(key)
                    values[key] = value
                    if (comment or hangover_comment) and key not in comments:
                        comments[key] = {}
                    if comment:
                        comments[key]['on'] = comment
                    if hangover_comment:
                        comments[key]['before'] = hangover_comment.replace('\n','',1)
                        hangover_comment = ''
                elif command:
                    raise ValueError("Uncommented text not matched: {}".format(repr(line)))
                else:
                    hangover_comment += '\n' + comment
            if hangover_comment:
                try:
                    comments[tail_comment]['after'] = hangover_comment.replace('\n','',1)
                except KeyError:
                    comments[tail_comment] = {'after': hangover_comment.replace('\n','',1)}
                except NameError:
                    raise ValueError("Can't handle a comment without an actual MDP entry")

        for key,value in values.items():
            self.__setattr__(key, value)
        for key,comment_dict in comments.items():
            for placement,comment in comment_dict.items():
                self.comment(key, comment, placement)

    def copy(self):
        return self.__class__(self.write())

    @classmethod
    def add_obsolete(cls, name, newname=None, docstring='Obsolete'):
        if newname:
            _, obs_dict = cls._get_mdp_entry(newname)
        else:
            obs_dict = {
                "docname": name,
                "docstring": docstring,
                "options": OrderedDict(),
                "default": 'Obsolete',
                "units": None
            }
        try:
            cls._get_mdp_entry(name)
        except AttributeError:
            name = process_name(name)
            cls.obsoletes.add(name)
            cls.options[name] = obs_dict
        else:
            raise ValueError('MDP option {} already exists.'.format(name))

    @classmethod
    def add_method(cls, func):
        setattr(cls, func.__name__, func)
        return func

    @classmethod
    def def_series(cls, pattern, start=1, stop=10):
        '''Define a series of values

        Pattern should be a regex string with a single capture group.
        All mdp options matching the regex pattern will have new copies
        added replacing the capture group with numbers taken from
        range(start, stop).'''
        pattern = re.compile(pattern)
        if pattern.groups != 1:
            raise ValueError('pattern should have exactly 1 capture group')
        new_options = {}
        for opt, opt_dict in cls.options.items():
            m = pattern.match(opt)
            if m:
                idx1, idx2 = m.span(1)
                bef_grp = m.string[:idx1]
                aft_grp = m.string[idx2:]
                for i in range(start, stop):
                    name = bef_grp + str(i) + aft_grp
                    docname_iter = iter(opt_dict['docname'])
                    docname = (
                        insert_breaks(docname_iter, bef_grp)
                        + str(i)
                        + insert_breaks(docname_iter, aft_grp)
                    )
                    new_dict = dict(opt_dict)
                    new_dict['docname'] = docname
                    new_options[name] = new_dict
        for k in new_options.keys():
            try:
                del cls.options[k]
            except KeyError:
                pass
        cls.options.update(new_options)

def insert_breaks(source, dest):
    """Insert - and _ from source into the same position in dest"""
    dest = iter(dest)
    out = ''
    for c in source:
        if c in '-_':
            out += c
        else:
            try:
                out += next(dest)
            except StopIteration:
                break
    return out


gmxdocs_path = os.path.dirname(__file__) + "/gmxdocs"

MDP20181 = MDPBase.create_subclass_from_doc(
    'MDP20181',
    f'{gmxdocs_path}/2018.1/mdp-options.rst'
)

MDP20183 = MDPBase.create_subclass_from_doc(
    'MDP20183',
    f'{gmxdocs_path}/2018.3/mdp-options.rst'
)

gmx_2018_series = [
    MDP20181,
    MDP20183
]

for cls in gmx_2018_series:
    cls.add_obsolete('title')
    cls.add_obsolete('nstxtcout', 'nstxout-compressed')
    cls.add_obsolete('xtc-precision', 'compressed-x-precision')
    cls.add_obsolete('xtc-grps', 'compressed-x-grps')

    cls.def_series('pullgroup(1)', start=1, stop=20)
    cls.def_series('pullcoord(1)', start=1, stop=10)
    cls.def_series('awh(1)', start=1, stop=10)
    cls.def_series('awh[1-9]+dim(1)', start=1, stop=5)
    cls.def_series('rot[a-z]+(0)', start=1, stop=10)

    MDP20183.options['awh1sharegroup']['options'] = OrderedDict()

    @cls.add_method
    def remove_pcouple(self):
        self.pcoupl = 'no'
        del self.pcoupltype
        del self.nstpcouple
        del self.tau_p
        del self.compressibility
        del self.ref_p
        self.refcoord_scaling = 'no'
        return self

    @cls.add_method
    def align_temps(self):
        if isinstance(self.tc_grps, Default):
            splits = [['🏳️‍🌈', self.ref_t, self.tau_t]]
        else:
            splits = zip(
                self.tc_grps.split(),
                self.ref_t.split(),
                self.tau_t.split()
            )

        maxlens = ((max(len(s) for s in tup), tup) for tup in splits)

        paddeds = ([s.ljust(l) for s in tup] for l, tup in maxlens)

        (
            self.tc_grps,
            self.ref_t,
            self.tau_t
        ) = ('  '.join(t).strip() for t in zip(*paddeds))

        if self.tc_grps == '🏳️‍🌈':
            del self.tc_grps

        return self


    @cls.add_method
    def set_temperature(self, temp):
        if self.annealing and self.annealing.lower() != 'no':
            raise ValueError(
                "This is an annealing simulation, set temps by hand"
            )
        if self.temperature_lambdas:
            raise ValueError(
                "Option 'temperature-lambdas' is set, set temps by hand"
            )
        if self.mc_temperature:
            raise ValueError(
                "Option 'mc-temperature' is set, set temps by hand"
            )
        if self.simulated_tempering and self.simulated_tempering.lower() != 'no':
            raise ValueError(
                "This is an simulated tempering simulation, set temps by hand"
            )

        if isinstance(temp, str):
            temp_str = temp
        else:
            temp_str = f'{temp:.2f}'

        if isinstance(self.tc_grps, Default):
            num_groups = 1
        else:
            num_groups = len(self.tc_grps.split())
        self.ref_t = ' '.join([temp_str] * num_groups)
        self.align_temps()

        if self.gen_vel.lower() == 'yes' or self.gen_temp:
            self.gen_temp = temp_str

        return self


    @cls.add_method
    def comment_total_time(self):
        total_ns = self.nsteps * self.dt / 1000
        self.comment('nsteps', "{} ns".format(total_ns), 'on')


    @cls.add_method
    def set_time(self, dt_fs, time_ns):
        dt_ps = dt_fs / 1000
        nsteps = int(time_ns / dt_ps * 1000)

        output_freqs = [
            ('nstxout', self.nstxout),
            ('nstvout', self.nstvout),
            ('nstfout', self.nstfout),
            ('nstlog', self.nstlog),
            ('nstcalcenergy', self.nstcalcenergy),
            ('nstenergy', self.nstenergy),
            ('nstxout_compressed', self.nstxout_compressed)
        ]
        for s, o in output_freqs:
            o = int(Default.unwrap(o))
            if o and nsteps % o:
                raise ValueError(f'nsteps={nsteps} is not divisible by {s}={o}')

        self.dt = dt_ps
        self.nsteps = nsteps
        self.comment_total_time()


    @cls.add_method
    def set_output_freq(self, freq_ps):
        dt_ps = float(self.dt)
        nst_output = int(freq_ps / dt_ps)
        nsteps = int(self.nsteps)

        if nsteps % nst_output:
            raise ValueError(
                f'nsteps={nsteps} is not divisible by {nst_output}'
            )

        print(f"Setting nstlog, nstenergy and nstxout_compressed to {nst_output}")

        self.nstlog = nst_output
        self.comment(
            'nstlog',
            f'{freq_ps} ps ; Steps between writing to .log file'
        )
        self.nstenergy = nst_output
        self.comment(
            'nstenergy',
            f'{freq_ps} ps ; Steps between writing energies to .edr file'
        )
        self.nstxout_compressed = nst_output
        self.comment(
            'nstxout_compressed',
            f'{freq_ps} ps ; Steps between writing compressed positions and velocities to .xtc file'
        )


