from collections import OrderedDict, Mapping
from .rstparser import parse_rst_file
import os
import re
from numbers import Number



class Default():
    """Class to signal that a setting hasn't been specified"""
    def __init__(self, object):
        self.obj = object

    def __repr__(self):
        return "Default({})".format(repr(self.obj))

class MDPBase():
    """Base class for creating MDP files"""
    @classmethod
    def create_subclass_from_doc(cls, name, docpath):
        """Create an MDP class from a mdp-options.rst file distributed with GROMACS"""
        mdp_entries, doc = parse_rst_file(docpath)
        attr = {
            'options': mdp_entries,
            'doc': doc
        }
        return type(name, (cls,), attr)

    def __init__(self, data=None):
        super().__init__()
        self.__dict__['values'] = {}
        self.__dict__['comments'] = {}
        if data is None: return

        try: 
            with open(data) as f:
                data = f.read()
        except (FileNotFoundError, TypeError, OSError):
            pass
        self.read_mdp(data)



    def __getattr__(self, name):
        """Return the value of the setting, or a Default object"""
        name, mdp_entry = self._get_mdp_entry(name)
        try:
            return self.values[name]
        except KeyError:
            return Default(mdp_entry['default'])

    def __setattr__(self, name, value):
        """Set the value of the setting, as long as its appropriate"""
        name, mdp_entry = self._get_mdp_entry(name)

        if mdp_entry['options']:
            if value.lower() in mdp_entry['options']:
                self.values[name] = value
            else: 
                raise AttributeError("Acceptable values for {}.{} are {}, not '{}'".format(self, name, list(mdp_entry['options']), value))
        elif mdp_entry['units']:
            if isinstance(value, Number):
                self.values[name] = value
            else:
                raise AttributeError("{}.{} should be a Number expressed in units of {}".format(self, name, mdp_entry['units']))
        else:
            self.values[name] = value

    def __delattr__(self, name):
        """Return a setting to its default value"""
        name, mdp_entry = self._get_mdp_entry(name)
        
        try:
            del self.values[name]
        except KeyError:
            pass

    @classmethod
    def _get_mdp_entry(cls, name):
        """Return the options dict entry for name"""
        name = name.replace('_', '-')
        try:
            return name, cls.options[name]
        except KeyError:
            raise AttributeError("{} object has no attribute '{}'".format(cls, name))


    def help(self, option, value=None):
        """Return the docstring for the option, or its value"""
        option, mdp_entry = self._get_mdp_entry(option)
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

    def write(self):
        """Return all the changed settings and comments in mdp format"""
        out = []
        keypad = 0
        for key in self.options:
            try:
                out.append((key, self.values[key]))
            except KeyError:
                pass
            else:
                keypad = max(keypad, len(key))

        outstr = ""
        fmt = "{{: <{}}} = {{}}".format(keypad)
        for key,value in out:
                comments = self.get_comments(key) 
                if 'before' in comments:
                    outstr += '; ' + str(comments['before']).replace('\n', '\n; ') + "\n"
                outstr += fmt.format(key, value)
                if 'on' in comments:
                    outstr += ' ; ' + str(comments['on']).replace('\n', '\n; ')
                outstr += "\n"
                if 'after' in comments:
                    outstr += '; ' + str(comments['after']).replace('\n', '\n; ') + "\n"
        outstr = outstr.replace('; \n', '\n')
        return outstr

    def comment(self, option, comment, placement="on"):
        """Add a comment to the option"""
        option, mdp_entry = self._get_mdp_entry(option)

        try:
            comment_dict = self.comments[option]
        except KeyError:
            comment_dict = {}
            self.comments[option] = comment_dict

        if placement not in ['before', 'on', 'after']:
            raise ValueError("placement should be 'before', 'on' or 'after'")

        comment_dict[placement] = str(comment)

    def get_comments(self, option):
        option, mdp_entry = self._get_mdp_entry(option)

        try:
            return self.comments[option]
        except KeyError:
            return {}

    _re_mdp_line_match = re.compile(r'\s*(?P<key>\S+?)\s*=\s*(?P<value>\S+)')

    def read_mdp(self, mdp):
        try:
            values = mdp.values
            comments = mdp.comments
        except AttributeError:
            hangover_comment = ''
            values = {}
            comments = {}
            for line in mdp:
                try:
                    command, comment = line.split(';', maxsplit=1)
                except ValueError:
                    command, comment = line, ''
                else:
                    comment += ' '
                match = self._re_mdp_line_match.match(command)
                if match:
                    key = match.group('key')
                    value = match.group('value')
                    values[key] = value
                    if (comment or hangover_comment) and key not in comments:
                        comments[key] = {}
                    if comment: 
                        comments[key]['on'] = comment
                    if hangover_comment: 
                        comments[key]['before'] = hangover_comment.replace('\n','',1)
                        hangover_comment = ''
                elif command:
                    raise ValueError("Uncommented text not matched!")
                else:
                    hangover_comment += '\n' + comment
            if hangover_comment:
                try:
                    comments[key]['after'] = hangover_comment.replace('\n','',1)
                except KeyError:
                    comments[key] = {'after': hangover_comment.replace('\n','',1)}
                except NameError:
                    raise ValueError("Can't handle a comment without an actual MDP entry")

        for key,value in values.items():
            self.__setattr__(key, value)
        for key,comment_dict in comments.items():
            for placement,comment in comment_dict.items():
                self.comment(key, comment, placement)





gmxdocs_path = os.path.dirname(__file__) + "/gmxdocs"

MDP20181 = MDPBase.create_subclass_from_doc('MDP20181', '{}/2018.1/mdp-options.rst'.format(gmxdocs_path))
