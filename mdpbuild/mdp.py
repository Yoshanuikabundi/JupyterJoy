from collections import OrderedDict
from .rstparser import parse_rst_file
import os


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

    def __getattr__(self, name):
        """Return the value of the setting, or a Default object"""
        mdp_entry = self._get_mdp_entry(name)
        try:
            return mdp_entry['value']
        except KeyError:
            return Default(mdp_entry['default'])

    def __setattr__(self, name, value):
        """Set the value of the setting, as long as its appropriate"""
        mdp_entry = self._get_mdp_entry(name)

        if mdp_entry['options']:
            if value in mdp_entry['options']:
                mdp_entry['value'] = value
            else: 
                raise AttributeError("Acceptable values for {}.{} are {}, not '{}'".format(self, name, list(mdp_entry['options']), value))
        elif mdp_entry['units']:
            if isinstance(value, Number):
                mdp_entry['value'] = value
            else:
                raise AttributeError("{}.{} should be a Number expressed in units of {}".format(self, name, mdp_entry['units']))
        else:
            mdp_entry['value'] = value

    def __delattr__(self, name):
        """Return a setting to its default value"""
        mdp_entry = self._get_mdp_entry(name)
        
        try:
            del mdp_entry['value']
        except KeyError:
            pass

    def _get_mdp_entry(self, name):
        """Return the options dict entry for name"""
        try:
            return self.options[name]
        except KeyError:
            raise AttributeError("{} object '{}' has no attribute '{}'".format(self.__class__, self, name))


    def help(self, option, value=None):
        """Return the docstring for the option, or its value"""
        mdp_entry = self._get_mdp_entry(option)
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
        for key, entry in self.options.items():
            if 'value' in entry:
                out.append((key, entry['value']))
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
        return outstr

    def comment(self, option, comment, placement="on"):
        """Add a comment to the option"""
        mdp_entry = self._get_mdp_entry(option)

        try:
            comment_dict = mdp_entry['comments']
        except KeyError:
            comment_dict = {}
            mdp_entry['comments'] = comment_dict

        if placement not in ['before', 'on', 'after']:
            raise ValueError("placement should be 'before', 'on' or 'after'")

        comment_dict[placement] = str(comment)

    def get_comments(self, option):
        mdp_entry = self._get_mdp_entry(option)

        try:
            return mdp_entry['comments']
        except KeyError:
            return {}

gmxdocs_path = os.path.dirname(__file__) + "/gmxdocs"

MDP20181 = MDPBase.create_subclass_from_doc('MDP20181', '{}/2018.1/mdp-options.rst'.format(gmxdocs_path))
