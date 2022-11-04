#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from . import util
from  .util import eprint
import configparser
import collections.abc

class Settings():

    def __init__(self):
        self.__dict__ = {}

    def print(self, indent = 0, fileptr=sys.stdout):
        if indent == 0:
            print("\nSettings:")
            self.print(1)
            print("\n")
        else:
 
            settings = {k: v for k, v in self.__dict__.items() if isinstance(v, Settings)}
            values = {k: v for k, v in self.__dict__.items() if not isinstance(v, Settings)}

            for key, value in sorted(settings.items()):
                print((" " * (indent*4)) + "[" + key + "]", file=fileptr)
                value.print(indent + 1, fileptr)

            if settings:
                print("")

            for key, value in sorted(values.items()):
                print((" " * (indent*4)) + "{:s} = {}".format(key, value), file=fileptr)

        return 0

    def __resolve(self, unresolved, resolved):
        if len(unresolved) == 0:
            return None
        elif len(unresolved) == 1:
            attribute = unresolved[0]
            if self.has(attribute):
                if isinstance(attribute, Settings):
                    raise RuntimeError("Settings path error: " + ".".join(resolved) + "." + attribute + " is expected to be an attribute")
            
            return self.__dict__.get(attribute)
        else:
            attribute = unresolved[0]
            if self.has(attribute):
                settings = self.value(attribute)
                if isinstance(settings, Settings):
                    resolved.append(unresolved.pop(0))
                    return settings.__resolve(unresolved, resolved)
                else:
                    raise RuntimeError("Settings path error: " + ".".join(resolved) + "." + attribute + " is not suppose to be an attribute")
            else:
                return None

    def setu(self, path: str, value):
        val = self.get(path)
        if val == None:
            self.set(path, value)

    def set(self, path: str, value): 
        attributes = path.split(".")   
        if len(path) == 0:
            raise RuntimeError("Cannot add to settings without a key.")

        def to_dictionary(attr):
            if (len(attr)) == 1: 
                return { attr[0]: value }
            else:
                attribute =  attr.pop(0)
                return { attribute: to_dictionary(attr) }

        self.update(to_dictionary(attributes))

    def get(self, path: str):
        attributes = path.split(".")
        res = self.__resolve(attributes, [])
        return res

    def to(self, path: str):
        value = self.get(path)
        if value == None:
            self.print()
            raise RuntimeError("No value is set for setting: " + path)

        return value

    def to_dir(self, path: str):
        value = self.to(path)
        if os.path.isdir(value):
            return value
        else:
            raise RuntimeError("Setting " + path + " contains a directory that does not exist: " + value)

    def to_file(self, path: str):
        value = self.to(path)
        if os.path.isfile(value):
            return value
        else:
            raise RuntimeError("Setting " + path + " contains a file that does not exist: " + value)

    def value(self, key: str):
        return self.__dict__[key]

    def require(self, key: str):
        if key not in self.__dict__:
            eprint()
            raise RuntimeError("Key '{}' in settings is required to be set".format(key))

    def has(self, key: str):
        return key in self.__dict__

    def remove(self, key: str):
        del self.__dict__

    def items(self):
        return self.__dict__.items()

    def keys(self):
        return self.__dict__.keys()

    def update(self, mapping) -> 'Settings':
        for key, value in mapping.items():
            if isinstance(value, collections.abc.Mapping) or isinstance(value, Settings):
                settings = self.__dict__.get(key, Settings())
                if not isinstance(settings, Settings):
                    settings = Settings()
                self.__dict__[key] = settings.update(value)
            else:
                self.__dict__[key] = value

        return self

    def to_dict(self) -> dict:
        dictionary = self.__dict__
        for key, value in dictionary.items():
            if isinstance(value, Settings):
                dictionary[key] = value.to_dict()

        return dictionary

    def set_config(self):
        # set the config file and locations
        platform = util.get_platform()
        if platform == "windows":
            self.config_filename = "paragnosis.ini"
            try:
                windir = os.path.join(os.environ['APPDATA'],'paragnosis')
                if not os.path.exists(windir):
                    os.mkdir(windir)
                self.config_dirs = [ windir ]
            except:
                eprint("Could not set windows config dir to %APPDATA%/paragnosis.")
                self.config_dirs = [ ]
        else:
            self.config_filename = ".paragnosisrc"
            self.config_dirs = [ os.path.expanduser('~') ]

        if self.has('config_dir'):
            self.config_dirs.extend(self.value('config_dir'))

        self.config_filenames = [ os.path.join(dir, self.config_filename) for dir in self.config_dirs ]
        self.config_filenames.extend( [ self.config_filename ] )

    def from_rc(self):
        self.set_config()
        read_files = []

        parser = configparser.ConfigParser()
        for filename in self.config_filenames:
            try:
                with open(filename) as stream:
                    # insert default section (for section-less config files)
                    parser.read_string("[root]\n" + stream.read())
                    sections = dict(parser._sections)

                    # add root section to root of settings
                    self.update(sections['root'])
                    del sections['root']
                    self.update(sections)

                    read_files.append(filename)
            except FileNotFoundError as error: 
                pass

        return read_files

    def to_config(self):
        parser = configparser.ConfigParser()
        return self._to_config(parser)

    def _to_config(self, parser):
        for key, value in self.__dict__.items():
            if isinstance(value, Settings):
                parser[key] = value._to_config({})
            elif isinstance(parser, dict):
                parser[key] = value

        return parser

    def write_rc(self):
        self.set_config()
        parser = self.to_config()

        dir = self.to('config_dir')
        self.to("dry_run") or util.require_dir(dir)

        filename = os.path.join(dir,self.config_filename)
        if self.to("dry_run"):
            print("Writing config to: " + filename + " (dry-run)")
        else:
            print("Writing config to: " + filename)

        try:
            configfile = open(filename, 'w')
        except IOError:
            raise RuntimeError("Could not write to config file: {}".format(filename))
        else:
            with configfile:
                parser.write(configfile)
            
        return 0


