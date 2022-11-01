#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
import subprocess
import os
import sys
import re
import platform
from distutils.util import strtobool

def get_platform():
    pf = platform.system()

    if re.compile(r'darwin', re.IGNORECASE).search(pf):
        return "osx"
    elif re.compile(r'win', re.IGNORECASE).search(pf):
        return "windows"
    elif re.compile(r'linux', re.IGNORECASE).search(pf):
        return "linux"
    elif re.compile(r'cygwin', re.IGNORECASE).search(pf):
        return "linux"
    else:
        return  pf

def select(lst: list, pattern: str, nr: int = 0):
    r = re.compile(pattern)
    results = list(filter(r.match, map(str, lst)))
    if nr == 0:
        return results
    elif len(results) == nr:
        if nr == 1:
            return results[0]
        else:
            return results 
    else:
        raise RuntimeError("Expected to find {} result(s) instead of {} using pattern {} in list:\n{}"
            .format(nr, len(results), pattern, "\n    ".join(lst)))

def call(command, verbose: bool, dryrun: bool, cwd = None, stdout = sys.stdout.fileno(), stderr = sys.stderr.fileno(), shell = False):
    if verbose or dryrun:
        if not shell:
            print_command = map(lambda w: w if ' ' not in w else "'" + w + "'", command.copy())
            print_command = " ".join(print_command)
        else: 
            print_command = command

        if cwd != None and cwd != ".": 
            print("{}>".format(cwd),print_command)
        else:
            print("> ",print_command)

    if dryrun:
        return 0
    else:
        return subprocess.call(command, cwd=cwd, stdout = stdout, stderr = stderr, shell = shell)

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def strip_margin(text):
    def remove_margin(line):
        return re.sub(r'^[ ]*\|', '', line)

    lines = [ remove_margin(line) for line in text.splitlines() ]
    return '\n'.join(lines)

def get_file_directory(path, *dirs):
    path = os.path.dirname(os.path.realpath(path))
    for directory in dirs:
        path = os.path.join(path, directory)

    return os.path.abspath(path)

def get_abs_dir(path, *dirs):
    path = os.path.realpath(path)
    for directory in dirs:
        path = os.path.join(path, directory)

    return os.path.abspath(path)

def find_ancestor_directory(directory):
    cdir = os.path.abspath(os.path.join(SCRIPT_DIR,".."))
    while not os.path.exists(os.path.join(cdir,directory)) and cdir != '/':
        cdir = os.path.abspath(os.path.join(cdir,".."))

    if cdir == '/':
        raise Exception("Could not find directory: {:s}".format(directory))
    else:
        return cdir

def require_dir(path, msg = ""):
    if not os.path.isdir(path):
        raise RuntimeError("Required directory '" + path + "' does not exist. " + msg)

def require_file(path, msg = ""):
    if not os.path.isfile(path):
        raise RuntimeError("Required file '" + path + "' does not exist. " + msg)

def require(file):
    try:
        os.stat(file)
    except:
        raise RuntimeError("required file '{:s}' not found".format(file))

def query_yes_no(question, default='no'):
    if default is None:
        prompt = " [y/n] "
    elif default == 'yes':
        prompt = " [Y/n] "
    elif default == 'no':
        prompt = " [y/N] "
    else:
        raise ValueError(f"Unknown setting '{default}' for default.")

    while True:
        try:
            resp = input(question + prompt).strip().lower()
            if default is not None and resp == '':
                return default == 'yes'
            else:
                return strtobool(resp)
        except ValueError:
            print("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

