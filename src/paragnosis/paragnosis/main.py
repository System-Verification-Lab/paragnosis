#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys

from sympy import O
from paragnosis.misc import *
from paragnosis.bdd import *
from paragnosis.experiments import *
import multiprocessing
from .defaults import defaults
from .settings import Settings
from .version import version

def process(settings):

    # print settings
    if settings.verbose:
        settings.print()

    if not settings.has('location'):
        sys.stderr.write("FATAL: location of ParaGnosis toolset unknown. Please use --source to specify it or define the 'location' variable in your ~/.paragnosisrc.\n")
        sys.exit(1)
    elif settings.list:
        list_bayesian_networks(settings)
        sys.exit(0)
    elif settings.version:
        print("Version {}".format(version))
        sys.exit(0)
    elif settings.command == None:
        settings.print_help()
        sys.exit(1)
    else:
        if settings.command == "inference":
            compare_inference(settings)
        elif settings.command == "encode":
            compare_encoding(settings)
        else:
            compare_compilation(settings)


def main():
    MAX_CORES = multiprocessing.cpu_count()
    CORES     = [2**exp for exp in range(0,10) if 2**exp <= MAX_CORES]

    # load the settings from default and from
    settings = Settings()
    settings.update(defaults)
    settings.from_rc()

    parser = argparse.ArgumentParser(description='WMC testing suite',add_help=False)

    parser.add_argument('--help',       action='help',                             help='Show this help message and exit')
    parser.add_argument('--version',    action='store_true',                       help='Print version and exit')
    parser.add_argument('--dry-run',    action='store_true',                       help='Print the commands instead of running them')
    parser.add_argument('--debug',      dest='verbose',      action='store_true',  help='Print debug info')
    parser.add_argument('--list',       action='store_true',                       help='Print list of available Bayesian networks')
    parser.add_argument('--output-dir', dest='output_dir',   metavar="#DIRECTORY", help='Set the output directory (current: {})'.format(settings.output_dir), type=str)
    parser.add_argument('--source',     dest='location',     metavar="#DIRECTORY", help='Set the source directory of the Paragnosis toolset (current: {})'.format(settings.get("location", "unknown")), type=str, default=settings.get("location"))

    bdd_options = ["wpbdd","parallel-wpbdd","pwpbdd","parallel-pwpbdd","mg","tdmg"]
    bdd_default = ["tdmg"]

    subparsers = parser.add_subparsers(dest='command')

    # encoding parser
    encoding = subparsers.add_parser('encoding', help='Create Bayesian network encodings')
    encoding.add_argument(dest='network', type=str,  metavar="NETWORK", help="Provide a hugin file or a Bayesian network name provided by ./pg --list")
    encoding.add_argument('--encoding-help',dest='help',action='store_true',help='Show help for encoding')
    encoding.add_argument('--args',dest='args',nargs=argparse.REMAINDER,help="Encoding arguments",metavar='ARG')

    # compilation parser
    compilation = subparsers.add_parser('compile', help='Compile Bayesian networks')
    compilation.add_argument(                dest='network',                                         metavar="NETWORK",     type=str,             help="Provide a hugin file or a Bayesian network name provided by ./pg --list")
    compilation.add_argument('--cores',      dest='cores',      required=False, default=CORES,       metavar="#CORES",      type=int,  nargs='+', help='Number of cores to use during parallel execution', )
    compilation.add_argument('--method',     dest='bdds',       required=False, default=bdd_default, metavar='BDD',         nargs='+',            help='Type of BDD. Options are ' + ', '.join(bdd_options), choices=bdd_options)
    compilation.add_argument('--partitions', dest='partitions', required=False, default=2,           metavar='#PARTITIONS', type=int,             help='Set number of partitions')
    compilation.add_argument('--overwrite',  dest='overwrite',  required=False,                                             action='store_true',  help='Overwrite ordering, partitioning, etc.')
    compilation.add_argument('--repeat',     dest='repeat',     required=False, default=1,                                  type=int,             help='Set number of compilation repeats')

    # inference parser
    inference = subparsers.add_parser('inference', help='Perform inference')
    inference.add_argument(                dest='network',                                         type=str,  metavar="NETWORK",               help="Provide a hugin file or a Bayesian network name provided by ./pg --list")
    inference.add_argument('--evidence',   dest='evidence',   required=False, default=None,        type=str,  metavar='#EVIDENCE',             help='Provide evidence in the form \'var=value[, var2=value2][, var3=value3]\'')
    inference.add_argument('--posteriors', dest='posteriors', required=False, default=None,        type=str,  metavar='#VARIABLES',            help='Provide variables in the form \'var[, var2][, var3]\'', nargs='?')
    inference.add_argument('--verify',     dest='verify',     required=False,                                 action='store_true',             help='Verify inference answers')
    inference.add_argument('--compare',    dest='compare',    required=False, default="",          type=str,  metavar="#OBSERVED",  nargs='?', help='Limit number of observerd variables during inference', const="")
    inference.add_argument('--cases',      dest='cases',      required=False,                      type=str,  metavar='#FILE',                 help='Provide a file with inference test cases\'')
    inference.add_argument('--method',     dest='bdds',       required=False, default=bdd_default,            metavar='BDD',        nargs='+', help='Type of BDD. Options are ' + ', '.join(bdd_options), choices=bdd_options)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)

    # parse provided args
    options = parser.parse_args()

    # update the settings
    settings.update(options.__dict__)
    settings.print_help = lambda : parser.print_help()

    # handle provided input
    process(settings)

#if __name__ == "__main__":
#   main(sys.argv[1:])

