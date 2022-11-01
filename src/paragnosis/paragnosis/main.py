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
import collections.abc

def process(settings):
    if settings.list == True:
        list_bayesian_networks()
        sys.exit(0)
    else:
        if settings.test != None:
            if settings.test == "inference":
                compare_inference(settings)
            elif settings.test == "compilation":
                compare_compilation(settings)
            else:
                compare_encoding(settings)


def entry_to_dict(key, value) -> dict:
    dest = key.split('.',1)
    if len(dest) == 2:
        return { dest[0]: entry_to_dict(dest[1], value) }
    else:
        return { key: value }

def update(d: dict, mapping) -> dict:
    for key, value in mapping.items():
        if isinstance(value, collections.abc.Mapping):
            dd = d.get(key, {})
            if not isinstance(dd, dict):
                dd = {}

            d[key] = update(dd, value)
        else:
            d[key] = value

    return d

def to_dict(parsed):
    dictionary = {}
    for key, value in vars(parsed).items():
        entry = entry_to_dict(key,value)
        update(dictionary,entry)

    return dictionary


def main():
    MAX_CORES = multiprocessing.cpu_count()
    CORES = [2**exp for exp in range(0,10) if 2**exp <= MAX_CORES]

    parser = argparse.ArgumentParser(description='WMC testing suite',add_help=False)

    parser.add_argument('--help',action='help',help='Show this help message and exit')
    parser.add_argument('--list', dest='list', action='store_true', help='Print list of available Bayesian networks', required=False)
    parser.add_argument('--cores', dest='cores', nargs='+', help='Number of cores to use during parallel execution', required=False, default=CORES,metavar="#CORES",type=int)
    test_options = ["compilation","inference","encoding"]
    parser.add_argument('--test',dest='test',choices=test_options, help='Choose what to test. Options are ' + ', '.join(test_options), required=False,metavar='TEST')
    parser.add_argument('--network',dest='networks',nargs='+', help='Bayesian network(s) used for testing',metavar="#NETWORK")

    bdd_options = ["wpbdd","parallel-wpbdd","pwpbdd","parallel-pwpbdd","sdd","sddr","obdd","zbdd","lazy","shafershanoy","ve","dlib","ace","acei","mg","tdmg"]
    group = parser.add_argument_group('inference and compilation arguments')
    group.add_argument('--bdd',dest='bdds',nargs='+', help='Type of BDD. Options are ' + ', '.join(bdd_options), choices=bdd_options,required=False,default=None,metavar='BDD')
    group.add_argument('--partitions',dest='partitions',help='Set number of partitions',default=2,type=int,metavar='#PARTITIONS')
    group.add_argument('--overwrite',dest='overwrite',action='store_true', help='Overwrite ordering, partitioning, etc.')
    group.add_argument('--verify', dest='verify', action='store_true', help='Verify inference answers', required=False)
    group.add_argument('--compare', dest='compare', type=str, nargs='?', metavar="#OBSERVED", help='Limit number of observerd variables during inference', required=False, const="", default="")
    group.add_argument('--evidence',dest='evidence', type=str, help='Provide evidence in the form \'var=value[, var2=value2][, var3=value3]\'', required=False, metavar='#EVIDENCE')
    group.add_argument('--posteriors',dest='posteriors', type=str, help='Provide variables in the form \'var[, var2][, var3]\'', nargs='?', required=False, metavar='#VARIABLES')
    group.add_argument('--cases',dest='cases', type=str, help='Provide a file with inference test cases\'', required=False, metavar='#FILE')

    group.add_argument('--verbose', dest='verbose', action='store_true', help=argparse.SUPPRESS, required=False)
    group.add_argument('--repeat',dest='repeat',help='Set number of compilation repeats',default=3,type=int, required=False)

    group = parser.add_argument_group('encoding arguments')
    group.add_argument('--encoding-help',dest='help',action='store_true',help='Show help for encoding')
    group.add_argument('--args',dest='args',nargs=argparse.REMAINDER,help="Encoding arguments",metavar='ARG')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)

    options = parser.parse_args()

    if options.test == 'compilation' or options.test == 'inference':
        if options.networks == None or options.bdds == None:
            parser.error('{} requires --network and --bdd'.format(options.test))

    if options.test == 'encoding':
        if not options.help and (options.networks == None):
            parser.error('{} requires --network'.format(options.test))

    # get all settings
    settings = Settings()
    settings.update(defaults)
    settings.from_rc()
    settings.update(to_dict(options))
     
    process(settings)

#if __name__ == "__main__":
#   main(sys.argv[1:])

