#!/usr/bin/env python3

import paragnosis.misc as misc
import paragnosis.globals as g
from paragnosis.bdd import *

def compare_inference(settings):
    for network in settings.networks:
        bdd = Bdd(settings)
        bdd.set_cores(settings.cores)
        bdd.set_overwrite(settings.overwrite)
        bdd.set_verify(settings.verify)
        bdd.set_compare(settings.compare.split(','))
        bdd.set_evidence(settings.evidence)
        bdd.set_posteriors(settings.posteriors)
        bdd.set_bayesian_network(network)
        bdd.set_partitions(settings.partitions)
        bdd.set_verbose(settings.verbose)
        #bdd.set_timeout(5)
        bdd.run_inference(settings.bdds)
        #bdd.set_timeout(None)
        if not settings.evidence:
            bdd.print_inference_results()


def compare_compilation(settings):
    for network in settings.networks:
        bdd = Bdd(settings)
        bdd.set_cores(settings.cores)
        bdd.set_overwrite(settings.overwrite)
        bdd.set_bayesian_network(network)
        bdd.set_repeat(settings.repeat);
        bdd.set_verbose(settings.verbose)
        bdd.set_partitions(settings.partitions)
        bdd.run_compilation(settings.bdds)
        bdd.print_compilation_results()

def compare_encoding(settings):
    if settings.help:
        bdd = Bdd()
        bdd.help_encoding()
    else:
        for network in settings.networks:
            bdd = Bdd()
            bdd.set_bayesian_network(network)
            bdd.set_verbose(settings.verbose)
            bdd.run_encoding(settings.args)

