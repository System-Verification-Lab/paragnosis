#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
from paragnosis.version import version

setup(
    name='ParaGnosis',
    version=version,
    zip_safe=False,
    description='A tool for Weigted model counting.',
    keywords='ParaGnosis',
    packages=['paragnosis'],
    scripts = ['bin/pg'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Software Development',
        'Topic :: Software Development :: Build Tools',
        'Programming Language :: Python :: 3'
    ]
    )
