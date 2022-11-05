#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from . import util
from .settings import Settings
import multiprocessing

defaults = Settings()
defaults.package_dir = util.get_file_directory(__file__)
defaults.output_dir  = os.path.join(os.getcwd(), ".pg.output")
defaults.partitions  = 2
defaults.overwrite   = False
defaults.cores       = [multiprocessing.cpu_count()] # [2**exp for exp in range(0,10) if 2**exp <= MAX_CORES]
