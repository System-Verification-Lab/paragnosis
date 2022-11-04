#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from . import util
from .settings import Settings

defaults = Settings()
defaults.package_dir = util.get_file_directory(__file__)
defaults.output_dir  = os.path.join(os.getcwd(), "output")