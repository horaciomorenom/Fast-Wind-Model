#!/usr/bin/env python3

import numpy as np
import matplotlib as plt
import hissw

ssw = hissw.Environment(ssw_packages=['packages/chianti/'], ssw_paths=['chianti'])
ssw.run('run_wind')


