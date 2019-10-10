#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches
import numpy as np

from datetime import datetime, timedelta

from utilities import fio, plotting

df, dfa = fio.read_AWS_wagerup()

