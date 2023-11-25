#!/usr/bin/env python3 
# -*- coding: utf-8 -*-

#--------1---------2---------3---------4---------5---------6---------7---------

"""Module for basic data treatment of disordered materials data.



Module: d4treat.py



This module contains functions and variables that are useful for the basic 

data treatment of disordered material diffractograms at D4.



The functions are thougth to use D4 data, but in fact they could be used for

a more general case.



For a listing of the functions in the module:

    dir(d4treat)



For a detailed help of all functions:

    help(d4treat)



Created on Wed Dec 30 23:47:15 2020

@author: Gabriel Cuello

---------------------------------------

"""

print('Module d4treat (by Gabriel Cuello, ILL, Jan 2022)')

print('Imported...')



#--------1---------2---------3---------4---------5---------6---------7---------8

import sys

import numpy as np                   # NumPy library

from scipy import integrate          # For integration

from .PhysicalMagnitudes import *

from .ExperimentalSettings import *

from .MiscelaneousCalculations import *

from .Element import *

from .FittingModels import *

from .FourierTransforms import *

from .InputOutput import *

#--------1---------2---------3---------4---------5---------6---------7---------8
