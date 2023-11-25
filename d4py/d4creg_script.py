#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 18:00:00 2023.

@author: cuello, robledo
"""

# Required modules (in the working directory):
#   - readD4.py
#   - D4treat.py
# Required files (in the working directory):
#   - effd4c.eff (or other efficiency file)
#   - dec.dec (or other shift file)

import sys,os

#/MyHome/python/modules
# path_modules = '/nethome/dif/cuello/python/modules/'
path_modules = 'modules/'
if not (path_modules in sys.path):
    sys.path.append(path_modules)

# Defines the working directory
#/MyHome/exp/D4-6-06-419/d4py/
#/home/cuello/MyHome/exp/D4-6-06-419/d4py/
#/net4/serdon/illdata/163/d4/exp_6-03-419/processed/D4pdf
#/nethome/dif/cuello/exp/D4-6-03-419/d4py/

CurrentWorkingDirectory = os.getcwd()
WorkingDirectory = "/net4/serdon/illdata/233/d4/exp_6-05-1067/processed/oski"
if WorkingDirectory != CurrentWorkingDirectory:
    WorkingDirectory = CurrentWorkingDirectory
    os.chdir(WorkingDirectory)
print(60*"-") 
print ("The working directory is: {}".format(os.getcwd()))
print(60*"-") 

# Import the regrouping script from the module readD4
from readD4 import d4creg, readParam

# Agregar opcion para cargar desde un archivo de parametros.
# para correr desde la linea de comandos

# --------1---------2---------3---------4---------5---------6---------7---------
print()
print(80 * "-")
print("Start of d4creg.py")

# paraFile = "pueba.par"
paraFile = sys.argv[1]
#print(sys.argv)
#paraFile = sys.argv[1]

runInfo = {}
runInfo = readParam(paraFile)
# --------1---------2---------3---------4---------5---------6---------7---------
# --------1---------2---------3---------4---------5---------6---------7---------
print()
print(80 * "-")
print()
print("End of d4creg.py")
print()
# --------1---------2---------3---------4---------5---------6---------7---------
