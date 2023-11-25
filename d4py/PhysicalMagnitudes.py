#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------1---------2---------3---------4---------5---------6---------7---------
"""
Module for basic data treatment of disordered materials data.

Module: PhysicalMagnitudes.py

This module contains functions and variables that are useful for the basic 
data treatment of disordered material diffractograms at D4.

The functions are thougth to use D4 data, but in fact they could be used for
a more general case.

For a listing of the functions in the module:
    dir(PhysicalMagnitudes)

For a detailed help of all functions:
    help(PhysicalMagnitudes)

Author: Gabriel Cuello
Date: Dec 30, 2020
Modified: Feb 7, 2023
--------------------------------------------------------------------------------
"""
# Feb 7, 2023: I split the original d4treat.py in smaller modules regrouping
# functions by topics.

#--------1---------2---------3---------4---------5---------6---------7---------
import sys
import numpy as np                   # NumPy library
from scipy import integrate          # For integration
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------8
def getMassNumber(molarMass=6.0):
    """
Calculate the mass number (A).
    
    A = molar_mass / neutron_mass

Input:
  - molarMass: float, optional
        Molar mas in amu. The default is 6.0 (carbon molar mass).

Output: 
  - A : float
        The mass number

Author: Gabriel Cuello
Date: Jan 13, 2022
--------------------------------------------------------------------------------
    """
    neutronMass = 1.0086649 # in amu
    A = molarMass / neutronMass
    return A
#--------1---------2---------3---------4---------5---------6---------7---------8

#--------1---------2---------3---------4---------5---------6---------7---------8
def getFreeXS(BoundXS=5.08,A=50.9415):
    """
Calculate the free cross section.
    
The free cross section is calculated as:
        sigma_free = sigma_bound * A^2 /(1+A)^2
where sigma_bound is the bound cross section (in barns) and A is the mass 
number.

Input: BoundXS, A
  - BoundXS : float, optional
        Bound cross section. 
        The default is 5.08 barns, the vanadium bound cross section.
  - A : float, optional
        Mass number. 
        The default is 50.9415, the vanadium mass number.

Output:
  - FreeXS : float
        The free cross section (same units as the bound cross section)
    
Author: Gabriel Cuello
Date: Jan 13, 2022
--------------------------------------------------------------------------------
    """
    FreeXS = BoundXS * A**2 / (1.0+A)**2
    return FreeXS
#--------1---------2---------3---------4---------5---------6---------7---------8

#--------1---------2---------3---------4---------5---------6---------7---------8
def extractDict(dictionary,key):
    """
Extract a sub-dictionary from a bigger one.
    
This function extracts a transversal dictionary from a bigger dictionary 
containig dictionaries. The extracted dictionary will have the same keys that 
the bigger one, but the values will be only those corresponding to the input 
key.

Input: dictionary, key
  - dictionary : a dictionary
      The big dictionary containing sub-dictionaries.
  - key : string
      The key of the subdictionaries that will be extracted.

Output: 
  - A dictionary with the keys of the big one and containing only the elements 
    with the given key.

Example:
    For the following dictionary
        elements = {'Ti':{'CohXS': 1.485,'IncXS': 2.87,  'AbsXS':6.09 },
                    'Nb':{'CohXS': 6.253,'IncXS': 0.0024,'AbsXS':1.15 },
                    'P': {'CohXS': 3.307,'IncXS': 0.005, 'AbsXS':0.172}}
    The instruction: incohXS = extractDict(elements,'IncXS') will produce:
        incohXS = {'Ti':2.87,'Nb':0.0024,'P':0.005}

Author: Gabriel Cuello
Date: Oct 29, 2021
--------------------------------------------------------------------------------
    """
    result = {}
    for key1,value in dictionary.items():
        result[key1] = dictionary[key1][key]
    return result
#--------1---------2---------3---------4---------5---------6---------7---------8

#--------1---------2---------3---------4---------5---------6---------7---------8
def extractAttr(Dict,attr):
    """
Extract a an attribue from a dictionary containing classes.
    
This function extracts a transversal dictionary from a bigger dictionary 
containig dictionaries. The extracted dictionary will have the same keys that 
the bigger one, but the values will be only those corresponding to the input 
key.

Input: dictionary, key
  - dictionary : a dictionary
      The big dictionary containing sub-dictionaries.
  - key : string
      The key of the subdictionaries that will be extracted.

Output: 
  - A dictionary with the keys of the big one and containing only the elements 
    with the given key.

Example:
    For the following dictionary
        elements = {'Ti':{'CohXS': 1.485,'IncXS': 2.87,  'AbsXS':6.09 },
                    'Nb':{'CohXS': 6.253,'IncXS': 0.0024,'AbsXS':1.15 },
                    'P': {'CohXS': 3.307,'IncXS': 0.005, 'AbsXS':0.172}}
    The instruction: incohXS = extractDict(elements,'IncXS') will produce:
        incohXS = {'Ti':2.87,'Nb':0.0024,'P':0.005}

Author: Jose Robledo
Date: Sept 29, 2023
--------------------------------------------------------------------------------
    """
    result = {}
    for key, val in Dict.items():
        if hasattr(val, attr):
            result[key] = getattr(val, attr)
        else:
            print(f"Attribute {attr} not found in {val}.")
    return result
#--------1---------2---------3---------4---------5---------6---------7---------8

#--------1---------2---------3---------4---------5---------6---------7---------8
def getAtomicDensity(density=1.0,molarMass=6.0):
    """
Calculate the atomic density (atoms/volume).

This function returns the atomic density (at/A3) given the macroscopic density 
(in g/cm3) and the average atomic molar mass (g/mol).
    
    Atomic density [at/A3] = 
           density [g/cm3] * NA [at/mol] / molarM [g/mol] * 10^(-24)
where molarM is the average molar mass per atom and NA = 6.02214076 10^(23) 
at/mol is the Avogadro's number.

For example, a water molecule has 3 atoms and the average molar mass is 
18 g/mol/3 = 6 g/mol, then
    AtomicDensity = 1 g/cm3 * 0.602214076 at/mol / 6 g/mol = 0.1 at/A3
        
If density is less or equal to 0, the function returns water atomic density, 
i.e, 0.1 at/A3.
    
Input: density, molarMass
    density : float, optional
        The macroscopic density in g/cm3. The default is 1.0.
    molarMass : float, optional
        The average molar mass per atom, in g/mol/at. The default is 6.0.

Output: float
    AtomicDensity : The atomic density, in atoms/A3
    
Author: Gabriel Cuello
Date: Oct 25, 2021
--------------------------------------------------------------------------------
    """
# Avogadros's number times 10^(-24): 
#    6.02214076 10^(23) * 10^(-24) = 0.602214076
# NA = 0.602214076
    if density <= 0:
        density = 1.0
        molarMass = 6.0
        print ('Attention!')
        print ('    Using the water density in getAtomicDensity function.')
    return density * 0.602214076 / molarMass
#--------1---------2---------3---------4---------5---------6---------7---------8

#--------1---------2---------3---------4---------5---------6---------7---------8
def getDensity(atomic_density=0.1,molarM=6.0):
    """
Calculate the macroscopic density (g/cm3).

This function returns the atomic density (at/A3) given the atomic density 
(in at/A3) and the average atomic molar mass (g/mol).
    
    density [g/cm3] = 
           Atomic density [at/A3] / NA [at/mol] * molarM [g/mol] * 10^(24)
where molarM is the average molar mass per atom and 
NA = 6.02214076 10^(23) at/mol is the Avogadro's number.

For example, a water molecule has 3 atoms and the average molar mass is 
18 g/mol/3 = 6 g/mol/at
    Density = 0.1 at/A3 / 0.602214076 at/mol * 6 = 1 g/cm3 
        
If density is less or equal to 0, the function returns water density, 1.0 g/cm3

Input: atomic_density, molarMass
 -  atomic density : float, optional
        The atomic density in at/A3. The default is 0.1 at/A3.
 -  molarMass : float, optional
        The average molar mass per atom, in g/mol/atom. The default is 6.0.

Output: float
 -  density : float
        The macroscopic density, in g/cm3

Author: Gabriel Cuello
Date: Oct 25 2021    
--------------------------------------------------------------------------------
    """
    # Avogadros's number times 10^(-24): 
    # 6.02214076 10^(23) * 10^(-24) = 0.602214076
    # NA = 0.602214076
    result = atomic_density / 0.602214076 * molarM
    if atomic_density <= 0: 
        result = 1.0
        print ('Attention! Using the water density in getDensity function.')
    print('Density ',"{:.6f}".format(result),' g/cm3')
    return result 
#--------1---------2---------3---------4---------5---------6---------7---------8


#--------1---------2---------3---------4---------5---------6---------7---------8
def getNofAtoms(atoms):
    """
getNofAtoms function
    
This function calculates de number of atoms in basic unit.

Use: 
    n = getNofAtoms(atoms)

Input:
    - atoms: A dictionary with the number of atoms in the sample. The keys are the chemical symbols.

Output:
    - natoms: the number of atoms in a basic unit (a kind of molecule)

Created on Sat Oct 23 2021
@author: Gabriel Cuello
---------------------------------------
    """
    natoms = 0
    for key,value in atoms.items():
        natoms += value
    return natoms
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def getConcentrations(atoms):
    """
getConcentrations function

This function calculate the concentration of each atom in the sample.

Use: 
    c = getConcentrations(atoms)

Input:
    - atoms: A dictionary with the number of atoms in the sample. 
             The keys are the chemical symbols and values the corresponding number of atoms.

Output:
    - concentration: A dictionary with the same keys and the corresponding concentration values.

Requires:
    - getNofAtoms
    
Created on Sat Oct 23 2021
@author: Gabriel Cuello
---------------------------------------
    """
    concentration = {}
    natoms = getNofAtoms(atoms)
    print ('Atomic concentrations:')
    for key, value in atoms.items():
        concentration[key] = value/natoms
        print ('   ',key,'=',"{:.6f}".format(concentration[key]))
    return concentration
#--------1---------2---------3---------4---------5---------6---------7---------


def AtomicAvg(concentration,magnitude):
    """
AtomicAvg function
    
This is a function that makes an atomic average of a given magnitude.
    
Use:
    average = AtomicAvg(concentration,magnitude)

Input:
    - concentration: a dictionary with the concentration of each atom in the sample
    - magnitude: a dictionary with the magnitude to average
    Note: both dictionaries must have the same keys, i.e., the chemical symbols

Output:
    - The average of the magnitude

Created on Wed Dec 30 2020
@author: Gabriel Cuello
---------------------------------------
    """
    average = 0
    for key,value in concentration.items():
        average += float(value) * float(magnitude[key])
    return average




#--------1---------2---------3---------4---------5---------6---------7---------
#--------1---------2---------3---------4---------5---------6---------7---------
def XS_model(Eval, Eeff, composition_vec, bound_xs_vec, A_vec):
    """
    Calculates the neutron total cross section in the epithermal
    limit.
    
    $\sigma = \sum_i n_i \sigma_{b,i} \left(\frac{A_i}{A_i+1}\right)**2 \left( 1 + \frac{Eff}{2 A_i Eval}\right$

    Parameters
    ----------
    Eval : float
        Energy value
    Eeff : float
        Effective Energy = kB * Teff
    composition_vec : list
        ordered list of composition coefficients.
    bound_xs_vec : list
        ordered list of bound cross sections.
    A_vec : list
        ordered list of molar mass numbers.

    Returns
    -------
    float
        Total cross section for the given energy.
        
    Created on Wed Feb 23, 2022
    @author: José Robledo
    """
    composition_vec = np.array(composition_vec)
    bound_xs_vec = np.array(bound_xs_vec)
    A_vec = np.array(A_vec)
    
    a_vals = (A_vec/(A_vec + 1))**2
    
    result = np.sum(composition_vec * bound_xs_vec * a_vals * (1 + Eeff / (2 * A_vec * Eval)))
    return result


def scattering_probability(E0, scat_xs, composition_vec, abs_xs_vec):
    """
        Calculate the scattering probability for a given neutron incident energy.

    Parameters
    ----------
    E0 : float
        Incident energy in meV.
    bound_xs : float
        Scattering cross section
    composition_vec : list
        ordered list with compositions.
    abs_xs_vec : list
        ordered list with absorption scattering cross sections.

    Returns
    -------
    float
        Scattering probability for a neutron of energy E0.
    
    Created on Wed Feb 23, 2022
    @author: José Robledo
    """
    composition_vec = np.array(composition_vec)
    abs_xs_vec = np.array(abs_xs_vec)

    # calculate total absorption cross section 
    
    abs_xs = np.dot(composition_vec, abs_xs_vec) / np.sqrt(E0/25.3)
    
    # calculate total cross section
    total_xs = abs_xs + scat_xs
    
    return 1 - abs_xs/total_xs

#--------1---------2---------3---------4---------5---------6---------7---------