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

#--------1---------2---------3---------4---------5---------6---------7---------
import sys
import numpy as np                   # NumPy library
from scipy import integrate          # For integration
from datetime import datetime
import calendar       # To find the weekday


# --------1---------2---------3---------4---------5---------6---------7---------
def getDate():
    """This function produce a string with the current date.

    Produce a string with the current date in the format:
        Dayname DD/MM/YYYY HH:MM:SS

    Returns
    -------
    current_datetime: string

    author: Gabriel Cuello (ILL)
    date:   Jan 2022
    """
    # datetime object containing current date and time
    now = datetime.now()
    # dd/mm/YY H:M:S
    #    curr_date = datetime.today()
    #    day = ' '+calendar.day_name[curr_date.weekday()]+' '
    day = " " + calendar.day_name[now.weekday()] + " "
    current_datetime = day + now.strftime("%d/%m/%Y %H:%M:%S")
    return current_datetime


# --------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------


def sf_fcc(h,k,l):
    """ This function calculates the structure factor for a fcc lattice for
    a given plane hkl.
    
    Input
    -----
    The 3 indices of Miller corresponding to the plane to evaluate 
    """
    result = 1.0 + np.exp(-1j*np.pi*(k+l)) + np.exp(-1j*np.pi*(h+l)) + np.exp(-1j*np.pi*(h+k))
    return result

def sf_bcc(h,k,l):
    """ This function calculates the structure factor for a bcc lattice for
    a given plane hkl.
    
    Input
    -----
    The 3 indices of Miller corresponding to the plane to evaluate 
    """
    result = 1.0 + np.exp(-1j*np.pi*(h+k+l))
    return result

def sf_sc(h,k,l):
    """ This function calculates the structure factor for a sc lattice for
    a given plane hkl. Obviously, it always returns 1.
    
    Input
    -----
    The 3 indices of Miller corresponding to the plane to evaluate 
    """
    return 1.0

def reflections_fcc(wavelength=0.5, twotheta0=0.0, lattice=3.52):
    """ This function produce a list of the angular positions for
        fcc lattice, for a given wavelength, zero angle correction
        and a lattice parameter.
        
        This function will be used to fit the diffractogram of the
        Nickel powder sample, varying the wavelength and the zero
        angle correction.
        
        Note that nickel has a fcc structure with a=3.52024 Å.
    
    Input
    -----
    wavelength in Å, 
    twotheta0 in degrees
    lattice parameter in Å
    
    Output
    ------
    A list with the angular positions (in degrees) of all allowed 
    reflections for a fcc structure, in the range of 0 to 180 degrees.
    
    """
    moduleMax = int(2.0*lattice/wavelength)
    sinus = []
    for h in range(moduleMax):
        for k in range(moduleMax):
            for l in range(moduleMax):
                hkl = [h,k,l]
                module = np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2)
#                plane = "({} {} {})".format(*hkl)
                sf = sf_fcc(*hkl).real
                if sf*module > 0:
                    sintheta = module*wavelength/2.0/lattice
                    if sintheta < 1.0:
                        sinus.append(sintheta)
    rad_fcc = 2.0 * np.arcsin(sorted(set(sinus)))
    deg_fcc1 = twotheta0 + np.array(180.0/np.pi * rad_fcc)
    deg_fcc2 = -twotheta0 + np.array(sorted(180.0 - deg_fcc1))
    deg_fcc = np.concatenate([deg_fcc1, deg_fcc2])
    return deg_fcc

#--------1---------2---------3---------4---------5---------6---------7---------
def getMassNumber(molarMass=6.0):
    """Calculate the mass number (A).
    
    A = molar_mass / neutron_mass

    Parameters
    ----------
    molarMass : float, optional
        Molar mas in amu. The default is 6.0 (carbon molar mass).

    Returns
    -------
    A : float
        The mass number

    Date: Jan 13, 2022
    @author: Gabriel Cuello (ILL)
    ---------------------------------------
    """
    neutronMass = 1.0086649 # in amu
    A = molarMass / neutronMass
    return A
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def getFreeXS(BoundXS=5.08,A=50.9415):
    """Calculate the free cross section.
    
    The free cross section is calculated as:
        sigma_free = sigma_bound * A^2 /(1+A)^2
    where sigma_bound is the bound cross section (in barns) and A is the mass
    number.

    Parameters
    ----------
    BoundXS : float, optional
        Bound cross section. The default is 5.08.
    A : float, optional
        Mass number. The default is 50.9415.

    Returns
    -------
    FreeXS : float
        The free cross section (same units as the boun cross section)
    
    Date: Jan 13, 2022
    @author: Gabriel Cuello (ILL)
    ---------------------------------------
    """
    FreeXS = BoundXS * A**2 / (1.0+A)**2
    return FreeXS
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def extractDict(dictionary,key):
    """Extract a sub-dictionary from a bigger one.
    
    This function extracts a transversal dictionary from a bigger dictionary 
    containig dictionaries. The extracted dictionary will have the same keys 
    that the bigger one, but the values will be only those corresponding to 
    the input key.

    Parameters
    ----------
    dictionary : dictionary
        The big dictionary containing sub-dictionaries.
    key : string
        The key of the subdictionaries that will be extracted.

    Returns
    -------
    A dictionary with the keys of the big one and containing only the elements 
    with the given key.

    Example
    -------
    For the following dictionary
        elements = {'Ti':{'CohXS': 1.485,'IncXS': 2.87,  'AbsXS':6.09 },
                    'Nb':{'CohXS': 6.253,'IncXS': 0.0024,'AbsXS':1.15 },
                    'P': {'CohXS': 3.307,'IncXS': 0.005, 'AbsXS':0.172}}
    The instruction: incohXS = extractDict(elements,'IncXS') will produce:
        incohXS = {'Ti':2.87,'Nb':0.0024,'P':0.005}

    Usage
    -----
            dict = extractDict(dictionary,key)
    
    Date: Fri Oct 29 2021
    @author: Gabriel Cuello (ILL)
    ---------------------------------------
    """
    result = {}
    for key1,value in dictionary.items():
        result[key1] = dictionary[key1][key]
    return result
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def getCylVolume(diameter=5.0,height=50.0):
    """Calculate the volume of a cylinder.

    This function returns the volume of a cylinder, given its diameter and 
    height.

    The units of these two distances must be the same, and the result is given 
    in the same unit^3 (if distances are in mm, the volume will be in mm^3).

    Parameters
    ----------
    diameter : float, optional
        The diameter of the cylinder. The default is 5.0.
    height : float, optional
        The height of the cylinder. The default is 50.0.

    Returns
    -------
    float
        The volume of the cylinder.

    Date: Fri Oct 29 2021
    @author: Gabriel Cuello (ILL)
    ---------------------------------------
    """
    return np.pi * (diameter/2.0)**2 * height
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def getAtomicDensity(density=1.0,molarMass=6.0):
    """Calculate the atomic density (atoms/volume).

    This function returns the atomic density (at/A3) given the macroscopic 
    density (in g/cm3) and the average atomic molar mass (g/mol).
    
       Atomic density [at/A3] = 
           density [g/cm3] * NA [at/mol] / molarM [g/mol] * 10^(-24)
    
    where molarM is the average molar mass per atom and 
    NA = 6.02214076 10^(23) at/mol is the Avogadro's number.

    For example, a water molecule has 3 atoms and the average molar mass is 
    18/3 = 6, then
        AtomicDensity = 1 g/cm3 * 0.602214076 at/mol / 6 = 0.1 at/A3
        
    If density is less or equal to 0, the function returns water atomic 
    density, 0.1 at/A3.
    
    Parameters
    ----------
    density : float, optional
        The macroscopic density in g/cm3. The default is 1.0.
    molarMass : float, optional
        The average molar mass per atom, in g/mol/atom. The default is 6.0.

    Returns
    -------
    AtomicDensity : float
        The atomic density, in atoms/A3

    Usage
    -----
            AtomicDensity = getAtomicDensity(density,molarM)
    
    Date: Mon Oct 25 2021
    @author: Gabriel Cuello
    ---------------------------------------
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
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def getDensity(atomic_density=0.1,molarM=6.0):
    """Calculate the macroscopic density (g/cm3).

    This function returns the atomic density (at/A3) given the atomic 
    density (in at/A3) and the average atomic molar mass (g/mol).
    
       density [g/cm3] = 
           Atomic density [at/A3] / NA [at/mol] * molarM [g/mol] * 10^(24)
    
    where molarM is the average molar mass per atom and 
    NA = 6.02214076 10^(23) at/mol is the Avogadro's number.

    For example, a water molecule has 3 atoms and the average molar mass is 
    18/3 = 6
        Density = 0.1 at/A3 / 0.602214076 at/mol * 6 = 1 g/cm3 
        
    If density is less or equal to 0, the function returns water density, 
    1.0 g/cm3
    
    Parameters
    ----------
    atomic density : float, optional
        The atomic density in at/A3. The default is 0.1 at/A3.
    molarMass : float, optional
        The average molar mass per atom, in g/mol/atom. The default is 6.0.

    Returns
    -------
    density : float
        The macroscopic density, in g/cm3

    Usage
    -----
        density = getDensity(atomic_density,molarM)

    Date: Mon Oct 25 2021
    @author: Gabriel Cuello
    ---------------------------------------
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
#--------1---------2---------3---------4---------5---------6---------7---------


#--------1---------2---------3---------4---------5---------6---------7---------
def setExpInfo(Proposal="6-01-000",mainProposer="nobody",
               experimenters="nobody",LC="Cuello",otherUsers="nobody",
               startDate="01/01/2000",endDate="02/01/2000",
               environment="A",logBook=0,logPage=0,instrument="D4"):
    """Create a dictionary with information about the experiment.
    
    The function also print out a summary on the screen.
    
    Parameters
    ----------
    Proposal : string, optional
        Proposal number. The default is "6-01-000".
    mainProposer : string, optional
        Main proposer. The default is "nobody".
    experimenters : string, optional
        Experiments on site. The default is "nobody".
    LC : string, optional
        The local contact. The default is "Cuello".
    otherUsers : string, optional
        Users in the proposal but no on site. The default is "nobody".
    startDate : string, optional
        Date of the starting day. The default is "01/01/2000".
    endDate : string, optional
        Date of the ending day. The default is "02/01/2000".
    environment : string, optional
        Sample environment. The default is "A".
    logBook : integer, optional
        Number of the D4 logbook. The default is 0.
    logPage : inetger, optional
        Page in the D4 logbook. The default is 0.
    instrument : string, optional
        Instrument. The default is "D4".

    Returns
    -------
    experiment : dictionary
        A dictionary containing the information about the experiment.

    Date: Sun Dec 26 2021
    @author: Gabriel Cuello
    ---------------------------------------
    """
    experiment = {'Proposal': Proposal}
    experiment['MP'] = mainProposer
    experiment['experimenters'] = experimenters
    experiment['LC'] = LC
    experiment['otherUsers'] = otherUsers
    experiment['instr'] = instrument
    experiment['startDate'] = startDate
    experiment['endDate'] = endDate
    experiment['environment'] = environment
    experiment['logBook'] = logBook
    experiment['logPage'] = logPage
    print (30*'-')
    print ('Experiment')
    print (4*' ','Instrument: {}'.format(experiment['instr']))
    print (4*' ','Proposal: {}'.format(experiment['Proposal']))
    print (4*' ','Main proposer: {}'.format(experiment['MP']))
    print (4*' ','Other users: {}'.format(experiment['otherUsers']))
    print (4*' ','On-site users: {}'.format(experiment['experimenters']))
    print (4*' ','Local contact: {}'.format(experiment['LC']))
    print (4*' ','Starting date: {0} ---->  Ending date: {1}'.
           format(experiment['startDate'],experiment['endDate']))
    print (4*' ','Sample environment: {}'.format(experiment['environment']))
    print (4*' ','D4 notebook: {0}  Page: {1}'.
           format(experiment['logBook'],experiment['logPage']))
    print ()
    return experiment
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def setBeamInfo(zeroAngle=0,wavelength=0.5,
                LohenSlit=2.5,GamsSlit=-2.5,topFlag=25.0,bottomFlag=-25.0):
    """Create a dictionary with information about the beam.
    
    The function also print out a summary on the screen.

    Parameters
    ----------
    zeroAngle : float, optional
        The zero-angle correction, as obtained from the calibration with Ni
        powder. The default is 0 deg.
    wavelength : float, optional
        The incident wavelength, as obtained from the calibration with Ni
        powder. The default is 0.5 A.
    LohenSlit : float, optional
        The position of the vertical slit (left in downstream direction). 
        This is the slit at Lohengrin side. The default is 2.5 mm.
    GamsSlit : float, optional
        The position of the vertical slit (right in downstream direction). 
        This is the slit at GAMS side. The default is -2.5 mm.
    topFlag : float, optional
        The position of the top horizontal slit. The default is 25.0 mm.
    bottomFlag : float, optional
        The position of the top horizontal slit. The default is -25.0 mm.

    Returns
    -------
    beam : dictionary
        A dictionary containing the information about the beam.

    Date: Sun Dec 26 2021
    @author: Gabriel Cuello
    ---------------------------------------
    """
    beam = {'zero': zeroAngle}    # zero angle correction
    beam['wlength'] = wavelength   # neutron wavelength
    # Vertical slits defining the horizontal size (width) of the beam
    # Lohengrin side, i.e., on the left in downstream direction
    # Gams side, i.e., on the right in downstream direction
    beam['LohenSlit'] = LohenSlit   # in mm
    beam['GamsSlit'] = GamsSlit     # in mm
    beam['width'] = beam['LohenSlit'] - beam['GamsSlit'] # Beam width (in mm)
    # Horizontal flags defining the vertical size (height) of the beam
    beam['topFlag'] = topFlag       # in mm
    beam['botFlag'] = bottomFlag    # in mm
    beam['height'] = beam['topFlag']-beam['botFlag'] # Beam height (in mm) 
    print (30*'-')
    print ('Beam characteristics')
    print (4*' ','Wavelength = {:8.6g} A'.format(beam['wlength']))
    print (4*' ','Zero angle = {:8.6g} deg'.format(beam['zero']))
    print (4*' ','Dimensions: (Width = {0:.6g} mm)x(Height = {1:.6g} mm)'.
           format(beam['width'],beam['height']))
    print ()
    return beam
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def setCanInfo(material='Vanadium',shape='Cylinder',
               outerDiam=5,innerDiam=4.8,height=60.0):
    can = {'material': material}
    can['shape'] = shape
    can['outerDiam'] = outerDiam   # outer diameter (in mm)
    can['innerDiam'] = innerDiam   # inner diameter (in mm)
    can['height'] = height     # height (in mm)
    can['wallThickness'] = (can['outerDiam']-can['innerDiam'])/2.0
    print (30*'-')
    print ('Container')
    print (4*' ','Type: {0} {1}'.format(can['material'],can['shape']))
    print (4*' ','Outer diameter = {} mm'.format(can['outerDiam']))
    print (4*' ','Inner diameter = {} mm'.format(can['innerDiam']))
    print (4*' ','Wall thickness = {:.3g} mm'.format(can['wallThickness']))
    print (4*' ','Height = {} mm'.format(can['height']))
    print ()
    return can
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def setBinInfo(AngularResolution=0.125,
               AMin=0.0, AMax=140.0, AStep=0.125,
               QMin=0.0, QMax=23.5,  QStep=0.02,
               RMin=0.0, RMax=20.0,  RStep=0.01):
    # Experimental width of an angular channel (in degrees)
    # 0.125 deg (January 2022)
    binning = {'Ares': AngularResolution}
    # Tuples with initial, final and step for each scale
    #   in degrees for angular scale
    #   in 1/A for Q-scale
    #   in A for R-scale
    binning['Abin'] = (AMin,AMax,AStep)
    binning['Qbin'] = (QMin,QMax,QStep)
    binning['Rbin'] = (RMin,RMax,RStep)
    
    binning['NbrPointsA'] = int((binning['Abin'][1]-binning['Abin'][0])
                                /binning['Abin'][2])
    binning['NbrPointsQ'] = int((binning['Qbin'][1]-binning['Qbin'][0])
                                /binning['Qbin'][2])
    binning['NbrPointsR'] = int((binning['Rbin'][1]-binning['Rbin'][0])
                                /binning['Rbin'][2])
    print (30*'-')
    print ('Binning')
    print (4*' ','Angular channel width = {:.3g} deg'.format(binning['Ares']))
    print (4*' ','In angle: from {0:.3g} deg to {1:.3g} deg, in steps of {2:.3f} deg, thus {3:.5g} points.'.
           format(binning['Abin'][0],binning['Abin'][1],binning['Abin'][2],binning['NbrPointsA']))
    print (4*' ','In Q: from {0:.3g} 1/A to {1:.3g} 1/A, in steps of {2:.3f} 1/A, thus {3:.5g} points.'.
           format(binning['Qbin'][0],binning['Qbin'][1],binning['Qbin'][2],binning['NbrPointsQ']))
    print (4*' ','In R: from {0:.3g} A to {1:.3g} A, in steps of {2:.3f} A, thus {3:.5g} points.'.
           format(binning['Rbin'][0],binning['Rbin'][1],binning['Rbin'][2],binning['NbrPointsR']))
    print ()
    return binning
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def setNumorInfo(total=('Exp',0,1),container=('Can',0,1),
            environment=('Env',0,1),
            nickel=('Ni5',0,1),vanadium=('Van',0,1),absorber=('Abs',0,1),
            sample01 =('S01',0,1),sample02 =('S02',0,1),sample03 =('S03',0,1),
            sample04 =('S04',0,1),sample05 =('S05',0,1),sample06 =('S06',0,1),
            sample07 =('S07',0,1),sample08 =('S08',0,1),sample09 =('S09',0,1),
            sample10 =('S10',0,1),sample11 =('S11',0,1),sample12 =('S12',0,1),
            sample13 =('S13',0,1),sample14 =('S14',0,1),sample15 =('S15',0,1),
            sample16 =('S16',0,1),sample17 =('S17',0,1),sample18 =('S18',0,1),
            sample19 =('S19',0,1),sample20 =('S20',0,1),sample21 =('S21',0,1)):


    numors = {'experiment': total,'container': container,
              'environment': environment,
              'nickel': nickel,'vanadium': vanadium,'absorber': absorber,
              'sample01': sample01,'sample02': sample02,'sample03': sample03,
              'sample04': sample04,'sample05': sample05,'sample06': sample06,
              'sample07': sample07,'sample08': sample08,'sample09': sample09,
              'sample10': sample10,'sample11': sample11,'sample12': sample12,
              'sample13': sample13,'sample14': sample14,'sample15': sample15,
              'sample16': sample16,'sample17': sample17,'sample18': sample18,
              'sample19': sample19,'sample20': sample20,'sample21': sample21}

    print ('Numors')

    if total[1]!=0:
        print ('{}Total {}: {}-{}, {} numors.'.format(
            4*' ',*numors['experiment'],numors['experiment'][2]-numors['experiment'][1]+1))
    if container[1]!=0:
        print ('{}{}: {}-{}, {} numors.'.format(
            4*' ',*numors['container'],numors['container'][2]-numors['container'][1]+1))
    if environment[1]!=0:
        print ('{}{}: {}-{}, {} numors.'.format(
            4*' ',*numors['environment'],numors['environment'][2]-numors['environment'][1]+1))
    if nickel[1]!=0:
        print ('{}{}: {}-{}, {} numors.'.format(
            4*' ',*numors['nickel'],numors['nickel'][2]-numors['nickel'][1]+1))
    if vanadium[1]!=0:
        print ('{}{}: {}-{}, {} numors.'.format(
            4*' ',*numors['vanadium'],numors['vanadium'][2]-numors['vanadium'][1]+1))
    if absorber[1]!=0:
        print ('{}{}: {}-{}, {} numors.'.format(
            4*' ',*numors['absorber'],numors['absorber'][2]-numors['absorber'][1]+1))

    if sample01[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample01'],numors['sample01'][2]-numors['sample01'][1]+1))
    if sample02[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample02'],numors['sample02'][2]-numors['sample02'][1]+1))
    if sample03[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample03'],numors['sample03'][2]-numors['sample03'][1]+1))
    if sample04[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample04'],numors['sample04'][2]-numors['sample04'][1]+1))
    if sample05[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample05'],numors['sample05'][2]-numors['sample05'][1]+1))
    if sample06[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample06'],numors['sample06'][2]-numors['sample06'][1]+1))
    if sample07[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample07'],numors['sample07'][2]-numors['sample07'][1]+1))
    if sample08[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample08'],numors['sample08'][2]-numors['sample08'][1]+1))
    if sample09[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample09'],numors['sample09'][2]-numors['sample09'][1]+1))
    if sample10[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample10'],numors['sample10'][2]-numors['sample10'][1]+1))
    if sample11[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample11'],numors['sample11'][2]-numors['sample11'][1]+1))
    if sample12[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample12'],numors['sample12'][2]-numors['sample12'][1]+1))
    if sample13[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample13'],numors['sample13'][2]-numors['sample13'][1]+1))
    if sample14[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample14'],numors['sample14'][2]-numors['sample14'][1]+1))
    if sample15[1]!=0:
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample15'],numors['sample15'][2]-numors['sample15'][1]+1))
    if sample16[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample16'],numors['sample16'][2]-numors['sample16'][1]+1))
    if sample17[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample17'],numors['sample17'][2]-numors['sample17'][1]+1))
    if sample18[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample18'],numors['sample18'][2]-numors['sample18'][1]+1))
    if sample19[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample19'],numors['sample19'][2]-numors['sample19'][1]+1))
    if sample20[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample20'],numors['sample20'][2]-numors['sample20'][1]+1))
    if sample21[1]!=0:    
        print ('{}Sample {}: {}-{}, {} numors.'.format(
            4*' ',*numors['sample20'],numors['sample21'][2]-numors['sample21'][1]+1))
    return numors
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def setVanaInfo(IncXS=5.08,CohXS=0.0184,ScaXS=5.1,AbsXS=5.08,
                CohSL=-0.3824,molarM=50.9415,NAtoms=1.0,
                Diam=6.08,Height=50.0,density=6.51,wavelength=0.5):
    vana = {}
    vana['IncXS'] = IncXS # Incoherent scattering cross section (sigma bound), in barns
    vana['CohXS'] = CohXS # Coherent scattering cross section, in barns
    vana['ScaXS'] = ScaXS # Coherent scattering cross section, in barns
    vana['AbsXS'] = AbsXS # Absorption cross section, in barns
    vana['AbsWW'] = vana['AbsXS'] * wavelength/1.8 
    # Coherent scattering length (fm)
    vana['CohSL'] = CohSL # Coherent scattering length, in fm
    vana['SelfQ0'] = vana['IncXS']/4./np.pi
    vana['diam'] = Diam # diameter (mm)
    vana['NAtoms'] = NAtoms # Number of atoms in a unit
    vana['molarM'] = molarM # Molar mass of vanadium, in g/mol
    vana['den_gcc'] =density # Macroscopic density, in g/cm3
    vana['den_aac'] = getAtomicDensity(density=vana['den_gcc'],molarMass=vana['molarM'])
    vana['MassNbr'] = getMassNumber(molarMass=vana['molarM']) # Ratio mass to neutron mass
    vana['FreeIncXS'] = getFreeXS(BoundXS=vana['IncXS'],A=vana['MassNbr']) 
    vana['FreeCohXS'] = getFreeXS(BoundXS=vana['CohXS'],A=vana['MassNbr'])
    vana['FreeScaXS'] = vana['FreeIncXS'] + vana['FreeCohXS']
    vana['volume'] = getCylVolume(diameter=vana['diam'],height=Height)/1000.0
    vana['MacroXS'] = vana['den_aac']*(vana['ScaXS']+vana['AbsWW']) # Units: 1/cm
    # pi*D/4 is the equivalent thickness of a plate corresponding to a cylinder on diameter D
    # The factor 0.1 is to convert the diameter in mm to cm. 
    vana['Transm'] = np.exp(-vana['MacroXS']*np.pi*vana['diam']/4.0*0.1)
    print (30*'-')
    print ('Standard of vanadium')
    print(4*' ','The standard is a cylinder of {:.3g} mm of diameter.'.format(vana['diam']))
    print(4*' ','Volume in the beam = {:.3g} cm3.'.format(vana['volume']))
    print ()
    print (4*' ','Bound coherent cross section: {:.6g} barns/atom'.format(vana['CohXS']))
    print (4*' ','Bound incoherent cross section: {:.6g} barns/atom'.format(vana['IncXS']))
    print (4*' ','Bound scattering cross section: {:.6g} barns/atom'.format(vana['ScaXS']))
    print ()
    print (4*' ','Free coherent cross section: {:.6g} barns/atom'.format(vana['FreeCohXS']))
    print (4*' ','Free incoherent cross section: {:.6g} barns/atom'.format(vana['FreeIncXS']))
    print (4*' ','Free scattering cross section: {:.6g} barns/atom'.format(vana['FreeScaXS']))
    print ()
    print (4*' ','Absorption cross section: {:.6g} barns/atom'.format(vana['AbsXS']))
    print (4*' ','Absorption cross section at {:.6g} A: {:.6g} barns/atom'.format(wavelength,vana['AbsWW']))
    print ()
    print(4*' ','b**2 = sigma/4/pi at Q=0 (self) is',"{:.6g} barns/sterad/atom.".format(vana['SelfQ0']))
    print(4*' ','b**2 = sigma/4/pi at Q=infty is',"{:.6g} barns/sterad/atom".format(vana['FreeIncXS']/4./np.pi))
    print ()
    print(4*' ','The molar mass is',"{:.6g} g/mol.".format(vana['molarM']))
    print(4*' ','Mass number, A = ',"{:.6g} (= mass/neutron_mass)".format(vana['MassNbr']))
    print ()
    print(4*' ','Density = ',"{:.6g} g/cm3.".format(vana['den_gcc']))
    print(4*' ','Atomic density = {:.6g} atoms/A3.'.format(vana['den_aac']))
    print(4*' ','Macroscopic cross section = {:.6g} 1/cm.'.format(vana['MacroXS']))
    print(4*' ','Transmission = {:.6g} .'.format(vana['Transm']))
    print ()
    return vana
#--------1---------2---------3---------4---------5---------6---------7---------


#--------1---------2---------3---------4---------5---------6---------7---------
def setSampleInfo(comp_sample,atoms_sample,natoms_sample,
                          elements,c_sample,wavelength=0.5,beamH=50.0,
                          vanavol=1.0,vanadens=0.0769591,
                          height=60.0,diameter=5.0,mass=1.0,density=1.0,
                          title='sample'):
    # Creating a new dictionary with the sample basic information
    sample ={'Title':title}
    
    sample['height'] = height # sample height (in mm)
    # This is the height of the sample seen by the beam
    sample['heightInBeam'] = min(sample['height'],beamH)
    # sample diameter, i.e., container inner diameter
    # If it is a delf-contained sample (no container), put the sample diameter (in mm)
    sample['diam'] = diameter 
    sample['mass'] = mass # sample mass (in g)
    sample['NAtoms'] = natoms_sample # number of atoms in 1 unit
    
    # average molar mass per atom (in g/mol/atom)
    sample['molarM'] = AtomicAvg(c_sample,extractDict(elements,'molarM'))
    
    # Volume = pi*radius^2*height (in cm3)
    sample['volume'] = getCylVolume(sample['diam'],sample['height'])/1000.0 
    # Volume of the sample in the beam (in cm3)
    sample['volumeInBeam'] = getCylVolume(diameter=sample['diam'],height=sample['heightInBeam'])/1000.0 
    
    sample['den_gcc'] = density # macroscopic density of the sample (in g/cm3)
    
    # Effective density, simply mass/volume (g/cm3)
    sample['effden_gcc'] = sample['mass']/sample['volume'] 
    
    # Packing fraction
    sample['packing'] = sample['effden_gcc']/sample['den_gcc']
    # Filling fraction: 
    #    1 if completely-filled container, less than 1 for a partially-filled container
    sample['filling'] = sample['heightInBeam']/beamH
    
    # atomic density (atoms/A3)
    sample['den_aac'] = getAtomicDensity(density=sample['den_gcc'],molarMass=sample['molarM'])
    
    # Effective atomic density (atoms/A3)
    sample['effden_aac'] = getAtomicDensity(density=sample['effden_gcc'],molarMass=sample['molarM'])
    
    # Coherent scattering cross section (barns)
    sample['CohXS'] = AtomicAvg(c_sample,extractDict(elements,'CohXS')) 
    # Incoherent scattering cross section (barns)
    sample['IncXS'] = AtomicAvg(c_sample,extractDict(elements,'IncXS')) 
    # Scattering cross section (barns)
    sample['ScaXS'] = AtomicAvg(c_sample,extractDict(elements,'ScaXS')) 
    # Absorption cross section (barns)
    sample['AbsXS'] = AtomicAvg(c_sample,extractDict(elements,'AbsXS')) 
    # Absorption cross section at working wavelength (barns)
    sample['AbsWW'] = sample['AbsXS'] * wavelength/1.8 
    # Coherent scattering length (fm)
    sample['CohSL'] = AtomicAvg(c_sample,extractDict(elements,'CohSL')) 
    # Atomic number (amu)
    sample['MassNbr'] = AtomicAvg(c_sample,extractDict(elements,'MassNbr')) 
    # Free coherent scattering cross section (barns)
    sample['FreeCohXS'] = getFreeXS(BoundXS=sample['CohXS'],A=sample['MassNbr'])
    # Free incoherent scattering cross section (barns)
    sample['FreeIncXS'] = getFreeXS(BoundXS=sample['IncXS'],A=sample['MassNbr'])
    # Free scattering cross section (barns)
    sample['FreeScaXS'] = getFreeXS(BoundXS=sample['ScaXS'],A=sample['MassNbr'])
    # Self incoherent scattering cross section at Q=0 (barns), selfQ0 = IncXS/4/pi
    sample['SelfQ0'] = sample['IncXS']/4.0/np.pi
    
    # Ratio of sample volume to vanadium volume
    sample['VolRatioSV'] = sample['volumeInBeam']/vanavol
    # Ratio of sample density to vanadium desity
    sample['DenRatioSV'] = sample['effden_aac']/vanadens
    
    print (80*'-')
    print ('Sample',sample['Title'])
    
    text = ' '
    print ()
    print ('Composition:')
    for key,value in comp_sample.items():
        text += '+ '+str(value)+' of '+key+' '
    #    print (8*' ',key,value)
    print (4*' ','A "unit" is',text[3:])
    
    text = ' '
    print ()
    print ('Atomic composition:')
    for key,value in atoms_sample.items():
        text += '+ '+str(value)+' '+key+' atoms '
    #    print (8*' ',key,value)
    print (4*' ','A "unit" has',text[3:])
    print (4*' ','Number of atoms in one unit = {} atoms'.format(sample['NAtoms']))
    
        
    print ()
    print ('Atomic concentrations:')
    for key,value in c_sample.items():
        print (8*' ',key,"{:.6f}".format(value))
    
    print ()
    print ('Average molar mass: {:.6g} g/mol/atom'.format(sample['molarM']))
    print ('Mass number: {:.6g} amu'.format(sample['MassNbr']))
    print ()
    print ('Coherent cross section: {:.6g} barns/atom'.format(sample['CohXS']))
    print ('Incoherent cross section: {:.6g} barns/atom'.format(sample['IncXS']))
    print ('Scattering cross section: {:.6g} barns/atom'.format(sample['ScaXS']))
    print ()
    print ('Free coherent cross section: {:.6g} barns/atom'.format(sample['FreeCohXS']))
    print ('Free incoherent cross section: {:.6g} barns/atom'.format(sample['FreeIncXS']))
    print ('Free scattering cross section: {:.6g} barns/atom'.format(sample['FreeScaXS']))
    print ()
    print ('Absorption cross section: {:.6g} barns/atom'.format(sample['AbsXS']))
    print ('Absorption cross section at {:.6g} A: {:.6g} barns/atom'.format(wavelength,sample['AbsWW']))
    print ()
    print ('Coherent scattering length: {:.6g} fm'.format(sample['CohSL']))
    print ('Self scattering cross section at Q=0: {:.6g} barns/steradian/atom'.format(sample['SelfQ0']))
    print ()
    print ('Density: {:.6g} g/cm3'.format(sample['den_gcc']))
    print ('Atomic density: {:.6g} atoms/A3'.format(sample['den_aac']))
    print ()
    print ('Cylindrical sample')
    print (4*' ','Diameter: {:.6g} mm'.format(sample['diam']))
    print (4*' ','Height: {:.6g} mm'.format(sample['height']))
    print (4*' ','Volume: {:.6g} cm3'.format(sample['volume']))
    print (4*' ','Mass: {:.6g} g'.format(sample['mass']))
    print (4*' ','Effective density: {:.6g} g/cm3'.format(sample['effden_gcc']))
    print (4*' ','Effective atomic density: {:.6g} atoms/A3'.format(sample['effden_aac']))
    print ()
    print (4*' ','Sample/Vanadium volume fraction: {:.6g}'.format(sample['VolRatioSV']))
    print (4*' ','Sample/Vanadium density fraction: {:.6g}'.format(sample['DenRatioSV']))
    print (4*' ','Packing fraction: {:.6g}'.format(sample['packing']))
    print (4*' ','Filling fraction: {:.6g}'.format(sample['filling']))
    
    return sample
#--------1---------2---------3---------4---------5---------6---------7---------


#--------1---------2---------3---------4---------5---------6---------7---------
def stop(message):
    """
This function allows stopping the program at this point.
Requires an input from the console:
    y: continue running
    n: stops the program
    
    Parameters
    ----------
    message : STRING
        This a message to be printed just before asking whether continuing or not.

    Returns
    -------
    None

    """
    print ()
    print (30*'<>')
    print (message)
    answer = input("Do you want to continue? (y/n) :")
    if answer == 'n':
        sys.exit("You stopped the program. Bye!")
    return
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def ang2q(x,wlength=0.5):
    """
ang2q function

This function converts a list of angles in the corresponding Q values, using the Bragg law
         Q = 4 pi / lambda sin(angle/2 / 180 * pi)
where angle is the scattering angle (2theta).
    
Use: 
    Q = ang2q(angles,wavelength)

Input:
    - x: A list with the angles (in degrees).
    - wlenght: The wavelength (in A).

Output:
    - q: A list with the Q values (in 1/A).

Created on Sun Oct 24 2021
@author: Gabriel Cuello
---------------------------------------
    """
    q = 4.0*np.pi/wlength * np.sin(np.array(x)/2.0*np.pi/180.0)
    return q
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def q2ang (x,wlength=0.5):
    """
q2ang function

    This function converts a list of Q values in the corresponding angles, 
    using the Bragg law
             angle = 2 arcsin(Q lambda / 4 pi) / pi * 180
    where angle is the scattering angle (2theta).

    Usage
    ----- 
        angles = ang2q(Q,wavelength)

Input:
    - x: A list with the Q values (in 1/A).
    - wlenght: The wavelength (in A).

Output:
    - angles: A list with the angles (in degrees).

    Created on Sun Oct 24 2021
    @author: Gabriel Cuello
    ---------------------------------------
    """
    ang = 360.0/np.pi * np.arcsin(np.array(x)*wlength/4.0/np.pi)
    return ang
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
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

#--------1---------2---------3---------4---------5---------6---------7---------
def inelastic(x,A=51.0,lowQ=0.4,Q0=7.0,dQ=2.4):
    """
Inelastic behaviour function

Use: 
    y = inelastic(x,A,lowQ,Q0,dQ)

This function takes into account the inelastic behaviour of the incoherent signal.

        inelastic(x) = lowQ * [(1 + A**2/(1+A**2) * delta ] / (1+delta)

where delta = exp[(x-q0)/dq]

Note that inelastic(infty) = lowQ * A**2 / (1+A)**2
    i.e. sigma_free = sigma_bound * A**2 / (1+A)**2
    or   b_free = b_bound * A / (1+A)

Input:
    - x: a range of x
    - A: mass number (M/m, where M is the atomic mass and m is the neutron mass)
    - lowQ: limiting value at Q=0, i.e., sigma_bound/4/pi
    - Q0: inflection point of the sigmoidal function
    - dQ: width of the sigmoidal function

Output:
    - The values of the function at point x

Created on Wed Dec 30 2020
Modified on Sat Oct 09 2021
@author: Gabriel Cuello
---------------------------------------
    """
    delta = np.exp((np.array(x)-Q0)/dQ)
    return lowQ*(1+A**2/(1+A)**2*delta)/(1+delta)
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def Lorch(q,qmax=23.5):
    """
Lorch function

Use: 
    y = Lorch(Q,Qmax)

This is a window function.
       If Q<0 or Q>=Qmax, Lorch(Q)=0.0
       If 0<Q<Qmax, Lorch(Q) = sin(pi*Q/Qmax) / (pi*Q/Qmax)
       If Q=0, Lorch(0) = 1.0
    
Created on Wed Dec 30 2020
Modified on Sat Oct 09 2021
@author: Gabriel Cuello
---------------------------------------
    """
    if qmax<=0: # qmax must be positive, otherwise the function
                # produces an error message and returns -999
        print ('ERROR: Non positive value for qmax in Lorch function.')
        return -999

    if ((q<0.0) or (q>=qmax)): # Outside the range [0,qmax) returns 0
        return 0.0

    elif q == 0.0: # Special case sin(x)/x = 1 for x = 0
        return 1.0
    else:
        a = q*np.pi/qmax # sin(x)/x
    return np.sin(a)/a
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def getSineIntegral(x):
    """
    This function evaluates the Sine Integral function:
        Si(x) = int_0^x sin(t)/t dt

    Returns
    -------
    result : float
        DESCRIPTION.

    """
    npoints = 10000
    si = [] # A list containing the values of the window function
    t = []
    t.append(0.0)
    si.append(1.0) # A list containing the values of the window function
    for i in range(1,npoints):
        t.append(i*x/float(npoints))
        si.append(np.sin(t[i])/t[i])
    result = integrate.simps(si,t)
    return result
#--------1---------2---------3---------4---------5---------6---------7---------
    
    
#--------1---------2---------3---------4---------5---------6---------7---------
def step(q,qmax=23.5):
    """
Step function

Use: 
    y = step(Q,Qmax)

This is a rectangular window function
        returns 1 for 0<=Q<Qmax
        returns 0 otherwise
    
Created on Sat Oct 09 2021
@author: Gabriel Cuello
---------------------------------------
    """
    if qmax<=0: # qmax must be positive, otherwise the function
                # produces an error message and returns -999
        print ('ERROR: Non positive value for qmax in step function.')
        return -999

    if ((q<0.0) or (q>=qmax)): # Outside the range [0,qmax) returns 0
        return 0.0
    else:
        return 1.0
#--------1---------2---------3---------4---------5---------6---------7---------


#--------1---------2---------3---------4---------5---------6---------7-----
def LorchN(x,xmax=0):
#   The default value of xmax is 0
#   If xmax less or equal to 0 or greater than the last element
#   then xmax is the last element of the input x array.
    if (xmax <= 0) or (xmax > x[-1]):
        xmax = x[-1]

    lorch = np.ones(len(x))
    for i in range(len(x)):
#       Outside the range [0,xmax), lorch=0
        if ((x[i]<0.0) or (x[i]>=xmax)): 
            lorch[i] = 0.0
        elif x[i] != 0.0:
            a = x[i]*np.pi/xmax # sin(x)/x
            lorch[i] = np.sin(a)/a
        else:
            lorch[i] = 1.0

#   Integral of the Lorch function for normalisation
    integralLorch = integrate.simps(lorch,x)
    for i in range(len(x)):
        lorch[i] *= xmax / integralLorch
    return lorch
#--------1---------2---------3---------4---------5---------6---------7-----

#--------1---------2---------3---------4---------5---------6---------7-----
def getSineFT(vvec,yvec,uvec,vmax=0,c=1.0,s=1.0,w=0):
# Number of points in the function to be transformed
    nbr_v = len(vvec)  # abcissa
    nbr_y = len(yvec)  # ordinate
    # These numbers must be equal, otherwise print error message
    if nbr_v != nbr_y:
        print ('ERROR')
        return None

#   The default value of vmax is 0
#   If vmax less or equal to 0 or greater than the last element
#   then vmax is the last element of the input v array.
    if (vmax <= 0) or (vmax > vvec[-1]):
        vmax = vvec[-1]

    win=np.ones(nbr_v)
    if w == 1:
        win = LorchN(vvec,vmax)

    stf = np.zeros(len(uvec))
    for i in range(len(uvec)):
#        integ = soq*q*np.sin(q*r)
        integ = (yvec-s)*vvec*np.sin(vvec*uvec[i])*win
#        result = integrate.simps(integ,q)
#        pcf.append(result)
#        integ = soq*q*np.sin(q*r)*win
        stf[i] = 2.0 / np.pi * c * integrate.simps(integ,vvec)
    # Pair correlation function or G(r)
    return stf
#--------1---------2---------3---------4---------5---------6---------7-----

#--------1---------2---------3---------4---------5---------6---------7---------
def sineFT(StrFactor,nr,qmax=23.5,density=0.1,constant=1.0,selfsca=1.0,window=0):
    """
Function sineFT

This function performs the sinus Fourier tranform of a given function.
    
Use: 
    sinusFT(density,Q,SoQ,Qmax,constant,Rrange,window)
    
The integral is performed from 0 to Qmax, and the integrand is: 
        constant * (S(Q)-self) * Q * sin(Q*R) * window(Q)
        
Input:
    - density: the atomic density of the sample (in atoms/A3)
    - Q: list containing the Q range (in 1/A)
    - SoQ: list containing the structure factor
    - Qmax: maximum value of the experimental Q-range (in 1/A)
    - constant: a multiplicative constant (it should be 1)
    - constant: an additve constant (it should be 1)
    - Rrange: list containing the R range (in A)
    - window: 0 or 1, for step or Lorch window function, respectively

Output: 
    pcf,pdf,rdf,tor,run

The output are lists with the different real space functions.
    - pcf: pair correlation function
    - pdf: pair distribution function
    - rdf: radial distribution function or linear density
    - tor: rdf/R, which produces more symmetrical peaks
    - run: running integral of the rdf, i.e., integral from 0 to R

Created on Wed Dec 30 2020
Modified on Sat Oct 09 2021
@author: Gabriel Cuello
---------------------------------------
    """
    soq = StrFactor[:,1]
    q = StrFactor[:,0]
    pcf = []
    pdf = []
    rdf = []
    tor = []
    run = []
#    infile = np.genfromtxt(filename, skip_header= 7) #creates an array from xA, yA, zA 
#    q  = infile[:, 0]
#    soq  = infile[:, 1]  
#    #err  = infile[:, 2]  

    win = [] # A list containing the values of the window function

    if window ==0:
        for qu in q:
            win.append(step(qu,qmax))
    else:
        for qu in q:
            win.append(Lorch(qu,qmax))
#       Integral of the Lorch function for normalisation
        integralLorch = integrate.simps(win,q)
#        print ('Normalisation of Lorch function: ',integralLorch/qmax)
        for i in range(len(win)):
            win[i] = win[i] * qmax / integralLorch

    deltaR = nr[1]-nr[0]
    for r in nr:
#        integ = soq*q*np.sin(q*r)
        integ = (soq-selfsca)*q*np.sin(q*r)*win
#        result = integrate.simps(integ,q)
#        pcf.append(result)
#        integ = soq*q*np.sin(q*r)*win
        result = 2.0 / np.pi * constant * integrate.simps(integ,q)
    # Pair correlation function or G(r)
        pcf.append(result)
    # Pair distribution function or g(r)
        if r <= 0:
            pdf.append(0.0)
        else:
            pdf.append(result / 4.0 / np.pi / r /density + 1.0)
    # Radial distribution function or RDF(r)
        rdf.append(result * r + 4.0 * np.pi * density * r**2)
    # T(r) = RDF(r)/r; symmetric peaks for fitting
        tor.append(result + 4.0 * np.pi * density * r)
    # Running integral of the RDF(r), i.e., integral from 0 to r
    for r in nr:
        ind = 1+int(r/deltaR)
        xr = nr[0:ind]
        integ = np.array(rdf)[0:ind]
        result = integrate.simps(integ,xr)
        run.append(result)
    rrr = np.array(nr)
    pcf = np.array(pcf)
    pdf = np.array(pdf)
    rdf = np.array(rdf)
    tor = np.array(tor)
    run = np.array(run)
    rrr = rrr.reshape(rrr.shape[0],1)
    pcf = pcf.reshape(pcf.shape[0],1)
    pdf = pdf.reshape(pdf.shape[0],1)
    rdf = rdf.reshape(rdf.shape[0],1)
    tor = tor.reshape(tor.shape[0],1)
    run = run.reshape(run.shape[0],1)
    fou = np.concatenate((rrr,pcf),axis=1)
    fou = np.concatenate((fou,pdf),axis=1)
    fou = np.concatenate((fou,rdf),axis=1)
    fou = np.concatenate((fou,tor),axis=1)
    fou = np.concatenate((fou,run),axis=1)
    return fou
#--------1---------2---------3---------4---------5---------6---------7---------


def backFT(qu,ar,pdf,density,cut=0):
    pdf_cut  = pdf.copy()
    #pdf_lorch_cut = pdf_lorch.copy()

#     for i in range(len(pdf)-1,0,-1):
#         if pdf[i] < 0:
#             break
#     for j in range(i+1):
#         pdf_cut[j] = 0.0
    if cut == -1:
        for i in range(len(pdf)-1,0,-1):
            if pdf[i] < 0:
                break
        for j in range(i+1):
            pdf_cut[j] = 0.0
    elif cut > 0.0:
        for i in range(len(pdf)):
            if qu[i] < cut:
                pdf_cut[i] = 0.0
                break

#--------1---------2---------3---------4---------5---------6---------7-----
# Calculation of the structure factor as the back Fourier transform of the
# given pdf
    soq = 1+2.0*np.pi**2*density/qu*getSineFT(ar,pdf_cut,qu,w=0)
#--------1---------2---------3---------4---------5---------6---------7-----

#--------1---------2---------3---------4---------5---------6---------7-----
# Pair correlation function, G(R)
    pcf  = getSineFT(qu,soq,ar,w=0)
#--------1---------2---------3---------4---------5---------6---------7-----

# Pair distribution function, g(R)
    pdf  = 1+pcf /4.0/np.pi/density/ar
# Radial distribution function, RDF(R)
    rdf  = 4.0*np.pi*density*ar*ar*pdf
# Linearised radial distribution function, T(R) = RDF(R)/R
    tor  = rdf /ar

# Running integral of the radial distribution function
    run  = np.zeros(len(ar))
    for i in range(1,len(ar)):
        run[i]  = integrate.simps(rdf[0:i], ar[0:i])
    return soq,pcf,pdf,rdf,tor,run



#--------1---------2---------3---------4---------5---------6---------7---------
def Lorentzian(x,A,q0,gamma,bckg):
    """
Lorentzian function
    
This is a Lorentzian function.

    Lor(x) = bckg + A * gamma**2/4/((x-q0)**2+gamma**2/4)

Use:
    array = Lorentzian(x,A,q0,gamma,bckg)

Input:
    - x: a range of x
    - A: the scale of the Lorentzian
    - q0: the centre of the Lorentizan
    - gamma: the width of the Lorentzian
    - bckg: an additive constant

Output:
    - A list with the values of the function
    
Created on Wed Dec 30 2020
@author: Gabriel Cuello
---------------------------------------
    """
    q = np.array(x)
    lor = A * gamma*gamma /4.0 /((q-q0)*(q-q0)+gamma*gamma/4.0) + bckg
    return lor
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def Gaussian(x,A,x0,sigma,bckg):
    """
Gaussian function
    
This is a Gaussian function.

    Gau(x) = bckg + A /sqrt(2pi)/sigma * exp(-(x-x0)^2/2/sigma^2)

Use:
    array = Gaussian(x,A,x0,sigma,bckg)

Input:
    - x: a range of x
    - A: the area of the Gaussian
    - x0: the centre of the Gaussian
    - sigma: the width of the Gaussian
    - bckg: an additive constant

Output:
    - A list with the values of the function
    
Created on Sat Dec 15 2022
@author: Gabriel Cuello
---------------------------------------
    """
    s = np.array(x)
    gau = A/np.sqrt(2*np.pi)/sigma * np.exp(-(s-x0)**2/2.0/sigma**2) + bckg
    return gau
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def LorGau(x,f0=1.0,eta=0.5,sigma=2.0,gamma=2.0,bckg=0.0):
    """
LorGau function
    
This is bell function centred in 0 and obtained as a linear combination of 
Gaussian and Lorentzian functions.

    LorGau(x) = background + factor * ( eta * Gau(x) + (1-eta) * Lor(x) )
    
        Gau(x) = exp(-x**2/2/sigma**2)
        Lor(x) = gamma**2/4/(x**2+gamma**2/4)
    
Note that LorGau(0) = background + factor

Use:
    array = LorGau(x,factor,eta,sigma,gamma,background)

Input:
    - x: a range of x
    - factor: a multiplicative constant. It is the value of the function 
      at the origin if the background is 0.
    - eta: the weight of the Gaussian contribution in [0,1]
    - sigma: the width of the Gaussian
    - gamma: the width of the Lorentzian
    - background: an additive constant

Output:
    - A list with the values of the function
    
Created on Wed Dec 30 2020
@author: Gabriel Cuello
---------------------------------------
    """
    q = np.array(x)
    lor = gamma*gamma /4.0 /(q*q+gamma*gamma/4.0)
    gau = np.exp(-1.0*q*q/2.0/sigma/sigma)
#    free = A*A/(1.0+A)**2 
#    return (f0 - free) * (eta * gau + (1.0-eta) * lor) + f0*free
    return f0 * (eta * gau + (1.0-eta) * lor) + bckg
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def GaussianA(x,A,x0,sigma,asym=0):
    """
Asymmetric Gaussian function
    
This is an asymmetric Gaussian function with a multiplying exponential

    Gau(x) = A /sqrt(2pi)/sigma * exp(-((x-x0)**2)/2.0/sigma**2) * exp(-asym*(x-x0))

Use:
    array = GaussianA(x,A,x0,sigma,asym)

Input:
    - x: a range of x
    - A: the area of the Gaussian
    - x0: the centre of the Gaussian
    - sigma: the width of the Gaussian
    - asym: exponential rate to take into account some asymmetry

Output:
    - A list with the values of the function

Remarks:
    Here the asymmetry is accounted by a multiplying exponential 'centred' in x0
    with a single parameter, i.e., exp(-asym * (x-x0)).
    
    This does not work very well, so finally it is better to use asym=0, converting
    the function in a classical Gaussian. Not that asym=0 is the default value.
    
    
Created on Dec 27 2022
@author: Gabriel Cuello
---------------------------------------
    """
    s = np.array(x)
    gau = A/np.sqrt(2*np.pi)/sigma * np.exp(-((s-x0)**2)/2.0/sigma**2) *np.exp(-asym*(s-x0))
    return gau
#--------1---------2---------3---------4---------5---------6---------7---------


#--------1---------2---------3---------4---------5---------6---------7---------
def niPeaks10(x,I0,slope,quad,wavelength,twotheta0,
            A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,
            G0,G1,G2,G3,G4,G5,G6,G7,G8,G9,
            S0,S1,S2,S3,S4,S5,S6,S7,S8,S9):
    """
    This function produces a modelised diffractogram for a nickel powder with
    the firts 10 reflections of the fcc structure.
    """
# The centroids of each peak is taken from the list of reflections of a fcc lattice
# Note the value of the lattice constant (a=3.52024 Å) that sould not be changed unless
# having very good reasons to do so.
    C0,C1,C2,C3,C4,C5,C6,C7,C8,C9 = reflections_fcc(wavelength, twotheta0, lattice=3.52024)[:10]
#    C0,C1,C2,C3,C4,C5,C6,C7,C8,C9 = (centre[0],centre[1],centre[2],centre[3],centre[4],
#                                     centre[5],centre[6],centre[7],centre[8],centre[9])

# Here the modelised diffractogram is created adding to the background 
# 10 asymmetric Gaussian functions
    diff = I0+slope*x+quad*x*x+(GaussianA(x,A0,C0,S0,G0)+GaussianA(x,A1,C1,S1,G1)+
                       GaussianA(x,A2,C2,S2,G2)+GaussianA(x,A3,C3,S3,G3)+
                       GaussianA(x,A4,C4,S4,G4)+GaussianA(x,A5,C5,S5,G5)+
                       GaussianA(x,A6,C6,S6,G6)+GaussianA(x,A7,C7,S7,G7)+
                       GaussianA(x,A8,C8,S8,G8)+GaussianA(x,A9,C9,S9,G9))
    return diff
#--------1---------2---------3---------4---------5---------6---------7---------


#--------1---------2---------3---------4---------5---------6---------7---------
def polyQ2(x,a0,a1,a2):
    '''
polyQ2 function
    
This is a 2nd-degree polynomial.

The first argument is a list with the independent variable x.
The other arguments are the coefficients of the different powers of x.

Usage:
    y = polyQ2(x,a0,a1,a2)

Input:
    - x: a range of x
    - a_i: polynomial coefficients

Output:
    - A list with the values of the function
    
    
Created on Wed Dec 30 23:47:15 2020
@author: Gabriel Cuello
---------------------------------------
    '''
    q = np.array(x)
    return (a2 * q + a1 ) * q + a0
#--------1---------2---------3---------4---------5---------6---------7---------
    
#--------1---------2---------3---------4---------5---------6---------7---------
def polyQ4(x,a0,a1,a2,a3,a4):
    """
polyQ4 function
    
This is 4th-degree polynomial.

    polyQ4 = a0 + (a1 + (a2 + (a3 + a4 * q) * q) * q) * q
    
Use:
    y = polyQ4(x,a0,a1,a2,a3,a4)

Input:
    - x: a range of x (a list)
    - a_i: polynomial coefficients

Output:
    - A list with the values of the function
    
Created on Wed Dec 30 2020
@author: Gabriel Cuello
---------------------------------------
    """
    q = np.array(x)
#    a0 + a1 * q + a2 * q*q + a3 * q*q*q + a4 * q*q*q*q
    return a0 + (a1 + (a2 + (a3 + a4 * q) * q) * q) * q
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def vanaQdep(x,a0,a1,a2,A=51.0,lowQ=0.4,Q0=7.4,dQ=2.4):
    """
vanaQdep function
    
This is a function that describes the general behaviour of the incoherent 
Vanadium signal as function of the momentum transfer (Q).

    vanadium(Q) = polyQ2 (Q,a0,a1,a2) * inelastic(Q,A,lowQ,Q0,dQ)
    
The first factor takes into account the instrument related effects (resolution) 
and the second one accounts for the Q dependence of the inelasticity.
    
Use:
    y = vanadium(x,A,lowQ,Q0,dQ,a0,a1,a2)

Input:
    - x: a range of x
    - a_i: polynomial coefficients
    - A,lowQ,Q0,dQ: parameters of the sigmoidal function (see help for the 
      inelastic function)

Output:
    - A list with the values of the function

Requires:
    - polyQ2
    - inelastic

Created on Wed Dec 30 2020
@author: Gabriel Cuello
---------------------------------------
    """
    q = np.array(x)
    polynomial = polyQ2(q,a0,a1,a2)        # Polynomial contribution
    sigmoidal = inelastic(q,A=A,lowQ=lowQ,Q0=Q0,dQ=dQ) # Inelasticity effect
    return sigmoidal*polynomial
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
        average += value * magnitude[key]
    return average

#--------1---------2---------3---------4---------5---------6---------7---------
def ratio(y1,e1,y2,e2):
    y_rat = []
    e_rat = []
    for i in range(len(y1)):
        if (y2[i] != 0):
            y_rat.append(y1[i]/y2[i])
            e_rat.append(np.sqrt(y2[i]**2 * e1[i]**2 + y1[i]**2 * e2[i]**2)/y2[i]**2)
        else:
            y_rat.append(0.0)
            e_rat.append(0.0)
    return y_rat,e_rat
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def saveFile_3col(filename,data,heading):
    """
saveFile_3col function

This function creates a 3 column ASCII file, with x,y,error.

As first line of the heading, it writes the filename preceded by "# "
Then, it prints as many lines as elements contains the list heading.
Finally, it writes a line for each point with (x,y,error)
    
Input:
    - filename: a string containing the output filename
    - x,y,e: 3 lists with the same number of elements, containing abcissa, ordinate and error
    - heading: A list where each element will be a line in the heading of the output file
    
Use:
    saveFile_3col(filename,x,y,e,heading)
    
Created on Thu Oct 19 2021
@author: Gabriel Cuello
---------------------------------------
    """
    x = data[:,0]
    y = data[:,1]
    e = data[:,2]
    with open(filename,'w') as datafile:
        datafile.write("# "+filename+'\n')
        for i in range(len(heading)):
            datafile.write("# "+heading[i]+'\n')
        for i in range(len(x)):
            datafile.write("{: 9.3f}".format(x[i])+' '+
                           "{:18.6f}".format(y[i])+' '+
                           "{:18.6f}".format(e[i])+'\n')
        print('File '+filename+' saved')
    return
#--------1---------2---------3---------4---------5---------6---------7---------


#--------1---------2---------3---------4---------5---------6---------7---------
def saveFile_2col(filename,x,y,heading):
    """
saveFile_3col function

This function creates a 2 column ASCII file, with x,y.

As first line of the heading, it writes the filename preceded by "# "
Then, it prints as many lines as elements contains the list heading.
Finally, it writes a line for each point with (x,y)
    
Input:
    - filename: a string containing the output filename
    - x,y: 2 lists with the same number of elements, containing abcissa, ordinate and error
    - heading: A list where each element will be a line in the heading of the output file
    
Use:
    saveFile_3col(filename,x,y,e,heading)
    
Created on Sat Dec 24 2022
@author: Gabriel Cuello
---------------------------------------
    """
    with open(filename,'w') as datafile:
        datafile.write("# "+filename+'\n')
        for i in range(len(heading)):
            datafile.write("# "+heading[i]+'\n')
        for i in range(len(x)):
            datafile.write("{: 9.3f}".format(x[i])+' '+
                           "{:18.6f}".format(y[i])+'\n')
        print('File '+filename+' saved')
    return
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def saveCorrelations(filename,ar,pcr,pdf,rdf,tor,run,heading):
    with open(filename,'w') as datafile:
        datafile.write("# "+filename+'\n')
        for i in range(len(heading)):
            datafile.write("# "+heading[i]+'\n')
        for i in range(len(ar)):
            datafile.write("{: 9.3f}".format(ar[i])+' '+
                           "{:12.6f}".format(pcr[i])+' '+
                           "{:12.6f}".format(pdf[i])+' '+
                           "{:12.6f}".format(rdf[i])+' '+
                           "{:12.6f}".format(tor[i])+' '+
                           "{:12.6f}".format(run[i])+'\n')
        print('File '+filename+' saved')
    return
#--------1---------2---------3---------4---------5---------6---------7---------


#--------1---------2---------3---------4---------5---------6---------7---------
def saveRSCF(filename,fou,heading):
    """
saveRSF function

This function creates a 6 column ASCII file, with the real space correlation functions

As first line of the heading, it writes the filename preceded by "# "
Then, it prints as many lines as elements contains the list heading.
Finally, it writes a line for each point with (x,y,error)
    
Input:
    - filename: a string containing the output filename
    - x,y1,y2,y3,y4,y5: 6 lists with the same number of elements, containing r and the functions
    - heading: A list where each element will be a line in the heading of the output file
    
Use:
    saveRSCF(filename,x,y1,y2,y3,y4,y5,heading)
    
Created on Mon Dec 13 2021
@author: Gabriel Cuello
---------------------------------------
    """
    with open(filename,'w') as datafile:
        datafile.write("# "+filename+'\n')
        for i in range(len(heading)):
            datafile.write("# "+heading[i]+'\n')
        for i in range(len(fou)):
            datafile.write("{: 9.3f}".format(fou[i][0])+' '+
                           "{:12.6f}".format(fou[i][1])+' '+
                           "{:12.6f}".format(fou[i][2])+' '+
                           "{:12.6f}".format(fou[i][3])+' '+
                           "{:12.6f}".format(fou[i][4])+' '+
                           "{:12.6f}".format(fou[i][5])+'\n')
        print('File '+filename+' saved')
    return
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def wsum2(w1,data1,w2,data2):
    '''
wsum2 function
    
This function sums (or subtracts) two set of values (y and error).
Depending on the parameters w1 and w2, different operations are done.
w1 and w2 are the weights for the sets 1 (y1,e1) and 2 (y2,e2).
    
    Case (1): both weights are zero (w1=w2=0)
        The function returns the sum weighted with the inveres of the square errors.
    Case (2): w1=0 and w2!=0
        The function returns the set 2 multiplied by w2
    Case (3): w1!=0 and w2=0
        If w1 is not in the range (0,1], error.
        Otherwise, the function returns the weighted sum with w2=1-w1
    Case (4): Both weights are different from zero (w1!=0 and w2!=0)
        The functions returns simply w1*y1+w2*y2.
        Note that the weights can also be negative, so this case includes the subtraction.
    
Use:
    ysum,esum = wsum2(w1,y1,e1,w2,y2,e2)

Requires:   numpy

Created on Thu Jan 01 2021
@author: Gabriel Cuello
---------------------------------------
    '''
    x1 = data1[:,0]
    y1 = data1[:,1]
    e1 = data1[:,2]
#    x2 = data1[:,0]
    y2 = data2[:,1]
    e2 = data2[:,2]
    if (len(y1)!=len(y2)):
        print('--- Error in the binary sum (wsum2).')
        print('--- The input vectors have not the same length.')
        return
    length = len(y1)
    ysum = []
    esum = []
    if ((w1==0) and (w2==0)):
        for i in range(length):
            w = 0
            sqerr = e1[i]**2+e2[i]**2
            if (sqerr!=0): w = e2[i]**2/sqerr
            ysum.append(w*y1[i]+(1-w)*y2[i])
            esum.append(np.sqrt(w**2*e1[i]**2+(1.0-w)**2*e2[i]**2))
    elif (w1==0):
        for i in range(length):
            ysum.append(w2*y2[i])
            esum.append(np.sqrt(w2**2*e2[i]**2))
    elif (w2==0):
        if ((w1>0) and (w1<=1)):
            for i in range(length):
                ysum.append(w1*y1[i]+(1-w1)*y2[i])
                esum.append(np.sqrt(w1**2*e1[i]**2+(1.0-w1)**2*e2[i]**2))
        else:
            print('--- Error in the binary sum (wsum2).')
            print('--- The the weight of first set of data should be between 0 and 1.')
    else: # both weights are different of zero and are simply considered as factors
        for i in range(length):
            ysum.append(w1*y1[i]+w2*y2[i])
            esum.append(np.sqrt(w1**2*e1[i]**2+w2**2*e2[i]**2))
    ysum = np.array(ysum)
    esum = np.array(esum)
    x1 = x1.reshape(x1.shape[0],1)
    ysum = ysum.reshape(ysum.shape[0],1)
    esum = esum.reshape(esum.shape[0],1)
    summed = np.concatenate((x1,ysum),axis=1)
    summed = np.concatenate((summed,esum),axis=1)
    return summed
# End of wsum2
#--------1---------2---------3---------4---------5---------6---------7---------


def fittingRange(xmin,xmax,ymin,ymax,abcissas,ordinates,errors):
    if (len(abcissas) != len(ordinates)):
        print ('ERROR: Ordinate and abcissa must have the same length in function fittingRange.')
    if (xmax<xmin) or (ymax<ymin):
        print ('ERROR: The limits of the rectangle are wrong in the function fittingRange.')
    x_range = []
    y_range = []
    e_range = []
    for i in range(len(abcissas)):
        if (abcissas[i] >= xmin) and (abcissas[i] <= xmax):
            if (ordinates[i] >= ymin) and (ordinates[i] <= ymax):
                x_range.append(abcissas[i])
                y_range.append(ordinates[i])
                e_range.append(errors[i])
    # np.array converts a list in a numpy array
    if (len(x_range)<10):
        print ('WARNING: Less than 10 points for the fitting in fittingRange.')
    if (len(x_range)==0):
        print ('ERROR: No points for the fitting in fittingRange.')

    x = np.array(x_range)
    y = np.array(y_range)
    e = np.array(e_range)
    x = x.reshape(x.shape[0],1)
    y = y.reshape(y.shape[0],1)
    e = e.reshape(e.shape[0],1)
#    limited = np.concatenate((x,y),axis=1)
#    limited = np.concatenate((limited,e),axis=1)
    return x,y,e


#--------1---------2---------3---------4---------5---------6---------7---------
def fit_range(xmin,xmax,ymin,ymax,data):
    """
fit_range function

This function defines a rectangle in the x,y plane. All points inside this rectangle will be used
in the fitting procedure.

Input:
    - xmin,xmax,ymin,ymax: values defining a rectangle in the x,y plane.
    - x_dat,y_dat: data points

Output:
    - x_range,y_range: two lists containing the points falling in the rectangle


Use:
    x_range,y_range = fit_range(xmin,xmax,ymin,ymax,x_dat,y_dat)

Created on Thu Jan 01 2021
@author: Gabriel Cuello
---------------------------------------
    """
    x_dat = data[:,0]
    y_dat = data[:,1]
    e_dat = data[:,2]
    if (len(x_dat) != len(y_dat)):
        print ('ERROR: Ordinate and abcissa must have the same length in the function fit_range.')
    if (xmax<xmin) or (ymax<ymin):
        print ('ERROR: The limits of the rectangle are wrong in the function fit_range.')
    x_range = []
    y_range = []
    e_range = []
    for i in range(len(x_dat)):
        if (x_dat[i] >= xmin) and (x_dat[i] <= xmax):
            if (y_dat[i] >= ymin) and (y_dat[i] <= ymax):
                x_range.append(x_dat[i])
                y_range.append(y_dat[i])
                e_range.append(e_dat[i])
    # np.array converts a list in a numpy array
    if (len(x_range)<10 or len(y_range)<0):
        print ('WARNING: Less than 10 points for the fitting in the function fit_range.')
    if (len(x_range)==0 or len(y_range)==0):
        print ('ERROR: No points for the fitting in the function fit_range.')

    x = np.array(x_range)
    y = np.array(y_range)
    e = np.array(e_range)
    x = x.reshape(x.shape[0],1)
    y = y.reshape(y.shape[0],1)
    e = e.reshape(e.shape[0],1)
    limited = np.concatenate((x,y),axis=1)
    limited = np.concatenate((limited,e),axis=1)
    return limited
# End of fit_range
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def get_xlim(xmin,xmax,dbin):
    '''
get_xlim
    
Given minimum, maximum and step values, this function produces a binning.
The function returns an integer number (nbins) with the number of bins, and two lists:
    (1) A list with the limiting values for each bin (dimension 1 x (nbins+1)).
    (2) A list with the value corresponding to each bin (dimension 1 x nbins). 
        These values correspond to the center of the bins.
    
Use:     
    number_of_bins,x_lim,x_bin = get_xlim(xmin,xmax,dbin)

Called by:  rebin

Created on Wed Dec 30, 2020
@author: Gabriel Cuello
---------------------------------------
    '''
    xini = xmin-dbin/2.0 # coordinate of the left side of the first bin
    xfin = xmax+dbin/2.0 # coordinate of the right side of the last bin
    nb = int((xfin-xini)/dbin + 1) # number of bins
    # Initialising the two output lists
    x_lim = []
    x_bin = []
    for i in range(nb+1):
        x_lim.append(xmin-dbin/2.0+i*dbin)
        x_bin.append(xmin+i*dbin)
#    x_bin.pop() # Removing the last element of this list
    return nb,x_lim,x_bin
# End of det_xlim
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def get_bins(x,xdel,x_lim):
    '''
get_bins function
    
Given an experimental value (x), the width that each point represents (xdel) and the list with the 
limiting values for each bin (x_lim) this function returns two lists:
   (1) A list of the bins covered (completely or partially) by the rectangle
   (2) A list containing the fraction of the rectangle covering each bin
Both lists have the same variable dimension.
    
Use:
    covered_bins,covered_fractions = get_bins(x,xdel,x_lim)

Called by:  rebin

Created on Wed Dec 30, 2020
@author: Gabriel Cuello
---------------------------------------
    '''
    # Initialising the two output lists
    bins = []
    frac = []
    # Determines the size of a bin, which is constant
    dbin = x_lim[1]-x_lim[0]

    # The rectangle corresponding to the x value has a width of xdel, 
    # half on each side of x
    x1 = x-xdel/2.0 # coordinate of the left side  of the rectangle
    x2 = x+xdel/2.0 # coordinate of the right side of the rectangle

    # Determining the bins where left and right sides of the rectangle fall
    b1 = int((x1-x_lim[0])/dbin)
    b2 = int((x2-x_lim[0])/dbin)
    if b1 < 0 or b2 >= len(x_lim):
        return bins,frac
    deltab = b2-b1
    # There are 3 possible cases:
    #    (1) deltab = 0
    #        The rectangle falls completely in a single bin.
    #    (2) deltab = 1
    #        The rectangle falls on 2 bins, covering partially each one.
    #    (3) deltab > 1
    #        The rectangle falls on more than 2 bins, covering partially the 
    #        first and the last ones, and completely the bins in between.
    if deltab == 0:         # Case (1)
    #   b1 (=b2) is the single bin where the rectangle falls.
    #   Then the fraction is equal to 1.0
        bins.append(b1)
        frac.append(1.0)
    elif deltab == 1:       # Case (2)
        f1 = (x_lim[b1+1]-x1)/xdel
        bins.append(b1)
        frac.append(f1)
        bins.append(b1+1)
        frac.append(1.0-f1)
    elif deltab > 1:        # Case (3)
        f1 = (x_lim[b1+1]-x1)/xdel # First bin
        bins.append(b1)
        frac.append(f1)
        for i in range(1,deltab): # Intermediate bins
            bins.append(b1+i)
            frac.append(dbin/xdel)
        f2 = (x2-x_lim[b2])/xdel # Last bin
        bins.append(b2)
        frac.append(f2)
    else:
        print ('ERROR in get_bins')

    # # There are 3 possible cases:
    # #    (1) The rectangle falls completely in a single bin
    # #    (2) The rectangle falls on 2 bins, covering partially each one
    # #    (3) The rectangle falls more than 2 bins, covering partially the 
    # #        first and the last one, and completely the bins in between.
    # if (b1 == b2):        # Case (1)
    #     bins.append(b1)   # Bin where the rectangle falls
    #     frac.append(1.0)  # Obviously, the fraction of the rectangle on this bin is 1
    # else:                 # Cases (2) and (3)
    #     fb1 = (x_lim[b1+1] - x1)/xdel   # Fraction of the rectangle falling on the first bin 
    #     bins.append(b1)
    #     frac.append(fb1)
    #     if (b2-b1 > 1):   # Only case (3) 
    #         # Loop running over all the bins completely covered by the rectangle
    #         for i in range(1,b2-b1):    
    #             bins.append(b1+i)
    #             # Obviously, the fraction is the fraction of the bin width to the rectangle width
    #             frac.append(dbin/xdel)  
    #     fb2 = (x2 - x_lim[b2])/xdel     # Fraction of the rectangle falling on the last bin 
    #     bins.append(b2)
    #     frac.append(fb2)
    return bins,frac
# End of get_bins
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def rebin(xdel,wlength,data,xmin,xmax,dbin):
    '''
rebin function

This function makes a rebinning of the experimetal data.

The rebinning can be made in angular or Q-scale, depending on the input parameter wavelength.
    A non-positive wavelength produces a binning in scattering angle.
    A positive wavelength produces a binning in Q-scale. The wavelength must be in Å.
    
The parameter xdel is the width of a channel in scattering angle, assumed constant for the 
whole diffractogram. 
At D4, having 64 channels covering 8 degrees, this value is 0.125 degrees.
    
Then, 3 list containing the experimental data follow: x_dat, y_dat, e_dat
    x_dat: abcissa in degrees if angular scale or in 1/Å if Q-scale
    y_dat: intensity in arbitraty units
    e_dat: intensity error in the same units as y_dat
    
Finally, 3 floats corresponding to the binninig: minimum, maximum and step
The units depends on the scale of the experimental data: 
    degrees for angular scale and 1/Å for Q-scale
    
The output are 3 lists containing the binned data: x_bin, y_bin, e_bin
    
Use:
    x_bin, y_bin, e_bin = rebin(xdel,wavelength,x_dat,y_dat,e_dat,xmin,xmax,dbin)

Requires:
    get_xlim, get_bins, numpy

Called by:
    rebin
    
Created on Wed Dec 30, 2020
@author: Gabriel Cuello
---------------------------------------
    '''
    x_dat = data[:,0]
    y_dat = data[:,1]
    e_dat = data[:,2]
    # Call get_xlim to obtain the number of bins, the limiting values for the bins and the values 
    # of the bins
    nbins,x_lim,x_bin = get_xlim(xmin,xmax,dbin)
    
    # Creates lists for storing the new y, new error and the fraction
    y_bin = []
    e_bin = []
    f_bin = []
    
    # Initialise the lists that will serve as accumulators
    for i in range(nbins+1): 
        y_bin.append(0.0)
        e_bin.append(0.0)
        f_bin.append(0.0)

    if (wlength <= 0):
        # The binning is in angular scale
        for i in range(len(x_dat)):
            if (np.isnan(y_dat[i]) == False) and (np.isnan(e_dat[i]) == False):
                # For each experimental x value, calls get_bins, which returns the bins covered by that 
                # point and the corresponding fractions
                bins,frac = get_bins(x_dat[i],xdel,x_lim)
                # For each of these bins we add the corresponding fraction of y and error.
                # Because these fractions act as weighting factors, the fractions are accumulated for 
                # normalisation purposes
                for j in range(len(bins)):
                    y_bin[bins[j]] += frac[j]*y_dat[i]
                    e_bin[bins[j]] += frac[j]*e_dat[i]
                    f_bin[bins[j]] += frac[j]
    else:
        # The binning is in Q-scale
        for i in range(len(x_dat)):
            if (np.isnan(y_dat[i]) == False) and (np.isnan(e_dat[i]) == False):
                # ${\rm d}Q = \frac{2\pi}{\lambda} \sqrt{1-\frac{Q\lambda}{4 \pi}}
                #  {\rm d}2\theta \frac{pi}{180}$
                qdel = 2.0*np.pi/wlength * np.sqrt(1.0-(x_dat[i]*wlength/4.0/np.pi)**2)
                qdel *= xdel * np.pi/180.0 
                # For each experimental x value, calls get_bins, which returns the bins covered by that 
                # point and the corresponding fractions
                bins,frac = get_bins(x_dat[i],qdel,x_lim)
    #            print (bins,frac)
                # For each of these bins we add the corresponding fraction of y and error.
                # Because these fractions act as weighting factors, the fractions are accumulated for 
                # normalisation purposes
    #            print ('--->', len(bins),y_bin)
                for j in range(len(bins)):
                    # print ('i=',i,'j=',j,'len(bins)=',len(bins),
                    #        'len(y_bin)=',len(y_bin),'y_bin[bins[j]]=',y_bin[bins[j]],'bins[j]=',bins[j],
                    #        'frac[j]=',frac[j],
                    #        'y_dat[i]=',y_dat[i],'len(y_dat)=',len(y_dat),
                    #        'len(x_dat)=',len(x_dat)
                    #        )
                    y_bin[bins[j]] += frac[j]*y_dat[i]
                    e_bin[bins[j]] += frac[j]*e_dat[i]
                    f_bin[bins[j]] += frac[j]
        
    # Normalisation of y and errors. If fraction is 0, that bin has no data.
    for i in range(nbins+1):
        if (f_bin[i] != 0):
            y_bin[i] /= f_bin[i]
            e_bin[i] /= f_bin[i]
#            print(x_bin[i],y_bin[i],e_bin[i])
    x_bin = np.array(x_bin)
    y_bin = np.array(y_bin)
    e_bin = np.array(e_bin)
#   Reshapes the arrays as 2D arrays, but with 1 column
    x_bin=x_bin.reshape(x_bin.shape[0],1)
    y_bin=y_bin.reshape(y_bin.shape[0],1)
    e_bin=e_bin.reshape(e_bin.shape[0],1)
#   Concatenates the arrays to have a 3-column matrix 
    binned = np.concatenate((x_bin,y_bin),axis=1)
    binned = np.concatenate((binned,e_bin),axis=1)
    return binned
    
#    return x_bin,y_bin,e_bin
# End of rebin
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def read_3col(filename):
    '''
read_3col function

Opens a file and read it line by line. This function assumes it is a 3-columns file.
The lines containig the symbol # are ignored. Be careful, # could be anywhere in the line!
The empty lines are also ignored.
As output it produces 3 lists with abcissas, ordinates and errors for the ordinates.
    
Use:
    x_dat, y_dat, e_dat = readD4_3col('mydatafile.dat')
    
Created on Wed Dec 30, 2020
@author: Gabriel Cuello
---------------------------------------
    '''
    data = open(filename,'r') # Opens the data file in read only mode 

    # Creating the lists that will contain abcissas, ordinates and errors
    x = []
    y = []
    e = []
    
    # Reading the file line by line
    for dataline in data.readlines():
    # dataline is a string that contains each line of the file in data.
    # note that the last character of the string is a 'carriage return', \n
        if '#' not in dataline:  # Only the lines without # are treated.
            # the method .strip(' ') removes blanks at the beginning of the string
            row = dataline.strip(' ')[:-1]
            if len(row)>0:  # Only the no-empty lines are treated
                columns = row.split()   # This method split the line using the spaces
                x.append(float(columns[0]))
                y.append(float(columns[1]))
                if (len(columns)==3):
                    if (columns[2] == 'i') or (columns[2] == 'o'):
                        e.append(float(0.0))
                    else:
                        e.append(float(columns[2]))
                else:
                    e.append(float(0.0))
    data.close()            
    print ('The data file {} read with no errors. Number of data = {}'.
           format(filename,len(x)))
#    return np.array(x,y,e)
#   Converts the lists in arrays
    xa=np.array(x)
    ya=np.array(y)
    ea=np.array(e)
#   Reshapes the arrays as 2D arrays, but with 1 column
    xa=xa.reshape(xa.shape[0],1)
    ya=ya.reshape(ya.shape[0],1)
    ea=ea.reshape(ea.shape[0],1)
#   Concatenates the arrays to have a 3-column matrix 
    data = np.concatenate((xa,ya),axis=1)
    data = np.concatenate((data,ea),axis=1)
    return data
# End of read_3col
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------

