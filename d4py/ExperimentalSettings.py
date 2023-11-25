import numpy as np                   # NumPy library
from .PhysicalMagnitudes import *
from .MiscelaneousCalculations import *

#--------1---------2---------3---------4---------5---------6---------7---------8
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
--------------------------------------------------------------------------------
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
def setNumorInfo(totalNumors=('Exp',0,1),containerNumors=('Can',0,1),
            environmentNumors=('Env',0,1),
            nickelNumors=('Ni5',0,1),vanadiumNumors=('Van',0,1),absorberNumors=('Abs',0,1),
            sampleNumors =('S01',0,1)):


    numors = {'experiment': totalNumors,'container': containerNumors,
              'environment': environmentNumors,
              'nickel': nickelNumors,'vanadium': vanadiumNumors,'absorber': absorberNumors,
              'sample': sampleNumors}

    print ('Numors')

    if totalNumors[1]!=0:
        print ('{}Total {}: {} -{}, {} numors.'.format(
            4*' ',*numors['experiment'],int(numors['experiment'][2])-int(numors['experiment'][1])+1))
    if containerNumors[1]!=0:
        print ('{}{}: {} -{}, {} numors.'.format(
            4*' ',*numors['container'],int(numors['container'][2])-int(numors['container'][1])+1))
    if environmentNumors[1]!=0:
        print ('{}{}: {} -{}, {} numors.'.format(
            4*' ',*numors['environment'],int(numors['environment'][2])-int(numors['environment'][1])+1))
    if nickelNumors[1]!=0:
        print ('{}{}: {} -{}, {} numors.'.format(
            4*' ',*numors['nickel'],int(numors['nickel'][2])-int(numors['nickel'][1])+1))
    if vanadiumNumors[1]!=0:
        print ('{}{}: {} -{}, {} numors.'.format(
            4*' ',*numors['vanadium'],int(numors['vanadium'][2])-int(numors['vanadium'][1])+1))
    if absorberNumors[1]!=0:
        print ('{}{}: {} -{}, {} numors.'.format(
            4*' ',*numors['absorber'],int(numors['absorber'][2])-int(numors['absorber'][1])+1))
    if sampleNumors[1]!=0:    
        print ('{}Sample {}: {} -{}, {} numors.'.format(
            4*' ',*numors['sample'],int(numors['sample'][2])-int(numors['sample'][1])+1))
    return numors
#--------1---------2---------3---------4---------5---------6---------7---------

#--------1---------2---------3---------4---------5---------6---------7---------
def setVanaInfo(IncXS=5.08,CohXS=0.0184,ScaXS=5.1,AbsXS=5.08,
                CohSL=-0.3824,molarM=50.9415,NAtoms=1.0,
                Diam=6.08,Height=50.0,density=6.51):
    vana = {}
    vana['IncXS'] = IncXS # Incoherent scattering cross section (sigma bound), in barns
    vana['CohXS'] = CohXS # Coherent scattering cross section, in barns
    vana['ScaXS'] = ScaXS # Coherent scattering cross section, in barns
    vana['AbsXS'] = AbsXS # Absorption cross section, in barns
    vana['CohSL'] = CohSL # Coherent scattering length, in fm
    vana['SelfQ0'] = vana['IncXS']/4./np.pi
    vana['diam'] = Diam # diameter (mm)
    vana['NAtoms'] = NAtoms # Number of atoms in a unit
    vana['molarM'] = molarM # Molar mass of vanadium, in g/mol
    vana['den_gcc'] =density # Macroscopic density, in g/cm3
    vana['den_aac'] = getAtomicDensity(density=vana['den_gcc'],molarMass=vana['molarM'])
    vana['A'] = getMassNumber(molarMass=vana['molarM']) # Ratio mass to neutron mass
    vana['FreeIncXS'] = getFreeXS(BoundXS=vana['IncXS'],A=vana['A']) 
    vana['volume'] = getCylVolume(diameter=vana['diam'],height=Height)/1000.0
    print (30*'-')
    print ('Standard of vanadium')
    print(4*' ','The standard is a cylinder of {:.3g} mm of diameter.'.format(vana['diam']))
    print(4*' ','Volume in the beam = {:.3g} cm3.'.format(vana['volume']))
    print ()
    print(4*' ','Bound incoherent cross section {:.6g} barns/atom.'.format(vana['IncXS']))
    print(4*' ','Free incoherent cross section',"{:.6g} barns/atom.".format(vana['FreeIncXS']))
    print ()
    print(4*' ','b**2 = sigma/4/pi at Q=0 (self) is',"{:.6g} barns/sterad/atom.".format(vana['SelfQ0']))
    print(4*' ','b**2 = sigma/4/pi at Q=infty is',"{:.6g} barns/sterad/atom".format(vana['FreeIncXS']/4./np.pi))
    print ()
    print(4*' ','The molar mass is',"{:.6g} g/mol.".format(vana['molarM']))
    print(4*' ','Mass number, A = ',"{:.6g} (= mass/neutron_mass)".format(vana['A']))
    print ()
    print(4*' ','Density = ',"{:.6g} g/cm3.".format(vana['den_gcc']))
    print(4*' ','Atomic density = {:.6g} atoms/A3.'.format(vana['den_aac']))
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
    sample['molarM'] = AtomicAvg(c_sample,extractAttr(elements,'weight'))
    
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
    sample['CohXS'] = AtomicAvg(c_sample,extractAttr(elements,'sig_coh')) 
    # Incoherent scattering cross section (barns)
    sample['IncXS'] = AtomicAvg(c_sample,extractAttr(elements,'sig_inc')) 
    # Scattering cross section (barns)
    sample['ScaXS'] = AtomicAvg(c_sample,extractAttr(elements,'sig_sca')) 
    # Absorption cross section (barns)
    sample['AbsXS'] = AtomicAvg(c_sample,extractAttr(elements,'sig_abs')) 
    # Absorption cross section at working wavelength (barns)
    sample['AbsWW'] = sample['AbsXS'] * wavelength/1.8 
    # Coherent scattering length (fm)
    sample['CohSL'] = AtomicAvg(c_sample,extractAttr(elements,'re_bcoh')) 
    # Atomic number (amu)
    sample['A'] = AtomicAvg(c_sample,extractAttr(elements,'A')) 
    # Free coherent scattering cross section (barns)
    sample['FreeCohXS'] = getFreeXS(BoundXS=sample['CohXS'],A=sample['A'])
    # Free incoherent scattering cross section (barns)
    sample['FreeIncXS'] = getFreeXS(BoundXS=sample['IncXS'],A=sample['A'])
    # Free scattering cross section (barns)
    sample['FreeScaXS'] = getFreeXS(BoundXS=sample['ScaXS'],A=sample['A'])
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
    print ('Mass number: {:.6g} amu'.format(sample['A']))
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
