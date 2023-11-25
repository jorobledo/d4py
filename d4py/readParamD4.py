#--------1---------2---------3---------4---------5---------6---------7---------8    

print('Module readParamD4 (Matias Ilarri)')
print('Imported...')
import sys
from Element import elemento
import d4treat as d4
import numpy as np
from scipy.interpolate import interp1d

#--------1---------2---------3---------4---------5---------6---------7---------8
def generate_correction_factor(xvals, qvals, single, s_aten, total):
    '''    
This function creates the multiplying factor that corrects for attenuation and 
multiple scattering contributions.

The input arrays are produced by the Monte Carlo simulation using fewer points
that the actual number of points in the experimental Q-scale. This is possible
because these curves are smooth. Reducing the number of points in the Q-range
allows to reduce the calculation time.

Input: xvals, qvals, single, s_aten, total
    - xvals : (array),
        Experimental Q-scale.
    - qvals : (array),
        Q-scale for the MC simulation.
    - single : (array),
        Single scattering with attenuation correction from MC.
    - s_aten : (array),
        Single scattering without attenuation correction from MC.
    - total : (array),
        Total scattering result from Monte Carlo simulation.

Output:
    - correction_factor: (array) contaning the correction factor to apply to
      de total scatteing signal. The values correspond to the experimental
      Q-scale.

Author: Gabriel Cuello
Date: 20/02/2023
--------------------------------------------------------------------------------
    '''
# Multiple Scattering factor. This is the factor that multiplies the total
# scattering to get the single scattering contribution.
    MS_factor = single / total
# Attenuation factor. This is the factor that multiplies single scattering
# without any attenuation correction to get the corrected signle scattering.
    atten_factor = single / s_aten
# Here we calculated a single global factor that perform both corrections.
    global_factor = MS_factor * atten_factor
# If the value is NaN, the factor is replaced by 0
    global_factor[np.isnan(global_factor)]=0
# The function interp1d creates a function (correction_factor), which
# interpolates the factor obtained from the Monte Carlo simulation
#
    correction_factor = interp1d(qvals, global_factor,fill_value=0)
    return correction_factor(xvals)
#--------1---------2---------3---------4---------5---------6---------7---------8    


#--------1---------2---------3---------4---------5---------6---------7---------8    
def atomDisarm(comp_sample):
    '''
This function calculates the number of atoms of each element in the sample.

The sample is composed by several compounds, each one with a given chemical 
formula. Starting with the list of these compounds and their concentrations in
the sample, this function uses the chemical formulas and ddetermines the
concentration of each atom in the sample. 
    
Input:
  - comp_sample(dict)
    This is a dictionary where the keys are the chemical formulas of compounds
    and the values are the concentrations of those compounds in the sample.
        
Output:
  - totalAtomSample (dict)
    This is a dictionary where the keys are the chemical symbol of the elements
    and the values are the concentrations of those elements in the sample.

Author: Gabriel Cuello
Date: 20/02/2023
--------------------------------------------------------------------------------
    '''
# List of compounds (or components) in the sample
    component = list(comp_sample.keys())   
    atoms_sample = []
    # chemform is the chemical formula of each compound
    for chemform in component:
        element="-"
        aux="-"
        nElement=1
        for j in chemform:   
            if(j.isupper()):
                if(aux.isupper()):
                    element=aux
                    atoms_sample.append([chemform, element, nElement])
                    element=j
                    aux=j
                    nElement=1
                elif(element=="-"):
                    element=j
                else:
                    atoms_sample.append([chemform, element, nElement])
                    element=j
                    aux=j
                    nElement=1
            elif(j.islower()):
                element=element+j
                aux=j
            else:
                if(aux.isnumeric()):
                    nElement=int(aux)*10+int(j)
                else:
                    nElement=int(j)
                    aux=j
        atoms_sample.append([chemform, element, nElement])
    
    # print("----------------")
    # print(atoms_sample)   
    # print("----------------")
    
    totalAtomSample = {}
    for i in range(len(atoms_sample)):
        if atoms_sample[i][1] in totalAtomSample:
            totalAtomSample[atoms_sample[i][1]] += comp_sample[atoms_sample[i][0]]*int(atoms_sample[i][2])
        else:
            totalAtomSample[atoms_sample[i][1]] = comp_sample[atoms_sample[i][0]]*int(atoms_sample[i][2])
    # print(totalAtomSample)
    
    return totalAtomSample
#--------1---------2---------3---------4---------5---------6---------7---------8
        
#--------1---------2---------3---------4---------5---------6---------7---------8
class readParamD4:
    '''
This class read the parameters file for the data treatment.
    
Input:
  - parametersFile: a string with the filename
       
Selfs:
    'Proposal' 'Proposer' 'experimenters' 'otherUsers' 'instrument' 'LC'
    'startDate' 'endDate' 'environmentInfo' 'logBook' 'logPage'
    'zeroAngle' 'wavelength' 'LohenSlit' 'GamsSlit' 'topFlag' 'bottomFlag'
    'material' 'shape' 'outerDiam' 'innerDiam' 'height' 'AngularResolution'
    'AMin' 'AMax' 'AStep' 'QMin' 'QMax' 'QStep' 'RMin' 'RMax' 'RStep'
    'NAtoms' 'Diam' 'density' 'heightSample' 'massSample' 'densitySample'
    'titleSample' 'sample' 'container' 'environment' 'vanadium'
        
    comp_sample
        contains the proportions of each sample component
    atoms_sample
        contains the atoms of the components with their total amounts
    elements
        contains the information of the elements of the sample
    natoms_sample
        Number of atoms in a unit
    c_sample
        
Functions or methods:
  - __init__(self, parametersFile)
  - noLoad(self,var):

Author: Gabriel Cuello
Date: 20/02/2023
--------------------------------------------------------------------------------
    '''
#--------1---------2---------3---------4---------5---------6---------7---------8    
#---5----1---------2---------3---------4---------5---------6---------7---------8    
    def __init__(self, parametersFile):
        # Initialise the dictionary that will contain the parameters
        params={}
        # print(parametersFile)
            
        with open(parametersFile, "r") as parfile:
            lines = parfile.readlines()

# Here verifies that any line with more than 1 character starts with one of
# these characters: #, ! or <
        for i in range(len(lines)):
            if (len(lines[i]) > 1):
                first = ((lines[i][0] != "#") and (lines[i][0] != "!") and 
                         (lines[i][0] != "<"))
                if first:
                    print ("Wrong input in line: ",i+1," file: ",parametersFile)
                    sys.exit()
        
# Here the blank lines or those starting with # or ! are considered as comments
# i.e., they are ignored.
        for i in range(len(lines)):
            if (lines[i][0] == "#") or (lines[i][0] == "!"):
                pass
            elif len(lines[i]) == 1:
                pass
            elif (lines[i][0] == "<" ):
# If the line starts with <, it is a line with useful information
                line = lines[i].split(">")
                if(line[1].find(";")== -1):
                    print("Missing value for parameter {} \n".format(line[0][1:]))
                line2=line[1][1:].split(";")
                # print(line[0][1:]) # aca tengo el nombre de la variable       
                # print(line2[0]) # aca tengo el parametro
                params[line[0][1:]]=line2[0]
                
        comp_sample={}
        #Codigo para cargar muestras
        aux=0
        while(("comp{}".format(aux)) in params):
            auxComp=params["comp{}".format(aux)].split("=")
            comp_sample[auxComp[0]]=float(auxComp[1])
            aux=aux+1
        if aux==0: self.noLoad("comp") # component load check
            
        atoms_sample = atomDisarm(comp_sample)
        elements={}
        
        for i in atoms_sample.keys():
            modif={}
            if i in params.keys():
                line = params[i].split(" ")
                for j in line:
                    line2=j.split("=")
                    modif[line2[0]]=line2[1]
                if "isot" in modif.keys(): #change isotope
                    i=modif["isot"]
                    modif.pop("isot")
            auxElement = elemento(i, table=0, **modif)        #get the element from the element class 
            elements[auxElement.symbol]=auxElement.__dict__   #with the modifications written in the parameter file

        
        natoms_sample = d4.getNofAtoms(atoms_sample) # Number of atoms in a unit
        c_sample = {}
        for key, value in atoms_sample.items():
            elements[key]['MassNbr'] = d4.getMassNumber(float(elements[key]['weight']))
            c_sample[key] = value/natoms_sample
        
        def stringValue(key,params,defaultValue):
            value = defaultValue[key]
            if key in params.keys():
                value = params[key]
            return value
        
        self.Proposal        = stringValue("Proposal"       ,params,defaultValue)
        self.mainProposer    = stringValue("mainProposer"   ,params,defaultValue)
        self.experimenters   = stringValue("experimenters"  ,params,defaultValue)
        self.otherUsers      = stringValue("otherUsers"     ,params,defaultValue)
        self.instrument      = stringValue("instrument"     ,params,defaultValue)
        self.LC              = stringValue("LC"             ,params,defaultValue)
        self.startDate       = stringValue("startDate"      ,params,defaultValue)
        self.endDate         = stringValue("endDate"        ,params,defaultValue)
        self.environmentInfo = stringValue("environmentInfo",params,defaultValue)
        self.logBook         = stringValue("logBook"        ,params,defaultValue)
        self.logPage         = stringValue("logPage"        ,params,defaultValue)

        self.zeroAngle   = float(stringValue("zeroAngle" ,params,defaultValue))
        self.wavelength  = float(stringValue("wavelength",params,defaultValue))
        self.LohenSlit   = float(stringValue("LohenSlit" ,params,defaultValue))
        self.GamsSlit    = float(stringValue("GamsSlit"  ,params,defaultValue))
        self.topFlag     = float(stringValue("topFlag"   ,params,defaultValue))
        self.bottomFlag  = float(stringValue("bottomFlag",params,defaultValue))

        self.material    = stringValue("material" ,params,defaultValue)
        self.shape       = stringValue("shape"    ,params,defaultValue)

        self.outerDiam         = float(stringValue("outerDiam"        ,params,defaultValue))
        self.innerDiam         = float(stringValue("innerDiam"        ,params,defaultValue))
        self.height            = float(stringValue("height"           ,params,defaultValue))
        self.AngularResolution = float(stringValue("AngularResolution",params,defaultValue))

        self.AMin  = float(stringValue("AMin" ,params,defaultValue))
        self.AMax  = float(stringValue("AMax" ,params,defaultValue))
        self.AStep = float(stringValue("AStep",params,defaultValue))
        self.QMin  = float(stringValue("QMin" ,params,defaultValue))
        self.QMax  = float(stringValue("QMax" ,params,defaultValue))
        self.QStep = float(stringValue("QStep",params,defaultValue))
        self.RMin  = float(stringValue("RMin" ,params,defaultValue))
        self.RMax  = float(stringValue("RMax" ,params,defaultValue))
        self.RStep = float(stringValue("RStep",params,defaultValue))

        self.NAtoms  = float(stringValue("NAtoms" ,params,defaultValue))
        self.Diam    = float(stringValue("Diam"   ,params,defaultValue))
        self.density = float(stringValue("density",params,defaultValue))

        self.comp_sample   = comp_sample
        self.atoms_sample  = atoms_sample
        self.elements      = elements
        self.natoms_sample = natoms_sample
        self.c_sample      = c_sample
        
        self.heightSample  = float(stringValue("heightSample" ,params,defaultValue))
        self.massSample    = float(stringValue("massSample"   ,params,defaultValue))
        self.densitySample = float(stringValue("densitySample",params,defaultValue))
        self.densitySample = float(stringValue("densitySample",params,defaultValue))

        self.titleSample    = stringValue("titleSample" ,params,defaultValue)

        if "sample" in params.keys():
            self.sample=params["sample"]
        else: 
            print(" Missing sample data file ")
        if "environment" in params.keys():
            self.environment=params["environment"]
        else: 
            print(" Missing environment data file ")
        if "container" in params.keys():
            self.container=params["container"]
        else: 
            print(" Missing container data file ")
        if "vanadium" in params.keys():
            self.vanadium=params["vanadium"]
        else: 
            print(" Missing vanadium data file ")       
        
    defaultValue={'Proposal'      : (" Empty "),
                  'Proposer'      : (" Empty "),
                  'experimenters' : (" Empty "),
                  'otherUsers'    : (" Empty "),
                  'instrument'    : (" Empty "),
                  'LC'            : (" Empty "),
                  'startDate'     : (" Empty "),
                  'endDate'       : (" Empty "),
                  'environment'   : (" Empty "),
                  'logBook'       : (" Empty "),
                  'logPage'       : (" Empty "),
                  'zeroAngle' : (0.0),  'wavelength': (0.5),
                  'LohenSlit' : (6.0),  'GamsSlit'  : (-6.0),
                  'topFlag'   : (25.0), 'bottomFlag': (-25.0),
                  'material'  : ("V"),  'shape'     : ("cylinder"),
                  'outerDiam' : (5.0),  'innerDiam' : (4.8),
                  'height' : (60.0),
                  'AngularResolution' : (0.125),
                  'AMin' : (0.0), 'AMax' : (140.0), 'AStep' : (0.125),
                  'QMin' : (0.0), 'QMax' : (23.5),  'QStep' : (0.02),
                  'RMin' : (0.0), 'RMax' : (20.0),  'RStep' : (0.01),
                  'NAtoms' : (1.0), 'Diam' : (6.08), 'density' : (6.51),
                  'heightSample' : (59.5), 'massSample' : (2.1271),
                  'densitySample' : (3.708),
                  'titleSample' : (" None ")}
        return
#---5----1---------2---------3---------4---------5---------6---------7---------8    
#--------1---------2---------3---------4---------5---------6---------7---------8    

#--------1---------2---------3---------4---------5---------6---------7---------8    
    def noLoad(self,var):
        print('-'*70)
        print(" You must load the {} of the sample to be processed ".format(var))
        print('-'*70)
        return
#--------1---------2---------3---------4---------5---------6---------7---------8    


