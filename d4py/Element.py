import math
import numpy as np
print('Module Element (Maria Sanchez)')
print('Imported...')


class elemento:
    '''
    class elemento:
    
    
        Class that saves by default every nuclear value given by the "Bound Scattering Lengths and Cross Sections of the Elements and Their      Isotopes" in Appendix 2 of Neutron Scattering Fundamentals, volume 44 (dictionary at the bottom called Isotope_dict) inside an object.
    
    The user can also replace any default value with one of his choice or use this class to consult the value of any isotope variable.
    
    object = elemento( Iso_name, table, changes)
    
    Input: 
         - Iso_name : (string) Name of the isotope, written as 'Li', 'Li6', 'Li7'... If the isotope does not exist, a warning will be                printed
         
         - Table : (int) This is a optional input with a default value of 0, which corresponds with no display. Type 1 if the user wants              to have a table with all the values linked to the isotope used displayed on screen. Type 2 if you want all the tables from all            the isotope of the same element displayed. 
         
         - Changes : (string = float) This is an optional input, you can write as many as you want. If you want to use a different value              than the ones from the Isotope_dict dictionary, type the name of the variable, = , the new value. If this new value is somehow            incoherent with its physical meaning, the user will recieve a warning as a print.
         
   Titanium = elemento('Ti')         Creates an object containing the default values of natural titanium. No output
   Titanium = elemento('Ti', 1)      Creates an object containing the default values of natural titanium and displays them on screen
   Titanium = elemento('Ti', 2)      Creates an object containing the default values of natural titanium and displays one tables of                                            Titanium isotopes on screen
   
   Titanium = elemento('Ti', 1, A=5, weight = 3)   Creates an object containing the default values of natural titanium except for the                                                        atomic number A that is assigned the value 5, and the weight, replaced by 3. No output
    
   The keywords used are:
   
        symbol : (str) symbol of the isotope
        iso: (str) isotope (‘nat’, 1,2,3…)
        name: (str) name of the element
        A : (int)  mass number/ number of neutrons and protons
        Z : (int) Atomic number/ number of protons
        weight : (float) number of protons and neutrons times their masses
        excess : (float) 
        spin : (str) nuclear spin
        parity : (str) nuclear parity
        abundance: (float) nuclear abundance (%) for a stable nuclei
        life: (str) half-life for a unstable nuclei followed by its units
        re_bcoh: (float) real part of the bound coherent scattering length
        im_bcoh: (float) imaginary part of the bound coherent scattering length
        bplus: (float) bound scattering length for the (I + ½) state
        bminus: (float) bound scattering length for the (I - ½) state
        re_binc: : (float) real part of the bound incoherent scattering length
        im_binc: (float) imaginary part of the bound incoherent scattering length
        sig_coh: (float) module of bound coherent cross section (σ)
        sig_inc: (float)  module of bound incoherent cross section
        sig_sca: (float) total bound scattering cross section
        sig_abs: (float) absorption cross section for thermal (ν=2200 m/s) neutrons
        contrast: 
        neut: (int) number of neutrons, A-Z
        bcoh: (complex) total coherent scattering length 
        binc: (complex) total incoherent scattering length
        scoh_bound: (complex) bound coherent cross section

                𝑠𝑐𝑜ℎ_𝑏𝑜𝑢𝑛𝑑=4𝜋 𝑏𝑐𝑜ℎ∙𝑏𝑐𝑜ℎ'        bcoh' = conjugate

        sinc_bound: (complex) bound incoherent cross section

                𝑠𝑖𝑛𝑐_𝑏𝑜𝑢𝑛𝑑=4𝜋 𝑏𝑖𝑛𝑐∙𝑏𝑖𝑛𝑐'         binc' = conjugate
    '''

    def __init__(self, Iso_name , table =0, **kwargs):
        if ((Iso_name in isotope_dict) == False):
            print("\n The isotope you wrote does not exist in this dictionary")
            
            # dicSearch = dict(list(zip(nombres, isotope_dict[Iso_name[0]])))# take parameters from table
            # if (table ==1):
            #     # Print the names of the columns.
            #     print("{:<15}".format(Iso_name))
            #     print("{:<10} {:<10} {:<10}".format('PARAMETER', 'VALUE', 'UNITS'))
            #     # print each data item.
            #     for key in dicSearch:
            #         print("{:<10} {:<10} {:<10}".format(key, dicSearch[key], units[key]))
            # if (table==2):
            #     for key in isotope_dict:
            #         if isotope_dict[key][0] == isotope_dict[Iso_name][0]:
            #             doc = dict(list(zip(nombres, isotope_dict[key])))
            #             print("\n")
            #             print("{:<15}".format(key.upper()))
            #             print("\n")
            #             print("{:<10} {:<10} {:<10}".format('PARAMETER', 'VALUE', 'UNITS'))
            #             for key in doc:
            #                 print("{:<10} {:<10} {:<10}".format(key, dicSearch[key], units[key]))

        else:
            dic = dict(list(zip(nombres, isotope_dict[Iso_name])))# take parameters from table
            for key in dic:
                if dic[key] == 'NULL': # replaces null with zero in all fields 
                    dic[key] = 0
            dic['bcoh'] = 0
            dic['binc'] = 0
            dic['scoh_bound'] = 0
            dic['sinc_bound'] = 0
            dic['neut']= 0
            for key, value in kwargs.items():
                if (key != 'b_length' and key!= 's_bound' and key!= 's_free'):
                    # print("\n")
                    for element in myCorr[key]:
                        # print("\n")
                        element(key, float(value))
                    # print("\n")
                    dic[key] = value

            bcoh = complex(float(dic['re_bcoh']), float(dic['im_bcoh']))
            binc = complex(float(dic['re_binc']), float(dic['im_binc']))
            dic['bcoh'] = bcoh
            dic['binc'] = binc

            scoh_bound = 4 * np.pi * bcoh * bcoh.conjugate()
            sinc_bound = 4 * np.pi * binc * binc.conjugate()

            dic['scoh_bound'] = complex(round(scoh_bound.real, 4),round(scoh_bound.imag, 4))
            dic['sinc_bound'] = complex(round(sinc_bound.real, 4), round(sinc_bound.imag, 4))
            dic['neut'] = dic['A']-dic['Z']

            for key, value in kwargs.items():
                if (key == ('b_length' or 's_bound' or 's_free')) :
                    for element in myCorr[key]:
                        element(key, float(value))
                    dic[key]= value
            if (table ==1):
                # Print the names of the columns.
                print("{:<15}".format(Iso_name))
                print("{:<10} {:<10} {:<10}".format('PARAMETER', 'VALUE', 'UNITS'))
                # print each data item.
                for key in dic:
                    print("{:<10} {:<10} {:<10}".format(key, dic[key], units[key]))
            if (table==2):
                for key in isotope_dict:
                    if isotope_dict[key][0] == isotope_dict[Iso_name][0]:
                        doc = dict(list(zip(nombres, isotope_dict[key])))
                        print("\n")
                        print("{:<15}".format(key.upper()))
                        print("\n")
                        print("{:<10} {:<10} {:<10}".format('PARAMETER', 'VALUE', 'UNITS'))
                        for key in doc:
                            print("{:<10} {:<10} {:<10}".format(key, dic[key], units[key]))

            self.symbol = dic['symbol']
            self.iso = dic['iso']
            self.name = dic['name']
            self.A = dic['A']
            self.Z = dic['Z']
            self.weight = dic['weight']
            self.excess = dic['excess']
            self.spin = dic['spin']
            self.parity = dic['parity']
            self.abundance = dic['abundance']
            self.life = dic['life']
            self.re_bcoh = dic['re_bcoh']
            self.im_bcoh = dic['im_bcoh']
            self.bplus = dic['bplus']
            self.bminus = dic['bminus']
            self.re_binc = dic['re_binc']
            self.im_binc = dic['im_binc']
            self.sig_coh = dic['sig_coh']
            self.sig_inc = dic['sig_inc']
            self.sig_sca = dic['sig_sca']
            self.sig_abs = dic['sig_abs']
            self.contrast = dic['contrast']
            self.neut = dic['neut']
            self.scoh_bound = dic['scoh_bound']
            self.sinc_bound = dic['sinc_bound']
            self.bcoh = dic['bcoh']
            self.binc = dic['binc']

    @classmethod
    def getIsotope(self, Iso_name):
        lista=[]
        
        # print(name)
        if Iso_name=="All":
            xName=""
            for key in isotope_dict:                
                if isotope_dict[key][2] != xName:
                    lista.append(isotope_dict[key][2])
                xName=isotope_dict[key][2]
        else:
            for key in isotope_dict:
                if isotope_dict[key][2] == Iso_name:
                    if(isotope_dict[key][1]=="nat"):
                        lista.append(isotope_dict[key][0])
                    else:
                        lista.append(isotope_dict[key][0]+isotope_dict[key][1])
        return lista
        
def Positive(key, number):
    if((float(number) > 0) == False):
        print(" WARNING. ", key, " parameter must have a positive value")
def Negative(key, number):
    if ((float(number) < 0) == False):
        print(" WARNING.", key, "parameter must have a negative value")
def SemiWhole(key, number):
    if((float(number)%(1/2)) != 0):
        print(" WARNING.", key, "parameter must be semi-whole")
def Top(key, value):
    if (key == 'weight'): top = 248.072
    if (key == 'abundance'): top = 100
    if (key == 'Z'): top = 96
    if(value > top):
        print(" WARNING.", key, "parameter can't be bigger than", top)
def Int(key, value):
    if (value.is_integer() == False):
        print(" WARNING.", key, "parameter has to be a whole number")
def String(key, value):
    if (isinstance(value, str) == False):
        print(" WARNING.", key, "parameter has to be a string")
def Complex(key, value):
    if (isinstance(value, complex) == False):
        print(" WARNING.", key, "parameter has to be a complex number")
def getFree(isot,key):
    s = complex(round(math.pow((isot.A/(isot.A +1)), 2) * key.real, 4), round(math.pow((isot.A/(isot.A +1)), 2) * key.imag, 4))
    return s
def dummyOk(key, value):
    return

myCorr = {
    "symbol": [String],
    "iso": [String],
    "name": [String],
    "A": [Int, Positive],
    "Z": [Int,Positive],
    "weight": [Top, Positive],
    #"excess": ,
    "spin": [SemiWhole],
    #"parity": ,
    "abundance": [Top, Positive],
    "life": [Positive, String],
    "re_bcoh": [dummyOk],
    #"im_bcoh": ,
    #"bplus": ,
    #"bminus": ,
    #"re_binc": ,
    #"im_binc": ,
    #"sig_coh": ,
    #"sig_inc": ,
    #"sig_sca": ,
    #"sig_abs": ,
    #"contrast": ,
    "neut": [Int, Positive],
    "bcoh": [Complex],
    "binc": [Complex],
    "scoh_bound": [Complex] ,
    "sinc_bound":[Complex]
}
units = {
    "symbol":'-'  ,
    "iso": '-',
    "name":'-' ,
    "A": 'uma',
    "Z": '-',
    "weight": 'uma',
    "excess":'-' ,
    "spin":'-',
    "parity":'-' ,
    "abundance": '%',
    "life":'-' ,
    "re_bcoh": 'fm',
    "im_bcoh":'fm' ,
    "bplus": 'fm',
    "bminus": 'fm',
    "re_binc": 'fm',
    "im_binc": 'fm',
    "sig_coh": 'barns',
    "sig_inc": 'barns',
    "sig_sca": 'barns',
    "sig_abs": 'barns',
    "contrast":'-' ,
    "neut":'-' ,
    "bcoh":'fm' ,
    "binc":'fm',
    "scoh_bound":'barns',
    "sinc_bound":'barns'
}

nombres = ["symbol", "iso", "name", "A", "Z", "weight", "excess", "spin", "parity", "abundance", "life", "re_bcoh","im_bcoh","bplus", "bminus", "re_binc", "im_binc", "sig_coh", "sig_inc", "sig_sca", "sig_abs", "contrast", "neut", "bcoh", "binc", "scoh_bound", "sinc_bound"]
isotope_dict = {
'n1': ('n', '1', 'neutron', 1, 0, 1.00866, 8.07132, '1/2', '+', 'NULL', '10.24m', -37.8, 0, 'NULL', 'NULL', 'NULL', 'NULL', 44.89, 0, 44.89, 0, 'NULL'),
'H1': ('H', '1', 'proton', 1, 1, 1.00728, 6.77799, '1/2', '+', 99.985, 'inf', -37.7423, 0, 10.817, -47.42, 25.217, 0, 1.7589, 79.91, 81.6689, 0.3326, 100.893),
'H2': ('H', '2', 'deuteron', 2, 1, 2.01355, 12.6247, '1', '+', 0.015, 'inf', 6.674, 0, 9.53, 0.975, 4.03, 0, 5.597, 2.04, 7.637, 0.000519, 2.18612),
'H3': ('H', '3', 'triton', 3, 1, 3.0155, 14.4388, '1/2', '+', 0, '12.312y', 4.792, 0, 4.18, 6.56, -1.04, 0, 2.89, 0.14, 3.03, 0.000006, 0.642565),
'H': ('H', 'nat', 'Hydrogen', 1, 1, 1.00794, 0, 'NULL', 'NULL', 100, 'inf', -3.739, 0, 'NULL', 'NULL', 'NULL', 'NULL', 1.7568, 80.26, 82.0168, 0.3326, 0),
'He3': ('He', '3', 'Helium', 3, 2, 'NULL', 14.9312, '1/2', '+', 0.000137, 'inf', 5.74, 1.483, 4.5, 9.3, -2.1, 2.568, 4.42, 1.38, 5.8, 5333, 2.30713),
'He4': ('He', '4', 'Helium', 4, 2, 'NULL', 2.4249, '0', '+', 99.9999, 'inf', 3.26, 0, 'NULL', 'NULL', 0, 0, 1.34, 0, 1.34, 0, 0),
'He': ('He', 'nat', 'Helium', 4, 2, 4.0026, 0, 'NULL', 'NULL', 100, 'inf', 3.26, 0, 'NULL', 'NULL', 'NULL', 'NULL', 1.34, 0, 1.34, 0.00747, 0),
'Li6': ('Li', '6', 'Lithium', 6, 3, 'NULL', 14.0868, '1', '+', 7.59, 'inf', 2, 0.261, 0.67, 4.67, -1.89, 0.26, 0.51, 0.46, 0.97, 940, 0.126903),
'Li7': ('Li', '7', 'Lithium', 7, 3, 'NULL', 14.9081, '3/2', '-', 92.41, 'inf', -2.22, 0, -4.15, 1, -2.49, 0, 0.619, 0.78, 1.399, 0.0454, 0.365208),
'Li': ('Li', 'nat', 'Lithium', 7, 3, 6.941, 0, 'NULL', 'NULL', 100, 'inf', -1.9, 0, 'NULL', 'NULL', 'NULL', 'NULL', 0.454, 0.92, 1.374, 70.5, 0),
'Be': ('Be', 'nat', 'Beryllium', 9, 4, 9.01218, 11.3476, '3/2', '-', 100, 'inf', 7.79, 0, 'NULL', 'NULL', 0.12, 0, 7.63, 0.0018, 7.6318, 0.0076, 0),
'B10': ('B', '10', 'Boron', 10, 5, 'NULL', 12.0507, '3', '+', 19.8, 'inf', 0.1, -1.066, -4.2, 5.2, -4.7, 1.231, 0.144, 3, 3.144, 3835, 0.959256),
'B11': ('B', '11', 'Boron', 11, 5, 'NULL', 8.6679, '3/2', '-', 80.2, 'inf', 6.65, 0, 5.6, 8.3, -1.3, 0, 5.56, 0.21, 5.77, 20.0055, 0.571776),
'B': ('B', 'nat', 'Boron', 11, 5, 10.811, 0, 'NULL', 'NULL', 100, 'inf', 5.3, -0.213, 'NULL', 'NULL', 'NULL', 'NULL', 3.54, 1.7, 5.24, 767, 0),
'C12': ('C', '12', 'Carbon', 12, 6, 'NULL', 0, '0', '+', 98.89, 'inf', 6.6535, 0, 'NULL', 'NULL', 0, 0, 5.559, 0, 5.559, 0.00353, 0.00153479),
'C13': ('C', '13', 'Carbon', 13, 6, 'NULL', 3.125, '1/2', '-', 1.11, 'inf', 6.19, 0, 5.6, 6.2, -0.25, 0, 4.81, 0.034, 4.844, 0.00137, 0.133144),
'C': ('C', 'nat', 'Carbon', 12, 6, 12.0107, 0, 'NULL', 'NULL', 100, 'inf', 6.6484, 0, 'NULL', 'NULL', 'NULL', 'NULL', 5.551, 0.001, 5.552, 0.0035, 0),
'N14': ('N', '14', 'Nitrogen', 14, 7, 'NULL', 2.8634, '1', '+', 99.634, 'inf', 9.37, 0, 10.7, 6.2, 2.1, 0, 11.03, 0.5, 11.53, 1.91, 0.00213789),
'N15': ('N', '15', 'Nitrogen', 15, 7, 'NULL', 0.1014, '1/2', '-', 0.366, 'inf', 6.44, 0, 6.77, 6.21, 0.24, 0, 5.21, 0.00005, 5.21005, 0.000024, 0.526609),
'N': ('N', 'nat', 'Nitrogen', 14, 7, 14.0067, 0, 'NULL', 'NULL', 100, 'inf', 9.36, 0, 'NULL', 'NULL', 'NULL', 'NULL', 11.01, 0.5, 11.51, 1.9, 0),
'O16': ('O', '16', 'Oxygen', 16, 8, 'NULL', -4.737, '0', '+', 99.762, 'inf', 5.805, 0, 'NULL', 'NULL', 0, 0, 4.232, 0, 4.232, 0.0001, 0),
'O17': ('O', '17', 'Oxygen', 17, 8, 'NULL', -0.8088, '5/2', '+', 0.038, 'inf', 5.66, 0, 5.86, 5.41, 0.17, 0, 4.2, 0.004, 4.204, 0.236, 0.049333),
'O18': ('O', '18', 'Oxygen', 18, 8, 'NULL', -0.7815, '0', '+', 0.2, 'inf', 5.84, 0, 'NULL', 'NULL', 0, 0, 4.29, 0, 4.29, 0.00016, 0.0120949),
'O': ('O', 'nat', 'Oxygen', 16, 8, 15.9994, 0, 'NULL', 'NULL', 100, 'inf', 5.805, 0, 'NULL', 'NULL', 'NULL', 'NULL', 4.232, 0, 4.232, 0.00019, 0),
'F': ('F', 'nat', 'Fluorine', 19, 9, 18.9984, -1.4874, '1/2', '+', 100, 'inf', 5.654, 0, 5.632, 5.767, -0.082, 0, 4.017, 0.0008, 4.0178, 0.0096, 0),
'Ne20': ('Ne', '20', 'Neon', 20, 10, 'NULL', -7.0419, '0', '+', 90.48, 'inf', 4.631, 0, 'NULL', 'NULL', 0, 0, 2.695, 0, 2.695, 0.036, 0.028674),
'Ne21': ('Ne', '21', 'Neon', 21, 10, 'NULL', -5.7318, '3/2', '+', 0.27, 'inf', 6.66, 0, 'NULL', 'NULL', 0.6, 0, 5.6, 0.05, 5.65, 0.67, 1.12753),
'Ne22': ('Ne', '22', 'Neon', 22, 10, 'NULL', -8.0247, '0', '+', 9.25, 'inf', 3.87, 0, 'NULL', 'NULL', 0, 0, 1.88, 0, 1.88, 0.046, 0.281627),
'Ne': ('Ne', 'nat', 'Neon', 20, 10, 20.1797, 0, 'NULL', 'NULL', 100, 'inf', 4.566, 0, 'NULL', 'NULL', 'NULL', 'NULL', 2.62, 0.008, 2.628, 0.039, 0),
'Na': ('Na', 'nat', 'Sodium', 23, 11, 22.9898, -9.5299, '3/2', '+', 100, 'inf', 3.63, 0, 6.42, -1, 3.59, 0, 1.66, 1.62, 3.28, 0.53, 0),
'Mg24': ('Mg', '24', 'Magnesium', 24, 12, 'NULL', -13.9336, '0', '+', 78.99, 'inf', 5.49, 0, 'NULL', 'NULL', 0, 0, 4.03, 0, 4.03, 0.05, 0.0432485),
'Mg25': ('Mg', '25', 'Magnesium', 25, 12, 'NULL', -13.1928, '5/2', '+', 10, 'inf', 3.62, 0, 4.73, 1.76, 1.48, 0, 1.65, 0.28, 1.93, 0.19, 0.546413),
'Mg26': ('Mg', '26', 'Magnesium', 26, 12, 'NULL', -16.2146, '0', '+', 11.01, 'inf', 4.89, 0, 'NULL', 'NULL', 0, 0, 3, 0, 3, 0.0382, 0.172323),
'Mg': ('Mg', 'nat', 'Magnesium', 24, 12, 24.305, 0, 'NULL', 'NULL', 100, 'inf', 5.375, 0, 'NULL', 'NULL', 'NULL', 'NULL', 3.631, 0.08, 3.711, 0.063, 0),
'Al': ('Al', 'nat', 'Aluminium', 27, 13, 26.9815, -17.1967, '5/2', '+', 100, 'inf', 3.449, 0, 3.67, 3.15, 0.256, 0, 1.495, 0.0082, 1.5032, 0.231, 0),
'Si28': ('Si', '28', 'Silicon', 28, 14, 27.9769, -21.4928, '0', '+', 92.22, 'inf', 4.106, 0, 'NULL', 'NULL', 0, 0, 2.12, 0, 2.12, 0.177, 0.0214273),
'Si29': ('Si', '29', 'Silicon', 29, 14, 28.9765, -21.895, '1/2', '+', 4.69, 'inf', 4.7, 0, 4.5, 4.7, -1.08, 0, 2.78, 0.79, 3.57, 0.101, 0.282186),
'Si30': ('Si', '30', 'Silicon', 30, 14, 29.9738, -24.4329, '0', '+', 3.09, 'inf', 4.58, 0, 'NULL', 'NULL', 0, 0, 2.64, 0, 2.64, 0.107, 0.217548),
'Si': ('Si', 'nat', 'Silicon', 28, 14, 28.085, 'NULL', 'NULL', 'NULL', 100, 'inf', 4.15071, 0, 'NULL', 'NULL', 'NULL', 'NULL', 2.1633, 0.004, 2.167, 0.171, 0),
'P': ('P', 'nat', 'Phosphorus', 31, 15, 30.974, -24.4409, '1/2', '+', 100, 'inf', 5.13, 0, 'NULL', 'NULL', 0.3, 0, 3.307, 0.005, 3.312, 0.172, 0),
'S32': ('S', '32', 'Sulfur', 32, 16, 31.9721, -26.0157, '0', '+', 95.02, 'inf', 2.804, 0, 'NULL', 'NULL', 0, 0, 0.988, 0, 0.988, 0.54, 0.0299791),
'S33': ('S', '33', 'Sulfur', 33, 16, 32.9715, -26.586, '3/2', '+', 0.75, 'inf', 4.74, 0, 'NULL', 'NULL', 1.5, 0, 2.8, 0.3, 3.1, 0.54, 1.77193),
'S34': ('S', '34', 'Sulfur', 34, 16, 33.9679, -29.9318, '0', '+', 4.21, 'inf', 3.48, 0, 'NULL', 'NULL', 0, 0, 1.52, 0, 1.52, 0.227, 0.494113),
'S36': ('S', '36', 'Sulfur', 36, 16, 35.9671, -28.8464, '0', '+', 0.02, 'inf', 3, 0, 'NULL', 'NULL', 0, 0, 1.1, 0, 1.1, 0.15, 0.11037),
'S': ('S', 'nat', 'Sulfur', 32, 16, 32.065, 0, 'NULL', 'NULL', 100, 'inf', 2.847, 0, 'NULL', 'NULL', 'NULL', 'NULL', 1.0186, 0.007, 1.0256, 0.53, 0),
'Cl35': ('Cl', '35', 'Chlorine', 35, 17, 34.9689, -29.0135, '3/2', '+', 75.78, 'inf', 11.709, 0, 16.3, 4, 0, 0, 17.06, 4.7, 21.8, 44.1, 0.494105),
'Cl37': ('Cl', '37', 'Chlorine', 37, 17, 36.9659, -31.7615, '3/2', '+', 24.22, 'inf', 3.08, 0, 3.1, 3.05, 0.02, 0, 1.19, 0.001, 1.19, 0.433, 0.896618),
'Cl': ('Cl', 'nat', 'Chlorine', 35, 17, 35.453, 0, 'NULL', 'NULL', 100, 'inf', 9.5792, 0, 'NULL', 'NULL', 'NULL', 'NULL', 11.528, 5.3, 16.828, 33.5, 0),
'Ar36': ('Ar', '36', 'Argon', 36, 18, 35.6975, -30.2315, '0', '+', 0.337, 'inf', 24.9, 0, 'NULL', 'NULL', 0, 0, 77.9, 0, 77.9, 5.2, 169.132),
'Ar38': ('Ar', '38', 'Argon', 38, 18, 37.9627, -34.7146, '0', '+', 0.063, 'inf', 3.5, 0, 'NULL', 'NULL', 0, 0, 1.5, 0, 1.5, 0.8, 2.36143),
'Ar40': ('Ar', '40', 'Argon', 40, 18, 39.9624, -35.0399, '0', '+', 99.6, 'inf', 1.84, 0, 'NULL', 'NULL', 0, 0, 0.421, 0, 0.421, 0.66, 0.0709827),
'Ar': ('Ar', 'nat', 'Argon', 40, 18, 39.948, 'NULL', '0', '+', 100, 'inf', 1.909, 0, 'NULL', 'NULL', 'NULL', 'NULL', 0.458, 0.225, 0.683, 0.675, 0),
'K39': ('K', '39', 'Potassium', 39, 19, 38.9637, -33.807, '3/2', '+', 93.2581, 'inf', 3.72, 0, 5.15, 5.15, 1.43, 0, 1.76, 0.25, 2.01, 2.1, 0.0274336),
'K40': ('K', '40', 'Potassium', 40, 19, 39.964, -33.5352, '4', '-', 0.0117, 'inf', 3.1, 0, 'NULL', 'NULL', 'NULL', 'NULL', 1.1, 0.5, 1.6, 35, 0.286504),
'K41': ('K', '41', 'Potassium', 41, 19, 40.9618, -35.5591, '3/2', '+', 6.7302, 'inf', 2.69, 0, 'NULL', 'NULL', 1.51, 0, 0.91, 0.3, 1.2, 1.46, 0.462755),
'K': ('K', 'nat', 'Potassium', 39, 19, 39.0983, 'NULL', 'NULL', 'NULL', 100, 'inf', 3.67, 0, 'NULL', 'NULL', 'NULL', 'NULL', 1.69, 0.27, 1.96, 2.1, 0),
'Ca40': ('Ca', '40', 'Calcium', 40, 20, 39.9626, -34.8463, '0', '+', 96.941, 'inf', 4.78, 0, 'NULL', 'NULL', 0, 0, 2.9, 0, 2.9, 0.41, 0.0343323),
'Ca42': ('Ca', '42', 'Calcium', 42, 20, 41.9586, -38.5471, '0', '+', 0.647, 'inf', 3.36, 0, 'NULL', 'NULL', 0, 0, 1.42, 0, 1.42, 0.68, 0.488927),
'Ca43': ('Ca', '43', 'Calcium', 43, 20, 42.9588, -38.4086, '7/2', '-', 0.135, 'inf', -1.56, 0, 'NULL', 'NULL', 'NULL', 'NULL', 0.31, 0.5, 0.8, 6.2, 0.889832),
'Ca44': ('Ca', '44', 'Calcium', 44, 20, 43.9555, -41.4685, '0', '+', 2.086, 'inf', 1.42, 0, 'NULL', 'NULL', 0, 0, 0.25, 0, 0.25, 0.88, 0.908719),
'Ca46': ('Ca', '46', 'Calcium', 46, 20, 45.9537, -43.1351, '0', '+', 0.004, 'inf', 3.55, 0, 'NULL', 'NULL', 0, 0, 1.6, 0, 1.6, 0.7, 0.429493),
'Ca48': ('Ca', '48', 'Calcium', 48, 20, 47.9525, -44.2141, '0', '+', 0.187, 'inf', 0.39, 0, 'NULL', 'NULL', 0, 0, 0.019, 0, 0.019, 1.09, 0.993115),
'Ca': ('Ca', 'nat', 'Calcium', 40, 20, 40.078, 'NULL', 'NULL', 'NULL', 100, 'inf', 4.7, 0, 'NULL', 'NULL', 'NULL', 'NULL', 2.78, 0.05, 2.83, 0.43, 0),
'Sc': ('Sc', 'nat', 'Scandium', 35, 21, 44.9559, -41.0678, '7/2', '-', 100, 'inf', 12.1, 0, 6.91, 18.99, -6.02, 0, 19.03, 4.5, 23.5, 27.5, 0),
'Ti46': ('Ti', '46', 'Titanium', 46, 22, 'NULL', -44.1234, '0', '+', 8.25, 'inf', 4.72, 0, 'NULL', 'NULL', 0, 0, 3.05, 0, 3.05, 0.59, 0.961662),
'Ti47': ('Ti', '47', 'Titanium', 47, 22, 'NULL', -44.9324, '5/2', '-', 7.44, 'inf', 3.53, 0, 0.46, 7.64, -3.5, 0, 1.66, 1.5, 3.16, 1.7, 0.0972096),
'Ti48': ('Ti', '48', 'Titanium', 48, 22, 'NULL', -48.4877, '0', '+', 73.72, 'inf', -5.86, 0, 'NULL', 'NULL', 0, 0, 4.65, 0, 4.65, 7.84, 2.02368),
'Ti49': ('Ti', '49', 'Titanium', 49, 22, 'NULL', -48.5588, '7/2', '-', 5.41, 'inf', 0.98, 0, 2.6, -1.2, 1.9, 0, 0.14, 3.3, 3.44, 2.2, 0.915435),
'Ti50': ('Ti', '50', 'Titanium', 50, 22, 'NULL', -51.4267, '0', '+', 5.18, 'inf', 5.88, 0, 'NULL', 'NULL', 0, 0, 4.8, 0, 4.8, 0.179, 2.04435),
'Ti': ('Ti', 'nat', 'Titanium', 48, 22, 47.867, 0, 'NULL', 'NULL', 100, 'inf', -3.37, 0, 'NULL', 'NULL', 'NULL', 'NULL', 1.485, 2.87, 4.355, 6.09, 0),
'V50': ('V', '50', 'Vanadium', 50, 23, 'NULL', -49.2216, '6', '+', 0.25, 'inf', 7.6, 0, 'NULL', 'NULL', 'NULL', 'NULL', 7.3, 0.5, 7.8, 60, 293.32),
'V51': ('V', '51', 'Vanadium', 51, 23, 'NULL', -52.2014, '7/2', '-', 99.75, 'inf', -0.402, 0, 4.93, -7.58, 6.35, 0, 0.0203, 5.07, 5.0903, 4.9, 0.176536),
'V': ('V', 'nat', 'Vanadium', 51, 23, 50.9415, 0, 'NULL', 'NULL', 100, 'inf', -0.443, 0, 'NULL', 'NULL', 'NULL', 'NULL', 0.01838, 5.08, 5.09838, 5.08, 0),
'Cr50': ('Cr', '50', 'Chromium', 50, 24, 49.946, -50.2595, '0', '+', 4.345, 'inf', -4.5, 0, 'NULL', 'NULL', 0, 0, 2.54, 0, 2.54, 15.8, 0.532555),
'Cr52': ('Cr', '52', 'Chromium', 52, 24, 51.9405, -55.4169, '0', '+', 83.789, 'inf', 4.914, 0, 'NULL', 'NULL', 0, 0, 3.042, 0, 3.042, 0.76, 0.827517),
'Cr53': ('Cr', '53', 'Chromium', 53, 24, 52.9407, -55.2847, '3/2', '-', 9.501, 'inf', -4.2, 0, 1.16, -13, 6.9, 0, 2.22, 5.93, 8.15, 18.1, 0.335026),
'Cr54': ('Cr', '54', 'Chromium', 54, 24, 53.9389, -56.9325, '0', '+', 2.365, 'inf', 4.55, 0, 'NULL', 'NULL', 0, 0, 2.6, 0, 2.6, 0.36, 0.566801),
'Cr': ('Cr', 'nat', 'Chromium', 52, 24, 51.9961, 'NULL', 'NULL', 'NULL', 100, 'inf', 3.635, 0, 'NULL', 'NULL', 'NULL', 'NULL', 1.66, 1.83, 3.49, 3.05, 0),
'Mn': ('Mn', 'nat', 'Manganese', 55, 25, 54.938, -57.7106, '5/2', '-', 100, 'inf', -3.75, 0, -4.93, -1.46, -1.71, 0, 1.75, 0.4, 2.15, 13.3, 0),
'Fe54': ('Fe', '54', 'Iron', 54, 26, 53.9396, -56.2525, '0', '+', 5.845, 'inf', 4.2, 0, 'NULL', 'NULL', 0, 0, 2.2, 0, 2.2, 2.25, 0.802469),
'Fe56': ('Fe', '56', 'Iron', 56, 26, 55.9349, -60.6054, '0', '+', 91.754, 'inf', 10.1, 0, 'NULL', 'NULL', 0, 0, 12.42, 0, 12.42, 2.59, 0.142297),
'Fe57': ('Fe', '57', 'Iron', 57, 26, 56.9354, -60.1801, '1/2', '-', 2.119, 'inf', 2.3, 0, 'NULL', 'NULL', 'NULL', 'NULL', 0.66, 0.3, 1.03, 2.48, 0.940763),
'Fe58': ('Fe', '58', 'Iron', 58, 26, 57.9333, -62.1534, '0', '+', 0.282, 'inf', 15, 0, 'NULL', 'NULL', 0, 0, 28, 0, 28, 1.28, 1.51953),
'Fe': ('Fe', 'nat', 'Iron', 56, 26, 55.845, 'NULL', 'NULL', 'NULL', 100, 'inf', 9.45, 0, 'NULL', 'NULL', 'NULL', 'NULL', 11.22, 0.4, 11.62, 2.56, 0),
'Co': ('Co', 'nat', 'Cobalt', 59, 27, 58.9332, -62.2284, 'NULL', 'NULL', 100, 'inf', 2.49, 0, -9.21, 3.58, -6.43, 0, 0.779, 4.8, 5.6, 37.18, 0),
'Ni58': ('Ni', '58', 'Nickel', 58, 28, 57.9353, -60.2277, '0', '+', 68.0769, 'inf', 68.077, 0, 'NULL', 'NULL', 0, 0, 26.1, 0, 26.1, 4.6, 42.6844),
'Ni60': ('Ni', '60', 'Nickel', 60, 28, 59.9308, -64.4721, '0', '+', 26.2231, 'inf', 2.8, 0, 'NULL', 'NULL', 0, 0, 0.99, 0, 0.99, 2.9, 0.9261),
'Ni61': ('Ni', '61', 'Nickel', 61, 28, 60.9311, -64.2209, '3/2', '+', 1.1399, 'inf', 7.6, 0, 'NULL', 'NULL', -3.9, 0, 7.26, 1.9, 5.6, 2.5, 0.455557),
'Ni62': ('Ni', '62', 'Nickel', 62, 28, 61.9283, -66.7461, '0', '+', 3.6345, 'inf', -8.7, 0, 'NULL', 'NULL', 0, 0, 9.5, 0, 9.5, 14.5, 0.286549),
'Ni64': ('Ni', '64', 'Nickel', 64, 28, 63.928, -67.0993, '0', '+', 0.9256, 'inf', -0.37, 0, 'NULL', 'NULL', 0, 0, 0.017, 0, 0.017, 1.52, 0.99871),
'Ni': ('Ni', 'nat', 'Nickel', 59, 28, 58.6934, 'NULL', 'NULL', 'NULL', 100, 'inf', 10.3, 0, 'NULL', 'NULL', 'NULL', 'NULL', 13.3, 5.2, 18.5, 4.49, 0),
'Cu63': ('Cu', '63', 'Copper', 63, 29, 'NULL', -65.5795, '3/2', '-', 69.15, 'inf', 6.477, 0, 'NULL', 'NULL', 0.22, 0, 5.2, 0.006, 5.206, 4.5, 0.295732),
'Cu65': ('Cu', '65', 'Copper', 65, 29, 'NULL', -67.2637, '3/2', '-', 30.85, 'inf', 10.204, 0, 'NULL', 'NULL', 1.82, 0, 14.1, 0.4, 14.5, 2.17, 0.747959),
'Cu': ('Cu', 'nat', 'Copper', 64, 29, 63.546, 0, 'NULL', 'NULL', 100, 'inf', 7.718, 0, 'NULL', 'NULL', 'NULL', 'NULL', 7.485, 0.55, 8.035, 3.78, 0),
'Zn64': ('Zn', '64', 'Zinc', 64, 30, 63.9291, -66.0036, '0', '+', 48.63, 'inf', 5.23, 0, 'NULL', 'NULL', 0, 0, 3.42, 0, 3.42, 0.93, 0.152174),
'Zn66': ('Zn', '66', 'Zinc', 66, 30, 65.926, -68.8994, '0', '+', 27.9, 'inf', 5.98, 0, 'NULL', 'NULL', 0, 0, 4.48, 0, 4.48, 0.62, 0.108423),
'Zn67': ('Zn', '67', 'Zinc', 67, 30, 66.9271, -67.8804, '5/2', '-', 4.1, 'inf', 7.58, 0, 5.8, 10.1, -1.5, 0, 7.18, 0.28, 7.46, 6.8, 0.780909),
'Zn68': ('Zn', '68', 'Zinc', 68, 30, 67.9248, -70.0072, '0', '+', 18.75, 'inf', 6.04, 0, 'NULL', 'NULL', 0, 0, 4.57, 0, 4.57, 1.1, 0.130778),
'Zn70': ('Zn', '70', 'Zinc', 70, 30, 69.9253, -69.5647, '0', '+', 0.62, 'inf', 6, 0, 'NULL', 'NULL', 0, 0, 4.5, 0, 4.5, 0.092, 0.11585),
'Zn': ('Zn', 'nat', 'Zinc', 65, 30, 65.38, 0, 'NULL', 'NULL', 100, 'inf', 5.68, 0, 'NULL', 'NULL', 'NULL', 'NULL', 4.054, 0.077, 4.131, 1.11, 0),
'Ga69': ('Ga', '69', 'Gallium', 69, 31, 68.9256, -69.3278, '3/2', '-', 60.108, 'inf', 8.0403, 0, 6.3, 10.5, -0.85, 0, 7.8, 0.091, 7.89, 2.18, 0.217104),
'Ga71': ('Ga', '71', 'Gallium', 71, 31, 70.9247, -70.1402, '3/2', '-', 39.892, 'inf', 6.17, 0, 5.5, 7.8, -0.82, 0, 5.15, 0.084, 5.23, 3.61, 0.283273),
'Ga': ('Ga', 'nat', 'Gallium', 70, 31, 69.723, 0, 'NULL', 'NULL', 100, 'inf', 7.288, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.675, 0.16, 6.835, 2.75, 0),
'Ge70': ('Ge', '70', 'Germanium', 70, 32, 69.9242, -70.5631, '0', '+', 20.84, 'inf', 10, 0, 'NULL', 'NULL', 0, 0, 12.6, 0, 12.6, 3, 0.492666),
'Ge72': ('Ge', '72', 'Germanium', 72, 32, 71.9221, -72.5859, '0', '+', 27.54, 'inf', 8.51, 0, 'NULL', 'NULL', 0, 0, 9.1, 0, 9.1, 0.8, 0.0809902),
'Ge73': ('Ge', '73', 'Germanium', 73, 32, 72.9235, -71.2975, '9/2', '+', 7.73, 'inf', 5.02, 0, 5.5, 7.8, 3.43, 0, 3.17, 1.5, 4.7, 15.1, 0.623842),
'Ge74': ('Ge', '74', 'Germanium', 74, 32, 73.9212, -73.4224, '0', '+', 36.28, 'inf', 7.58, 0, 'NULL', 'NULL', 0, 0, 7.2, 0, 7.2, 0.4, 0.142368),
'Ge76': ('Ge', '76', 'Germanium', 76, 32, 75.9214, -73.213, '0', '+', 7.61, 'inf', 8.2, 0, 'NULL', 'NULL', 0, 0, 8, 0, 8, 0.16, 0.0036686),
'Ge': ('Ge', 'nat', 'Germanium', 73, 32, 72.64, 0, 'NULL', 'NULL', 100, 'inf', 8.185, 0, 'NULL', 'NULL', 'NULL', 'NULL', 8.42, 0.18, 8.6, 2.2, 0),
'As': ('As', 'nat', 'Arsenic', 75, 33, 74.9216, -73.0324, '3/2', '-', 100, 'inf', 6.58, 0, 6.04, 7.47, -0.69, 0, 5.44, 0.06, 5.5, 4.5, 0),
'Se74': ('Se', '74', 'Selenium', 74, 34, 73.9225, -72.2127, '0', '+', 0.89, 'inf', 0.8, 0, 'NULL', 'NULL', 0, 0, 0.1, 0, 0.1, 51.8, 0.989925),
'Se76': ('Se', '76', 'Selenium', 76, 34, 75.9192, -75.252, '0', '+', 9.37, 'inf', 12.2, 0, 'NULL', 'NULL', 0, 0, 18.7, 0, 18.7, 85, 1.34317),
'Se77': ('Se', '77', 'Selenium', 77, 34, 76.9199, -74.5996, '1/2', '-', 7.63, 'inf', 8.25, 0, 'NULL', 'NULL', -0.6, 0, 8.6, 0.05, 8.65, 42, 0.0714977),
'Se78': ('Se', '78', 'Selenium', 78, 34, 77.9173, -77.0261, '0', '+', 23.77, 'inf', 8.24, 0, 'NULL', 'NULL', 0, 0, 8.5, 0, 8.5, 0.43, 0.0689017),
'Se80': ('Se', '80', 'Selenium', 80, 34, 79.9165, -77.7599, '0', '+', 49.61, 'inf', 7.48, 0, 'NULL', 'NULL', 0, 0, 7.03, 0, 7.03, 0.61, 0.119181),
'Se82': ('Se', '82', 'Selenium', 82, 34, 81.9167, -77.594, '0', '+', 8.73, 'inf', 6.43, 0, 'NULL', 'NULL', 0, 0, 5.05, 0, 5.05, 0.044, 0.349113),
'Se': ('Se', 'nat', 'Selenium', 79, 34, 78.96, 0, 'NULL', 'NULL', 100, 'inf', 7.97, 0, 'NULL', 'NULL', 'NULL', 'NULL', 7.98, 0.32, 8.3, 11.7, 0),
'Br79': ('Br', '79', 'Bromine', 79, 35, 78.9183, -76.0685, '3/2', '-', 50.69, 'inf', 6.79, 0, 'NULL', 'NULL', -1.1, 0, 5.81, 0.15, 5.96, 11, 0.240261),
'Br81': ('Br', '81', 'Bromine', 81, 35, 80.9163, -77.9748, '3/2', '-', 49.31, 'inf', 6.78, 0, 'NULL', 'NULL', -0.6, 0, 5.79, 0.05, 5.84, 2.7, 0.242497),
'Br': ('Br', 'nat', 'Bromine', 80, 35, 79.904, 'NULL', 'NULL', 'NULL', 100, 'inf', 7.79, 0, 'NULL', 'NULL', 'NULL', 'NULL', 5.8, 0.1, 5.9, 6.9, 0),
'Kr78': ('Kr', '78', 'Krypton', 78, 36, 77.9204, -74.1797, '0', '+', 0.355, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 6.4, 1),
'Kr80': ('Kr', '80', 'Krypton', 80, 36, 79.9164, -77.8925, '0', '+', 2.2286, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 11.8, 1),
'Kr82': ('Kr', '82', 'Krypton', 82, 36, 81.9135, -80.5895, '0', '+', 11.593, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 29, 1),
'Kr83': ('Kr', '83', 'Krypton', 83, 36, 82.9141, -79.9817, '9/2', '+', 11.5, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0, 'NULL', 185, 1),
'Kr84': ('Kr', '84', 'Krypton', 84, 36, 83.9115, -82.431, '0', '+', 56.987, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 6.6, 0.113, 1),
'Kr86': ('Kr', '86', 'Krypton', 86, 36, 85.9106, -83.2656, '0', '+', 17.279, 'inf', 8.07, 0, 'NULL', 'NULL', 0, 0, 8.2, 0, 8.2, 0.003, 0.0676896),
'Kr': ('Kr', 'nat', 'Krypton', 84, 36, 83.798, 'NULL', 'NULL', 'NULL', 100, 'inf', 7.81, 0, 'NULL', 'NULL', 'NULL', 'NULL', 7.67, 0.01, 7.68, 25, 0),
'Rb85': ('Rb', '85', 'Rubidium', 85, 37, 84.9118, -82.1673, '5/2', '-', 72.17, 'inf', 7.07, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.2, 0.5, 6.7, 0.48, 0.00282286),
'Rb87': ('Rb', '87', 'Rubidium', 87, 37, 86.9092, -84.5978, '3/2', '-', 27.83, 'inf', 7.27, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.6, 0.5, 7.1, 0.12, 0.0543925),
'Rb': ('Rb', 'nat', 'Rubidium', 85, 37, 85.4678, 'NULL', 'NULL', 'NULL', 100, 'inf', 7.08, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.32, 0.5, 6.8, 0.38, 0),
'Sr84': ('Sr', '84', 'Strontium', 84, 38, 83.9134, -80.6438, '0', '+', 0.56, 'inf', 5, 0, 'NULL', 'NULL', 0, 0, 6, 0, 6, 0.87, 0.492699),
'Sr86': ('Sr', '86', 'Strontium', 86, 38, 85.9093, -84.5236, '0', '+', 9.86, 'inf', 5.68, 0, 'NULL', 'NULL', 0, 0, 4.04, 0, 4.04, 1.04, 0.34533),
'Sr87': ('Sr', '87', 'Strontium', 87, 38, 86.9089, -84.8804, '9/2', '+', 7, 'inf', 7.41, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.88, 0.5, 7.4, 16, 0.114198),
'Sr88': ('Sr', '88', 'Strontium', 88, 38, 87.9056, -87.9217, '0', '+', 82.58, 'inf', 7.16, 0, 'NULL', 'NULL', 0, 0, 6.42, 0, 6.42, 0.058, 0.0402838),
'Sr': ('Sr', 'nat', 'Strontium', 87, 38, 87.62, 'NULL', 'NULL', 'NULL', 100, 'inf', 7.02, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.19, 0.06, 6.25, 1.28, 0),
'Y': ('Y', 'nat', 'Yttrium', 89, 39, 88.9059, -87.7018, '1/2', '-', 100, 'inf', 7.75, 0, 8.4, 5.8, 1.1, 0, 7.55, 0.15, 7.7, 1.28, 0),
'Zr90': ('Zr', '90', 'Zirconium', 90, 40, 'NULL', -88.7673, '0', '+', 51.45, 'inf', 6.5, 0, 'NULL', 'NULL', 0, 0, 5.1, 0, 5.1, 0.011, 0.175861),
'Zr91': ('Zr', '91', 'Zirconium', 91, 40, 'NULL', -87.8904, '5/2', '+', 11.22, 'inf', 8.8, 0, 7.9, 10.1, -1.08, 0, 9.5, 0.15, 9.65, 1.17, 0.510565),
'Zr92': ('Zr', '92', 'Zirconium', 92, 40, 'NULL', -88.4539, '0', '+', 17.15, 'inf', 7.52, 0, 'NULL', 'NULL', 0, 0, 6.9, 0, 6.9, 0.22, 0.103087),
'Zr94': ('Zr', '94', 'Zirconium', 94, 40, 'NULL', -87.2668, '0', '+', 17.38, 'inf', 8.3, 0, 'NULL', 'NULL', 0, 0, 8.4, 0, 8.4, 0.0499, 0.343786),
'Zr96': ('Zr', '96', 'Zirconium', 96, 40, 'NULL', -85.4428, '0', '+', 2.8, 'inf', 5.5, 0, 'NULL', 'NULL', 0, 0, 3.8, 0, 3.8, 0.0229, 0.409936),
'Zr': ('Zr', 'nat', 'Zirconium', 91, 40, 91.224, 0, 'NULL', 'NULL', 100, 'inf', 7.16, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.44, 0.02, 6.46, 0.185, 0),
'Nb': ('Nb', 'nat', 'Niobium', 93, 41, 92.9064, -87.2083, '9/2', '+', 100, 'inf', 7.054, 0, 7.06, 7.35, -0.139, 0, 6.253, 0.0024, 6.2554, 1.15, 0),
'Mo92': ('Mo', '92', 'Molybdenum', 92, 42, 91.9068, -86.805, '0', '+', 14.84, 'inf', 6.93, 0, 'NULL', 'NULL', 0, 0, 6, 0, 6, 0.019, 0.0666488),
'Mo94': ('Mo', '94', 'Molybdenum', 94, 42, 93.9051, -88.4097, '0', '+', 9.25, 'inf', 6.82, 0, 'NULL', 'NULL', 0, 0, 5.81, 0, 5.81, 0.015, 0.0330556),
'Mo95': ('Mo', '95', 'Molybdenum', 95, 42, 94.9058, -87.7075, '5/2', '+', 15.92, 'inf', 6.93, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6, 0.5, 6.5, 13.1, 0.0666488),
'Mo96': ('Mo', '96', 'Molybdenum', 96, 42, 95.9047, -88.7905, '0', '+', 16.68, 'inf', 6.22, 0, 'NULL', 'NULL', 0, 0, 4.83, 0, 4.83, 0.5, 0.140718),
'Mo97': ('Mo', '97', 'Molybdenum', 97, 42, 96.906, -87.5404, '5/2', '+', 9.55, 'inf', 7.26, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.59, 0.5, 7.1, 2.5, 0.170653),
'Mo98': ('Mo', '98', 'Molybdenum', 98, 42, 97.9054, -88.1117, '0', '+', 24.13, 'inf', 6.6, 0, 'NULL', 'NULL', 0, 0, 5.44, 0, 5.44, 0.127, 0.0325181),
'Mo100': ('Mo', '100', 'Molybdenum', 100, 42, 99.9075, -86.1843, '0', '+', 9.63, 'inf', 6.75, 0, 'NULL', 'NULL', 0, 0, 5.69, 0, 5.69, 0.4, 0.011958),
'Mo': ('Mo', 'nat', 'Molybdenum', 96, 42, 95.96, 'NULL', 'NULL', 'NULL', 100, 'inf', 6.71, 0, 'NULL', 'NULL', 'NULL', 'NULL', 5.67, 0.04, 5.71, 2.48, 0),
'Tc': ('Tc', 'nat', 'Technecium', 99, 43, 98, -87.3231, '9/2', '+', 100, '2.111E+5y', 6.8, 0, 'NULL', 'NULL', 'NULL', 'NULL', 5.8, 0.5, 6.3, 20, 0),
'Ru96': ('Ru', '96', 'Ruthenium', 96, 44, 95.9076, -86.0721, '0', '+', 5.54, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 0.28, 1),
'Ru98': ('Ru', '98', 'Ruthenium', 98, 44, 97.9053, -88.2245, '0', '+', 1.87, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 0.8, 1),
'Ru99': ('Ru', '99', 'Ruthenium', 99, 44, 98.9059, -87.617, '5/2', '+', 12.76, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 6.9, 1),
'Ru100': ('Ru', '100', 'Ruthenium', 100, 44, 99.9042, -89.219, '0', '+', 12.6, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 4.8, 1),
'Ru101': ('Ru', '101', 'Ruthenium', 101, 44, 100.906, -87.9497, '5/2', '+', 17.06, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 3.3, 1),
'Ru102': ('Ru', '102', 'Ruthenium', 102, 44, 101.904, -89.098, '0', '+', 31.55, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 1.17, 1),
'Ru104': ('Ru', '104', 'Ruthenium', 104, 44, 103.905, -88.0889, '0', '+', 18.62, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 0.31, 1),
'Ru': ('Ru', 'nat', 'Ruthenium', 101, 44, 101.07, 'NULL', 'NULL', 'NULL', 100, 'inf', 7.02, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.21, 0.4, 6.6, 2.56, 0),
'Rh': ('Rh', 'nat', 'Rhodium', 103, 45, 102.906, -88.0222, '1/2', '-', 100, 'inf', 5.9, 0, 8.15, 6.74, 0.614, 0, 4.34, 0.047, 4.39, 144.8, 0),
'Pd102': ('Pd', '102', 'Palladium', 102, 46, 101.906, -87.9251, '0', '+', 1.02, 'inf', 7.7, 0, 'NULL', 'NULL', 0, 0, 7.5, 0, 7.5, 3.4, 0.697487),
'Pd104': ('Pd', '104', 'Palladium', 104, 46, 103.904, -89.39, '0', '+', 11.14, 'inf', 7.7, 0, 'NULL', 'NULL', 0, 0, 7.5, 0, 7.5, 0.6, 0.697487),
'Pd105': ('Pd', '105', 'Palladium', 105, 46, 104.905, -88.4128, '5/2', '+', 22.33, 'inf', 5.5, 0, 'NULL', 'NULL', -2.6, 0, 3.8, 0.8, 4.6, 20, 0.133935),
'Pd106': ('Pd', '106', 'Palladium', 106, 46, 105.903, -89.9025, '0', '+', 27.33, 'inf', 6.4, 0, 'NULL', 'NULL', 0, 0, 5.1, 0, 5.1, 0.304, 0.172695),
'Pd108': ('Pd', '108', 'Palladium', 108, 46, 107.904, -89.5243, '0', '+', 26.46, 'inf', 4.1, 0, 'NULL', 'NULL', 0, 0, 2.1, 0, 2.1, 8.5, 0.518726),
'Pd110': ('Pd', '110', 'Palladium', 110, 46, 109.905, -88.3492, '0', '+', 11.72, 'inf', 7.7, 0, 'NULL', 'NULL', 0, 0, 7.5, 0, 7.5, 0.226, 0.697487),
'Pd': ('Pd', 'nat', 'Palladium', 106, 46, 106.42, 'NULL', 'NULL', 'NULL', 100, 'inf', 5.91, 0, 'NULL', 'NULL', 'NULL', 'NULL', 4.39, 0.09, 4.48, 6.9, 0),
'Ag107': ('Ag', '107', 'Silver', 107, 47, 'NULL', -88.4017, '1/2', '-', 51.839, 'inf', 7.555, 0, 8.14, 5.8, 1, 0, 7.17, 0.13, 7.3, 37.6, 0.627542),
'Ag109': ('Ag', '109', 'Silver', 109, 47, 'NULL', -88.7227, '1/2', '-', 48.161, 'inf', 4.165, 0, 3.24, 6.9, -1.6, 0, 2.18, 0.32, 2.5, 91, 0.505355),
'Ag': ('Ag', 'nat', 'Silver', 108, 47, 107.868, 0, 'NULL', 'NULL', 100, 'inf', 5.922, 0, 'NULL', 'NULL', 'NULL', 'NULL', 4.407, 0.58, 4.987, 63.3, 0),
'Cd106': ('Cd', '106', 'Cadmium', 106, 48, 105.906, -87.1325, '0', '+', 1.25, 'inf', 5, 0, 'NULL', 'NULL', 0, 0, 3.1, 0, 3.1, 1, 0.0495867),
'Cd108': ('Cd', '108', 'Cadmium', 108, 48, 107.904, -89.2523, '0', '+', 0.89, 'inf', 5.31, 0, 'NULL', 'NULL', 0, 0, 3.7, 0, 3.7, 1.1, 0.18377),
'Cd110': ('Cd', '110', 'Cadmium', 110, 48, 109.903, -90.353, '0', '+', 12.49, 'inf', 5.78, 0, 'NULL', 'NULL', 0, 0, 4.4, 0, 4.4, 11, 0.4026),
'Cd111': ('Cd', '111', 'Cadmium', 111, 48, 110.904, -89.2575, '1/2', '+', 12.8, 'inf', 6.47, 0, 'NULL', 'NULL', 'NULL', 'NULL', 5.3, 0.3, 5.6, 24, 0.757466),
'Cd112': ('Cd', '112', 'Cadmium', 112, 48, 111.903, -90.5805, '0', '+', 24.13, 'inf', 6.34, 0, 'NULL', 'NULL', 0, 0, 5.1, 0, 5.1, 2.2, 0.687551),
'Cd113': ('Cd', '113', 'Cadmium', 113, 48, 112.904, -89.0493, '1/2', '+', 12.22, 'inf', -8, 5.73, 'NULL', 'NULL', 'NULL', 'NULL', 12.1, 0.3, 12.4, 20600, 3.06538),
'Cd114': ('Cd', '114', 'Cadmium', 114, 48, 113.903, -90.0209, '0', '+', 28.73, 'inf', 7.48, 0, 'NULL', 'NULL', 0, 0, 7.1, 0, 7.1, 0.34, 1.34899),
'Cd116': ('Cd', '116', 'Cadmium', 116, 48, 115.905, -88.7194, '0', '+', 7.49, 'inf', 6.26, 0, 'NULL', 'NULL', 0, 0, 5, 0, 5, 0.075, 0.645231),
'Cd': ('Cd', 'nat', 'Cadmium', 112, 48, 112.411, 'NULL', 'NULL', 'NULL', 100, 'inf', 4.83, -0.7, 'NULL', 'NULL', 'NULL', 'NULL', 3.04, 3.46, 6.5, 2520, 0),
'In113': ('In', '113', 'Indium', 113, 49, 112.904, -89.3696, '9/2', '+', 4.29, 'inf', 5.39, 0, 'NULL', 'NULL', 'NULL', 'NULL', 3.65, 0.000037, 3.65, 12, 0.757843),
'In115': ('In', '115', 'Indium', 115, 49, 114.904, -89.5366, '9/2', '+', 95.71, 'inf', 4, -0.0562, 2.1, 6.4, -2.1, 0, 2.02, 0.55, 2.57, 202, 0.0317037),
'In': ('In', 'nat', 'Indium', 115, 49, 114.818, 'NULL', 'NULL', 'NULL', 100, 'inf', 4.065, -0.0539, 'NULL', 'NULL', 'NULL', 'NULL', 2.08, 0.54, 2.62, 193.8, 0),
'Sn112': ('Sn', '112', 'Tin', 112, 50, 111.905, -88.6613, '0', '+', 0.97, 'inf', 6, 0, 'NULL', 'NULL', 0, 0, 4.5, 0, 4.5, 1, 0.0709827),
'Sn114': ('Sn', '114', 'Tin', 114, 50, 113.903, -90.5609, '0', '+', 0.66, 'inf', 6, 0, 'NULL', 'NULL', 0, 0, 4.8, 0, 4.8, 0.114, 0.0709827),
'Sn115': ('Sn', '115', 'Tin', 115, 50, 114.903, -90.036, '1/2', '+', 0.34, 'inf', 6, 0, 'NULL', 'NULL', 'NULL', 'NULL', 4.5, 0.3, 4.8, 30, 0.0709827),
'Sn116': ('Sn', '116', 'Tin', 116, 50, 115.902, -91.5281, '0', '+', 14.54, 'inf', 6.1, 0, 'NULL', 'NULL', 0, 0, 4.42, 0, 4.42, 0.14, 0.0397574),
'Sn117': ('Sn', '117', 'Tin', 117, 50, 116.903, -90.4, '1/2', '+', 7.68, 'inf', 6.59, 0, 0.22, -0.23, 0.19, 0, 5.28, 0.3, 5.6, 2.3, 0.120707),
'Sn118': ('Sn', '118', 'Tin', 118, 50, 117.902, -91.6561, '0', '+', 24.22, 'inf', 6.23, 0, 'NULL', 'NULL', 0, 0, 4.63, 0, 4.63, 0.22, 0.00160707),
'Sn119': ('Sn', '119', 'Tin', 119, 50, 118.903, -90.0684, '1/2', '+', 8.59, 'inf', 6.28, 0, 0.14, 0, 0.06, 0, 4.71, 0.3, 5, 2.2, 0.0177487),
'Sn120': ('Sn', '120', 'Tin', 120, 50, 119.902, -91.1051, '0', '+', 32.58, 'inf', 6.67, 0, 'NULL', 'NULL', 0, 0, 5.29, 0, 5.29, 0.14, 0.148082),
'Sn122': ('Sn', '122', 'Tin', 122, 50, 121.903, -89.946, '0', '+', 4.63, 'inf', 5.93, 0, 'NULL', 'NULL', 0, 0, 4.14, 0, 4.14, 0.18, 0.0925333),
'Sn124': ('Sn', '124', 'Tin', 124, 50, 123.905, -88.2367, '0', '+', 5.79, 'inf', 5.79, 0, 'NULL', 'NULL', 0, 0, 4.48, 0, 4.48, 0.133, 0.134876),
'Sn': ('Sn', 'nat', 'Tin', 119, 50, 118.71, 'NULL', 'NULL', 'NULL', 100, 'inf', 6.225, 0, 'NULL', 'NULL', 'NULL', 'NULL', 4.871, 0.022, 4.892, 0.626, 0),
'Sb121': ('Sb', '121', 'Antimony', 121, 51, 'NULL', -89.5951, '5/2', '+', 57.21, 'inf', 5.71, 0, 5.7, 5.8, -0.05, 'NULL', 4.1, 0.0003, 4.1003, 5.75, 0.0509011),
'Sb123': ('Sb', '123', 'Antimony', 123, 51, 'NULL', -89.2241, '7/2', '+', 42.79, 'inf', 5.38, 0, 5.2, 5.4, -0.1, 0, 3.64, 0.001, 3.641, 3.8, 0.067059),
'Sb': ('Sb', 'nat', 'Antimony', 122, 51, 121.76, 0, 'NULL', 'NULL', 100, 'inf', 5.57, 0, 'NULL', 'NULL', 'NULL', 'NULL', 3.9, 0, 3.9, 4.91, 0),
'Te120': ('Te', '120', 'Tellurium', 120, 52, 119.904, -89.4046, '0', '+', 0.09, 'inf', 5.3, 0, 'NULL', 'NULL', 0, 0, 3.5, 0, 3.5, 2.3, 0.129327),
'Te122': ('Te', '122', 'Tellurium', 122, 52, 121.903, -90.314, '0', '+', 2.55, 'inf', 3.8, 0, 'NULL', 'NULL', 0, 0, 1.8, 0, 1.8, 3.4, 0.55242),
'Te123': ('Te', '123', 'Tellurium', 123, 52, 122.904, -89.1719, '1/2', '+', 0.89, 'inf', -0.05, -0.116, -1.2, 3.5, -2.04, 0, 0.002, 0.52, 0.52, 418, 0.999505),
'Te124': ('Te', '124', 'Tellurium', 124, 52, 123.903, -90.5245, '0', '+', 4.74, 'inf', 7.95, 0, 'NULL', 'NULL', 0, 0, 8, 0, 8, 6.8, 0.959014),
'Te125': ('Te', '125', 'Tellurium', 125, 52, 124.904, -89.0222, '1/2', '+', 7.07, 'inf', 5.01, 0, 4.9, 5.2, -0.26, 0, 3.17, 0.008, 3.18, 1.55, 0.222001),
'Te126': ('Te', '126', 'Tellurium', 126, 52, 125.903, -90.0646, '0', '+', 18.84, 'inf', 5.55, 0, 'NULL', 'NULL', 0, 0, 3.88, 0, 3.88, 1.04, 0.0452508),
'Te128': ('Te', '128', 'Tellurium', 128, 52, 127.904, -88.9921, '0', '+', 31.74, 'inf', 5.88, 0, 'NULL', 'NULL', 0, 0, 4.36, 0, 4.36, 0.215, 0.0716624),
'Te130': ('Te', '130', 'Tellurium', 130, 52, 129.906, -87.3514, '0', '+', 34.08, 'inf', 6.01, 0, 'NULL', 'NULL', 0, 0, 4.55, 0, 4.55, 0.29, 0.119573),
'Te': ('Te', 'nat', 'Tellurium', 128, 52, 127.6, 0, 'NULL', 'NULL', 100, 'inf', 5.68, 0, 'NULL', 'NULL', 'NULL', 'NULL', 4.23, 0.09, 4.32, 4.7, 0),
'I': ('I', 'nat', 'Iodine', 127, 53, 126.904, -88.9831, '5/2', '+', 100, 'inf', 5.28, 0, 6.6, 3.4, 1.58, 0, 3.5, 0.31, 3.81, 6.15, 0),
'Xe124': ('Xe', '124', 'Xenon', 124, 54, 123.906, -87.6601, '0', '+', 0.0952, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 165, 1),
'Xe126': ('Xe', '126', 'Xenon', 126, 54, 125.904, -89.1685, '0', '+', 0.089, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 3.5, 1),
'Xe128': ('Xe', '128', 'Xenon', 128, 54, 127.904, -89.86, '0', '+', 1.9102, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 8, 1),
'Xe129': ('Xe', '129', 'Xenon', 129, 54, 128.905, -88.6974, '1/2', '+', 26.4006, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 21, 1),
'Xe130': ('Xe', '130', 'Xenon', 130, 54, 129.904, -89.8817, '0', '+', 4.071, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 26, 1),
'Xe131': ('Xe', '131', 'Xenon', 131, 54, 130.905, -88.4152, '3/2', '+', 21.332, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 85, 1),
'Xe132': ('Xe', '132', 'Xenon', 132, 54, 131.904, -89.2805, '0', '+', 26.9086, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 0.45, 1),
'Xe134': ('Xe', '134', 'Xenon', 134, 54, 133.905, -88.1245, '0', '+', 10.4357, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 0.265, 1),
'Xe136': ('Xe', '136', 'Xenon', 136, 54, 135.907, -86.4251, '0', '+', 8.8573, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 0.26, 1),
'Xe': ('Xe', 'nat', 'Xenon', 131, 54, 131.293, 'NULL', 'NULL', 'NULL', 100, 'inf', 4.69, 0, 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0, 'NULL', 165, 0),
'Cs': ('Cs', 'nat', 'Caesium', 133, 55, 132.905, -88.071, '7/2', '+', 100, 'inf', 5.42, 0, 'NULL', 'NULL', 1.29, 0, 3.69, 0.21, 3.9, 29, 0),
'Ba130': ('Ba', '130', 'Barium', 130, 56, 129.906, -87.2616, '0', '+', 0.106, 'NULL', -3.6, 0, 'NULL', 'NULL', 0, 0, 1.6, 0, 1.6, 30, 'NULL'),
'Ba132': ('Ba', '132', 'Barium', 132, 56, 131.905, -88.4348, '0', '+', 0.101, 'inf', 7.8, 0, 'NULL', 'NULL', 0, 0, 7.6, 0, 7.6, 7, 'NULL'),
'Ba134': ('Ba', '134', 'Barium', 134, 56, 133.904, -88.9499, '0', '+', 2.417, 'inf', 5.7, 0, 'NULL', 'NULL', 0, 0, 4.08, 0, 4.08, 2, 'NULL'),
'Ba135': ('Ba', '135', 'Barium', 135, 56, 134.906, -87.8505, '3/2', '+', 6.592, 'inf', 4.66, 0, 'NULL', 'NULL', 'NULL', 'NULL', 2.74, 0.5, 3.2, 5.8, 'NULL'),
'Ba136': ('Ba', '136', 'Barium', 136, 56, 135.905, -88.8869, '0', '+', 7.854, 'inf', 4.9, 0, 'NULL', 'NULL', 0, 0, 3.03, 0, 3.03, 0.68, 'NULL'),
'Ba137': ('Ba', '137', 'Barium', 137, 56, 136.906, -87.7212, '3/2', '+', 11.232, 'inf', 6.82, 0, 'NULL', 'NULL', 'NULL', 'NULL', 5.86, 0.5, 6.4, 3.6, 'NULL'),
'Ba138': ('Ba', '138', 'Barium', 138, 56, 137.905, -88.2616, '0', '+', 1.6987, 'inf', 4.83, 0, 'NULL', 'NULL', 0, 0, 2.94, 0, 2.94, 0.27, 'NULL'),
'Ba': ('Ba', 'nat', 'Barium', 137, 56, 137.327, 'NULL', 'NULL', 'NULL', 100, 'inf', 5.07, 0, 'NULL', 'NULL', 0, 0, 3.23, 0.15, 3.38, 1.1, 0),
'La138': ('La', '138', 'Lanthanum', 138, 57, 137.907, -86.5247, '5', '+', 0.08881, 'inf', 8, 0, 'NULL', 'NULL', 'NULL', 'NULL', 8, 0.5, 8.5, 57, 0.0574041),
'La139': ('La', '139', 'Lanthanum', 139, 57, 138.906, -87.2314, '7/2', '+', 99.9119, 'inf', 8.24, 0, 11.4, 4.5, 3, 0, 8.53, 11.13, 9.66, 8.93, 0),
'La': ('La', 'nat', 'Lanthanum', 139, 57, 138.905, 'NULL', 'NULL', 'NULL', 100, 'inf', 8.24, 0, 'NULL', 'NULL', 'NULL', 'NULL', 8.53, 1.13, 9.66, 8.97, 0),
'Ce136': ('Ce', '136', 'Cerium', 136, 58, 135.907, -86.4683, '0', '+', 0.185, 'inf', 5.76, 0, 'NULL', 'NULL', 0, 0, 4.23, 0, 4.23, 7.3, 0.416297),
'Ce138': ('Ce', '138', 'Cerium', 138, 58, 137.906, -87.5685, '0', '+', 0.251, 'inf', 6.65, 0, 'NULL', 'NULL', 0, 0, 5.64, 0, 5.64, 1.1, 0.887785),
'Ce140': ('Ce', '140', 'Cerium', 140, 58, 139.905, -88.0833, '0', '+', 88.45, 'inf', 4.81, 0, 'NULL', 'NULL', 0, 0, 2.94, 0, 2.94, 0.57, 0.0123583),
'Ce142': ('Ce', '142', 'Cerium', 142, 58, 141.909, -84.5385, '0', '+', 11.114, 'inf', 4.72, 0, 'NULL', 'NULL', 0, 0, 2.84, 0, 2.84, 0.95, 0.0489721),
'Ce': ('Ce', 'nat', 'Cerium', 140, 58, 140.116, 'NULL', 'NULL', 'NULL', 100, 'inf', 4.84, 0, 'NULL', 'NULL', 'NULL', 'NULL', 2.94, 0, 2.94, 0.63, 0),
'Pr': ('Pr', 'nat', 'Praseodymium', 141, 59, 140.908, -86.0209, '5/2', '+', 100, 'inf', 4.58, 0, 'NULL', 'NULL', -0.055, 0, 2.64, 0.01, 2.66, 11.5, 0),
'Nd142': ('Nd', '142', 'Neodymium', 142, 60, 141.908, -85.9552, '0', '+', 27.152, 'inf', 7.7, 0, 'NULL', 'NULL', 0, 0, 7.5, 0, 7.5, 18.7, 0.00260247),
'Nd143': ('Nd', '143', 'Neodymium', 143, 60, 142.91, -84.0074, '7/2', '+', 12.174, 'inf', 14, 0, 'NULL', 'NULL', -21, 0, 25, 55, 80, 337, 2.31439),
'Nd144': ('Nd', '144', 'Neodymium', 144, 60, 143.91, -83.7532, '0', '+', 23.798, 'inf', 2.8, 0, 'NULL', 'NULL', 0, 0, 1, 0, 1, 3.6, 0.867424),
'Nd145': ('Nd', '145', 'Neodymium', 145, 60, 144.913, -81.4371, '7/2', '+', 8.293, 'inf', 14, 0, 'NULL', 'NULL', 'NULL', 'NULL', 25, 5, 30, 42, 2.31439),
'Nd146': ('Nd', '146', 'Neodymium', 146, 60, 145.913, -80.931, '0', '+', 17.189, 'inf', 8.7, 0, 'NULL', 'NULL', 0, 0, 9.5, 0, 9.5, 1.4, 0.279929),
'Nd148': ('Nd', '148', 'Neodymium', 148, 60, 147.917, -77.4134, '0', '+', 5.756, 'inf', 5.7, 0, 'NULL', 'NULL', 0, 0, 4.1, 0, 4.1, 2.5, 0.450589),
'Nd150': ('Nd', '150', 'Neodymium', 150, 60, 149.921, -73.6897, '0', '+', 5.638, 'inf', 5.28, 5.28, 'NULL', 'NULL', 0, 0, 3.5, 0, 3.5, 1.2, 0.0571445),
'Nd': ('Nd', 'nat', 'Neodymium', 144, 60, 144.242, 'NULL', 'NULL', 'NULL', 100, 'inf', 7.69, 0, 'NULL', 'NULL', 'NULL', 'NULL', 7.43, 9.2, 16.6, 50.5, 0),
'Pm': ('Pm', 'nat', 'Promethium', 147, 61, 'NULL', -79.0479, '7/2', '+', 100, '2.6234y', 12.6, 0, 'NULL', 'NULL', 'NULL', 'NULL', 20, 1.3, 21.3, 168.4, 0),
'Sm144': ('Sm', '144', 'Samarium', 144, 62, 143.912, -81.972, '0', '+', 3.07, 'inf', -3, 0, 'NULL', 'NULL', 0, 0, 1, 0, 1, 0.7, 2.30579),
'Sm147': ('Sm', '147', 'Samarium', 147, 62, 146.915, -79.2721, '7/2', '-', 14.99, 'inf', 14, 0, 'NULL', 'NULL', -11, 0, 25, 14, 39, 57, 70.9927),
'Sm148': ('Sm', '148', 'Samarium', 148, 62, 147.915, -79.3422, '0', '+', 11.24, 'inf', -3, 0, 'NULL', 'NULL', 0, 0, 1, 0, 1, 2.4, 2.30579),
'Sm149': ('Sm', '149', 'Samarium', 149, 62, 148.917, -77.1419, '7/2', '-', 13.82, 'inf', 18.7, -11.7, 'NULL', 'NULL', -31.4, -10.3, 63.5, 137, 200, 42080, 177.725),
'Sm150': ('Sm', '150', 'Samarium', 150, 62, 149.917, -77.0573, '0', '+', 7.38, 'inf', 14, 0, 'NULL', 'NULL', 0, 0, 25, 0, 25, 104, 70.9927),
'Sm152': ('Sm', '152', 'Samarium', 152, 62, 151.92, -74.7688, '0', '+', 26.75, 'inf', -5, 0, 'NULL', 'NULL', 0, 0, 3.1, 0, 3.1, 206, 8.18274),
'Sm154': ('Sm', '154', 'Samarium', 154, 62, 153.922, -72.4616, '0', '+', 22.75, 'inf', -8, 0, 'NULL', 'NULL', 0, 0, 11, 0, 11, 8.4, 22.5078),
'Sm': ('Sm', 'nat', 'Samarium', 150, 62, 150.36, 'NULL', 'NULL', 'NULL', 100, 'inf', 0, -1.65, 'NULL', 'NULL', 'NULL', 'NULL', 0.422, 39, 39.4, 5922, 0),
'Eu151': ('Eu', '151', 'Europium', 151, 63, 150.92, -74.6591, '5/2', '+', 47.81, 'inf', 6.92, 2.53, 'NULL', 'NULL', -4.5, 2.14, 5.5, 3.1, 8.6, 9100, 0.829235),
'Eu153': ('Eu', '153', 'Europium', 153, 63, 152.921, -73.3735, '5/2', '+', 52.19, 'inf', 8.22, 0, 'NULL', 'NULL', -3.2, 0, 8.5, 1.3, 9.8, 312, 1.27675),
'Eu': ('Eu', 'nat', 'Europium', 152, 63, 151.964, 'NULL', 'NULL', 'NULL', 100, 'inf', 5.3, -1.26, 'NULL', 'NULL', 'NULL', 'NULL', 6.57, 2.5, 9.2, 4530, 0),
'Gd152': ('Gd', '152', 'Gadolinium', 152, 64, 151.92, -74.7142, '0', '+', 0.2, 'inf', 10, 0, 'NULL', 'NULL', 0, 0, 13, 0, 13, 735, 0.644435),
'Gd154': ('Gd', '154', 'Gadolinium', 154, 64, 153.921, -73.7132, '0', '+', 2.18, 'inf', 10, 0, 'NULL', 'NULL', 0, 0, 13, 0, 13, 85, 0.644435),
'Gd155': ('Gd', '155', 'Gadolinium', 155, 64, 154.923, -72.0771, '3/2', '-', 14.8, 'inf', 13.8, -17, 'NULL', 'NULL', -5, -13.16, 40.8, 25, 66, 61100, 0.704722),
'Gd156': ('Gd', '156', 'Gadolinium', 156, 64, 155.922, -72.5422, '0', '+', 20.47, 'inf', 6.3, 0, 'NULL', 'NULL', 0, 0, 5, 0, 5, 1.5, 0.858876),
'Gd157': ('Gd', '157', 'Gadolinium', 157, 64, 156.924, -70.8307, '3/2', '-', 15.65, 'inf', -1.14, -72, 'NULL', 'NULL', 5, -55.8, 650, 349, 104.4, 259000, 17.4371),
'Gd158': ('Gd', '158', 'Gadolinium', 158, 64, 157.924, -70.6967, '0', '+', 24.84, 'inf', 9, 0, 'NULL', 'NULL', 0, 0, 10, 0, 10, 2.2, 0.711992),
'Gd160': ('Gd', '160', 'Gadolinium', 160, 64, 159.927, -67.9486, '0', '+', 21.58, 'inf', 9.15, 0, 'NULL', 'NULL', 0, 0, 10.52, 0, 10.52, 0.77, 0.702312),
'Gd': ('Gd', 'nat', 'Gadolinium', 157, 64, 157.25, 'NULL', 'NULL', 'NULL', 100, 'inf', 9.5, -13.82, 'NULL', 'NULL', 'NULL', 'NULL', 29.3, 151, 180, 59700, 0),
'Tb': ('Tb', 'nat', 'Terbium', 159, 65, 158.925, -69.539, '3/2', '+', 100, 'inf', 7.34, 0, 6.8, 8.1, -0.17, 0, 6.48, 0.004, 6.48, 23.4, 0),
'Dy156': ('Dy', '156', 'Dysprosium', 156, 66, 155.924, -70.5298, '0', '+', 0.056, 'inf', 6.1, 0, 'NULL', 'NULL', 0, 0, 4.7, 0, 4.7, 33, 0.869752),
'Dy158': ('Dy', '158', 'Dysprosium', 158, 66, 157.924, -70.4121, '0', '+', 0.095, 'inf', 6.7, 0, 'NULL', 'NULL', 0, 0, 5, 0, 5, 43, 0.84287),
'Dy160': ('Dy', '160', 'Dysprosium', 160, 66, 159.925, -69.6781, '0', '+', 2.329, 'inf', 6.7, 0, 'NULL', 'NULL', 0, 0, 5.6, 0, 5.6, 56, 0.84287),
'Dy161': ('Dy', '161', 'Dysprosium', 161, 66, 160.927, -68.0611, '5/2', '+', 18.889, 'inf', 10.3, 0, 14.5, 4.2, -0.17, 0, 13.3, 3, 16, 600, 0.628648),
'Dy162': ('Dy', '162', 'Dysprosium', 162, 66, 161.927, -68.1868, '0', '+', 25.475, 'inf', -1.4, 0, 'NULL', 'NULL', 0, 0, 0.25, 0, 0.25, 194, 0.993139),
'Dy163': ('Dy', '163', 'Dysprosium', 163, 66, 162.929, -66.3865, '5/2', '-', 24.896, 'inf', 5, 0, 6.1, 3.5, 1.3, 0, 3.1, 0.21, 3.3, 124, 0.912491),
'Dy164': ('Dy', '164', 'Dysprosium', 164, 66, 163.929, -65.9733, '0', '+', 28.26, 'inf', 49.4, -0.79, 'NULL', 'NULL', 0, 0, 307, 0, 307, 2840, 7.54428),
'Dy': ('Dy', 'nat', 'Dysprosium', 163, 66, 162.5, 'NULL', 'NULL', 'NULL', 100, 'inf', 16.9, 0.276, 'NULL', 'NULL', 'NULL', 'NULL', 35.9, 54.4, 90.3, 994, 0),
'Ho': ('Ho', 'nat', 'Holmium', 165, 67, 164.93, -64.9046, '7/2', '-', 100, 'inf', 8.44, 0, 6.9, 10.3, -1.69, 0, 8.06, 0.36, 8.42, 64.7, 0),
'Er162': ('Er', '162', 'Erbium', 162, 68, 161.929, -66.3426, '0', '+', 0.139, 'inf', 9.01, 0, 'NULL', 'NULL', 0, 0, 9.7, 0, 9.7, 19, 0.337749),
'Er164': ('Er', '164', 'Erbium', 164, 68, 163.929, -65.9496, '0', '+', 1.601, 'inf', 7.95, 0, 'NULL', 'NULL', 0, 0, 8.4, 0, 8.4, 13, 0.0415002),
'Er166': ('Er', '166', 'Erbium', 166, 68, 165.93, -64.9316, '0', '+', 33.503, 'inf', 10.51, 0, 'NULL', 'NULL', 0, 0, 14.1, 0, 14.1, 19.6, 0.820248),
'Er167': ('Er', '167', 'Erbium', 167, 68, 166.932, -63.2967, '7/2', '+', 22.869, 'inf', 3.06, 0, 5.3, 0, 2.6, 0, 1.1, 0.13, 1.2, 659, 0.845699),
'Er169': ('Er', '168', 'Erbium', 168, 68, 167.932, -62.9967, '0', '+', 26.978, 'inf', 7.43, 0, 'NULL', 'NULL', 0, 0, 6.9, 0, 6.9, 2.74, 0.0902905),
'Er170': ('Er', '170', 'Erbium', 170, 68, 169.935, -60.1146, '0', '+', 14.91, 'inf', 9.61, 0, 'NULL', 'NULL', 0, 0, 11.6, 0, 11.6, 5.8, 0.52185),
'Er': ('Er', 'nat', 'Erbium', 167, 68, 167.259, 'NULL', 'NULL', 'NULL', 100, 'inf', 7.79, 0, 'NULL', 'NULL', 'NULL', 'NULL', 7.63, 1.1, 8.7, 159, 0),
'Tm': ('Tm', 'nat', 'Thulium', 169, 69, 168.934, -61.28, '1/2', '+', 100, 'inf', 9.61, 0, 'NULL', 'NULL', 0.9, 0, 6.28, 0.1, 6.38, 100, 0),
'Yb168': ('Yb', '168', 'Ytterbium', 168, 70, 167.934, -61.5746, '0', '+', 0.123, 'inf', -4.07, -0.62, 'NULL', 'NULL', 0, 0, 2.13, 0, 2.13, 2230, 0.889945),
'Yb170': ('Yb', '170', 'Ytterbium', 170, 70, 169.935, -60.769, '0', '+', 2.982, 'inf', 6.8, 0, 'NULL', 'NULL', 0, 0, 5.8, 0, 5.8, 11.4, 0.699756),
'Yb171': ('Yb', '171', 'Ytterbium', 171, 70, 170.936, -59.3121, '1/2', '-', 14.09, 'inf', 9.7, 0, 9.5, 19.4, -5.59, 0, 11.7, 3.9, 15.6, 48.6, 0.389058),
'Yb172': ('Yb', '172', 'Ytterbium', 172, 70, 171.936, -59.2603, '0', '+', 21.68, 'inf', 9.5, 0, 'NULL', 'NULL', 0, 0, 11.2, 0, 11.2, 0.8, 0.413992),
'Yb173': ('Yb', '173', 'Ytterbium', 173, 70, 172.938, -57.5563, '5/2', '-', 16.103, 'inf', 9.56, 0, 2.5, 13.3, -5.3, 0, 11.5, 3.5, 15, 17.1, 0.406566),
'Yb174': ('Yb', '174', 'Ytterbium', 174, 70, 173.939, -56.9496, '0', '+', 32.026, 'inf', 19.2, 0, 'NULL', 'NULL', 0, 0, 46.8, 0, 46.8, 69.4, 1.39364),
'Yb176': ('Yb', '176', 'Ytterbium', 176, 70, 175.943, -53.4941, '0', '+', 12.996, 'inf', 8.7, 0, 'NULL', 'NULL', 0, 0, 9.6, 0, 9.6, 2.85, 0.508532),
'Yb': ('Yb', 'nat', 'Ytterbium', 173, 70, 173.054, 'NULL', 'NULL', 'NULL', 100, 'inf', 12.41, 0, 'NULL', 'NULL', 'NULL', 'NULL', 19.42, 4, 23.4, 34.8, 0),
'Lu175': ('Lu', '175', 'Lutetium', 175, 71, 174.941, -55.1707, '7/2', '+', 97.401, 'inf', 7.28, 0, 'NULL', 'NULL', -2.2, 0, 6.59, 0.6, 7.2, 21, 0.0195117),
'Lu176': ('Lu', '176', 'Lutetium', 176, 71, 175.943, -53.3874, '7', '-', 2.599, 'inf', 6.1, -0.57, 'NULL', 'NULL', -3, 0.61, 4.7, 1.2, 5.9, 2065, 0.277954),
'Lu': ('Lu', 'nat', 'Lutetium', 175, 71, 174.967, 'NULL', 'NULL', 'NULL', 100, 'inf', 7.21, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.53, 0.7, 7.2, 74, 0),
'Hf174': ('Hf', '174', 'Hafnium', 174, 72, 173.94, -55.8466, '0', '+', 0.16, 'inf', 10.9, 0, 'NULL', 'NULL', 0, 0, 15, 0, 15, 561, 0.967936),
'Hf176': ('Hf', '176', 'Hafnium', 176, 72, 175.941, -54.5775, '0', '+', 5.26, 'inf', 6.61, 0, 'NULL', 'NULL', 0, 0, 5.5, 0, 5.5, 23.5, 0.276296),
'Hf177': ('Hf', '177', 'Hafnium', 177, 72, 176.943, -52.8896, '7/2', '+', 18.6, 'inf', 0.8, 0, 'NULL', 'NULL', -0.9, 0, 0.1, 0.1, 0.2, 373, 0.989399),
'Hf178': ('Hf', '178', 'Hafnium', 178, 72, 177.944, -52.4443, '0', '+', 27.28, 'inf', 5.9, 0, 'NULL', 'NULL', 0, 0, 4.4, 0, 4.4, 84, 0.423417),
'Hf179': ('Hf', '179', 'Hafnium', 179, 72, 178.946, -50.4719, '9/2', '+', 13.62, 'inf', 7.46, 0, 'NULL', 'NULL', -1.06, 0, 7, 0.14, 7.1, 41, 0.0782023),
'Hf180': ('Hf', '180', 'Hafnium', 180, 72, 179.947, -49.7884, '0', '+', 35.08, 'inf', 13.2, 0, 'NULL', 'NULL', 0, 0, 21.9, 0, 21.9, 13.04, 1.88606),
'Hf': ('Hf', 'nat', 'Hafnium', 180, 72, 178.49, 'NULL', 'NULL', 'NULL', 100, 'inf', 7.77, 0, 'NULL', 'NULL', 'NULL', 'NULL', 7.6, 2.6, 10.2, 104.1, 0),
'Ta180': ('Ta', '180', 'Tantalum', 180, 73, 179.947, -48.8591, '9', '-', 0.01201, 'inf', 7, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6.2, 0.5, 7, 563, 0.0262188),
'Ta181': ('Ta', '181', 'Tantalum', 181, 73, 180.948, -48.4416, '7/2', '+', 99.988, 'inf', 6.91, 0, 'NULL', 'NULL', -0.29, 0, 6, 0.01, 6.01, 20.5, 0),
'Ta': ('Ta', 'nat', 'Tantalum', 181, 73, 180.948, 'NULL', 'NULL', 'NULL', 100, 'inf', 6.91, 0, 'NULL', 'NULL', 'NULL', 'NULL', 6, 0.01, 6.01, 20.6, 0),
'W180': ('W', '180', 'Tungsten', 180, 74, 179.947, -49.6445, '0', '+', 0.12, 'inf', 5, 0, 'NULL', 'NULL', 0, 0, 3, 0, 3, 30, 0.105704),
'W182': ('W', '182', 'Tungsten', 182, 74, 181.948, -48.2475, '0', '+', 26.5, 'inf', 7.04, 0, 'NULL', 'NULL', 0, 0, 6.1, 0, 6.1, 20.75, 1.19202),
'W183': ('W', '183', 'Tungsten', 183, 74, 182.95, -46.367, '1/2', '-', 14.31, 'inf', 6.59, 0, 6.3, 7, -0.3, 0, 5.36, 0.3, 5.7, 10.1, 0.920745),
'W184': ('W', '184', 'Tungsten', 184, 74, 183.951, -45.7073, '0', '+', 30.64, 'inf', 7.55, 0, 'NULL', 'NULL', 0, 0, 7.03, 0, 7.03, 1.7, 1.52112),
'W186': ('W', '186', 'Tungsten', 186, 74, 185.954, -42.5095, '0', '+', 28.43, 'inf', -0.73, 0, 'NULL', 'NULL', 0, 0, 0.065, 0, 0.065, 87.9, 0.976431),
'W': ('W', 'nat', 'Tungsten', 184, 74, 183.84, 'NULL', 'NULL', 'NULL', 100, 'inf', 4.755, 0, 'NULL', 'NULL', 'NULL', 'NULL', 2.97, 1.63, 4.6, 18.3, 0),
'Re185': ('Re', '185', 'Rhenium', 185, 75, 184.953, -43.8222, '5/2', '+', 37.4, 'inf', 9, 0, 'NULL', 'NULL', -2, 0, 10.2, 0.5, 10.7, 112, 0.0430057),
'Re187': ('Re', '187', 'Rhenium', 187, 75, 186.956, -41.2157, '5/2', '+', 62.6, 'inf', 9.3, 0, 'NULL', 'NULL', -2.8, 0, 10.9, 1, 11.9, 76.4, 0.0218573),
'Re': ('Re', 'nat', 'Rhenium', 186, 75, 186.207, 'NULL', 'NULL', 'NULL', 100, 'inf', 9.2, 0, 'NULL', 'NULL', 'NULL', 'NULL', 10.6, 0.9, 11.9, 89.7, 0),
'Os184': ('Os', '184', 'Osmium', 184, 76, 183.952, -44.2561, '0', '+', 0.02, 'inf', 10.2, 0, 'NULL', 'NULL', 0, 0, 13, 0, 13, 3000, 0.0912743),
'Os186': ('Os', '186', 'Osmium', 186, 76, 185.954, -42.9995, '0', '+', 1.59, 'inf', 12, 0, 'NULL', 'NULL', 0, 0, 17, 0, 17, 80, 0.257752),
'Os187': ('Os', '187', 'Osmium', 187, 76, 186.956, -41.2182, '1/2', '-', 1.96, 'inf', 10, 0, 'NULL', 'NULL', 'NULL', 'NULL', 13, 0.3, 13, 320, 0.126561),
'Os188': ('Os', '188', 'Osmium', 188, 76, 187.956, -41.1364, '0', '+', 13.24, 'inf', 7.8, 0, 'NULL', 'NULL', 0, 0, 7.3, 0, 7.3, 4.7, 0.4686),
'Os189': ('Os', '189', 'Osmium', 189, 76, 188.958, -38.9854, '3/2', '-', 16.15, 'inf', 11, 0, 'NULL', 'NULL', 'NULL', 'NULL', 14.4, 0.5, 14.9, 25, 0.0568609),
'Os190': ('Os', '190', 'Osmium', 190, 76, 189.958, -38.7063, '0', '+', 26.26, 'inf', 11.4, 0, 'NULL', 'NULL', 0, 0, 15.2, 0, 15.2, 25, 0.135121),
'Os192': ('Os', '192', 'Osmium', 192, 76, 191.961, -35.8805, '0', '+', 40.78, 'inf', 11.9, 0, 'NULL', 'NULL', 0, 0, 16.6, 0, 16.6, 2, 0.236877),
'Os': ('Os', 'nat', 'Osmium', 190, 76, 190.23, 'NULL', 'NULL', 'NULL', 100, 'inf', 10.7, 0, 'NULL', 'NULL', 'NULL', 'NULL', 14.4, 0.3, 14.7, 16, 0),
'Ir191': ('Ir', '191', 'Iridium', 191, 77, 190.961, -36.7064, '3/2', '+', 37.3, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 954, 1),
'Ir193': ('Ir', '193', 'Iridium', 193, 77, 192.963, -34.5338, '3/2', '+', 62.7, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 111, 1),
'Ir': ('Ir', 'nat', 'Iridium', 192, 77, 192.217, 'NULL', 'NULL', 'NULL', 100, 'inf', 10.6, 0, 'NULL', 'NULL', 'NULL', 'NULL', 14.1, 0, 14, 425, 0),
'Pt190': ('Pt', '190', 'Platinium', 190, 78, 189.96, -37.3234, '0', '+', 0.012, 'inf', 9, 0, 'NULL', 'NULL', 0, 0, 10, 0, 10, 152, 0.121094),
'Pt192': ('Pt', '192', 'Platinium', 192, 78, 191.961, -36.2929, '0', '+', 0.782, 'inf', 9.9, 0, 'NULL', 'NULL', 0, 0, 12.3, 0, 12.3, 10, 0.0634766),
'Pt194': ('Pt', '194', 'Platinium', 194, 78, 193.963, -34.7631, '0', '+', 32.86, 'inf', 10.55, 0, 'NULL', 'NULL', 0, 0, 14, 0, 14, 1.44, 0.207709),
'Pt195': ('Pt', '195', 'Platinium', 195, 78, 194.965, -32.7968, '1/2', '-', 33.78, 'inf', 8.91, 0, 'NULL', 'NULL', 1, 0, 9.8, 0.13, 9.9, 27.5, 0.138584),
'Pt196': ('Pt', '196', 'Platinium', 196, 78, 195.965, -32.6474, '0', '+', 25.21, 'inf', 9.89, 0, 'NULL', 'NULL', 0, 0, 12.3, 0, 12.3, 0.72, 0.0613292),
'Pt198': ('Pt', '198', 'Platinium', 198, 78, 197.968, -29.9077, '0', '+', 7.36, 'inf', 7.8, 0, 'NULL', 'NULL', 0, 0, 7.6, 0, 7.6, 3.66, 0.339844),
'Pt1': ('Pt', 'nat', 'Platinium', 195, 78, 195.084, 'NULL', 'NULL', 'NULL', 100, 'inf', 9.6, 0, 'NULL', 'NULL', 'NULL', 'NULL', 11.58, 0.13, 11.71, 10.3, 0),
'Au': ('Au', 'nat', 'Gold', 197, 79, 196.967, -31.1411, '3/2', '+', 100, 'inf', 7.9, 0, 6.26, 9.9, -1.76, 0, 7.32, 0.43, 7.75, 98.65, 0),
'Hg196': ('Hg', '196', 'Mercury', 196, 80, 195.966, -31.8267, '0', '+', 0.15, 'inf', 30.3, -0.85, 'NULL', 'NULL', 0, 0, 115, 0, 115, 3080, 4.79203),
'Hg198': ('Hg', '198', 'Mercury', 198, 80, 197.967, -30.9544, '0', '+', 9.97, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 2.03, 1),
'Hg199': ('Hg', '199', 'Mercury', 199, 80, 198.968, -29.5471, '1/2', '-', 16.87, 'inf', 16.9, -0.6, 'NULL', 'NULL', -15.5, 0, 36, 30, 66, 2150, 0.802703),
'Hg200': ('Hg', '200', 'Mercury', 200, 80, 199.968, -29.5041, '0', '+', 23.1, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 60, 1),
'Hg201': ('Hg', '201', 'Mercury', 201, 80, 200.97, -27.6633, '3/2', '-', 13.18, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 7.8, 1),
'Hg202': ('Hg', '202', 'Mercury', 202, 80, 201.971, -27.3459, '0', '+', 29.86, 'inf', 11.002, 0, 'NULL', 'NULL', 0, 0, 15.2108, 0, 15.2108, 4.89, 0.236961),
'Hg204': ('Hg', '204', 'Mercury', 204, 80, 203.973, -24.6902, '0', '+', 6.87, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 0, 0, 'NULL', 0, 'NULL', 0.43, 1),
'Hg': ('Hg', 'nat', 'Mercury', 200, 80, 200.59, 'NULL', 'NULL', 'NULL', 100, 'inf', 12.595, 0, 'NULL', 'NULL', 'NULL', 'NULL', 20.24, 6.6, 26.84, 372.3, 0),
'Tl203': ('Tl', '203', 'Thallium', 203, 81, 202.972, -25.7612, '1/2', '+', 29.524, 'inf', 8.51, 0, 9.08, 6.62, 1.06, 0, 6.14, 0.14, 6.28, 11.4, 0.0597012),
'Tl205': ('Tl', '205', 'Thallium', 205, 81, 204.974, -23.8206, '1/2', '+', 70.48, 'inf', 8.87, 0, 5.15, 9.43, -0.242, 0, 11.39, 0.007, 11.4, 0.104, 0.0215368),
'Tl': ('Tl', 'nat', 'Thallium', 204, 81, 204.383, 'NULL', 'NULL', 'NULL', 100, 'inf', 8.776, 0, 'NULL', 'NULL', 'NULL', 'NULL', 9.678, 0.21, 9.89, 3.43, 0),
'Pb204': ('Pb', '204', 'Lead', 204, 82, 203.973, -25.1097, '0', '+', 1.4, 'inf', 10.893, 0, 'NULL', 'NULL', 0, 0, 12.3, 0, 12.3, 0.65, 0.342601),
'Pb206': ('Pb', '206', 'Lead', 206, 82, 205.974, -23.7854, '0', '+', 24.1, 'inf', 9.221, 0, 'NULL', 'NULL', 0, 0, 10.68, 0, 10.68, 0.03, 0.0379272),
'Pb207': ('Pb', '207', 'Lead', 207, 82, 206.976, -22.4519, '0', '+', 22.1, 'inf', 9.286, 0, 'NULL', 'NULL', 0.14, 0, 10.82, 0.002, 10.82, 0.699, 0.0243158),
'Pb208': ('Pb', '208', 'Lead', 208, 82, 207.977, -21.7485, '0', '+', 52.4, 'inf', 9.494, 0, 'NULL', 'NULL', 0, 0, 11.34, 0, 11.34, 0.00048, 0.019883),
'Pb': ('Pb', 'nat', 'Lead', 207, 82, 207.2, 'NULL', 'NULL', 'NULL', 100, 'inf', 9.401, 0, 'NULL', 'NULL', 'NULL', 'NULL', 11.115, 0.003, 11.118, 0.171, 0),
'Bi': ('Bi', 'nat', 'Bismuth', 209, 83, 208.98, -18.2585, '9/2', '-', 100, 'inf', 8.53, 0, 8.26, 8.74, 0.22, 0.22, 9.148, 0.0084, 9.156, 0.0338, 0),
'Po': ('Po', 'nat', 'Polonium', 209, 84, 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0),
'At': ('At', 'nat', 'Astatine', 210, 85, 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0),
'Rn': ('Rn', 'nat', 'Radon', 222, 86, 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0),
'Fr': ('Fr', 'nat', 'Francium', 223, 87, 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0),
'Ra': ('Ra', 'nat', 'Radium', 226, 88, 'NULL', 'NULL', '0', '+', 100, '1600y', 10, 0, 'NULL', 'NULL', 0, 0, 13, 0, 13, 12.8, 0),
'Ac': ('Ac', 'nat', 'Actinium', 227, 89, 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0),
'Th': ('Th', 'nat', 'Thorium', 232, 90, 232.038, 35.4483, '0', '+', 100, 'inf', 10.31, 0, 'NULL', 'NULL', 0, 0, 13.36, 0, 13.36, 7.37, 0),
'Pa': ('Pa', 'nat', 'Protactinium', 231, 91, 231.036, 'NULL', '3/2', '-', 100, '159200y', 9.1, 0, 'NULL', 'NULL', 'NULL', 'NULL', 10.4, 0.1, 10.5, 200.6, 0),
'U233': ('U', '233', 'Uranium', 233, 92, 233.04, 36.92, '5/2', '+', 'NULL', '159200y', 10.1, 0, 'NULL', 'NULL', -1, 0, 12.8, 0.1, 12.9, 574.7, 0.439886),
'U234': ('U', '234', 'Uranium', 234, 92, 234.041, 38.1466, '0', '+', 0.0054, 'inf', 12.4, 0, 'NULL', 'NULL', 0, 0, 19.3, 0, 19.3, 100.1, 1.17034),
'U235': ('U', '235', 'Uranium', 235, 92, 235.044, 40.9205, '7/2', '-', 0.7204, 'inf', 10.5, 0, 'NULL', 'NULL', -1.3, 0, 13.78, 0.2, 14, 680.9, 0.556195),
'U238': ('U', '238', 'Uranium', 238, 92, 238.051, 47.3089, '0', '+', 99.2742, 'inf', 0, 0, 'NULL', 'NULL', 0, 0, 8.871, 0, 8.871, 2.68, 1),
'U': ('U', 'nat', 'Uranium', 238, 92, 238.029, 'NULL', 'NULL', 'NULL', 100, 'inf', 8.417, 0, 'NULL', 'NULL', 'NULL', 'NULL', 8.903, 0.005, 8.908, 7.57, 0),
'Np': ('Np', 'nat', 'Neptunium', 237, 93, 'NULL', 44.8733, '5/2', '+', 100, '2144000y', 10.55, 0, 'NULL', 'NULL', 'NULL', 'NULL', 14, 0.5, 14.5, 175.9, 0),
'Pu238': ('Pu', '238', 'Plutonium', 238, 94, 238.05, 46.1647, '0', '+', 'NULL', '87.74y', 14.1, 0, 'NULL', 'NULL', 0, 0, 25, 0, 25, 558, 0),
'Pu239': ('Pu', '239', 'Plutonium', 239, 94, 239.052, 48.5899, '1/2', 'NULL', 'NULL', '24110y', 7.7, 0, 'NULL', 'NULL', -1.3, 0, 7.5, 0.2, 7.7, 1017.3, 0),
'Pu240': ('Pu', '240', 'Plutonium', 240, 94, 240.054, 50.127, '0', '+', 'NULL', '6561y', 3.5, 0, 'NULL', 'NULL', 0, 0, 1.54, 0, 1.54, 289.6, 0),
'Pu242': ('Pu', '242', 'Plutonium', 242, 94, 242.059, 54.7184, '0', '+', 'NULL', '375000y', 8.1, 0, 'NULL', 'NULL', 0, 0, 8.2, 0, 8.2, 18.5, 0),
'Pu': ('Pu', 'nat', 'Plutonium', 244, 94, 'NULL', 'NULL', 'NULL', 'NULL', 100, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0),
'Am': ('Am', 'nat', 'Americium', 243, 95, 'NULL', 'NULL', '5/2', '-', 100, '7370y', 8.3, 0, 'NULL', 'NULL', -2, 0, 8.7, 0.3, 9, 75.3, 0),
'Cm244': ('Cm', '244', 'Curium', 244, 96, 244.063, 58.4537, '0', '+', 'NULL', '18.1y', 9.5, 0, 'NULL', 'NULL', 0, 0, 11.3, 0, 11.3, 16.2, 0),
'Cm246': ('Cm', '246', 'Curium', 246, 96, 246.067, 62.6184, '0', '+', 'NULL', '4706y', 9.3, 0, 'NULL', 'NULL', 0, 0, 10.9, 0, 10.9, 1.36, 0),
'Cm248': ('Cm', '248', 'Curium', 248, 96, 248.072, 67.3922, '0', '+', 'NULL', '348000y', 7.7, 0, 'NULL', 'NULL', 0, 0, 7.5, 0, 7.5, 3, 0),
'Cm': ('Cm', 'nat', 'Curium', 247, 96, 'NULL', 'NULL', 'NULL', 'NULL', 100, 'inf', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 0) }


