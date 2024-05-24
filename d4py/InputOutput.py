import numpy as np  # NumPy library
import sys
import os
import statistics as st
import matplotlib.pyplot as plt

from .d4nifit import find_peaks_in_range, find_minimum_within_range

# --------1---------2---------3---------4---------5---------6---------7---------8
def getRunningParams(parfile):
    """
    This function reads the parameters file for data treatment notebooks.

    As first line of the heading, it writes the filename preceded by "# "
    Then, it prints as many lines as elements contains the list heading.
    Finally, it writes a line for each point with (x,y,error)

    Input:
        - parfile: a string containing the name of the parameters file.

    Output:
      - inputPar (dictionary): Contains all the parameters read from the
        parameters file

    Date: 22/02/2023
    Author: Gabriel Cuello
    --------------------------------------------------------------------------------
    """
    shapes = ["cylinder"]
    containers = ["SiO2", "V", "TiZr", "Nb"]
    absorbers = ["B"]
    environments = ["V"]
    Nelements = 0

    inputPar = {}  # Initialising the dictionary
    with open(parfile, "r") as f:
        lines = f.readlines()

    # Testing the first character of each line non-blank line.
    # If this character is not #, ! or <, stops the program.
    for i in range(len(lines)):
        if len(lines[i]) > 1:
            first = (
                (lines[i][0] != "#") and (lines[i][0] != "!") and (lines[i][0] != "<")
            )
            if first:
                print("Wrong input in line: ", i + 1, " file: ", parfile)
                sys.exit()

    # Loop over all lines in the input file
    for i in range(len(lines)):
        #    print (lines[i][0:-1])
        if (lines[i][0] == "#") or (lines[i][0] == "!"):
            #        print ("Comment line",len(lines[i]))
            pass
        elif len(lines[i]) == 1:
            #        print("Blank line")
            pass
        elif lines[i][0] == "<":
            line = lines[i].split(" ")
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<ins>":
                inputPar[line[0]] = ("Instrument   ", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<pro>":
                inputPar[line[0]] = ("Proposal     ", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<mpr>":
                inputPar[line[0]] = ("Main Proposer", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<lco>":
                inputPar[line[0]] = ("Local Contact", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<dte>":
                hyphen = line[1].find("-")
                inputPar["<sta>"] = ("Starting date", line[1][0:hyphen])
                inputPar["<end>"] = ("Ending date  ", line[1][hyphen + 1 :])
            #                print ("{}: {}"
            #                       .format(inputPar['<sta>'][0],inputPar['<sta>'][1]))
            #                print ("{}: {}"
            #                       .format(inputPar['<end>'][0],inputPar['<end>'][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<cyc>":
                inputPar[line[0]] = ("Cycle        ", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<log>":
                hyphen = line[1].find("-")
                inputPar["<nbk>"] = ("Notebook nbr.", line[1][0:hyphen])
                inputPar["<pag>"] = ("Page         ", line[1][hyphen + 1 :])
            #                print ("{}: {}"
            #                       .format(inputPar['<nbk>'][0],inputPar['<nbk>'][1]))
            #                print ("{}: {}"
            #                       .format(inputPar['<pag>'][0],inputPar['<pag>'][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<env>":
                inputPar[line[0]] = ("Environment  ", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<zac>":
                inputPar[line[0]] = ("Zero angle   ", float(line[1]))
            #                print ("{}: {} deg"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<wle>":
                inputPar[line[0]] = ("Wavelength   ", float(line[1]))
            #                print ("{}: {} Å"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<sli>":
                inputPar["<lft>"] = ("Left (Lohen) ", float(line[1]))
                inputPar["<rgt>"] = ("Right (Gams) ", float(line[2]))
            #                print ("{}: {} mm"
            #                       .format(inputPar['<lft>'][0],inputPar['<lft>'][1]))
            #                print ("{}: {} mm"
            #                       .format(inputPar['<rgt>'][0],inputPar['<rgt>'][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<fla>":
                inputPar["<top>"] = ("Top          ", float(line[1]))
                inputPar["<bot>"] = ("Bottom       ", float(line[2]))
            #                print ("{}: {} mm"
            #                       .format(inputPar['<top>'][0],inputPar['<top>'][1]))
            #                print ("{}: {} mm"
            #                       .format(inputPar['<bot>'][0],inputPar['<bot>'][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<apr>":
                inputPar[line[0]] = ("Angular prec.", float(line[1]))
            #                print ("{}: {} deg"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<asc>":
                inputPar[line[0]] = (
                    "Angular scale",
                    float(line[1]),
                    float(line[2]),
                    float(line[3]),
                )
            #                print ("{}: from {} deg to {} deg in steps of {} deg"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1],
            #                               inputPar[line[0]][2],inputPar[line[0]][3]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<qsc>":
                inputPar[line[0]] = (
                    "Q       scale",
                    float(line[1]),
                    float(line[2]),
                    float(line[3]),
                )
            #                print ("{}: from {} 1/Å to {} 1/Å in steps of {} 1/Å"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1],
            #                               inputPar[line[0]][2],inputPar[line[0]][3]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<rsc>":
                inputPar[line[0]] = (
                    "r       scale",
                    float(line[1]),
                    float(line[2]),
                    float(line[3]),
                )
            #                print ("{}: from {} Å to {} Å in steps of {} Å"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1],
            #                               inputPar[line[0]][2],inputPar[line[0]][3]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Cmater>":
                inputPar[line[0]] = ("Container material  ", line[1])
                if line[1] in containers:
                    #                    print ("{}: {}"
                    #                           .format(inputPar[line[0]][0],inputPar[line[0]][1]))
                    pass
                else:
                    print(
                        "The container {} is not available".format(inputPar[line[0]][1])
                    )
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Cshape>":
                inputPar[line[0]] = (
                    "Container shape     ",
                    line[1],
                    line[2],
                    line[3],
                    line[4],
                )
                if line[1] in shapes:
                    #                    print ("{}: A {} with {} mm (id) {} mm (od) and {} mm height"
                    #                           .format(inputPar[line[0]][0],inputPar[line[0]][1],
                    #                                   inputPar[line[0]][2],inputPar[line[0]][3],
                    #                                   inputPar[line[0]][4]))
                    pass
                else:
                    print("The shape {} is not available".format(inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Amater>":
                inputPar[line[0]] = ("Absorber material   ", line[1])
                if line[1] in absorbers:
                    #                    print ("{}: {}"
                    #                           .format(inputPar[line[0]][0],inputPar[line[0]][1]))
                    pass
                else:
                    print(
                        "The absorber {} is not available".format(inputPar[line[0]][1])
                    )
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Ashape>":
                inputPar[line[0]] = (
                    "Absorber shape      ",
                    line[1],
                    line[2],
                    line[3],
                    line[4],
                )
                if line[1] in shapes:
                    #                    print ("{}: A {} with {} mm (id) {} mm (od) and {} mm height"
                    #                           .format(inputPar[line[0]][0],inputPar[line[0]][1],
                    #                                   inputPar[line[0]][2],inputPar[line[0]][3],
                    #                                   inputPar[line[0]][4]))
                    pass
                else:
                    print("The shape {} is not available".format(inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Vshape>":
                inputPar[line[0]] = (
                    "Vanadium shape      ",
                    line[1],
                    line[2],
                    line[3],
                    line[4],
                )
                if line[1] in shapes:
                    #                    print ("{}: A {} with {} mm (id) {} mm (od) and {} mm height"
                    #                           .format(inputPar[line[0]][0],inputPar[line[0]][1],
                    #                                   inputPar[line[0]][2],inputPar[line[0]][3],
                    #                                   inputPar[line[0]][4]))
                    pass
                else:
                    print("The shape {} is not available".format(inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Emater>":
                inputPar[line[0]] = ("Environment material", line[1])
                if line[1] in environments:
                    #                    print ("{}: {}"
                    #                           .format(inputPar[line[0]][0],inputPar[line[0]][1]))
                    pass
                else:
                    print(
                        "The environment {} is not available".format(
                            inputPar[line[0]][1]
                        )
                    )
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Eshape>":
                inputPar[line[0]] = (
                    "Environment shape   ",
                    line[1],
                    line[2],
                    line[3],
                    line[4],
                )
                if line[1] in shapes:
                    #                    print ("{}: A {} with {} mm (id) {} mm (od) and {} mm height"
                    #                           .format(inputPar[line[0]][0],inputPar[line[0]][1],
                    #                                   inputPar[line[0]][2],inputPar[line[0]][3],
                    #                                   inputPar[line[0]][4]))
                    pass
                else:
                    print(
                        "The Environment {} is not available".format(
                            inputPar[line[0]][1]
                        )
                    )
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Sdescr>":
                description = lines[i][9:]
                hash = description.find("#")
                inputPar[line[0]] = ("Sample description  ", description[0:hash])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<StempK>":
                inputPar[line[0]] = ("Sample temperature  ", float(line[1]))
            #                print ("{}: {} K"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0][0:5] == "<Sato":
                Nelements += 1
                usefulinfo = lines[i][9:]
                hash = usefulinfo.find("#")
                # remove blank before and after the useful info, and split it
                atomdat = (usefulinfo[0:hash].strip()).split(" ")
                # new list without empty strings
                atompar = [x for x in atomdat if x != ""]
                inputPar[line[0]] = ("Atom " + str(Nelements).zfill(2), *atompar)
                inputPar["<SNelem>"] = ("Sample nbr. elements", Nelements)
            #                print (*inputPar[line[0]])
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Sdegcc>":
                inputPar[line[0]] = ("Sample density      ", float(line[1]))
            #                print ("{}: {} g/cm3"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Smassc>":
                inputPar[line[0]] = ("Sample mass in can  ", float(line[1]))
            #                print ("{}: {} g"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Sshape>":
                inputPar[line[0]] = (
                    "Sample shape        ",
                    line[1],
                    line[2],
                    line[3],
                    line[4],
                )
                if line[1] in shapes:
                    #                    print ("{}: A {} with {} mm (id) {} mm (od) and {} mm height"
                    #                           .format(inputPar[line[0]][0],inputPar[line[0]][1],
                    #                                   inputPar[line[0]][2],inputPar[line[0]][3],
                    #                                   inputPar[line[0]][4]))
                    pass
                else:
                    print("The shape {} is not available".format(inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Sfulln>":
                inputPar[line[0]] = ("Sample fullness     ", float(line[1]))
            #                print ("{}: {} g"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<Stitle>":
                inputPar[line[0]] = ("Sample title        ", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<CoFile>":
                inputPar[line[0]] = ("Container file      ", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<AbFile>":
                inputPar[line[0]] = ("Absorber file       ", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<VaFile>":
                inputPar[line[0]] = ("Vanadium file       ", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<EnFile>":
                inputPar[line[0]] = ("Environment file    ", line[1])
            #                print ("{}: {}"
            #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
            # --------1---------2---------3---------4---------5---------6---------7---------8
            if line[0] == "<SaFile>":
                inputPar[line[0]] = ("Sample file         ", line[1])
        #                print ("{}: {}"
        #                       .format(inputPar[line[0]][0],inputPar[line[0]][1]))
        # --------1---------2---------3---------4---------5---------6---------7---------8
        else:
            print("Input error in line: ", i + 1, " file: ", parafile)
            sys.exit()
    return inputPar


# --------1---------2---------3---------4---------5---------6---------7---------8


class Measurement:
    """
    An object of this class contains all input data for a given sample.

        Input:
        ---------
        -parametersFile (str)
           contains the name of the parameter file

    Date: 24/02/2023
    Author: Gabriel Cuello
    """

    def __init__(self, inputPar):
        self.instrument = inputPar["<ins>"][1]
        self.proposal = inputPar["<pro>"][1]
        self.mainProposer = inputPar["<mpr>"][1]
        self.localContact = inputPar["<lco>"][1]
        self.startDate = inputPar["<sta>"][1]
        self.endDate = inputPar["<end>"][1]
        self.cycle = inputPar["<cyc>"][1]
        self.notebook = inputPar["<nbk>"][1]
        self.page = inputPar["<pag>"][1]
        self.envCode = inputPar["<env>"][1]
        self.zeroAngle = float(inputPar["<zac>"][1])
        self.wavelength = float(inputPar["<wle>"][1])
        self.vslits = (float(inputPar["<lft>"][1]), float(inputPar["<rgt>"][1]))
        self.hslits = (float(inputPar["<top>"][1]), float(inputPar["<bot>"][1]))
        self.beamHeight = self.hslits[0] - self.hslits[1]
        self.beamWidth = self.vslits[0] - self.vslits[1]
        self.angularPrecision = float(inputPar["<apr>"][1])
        self.aScale = (
            float(inputPar["<asc>"][1]),
            float(inputPar["<asc>"][2]),
            float(inputPar["<asc>"][3]),
        )
        self.qScale = (
            float(inputPar["<qsc>"][1]),
            float(inputPar["<qsc>"][2]),
            float(inputPar["<qsc>"][3]),
        )
        self.rScale = (
            float(inputPar["<rsc>"][1]),
            float(inputPar["<rsc>"][2]),
            float(inputPar["<rsc>"][3]),
        )
        self.container = (
            inputPar["<Cmater>"][1],
            inputPar["<Cshape>"][1],
            float(inputPar["<Cshape>"][2]),
            float(inputPar["<Cshape>"][3]),
            float(inputPar["<Cshape>"][4]),
        )
        self.vanadium = (
            "Vanadium",
            inputPar["<Vshape>"][1],
            float(inputPar["<Vshape>"][2]),
            float(inputPar["<Vshape>"][3]),
            float(inputPar["<Vshape>"][4]),
        )
        self.absorber = (
            inputPar["<Amater>"][1],
            inputPar["<Ashape>"][1],
            float(inputPar["<Ashape>"][2]),
            float(inputPar["<Ashape>"][3]),
            float(inputPar["<Ashape>"][4]),
        )
        self.environment = (
            inputPar["<Emater>"][1],
            inputPar["<Eshape>"][1],
            float(inputPar["<Eshape>"][2]),
            float(inputPar["<Eshape>"][3]),
            float(inputPar["<Eshape>"][4]),
        )
        self.sampleDescription = inputPar["<Sdescr>"][1]
        self.sampleTemperature = float(inputPar["<StempK>"][1])
        self.sampleNbrElements = int(inputPar["<SNelem>"][1])
        self.sampleMass = float(inputPar["<Smassc>"][1])
        self.sampleDensity = float(inputPar["<Sdegcc>"][1])
        self.sampleShape = (
            inputPar["<Sshape>"][1],
            float(inputPar["<Sshape>"][2]),
            float(inputPar["<Sshape>"][3]),
            float(inputPar["<Sshape>"][4]),
        )
        self.sampleFullness = float(inputPar["<Sfulln>"][1])
        self.sampleHeight = self.sampleShape[3]
        self.sampleInnerRadius = float(inputPar["<Sshape>"][2]) / 2.0
        self.sampleOuterRadius = float(inputPar["<Sshape>"][3]) / 2.0
        self.sampleHeight = float(inputPar["<Sshape>"][4])

        # This is just a check to avoid mistyping
        if float(self.container[2]) != float(self.sampleShape[2]):
            print(
                "WARNING! Container inner diameter should be equal",
                "to sample outer diameter",
            )
            msg = "    Container id = {} mm != Sample od = {} mm"
            print(msg.format(float(self.container[2]), float(self.sampleShape[2])))
        self.sampleTitle = inputPar["<Stitle>"][1]

        self.sampleVolume = (
            np.pi
            * self.sampleHeight
            / 1000.0
            * (self.sampleOuterRadius - self.sampleInnerRadius) ** 2
        )
        self.sampleEffDensity = self.sampleMass / self.sampleVolume
        self.samplePackingFraction = self.sampleEffDensity / self.sampleDensity

        self.containerFile = inputPar["<CoFile>"][1].strip()
        self.vanadiumFile = inputPar["<VaFile>"][1].strip()
        self.environmentFile = inputPar["<EnFile>"][1].strip()
        self.sampleFile = inputPar["<SaFile>"][1].strip()
        self.absorberFile = inputPar["<AbFile>"][1].strip()

        #         print (self.sampleFile)
        #         print (self.environmentFile)
        #         print (self.containerFile)
        #         print (self.vanadiumFile)
        #         print (self.absorberFile)

        # Reading scattering data from the 3-column ASCII files
        #        print('Opening',self.sampleFile)
        self.sampleData = read_xye(self.sampleFile)
        self.environmentData = read_xye(self.environmentFile)
        self.containerData = read_xye(self.containerFile)
        self.vanadiumData = read_xye(self.vanadiumFile)
        self.absorberData = read_xye(self.absorberFile)
        # self.sampleData      = read_3col(self.sampleFile)
        # self.environmentData = read_3col(self.environmentFile)
        # self.containerData   = read_3col(self.containerFile)
        # self.vanadiumData    = read_3col(self.vanadiumFile)
        # self.absorberData    = read_3col(self.absorberFile)

        atoms = {}
        natoms = 0.0
        for elem in range(self.sampleNbrElements):
            ordatom = str(elem + 1).zfill(2)
            label = "<Sato" + ordatom + ">"
            #            print (elem,ordatom,inputPar[label][0],inputPar[label][1],inputPar[label][2])
            natoms += float(inputPar[label][2])
            atoms[inputPar[label][1]] = float(inputPar[label][2])
        #        print (atoms)
        self.natoms = natoms
        self.atoms = atoms


# --------1---------2---------3---------4---------5---------6---------7---------8
def saveCORRECT(expt):
    parfile = expt.sampleTitle + ".com"
    with open(parfile, "w") as f:
        f.write("! " + parfile + "\n")
        f.write("!\n")
        f.write("instr " + expt.instrument + "\n")
        line = 'sample "{}.adat" {} /temperature={} /density={} /packing={:8.6f} /fullness={}'.format(
            expt.sampleTitle,
            expt.sampleOuterRadius / 10.0,
            expt.sampleTemperature,
            expt.sampleDensity,
            expt.samplePackingFraction,
            expt.sampleFullness,
        )
        f.write(line + "\n")
        # for i in range(expt.sampleNbrElements):
        #     keyatom = '<Sato'+str(i+1).zfill(2)+'>'
        #     f.write('component {}\n'.format(inputPar[keyatom][1]))
        #     print('component {}'.format(inputPar[keyatom][1]))
        can = expt.containerFile[0:-5] + ".adat"
        f.write('container "{}" {}\n'.format(can, expt.container[3] / 20.0))
        bckg = expt.environmentFile[0:-5] + ".adat"
        f.write('background "{}" 0.8\n'.format(bckg))
        f.write('! black "absorber.adat" 0.93\n')
        vana = expt.vanadiumFile[0:-5] + ".adat"
        f.write(
            'vanadium "{}" {} /smoothing=1 /multiplier=1.02\n'.format(
                vana, expt.vanadium[3] / 20.0
            )
        )
        f.write('background /vanadium "{}" 0.85\n'.format(bckg))
        f.write("wavelenght {}\n".format(expt.wavelength))
        f.write("! zeroangle = {} already subtracted\n".format(expt.zeroAngle))
        f.write("zeroangle {}\n".format(0.0))
        f.write("beam {} {}\n".format(expt.beamHeight / 10.0, expt.beamWidth / 10.0))
        f.write("! xout angle\n")
        f.write("! output " + expt.sampleTitle + ".corr\n")
        f.write('! title "' + expt.sampleTitle + '.corr (after correct)"\n')
        f.write("xout q\n")
        f.write("output " + expt.sampleTitle + ".corr.q\n")
        f.write('title "' + expt.sampleTitle + '.corr.q (after correct)"\n')
        f.write("spectrum 1\n")
        f.write("execute/nopause\n")
        f.write("quit\n")
    return


# --------1---------2---------3---------4---------5---------6---------7--------8
def saveFile_xye(filename, x, y, e, heading):
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
    with open(filename, "w") as datafile:
        datafile.write("# " + filename + "\n")
        for i in range(len(heading)):
            datafile.write("# " + heading[i] + "\n")
        for i in range(len(x)):
            datafile.write(
                "{: 9.3f}".format(x[i])
                + " "
                + "{:18.6f}".format(y[i])
                + " "
                + "{:18.6f}".format(e[i])
                + "\n"
            )
        print("File " + filename + " saved")
    return


# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7--------8
def saveFile_3col(filename, data, heading):
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
    x = data[:, 0]
    y = data[:, 1]
    e = data[:, 2]
    with open(filename, "w") as datafile:
        datafile.write("# " + filename + "\n")
        for i in range(len(heading)):
            datafile.write("# " + heading[i] + "\n")
        for i in range(len(x)):
            datafile.write(
                "{: 9.3f}".format(x[i])
                + " "
                + "{:18.6f}".format(y[i])
                + " "
                + "{:18.6f}".format(e[i])
                + "\n"
            )
        print("File " + filename + " saved")
    return


# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7---------
def saveCorrelations(filename, ar, pcr, pdf, rdf, tor, run, heading):
    with open(filename, "w") as datafile:
        datafile.write("# " + filename + "\n")
        for i in range(len(heading)):
            datafile.write("# " + heading[i] + "\n")
        for i in range(len(ar)):
            datafile.write(
                "{: 9.3f}".format(ar[i])
                + " "
                + "{:12.6f}".format(pcr[i])
                + " "
                + "{:12.6f}".format(pdf[i])
                + " "
                + "{:12.6f}".format(rdf[i])
                + " "
                + "{:12.6f}".format(tor[i])
                + " "
                + "{:12.6f}".format(run[i])
                + "\n"
            )
        print("File " + filename + " saved")
    return


# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7---------
def saveRSCF(filename, fou, heading):
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
    with open(filename, "w") as datafile:
        datafile.write("# " + filename + "\n")
        for i in range(len(heading)):
            datafile.write("# " + heading[i] + "\n")
        for i in range(len(fou)):
            datafile.write(
                "{: 9.3f}".format(fou[i][0])
                + " "
                + "{:12.6f}".format(fou[i][1])
                + " "
                + "{:12.6f}".format(fou[i][2])
                + " "
                + "{:12.6f}".format(fou[i][3])
                + " "
                + "{:12.6f}".format(fou[i][4])
                + " "
                + "{:12.6f}".format(fou[i][5])
                + "\n"
            )
        print("File " + filename + " saved")
    return


# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7---------
def read_xye(filename):
    """
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
    """
    data = open(filename, "r")  # Opens the data file in read only mode

    # Creating the lists that will contain abcissas, ordinates and errors
    x = []
    y = []
    e = []

    # Reading the file line by line
    for dataline in data.readlines():
        # dataline is a string that contains each line of the file in data.
        # note that the last character of the string is a 'carriage return', \n
        if "#" not in dataline:  # Only the lines without # are treated.
            # the method .strip(' ') removes blanks at the beginning of the string
            row = dataline.strip(" ")[:-1]
            if len(row) > 0:  # Only the no-empty lines are treated
                columns = row.split()  # This method split the line using the spaces
                x.append(float(columns[0]))
                y.append(float(columns[1]))
                if len(columns) == 3:
                    if (columns[2] == "i") or (columns[2] == "o"):
                        e.append(float(0.0))
                    else:
                        e.append(float(columns[2]))
                else:
                    e.append(float(0.0))
    data.close()
    print(
        "The data file {} read with no errors. Number of data = {}".format(
            filename, len(x)
        )
    )
    #    return np.array(x,y,e)
    #   Converts the lists in arrays
    xa = np.array(x)
    ya = np.array(y)
    ea = np.array(e)
    return xa, ya, ea


# End of read_xye
# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7---------
def read_3col(filename):
    """
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
    """
    data = open(filename, "r")  # Opens the data file in read only mode

    # Creating the lists that will contain abcissas, ordinates and errors
    x = []
    y = []
    e = []

    # Reading the file line by line
    for dataline in data.readlines():
        # dataline is a string that contains each line of the file in data.
        # note that the last character of the string is a 'carriage return', \n
        if "#" not in dataline:  # Only the lines without # are treated.
            # the method .strip(' ') removes blanks at the beginning of the string
            row = dataline.strip(" ")[:-1]
            if len(row) > 0:  # Only the no-empty lines are treated
                columns = row.split()  # This method split the line using the spaces
                x.append(float(columns[0]))
                y.append(float(columns[1]))
                if len(columns) == 3:
                    if (columns[2] == "i") or (columns[2] == "o"):
                        e.append(float(0.0))
                    else:
                        e.append(float(columns[2]))
                else:
                    e.append(float(0.0))
    data.close()
    print(
        "The data file {} read with no errors. Number of data = {}".format(
            filename, len(x)
        )
    )
    #    return np.array(x,y,e)
    #   Converts the lists in arrays
    xa = np.array(x)
    ya = np.array(y)
    ea = np.array(e)
    #   Reshapes the arrays as 2D arrays, but with 1 column
    xa = xa.reshape(xa.shape[0], 1)
    ya = ya.reshape(ya.shape[0], 1)
    ea = ea.reshape(ea.shape[0], 1)
    #   Concatenates the arrays to have a 3-column matrix
    data = np.concatenate((xa, ya), axis=1)
    data = np.concatenate((data, ea), axis=1)
    return data


# End of read_3col
# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7---------
class DataXYE:
    """
    Type: Class

    Object:
        To read data files in ASCII format, two (X Y) or three (X Y E) columns

    Input:
        filename: (string) Filename of the file containing the data.

    Output:
        An instance created with the following attributes and methods.

        self.filename: (string) Filename of the input file
        self.basename: (string) File basename (what it is before extention's dot)
        self.ext: (string) File extension (without the dot)
        self.x: Abscissas
        self.y: Ordinates
        self.e: Errors (or 3rd coordinate). Returns -1 for 2-column files.
        self.head: List of strings with each line in the file header.

        self.xave: Mean value of the abscissas
        self.yave: Mean value of the ordinates
        self.eave: Mean value of the errors

        self.xmin: Minimum value of the abscissas
        self.ymin: Minimum value of the ordinates
        self.emin: Minimum value of the errors

        self.xmax: Maximum value of the abscissas
        self.ymax: Maximum value of the ordinates
        self.emax: Maximum value of the errors

        self.plot(): Makes a simple plot of y coordinate
        self.show(): Shows the data on screen (as a 3-column table)
        self.header(): Prints the header of the file

    Author: Gabriel Cuello
    Created: 29/12/2022
    Modified:
    #--------1---------2---------3---------4---------5---------6---------7---------
    """

    def __init__(self, filename):
        """
        Type: Main function of the Class DataXYE
            The file is read and the attributes are defined here.

        Input:
            filename: (string) The filename of the file containing the data.

        Output:
            The attributes that can be accessed by the instances.
            See the help of this Class for a complete list of attributes.

        Author: Gabriel Cuello
        Created: 29/12/2022
        Modified:
        #--------1---------2---------3---------4---------5---------6---------7---------
        """
        self.filename = filename
        self.basename = os.path.splitext(filename)[0]
        self.ext = os.path.splitext(filename)[1][
            1:
        ]  # Exclude 1st character to avoid the dot
        self.x = []
        self.y = []
        self.e = []
        self.head = []

        data = open(self.filename, "r")
        lines = data.readlines()
        for dataline in lines:
            row = dataline.strip(" ")[:-1]
            if len(row) > 0:  # Only the non empty lines are treated
                if row[0] == "#" or row[0] == "!":
                    self.head.append(row)
                else:
                    columns = row.split()  # This method split the line using the spaces
                    if len(columns) == 2:
                        self.x.append(float(columns[0]))
                        self.y.append(float(columns[1]))
                        self.e.append(-1.0)
                    elif len(columns) == 3:
                        self.x.append(float(columns[0]))
                        self.y.append(float(columns[1]))
                        self.e.append(float(columns[2]))
                    else:
                        print("Wrong file format")
                        sys.exit()
        data.close()
        self.x = np.array(self.x)
        self.y = np.array(self.y)
        self.e = np.array(self.e)
        self.xave = st.mean(self.x)
        self.xmin = min(self.x)
        self.xmax = max(self.x)
        self.yave = st.mean(self.y)
        self.ymin = min(self.y)
        self.ymax = max(self.y)
        self.eave = st.mean(self.e)
        self.emin = min(self.e)
        self.emax = max(self.e)
        self.peaks_x, self.peaks_y = find_peaks_in_range(self.x, self.y, 5.0, 40.0)
        self.xminr, self.yminr = find_minimum_within_range(self.x, self.y, 5.0, 40.0)

    def plot(self, file_format=0, xmin=None, xmax=None, ymin=None, ymax=None):
        """
        Type: Method in DataXYE class

        Object:
            To make a simple plot of the ordinates as function of abscissas.
            To produce a file with the plot.

        Input:
            xmin,xmax: Minimum and maximum values of the x-axis (float, optional)
            ymin,ymax: Minimum and maximum values of the y-axis (float, optional)
            file_format: A string that defines the format (and extension) of the ouput
                         file (string, optional)

        Output:
            A simple plot on the screen.
            A file with the plot in a graphical file.

        Remarks:
          * Several formats are possible for the output file. The kind of file is
            defined by the input parameter file_format, which can must take one
            of the following values: 'png','pdf','jpg','tiff','svg','jpeg','ps','eps'.
            If this paramteter is not present, it takes the default value 0 and no
            output file is created.

          * The output file has the same basename as the input file, but the extension
            corresponding to chosen format.

          * The limits of the axes are optional. Their default value is None, which
            will produce a plot with automatic limits.

        Author: Gabriel Cuello
        Created: 29/12/2022
        Modified:
        #--------1---------2---------3---------4---------5---------6---------7---------
        """
        A_ext = ["adat", "Adat", "reg"]
        Q_ext = ["qdat", "Qdat", "Qreg", "soq", "SoQ"]
        R_ext = ["pcf", "pdf", "tor", "rdf", "run"]

        plt.figure(figsize=(9, 6))

        plt.plot(self.x, self.y, "r-+", label=self.filename)

        plt.legend(loc="best")
        plt.title("Data in " + self.filename)
        plt.xlabel("Abscissa")
        if self.ext in A_ext:
            plt.xlabel(r"$2\theta$ (˚)")
        elif self.ext in Q_ext:
            plt.xlabel(r"$Q$ (Å${-1}$)")
        elif self.ext in R_ext:
            plt.xlabel("$R$ (Å)")
        plt.ylabel("Intensity (arb. units)")
        plt.axis([xmin, xmax, ymin, ymax])
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        if file_format in ["png", "pdf", "jpg", "tiff", "svg", "jpeg", "ps", "eps"]:
            file_fig = "../regdata/" + self.basename + "." + file_format
            plt.savefig(file_fig, format=file_format)
            print("Figure saved on {}".format(file_fig))

    def show(self):
        """
        Type: Method in DataXYE class

        Object:
            To show the data on the screen.

        Input: None

        Output:
            Print out of data on the screen in a 3-column table.

        Author: Gabriel Cuello
        Created: 29/12/2022
        Modified:
        #--------1---------2---------3---------4---------5---------6---------7---------
        """
        A_ext = ["adat", "Adat", "reg"]
        Q_ext = ["qdat", "Qdat", "Qreg", "soq", "SoQ"]
        R_ext = ["pcf", "pdf", "tor", "rdf", "run"]

        if self.ext in A_ext:
            print(f"{'Ang':>6}{'Intensity':>12}{'Error':>12}")
        elif self.ext in Q_ext:
            print(f"{'Q':>6}{'Intensity':>12}{'Error':>12}")
        elif self.ext in R_ext:
            print(f"{'R':>6}{'Intensity':>12}{'Ignore':>12}")
        else:
            print(f"{'x':>6}{'y':>12}{'e':>12}")
        for i in range(len(self.x)):
            print(f"{self.x[i]:>6}{self.y[i]:>12}{self.e[i]:>12}")

    def header(self, lines=1000):
        """
        Type: Method in DataXYE class

        Object:
            To print the file header on the screen.

        Input:
            lines: number of lines to be printed

        Output:
            Print out of file heading on the screen.

        Remarks:
          * The default value for the input parameter lines is 1000. But this method
            will print only the number of lines in the header, unless this number is
            greater than 1000.

        Author: Gabriel Cuello
        Created: 29/12/2022
        Modified:
        #--------1---------2---------3---------4---------5---------6---------7---------
        """
        for i in range(min(lines, len(self.head))):
            print(f"{i+1:>3} {self.head[i]}")


# --------1---------2---------3---------4---------5---------6---------7---------
