#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 16:10:50 2021.

18/01/2022 Modified by José Robledo to include the reading of Nexus files
 He added the function load_nxs()
@author: cuello
"""
import numpy as np
import matplotlib.pyplot as plt  # For plotting
import matplotlib.colors as mcolors
from datetime import datetime
import calendar
import h5py
import sys, os
from .MiscelaneousCalculations import ang2q, rebin

# --------1---------2---------3---------4---------5---------6---------7---------8
def getDate():
    """
    Produce a string with the current date in the format:
        Dayname DD/MM/YYYY HH:MM:SS

    Input: nothing

    Output: string, current_datetime

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # datetime object containing current date and time
    now = datetime.now()
    # dd/mm/YY H:M:S
    #    curr_date = datetime.today()
    #    day = ' '+calendar.day_name[curr_date.weekday()]+' '
    day = " " + calendar.day_name[now.weekday()] + " "
    current_datetime = day + now.strftime("%d/%m/%Y %H:%M:%S")
    return current_datetime


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getDTC_mcounts(monitor, time, dead_time=2.4e-6):
    """
    Calculate the dead-time corrected monitor counts.

    The measured counting rate is
            nu_meas = monitor / time
    then, the dead-time corrected rate is
            nu = nu_meas / (1 - nu_meas * tau)
    where tau is the dead-time in seconds and nu is given in counts/second.

    Thus, the dead-time corrected counts are: nu * time.

    Input: count, time, dead_time
        count : measured counts on the monitor
        time : measuring time (in seconds) (if negative, it changes to 1e-6)
        dead_time : monitor dead-time (in seconds) (default = 2.4e-6)

    Output: Dead-time corrected monitor

    Notes:
        (1) For D4, in 2021 the program d4creg uses tau=2.4E-06 s for the monitor.
        (2) The program d4creg uses the following expresion:
                 nu = (1-sqrt(1-4*nu_meas*tau))/2/tau
        (3) In the Internal ILL Report (ILL98FI15T) by H.E. Fischer, P. Palleau and
            D. Feltin the tau for the monitor is 7.2E-06 s.
        (4) At some point Henry Fischer suggested using tau=12E-6 s for the monitor.

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # If the input counting time is negative, it is forced to be 1e-6 sec, to avoid
    # a division by 0.
    if time <= 0:
        time = 1e-6
    nu_meas = monitor / time  # counting rate
    correction = 1.0 - nu_meas * dead_time
    if correction > 0:
        factor = 1.0 / correction
        nu = nu_meas * factor

        with open("../logfiles/d4creg.log", "a") as file:
            file.write(
                "  Monitor DT = {:4.2e} s. Correction factor = {:10.6f} \n".format(
                    dead_time, factor
                )
            )

    #        print(
    #            "  Monitor DT = {:4.2e} s. Correction factor = {:10.6f}".format(
    #                dead_time, factor))
    else:
        with open("../logfiles/d4creg.log", "a") as file:
            file.write("--- ERROR in Function getDTC_mcounts: correction factor < 0 \n")
    #        print("--- ERROR in Function getDTC_mcounts: correction factor < 0")
    return nu * time


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getDTC_dcounts(counts, total_counts, time, dead_time=7e-06):
    """
    Calculate the dead-time corrected detector counts.

    The measured counting rate is
            nu_meas = counts / time
    then, the dead-time corrected rate is
            nu = nu_meas / (1 - nu_meas * tau)
    where tau is the dead-time in seconds and nu is given in counts/second.

    Thus, the dead-time corrected counts are:  nu * time.

    Input: count, time, dead_time
        count : measured counts on the detector (this a matrix, 9x64 or 9x128)
        time : measuring time (in seconds) (if negative, it changes to 1e-6)
        dead_time : monitor dead-time (in seconds) (default = 7e-6)

    Output: Dead-time corrected detector counts (a matrix, 9x64 or 9x128)

    Notes:
        (1) For D4, in 2021 the program d4creg uses tau=7E-06 for the microstrip
            detector.
        (2) The program d4creg uses the following expresion:
                 nu = (1-sqrt(1-4*nu_m*tau))/2/tau
        (3) In the Internal ILL Report (ILL98FI15T) by H.E. Fischer, P. Palleau and
            D. Feltin the tau for the detector is 7.17E-06 s.
        (4) At some point Henry Fischer suggested using tau=0 s for the detector.

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # If the input counting time is negative, it is forced to be 1e-6 sec, to avoid
    # a division by 0.
    if time <= 0:
        time = 1e-6
    nu_meas = counts / time
    correction = 1.0 - nu_meas * dead_time
    if (correction.all()) > 0:
        factor = 1.0 / correction
        nu = nu_meas * factor
        with open("../logfiles/d4creg.log", "a") as file:
            file.write(
                "  Monitor DT = {:4.2e} s. Correction factor = {:10.6f}".format(
                    dead_time, factor[0][0]
                )
                + " (1st detector cell)\n"
            )
    #        print(
    #            "  Monitor DT = {:4.2e} s. Correction factor = {:10.6f}".format(
    #                dead_time, factor[0][0]), "(1st detector cell")
    else:
        with open("../logfiles/d4creg.log", "a") as file:
            file.write("--- ERROR in Function getDTC_dcounts: correction factor < 0 \n")
    #        print("--- ERROR in Function getDTC_dcounts: correction factor < 0")
    return nu * time, total_counts * factor[0][0]


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def printRates(detector, counts, time, reactor_power=52):
    """
    Print out the basic information about counting rates.

    Input: detector, counts, time, reactor_power
        detector : string, 'detector' or 'monitor'
        counts : float, total number of counts on detectors or on the monitor
        time : float, counting time in seconds (if negative, it changes to 1e-6)
        reactor_power : float, reactor power in MW (default value = 52 MW)

    Output: None

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # If the input counting time is negative, it is forced to be 1e-6 sec, to avoid
    # a division by 0.
    if time <= 0:
        time = 1e-6
    out = "  Counts on {} = {:10.8g}  ===> counting-rate = {:10.8g} c/s".format(
        detector, counts, counts / time
    )
    with open("../logfiles/d4creg.log", "a") as file:
        file.write(out + "\n")

    #    print(out)
    if reactor_power <= 0:
        reactor_power = 52
    out = "  Reactor-power normalised {} counting-rate = {:10.6g} c/s/MW".format(
        detector, counts / time / reactor_power
    )
    with open("../logfiles/d4creg.log", "a") as file:
        file.write(out + "\n")
    #    print(out)
    return


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getRefLines(Lines):
    """
    Determine the reference lines in the ASCII numor files (for D4).

    Input: Lines
        Lines : list of string, containing the all lines in the ASCII file
            Each line is an lement in the string

    Output: lineR, AS[0], AS[1], FS[0], FS[1], FS[2], SS[0]
        lineR: index of the last line with R's
        AS: list with the indexes of lines with A's
        FS: list with the indexes of lines with F's
        SS: list with the indexes of lines with S's

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    a = 0
    f = 0
    s = 0
    FS = []
    AS = []
    SS = []
    for i in range(len(Lines)):
        if Lines[i][0:10] == "RRRRRRRRRR":
            lineR = i
        if Lines[i][0:10] == "FFFFFFFFFF":
            FS.append(i)
            f += 1
        if Lines[i][0:10] == "AAAAAAAAAA":
            AS.append(i)
            a += 1
        if Lines[i][0:10] == "SSSSSSSSSS":
            SS.append(i)
            s += 1
    return lineR, AS[0], AS[1], FS[0], FS[1], FS[2], SS[0]


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def readD4numorASCII(filename):
    """
    Read one numor in ASCII format.

    Create a dictionary with the information in the file header and a 2D-matrix
    with the registered counts.

    Input: filename
        filename: string, the filename of the numor to be read.

    Output: data, counts
        data: A dictionary with the information in the file header.
        counts: A 2D-matrix (ndet x ncell) containing the counts of detection cells.

    Notes:
        (1) This function is already adapted for the electronically-doubled number
            of cells, i.e., 128 cell instead of 64 for each detection bank.
        (2) Inside the function there are 7 variables containing the linenumbers
            where the lines with R, A, F and S are found in the numor file.
        (3) The variable cdel contains the number of columns between two
            consecutive numbers in a line of floats.
        (4) In the case of future format modifications of the ASCII numors, these
            variables should be updated. From 2023 the ASCII format for D4 data will
            be discontinued, so no changes to this format are expected.
        (5) The labels of the variables given in the numor file are used as keys for
            the dictionary containing the header information.

    Author: Gabriel Cuello (ILL)
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # --------1---------2---------3---------4---------5---------6---------7---------8
    # Reading the file and creating the list Lines, for which each element is a
    # line in the file, stored as a string.
    with open(filename, "r") as rawdata:
        Lines = rawdata.readlines()

    # Automatically determines the reference lines in a numor
    lineR, lineAL, lineAI, lineFD, lineFS, lineFM, lineSI = getRefLines(Lines)
    # number of columns between 2 consecutive numbers in a line of floats
    cdel = 16

    # Creating a dictionary which will contain the information about the numor.
    # The keys are defined depending on the information contained in that line.
    # This definition is done in an arbitrary but natural way, here in the code.
    # The keys names could be modified if necessary.
    data = {}
    col = (Lines[lineR + 1].strip(" ")[:-1]).split()
    data["numor"] = str(col[0]).zfill(6)
    data["label"] = Lines[lineAL + 2][:-1]
    data["user"] = Lines[lineAI + 2].strip(" ")[20:-1]
    data["LC"] = Lines[lineAI + 3].strip(" ")[20:-1]
    data["title"] = Lines[lineAI + 4][20:-1]
    data["subtitle"] = Lines[lineAI + 5][20:-1]
    data["proposal"] = Lines[lineAI + 6][20:-1]
    data["startTime"] = Lines[lineAI + 7][20:38]
    data["endTime"] = Lines[lineAI + 7][60:79]

    # Here starts the reading of the floating numbers blocks.
    # Each block starts with a line of 80 'F'.
    # The keys for the dictionary are automatically assigned with the lables found
    # in the numor file. Those labels have been agreed between SCI and IR, and
    # cannot be changed. To change them a new agreement must be found with SCI.
    # Only the SCI team can change labels in a numor file. This is true in
    # general for any change in the format of numors, either ASCII or Nexus.

    # DETECTOR: Reading the detector block from the numor file.
    # This block contains a information about the instrument configuration.
    # Some of the most useful entries are: monitor couts, counting time, 2theta,
    # counts on each detection bank, total number of counts, etc.
    lref = lineFD + 3  # Reference line: the line where the 1st label appears.
    col = (Lines[lineFD + 1].strip(" ")[:-1]).split()
    ldel = int(col[1]) - 1  # Nbr of lines between the labels and data blocks.
    for i in range(ldel):
        for j in range(5):
            parameter = None
            while parameter is None:
                try:
                    parameter = float(
                        (Lines[lref + ldel + i][1 + j * cdel: (j + 1) * cdel]).lstrip(
                            " "
                        )
                    )
                except Exception:
                    parameter = float(0.0)
            data[(Lines[lref + i][1 + j * cdel: (j + 1) * cdel]).lstrip(" ")] = (
                parameter
            )

    # SAMPLE ENVIRONNEMENT AND REACTOR PARAMETERS
    # This block contains information about the sample environment (temperatures,
    # pressures, etc.).
    # There is also information about the reactor power and reactor cycle.
    # In fact, the reactor cycle is always 0 (at least at present, January 2022).
    lref = lineFS + 3  # Reference line: the line where the 1st label appears.
    col = (Lines[lineFS + 1].strip(" ")[:-1]).split()
    ldel = int(col[1]) - 1  # Nbr of lines between the labels and data blocks.
    for i in range(ldel):
        for j in range(5):
            parameter = None
            while parameter is None:
                try:
                    parameter = float(
                        (Lines[lref + ldel + i][1 + j * cdel: (j + 1) * cdel]).lstrip(
                            " "
                        )
                    )
                except Exception:
                    parameter = float(0.0)
            data[(Lines[lref + i][1 + j * cdel: (j + 1) * cdel]).lstrip(" ")] = (
                parameter
            )

    # MONOCHROMATOR
    # This block contains all the parameters from the monochromator.
    # Many parameters are yet properly updated (Jan 2022), so be careful with this.
    # An useful parameter is the incident energy and the d-spacing. These values
    # could be used to determine the default working wavelength.
    lref = lineFM + 3  # Reference line: the line where the 1st label appears.
    col = (Lines[lineFM + 1].strip(" ")[:-1]).split()
    ldel = int(col[1]) - 1  # Nbr of lines between the labels and data blocks.
    for i in range(ldel):
        for j in range(5):
            parameter = None
            while parameter is None:
                try:
                    parameter = float(
                        (Lines[lref + ldel + i][1 + j * cdel: (j + 1) * cdel]).lstrip(
                            " "
                        )
                    )
                except Exception:
                    parameter = float(0.0)
            data[(Lines[lref + i][1 + j * cdel: (j + 1) * cdel]).lstrip(" ")] = (
                parameter
            )

    # COUNTS
    # These are 9 consecutive blocks containing the counts for each detector bank.
    # Each block starts with a line of 80 'S'.
    # The next line contains the current bank, the remaining banks, the total
    # number number of banks (9) and the numor.
    # Then, there is a line with 80 'I', preceding the block of integer numbers.
    # The next line contains a single integer: number of cells in a detector (64).
    # Then, the individual counts for each cell are written (10 values per line).
    # Note that here the script determenines the number of detectors (ndet) and
    # the number of cells per dtector (ncell).
    # At presesent (Jan 2022), there are 9 detectors and 64 cells per detector.
    # This could be changed in a near future to 9 x 128. Thanks to the fact of
    # reading this information frmo the numor, this change should be transparent
    # for this script, provided that we keep the same format for recording the
    # counts, i.e., a single line with the number of cells per detctor (128) and
    # 10 values per line with the counts.
    lref = lineSI  # Reference line: the first line with 80 'S'
    col = (Lines[lineSI + 1].strip(" ")[:-1]).split()
    ndet = int(col[2])  # number of detection banks.
    data["ndet"] = ndet
    ncell = int(Lines[lref + 3].strip(" ")[:-1])
    data["ncell"] = ncell

    # Initialising the matrix for the counts: 9x65 (or 9x129)
    counts = np.zeros((ndet, ncell))

    # A first loop over the 9 detectors
    for det in range(ndet):
        # There are 10 numbers on each line.
        # Then this loop goes up to 6 (or 12), reading 10 numbers at each loop.
        for j in range(int(ncell / 10)):
            # col is a list of 10 elements containing the numbers
            col = (Lines[lref + 4 + j + 11 * det].strip(" ")[:-1]).split()
            # The numbers are assigned to the matrix counts
            for i in range(len(col)):
                counts[det, i + 10 * j] = col[i]
        # Now, the remaining line is read and data stored in the list col.
        # This line contains 4(8) numbers, depending on the nbr of cells, 64(128).
        col = (Lines[lref + 4 + int(ncell / 10) + 11 * det].strip(" ")[:-1]).split()
        # The numbers are assigned to the matrix counts
        for i in range(len(col)):
            counts[det, i + 10 * int(ncell / 10)] = col[i]

    print("Numor {} loaded without errors.".format(data["numor"]))
    print(2 * " ", data["subtitle"])
    print("  Starts: {}  Ends: {}".format(data["startTime"], data["endTime"]))
    print("  Counting time ={:10.6g} s".format(data["CntTime (sec)"]))
    print("  Angle (1st cell 1st det) ={:10.6g} degrees".format(data["Theta_lue(deg)"]))
    printRates(
        "monitor", data["MonitorCnts"], data["CntTime (sec)"], data["RtrPower (MW)"]
    )
    printRates(
        "detector", data["TotalCnts"], data["CntTime (sec)"], data["RtrPower (MW)"]
    )
    print()

    return data, counts


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def load_nxs(nxs_path, ncell=64):
    """
    Read one numor in Nexus format.

    Creates a dictionary with the information in the metadata and a 2D-matrix with
    the registered counts.

    Equivalent to readD4numorASCII() but for nexus files.

    Input: nxs_path, ncell
        nxs_path: string
            The filename and path of the nexus numor to be read (ending with ".nxs")
        ncell: int
            The number of detection cells of the instrument. To be modified to 128
            only when needed.

    Output: metadata, counts
        metadata: A dictionary with the information in the file metadata.
        counts: A 2D-matrix (ndet x ncell) containing the counts of detection cells.

    Notes:
        (1) This function is adapted for the electronically-doubled number of cells,
            i.e., 128 cell instead of 64 for each detection bank.
            In case of using 128 cells, modify ncell to 128.
        (2) This function uses h5py to read the nexus files. The nexus files still
            contains more information than the one loaded here. This function only
            load the information necessary to run all functions in this module.
        (3) The labels of the variables given in the numor file are used as keys for
            the dictionary containing the header information.

    Author: José Robledo (CAB-CNEA)
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # read from nexus
    with h5py.File(nxs_path, "r") as nxs:
        entry0 = nxs["entry0"]

        # create metadata dictionary
        metadata = {
            "label": f"d4 {entry0['user/name'][0].decode('utf-8')}{entry0['user/namelocalcontact'][0].decode('utf-8')}{entry0['start_time'][0].decode('utf-8')}",
            "user": entry0["user/name"][0].decode("utf-8"),
            "LC": entry0["user/namelocalcontact"][0].decode("utf-8"),
            "proposal": entry0["user/proposal"][0].decode("utf-8"),
            "ndet": 9,
            "ncell": ncell,
            "RtrPower (MW)": 1,
            "numor": str(entry0["run_number"][0]),
            "subtitle": entry0["experiment_identifier"][0].decode("utf-8"),
            "title": entry0["title"][0].decode("utf-8"),
            "startTime": entry0["start_time"][0].decode("utf-8"),
            "endTime": entry0["end_time"][0].decode("utf-8"),
            "MonitorCnts": entry0["monitor/data"][0, 0, 0],
            "CntTime (sec)": entry0["time"][0],
            "TotalSteps": entry0["data_scan/total_steps"][0],
            "Theta_lue(deg)": entry0["instrument/2theta/value"][0],
            "Theta_des(deg)": entry0["instrument/2theta/target_value"][0],
            "Omega(deg)": entry0["instrument/omega/value"][0],
            "A1": entry0["instrument/A1/value"][0],
        }

        for name, data in zip(
            entry0["data_scan/scanned_variables/variables_names/property"],
            entry0["data_scan/scanned_variables/data"],
        ):
            if name.decode("utf-8") == "TotalCount":
                name = "TotalCnts"
            else:
                name = name.decode("utf-8")
            metadata[name] = data[0]

        # create counts
        counts = (
            entry0["data_scan/detector_data/data"][0]
            .squeeze()
            .reshape(metadata["ndet"], metadata["ncell"])
        )
    return metadata, counts


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def readParam(parametersFile):
    runInfo = {}

    # --------1---------2---------3---------4---------5---------6---------7---------8
    # This block rarely changes
    #   Angular range covered by a detection bank. At present this value is 8 deg.
    runInfo["angular_range_bank"] = 8.0
    # Dead-time (in s) for monitor and detector.
    #   If these values are equal to 0, no dead-time correction is performed.
    runInfo["dead_time"] = (2.4e-6, 7.0e-6)
    # Tolerance (in cell units) for warning/error in angle positioning
    #   The requested angle is compared to the read angle.
    #   If the difference is bigger than first value, only a warning and continue.
    #   If this difference is bigger than second value, error message and stop.
    runInfo["cellTolerance"] = (0.1, 1.0)
    # Format of the numors: can be 'ASCII' or 'Nexus'
    runInfo["dataFormat"] = "Nexus"
    # --------1---------2---------3---------4---------5---------6---------7---------8

    # --------1---------2---------3---------4---------5---------6---------7---------8
    # This block can change for each experiment.
    # Path (relative to the working directory) to the rawdata
    runInfo["path_raw"] = "rawdata/"
    # Extension for the output data files:
    #   .reg:  Files in D4 format (angle, counts, error, Q)
    #   .adat: Files in CORRECT format (angle, counts, error)
    #   .qdat: Files in CORRECT format (Q, counts, error)
    #   .cdat: Files with single numor (angle, counts, error, Q, cell)
    #   .nxs:  Files in Nexus format, typically the rawdata
    #   .log:  Logfiles
    runInfo["ext"] = (".reg", ".adat", ".qdat", ".cdat", ".nxs", ".log")
    # Efficiency file, containing the relative efficiencies (I=counts/eff).
    #   'ones': replace all efficiencies by 1.
    runInfo["efffile"] = "effd4c.eff"
    # Shift file, containing the small angular correction for each detection bank.
    #   'zeros': replace all shifts by 0.
    #   'manual': shifts given in function getDec (in module readD4.py).
    runInfo["decfile"] = "dec.dec"
    runInfo["logfile"] = "d4creg" + runInfo["ext"][5]
    # --------1---------2---------3---------4---------5---------6---------7---------8

    # --------1---------2---------3---------4---------5---------6---------7---------8
    # This block can change for every experiment, and even for some samples.
    # Defines the normalisation mode and the values used to normalise.
    #   0: 'monitor' or 'time'
    #   1: normalising time in seconds. Normal reactor power, 80s approx 1E6 mon.
    #   2: normalising monitor counts. Usually 1E6 monitor counts.
    runInfo["normalisation"] = ("monitor", 80, 1000000)
    # Write individual files for each numor depending on 'writeNumor' value (0/1).
    #   0: no individual numor files are created.
    #   1: an individual numor file is created for each numor in the range.
    #   If only one numor in the range, 'writeNumor' is forced to 1.
    runInfo["writeNumor"] = 0
    runInfo["plotDiff"] = 0
    # Zero-angle correction (deg) and working wavelength (A).
    runInfo["twotheta0"] = -0.1013
    runInfo["wavelength"] = 0.49857
    # Binning in angular and Q scales.
    runInfo["angular_scale"] = (0.0, 140.0, 0.125)  # Angular scale (degrees)
    runInfo["q_scale"] = (0.0, 23.8, 0.02)  # Q-scale (1/A)
    # --------1---------2---------3---------4---------5---------6---------7---------8

    # --------1---------2---------3---------4---------5---------6---------7---------8
    # This block changes for every sample
    # Name identifying the sample. Root for the output filenames (root.ext).
    runInfo["file"] = "vanadium"
    # Numors to be included. Individual numors or range of numors are accepted.
    # Examples:
    #   - ['387229-387288']
    #   - ['387229','387231','387231','387235-387270','387275','387280-387288']
    #   - ['387229']
    # runInfo['numorLst'] = ['387229-387230']
    runInfo["numorLst"] = ["387229-387288"]

    with open("../logfiles/d4creg.log", "w") as file:
        file.write("# d4creg.log\n\n")

    with open(parametersFile, "r") as parfile:
        lines = parfile.readlines()

    for i in range(len(lines)):
        if len(lines[i]) > 1:
            first = (
                (lines[i][0] != "#") and (lines[i][0] != "!") and (lines[i][0] != "<")
            )
            if first:
                #            print (first)
                print("Wrong input in line: ", i + 1, " file: ", parametersFile)
                sys.exit()

    for i in range(len(lines)):
        #    print (lines[i][0:-1])
        if (lines[i][0] == "#") or (lines[i][0] == "!"):
            #        print ("Comment line",len(lines[i]))
            pass
        elif len(lines[i]) == 1:
            #        print("Blank line")
            pass
        elif (lines[i][0] == "<") and (lines[i][4] == ">"):
            line = lines[i].split(" ")
            if line[0] == "<tol>":
                runInfo["cellTolerance"] = (float(line[1]), float(line[2]))
                print(
                    "Cell tolerance: {} cells for warning, {} cells for error".format(
                        *runInfo["cellTolerance"]
                    )
                )
            if line[0] == "<tau>":
                runInfo["dead_time"] = (float(line[1]), float(line[2]))
                print(
                    "Dead times: {} s for monitor, {} s for detector".format(
                        *runInfo["dead_time"]
                    )
                )
            if line[0] == "<ext>":
                #            runInfo["ext"] = (line[j] for j in range(1,7))
                runInfo["ext"] = (line[1], line[2], line[3], line[4], line[5], line[6])
                print("Extensions: {}, {}, {}, {}, {}, {}".format(*runInfo["ext"]))
            if line[0] == "<eff>":
                runInfo["efffile"] = line[1]
                print("Efficiency file: {}".format(runInfo["efffile"]))
            if line[0] == "<dec>":
                runInfo["decfile"] = line[1]
                print("Shifts file: {}".format(runInfo["decfile"]))
            if line[0] == "<log>":
                runInfo["logfile"] = line[1]
                print("Shifts file: {}".format(runInfo["decfile"]))
            if line[0] == "<fmt>":
                runInfo["dataFormat"] = line[1]
                print("Data format: {}".format(runInfo["dataFormat"]))
            if line[0] == "<asc>":
                runInfo["angular_scale"] = (
                    float(line[1]),
                    float(line[2]),
                    float(line[3]),
                )
                print(
                    "Angular scale: from {} deg to {} deg in steps of {} deg".format(
                        *runInfo["angular_scale"]
                    )
                )
            if line[0] == "<qsc>":
                runInfo["q_scale"] = (float(line[1]), float(line[2]), float(line[3]))
                print(
                    "Q scale: from {} 1/A to {} 1/A in steps of {} 1/A".format(
                        *runInfo["q_scale"]
                    )
                )
            if line[0] == "<wri>":
                if line[1] == "True":
                    runInfo["writeNumor"] = 1
                else:
                    runInfo["writeNumor"] = 0
                print(
                    "Write individual files for each numor: {} {}".format(
                        line[1], runInfo["writeNumor"]
                    )
                )
            if line[0] == "<plo>":
                runInfo["plotDiff"] = 0
                if line[1] == "True":
                    runInfo["plotDiff"] = 1
                else:
                    runInfo["plotDiff"] = 0
                print("Plot diffractograms: {} {}".format(line[1], runInfo["plotDiff"]))
            if line[0] == "<wle>":
                runInfo["wavelength"] = float(line[1])
                print("Incident wavelength: {} A".format(runInfo["wavelength"]))
            if line[0] == "<zac>":
                runInfo["twotheta0"] = float(line[1])
                print("Zero-angle: {} deg".format(runInfo["twotheta0"]))
            if line[0] == "<rdp>":
                runInfo["path_raw"] = line[1]
                print("Relative path for rawdata: {}".format(runInfo["path_raw"]))
            if line[0] == "<nor>":
                runInfo["normalisation"] = (line[1], float(line[2]), float(line[3]))
                print("Normalisation mode: {}".format(runInfo["normalisation"][0]))
                print(
                    "Normalisation constants: {1} s or {2} monitor counts".format(
                        *runInfo["normalisation"]
                    )
                )
            if line[0] == "<out>":
                runInfo["file"] = line[1]
                print("Base name for output files: {}".format(runInfo["file"]))

            if line[0] == "<num>":
                runInfo["numorLst"] = [line[1]]
                print("List of numors: {}".format(runInfo["numorLst"]))

            if line[0] == "<add>":
                runInfo["numorLst"].append(line[1])
                print("List of numors: {}".format(runInfo["numorLst"]))

            elif line[0] == "<run>":
                print("Calling d4creg...")
                d4creg(runInfo)

        else:
            print("Input error in line: ", i + 1, " file: ", parametersFile)
            sys.exit()
    return runInfo


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getEff(filename="effd4c.eff", ndet=9, ncells=64):
    """
    Generate the relative efficiencies.

    Reads the relative efficiencies from the given file (default 'effd4c.eff').
    The reserved name 'ones' fixes all efficiencies to 1.

    Input: filename, ndet, ncells
        filename : string, optional
            The filename for the relative efficiencies. Default: 'effd4c.eff'.
            If the string is equal to 'ones': all the efficiencies are 1.
        ndet : integer, optional
            Number of detector banks. The default is 9.
        ncells : integer, optional
            Number of cells in one detector bank. The default is 64.

    Output:
        efficiencies: a 2D-matrix with the efficiencies
                      axis 0: detector bank
                      axis 1: cell

    Note:
        (1) The efficiency is the coefficient that multiplies the expected counting
            value to give the measured value, i.e., I_expected = I_meas/eff.

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # Initialising the matrix for the efficiencies: 9x64 (or 9x128)
    efficiencies = np.ones((ndet, ncells))  # Initialise the vector with 1's
    if filename == "ones":
        print("No relative efficiencies")
    else:
        with open(filename, "r") as shifts:
            Lines = shifts.readlines()
        for i in range(len(Lines)):
            if "#" not in Lines[i]:  # Only the lines without # are treated.
                row = (Lines[i].strip(" "))[0:-1]
                if len(row) > 0:  # Only the no-empty lines are treated
                    col = row.split()  # This method split the line
                    if float(col[2]) <= 0:
                        col[2] = "nan"
                    efficiencies[int(col[0]) - 1, int(col[1]) - 1] = col[2]
    return efficiencies


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getDec(filename="dec.dec", ndet=9):
    """
    Reads from a file the angular shifts for each detection bank.

    The shifts can be manually introduced within this function using the reserved
    name 'manual'. The reserved name 'zeros' fixes all shifts to 0.

    Input: filename, ndet
      - filename : string, optional
            The filename for the angular shifts. The default is 'dec.dec'.
            If the string is:
                - 'zeros': all the shifts are zero.
                - 'manual': shifts are fixed to the values inside this function.
      - ndet : integer, optional
            Number of detector banks. The default is 9.

    Output: zero
        zero : array with the 9 shift values (in deg).

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    if filename == "zeros":
        print("No angular shifts")
        zero = np.zeros(ndet)  # This initialise the vector with zeros
    elif filename == "manual":
        print("Built-in shifts:")
        zero = np.array(
            [0.000, -0.027, -0.026, -0.087, -0.125, -0.160, -0.220, -0.250, -0.340]
        )
    else:
        zero = np.zeros(ndet)
        with open(filename, "r") as shifts:
            Lines = shifts.readlines()
        for i in range(len(Lines)):
            if "#" not in Lines[i]:  # Only the lines without # are treated.
                row = (Lines[i].strip(" "))[0:-1]
                if len(row) > 0:  # Only the no-empty lines are treated
                    columns = row.split()  # This method split the line
                    zero[int(columns[0]) - 1] = float(columns[1])
    return zero


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getErrorsNumor(counts):
    """
    Calculates the errors corresponding to the detector counts.

    The errors are calculated assuming a Poisson statistics, i.e., as the square
    root of the number of counts.

    Input: counts
      - counts : 2D-matrix of floats (ndet x ncell)
            Contains the registered counts for each cell of each detector

        Returns
        -------
      - errors : 2D-matrix of floats (ndet x ncell)
            Contain the square root of the counts.

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # Initialising the matrix for the errors. Same dimensions as counts matrix.
    errors = np.zeros((counts.shape[0], counts.shape[1]))
    errors = np.sqrt(counts)
    return errors


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getNormalFactors(monitor, time):
    """
    Defines the normalisation factors.

    Definition of a 2x2 matrix containing the normalisation factors i.e., monitor
    counts or counting time
        - 1st row: monitor counts, error of monitor counts
        - 2nd row: counting time, error in counting time (in seconds)

    The row index is the normalisation method. Thus,
    normal[norm][0], normal[norm][1] = value, error
    where norm = 0 --> normalisation by monitor
          norm = 1 --> normalisation by counting time

    Input: monitor, time
        monitor : float
            monitor counts
        time : float
            counting time (in seconds)

    Output: normal
      - normal: normalisation factors and errors (2x2 matrix)

    Notes:
      - The error in counting time is supposed to be 0.01 seconds.

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    normal = np.zeros((2, 2))
    normal[0][0] = monitor  # Monitor counts
    normal[0][1] = np.sqrt(monitor)  # Error monitor counts
    normal[1][0] = time  # Counting time
    normal[1][1] = 0.01  # Error counting time (s)
    return normal


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def normalise(
    counts, errors, normal, norm_type="monitor", norm_time=120, norm_mon=1000000
):
    """
    Normalise the detector counts by monitor counts or counting time.

    Normalises the counts and errors by monitor or counting time.
        - monitor:  counts/monitor * norm_mon
        - time:  counts/time * norm_time

    Input: counts,errors,normal,norm_type=,norm_time,norm_mon
      - counts : 2D matrix (ndet x ncell)
            Matrix containing the counts of a single numor
      - errors : 2D matrix (ndet x ncell)
            Matrix containing the errors of a single numor
      - normal : 2D matrix (2 x 2)
            Matrix containing the normalisation factors and errors
      - norm_type : string, optional
            Defines the normalusaion mode. The default is 'monitor'.
            'monitor': normalisation by monitor
            'time': normalisation by counting time
      - norm_time : float, optional
            Counting time for normalisation (in seconds). The default is 120 s.
      - norm_mon : float, optional
            Monitor counts for normalisation. The default is 1000000.

    Output:
      - counts : 2D matrix (ndet x ncell)
            Matrix containing the normalised counts of a single numor
      - errors : 2D matrix (ndet x ncell)
            Matrix containing the normalised errors of a single numor

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    norm = [norm_mon, norm_time]
    ntype = 0  # normalisation by monitor
    if norm_type == "time":
        ntype = 1  # normalisation by time
    counts = norm[ntype] * counts / normal[ntype][0]
    relative_error = np.sqrt(
        (errors / counts) ** 2 + (normal[ntype][1] / normal[ntype][0]) ** 2
    )
    errors = counts * relative_error
    return counts, errors


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getAngles(counts, ttheta_ref, zeros, cell_step=0.125, det_step=15.0):
    """
    Determine the angles and cell labels fro a gove numor.

    For a given numor, i.e., a reference scattering angle, this function
    returns the list of angles for each cell of the 9 detection blocks.
    An identification of each cell is produced too, as:
        det * 1000 + cell
    in this way 7050 is the 50th cell of detector 7.

    Input:
      - counts : 2D matrix (ndet x ncell)
            Counts for a given numor.
            Used only to determine the total number of cells.
      - ttheta_ref : float
            Reference angle for the 1st cell of the 1st detector (in deg)
      - zeros : array of floats
            The angular shifts of each detection block (in deg)
      - cell_step : float, optional
            The angular range covered by 1 cell, in degrees. Default: 0.125 deg.
      - det_step : float, optional
            Angular distance between the 1st cells of 2 detection blocks, in deg.
            The default is 15.0 degrees.

    Output: angles, cells
      - angles : array of angles, float
            Angle of each cell
      - cells : array of identifiers, float
            Indentifier for each cell

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    angles = np.zeros((counts.shape[0], counts.shape[1]))
    cells = np.zeros((counts.shape[0], counts.shape[1]))
    ndet = counts.shape[0]
    ncel = counts.shape[1]
    for det in range(ndet):
        for cel in range(ncel):
            angles[det][cel] = (
                ttheta_ref + zeros[det] + (det * det_step + cel * cell_step)
            )
            cells[det][cel] = (det + 1) * 1000 + cel + 1
    return angles, cells


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def saveOneNumor(x, y, e, q, c, head, runInfo):
    """
    Save an individual file for a single numor.

    Create a file with the numor name and extension .cdat.
    The file contains a header with basic inofrmation, then five columns with
    angle (deg), counts, error, Q (1/Å), and cell.

    Input: x, y, e, q, c, head, runInfo
      - x : float
            Angle in degrees
      - y : float
            Counts
      - e : float
            Errors
      - q : float
            Momentum transfer (1/A)
      - c : integer
            Cell identifier
      - head : dictionary
            Contains the basic information about the numor.
      - runInfo : dictionaty
            Contains the basic information about the running parameters.

    Output: None
        A file with extemsion .cdat is created.

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    filename = head["numor"] + runInfo["ext"][3]
    with open(filename, "w") as datafile:
        datafile.write("# {} \n".format(filename))
        datafile.write("# Sample: {}\n".format(head["subtitle"]))
        datafile.write(
            "# Starts: {}  Ends: {} \n".format(head["startTime"], head["endTime"])
        )
        datafile.write("# Counting time: {:8.4f} s\n".format(head["CntTime (sec)"]))
        datafile.write(
            "# Monitor: {:10.2f} counts ({:8.2f} c/s)\n".format(
                head["MonitorCnts"], head["MonitorCnts"] / head["CntTime (sec)"]
            )
        )
        datafile.write(
            "# Detector: {:10.2f} counts ({:8.2f} c/s)\n".format(
                head["TotalCnts"], head["TotalCnts"] / head["CntTime (sec)"]
            )
        )
        datafile.write("#" + 80 * "-" + "\n")
        datafile.write(
            "# Angle(deg)"
            + 6 * " "
            + "Counts"
            + 14 * " "
            + "Error"
            + 8 * " "
            + "Q(1/Å)"
            + 6 * " "
            + "Cell \n"
        )
        for i in range(len(x)):
            datafile.write(
                "{: 9.3f} {:18.6f} {:18.6f} {:9.3f} {:9.0f}".format(
                    x[i], y[i], e[i], q[i], c[i]
                )
                + "\n"
            )
        print("File {} saved.".format(filename))

    return


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getOneNumor(numorfile, runInfo):
    """
    Reads a numor file and create the counting lists.

    This function puts together all the information coming from a single numor.
    But also:
        - Check the angular positioning, comparing the requested angle to the
          read one. If the difference is greater than a given number of cells
          (1 cell), produces an error message and stops the program.
          A second threshold (0.1 cell) produces only a warning on the screen.
          These two limits are controlled via runInfo dictionary.
        - Make the deadtime correcton on monitor and detector counts.
        - Divide the registered counts by the relative efficiency of cells.
        - Correct the angles by the angular shifts of each detection bank.

    Required functions:
        readD4numorASCII: Read one numor in ASCII format
        getDTC_mcounts:   Correct monitor counts by deat-time
        getDTC_dcounts:   Correct detector counts by deat-time
        printRates:       Print the counting rates on the screen
        getEff:           Read the relative efficiencies
        getErrorsNumor:   Calculate the errors in detector counts
        getNormalFactors: Determine the normalisation factors
        normalise:        Perform the data normalisation
        getDec:           Read the angular shift of each detection bank
        getAngles:        Determine the angles and cell ID for each cell
        ang2q:            Convert scattering angles in Q
        saveOneNumor:     Save on disk a single numor

    Input: numorfile, runInfo
      - numorfile : string
            The name identifying the numor.
      - runInfo : dictionary
            Contains basic information about the running parameters.

    Output: x, y, e, q, c, head
      - x : float
            Angle in degrees
      - y : float
            Counts
      - e : float
            Errors
      - q : float
            Momentum transfer (1/A)
      - c : integer
            Cell identifier
      - head : dictionary
            Contains the basic information about the numor.

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    write = runInfo["writeNumor"]  # write individual numors: 0 or 1
    efffile = runInfo["efffile"]  # efficiency file
    decfile = runInfo["decfile"]  # shifting file
    norm_type = runInfo["normalisation"][0]  # 'monitor' or 'time'
    norm_time = runInfo["normalisation"][1]  # time for normalisation
    norm_mon = runInfo["normalisation"][2]  # monitor counts for normalisation
    dtm = runInfo["dead_time"][0]  # monitor dead-time
    dtd = runInfo["dead_time"][1]  # detector dead-time
    cell_tol = runInfo["cellTolerance"]  # cell tolerance: warning/error

    #    print(80 * "-")
    if runInfo["dataFormat"] == "Nexus":
        # Reading the numor in Nexus format
        head, counts = load_nxs(numorfile + runInfo["ext"][4])
    else:
        # Reading the numor in ASCII format
        head, counts = readD4numorASCII(numorfile)

    # Checking the angular positioning
    # For 64 cells, 0.125 deg
    angular_range_cell = runInfo["angular_range_bank"] / head["ncell"]
    deltaTheta = np.abs(head["Theta_lue(deg)"] - head["Theta_des(deg)"])
    if deltaTheta > cell_tol[1] * angular_range_cell:
        with open("../logfiles/d4creg.log", "a") as file:
            file.write("   ERROR! The angle differs from the requested angle by\n")
            file.write("          more than {} cell\n".format(cell_tol[1]))
            file.write(
                "          Angle = {} deg, requested angle = {} deg\n".format(
                    head["Theta_lue(deg)"], head["Theta_des(deg)"]
                )
            )
        print("   ERROR! The angle differs from the requested angle by")
        print("          more than {} cell".format(cell_tol[1]))
        print(
            "          Angle = {} deg, requested angle = {} deg".format(
                head["Theta_lue(deg)"], head["Theta_des(deg)"]
            )
        )
        print(0 / 0)  # This terminates the program
    elif deltaTheta > cell_tol[0] * angular_range_cell:
        with open("../logfiles/d4creg.log", "a") as file:
            file.write("   WARNING! The angle differs from the requested angle by\n")
            file.write("            more than {} cell\n".format(cell_tol[0]))
            file.write(
                "            Angle = {} deg, requested angle = {} deg\n".format(
                    head["Theta_lue(deg)"], head["Theta_des(deg)"]
                )
            )
        print("   WARNING! The angle differs from the requested angle by ")
        print("            more than {} cell".format(cell_tol[0]))
        print(
            "            Angle = {} deg, requested angle = {} deg".format(
                head["Theta_lue(deg)"], head["Theta_des(deg)"]
            )
        )

    # Dead-time corrections
    with open("../logfiles/d4creg.log", "a") as file:
        file.write("Dead-time corrections\n")
    #    print("Dead-time corrections")
    head["MonitorCnts"] = getDTC_mcounts(
        head["MonitorCnts"], head["CntTime (sec)"], dead_time=dtm
    )
    counts, head["TotalCnts"] = getDTC_dcounts(
        counts, head["TotalCnts"], head["CntTime (sec)"], dead_time=dtd
    )

    printRates(
        "monitor", head["MonitorCnts"], head["CntTime (sec)"], head["RtrPower (MW)"]
    )
    printRates(
        "detector", head["TotalCnts"], head["CntTime (sec)"], head["RtrPower (MW)"]
    )

    # Relative efficiencies correction
    efficiency = getEff(filename=efffile, ndet=head["ndet"], ncells=head["ncell"])
    counts = counts / efficiency  # raw data are divided by the efficiencies
    # Calculation of the experimental errors
    errors = getErrorsNumor(counts)
    # Normalising the data
    normal = getNormalFactors(head["MonitorCnts"], head["CntTime (sec)"])
    counts, errors = normalise(
        counts,
        errors,
        normal,
        norm_type=norm_type,
        norm_time=norm_time,
        norm_mon=norm_mon,
    )
    # Calculating the angular coordinates
    zeros = getDec(filename=decfile, ndet=head["ndet"])
    #   2theta = 2theta.raw + zeros[det]
    angles, cells = getAngles(counts, head["Theta_lue(deg)"], zeros)

    with open("../logfiles/d4creg.log", "a") as file:
        file.write(80 * "-" + "\n")
    #    print(80 * "-")

    x, y, e, c = [], [], [], []
    for det in range(head["ndet"]):
        for cell in range(len(counts[0, :])):
            x.append(angles[det][cell])
            y.append(counts[det][cell])
            e.append(errors[det][cell])
            c.append(cells[det][cell])

    ang = np.array(x)
    q = ang2q(ang, wlength=runInfo["wavelength"])
    if write == 1:
        saveOneNumor(x, y, e, q, c, head, runInfo)

    return x, y, e, q, c, head


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getNumorFiles(numorLst):
    """
    Creates a detailed list of numors.

    Take a short syntax list of numors and returns a detailed list of numors
    (strings that are used as filenames for readung numors).

    Input: numorList
      - numorLst : List of numors in short syntax (strings)
            Each element can be a single numor (6 characters) or a range of numors
            in the format 111111-222222, where 111111 and 222222 are the first and
            last numors in the range. The last numor must be bigger or equal to the
            first numor.

    Output: numorFiles
      - numorFiles : Detailed list of numors

    author: Gabriel Cuello
    date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    numorFiles = []
    lenNumor = 6  # Default length of a numor name, i.e., 6 characters
    for i in range(len(numorLst)):
        if len(numorLst[i]) == lenNumor:
            numorFiles.append(numorLst[i])
        elif len(numorLst[i]) == 2 * lenNumor + 1:  # this corresponde to a range
            firstNum = int(numorLst[i][0:lenNumor])
            lastNum = int(numorLst[i][lenNumor + 1: 2 * lenNumor + 1])
            if firstNum > lastNum:
                with open("../logfiles/d4creg.log", "a") as file:
                    file.write("ERROR with the range \n", numorLst[i])
                    file.write("      The last numor must be less than or equal\n")
                    file.write("      to the first one.\n")
                print("ERROR with the range ", numorLst[i])
                print("      The last numor must be less than or equal")
                print("      to the first one.")
            for j in range(firstNum, lastNum + 1):
                numorFiles.append(str(j).zfill(lenNumor))
        else:
            with open("../logfiles/d4creg.log", "a") as file:
                file.write("ERROR in the list of numors (getNumorFiles function)\n")
                file.write("      It is likely a problem with a numor not having six\n")
                file.write("      characters.\n")
            print("ERROR in the list of numors (getNumorFiles function)")
            print("      It is likely a problem with a numor not having six")
            print("      characters.")
    #     with open('../logfiles/d4creg.log', 'a') as file:
    #         file.write("List of numors included in the diffractogram ({} in total):\n".format(
    #             len(numorFiles)))
    #         file.write(numorFiles)
    #         file.write("\n")

    #    print(
    #        "List of numors included in the diffractogram ({} in total):".format(
    #            len(numorFiles)))
    #    print(numorFiles)
    #    print()
    return numorFiles


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getNumors(runInfo):
    """
    Puts all numors together in the same matrix.

    This function creates a 3D matrix with all the experimental data.

    Required functions:
        getNumorFiles: Create the detailed list of numors in the range
        getOneNumor: Read a numor file and create the counting lists

    Input: runInfo
        runInfo : dictionary
            Contains basic information about the running parameters.

    Output: numor, head
        numor : 3D matrix (nbrNumors x 5 x ndata)
            The 3D matrix numor contains all the data:
            - The first dimension contains the order of the numor, starting from 0
            - The second dimension contains the variable indexes:
                  0: the angle (deg)
                  1: the counts
                  2: the errors
                  3: the Q values (1/A)
                  4: the cell (20xx corresponds to the cell xx of detector 2)
            - The third dimension contains the values of the corresponding variable
        head : dictionary
            Contains the basic information about the numor.

    author: Gabriel Cuello
    date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # --------1---------2---------3---------4---------5---------6---------7---------8
    numorFiles = getNumorFiles(runInfo["numorLst"])
    if len(numorFiles) == 1:
        runInfo["writeNumor"] = 1

    head = []

    # Reading the first numor
    numorfile = runInfo["path_raw"] + numorFiles[0]
    #    print("Numor 1/{}".format(len(numorFiles)),' ',numorFiles[0])

    with open("../logfiles/d4creg.log", "a") as file:
        file.write(runInfo["file"] + "\n")
        file.write("Numor 1/{}".format(len(numorFiles)) + " " + numorFiles[0] + "\n")

    angle, count, error, qval, cell, header = getOneNumor(numorfile, runInfo)
    #    print(header["ndet"])
    angle = np.array(angle)
    head.append(header)

    numor = np.zeros((len(numorFiles), 5, len(angle)))

    numor[0, 0, :] = angle - runInfo["twotheta0"]
    numor[0, 1, :] = count
    numor[0, 2, :] = error
    numor[0, 3, :] = qval
    numor[0, 4, :] = cell

    # Reading the other numors, if more than 1 (note that the loop starts at 1)
    for i in range(1, len(numorFiles)):
        numorfile = runInfo["path_raw"] + numorFiles[i]
        #        print("Numor {}/{}".format(i + 1, len(numorFiles)),' ',numorFiles[i])

        with open("../logfiles/d4creg.log", "a") as file:
            file.write(
                "Numor {}/{}".format(i + 1, len(numorFiles))
                + " "
                + numorFiles[i]
                + "\n"
            )

        angle, count, error, qval, cell, header = getOneNumor(numorfile, runInfo)
        angle = np.array(angle)
        head.append(header)

        numor[i, 0, :] = angle - runInfo["twotheta0"]
        numor[i, 1, :] = count
        numor[i, 2, :] = error
        numor[i, 3, :] = qval
        numor[i, 4, :] = cell

    return numor, head


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def getDiffA(numor, head, runInfo):
    """
    Creates the diffractograms in angular scale.

    Produces 10 diffractograms as function of the scattering angle (in deg).
    The diffractograms 1 to 9 correspond to the detector of the same index.
    The diffractogram 0 corresponds to a single diffractogram regrouping all
    detectors.

    Required functions:
      - ang2q: Convert scattering angles in Q

    Input: numor, head, runInfo
      - numor : 3D matrix (nbrNumors x 5 x nbrAngles)
            The 3D matrix numor contains all the data:
            - The first dimension contains the order of the numor, starting from 0
            - The second dimension contains the variable indexes:
                  0: the angle (deg)
                  1: the counts
                  2: the errors
                  3: the Q values (1/A)
                  4: the cell (20xx corresponds to the cell xx of detector 2)
            - The third dimension contains the values of the corresponding variable
      - head : dictionary
            Contains the header of each numor
      - runInfo : dictionary
            Contains basic information about the running parameters.

    Output: diff
       - diff: a 3D matrix (10 x 5 x NbrAngles)
            - The first dimension is the detector (0 for the total diffractogram)
            - The second dimension contains the variables:
                0: the angle (deg)
                1: the counts
                2: the errors
                3: the weighting factor for each bin
                4: the Q values (1/A)
            - The third dimension contains the values of the corresponding variable

    Author: Gabriel Cuello
    Date:   Jan 2022
    --------------------------------------------------------------------------------
    """
    # --------1---------2---------3---------4---------5---------6---------7---------8
    ndet = head[0]["ndet"]  # Number of detectors
    ncell = head[0]["ncell"]  # Number of cells in one detector
    angmin = runInfo["angular_scale"][0]
    angmax = runInfo["angular_scale"][1]
    angstep = runInfo["angular_scale"][2]
    angular_range = runInfo["angular_range_bank"]
    angular_step = angular_range / ncell  # Angular step for 1 cell (in deg)

    if angstep < angular_step:
        print("WARNING: In getDiffA function")
        print("         The angular step should not be smaller than")
        print("         {:6.3f} degrees".format(angular_step))
        print(
            "         Then, it has been changed to {:6.3f} degrees".format(angular_step)
        )
        angstep = angular_step

    # Number of angles in the angular binning
    nangles = int((angmax - angmin) / angstep)

    # Defining the matrix that will contain the 10 diffractograms
    diff = np.zeros((10, 5, nangles))

    # The angular scale is the same for all diffractograms
    for i in range(ndet + 1):
        diff[i, 0, :] = np.arange(angmin, angmax, angstep)

    # The Q-scale is the same for all diffractograms
    for i in range(ndet + 1):
        diff[i, 4, :] = ang2q(diff[i, 0, :], wlength=runInfo["wavelength"])

    # Number of numors to regroup and number of angles in one numor
    nbrNumors = len(numor[:, 0, 0])
    nbrData = len(numor[0, 0, :])

    # This variable will be used for normalisation
    # It is divided by 1000000 just to avoid large numbers
    totalmon = 0
    for num in range(nbrNumors):
        totalmon += head[num]["MonitorCnts"] / 1000000

    for num in range(nbrNumors):
        # Extract the values from the matrix numor using clearer names
        angle = numor[num, 0, :]
        count = numor[num, 1, :]
        error = numor[num, 2, :]
        cell = numor[num, 4, :]
        # Calculates the normalisation factor for the weighted sum
        mon = head[num]["MonitorCnts"] / 1000000 / totalmon
        for ang in range(nbrData):
            # Left side of the rectangle
            angle1 = angle[ang] - angular_step / 2.0
            ang1 = (angle1 - angmin) / angstep
            a1 = int(ang1)
            # Right side of the rectangle
            angle2 = angle[ang] + angular_step / 2.0
            ang2 = (angle2 - angmin) / angstep
            a2 = int(ang2)
            # Fraction that falls on each bin: left=1, right=2
            frac1 = 1.0 - (ang1 - int(ang1))
            frac2 = ang1 - int(ang1)
            # Checks if counts and errors are numbers
            # If NaN, it is not added to the histogram
            if (np.isnan(count[ang]) is False) and (np.isnan(error[ang]) is False):
                # The index 0 is used for the total diffractogram
                diff[0, 1, a1] += frac1 * count[ang] * mon
                diff[0, 1, a2] += frac2 * count[ang] * mon
                diff[0, 2, a1] += frac1 * error[ang] * mon
                diff[0, 2, a2] += frac2 * error[ang] * mon
                diff[0, 3, a1] += frac1 * mon
                diff[0, 3, a2] += frac2 * mon
                # Creates a histogram for individual detectors
                # det is an integer from 1 to 9
                det = int(cell[ang] / 1000)
                diff[det, 1, a1] += frac1 * count[ang] * mon
                diff[det, 1, a2] += frac2 * count[ang] * mon
                diff[det, 2, a1] += frac1 * error[ang] * mon
                diff[det, 2, a2] += frac2 * error[ang] * mon
                diff[det, 3, a1] += frac1 * mon
                diff[det, 3, a2] += frac2 * mon

    # Normalisation by the weighting factor of each bin
    for det in range(10):
        for i in range(len(diff[0, 1, :])):
            if (np.isnan(diff[0, 1, i]) is False) and (diff[det, 3, i] > 0):
                diff[det, 1, i] = diff[det, 1, i] / diff[det, 3, i]
                diff[det, 2, i] = diff[det, 2, i] / diff[det, 3, i]
            else:
                diff[det, 1, i] = "NaN"
                diff[det, 2, i] = "NaN"

    return diff


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def saveDiffAngle(diffA, head, runInfo):
    """
    --------------------------------------------------------------------------------
    """
    # --------1---------2---------3---------4---------5---------6---------7---------8
    # Loop to count the number of data
    ndata = 0
    for i in range(len(diffA[0, 0, :])):
        angle = diffA[0, 0, i]
        count = diffA[0, 1, i]
        error = diffA[0, 2, i]
        if (np.isnan(count * error) is False) and (count > 0):
            ndata += 1

    # Writing a file with CORRECT format
    file = runInfo["file"] + runInfo["ext"][1]
    with open(file, "w") as adatfile:
        adatfile.write("# " + file + "\n")
        adatfile.write("#Block  1\n")
        adatfile.write("#========\n")
        adatfile.write("#\n")
        adatfile.write("#Instrument:" + head[0]["label"][0:2] + "\n")
        adatfile.write(
            "#User      :"
            + head[0]["user"].strip(" ")
            + "    exp_"
            + head[0]["proposal"].strip(" ")
            + "/processed"
            + "\n"
        )
        adatfile.write("#Run number:            1\n")
        adatfile.write("#Spectrum  :            1\n")
        adatfile.write("#Title     :" + file + "\n")

        adatfile.write("#Run date  :" + runInfo["runDate"] + "\n")

        mini = runInfo["angular_scale"][0]
        maxi = runInfo["angular_scale"][1]
        step = runInfo["angular_scale"][2]
        bins = int((maxi - mini) / step) + 1
        adatfile.write(
            "#X caption : 2Theta (degrees) binned from {} to {} by {} ({} bins)".format(
                mini, maxi, step, bins
            )
            + "\n"
        )
        adatfile.write("#Y caption : Counts/monitor \n")
        adatfile.write("#Histogram :            F \n")
        adatfile.write("#Points    :         {}   ".format(ndata) + "\n")

        for i in range(len(diffA[0, 0, :])):
            angle = diffA[0, 0, i]
            count = diffA[0, 1, i]
            error = diffA[0, 2, i]
            if (np.isnan(count * error) is False) and (count > 0):
                adatfile.write(
                    "{:8.4f}     {:12.8f}     {:12.8f} \n".format(angle, count, error)
                )
    print(
        4 * " "
        + "File "
        + runInfo["file"]
        + runInfo["ext"][1]
        + " (CORRECT format, in angle scale)"
    )

    listNum = " "
    for i in range(len(runInfo["numorLst"])):
        listNum += runInfo["numorLst"][i] + " "

    # Writing the data as 9 individual detectors (format reg from D4)
    file = runInfo["file"] + runInfo["ext"][0]
    with open(file, "w") as regfile:
        regfile.write("# " + file + "\n")
        regfile.write("# Equivalent command line: \n")
        regfile.write("#     " + runInfo["d4creg_cl"] + "\n")
        regfile.write("# Efficiency file: " + runInfo["efffile"] + "\n")
        regfile.write("# Shift file: " + runInfo["decfile"] + "\n")
        regfile.write("# Sample: " + head[0]["subtitle"] + "\n")
        regfile.write("# User: " + head[0]["user"].strip(" ") + "\n")
        regfile.write("# Local contact: " + head[0]["LC"].strip(" ") + "\n")
        regfile.write("# Proposal: " + head[0]["proposal"].strip(" ") + "\n")
        regfile.write("# Run date  :" + runInfo["runDate"] + "\n")
        regfile.write(
            "# Requested numors: {}, {} numors included \n".format(
                listNum, runInfo["nbrNumors"]
            )
        )
        regfile.write("# Normalisation type: " + runInfo["normalisation"][0] + "\n")
        regfile.write("# Zero-angle = " + str(runInfo["twotheta0"]) + " deg\n")
        regfile.write("# Wavelength = " + str(runInfo["wavelength"]) + " deg\n")
        regfile.write("# Monitor dead-time = " + str(runInfo["dead_time"][0]) + " s\n")
        regfile.write("# Detector dead-time = " + str(runInfo["dead_time"][1]) + " s\n")
        regfile.write("# 2theta = 2theta.raw - Zero-angle + shifts \n")
        regfile.write(
            "# Binned on angle (from {} to {} in steps of {}, deg), but not on Q \n".format(
                *runInfo["angular_scale"]
            )
        )
        regfile.write(
            "#     Q = 4pi/{}Å * sin(angle/2) \n".format(str(runInfo["wavelength"]))
        )

        for det in range(1, 10):
            regfile.write("# ----------- \n")
            regfile.write("# Detector " + str(det) + "\n")
            regfile.write("# ----------- \n")
            regfile.write("# Angle(deg)     Counts           Error            Q(1/Å)\n")
            for i in range(len(diffA[0, 0, :])):
                angle = diffA[det, 0, i]
                count = diffA[det, 1, i]
                error = diffA[det, 2, i]
                qval = diffA[det, 4, i]
                if (np.isnan(count * error) is False) and (count > 0):
                    regfile.write(
                        "{:8.4f}     {:12.8f}     {:12.8f}     {:8.4f} \n".format(
                            angle, count, error, qval
                        )
                    )

    print(
        "    File "
        + runInfo["file"]
        + runInfo["ext"][0]
        + " (D4 format, in angle and Q-scale)"
    )
    return


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def saveDiffQ(diffQ, head, runInfo):
    """
    --------------------------------------------------------------------------------
    """
    # --------1---------2---------3---------4---------5---------6---------7---------8
    # Loop to count the number of data
    ndata = 0
    for i in range(len(diffQ)):
        #        q = diffQ[i,0]
        count = diffQ[i, 1]
        error = diffQ[i, 2]
        if (np.isnan(count * error) is False) and (count > 0):
            ndata += 1

    # Writing a file with CORRECT format
    file = runInfo["file"] + runInfo["ext"][2]
    with open(file, "w") as qdatfile:
        qdatfile.write("# " + file + "\n")
        qdatfile.write("#Block  1\n")
        qdatfile.write("#========\n")
        qdatfile.write("#\n")
        qdatfile.write("#Instrument:" + head[0]["label"][0:2] + "\n")
        qdatfile.write(
            "#User      :"
            + head[0]["user"].strip(" ")
            + "    exp_"
            + head[0]["proposal"].strip(" ")
            + "/processed"
            + "\n"
        )
        qdatfile.write("#Run number:            1\n")
        qdatfile.write("#Spectrum  :            1\n")
        qdatfile.write("#Title     :" + file + "\n")
        qdatfile.write("#Run date  :" + runInfo["runDate"] + "\n")
        mini = runInfo["q_scale"][0]
        maxi = runInfo["q_scale"][1]
        step = runInfo["q_scale"][2]
        bins = int((maxi - mini) / step) + 1
        qdatfile.write(
            "#X caption : Q (1/Å) binned from {:4.3f} to {:4.3f} by {:4.3f} ({} bins)".format(
                mini, maxi, step, bins
            )
            + "\n"
        )
        qdatfile.write("#Y caption : Counts/monitor\n")
        qdatfile.write("#Histogram :            F\n")
        qdatfile.write("#Points    :         {}   ".format(ndata) + "\n")

        for i in range(len(diffQ)):
            q = diffQ[i, 0]
            count = diffQ[i, 1]
            error = diffQ[i, 2]
            if (np.isnan(count * error) is False) and (count > 0):
                qdatfile.write(
                    "{:8.4f}     {:12.8f}     {:12.8f} \n".format(q, count, error)
                )
    print(
        "    File "
        + runInfo["file"]
        + runInfo["ext"][2]
        + " (CORRECT format, in Q scale)"
    )
    return


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def defColors():
    """
    --------------------------------------------------------------------------------
    """
    # names = list(mcolors.BASE_COLORS)
    # print (names)
    # names = list(mcolors.TABLEAU_COLORS)
    # print (names)
    # names = list(mcolors.CSS4_COLORS)
    # print (names)

    # --------1---------2---------3---------4---------5---------6---------7---------8
    colors = []
    colors.append(mcolors.CSS4_COLORS["red"])  #  0
    colors.append(mcolors.CSS4_COLORS["blue"])  #  1
    colors.append(mcolors.CSS4_COLORS["green"])  #  2
    colors.append(mcolors.CSS4_COLORS["black"])  #  3
    colors.append(mcolors.CSS4_COLORS["cyan"])  #  4
    colors.append(mcolors.CSS4_COLORS["magenta"])  #  5
    colors.append(mcolors.CSS4_COLORS["gold"])  #  6
    colors.append(mcolors.CSS4_COLORS["orange"])  #  7
    colors.append(mcolors.CSS4_COLORS["brown"])  #  8
    colors.append(mcolors.CSS4_COLORS["gray"])  #  9
    colors.append(mcolors.CSS4_COLORS["silver"])  # 10
    colors.append(mcolors.CSS4_COLORS["pink"])  # 11
    colors.append(mcolors.CSS4_COLORS["purple"])  # 12
    colors.append(mcolors.CSS4_COLORS["navy"])  # 13
    colors.append(mcolors.CSS4_COLORS["teal"])  # 14
    colors.append(mcolors.CSS4_COLORS["lime"])  # 15
    colors.append(mcolors.CSS4_COLORS["olive"])  # 16
    colors.append(mcolors.CSS4_COLORS["salmon"])  # 17
    colors.append(mcolors.CSS4_COLORS["maroon"])  # 18
    colors.append(mcolors.CSS4_COLORS["chartreuse"])  # 19
    return colors


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def makePlotsA(runInfo, diffA, head, numor):
    """
    --------------------------------------------------------------------------------
    """
    colors = defColors()
    # Plot the individual numors
    plt.figure(figsize=(9, 5))
    plt.title(
        "Individual numors ({} in total) for sample {}".format(
            runInfo["nbrNumors"], runInfo["file"]
        )
    )

    amin, amax = runInfo["angular_scale"][0], runInfo["angular_scale"][1]
    ymin, ymax = 0.9 * np.nanmin(diffA[0, 1, :]), 1.1 * np.nanmax(diffA[0, 1, :])
    plt.axis([amin, amax, ymin, ymax])
    plt.xlabel("Scattering angle (degrees)")
    if runInfo["normalisation"][0] == "monitor":
        plt.ylabel("Counts/(" + str(runInfo["normalisation"][2]) + " monitor counts)")
    else:
        plt.ylabel("Counts/(" + str(runInfo["normalisation"][1]) + " seconds)")
    for i in range(numor.shape[0]):
        plt.plot(numor[i, 0, :], numor[i, 1, :], colors[i % 20], label=head[i]["numor"])
    plt.grid(True)
    plt.show()

    # Plot the individual detectors
    plt.figure(figsize=(9, 6))
    plt.title("Individual detectors for " + runInfo["file"])
    amin, amax = runInfo["angular_scale"][0], runInfo["angular_scale"][1]
    ymin, ymax = 0.9 * np.nanmin(diffA[0, 1, :]), 1.1 * np.nanmax(diffA[0, 1, :])
    plt.axis([amin, amax, ymin, ymax])
    plt.xlabel("Scattering angle (degrees)")
    if runInfo["normalisation"][0] == "monitor":
        plt.ylabel("Counts/(" + str(runInfo["normalisation"][2]) + " monitor counts)")
    else:
        plt.ylabel("Counts/(" + str(runInfo["normalisation"][1]) + " seconds)")
    for i in range(1, 10):
        plt.plot(diffA[i, 0, :], diffA[i, 1, :], colors[i], label="Detector " + str(i))
    plt.legend(loc="best")
    plt.grid(True)
    plt.show()

    # Plot the final diffractogram
    plt.figure(figsize=(9, 6))
    plt.title("Diffractogram for " + runInfo["file"])
    amin, amax = runInfo["angular_scale"][0], runInfo["angular_scale"][1]
    ymin, ymax = 0.9 * np.nanmin(diffA[0, 1, :]), 1.1 * np.nanmax(diffA[0, 1, :])
    plt.axis([amin, amax, ymin, ymax])
    plt.xlabel("Scattering angle (degrees)")
    if runInfo["normalisation"][0] == "monitor":
        plt.ylabel("Counts/(" + str(runInfo["normalisation"][2]) + " monitor counts)")
    else:
        plt.ylabel("Counts/(" + str(runInfo["normalisation"][1]) + " seconds)")
    #    linestyle = 'solid','dotted','dashed','dashdot'
    #    marker = 'o', 'v', '^','+','x'
    plt.plot(
        diffA[0, 0, :],
        diffA[0, 1, :],
        color=colors[0],
        linestyle="solid",
        label=runInfo["file"],
    )
    plt.legend(loc="best")
    plt.grid(True)
    plt.show()
    return


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def makePlotsQ(runInfo, diffQbin, head, numor):
    """
    --------------------------------------------------------------------------------
    """
    colors = defColors()

    # Plot the final diffractogram
    plt.figure(figsize=(9, 6))
    plt.title("Diffractogram for " + runInfo["file"])
    qmin, qmax = runInfo["q_scale"][0], runInfo["q_scale"][1]
    ymin, ymax = 0.9 * np.nanmin(diffQbin[:, 1]), 1.1 * np.nanmax(diffQbin[:, 1])
    plt.axis([qmin, qmax, ymin, ymax])
    plt.xlabel("Momentum transfer (1/Å)")
    if runInfo["normalisation"][0] == "monitor":
        plt.ylabel("Counts/(" + str(runInfo["normalisation"][2]) + " monitor counts)")
    else:
        plt.ylabel("Counts/(" + str(runInfo["normalisation"][1]) + " seconds)")
    #    linestyle = 'solid','dotted','dashed','dashdot'
    #    marker = 'o', 'v', '^','+','x'
    plt.plot(
        diffQbin[:, 0],
        diffQbin[:, 1],
        color=colors[0],
        linestyle="solid",
        label=runInfo["file"],
    )
    plt.legend(loc="best")
    plt.grid(True)
    plt.show()
    return


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def d4creg(runInfo):
    """
    --------------------------------------------------------------------------------
    """
    runInfo["runDate"] = getDate()  # Today's date

    numor, head = getNumors(runInfo)

    runInfo["nbrNumors"] = len(head)

    # Creating and printing the equivalent d4creg command line
    d4creg_cl = "d4creg"
    d4creg_cl += " -x " + str(runInfo["angular_scale"][0]) + " "
    d4creg_cl += str(runInfo["angular_scale"][1]) + " "
    d4creg_cl += str(runInfo["angular_scale"][2])
    d4creg_cl += " -z " + str(runInfo["twotheta0"])
    d4creg_cl += " -w " + str(runInfo["wavelength"])
    d4creg_cl += " -o " + runInfo["file"] + ".reg "
    d4creg_cl += head[0]["numor"] + " " + head[-1]["numor"]
    runInfo["d4creg_cl"] = d4creg_cl

    print()
    print("Equivalent d4creg command line")
    print(d4creg_cl)
    print()
    print(" - Efficiency: ", runInfo["efffile"])
    print(" - Det shifts: ", runInfo["decfile"])
    if runInfo["normalisation"][0] == "monitor":
        print(
            " - Normalisation to {1} of {0} counts".format(
                runInfo["normalisation"][0], runInfo["normalisation"][2]
            )
        )
    else:
        print(
            " - Normalisation to {1} s of counting {0}".format(
                runInfo["normalisation"][0], runInfo["normalisation"][1]
            )
        )
    print(" - Number of numors: ", runInfo["nbrNumors"])

    # Creating the diffactogram in angular scale
    diffA = getDiffA(numor, head, runInfo)

    # Saving the output files in angular scale
    print()
    print("Output:")
    saveDiffAngle(diffA, head, runInfo)

    if runInfo["plotDiff"] == 1:
        makePlotsA(runInfo, diffA, head, numor)

    qmin = runInfo["q_scale"][0]
    qmax = runInfo["q_scale"][1]
    qstep = runInfo["q_scale"][2]

    dataQ = np.zeros((diffA.shape[2], 3))

    ang = diffA[0, 0, :]
    dataQ[:, 0] = ang2q(ang, wlength=runInfo["wavelength"])

    for i in range(diffA.shape[2]):
        for j in range(1, 3):
            dataQ[i, j] = diffA[0, j, i]

    diffQbin = rebin(0.125, runInfo["wavelength"], dataQ, qmin, qmax, qstep)

    for i in range(len(diffQbin)):
        if diffQbin[i, 1] <= 0.0:
            for j in range(3):
                diffQbin[i, j] = "NaN"
                diffQbin[i - 1, j] = "NaN"

    saveDiffQ(diffQbin, head, runInfo)

    if runInfo["plotDiff"] == 1:
        makePlotsQ(runInfo, diffQbin, head, numor)

    # Here change the name of the logfile

    # Specify the current path of the file you want to rename
    current_file_path = "../logfiles/d4creg.log"

    # Specify the new name and path for the file
    new_file_name = "../logfiles/" + runInfo["logfile"] + ".log"

    # Rename the file
    try:
        os.rename(current_file_path, new_file_name)
        print(f"File '{current_file_path}' renamed to '{new_file_name}' successfully.")
    except FileNotFoundError:
        print(f"File '{current_file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

    return


# --------1---------2---------3---------4---------5---------6---------7---------8
