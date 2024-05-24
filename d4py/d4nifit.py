#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 17:17:33 2022

@author: cuello
"""

import sys
import os
import numpy as np  # NumPy library
import matplotlib.pyplot as plt  # For plotting
import lmfit as lm  # For fitting
import statistics as st
from scipy.signal import find_peaks
from .readD4 import getDate
from .FittingModels import fittingRange


def find_minimum_within_range(x, y, x_min, x_max):
    # Mask the data within the specified x range
    mask = (x >= x_min) & (x <= x_max)
    x_range = x[mask]
    y_range = y[mask]

    # Find the index of the minimum y value within the specified x range
    min_index = np.argmin(y_range)

    # Get the x position and minimum y value
    min_x_position = x_range[min_index]
    min_y_value = y_range[min_index]

    return min_x_position, min_y_value


def find_peaks_in_range(x, y, x_min, x_max):
    # Mask the data within the specified x range
    mask = (x >= x_min) & (x <= x_max)
    x_range = x[mask]
    y_range = y[mask]

    # Find peaks within the specified x range
    peaks, _ = find_peaks(y_range)
    _, ymin = find_minimum_within_range(x, y, x_min, x_max)

    # Get the x positions and maximum values of the peaks within the range
    peak_x_positions = x_range[peaks]
    peak_maximum_values = y_range[peaks] / ymin

    mask = peak_maximum_values > 2.0
    peak_x = peak_x_positions[mask]
    peak_y = peak_maximum_values[mask]

    return peak_x, peak_y


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


def setting_model_d4nifit(nifile, zac_ini=0.0):

    # Reading the data file and creating an instance of the class DataXYE
    nickel = DataXYE(nifile)

    # Plot the not yet normalised data
    nickel.plot()

    # Because the structure of Ni is fcc (with a lattice parametr of a = 3.52024 Å), the
    # first peak corresponds to the reflection (111):
    #         2 * a * sin (2theta/2) = sqrt(h**2+k**2+l**2) * lambda
    # where 2theta is the scattering angle of the observed first peak.
    # Thus, a nominal wavelength (and initial value for the fitting process) can easily
    # be estimated:
    #     lambda = 2 * a * sin(2theta/2) / sqrt(3)
    #
    # The instance has a property called peaks_x, which is a list of the peak positions
    # found in the range 5-40 degrees. Thus, the angular position of the first peak is
    # peaks_x[0]

    lattice = 3.52024
    wl_ini = 2 * lattice / np.sqrt(3.0) * np.sin(nickel.peaks_x[0] * np.pi / 360.0)
    print(f"Automatically calculated initial wavelength = {wl_ini:.3f} Å")

    # There is a flat valley between reflections (200) and (220). For normalisation
    # purposes (with the idea of having "standard" plots), the intensity in that flat
    # valley is used to normalise the data. The middle point between the 2nd peak (200)
    # and 3rd peak (220) is chosen for that normalisation. The positions of these peaks
    # are peaks_x[1] and peaks_x[2], respectively.

    flatValley_x = int((nickel.peaks_x[1] + nickel.peaks_x[2]) / 2)

    # Find the index where x equals the flatValley_x
    index = np.where(nickel.x == flatValley_x)

    # Check if the flatValley_x exists in x
    if index[0].size > 0:
        # Get the corresponding value of y using the index
        flatValley_y = nickel.y[index][0]
        print(
            f"Normalised data using the intensity at {flatValley_x:.0f} degrees: {flatValley_y:.2f}"
        )
    else:
        print(f"2theta = {flatValley_x:.0f} is not in the array of angles.")

    # Normalise the y-scale
    nickel.y = nickel.y / flatValley_y
    nickel.ymax = nickel.ymax / flatValley_y

    # Plot the already normalised data
    nickel.plot()

    # Using the nominal wavelength, some wavelength dependent parameters are initialised:
    # axes:   these are the limits of the plots showing the detailed region of fitting
    # limits: these are the limits used in the fitting
    # Note that by chance 93 times the wavelength is the good upper limit in angle.
    axes = [0.5 * flatValley_x, 93 * wl_ini, -0.1 * nickel.ymax, 1.1 * nickel.ymax]
    limits = [0.5 * flatValley_x, 93 * wl_ini, 0, 1.1 * nickel.ymax]

    # Here the coordinates of text blocks in plots are defined.
    # res: block with the results (wavelength and zero)
    # res: block with the list of reflections
    # fit: block with the information about the fit (iterations and other)
    xy_res = [0.81 * flatValley_x, 0.90 * max(nickel.y)]
    xy_ref = [1.60 * flatValley_x, 0.40 * max(nickel.y)]
    xy_fit = [1.85 * flatValley_x, 0.75 * max(nickel.y)]

    # FPAarea is the area of the first peak. Approximating this area by
    # A = 2 * sigma * I_max, where sigma is about 0.25, then A = 0.5 * I_max
    FPArea = 0.5 * (nickel.ymax - 1)
    #    print ('Area =', FPArea)

    # Setting the initial values for all parameters:
    # - the 3 first parameters are for the polynomial
    # - then wl and zac
    # - then the 10 areas, the 10 asymmetric factors and the 10 sigmas (Gaussian functions)
    # Note that the initial centroids are not setted because they are calculated from the
    # structure of the fcc lattice, for the initial wavelength and zero angle correction, and
    # the known lattice constant for Ni (a = 3.52024 Å).
    initial = [
        1,
        0.0007 * FPArea,
        0.0,
        wl_ini,
        zac_ini,
        1.00 * FPArea,
        0.56 * FPArea,
        0.54 * FPArea,
        0.76 * FPArea,
        0.24 * FPArea,
        0.12 * FPArea,
        0.42 * FPArea,
        0.40 * FPArea,
        0.32 * FPArea,
        0.37 * FPArea,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.2,
        0.2,
        0.2,
        0.2,
        0.2,
        0.2,
        0.3,
        0.3,
        0.4,
        0.4,
    ]

    planes = [
        "(111)",
        "(200)",
        "(220)",
        "(311)",
        "(222)",
        "(400)",
        "(331)",
        "(420)",
        "(422)",
        "(333)",
    ]

    # Calculates the first ten reflections
    refl = d4.reflections_fcc(wl_ini, zac_ini, lattice=3.52024)[:10]
    model = d4.niPeaks10(nickel.x, *initial)
    # initial_text = "lambda = {:.6f} Å \n2theta0 = {:.6f} deg".format(wl_ini,zac_ini)

    # Plot the figure with the starting point for the fitting
    plt.figure(figsize=(9, 6))

    plt.plot(nickel.x, model, "r-", label="Initial model")
    plt.plot(nickel.x, nickel.y, "b-+", label="Ni powder data")
    plt.plot(nickel.x, nickel.y - model, "g-+", label="Residuals")

    plt.legend(loc="best")
    plt.title("Nickel powder diffractogram: Starting point")
    plt.xlabel(r"$2\theta$ (˚)")
    plt.ylabel("Intensity (arb. units)")

    wlength = r"$\lambda = {:.3f}$ Å".format(wl_ini)
    zac = r"$2\theta_0 = {:.3f}$˚".format(zac_ini)

    plt.text(
        xy_res[0],
        xy_res[1],
        s=wlength + "\n" + zac,
        bbox=dict(facecolor="white", alpha=1.0, edgecolor="black"),
    )

    planes_fit = ""
    for i in range(10):
        planes_fit += "{}: {:.1f}˚\n".format(planes[i], refl[i])
    plt.text(
        xy_ref[0],
        xy_ref[1],
        s=planes_fit[:-1],
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="black"),
    )

    plt.axis(axes)
    plt.grid(True)
    plt.tight_layout()
    return axes, limits, xy_res, xy_ref, xy_fit, planes, initial, nickel


def fitting_model_d4nifit(nickel, limits, initial):
    #         """
    # Type: function

    # Object:
    #     Defines the fitting model for the fit of Ni diffractogram for wavelength
    #     calibration nad perform the fitting

    # Input:
    #     nickel: An instance of DataXYE class with the data
    #     limits: List with for values corresponding to the fitting box
    #     initial: Initial values of the fitting parameters

    # Output:
    #     An object containing the results of the fitting process.

    # Remarks:
    #   * In thsi function you can choose the parameters to be fitted.

    # Author: Gabriel Cuello
    # Created: 29/12/2022
    # Modified:
    # #--------1---------2---------3---------4---------5---------6---------7---------
    #         """
    xfit = nickel.x
    yfit = nickel.y
    efit = nickel.e

    xfit, yfit, efit = fittingRange(
        limits[0], limits[1], limits[2], limits[3], nickel.x, nickel.y, nickel.e
    )

    gmodel = lm.Model(d4.niPeaks10)
    gmodel.set_param_hint("I0", vary=True)
    gmodel.set_param_hint("slope", vary=True)
    gmodel.set_param_hint("quad", vary=False)
    gmodel.set_param_hint("wavelength", vary=True)
    gmodel.set_param_hint("twotheta0", vary=True)
    gmodel.set_param_hint("A0", vary=True)
    gmodel.set_param_hint("A1", vary=True)
    gmodel.set_param_hint("A2", vary=True)
    gmodel.set_param_hint("A3", vary=True)
    gmodel.set_param_hint("A4", vary=True)
    gmodel.set_param_hint("A5", vary=True)
    gmodel.set_param_hint("A6", vary=True)
    gmodel.set_param_hint("A7", vary=True)
    gmodel.set_param_hint("A8", vary=True)
    gmodel.set_param_hint("A9", vary=True)
    gmodel.set_param_hint("G0", vary=False)
    gmodel.set_param_hint("G1", vary=False)
    gmodel.set_param_hint("G2", vary=False)
    gmodel.set_param_hint("G3", vary=False)
    gmodel.set_param_hint("G4", vary=False)
    gmodel.set_param_hint("G5", vary=False)
    gmodel.set_param_hint("G6", vary=False)
    gmodel.set_param_hint("G7", vary=False)
    gmodel.set_param_hint("G8", vary=False)
    gmodel.set_param_hint("G9", vary=False)
    gmodel.set_param_hint("S0", vary=True)
    gmodel.set_param_hint("S1", vary=True)
    gmodel.set_param_hint("S2", vary=True)
    gmodel.set_param_hint("S3", vary=True)
    gmodel.set_param_hint("S4", vary=True)
    gmodel.set_param_hint("S5", vary=True)
    gmodel.set_param_hint("S6", vary=True)
    gmodel.set_param_hint("S7", vary=True)
    gmodel.set_param_hint("S8", vary=True)
    gmodel.set_param_hint("S9", vary=True)

    # Here the fit is performed
    result = gmodel.fit(
        yfit,
        x=xfit,
        I0=initial[0],
        slope=initial[1],
        quad=initial[2],
        wavelength=initial[3],
        twotheta0=initial[4],
        A0=initial[5],
        A1=initial[6],
        A2=initial[6],
        A3=initial[7],
        A4=initial[8],
        A5=initial[10],
        A6=initial[11],
        A7=initial[12],
        A8=initial[13],
        A9=initial[14],
        G0=initial[15],
        G1=initial[16],
        G2=initial[17],
        G3=initial[18],
        G4=initial[19],
        G5=initial[20],
        G6=initial[21],
        G7=initial[22],
        G8=initial[23],
        G9=initial[24],
        S0=initial[25],
        S1=initial[26],
        S2=initial[27],
        S3=initial[28],
        S4=initial[29],
        S5=initial[30],
        S6=initial[31],
        S7=initial[32],
        S8=initial[33],
        S9=initial[34],
    )
    return result


def showing_results_d4nifit(
    result, nickel, axes, limits, xy_res, xy_ref, xy_fit, planes
):
    print()
    print("Results of fitting the nickel powder sample")

    nickel_table = {}
    nickel_par = []

    for param in result.params.values():
        if param.value != 0:
            relative = abs(param.stderr / param.value)
        else:
            relative = 0.0
        nickel_table[param.name] = [param.value, param.stderr, relative]
        nickel_par.append(param.value)
        if param.name == "wavelength":
            wlength = r"$\lambda = ({:.6f} \pm {:.6f})$ Å".format(
                param.value, param.stderr
            )
        if param.name == "twotheta0":
            zac = r"$2\theta_0 = ({:.6f} \pm {:.6f})$˚".format(
                param.value, param.stderr
            )

    model = d4.niPeaks10(nickel.x, *nickel_par)
    refl = d4.reflections_fcc(nickel_par[3], nickel_par[4], lattice=3.52024)[:10]

    print(
        f"{result.params['I0'].name:>10} = \
          {nickel_table['I0'][0]:>9f} ± {nickel_table['I0'][1]:>9f}         {nickel_table['I0'][2]:>7.3%}"
    )

    print(
        f"{result.params['slope'].name:>10} = \
          {nickel_table['slope'][0]:>9f} ± {nickel_table['slope'][1]:>9f}         {nickel_table['slope'][2]:>7.3%}"
    )

    print(
        f"{result.params['quad'].name:>10} = \
          {nickel_table['quad'][0]:>9f} ± {nickel_table['quad'][1]:>9f}         {nickel_table['quad'][2]:>7.3%}"
    )

    print(
        f"{result.params['wavelength'].name:>10} = \
          ({nickel_table['wavelength'][0]:>8f} ± {nickel_table['wavelength'][1]:>9f}) Å      {nickel_table['wavelength'][2]:>7.3%}"
    )

    print(
        f"{result.params['twotheta0'].name:>10} = \
          ({nickel_table['twotheta0'][0]:>8.5f} ± {nickel_table['twotheta0'][1]:>9f})˚       {nickel_table['twotheta0'][2]:>7.3%}"
    )

    print(103 * "-")
    print(
        f"|{'Plane':>6s} | {'Center':>7s}  | {'Area        ':>20s}  |  {'%   ':>8s} ||\
           {'Sigma        ':>19s}  |   {'%   ':>8s} |"
    )
    print(103 * "-")
    print(
        f"|{planes[0]:>6s} | {refl[0]:>8.5f} | {nickel_table['A0'][0]:>8.5f} ± {nickel_table['A0'][1]:>9f}  |  {nickel_table['A0'][2]:>8.3%} ||\
          {nickel_table['S0'][0]:>8.5f} ± {nickel_table['S0'][1]:>9f}  |  {nickel_table['S0'][2]:>9.3%} |"
    )

    print(
        f"|{planes[1]:>6s} | {refl[1]:>8.5f} | {nickel_table['A1'][0]:>8.5f} ± {nickel_table['A1'][1]:>9f}  |  {nickel_table['A1'][2]:>8.3%} ||\
          {nickel_table['S1'][0]:>8.5f} ± {nickel_table['S1'][1]:>9f}  |  {nickel_table['S1'][2]:>9.3%} |"
    )
    print(
        f"|{planes[2]:>6s} | {refl[2]:>8.5f} | {nickel_table['A2'][0]:>8.5f} ± {nickel_table['A2'][1]:>9f}  |  {nickel_table['A2'][2]:>8.3%} ||\
          {nickel_table['S2'][0]:>8.5f} ± {nickel_table['S2'][1]:>9f}  |  {nickel_table['S2'][2]:>9.3%} |"
    )
    print(
        f"|{planes[3]:>6s} | {refl[3]:>8.5f} | {nickel_table['A3'][0]:>8.5f} ± {nickel_table['A3'][1]:>9f}  |  {nickel_table['A3'][2]:>8.3%} ||\
          {nickel_table['S3'][0]:>8.5f} ± {nickel_table['S3'][1]:>9f}  |  {nickel_table['S3'][2]:>9.3%} |"
    )
    print(
        f"|{planes[4]:>6s} | {refl[4]:>8.5f} | {nickel_table['A4'][0]:>8.5f} ± {nickel_table['A4'][1]:>9f}  |  {nickel_table['A4'][2]:>8.3%} ||\
          {nickel_table['S4'][0]:>8.5f} ± {nickel_table['S4'][1]:>9f}  |  {nickel_table['S4'][2]:>9.3%} |"
    )
    print(
        f"|{planes[5]:>6s} | {refl[5]:>8.5f} | {nickel_table['A5'][0]:>8.5f} ± {nickel_table['A5'][1]:>9f}  |  {nickel_table['A5'][2]:>8.3%} ||\
          {nickel_table['S5'][0]:>8.5f} ± {nickel_table['S5'][1]:>9f}  |  {nickel_table['S5'][2]:>9.3%} |"
    )
    print(
        f"|{planes[6]:>6s} | {refl[6]:>8.5f} | {nickel_table['A6'][0]:>8.5f} ± {nickel_table['A6'][1]:>9f}  |  {nickel_table['A6'][2]:>8.3%} ||\
          {nickel_table['S6'][0]:>8.5f} ± {nickel_table['S6'][1]:>9f}  |  {nickel_table['S6'][2]:>9.3%} |"
    )
    print(
        f"|{planes[7]:>6s} | {refl[7]:>8.5f} | {nickel_table['A7'][0]:>8.5f} ± {nickel_table['A7'][1]:>9f}  |  {nickel_table['A7'][2]:>8.3%} ||\
          {nickel_table['S7'][0]:>8.5f} ± {nickel_table['S7'][1]:>9f}  |  {nickel_table['S7'][2]:>9.3%} |"
    )
    print(
        f"|{planes[8]:>6s} | {refl[8]:>8.5f} | {nickel_table['A8'][0]:>8.5f} ± {nickel_table['A8'][1]:>9f}  |  {nickel_table['A8'][2]:>8.3%} ||\
          {nickel_table['S8'][0]:>8.5f} ± {nickel_table['S8'][1]:>9f}  |  {nickel_table['S8'][2]:>9.3%} |"
    )
    print(
        f"|{planes[9]:>6s} | {refl[9]:>8.5f} | {nickel_table['A9'][0]:>8.5f} ± {nickel_table['A9'][1]:>9f}  |  {nickel_table['A9'][2]:>8.3%} ||\
          {nickel_table['S9'][0]:>8.5f} ± {nickel_table['S9'][1]:>9f}  |  {nickel_table['S9'][2]:>9.3%} |"
    )
    print(103 * "-")

    plt.figure(figsize=(9, 6))
    plt.plot(nickel.x, model, "r-", label="Fit")
    plt.plot(nickel.x, nickel.y, "b+", label="Ni powder data")
    plt.plot(nickel.x, nickel.y - model, "g-+", label="Residuals")

    plt.legend(loc="best")
    plt.title("Final fit for " + nickel.filename + ", " + getDate())
    plt.xlabel(r"$2\theta$" + " (˚)")
    plt.ylabel("Intensity (arb. units)")

    ite = result.fit_report().find("function evals")
    chi = result.fit_report().find("chi-square")
    red = result.fit_report().find("reduced chi-square")
    aka = result.fit_report().find("Akaike")
    ite = result.fit_report()[ite + 18 : chi - 62]
    chi = result.fit_report()[chi + 20 : red - 5]
    red = result.fit_report()[red + 20 : aka - 5]
    info_fit = "Iterations ={}\nchi-sq ={}\nreduced chi-sq ={}".format(ite, chi, red)

    signature = "d4nifit.ipynb (Dec 2022), G.J. Cuello (ILL), cuello@ill.fr"

    plt.text(
        xy_fit[0],
        xy_fit[1],
        s=info_fit,
        bbox=dict(facecolor="white", alpha=1.0, edgecolor="black"),
    )
    plt.text(
        xy_res[0],
        xy_res[1],
        s=wlength + "\n" + zac,
        bbox=dict(facecolor="white", alpha=1.0, edgecolor="black"),
    )
    plt.text(
        axes[0] - 0.05 * (axes[1] - axes[0]),
        axes[2] - 0.1 * (axes[3] - axes[2]),
        s=signature,
        fontsize=6,
    )

    planes_fit = ""
    for i in range(10):
        planes_fit += "{}: {:.1f}˚\n".format(planes[i], refl[i])
    plt.text(
        xy_ref[0],
        xy_ref[1],
        s=planes_fit[:-1],
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="black"),
    )

    plt.axis(axes)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(nickel.basename + ".png")
    print("Results saved as " + nickel.basename + ".png")
    return
