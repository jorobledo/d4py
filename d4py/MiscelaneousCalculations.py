import numpy as np  # NumPy library
import sys


# --------1---------2---------3---------4---------5---------6---------7---------8
def getCylVolume(diameter=5.0, height=50.0):
    """
    Calculate the volume of a cylinder.

    This function returns the volume of a cylinder, given its diameter and height.

    The units of these two distances must be the same, and the result is given in
    the same unit^3 (if distances are in mm, the volume will be in mm^3).

    Input: diameter, height
      - diameter : float, optional
            The diameter of the cylinder. The default is 5.0.
        height : float, optional
            The height of the cylinder. The default is 50.0.

    Output: float
      - The volume of the cylinder.

    Author: Gabriel Cuello
    Date: Oct 29 2021
    --------------------------------------------------------------------------------
    """
    return np.pi * (diameter / 2.0) ** 2 * height


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def stop(message):
    """
    Stop the program at this point.
    This is not really useful for the Jupyter versions of the scripts, but it is
    for the python versions.

    Inputs: answer (from the console), message
     -  answer: string, 'n' or anything else
     -  message: string containing a message to be printed just before asking
        whether continuing or not.

    Output: nothing

    Author: Gabriel Cuello
    Date: Oct 25 2021
    --------------------------------------------------------------------------------
    """
    print()
    print(30 * "<>")
    print(message)
    answer = input("Do you want to continue? (y/n) :")
    if answer == "n":
        sys.exit("You stopped the program. Bye!")
    return None


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def ang2q(x, wlength=0.5):
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

    Author: Gabriel Cuello
    Date: Oct 24, 2021
    --------------------------------------------------------------------------------
    """
    q = 4.0 * np.pi / wlength * np.sin(np.array(x) / 2.0 * np.pi / 180.0)
    return q


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------8
def q2ang(x, wlength=0.5):
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
        - x: A list with the Q values (in 1/Å).
        - wlenght: The wavelength (in Å).

    Output:
        - angles: A list with the angles (in degrees).

    Author: Gabriel Cuello
    Date: Oct 24, 2021
    --------------------------------------------------------------------------------
    """
    ang = 360.0 / np.pi * np.arcsin(np.array(x) * wlength / 4.0 / np.pi)
    return ang


# --------1---------2---------3---------4---------5---------6---------7---------8


# --------1---------2---------3---------4---------5---------6---------7---------
def ratio(y1, e1, y2, e2):
    y_rat = []
    e_rat = []
    for i in range(len(y1)):
        if y2[i] != 0:
            y_rat.append(y1[i] / y2[i])
            e_rat.append(
                np.sqrt(y2[i] ** 2 * e1[i] ** 2 + y1[i] ** 2 * e2[i] ** 2) / y2[i] ** 2
            )
        else:
            y_rat.append(0.0)
            e_rat.append(0.0)
    return y_rat, e_rat


# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7---------
def wsum2(w1, data1, w2, data2):
    """
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
    """
    x1 = data1[:, 0]
    y1 = data1[:, 1]
    e1 = data1[:, 2]
    #    x2 = data1[:,0]
    y2 = data2[:, 1]
    e2 = data2[:, 2]
    if len(y1) != len(y2):
        print("--- Error in the binary sum (wsum2).")
        print("--- The input vectors have not the same length.")
        return
    length = len(y1)
    ysum = []
    esum = []
    if (w1 == 0) and (w2 == 0):
        for i in range(length):
            w = 0
            sqerr = e1[i] ** 2 + e2[i] ** 2
            if sqerr != 0:
                w = e2[i] ** 2 / sqerr
            ysum.append(w * y1[i] + (1 - w) * y2[i])
            esum.append(np.sqrt(w**2 * e1[i] ** 2 + (1.0 - w) ** 2 * e2[i] ** 2))
    elif w1 == 0:
        for i in range(length):
            ysum.append(w2 * y2[i])
            esum.append(np.sqrt(w2**2 * e2[i] ** 2))
    elif w2 == 0:
        if (w1 > 0) and (w1 <= 1):
            for i in range(length):
                ysum.append(w1 * y1[i] + (1 - w1) * y2[i])
                esum.append(np.sqrt(w1**2 * e1[i] ** 2 + (1.0 - w1) ** 2 * e2[i] ** 2))
        else:
            print("--- Error in the binary sum (wsum2).")
            print("--- The the weight of first set of data should be between 0 and 1.")
    else:  # both weights are different of zero and are simply considered as factors
        for i in range(length):
            ysum.append(w1 * y1[i] + w2 * y2[i])
            esum.append(np.sqrt(w1**2 * e1[i] ** 2 + w2**2 * e2[i] ** 2))
    ysum = np.array(ysum)
    esum = np.array(esum)
    x1 = x1.reshape(x1.shape[0], 1)
    ysum = ysum.reshape(ysum.shape[0], 1)
    esum = esum.reshape(esum.shape[0], 1)
    summed = np.concatenate((x1, ysum), axis=1)
    summed = np.concatenate((summed, esum), axis=1)
    return summed


# End of wsum2
# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7---------
def get_xlim(xmin, xmax, dbin):
    """
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
    """
    xini = xmin - dbin / 2.0  # coordinate of the left side of the first bin
    xfin = xmax + dbin / 2.0  # coordinate of the right side of the last bin
    nb = int((xfin - xini) / dbin + 1)  # number of bins
    # Initialising the two output lists
    x_lim = []
    x_bin = []
    for i in range(nb + 1):
        x_lim.append(xmin - dbin / 2.0 + i * dbin)
        x_bin.append(xmin + i * dbin)
    #    x_bin.pop() # Removing the last element of this list
    return nb, x_lim, x_bin


# End of det_xlim
# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7---------
def get_bins(x, xdel, x_lim):
    """
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
    """
    # Initialising the two output lists
    bins = []
    frac = []
    # Determines the size of a bin, which is constant
    dbin = x_lim[1] - x_lim[0]

    # The rectangle corresponding to the x value has a width of xdel,
    # half on each side of x
    x1 = x - xdel / 2.0  # coordinate of the left side  of the rectangle
    x2 = x + xdel / 2.0  # coordinate of the right side of the rectangle

    # Determining the bins where left and right sides of the rectangle fall
    b1 = int((x1 - x_lim[0]) / dbin)
    b2 = int((x2 - x_lim[0]) / dbin)
    if b1 < 0 or b2 >= len(x_lim):
        return bins, frac
    deltab = b2 - b1
    # There are 3 possible cases:
    #    (1) deltab = 0
    #        The rectangle falls completely in a single bin.
    #    (2) deltab = 1
    #        The rectangle falls on 2 bins, covering partially each one.
    #    (3) deltab > 1
    #        The rectangle falls on more than 2 bins, covering partially the
    #        first and the last ones, and completely the bins in between.
    if deltab == 0:  # Case (1)
        #   b1 (=b2) is the single bin where the rectangle falls.
        #   Then the fraction is equal to 1.0
        bins.append(b1)
        frac.append(1.0)
    elif deltab == 1:  # Case (2)
        f1 = (x_lim[b1 + 1] - x1) / xdel
        bins.append(b1)
        frac.append(f1)
        bins.append(b1 + 1)
        frac.append(1.0 - f1)
    elif deltab > 1:  # Case (3)
        f1 = (x_lim[b1 + 1] - x1) / xdel  # First bin
        bins.append(b1)
        frac.append(f1)
        for i in range(1, deltab):  # Intermediate bins
            bins.append(b1 + i)
            frac.append(dbin / xdel)
        f2 = (x2 - x_lim[b2]) / xdel  # Last bin
        bins.append(b2)
        frac.append(f2)
    else:
        print("ERROR in get_bins")

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
    return bins, frac


# End of get_bins
# --------1---------2---------3---------4---------5---------6---------7---------


# --------1---------2---------3---------4---------5---------6---------7---------
def rebin(xdel, wlength, data, xmin, xmax, dbin):
    """
    rebin function

    This function makes a rebinning of the experimetal data.

    The rebinning can be made in angular or Q-scale, depending on the input parameter wavelength.
        A non-positive wavelength produces a binning in scattering angle.
        A positive wavelength produces a binning in Q-scale. The wavelength must be in Å.

    The parameter xdel is the witdh of a channel in scattering angle, assumed constant for the
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
    """
    x_dat = data[:, 0]
    y_dat = data[:, 1]
    e_dat = data[:, 2]
    # Call get_xlim to obtain the number of bins, the limiting values for the bins and the values
    # of the bins
    nbins, x_lim, x_bin = get_xlim(xmin, xmax, dbin)

    # Creates lists for storing the new y, new error and the fraction
    y_bin = []
    e_bin = []
    f_bin = []

    # Initialise the lists that will serve as accumulators
    for i in range(nbins + 1):
        y_bin.append(0.0)
        e_bin.append(0.0)
        f_bin.append(0.0)

    if wlength <= 0:
        # The binning is in angular scale
        for i in range(len(x_dat)):
            if (np.isnan(y_dat[i]) is False) and (np.isnan(e_dat[i]) is False):
                # For each experimental x value, calls get_bins, which returns the bins covered by that
                # point and the corresponding fractions
                bins, frac = get_bins(x_dat[i], xdel, x_lim)
                # For each of these bins we add the corresponding fraction of y and error.
                # Because these fractions act as weighting factors, the fractions are accumulated for
                # normalisation purposes
                for j in range(len(bins)):
                    y_bin[bins[j]] += frac[j] * y_dat[i]
                    e_bin[bins[j]] += frac[j] * e_dat[i]
                    f_bin[bins[j]] += frac[j]
    else:
        # The binning is in Q-scale
        for i in range(len(x_dat)):
            if (np.isnan(y_dat[i]) is False) and (np.isnan(e_dat[i]) is False):
                # ${\rm d}Q = \frac{2\pi}{\lambda} \sqrt{1-\frac{Q\lambda}{4 \pi}}
                #  {\rm d}2\theta \frac{pi}{180}$
                qdel = (
                    2.0
                    * np.pi
                    / wlength
                    * np.sqrt(1.0 - (x_dat[i] * wlength / 4.0 / np.pi) ** 2)
                )
                qdel *= xdel * np.pi / 180.0
                # For each experimental x value, calls get_bins, which returns the bins covered by that
                # point and the corresponding fractions
                bins, frac = get_bins(x_dat[i], qdel, x_lim)
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
                    y_bin[bins[j]] += frac[j] * y_dat[i]
                    e_bin[bins[j]] += frac[j] * e_dat[i]
                    f_bin[bins[j]] += frac[j]

    # Normalisation of y and errors. If fraction is 0, that bin has no data.
    for i in range(nbins + 1):
        if f_bin[i] != 0:
            y_bin[i] /= f_bin[i]
            e_bin[i] /= f_bin[i]
    #            print(x_bin[i],y_bin[i],e_bin[i])
    x_bin = np.array(x_bin)
    y_bin = np.array(y_bin)
    e_bin = np.array(e_bin)
    #   Reshapes the arrays as 2D arrays, but with 1 column
    x_bin = x_bin.reshape(x_bin.shape[0], 1)
    y_bin = y_bin.reshape(y_bin.shape[0], 1)
    e_bin = e_bin.reshape(e_bin.shape[0], 1)
    #   Concatenates the arrays to have a 3-column matrix
    binned = np.concatenate((x_bin, y_bin), axis=1)
    binned = np.concatenate((binned, e_bin), axis=1)
    return binned


#    return x_bin,y_bin,e_bin
# End of rebin
# --------1---------2---------3---------4---------5---------6---------7---------
