import numpy as np  # NumPy library
from scipy import integrate  # For integration
from .FittingModels import Lorch, step, LorchN


# --------1---------2---------3---------4---------5---------6---------7-----
def getSineFT(vvec, yvec, uvec, vmax=0, c=1.0, s=1.0, w=0):
    # Number of points in the function to be transformed
    nbr_v = len(vvec)  # abcissa
    nbr_y = len(yvec)  # ordinate
    # These numbers must be equal, otherwise print error message
    if nbr_v != nbr_y:
        print("ERROR")
        return None

    #   The default value of vmax is 0
    #   If vmax less or equal to 0 or greater than the last element
    #   then vmax is the last element of the input v array.
    if (vmax <= 0) or (vmax > vvec[-1]):
        vmax = vvec[-1]

    win = np.ones(nbr_v)
    if w == 1:
        win = LorchN(vvec, vmax)

    stf = np.zeros(len(uvec))
    for i in range(len(uvec)):
        #        integ = soq*q*np.sin(q*r)
        integ = (yvec - s) * vvec * np.sin(vvec * uvec[i]) * win
        #        result = integrate.simps(integ,q)
        #        pcf.append(result)
        #        integ = soq*q*np.sin(q*r)*win
        stf[i] = 2.0 / np.pi * c * integrate.simps(integ, vvec)
    # Pair correlation function or G(r)
    return stf


# --------1---------2---------3---------4---------5---------6---------7-----


# --------1---------2---------3---------4---------5---------6---------7---------
def sineFT(StrFactor, nr, qmax=23.5, density=0.1, constant=1.0, selfsca=1.0, window=0):
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
    soq = StrFactor[:, 1]
    q = StrFactor[:, 0]
    pcf = []
    pdf = []
    rdf = []
    tor = []
    run = []
    #    infile = np.genfromtxt(filename, skip_header= 7) #creates an array from xA, yA, zA
    #    q  = infile[:, 0]
    #    soq  = infile[:, 1]
    #    #err  = infile[:, 2]

    win = []  # A list containing the values of the window function

    if window == 0:
        for qu in q:
            win.append(step(qu, qmax))
    else:
        for qu in q:
            win.append(Lorch(qu, qmax))
        #       Integral of the Lorch function for normalisation
        integralLorch = integrate.simps(win, q)
        #        print ('Normalisation of Lorch function: ',integralLorch/qmax)
        for i in range(len(win)):
            win[i] = win[i] * qmax / integralLorch

    deltaR = nr[1] - nr[0]
    for r in nr:
        #        integ = soq*q*np.sin(q*r)
        integ = (soq - selfsca) * q * np.sin(q * r) * win
        #        result = integrate.simps(integ,q)
        #        pcf.append(result)
        #        integ = soq*q*np.sin(q*r)*win
        result = 2.0 / np.pi * constant * integrate.simps(integ, q)
        # Pair correlation function or G(r)
        pcf.append(result)
        # Pair distribution function or g(r)
        if r <= 0:
            pdf.append(0.0)
        else:
            pdf.append(result / 4.0 / np.pi / r / density + 1.0)
        # Radial distribution function or RDF(r)
        rdf.append(result * r + 4.0 * np.pi * density * r**2)
        # T(r) = RDF(r)/r; symmetric peaks for fitting
        tor.append(result + 4.0 * np.pi * density * r)
    # Running integral of the RDF(r), i.e., integral from 0 to r
    for r in nr:
        ind = 1 + int(r / deltaR)
        xr = nr[0:ind]
        integ = np.array(rdf)[0:ind]
        result = integrate.simps(integ, xr)
        run.append(result)
    rrr = np.array(nr)
    pcf = np.array(pcf)
    pdf = np.array(pdf)
    rdf = np.array(rdf)
    tor = np.array(tor)
    run = np.array(run)
    rrr = rrr.reshape(rrr.shape[0], 1)
    pcf = pcf.reshape(pcf.shape[0], 1)
    pdf = pdf.reshape(pdf.shape[0], 1)
    rdf = rdf.reshape(rdf.shape[0], 1)
    tor = tor.reshape(tor.shape[0], 1)
    run = run.reshape(run.shape[0], 1)
    fou = np.concatenate((rrr, pcf), axis=1)
    fou = np.concatenate((fou, pdf), axis=1)
    fou = np.concatenate((fou, rdf), axis=1)
    fou = np.concatenate((fou, tor), axis=1)
    fou = np.concatenate((fou, run), axis=1)
    return fou


# --------1---------2---------3---------4---------5---------6---------7---------


def backFT(qu, ar, pdf, density, cut):
    pdf_cut = pdf.copy()
    #    print ('Cut-off = ',cut)
    # pdf_lorch_cut = pdf_lorch.copy()

    # for i in range(len(pdf)-1,0,-1):
    #     if pdf[i] < 0:
    #         break
    # for j in range(i+1):
    #     pdf_cut[j] = 0.0
    if len(cut) == 1 and cut[0] == -1:
        for i in range(len(pdf) - 1, 0, -1):
            if pdf[i] < 0:
                break
        for j in range(i + 1):
            pdf_cut[j] = 0.0
    elif len(cut) == 1 and cut[0] > 0.0:
        for i in range(len(pdf)):
            if ar[i] < cut[0]:
                pdf_cut[i] = 0.0
    elif len(cut) == 3:
        for i in range(len(pdf)):
            if ar[i] < cut[0] or (cut[1] < ar[i] and ar[i] < cut[2]):
                pdf_cut[i] = 0.0
    #                break

    #     width_fig = 10.
    #     plt.figure(figsize=(width_fig,2./3.*width_fig))
    #     plt.plot(ar,pdf_cut,label='Cut pdf')
    #     plt.title('Pair Distribution Function forced to be 0 in the repulsion region')
    #     plt.xlabel(r'$R$ (Å)')
    #     plt.ylabel('$g(R)$')
    #     plt.axis([0,10,-1,5])
    #     plt.grid(True)
    #     #plt.legend(bbox_to_anchor=(0.86, 0.35))
    #     plt.legend(loc='best')
    #     plt.show()
    # --------1---------2---------3---------4---------5---------6---------7-----
    # Calculation of the structure factor as the back Fourier transform of the
    # given pdf
    soq = 1 + 2.0 * np.pi**2 * density / qu * getSineFT(ar, pdf_cut, qu, w=0)
    # --------1---------2---------3---------4---------5---------6---------7-----
    #     width_fig = 10.
    #     plt.figure(figsize=(width_fig,2./3.*width_fig))
    #     plt.plot(ar,BAS00.tor,label=BAS00.sampleTitle+' tor')
    #     plt.plot(ar,BAS18.tor,label=BAS18.sampleTitle+' tor')
    #     plt.title('Radial Distribution Function $RDF(R)$')
    #     plt.xlabel(r'$R$ (Å)')
    #     plt.ylabel('$T(R)$')
    #     plt.axis([0,10,-1,10])
    #     plt.grid(True)
    #     #plt.legend(bbox_to_anchor=(0.86, 0.35))
    #     plt.legend(loc='best')
    #     plt.show()
    # --------1---------2---------3---------4---------5---------6---------7-----
    # Pair correlation function, G(R)
    pcf = getSineFT(qu, soq, ar, w=0)
    # --------1---------2---------3---------4---------5---------6---------7-----

    # Pair distribution function, g(R)
    pdf = 1 + pcf / 4.0 / np.pi / density / ar
    # Radial distribution function, RDF(R)
    rdf = 4.0 * np.pi * density * ar * ar * pdf
    # Linearised radial distribution function, T(R) = RDF(R)/R
    tor = rdf / ar

    # Running integral of the radial distribution function
    run = np.zeros(len(ar))
    for i in range(1, len(ar)):
        run[i] = integrate.simps(rdf[0:i], ar[0:i])
    return soq, pcf, pdf, rdf, tor, run
