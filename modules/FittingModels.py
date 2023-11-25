import numpy as np                   # NumPy library


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
    limited = np.concatenate((x,y),axis=1)
    limited = np.concatenate((limited,e),axis=1)
    return limited


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

#--------1---------2---------3---------4---------5---------6---------7---------8
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
--------------------------------------------------------------------------------
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

