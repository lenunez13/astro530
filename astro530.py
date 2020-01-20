import numpy as np
from astropy import constants as const

def Constants():
    '''
    Purpose: Returns several physical constants (as stored by Astropy) in SI units.
    
    Calling sequence: constants()
    
    Inputs: none
    
    Outputs: list containing the following constants:
             Gravitational constant G
             Speed of light c
             Planck's constant h
             Boltzmann constant k_b
             Stefan-Boltzmann constant sigma_sb
             
    Restrictions:
    
    Modifications:
    
    '''
    G = const.G.value
    c = const.c.value
    h = const.h.value
    k_B = const.k_B.value
    sigma_sb = const.sigma_sb.value
    
    return [G,c,h,k_B,sigma_sb]

G,c,h,k_B,sigma_sb = Constants()


def Planck(T,nu=None,lam=None,nu_tilde=None):
    '''
    Puropose: Calculates the Planck function. Returns a scalar value of specific intensity 
    (erg/sec/cm2/Hz/ster) for a given temperature and frequency/wavelength in accordance with
    equation 1.51 in "Radiative Processes in Astrophysics" by Rybicki & Lightman.
    
    Calling sequence: Planck(T,nu=None,lam=None)
    
    Inputs: Temperature: scalar "T", units K  
            Frequency: scalar "nu", units Hz
            Wavelength: scalar "lam", units cm
    
    Output(s): Specific intensity (erg/sec/cm2/Hz/ster) for the given temperature and frequency/wavelength.
    
    Restrictions:
    
    Modifications:
    '''
    
    if nu is not None or nu_tilde is not None:
        if nu_tilde is not None: 
            nu = nu_tilde*c
        B_nu = 2.*h*nu**3./c**2./(np.exp(h*nu/(k_B*T))-1.) #R&L Eq. (1.51)
        return B_nu
    
    elif lam is not None:
        B_lam = 2.*h*c**2./lam**5./(np.exp(h*c/(lam*k_B*T))-1.)
        return B_lam
    
def Integrate(funct,a,b,n=1e3):
    '''
    Purpose: Integrates a function using The Midpoint Rule (Stewart p. 378)
    
    Calling sequence: integrate(funct,a,b,n=100)
    
    Inputs: Array pertaining to the function: array "funct"
            Lower bound of integration: float "a"
            Upper bound of integration: float "b"
            Number of samples from the function: int "n"

    Output(s): Area under the function from a to b
    
    Restrictions: This function may be inaccurate when the limits of integration are 0 or
                  infinity, so the user must accordingly represent zero or infinity
                  using sufficiently small or large floats. 
    
    Modifications:
    '''
    dx = (b-a)/n
    x = np.linspace(a,b,n+1)
    Sum = 0
    for i in range(len(x)-1):
        dx = (b-a)/n
        midpoint = 0.5*(x[i]+x[i+1])
        Sum += funct(midpoint)*dx
    return Sum


        


    
    
    
    
    
    
    
    
    