"""
departure_coeffs.py
Routine to call Fortran program to generate departure coefficients
for non-LTE situations.
Based off IDL code by A.A. Kepley

Trey Wenger May 2016
"""
import numpy as np
import subprocess

def generate_departure_coeffs(electron_temp,electron_dens,fit_lines):
    """
    Run bsbn program to calculate departure coefficients.
    
    Inputs:
      electron_temp = electron temperature in K
      electron_dens = electron density in cm-3
      fit_lines     = array of energy levels of interest

    Returns:
      (bn,betan)
      bn          = departure coefficient for each fit_line
      betan       = 1. - gamma for each fit_line

    N.B. bsbn.f is a Fortran routine for calculating the departure
    coefficients and betas to correct for non-LTE effects. It is
    very similar to the routine in Appendix E.1 of Radio Recombination
    Lines: Their Physics and Astronomical Applications by Gordon
    and Sorochenko (2009). This code is based on Salem and Brocklehurst
    (1979) and Brocklehurst and Salem (1977).

    According to Gorden and Sorochenko, this code is good to about
    H40. At lower energy levels, you should use Storey and Hummer (1995).

    Right now we are only changing the temperature and density. The
    file call_bsbn.f prompts the user for temperature and density. This
    python function calls call_bsbn.f and supplies it with temperature
    and density.

    The file bsbn.inp controls other parameters, including the min
    and max levels to calculate (third and fourth values on last line)
    and the increment (second value on last line).

    These codes should be compiled using
    ifort -noautomatic -o call_bsbn.o bsbn.f call_bsbn.f
    """
    #
    # Call Fortran routine
    #
    run = subprocess.run(["./call_bsbn.o",format(electron_temp),format(electron_dens)])
    #
    # Read bsbn.f output and get bn and betan for our transition
    #
    output = np.genfromtxt('bsbn.out',dtype=None,skip_header=4,
                           usecols=(0,1,2),names='n,bn,betan')
    bn = np.zeros(len(fit_lines))
    betan = np.zeros(len(fit_lines))
    for i,n in enumerate(fit_lines):
        myrow = np.where(output['n'] == n)[0]
        if len(myrow) != 1:
            raise ValueError("Problem with level {0} in bsbn.out".format(n))
        bn[i] = output['bn'][myrow[0]]
        betan[i] = output['betan'][myrow[0]]
    return (bn,betan)
