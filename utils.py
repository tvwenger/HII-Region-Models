"""
General utilities for eRRL project.
Based on IDL code by A.A. Kepley

Functions:
voigt                    - calculate the voigt profile
rrl_freq                 - calculate the RRL frequency
doppler_width            - calculate doppler width
tau_c                    - calculate continuum optical depth
tau_l                    - calculate line optical depth
peak_line_fluxden        - calculate the peak RRL flux density
line_flux                - calculate the integrated line flux
min_num_HII_regions      - calculate the minimum number of HII regions
max_num_HII_regions      - calculate the maximum number of HII regions
num_HII_regions          - calculate the number of HII regions
num_HII_regions_los      - calculate the number of HII regions along
                           a line of sight
filling_factor           - calculate the HII region filling factor
thermal_fluxden          - calculate the thermal continuum flux
nonthermal_fluxden       - calculate the non-thermal continuum flux
                           using the measured continuum flux
spectral_info            - calculate spectral index and scaling factor
nonthermal_fluxden_model - calculate the nonthermal fluxden using the
                           spectral info from previous function
mass_ion                 - calculate mass of ionized gas
num_ion                  - calculate the number of H-ionizing photons per sec

Jun 2016 - Trey Wenger - Creation

dsb 21Jun2016 - Modify opacity threshold in thermal_fluxden().
dsb 30Sep2016 - Calculate mass_ion, num_ion.
dsb 02Dec2016 - Add max_num_HII_regions.
"""

import numpy as np
import math
import astropy.units as u
import astropy.constants as c
from scipy.special import wofz

import pdb
#pdb.set_trace()


def voigt(a,u):
    """
    Calculate the voigt profile at frequency offset u.
    Gives same result as IDL voigt function.

    The Voigt function is the real part of the Faddeeva function
    (wofz) evaluated at the complex position (u + ai)

    Inputs:
      a = Voigt damping constant
      u = offset from line center

    Returns:
      v = value of voigt function at u
    """
    z = (u + 1.j*a)
    v = wofz(z).real
    return v

def rrl_freq(fit_lines,deltan=1.):
    """
    Calculate the frequency of a hydrogen electronic transition from
    state n+deltan to n.
    Equation from Roelfsema and Goss 1992
    
    Inputs:
      fit_lines = array of lower energy levels
      deltan    = energy level transition

    Returns:
      freq_GHz = frequency of transition in GHz for each line
    """
    #
    # Rydberg constant for hydrogen (nuclear mass = proton mass)
    #
    Rx = c.Ryd/(1. + c.m_e/c.m_p)
    #
    # Transition frequency
    #
    freq = Rx * c.c * (1./fit_lines**2. - 1./(fit_lines+deltan)**2.)
    #
    # Convert to GHz
    #
    freq_GHz = freq.to('GHz').value
    return freq_GHz

def doppler_width(electron_temp,v_turbulent,freq):
    """
    Calculate the FWHM doppler width of a single HII region RRL
    Equation 13 from Viner, Vallee, and Hughes 1979

    Inputs:
      electron_temp = electron temperature in K
      v_turbulent   = turbulent velocity in km/s
      freq          = array of observing frequencies in GHz

    Returns:
      (deltaV_doppler, deltaNu_doppler)
      deltaV_doppler  = doppler FWHM in km/s
      deltaNu_doppler = doppler FWHM in GHz
    """
    # Calculate Doppler width of line in km/s
    # eq 13 in Viner, Vallee, and Hughes 1979
    deltaV_doppler = 1.66511 * np.sqrt(0.016499*electron_temp +
                                       2./3.*v_turbulent**2.)
    # and frequency in GHz
    deltaNu_doppler = ((deltaV_doppler*u.km/u.s)*freq*u.GHz/c.c)
    return (deltaV_doppler, deltaNu_doppler.to('GHz').value)
    

def tau_c(electron_temp,electron_dens,size,freq_obs,
               heII_abundance=0.,heIII_abundance=0.):
    """
    Calculate the continuum optical depth
    Equations from Viner, Vallee, and Hughes 1979 (eqs. 9,10)
    N.B. This equation agrees with Oster (1961)
         The Mezger and Henderson (1967) equation and Roeflsema and
         Goss (1992) equations are based on the approximation of
         Altenhoff et al. (1960).

    Inputs:
      electron_temp = electron temperature in K
      electron_dens = electron density in cm-3
      size          = line of sight distance through region in pc
      freq_obs      = array of frequencyies of observation in GHz
      heII_abundance  = abundance of singly ionized helium relative to
                        singly ionized hydrogen.
      heIII_abundance = abundance of doubly ionized helium relative to
                        singly ionized hydrogen.

    N.B. Found bug in AAK's IDL code. Missing * on line 164

    Returns:
      tau_c = continuum optical depth
    """
    #
    # defining constants (eq 10)
    #
    A = np.log(0.04955 * electron_temp**1.5 / freq_obs)
    n_iT = 1. + heII_abundance + 4.*heIII_abundance*(1.-np.log(2.)/A)
    n_iT = electron_dens * n_iT / (1. + heII_abundance +
                                   2.*heIII_abundance)
    #
    # equation 9, opacity
    #
    kappa_c = 9.7699e-21 * electron_temp**-1.5 * freq_obs**-2.
    kappa_c = kappa_c * electron_dens * n_iT * A
    #
    # convert opacity to optical depth
    #
    tau_c = kappa_c * (size*u.pc).to('cm').value
    return tau_c

def tau_l(electron_temp,electron_dens,size,
               fit_lines,deltaV_doppler,
               deltan=1,
               v_turbulent=0.0,use_voigt=False,
               heII_abundance=0.,heIII_abundance=0.,
               bn=1., betan=1.):
    """
    Calculate the line optical depth at line center
    Equations from Viner, Vallee, and Hughes 1979

    Inputs:
      electron_temp = electron temperature in K
      electron_dens = electron density in cm-3
      size          = line of sight distance through region in pc
      fit_lines     = array of lower energy level
      deltaV_doppler = FWHM doppler width (km/s)
      deltan        = energy level transition
      v_turbulent   = turbulent velocity in km/s
      use_voigt     = if True, use voigt profile
                      if False, use Gaussian profile
      heII_abundance  = abundance of singly ionized helium relative to
                        singly ionized hydrogen.
      heIII_abundance = abundance of doubly ionized helium relative to
                        singly ionized hydrogen.
      bn            = departure coefficient
      betan         = 1. - kT/(h*freq)*d(ln bn)/dn = (1. - gamma)

      N.B. note that Amanda says
      betan = (1. - kT/(h*freq))*d(ln bn)/dn
      which does not agree with Viner, Vallee, and Hughes. I don't
      think this will make a difference, though, since betan is
      calculated in the Fortran routines.

      Returns:
        tau_l = line optical depth at line center
    """
    #
    # Calculate oscillator strength
    # Based off approximation of Menzel 1969, eq 7, tables 3 and 4
    #
    if deltan == 1:
        M1 = 1.907749e-1
        M2 = -0.1877268
    elif deltan == 2:
        M1 = 2.633210e-2
        M2 = -0.3739199
    elif deltan == 3:
        M1 = 8.10562e-3
        M2 = -0.5218660
    else:
        raise ValueError("Not a valid deltan for oscillator strength calculation: {0} -- must be 3 or less.".format(deltan))
    oscillator_strength = fit_lines*M1*(1. + 1.5*deltan/fit_lines +
                                        M2/fit_lines**2.)
    #
    # Calculate height of line at line center
    #
    if use_voigt:
        #
        # Using voigt profile, equations 16, 17, 20, 21, 32, 69 in
        # Viner, Vallee and Hughes.
        # N.B. Amanda says something is wrong with this calculation
        # in her IDL code.
        #
        # Constants equation 32
        #
        TD = electron_temp + 40.42 * v_turbulent**2.
        #
        # Ratio of Lorentzian line width to Doppler line width, eq 69
        #
        deltaV_lor2dop = 9.76e-16 * electron_dens * electron_temp**(-0.5)
        deltaV_lor2dop *= TD**(-0.5) * fit_lines**7.0 / deltan
        deltaV_lor2dop *= (-11.4 + np.log(electron_temp * fit_lines/deltan))
        #
        # More constants, equation 17, 20, 21
        #
        A2 = np.sqrt(np.log(2.)) * deltaV_lor2dop
        #
        # Calculate voigt at line center i.e. V1 = 0 in eq. 16.)
        #
        height_center = voigt(A2,0.)
    else:
        #
        # No pressure broadening, eq 22 in Viner, Valee, and Hughes
        # at line center (freq_obs = freq_RRL)
        #
        height_center = 1.0
    #
    # Constants eq 12-15.
    # N.B. Factor of ni missing in ni equation. Perhaps ni = nH+ is
    # implicit?
    #
    ni = electron_dens / (1. + heII_abundance + 2.*heIII_abundance)
    A1 = 1.
    E1 = 157803. * A1 / fit_lines**2. / electron_temp
    #
    # line opacity at line center (eq 11) with bn = 1. and gamma = 0.
    #
    kappa_l = 1.4854e-22 * electron_dens * ni * fit_lines**2.
    kappa_l *= oscillator_strength
    kappa_l *= electron_temp**(-2.5) * deltaV_doppler**(-1.)
    kappa_l *= np.exp(E1) * height_center * bn * betan
    #
    # Convert to line optical depth at line center
    #
    tau_l = kappa_l * (size*u.pc).to('cm').value
    return tau_l

def peak_line_fluxden(electron_temp,omega_HII_region,omega_region,
                            rrl_freq,tau_c,tau_l_star,tau_l,bn=1.,
                            background_fluxden=0.):
    """
    Calculate the peak line flux density from an HII region and background
    Equation 1 from Anantharamaiah+1993

    N.B. Bug in AAK's IDL code. Need to keep double precision in exponentials
    to get accurate results. Luckily python uses double precision by default!

    Inputs:
      electron_temp = electron temperature (K)
      omega_HII_region = solid angle of single HII region (sterradians)
      omega_region     = solid angle of CASA region
      rrl_freq      = array of frequency of RRLs (GHz)
      tau_c         = array of continuum optical depths at each RRL frequency
      tau_l_star    = array of LTE line optical depths at each RRL frequency
      tau_l         = array of non-LTE line optical depths at each RRL frequency
      bn            = array of departure coefficient for each transition
      background_fluxden = background flux density (W/m2/Hz)

    Returns:
      (peak_line_fluxden_HII_region, peak_line_fluxden_background)
      peak_line_fluxden_HII_region = peak line flux density from the HII region (W/m2/Hz)
      peak_line_fluxden_background = peak line flux density from background (W/m2/Hz)
    """
    # First term in equation 1
    electron_temp = electron_temp * u.K
    rrl_freq = rrl_freq * u.GHz
    background_fluxden = background_fluxden * u.W/u.m**2./u.Hz
    peak_line_fluxden_HII_region = (bn*tau_l_star + tau_c)/(tau_l + tau_c)
    peak_line_fluxden_HII_region *= (1. - np.exp(-1.*(tau_l + tau_c)))
    peak_line_fluxden_HII_region += -1.*(1. - np.exp(-1.*tau_c))
    peak_line_fluxden_HII_region = peak_line_fluxden_HII_region * 2.*c.k_B*electron_temp*rrl_freq**2./c.c**2.*omega_HII_region
    # Second term in equation 1, including beam dilution correction
    peak_line_fluxden_background = background_fluxden*0.5*omega_HII_region/omega_region
    peak_line_fluxden_background = peak_line_fluxden_background * np.exp(-1.*tau_c) * (np.exp(-1.*tau_l) - 1.)
    # convert to W/m2/Hz

    return (peak_line_fluxden_HII_region.to('W/(m2 Hz)').value,
            peak_line_fluxden_background.to('W/(m2 Hz)').value)

def line_flux(peak_line_fluxden,deltaNu_doppler):
    """
    Calculate the integrated line flux assuming the line is Gaussian
    Area under a Gaussian is
    Area = sqrt(2\pi) * A * sigma
       where A = amplitude, sigma = width
    Since, FWHM = 2 * sqrt(2. * log(2)) * sigma
    Area = 0.5*sqrt(pi/log(2)) * A * FWHM

    Inputs:
      peak_line_fluxden = array of amplitude of Gaussian lines flux in W/m2/Hz
      deltaNu_doppler = array of doppler width of lines in GHz

    Returns:
      line_flux = line flux in W/m2
    """
    peak_line_fluxden = peak_line_fluxden * u.W/u.m**2./u.Hz
    deltaNu_doppler = deltaNu_doppler * u.GHz
    const = 0.5*np.sqrt(np.pi/np.log(2.))
    line_flux = const * peak_line_fluxden * deltaNu_doppler
    return line_flux.to('W/m2').value

def min_num_HII_regions(measured_line_fwhm,deltaV_doppler,
                             measured_omegaB,omega_region):
    """
    Calculate the minimum number of HII regions given the observed
    line width and model line width

    From equations on page 592 in Anantharamaiah+1993

    Inputs:
      measured_line_fwhm  = measured line FWHM (km/s)
      deltaV_doppler   = FWHM doppler RRL width of model HII region (km/s)
      measured_omegaB  = measured beam solid angle for this line (sr)
      omega_region     = solid angle of the model HII region (sterradian)

    Returns:
      min_num_HII_regions = minimum number of HII regions
    """
    min_num_HII_regions = (measured_line_fwhm/deltaV_doppler)*(omega_region/measured_omegaB)
    return min_num_HII_regions


def max_num_HII_regions(region_size, HII_region_size, num_HII_regions):
    """
    Calculate the maximum number of HII regions

    Volume of region > Volume of N HII regions

    Inputs:
      region_size      = size of line emitting region (pc)
      HII_region_size  = size of singel HII region (pc)
      num_HII_regions  = number of HII regions

    Returns:
      max_num_HII_regions = maximum number of HII regions
    """
    max_num_HII_regions = (region_size/HII_region_size)**3/num_HII_regions
    return max_num_HII_regions

    
def num_HII_regions(measured_line_flux,line_flux,
                         line_flux_override=None):
    """
    Calculate the number of HII regions within the extracted CASA
    region.

    From equations on page 592 in Anantharamaiah+1993

    Inputs:
      measured_line_flux = measured line flux (W/m2)
      line_flux   = the integrated RRL flux for a single HII region (W/m2)
      line_flux_override = use this line flux instead of the data
                            to estimate number of HII regions if set (W/m2)

    Returns:
      num_HII_regions = number of HII regions in beam
    """
    if line_flux_override is None:
        target_flux = measured_line_flux
    else:
        target_flux = line_flux_override
    #
    # Calculate the number of HII regions
    #
    num_HII_regions = target_flux / line_flux
    return num_HII_regions

def num_HII_regions_los(num_HII_regions,
                             HII_region_size,region_size):
    """
    Calculate the number of HII regions along a line of sight.

    From equations on page 592 in Anantharamaiah+1993

    Inputs:
      num_HII_regions = number of HII regions in beam
      HII_region_size = the linear size of the model HII region (pc)
      region_size = the linear size of the CASA region (pc)

    Returns:
      num_HII_regions_los = number of HII regions along a line of sight
    """
    #
    # Calculate the number of HII regions along a line of sight
    #
    num_HII_regions_los = num_HII_regions * (HII_region_size/region_size)**2.
    return num_HII_regions_los

def filling_factor(num_HII_regions,
                        HII_region_size,region_size):
    """
    Calculate the HII region filling factor.

    From equations on page 592 in Anantharamaiah+1993

    Inputs:
      num_HII_regions = number of HII regions in beam
      HII_region_size = the linear size of the model HII region (pc)
      region_size = the linear size of the CASA region (pc)

    Returns:
      filling_factor = HII region filling factor
    """
    #
    # Calculate the filling factor
    # 
    filling_factor = num_HII_regions * (HII_region_size/region_size)**3.
    return filling_factor
    
def thermal_fluxden(contfreqs,electron_temp,omega_HII_region,omega_region,
                         tau_c_cont,num_HII_regions,num_HII_regions_los):
    """
    Calculate the thermal flux density from a collection of HII regions.
    Equations 2-4 in Anantharamaiah+1993

    Inputs:
      contfreqs       = array of continuum frequencies (GHz)
      electron_temp   = electron temperature (K)
      omega_HII_region = the solid angle of the model HII region (sr)
      omega_region     = the solid angle of the CASA region (sr)
      tau_c_cont      = array of continuum optical depths at each
                        continuum frequency
      num_HII_regions = number of HII regions in model
      num_HII_regions_los = number of HII regions along LOS

    Returns:
      thermal_fluxden = thermal flux density (W/m2/Hz)
    """
    thermal_fluxden = np.zeros(len(contfreqs))
    electron_temp *= u.K
    contfreqs = contfreqs * u.GHz
    #
    # Differentiate between the low and high cutoffs
    # N.B. The paper says tau_c << 1, but AAK is using tau_c_los << 1
    #
    tau_c_cont_los = tau_c_cont * num_HII_regions_los
    low = tau_c_cont < 0.1
    high = tau_c_cont > 0.1
    #low = tau_c_cont_los < 0.1
    #high = tau_c_cont_los > 0.1
    #
    # low optical depth case; equation 2
    #
    thermal_fluxden_low = 2.*c.k_B*electron_temp*contfreqs[low]**2./c.c**2.*omega_HII_region
    thermal_fluxden_low *= (1. - np.exp(-1.*tau_c_cont[low]))*num_HII_regions
    thermal_fluxden[low] = thermal_fluxden_low.to('W/(m2 Hz)').value
    #
    # high optical depth case; equation 4
    #
    thermal_fluxden_high = 2.*c.k_B*electron_temp*contfreqs[high]**2./c.c**2.*omega_region
    thermal_fluxden_high *= (1. - np.exp(-1.*tau_c_cont_los[high]))
    thermal_fluxden[high] = thermal_fluxden_high.to('W/(m2 Hz)').value
    return thermal_fluxden

def nonthermal_fluxden(contfluxden,thermal_fluxden,tau_c_cont,
                            num_HII_regions_los):
    """
    Calculate the non-thermal flux density and unattenuated
    non-thermal flux density from a collection of HII regions.
    Equations 5 and 6 in Anantharamaiah+1993

    Inputs:
      contfluxden     = array of observed continuum flux densities (mJy)
      thermal_fluxden = thermal flux density 
      tau_c_cont      = numpy array of continuum optical depths at each
                        continuum frequency
      num_HII_regions_los = number of HII regions in model along los

    Returns:
      (nonthermal_fluxden, unatten_nonthermal_fluxden)
      nonthermal_fluxden         = non-thermal flux density (W/m2/Hz)
      unatten_nonthermal_fluxden = unattenuated non-thermal flux density (W/m2/Hz)
    """
    contfluxden = contfluxden * u.mJy
    thermal_fluxden = thermal_fluxden * u.W/u.m**2./u.Hz
    #
    # non-thermal flux density (eq 5)
    #
    nonthermal_fluxden = contfluxden - thermal_fluxden
    #
    # unattenuated non-thermal flux density (eq 6)
    #
    unatten_nonthermal_fluxden = nonthermal_fluxden / (1. - 0.5*num_HII_regions_los*(1.-np.exp(-1.*tau_c_cont)))

    return (nonthermal_fluxden.to('W/(m2 Hz)').value,
            unatten_nonthermal_fluxden.to('W/(m2 Hz)').value)

def spectral_info(freq,fluxden):
    """
    Calculate the spectral index and scaling factor for flux density vs freq

    Inputs:
      freq    = frequencies of observation (GHz)
      fluxden = observed flux density (W/m2/Hz)

    Returns:
      (spectral_index,scaling_factor)
      spectral_index = log10(high_freq_fluxden/low_freq_fluxden)/log10(high_freq/low_freq)
      scaling_factor = low_freq_fluxden/(low_freq^spectral_index) (W/m2/Hz^(1+spectral_index))
    """
    min_ind = np.argmin(freq)
    max_ind = np.argmax(freq)
    spectral_index = np.log10(fluxden[max_ind]/fluxden[min_ind])/np.log10(freq[max_ind]/freq[min_ind])
    scaling_factor = fluxden[min_ind]/(freq[min_ind]**spectral_index)
    return (spectral_index, scaling_factor)

def nonthermal_fluxden_model(spectral_index,scaling_factor,freq,
                             tau_c_cont,num_HII_regions_los):
    """
    Calculate the non-thermal flux density and unattenuated
    non-thermal flux density from a collection of HII regions using
    the scaling factor and spectral index.
    Equations 5 and 6 in Anantharamaiah+1993

    Inputs:
      spectral_index  = nonthermal spectral index
      scaling_factor  = nonthermal scaling factor (W/m2/Hz^(1+spectral_index))
      freq            = frequencies of observation
      tau_c_cont      = numpy array of continuum optical depths at each
                        continuum frequency
      num_HII_regions_los = number of HII regions in model along los

    Returns:
      (nonthermal_fluxden, unatten_nonthermal_fluxden)
      nonthermal_fluxden         = non-thermal flux density (W/m2/Hz)
      unatten_nonthermal_fluxden = unattenuated non-thermal flux density (W/m2/Hz)
    """
    #
    # unattenuated non-thermal flux density (eq 6)
    #
    unatten_nonthermal_fluxden = scaling_factor * freq**spectral_index
    #
    # non-thermal flux density (eq 5)
    #
    nonthermal_fluxden = 1. - 0.5*num_HII_regions_los*(1.-np.exp(-1.*tau_c_cont))
    nonthermal_fluxden *= unatten_nonthermal_fluxden
    return (nonthermal_fluxden,unatten_nonthermal_fluxden)

def mass_ion(num_HII_regions,HII_region_size,electron_dens):
    """
    Calculate the mass of ionized gas

    Inputs:
      num_HII_regions = number of HII regions in beam
      HII_region_size = the linear size of the model HII region (pc)
      electron_dens = electron density (cm-3)

    Returns:
      mass_ion = mass of ionzed gas (solar masses)
    """
    HII_region_size *= u.parsec
    electron_dens *= u.cm**(-3)
    m_h = 1.67352e-24*u.gram  # did not see m_h in astropy constants
    mass_ion = num_HII_regions*(math.pi/6.0)*HII_region_size.to(u.cm)**3*electron_dens*m_h
    return mass_ion.to(u.Msun).value

def num_ion(num_HII_regions,HII_region_size,electron_dens,electron_temp):
    """
    Calculate the number of H-ionizing photons

    Inputs:
      num_HII_regions = number of HII regions in beam
      HII_region_size = the linear size of the model HII region (pc)
      electron_dens = electron density (cm-3)
      electron_temp = electron temperature (K)

    Returns:
      num_ion = number of H-ionizing photons (s-1)
    """
    electron_dens /= u.cm**3
    HII_region_size *= u.parsec
    # Assume case B recombination rate (cm^3 s-1)
    # Use Hui & Gnedin (1997) approximation
    lambda_hi = 315614.0/electron_temp
    alpha_B = (2.753e-14*(lambda_hi**1.5/(1.0 + (lambda_hi/2.740)**0.407)**2.242))*u.cm**3/u.second
    num_ion = num_HII_regions*(math.pi/6.0)*HII_region_size.to(u.cm)**3*electron_dens**2*alpha_B
    return num_ion.value

