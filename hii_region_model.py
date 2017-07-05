"""
hii_region_model.py
Class and code for HII region models
Based on IDL code by A.A. Kepley

Trey Wenger May 2016

dsb 21Jun2016 - Edit comment in RSS routine
dsb 11Jul2016 - Fixed bug in call to get non-thermal emission (use los value)
                Modfiy non-thermal emission (unattenuated) constraint
                Add constraint to check for nan's in model properties
dsb 30Sep2016 - Add mass_ion, num_ion, SFR
dsb 02Dec2016 - Add max_num_HII_regions.
"""
import numpy as np
import math
import utils
from departure_coeffs import generate_departure_coeffs
import astropy.units as u
import astropy.constants as c

import pdb
#pdb.set_trace()


class HIIRegionModel:
    """
    HII Region Model Class
    """
    def __init__(self,linedata,contdata,omega_region,distance,
                 electron_temp=5000.,electron_dens=10.,
                 HII_region_size=1.,num_HII_regions=None,
                 v_turbulent=0.,heII_abundance=0.,heIII_abundance=0.,
                 n_for_num=58,fit_lines=[],background_fluxden=0.,
                 line_flux_override=None,use_voigt=True):
        """
        Initialize this HII Region Model.

        Inputs:
          linedata        = numpy array containing line data
          contdata        = numpy array containing continuum data
          omega_region    = solid angle of CASA region used to extract (arcsec2)
          distance        = distance to object (Mpc)
          electron_temp   = electron temperature (K)
          electron_dens   = electron density (cm-3)
          HII_region_size = HII region size (pc)
          num_HII_regions = number of HII regions
                            if None, estimate the number of HII regions
                            from the data and the model
          v_turbulent     = assumed turbulent velocity (km/s)
          heII_abundance  = abundance of singly ionized helium relative
                            to singly ionized hydrogen.
          heIII_abundance = abundance of doubly ionized helium relative
                            to singly ionized hydrogen.
          n_for_num       = RRL alpha line to use to calculate num HII regions
          fit_lines       = all RRL alpha lines to model
          background_fluxden = assumed background continuum flux density (mJy)
          line_flux_override = use this line flux to estimate number of
                               HII regions if set (W/m2)
          use_voigt          = if True, use a Voigt profile, otherwise Gaussian

        Returns:
          Nothing
        """
        self.linedata = linedata
        self.contdata = contdata
        self.distance = distance
        self.electron_temp = electron_temp
        self.electron_dens = electron_dens
        self.HII_region_size = HII_region_size
        self.num_HII_regions = num_HII_regions
        self.v_turbulent = v_turbulent
        self.heII_abundance = heII_abundance
        self.heIII_abundance = heIII_abundance
        self.n_for_num = n_for_num
        self.use_voigt = use_voigt
        self.line_flux_override = line_flux_override
        #
        # if fit_lines is empty, use all RRLs to calculate residuals
        #
        if len(fit_lines) == 0:
            self.fit_lines = linedata['n']
        else:
            self.fit_lines = np.array(fit_lines)
        # find indicies of fit_lines in linedata
        self.fit_lines_ind = np.zeros(len(fit_lines))
        for i,n in enumerate(self.fit_lines):
            myrow = np.where(self.linedata['n'] == n)[0]
            if len(myrow) != 1:
                raise ValueError("Problem with level {0} in linedata.".format(n))
            self.fit_lines_ind = myrow[0]
        # Find index of n_for_num in linedata
        myrow = np.where(self.linedata['n'] == self.n_for_num)[0]
        if len(myrow) != 1:
            raise ValueError("Problem with level {0} in linedata".format(self.n_for_num))
        self.n_for_num_linedata_ind = myrow[0]
        # Find index of n_for_num in fitlines
        myrow = np.where(self.fit_lines == self.n_for_num)[0]
        if len(myrow) != 1:
            raise ValueError("Problem with level {0} in fit_lines".format(self.n_for_num))
        self.n_for_num_fit_lines_ind = myrow[0]
        #
        # Hard code some parameters
        #
        self.deltan = 1 # only using alpha transitions
        # convert background flux density from mJy to W/m2/Hz
        self.background_fluxden = (background_fluxden*u.mJy).to('W/(m2 Hz)').value
        # convert solid angle of CASA region to sterrarians
        self.omega_region = (omega_region*u.arcsec**2.).to('sr').value
        # linear size of CASA region in pc
        self.region_size = np.sqrt(4. * self.omega_region * (self.distance * 1.e6)**2./np.pi)

    def calc_HII_region_properties(self):
        """
        Calculate physical properties for this HII region model

        Inputs:
          None

        Returns:
          None
        """
        #
        # Solid angle of single HII region (sterradians)
        #
        self.omega_HII_region = np.pi * (0.5*self.HII_region_size/(self.distance*1.e6))**2.
        #
        # Departure coefficients
        #
        self.bn,self.betan = generate_departure_coeffs(self.electron_temp,self.electron_dens,self.fit_lines)
        #
        # RRL frequency
        #
        self.rrl_freq = utils.rrl_freq(self.fit_lines,self.deltan)
        #
        # doppler width of a single HII region
        #
        self.deltaV_doppler, self.deltaNu_doppler = \
          utils.doppler_width(self.electron_temp,self.v_turbulent,
                                   self.rrl_freq)
        #
        # Continuum optical depth at RRL frequency and continuum frequencies
        #
        self.tau_c_rrl = utils.tau_c(self.electron_temp,self.electron_dens,
                                          self.HII_region_size,self.rrl_freq,
                                          heII_abundance=self.heII_abundance,heIII_abundance=self.heIII_abundance)
        self.tau_c_cont = utils.tau_c(self.electron_temp,self.electron_dens,
                                           self.HII_region_size,self.contdata['freq_GHz'],
                                           heII_abundance=self.heII_abundance,heIII_abundance=self.heIII_abundance)
        #
        # LTE line opacity
        #
        self.tau_l_star = utils.tau_l(self.electron_temp,self.electron_dens,
                                            self.HII_region_size,self.fit_lines,self.deltaV_doppler,deltan=self.deltan,
                                            v_turbulent=self.v_turbulent,use_voigt=self.use_voigt,
                                            heII_abundance=self.heII_abundance,heIII_abundance=self.heIII_abundance,
                                            bn=1.,betan=1.)
        #
        # non-LTE line opacity
        #
        self.tau_l = utils.tau_l(self.electron_temp,self.electron_dens,
                                      self.HII_region_size,self.fit_lines,self.deltaV_doppler,deltan=self.deltan,
                                      v_turbulent=self.v_turbulent,use_voigt=self.use_voigt,
                                      heII_abundance=self.heII_abundance,heIII_abundance=self.heIII_abundance,
                                      bn=self.bn,betan=self.betan)
        #
        # peak line flux density from a single HII region
        #
        self.peak_line_fluxden_HII_region,self.peak_line_fluxden_background = \
          utils.peak_line_fluxden(self.electron_temp,self.omega_HII_region,
                                       self.omega_region,self.rrl_freq,
                                       self.tau_c_rrl,self.tau_l_star,self.tau_l,
                                       bn=self.bn,background_fluxden=self.background_fluxden)
        self.peak_line_fluxden = self.peak_line_fluxden_HII_region + self.peak_line_fluxden_background
        #
        # line flux from a single HII region
        #
        self.line_flux = \
          utils.line_flux(self.peak_line_fluxden,self.deltaNu_doppler)
        #
        # minimum number of HII regions
        #
        self.min_num_HII_regions = \
          utils.min_num_HII_regions(self.linedata['fwhm_kms'][self.n_for_num_linedata_ind],
                                         self.deltaV_doppler,
                                         self.linedata['omegaB_sr'][self.n_for_num_linedata_ind],
                                         self.omega_region)
        #
        # number of HII regions; only calculate this if we were not
        # specified with a number of HII regions
        #
        if self.num_HII_regions is None:
            self.num_HII_regions = \
              utils.num_HII_regions(self.linedata['intflux_Wm2'][self.n_for_num_linedata_ind],
                                         self.line_flux[self.n_for_num_fit_lines_ind],
                                         line_flux_override=self.line_flux_override)
        #
        # maximum number of HII regions
        # check that num_HII_regions is sensible
        #
        if math.isnan(self.num_HII_regions):
            self.max_num_HII_regions = 1.0e10
        else:
            self.max_num_HII_regions = \
              utils.max_num_HII_regions(self.region_size,  self.HII_region_size, self.num_HII_regions)
        #
        # Number of HII regions along LOS
        #
        self.num_HII_regions_los = \
          utils.num_HII_regions_los(self.num_HII_regions,self.HII_region_size,
                                         self.region_size)
        #
        # HII region filling factor
        #
        self.filling_factor = \
          utils.filling_factor(self.num_HII_regions,self.HII_region_size,
                                    self.region_size)
        #
        # Calculate the thermal flux density
        #
        self.thermal_fluxden_cont = \
          utils.thermal_fluxden(self.contdata['freq_GHz'],self.electron_temp,
                                     self.omega_HII_region,self.omega_region,
                                     self.tau_c_cont,self.num_HII_regions,self.num_HII_regions_los)
        self.thermal_fluxden_rrl = \
          utils.thermal_fluxden(self.rrl_freq,self.electron_temp,
                                self.omega_HII_region,self.omega_region,
                                self.tau_c_rrl,self.num_HII_regions,self.num_HII_regions_los)
        #
        # Calculate the non-thermal flux density
        #
        self.nonthermal_fluxden_cont,self.unatten_nonthermal_fluxden_cont = \
          utils.nonthermal_fluxden(self.contdata['S_C_mJy'],self.thermal_fluxden_cont,
                                        self.tau_c_cont,self.num_HII_regions_los)
        #
        # Calculate nonthermal spectral index and scaling factor
        #
        self.nonthermal_spectral_index, self.nonthermal_scaling_factor = \
          utils.spectral_info(self.contdata['freq_GHz'],
                                   self.unatten_nonthermal_fluxden_cont)
        #
        # Calculate the non-thermal flux density using the scaling
        # factor and spectral index
        #
        self.nonthermal_fluxden_rrl,self.unatten_nonthermal_fluxden_rrl = \
          utils.nonthermal_fluxden_model(self.nonthermal_spectral_index,self.nonthermal_scaling_factor,
                                         self.rrl_freq,
                                         self.tau_c_rrl,self.num_HII_regions_los)
        # 
        # Calculate the mass of ionized gas
        #
        self.mass_ion = \
          utils.mass_ion(self.num_HII_regions,self.HII_region_size,
                         self.electron_dens)
        # 
        # Calculate the number of H-ionizing photons per sec
        #
        self.num_ion = \
          utils.num_ion(self.num_HII_regions,self.HII_region_size,
                         self.electron_dens,self.electron_temp)
        #
        # Calculate the star formation rate (SFR) 
        #  (Murphy et al. 2011 using Starburst99 code)
        #
        self.sfr = 7.29e-54*self.num_ion


    def check(self,verbose=False):
        """
        Check this model to make sure physical parameters are sane

        Inputs:
          verbose = if True, output info about failed checks

        Returns:
          True if everything is good
          False if there is a problem
        """
        good = True
        #
        # Check that filling factor is less than 1
        #
        if self.filling_factor > 1.0:
            if verbose: print("Filling factor > 1.0")
            good = False
        #
        # check that peak line intensity of single HII region is less
        # than observed line intensity
        #
        peak_line_fluxden = self.peak_line_fluxden*u.W/u.m**2./u.Hz
        if np.any(peak_line_fluxden.to('mJy').value > self.linedata['peak_mJy'][self.fit_lines_ind]):
            if verbose: print("HII region peak line flux > observed peak line flux")
            good = False
        #
        # Check that the estimated number of HII regions is greater than
        # the minimum possible number of HII regions
        #
        if self.min_num_HII_regions > self.num_HII_regions:
            if verbose: print("minimum number of HII regions is greater than number of HII regions")
            good = False
        #
        # Check that the estimated number of HII regions is less than
        # the maximum possible number of HII regions
        #
        if self.max_num_HII_regions < self.num_HII_regions:
            if verbose: print("maximumx number of HII regions is less than number of HII regions")
            good = False
        #
        # Check that thermal emission is not overproduced
        # (i.e. check that nonthermal is positive)
        #
        if np.any(self.nonthermal_fluxden_cont < 0.):
            if verbose: print("overproduces thermal continuum emission (attenuated)")
            good = False
        #
        # Check that the approximations used to calculate the unattenuated
        # non-thermal emission are valid (N.B., AAK did not check this). 
        #
        if np.any(self.tau_c_rrl > 0.5) and (self.num_HII_regions_los > 1.0):
            if verbose: print("approximation used to calculate unattenuated non-thermal emission fails")
        else:
            #
            # Check that unattenuated nonthermal spectral index is greater
            # than -1.5
            #
            if self.nonthermal_spectral_index < -1.5:
                if verbose: print("nonthermal spectrum too steep")
                good = False

        # Check if we produce nan's in some of the model properties
        if np.any(np.isnan(self.line_flux)):
            if verbose: print("line flux is nan")
            good = False
        if math.isnan(self.num_HII_regions):
            if verbose: print("number of HII regions is nan")
            good = False
        return good


    def rss(self):
        """
        Calculate the residual sum of squares for this model and the data.
    
        N.B. you can not use a chi-squared test here. The model is
        non-linear with the parameters, the parameters are probably
        not independent. 
    
        Inputs:
          linedata   = RRL observation data
          fit_lines  = lines used in fit
          fit_fluxes = fluxes from model (W/m2)
    
        Returns:
          rss      = residual sum of squares

        FYI with chisq:
        N.B. AAK was doing (observed-theory)^2/observed_err^2
        it should be (observed-theory)^2/observed_err^2/num_dof
        """
        rss = np.sum((self.linedata['intflux_Wm2'][self.fit_lines_ind] - self.line_flux*self.num_HII_regions)**2.)
        return rss    

