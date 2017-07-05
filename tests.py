"""
Unit tests for eRRL project

The following functions test the functions in utils.py with values
calculated by hand using the equations from the references.

It also checks against results calculated using IDL code from AAK.

test_voigt 
test_rrl_line_freq
test_doppler_width
test_tau_c
test_tau_l
test_peak_line_fluxden
test_line_flux
test_num_HII_regions
test_thermal_fluxden
test_nonthermal_fluxden
test_spectral_info

Trey Wenger May 2016
"""
import unittest
import utils
import numpy as np
from fit_model import fit_model

class TestUtils(unittest.TestCase):

    def test_voigt_1(self):
        """
        Test the voigt function:
        H(a,u) = a/pi * integral(-inf,inf)(e^-y^2 dy/(a^2 + (u-y)^2))
        Note that these answers agree with the IDL function voigt()
        """
        a = 5.
        u = 0.
        voigt_answer = 0.11070463773306864
        voigt_idl = 0.11070464
        voigt_check = utils.voigt(a,u)
        self.assertLess(np.abs(voigt_check-voigt_answer)/voigt_answer,1.e-4)
        self.assertLess(np.abs(voigt_idl-voigt_answer)/voigt_answer,1.e-4)

    def test_voigt_2(self):
        a = 10.
        u = 3.
        voigt_answer = 0.05160191683088553
        voigt_idl = 0.051601917
        voigt_check = utils.voigt(a,u)
        self.assertLess(np.abs(voigt_check-voigt_answer)/voigt_answer,1.e-4)
        self.assertLess(np.abs(voigt_idl-voigt_answer)/voigt_answer,1.e-4)

    def test_voigt_3(self):
        a = 15.
        u=-2.
        voigt_answer = 0.03688101492949091
        voigt_idl = 0.036881015
        voigt_check = utils.voigt(a,u)
        self.assertLess(np.abs(voigt_check-voigt_answer)/voigt_answer,1.e-4)
        self.assertLess(np.abs(voigt_idl-voigt_answer)/voigt_answer,1.e-4)

    def test_rrl_freq_1(self):
        """
        Test the rrl_freq function
        Equation from Roelfsema and Goss 1992
        """
        # case 1
        fit_lines = np.array([58,102])
        deltan = 1
        rrl_freq_answer = np.array([32.85219765,6.1068555]) # GHz
        rrl_freq_idl = np.array([32.852149429526861,6.1068408043009752]) # GHz
        rrl_freq_check = utils.rrl_freq(fit_lines,deltan=deltan)
        for answer,check in zip(rrl_freq_answer,rrl_freq_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(rrl_freq_idl,rrl_freq_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_rrl_freq_2(self):
        fit_lines = np.array([65,112])
        deltan = 3
        rrl_freq_answer = np.array([67.15323404,13.49751944]) # GHz
        rrl_freq_idl = np.array([67.153215100231193,13.497524174590783]) # GHz
        rrl_freq_check = utils.rrl_freq(fit_lines,deltan=deltan)
        for answer,check in zip(rrl_freq_answer,rrl_freq_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(rrl_freq_idl,rrl_freq_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)
            
    def test_doppler_width_1(self):
        """
        Test the doppler_width function
        Equation 13 from Viner, Vallee, and Hughes 1979
        """
        electron_temp = 5000. # K 
        v_turbulent = 0. # km/s
        freq = np.array([10.,20.,30.]) # GHz
        deltaV_doppler_answer = 15.12365432994584 # km/s
        deltaV_doppler_idl = 15.123654 # km/s
        deltaNu_doppler_answer = np.array([0.00050447,0.00100894,0.00151341]) # GHz
        deltaNu_doppler_idl = np.array([0.00050447080209084809,0.0010089416041816962,0.0015134124062725443]) # GHz
        deltaV_doppler_check,deltaNu_doppler_check = utils.doppler_width(electron_temp,v_turbulent,freq)
        self.assertLess(np.abs(deltaV_doppler_check-deltaV_doppler_answer)/deltaV_doppler_answer,1.e-4)
        for answer,check in zip(deltaNu_doppler_answer,deltaNu_doppler_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        self.assertLess(np.abs(deltaV_doppler_idl-deltaV_doppler_answer)/deltaV_doppler_answer,1.e-4)
        for idl,answer in zip(deltaNu_doppler_idl,deltaNu_doppler_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)
        
    def test_doppler_width_2(self):
        electron_temp = 10000. # K
        v_turbulent = 15. # km/s
        freq = np.array([15.,25.,35.]) # GHz
        deltaV_doppler_answer = 29.552301727587636 # km/s
        deltaV_doppler_idl = 29.552301 # km/s
        deltaNu_doppler_answer = np.array([0.00147864,0.0024644,0.00345016]) # GHz
        deltaNu_doppler_idl = np.array([0.0014786379828836769,0.0024643967399351446,0.0034501553951909293]) # GHz
        deltaV_doppler_check,deltaNu_doppler_check = utils.doppler_width(electron_temp,v_turbulent,freq)
        self.assertLess(np.abs(deltaV_doppler_check-deltaV_doppler_answer)/deltaV_doppler_answer,1.e-4)
        for answer,check in zip(deltaNu_doppler_answer,deltaNu_doppler_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        self.assertLess(np.abs(deltaV_doppler_idl-deltaV_doppler_answer)/deltaV_doppler_answer,1.e-4)
        for idl,answer in zip(deltaNu_doppler_idl,deltaNu_doppler_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_tau_c_1(self):
        """
        Test tau_c function
        Equations from Viner, Vallee, and Hughes 1979 (eqs. 9,10)
        """
        electron_temp = 5000. # K
        electron_dens = 100. # cm-3
        size = 10. # pc
        freq_obs = np.array([10.,20.,30.]) # GHz
        heII_abundance = 0.
        heIII_abundance = 0.
        tau_c_answer = np.array([6.36817617e-04,1.44428600e-04,6.03490261e-05])
        tau_c_idl = np.array([0.00063682221000715272,0.00014442964131460524,6.0349463698091224e-05])
        tau_c_check = utils.tau_c(electron_temp,electron_dens,size,freq_obs,
                                       heII_abundance=heII_abundance,heIII_abundance=heIII_abundance)
        for answer,check in zip(tau_c_answer,tau_c_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(tau_c_idl,tau_c_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_tau_c_2(self):
        electron_temp = 10000. # K
        electron_dens = 10. # cm-3
        size = 1. # pc
        freq_obs = np.array([15.,25.,35.]) # GHz
        heII_abundance = 0.08
        heIII_abundance = 0.01
        tau_c_answer = np.array([1.10201330e-07,3.71637032e-08,1.81179668e-08])
        # N.B. Found bug in AAK IDL code. Missing * on line 164. This is the fixed answer:
        tau_c_idl = np.array([1.1020132913129986e-07,3.7163705904493164e-08,1.8117965318156966e-08])
        tau_c_check = utils.tau_c(electron_temp,electron_dens,size,freq_obs,
                                       heII_abundance=heII_abundance,heIII_abundance=heIII_abundance)
        for answer,check in zip(tau_c_answer,tau_c_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(tau_c_idl,tau_c_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_tau_l_1(self):
        """
        Test tau_l function
        Equations from Viner, Vallee, and Hughes 1979
        """
        electron_temp = 10000. # K
        electron_dens = 150. # cm-3
        size = 50. # pc
        fit_lines = np.array([58.,102.])
        deltan = 1.
        v_turbulent = 15. # km/s
        use_voigt = True
        heII_abundance = 0.
        heIII_abundance = 0.
        bn = 1.
        betan = 1.
        # calculate Doppler velocity (frequency doesn't matter)
        deltaV_doppler, deltaNu_doppler = utils.doppler_width(electron_temp,v_turbulent,10.)
        tau_l_answer = np.array([6.69337456e-05,3.57980909e-04])
        tau_l_idl = np.array([6.6934236797953418e-05,0.00035798352477698156])
        tau_l_check = utils.tau_l(electron_temp,electron_dens,size,fit_lines,deltaV_doppler,
                                       deltan=deltan,v_turbulent=v_turbulent,
                                       use_voigt=use_voigt,heII_abundance=heII_abundance,
                                       heIII_abundance=heIII_abundance,bn=bn,betan=betan)
        for answer,check in zip(tau_l_answer,tau_l_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(tau_l_idl,tau_l_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_tau_l_2(self):
        electron_temp = 5000. # K
        electron_dens = 100. # cm-3
        size = 100. # pc
        fit_lines = np.array([65,115])
        deltan = 1
        v_turbulent = 0. # km/s
        use_voigt = True
        heII_abundance = 0.08
        heIII_abundance = 0.01
        bn = 0.9
        betan = -45.
        # calculate Doppler velocity (frequency doesn't matter)
        deltaV_doppler, deltaNu_doppler = utils.doppler_width(electron_temp,v_turbulent,10.)
        tau_l_answer = np.array([-0.03408153,-0.18429481])
        tau_l_idl = np.array([-0.034081771945503062,-0.18429615935433580])
        tau_l_check = utils.tau_l(electron_temp,electron_dens,size,fit_lines,deltaV_doppler,
                                       deltan=deltan,v_turbulent=v_turbulent,
                                       use_voigt=use_voigt,heII_abundance=heII_abundance,
                                       heIII_abundance=heIII_abundance,bn=bn,betan=betan)
        for answer,check in zip(tau_l_answer,tau_l_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(tau_l_idl,tau_l_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_tau_l_3(self):
        electron_temp = 5000. # K
        electron_dens = 100. # cm-3
        size = 100. # pc
        fit_lines = np.array([65.,115.])
        deltan = 1
        v_turbulent = 0. # km/s
        use_voigt = False
        heII_abundance = 0.08
        heIII_abundance = 0.01
        bn = 0.9
        betan = -45.
        # calculate Doppler velocity (frequency doesn't matter)
        deltaV_doppler, deltaNu_doppler = utils.doppler_width(electron_temp,v_turbulent,10.)
        tau_l_answer = np.array([-0.03408548,-0.18597207])
        tau_l_idl = np.array([-0.034085729216465736,-0.18597342453371252])
        tau_l_check = utils.tau_l(electron_temp,electron_dens,size,fit_lines,deltaV_doppler,
                                       deltan=deltan,v_turbulent=v_turbulent,
                                       use_voigt=use_voigt,heII_abundance=heII_abundance,
                                       heIII_abundance=heIII_abundance,bn=bn,betan=betan)
        for answer,check in zip(tau_l_answer,tau_l_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(tau_l_idl,tau_l_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_peak_line_fluxden_1(self):
        """
        Test peak_line_fluxden function
        Equation 1 from Anantharamaiah+1993
        """
        electron_temp = 5000. # K
        electron_dens = 100. # cm-3
        HII_region_size = 1 # pc
        distance = 3. # Mpc
        omega_HII_region = np.pi * (0.5*HII_region_size/(distance*1.e6))**2.
        omega_region = 60. # sq arcsec
        omega_region = omega_region * (np.pi/(180.*60.*60.))**2. # sr
        fit_lines = np.array([58,102])
        rrl_freq = utils.rrl_freq(fit_lines)
        heII_abundance = 0.
        heIII_abundance = 0.
        tau_c = utils.tau_c(electron_temp,electron_dens,HII_region_size,
                                 rrl_freq,heII_abundance=heII_abundance,heIII_abundance=heIII_abundance)
        v_turbulent=15.
        bn = 1.
        betan = 1.
        deltaV_doppler, deltaNu_doppler = utils.doppler_width(electron_temp,v_turbulent,rrl_freq)
        tau_l_star = utils.tau_l(electron_temp,electron_dens,HII_region_size,
                                      fit_lines,deltaV_doppler,deltan=1,v_turbulent=v_turbulent,
                                      heII_abundance=heII_abundance,heIII_abundance=heIII_abundance,
                                      bn=1.,betan=1.)
        tau_l = utils.tau_l(electron_temp,electron_dens,HII_region_size,
                                 fit_lines,deltaV_doppler,deltan=1,v_turbulent=v_turbulent,
                                 heII_abundance=heII_abundance,heIII_abundance=heIII_abundance,
                                 bn=bn,betan=betan)
        background_fluxden = 0.
        peak_line_fluxden_HII_region_answer = np.array([5.69481083e-34,1.05179646e-34]) # W/m2/Hz
        peak_line_fluxden_HII_region_idl = np.array([5.6917396965089532e-31,1.0519186812853735e-31])*1.e-3 # W/m2/Hz
        peak_line_fluxden_background_answer = np.array([0.,0.]) # W/m2/Hz
        peak_line_fluxden_background_idl = np.array([0.,0.])*1.e-3 # W/m2/Hz
        peak_line_fluxden_HII_region_check,peak_line_fluxden_background_check = \
          utils.peak_line_fluxden(electron_temp,omega_HII_region,
                                       omega_region,rrl_freq,tau_c,
                                       tau_l_star,tau_l,bn=bn,
                                       background_fluxden=background_fluxden)
        for answer,check in zip(peak_line_fluxden_HII_region_answer,peak_line_fluxden_HII_region_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for answer,check in zip(peak_line_fluxden_background_answer,peak_line_fluxden_background_check):
            self.assertEqual(check,answer)
        for idl,answer in zip(peak_line_fluxden_HII_region_idl,peak_line_fluxden_HII_region_answer):
            # difference is 0.000539. I think this is a physical constants problem
            self.assertLess(np.abs(idl-answer)/answer,1.e-3)
        for idl,answer in zip(peak_line_fluxden_background_idl,peak_line_fluxden_background_answer):
            self.assertEqual(idl,answer)

    def test_peak_line_fluxden_2(self):
        electron_temp = 10000. # K
        electron_dens = 10. # cm-3
        HII_region_size = 100. # pc
        distance = 3. # Mpc
        omega_HII_region = np.pi * (0.5*HII_region_size/(distance*1.e6))**2.
        omega_region = 30. # sq arcsec
        omega_region = omega_region * (np.pi/(180.*60.*60.))**2. # sr
        fit_lines = np.array([80.,115])
        rrl_freq = utils.rrl_freq(fit_lines)
        heII_abundance = 0.08
        heIII_abundance = 0.01
        tau_c = utils.tau_c(electron_temp,electron_dens,HII_region_size,
                                 rrl_freq,heII_abundance=heII_abundance,heIII_abundance=heIII_abundance)
        v_turbulent=0.
        bn = 0.9
        betan = -45.
        deltaV_doppler, deltaNu_doppler = utils.doppler_width(electron_temp,v_turbulent,rrl_freq)
        tau_l_star = utils.tau_l(electron_temp,electron_dens,HII_region_size,
                                      fit_lines,deltaV_doppler,deltan=1,v_turbulent=v_turbulent,
                                      heII_abundance=heII_abundance,heIII_abundance=heIII_abundance,
                                      bn=1.,betan=1.)
        tau_l = utils.tau_l(electron_temp,electron_dens,HII_region_size,
                                 fit_lines,deltaV_doppler,deltan=1,v_turbulent=v_turbulent,
                                 heII_abundance=heII_abundance,heIII_abundance=heIII_abundance,
                                 bn=bn,betan=betan)
        background_fluxden = 30. * 1.e-29 # W/m2/Hz
        peak_line_fluxden_HII_region_answer = np.array([7.45601155e-31,2.52914885e-31]) # W/m2/Hz
        peak_line_fluxden_HII_region_idl = np.array([7.4560553782987313e-28,2.5291638366348540e-28])*1.e-3 # W/m2/Hz
        peak_line_fluxden_background_answer = np.array([1.46109857e-32,4.31011918e-32]) # W/m2/Hz
        peak_line_fluxden_background_idl = np.array([1.4605445e-32,4.3102095e-32]) # W/m2/Hz
        peak_line_fluxden_HII_region_check,peak_line_fluxden_background_check = \
          utils.peak_line_fluxden(electron_temp,omega_HII_region,
                                       omega_region,rrl_freq,tau_c,
                                       tau_l_star,tau_l,bn=bn,
                                       background_fluxden=background_fluxden) 
        for answer,check in zip(peak_line_fluxden_HII_region_answer,peak_line_fluxden_HII_region_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for answer,check in zip(peak_line_fluxden_background_answer,peak_line_fluxden_background_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(peak_line_fluxden_HII_region_idl,peak_line_fluxden_HII_region_answer):            
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)
        for idl,answer in zip(peak_line_fluxden_background_idl,peak_line_fluxden_background_answer):
            # difference is 0.000379. I think this is a physical constants problem
            self.assertLess(np.abs(idl-answer)/answer,1.e-3)

    def test_line_flux_1(self):
        """
        Test line_flux function
        Integral of Gaussian with peak_line_fluxden amplitude and
        deltaNu_doppler width
        """
        peak_line_fluxden = np.array([20.0*1.e-29,100.*1.e-29]) # W/m2/Hz
        deltaNu_doppler = np.array([0.0005,0.005]) # GHz
        line_flux_answer = np.array([1.06446702e-22,5.32233510e-21]) # W/m2
        line_flux_idl = np.array([1.0644671e-22,5.3223353e-21]) # W/m2
        line_flux_check = utils.line_flux(peak_line_fluxden,deltaNu_doppler)
        for answer,check in zip(line_flux_answer,line_flux_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(line_flux_idl,line_flux_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_min_num_HII_regions_1(self):
        """
        Test min_num_HII_regions function
        Equations on page 592 of Anantharamaiah+1993
        """
        measured_line_fwhm = 50. # km/s
        deltaV_doppler = 35. # km/s
        omega_region = 30. # sq arcsec
        omega_region = omega_region * (np.pi/(180.*60.*60.))**2. # sr
        measured_omegaB = np.pi*0.1**2./4. * (np.pi/(180.*60.*60.))**2. # sr
        min_num_HII_regions_answer = 5456.74090600784
        min_num_HII_regions_idl = 5456.7412
        min_num_HII_regions_check = utils.min_num_HII_regions(measured_line_fwhm,deltaV_doppler,
                                                                   measured_omegaB,omega_region)
        self.assertLess(np.abs(min_num_HII_regions_check-min_num_HII_regions_answer)/min_num_HII_regions_answer,1.e-4)
        self.assertLess(np.abs(min_num_HII_regions_idl-min_num_HII_regions_answer)/min_num_HII_regions_answer,1.e-4)

    def test_min_num_HII_regions_2(self):
        measured_line_fwhm = 80. # km/s
        deltaV_doppler = 30. # km/s
        omega_region = 50. # sq arcsec
        omega_region = omega_region * (np.pi/(180.*60.*60.))**2. # sr
        measured_omegaB = np.pi*1.0**2./4. * (np.pi/(180.*60.*60.))**2. # sr
        min_num_HII_regions_answer = 169.76527263135503
        min_num_HII_regions_idl = 169.76526
        min_num_HII_regions_check = utils.min_num_HII_regions(measured_line_fwhm,deltaV_doppler,
                                                                   measured_omegaB,omega_region)
        self.assertLess(np.abs(min_num_HII_regions_check-min_num_HII_regions_answer)/min_num_HII_regions_answer,1.e-4)
        self.assertLess(np.abs(min_num_HII_regions_idl-min_num_HII_regions_answer)/min_num_HII_regions_answer,1.e-4)

    def test_num_HII_regions_1(self):
        """
        Test num_HII_regions function
        Equations on page 592 of Anantharamaiah+1993
        """
        measured_line_flux = 100.*1.e-30 # W/m2
        line_flux = 1.5e-29 # W/m2
        line_flux_override = None
        num_HII_regions_answer = 6.666666666666667
        num_HII_regions_idl = 6.6666665
        num_HII_regions_check = utils.num_HII_regions(measured_line_flux,line_flux,
                                                           line_flux_override=line_flux_override)
        self.assertLess(np.abs(num_HII_regions_check-num_HII_regions_answer)/num_HII_regions_answer,1.e-4)
        self.assertLess(np.abs(num_HII_regions_idl-num_HII_regions_answer)/num_HII_regions_answer,1.e-4)

    def test_num_HII_regions_2(self):
        measured_line_flux = 50.*1.e-30 # W/m2
        line_flux = 1.e-30 # W/m2
        line_flux_override = 8.e-30 # W/m2
        num_HII_regions_answer = 8.0
        num_HII_regions_idl = 8.0000000
        num_HII_regions_check = utils.num_HII_regions(measured_line_flux,line_flux,
                                                           line_flux_override=line_flux_override)
        self.assertLess(np.abs(num_HII_regions_check-num_HII_regions_answer)/num_HII_regions_answer,1.e-4)
        self.assertLess(np.abs(num_HII_regions_idl-num_HII_regions_answer)/num_HII_regions_answer,1.e-4)

    def test_num_HII_regions_los_1(self):
        """
        Test num_HII_regions function
        Equations on page 592 of Anantharamaiah+1993
        """
        num_HII_regions = 100.
        HII_region_size = 1. # pc
        region_size = 1000. # pc
        num_HII_regions_los_answer = 0.0001
        num_HII_regions_los_idl = 0.00010000001
        num_HII_regions_los_check = utils.num_HII_regions_los(num_HII_regions,
                                                                   HII_region_size,region_size)
        self.assertLess(np.abs(num_HII_regions_los_check-num_HII_regions_los_answer)/num_HII_regions_los_answer,1.e-4)
        self.assertLess(np.abs(num_HII_regions_los_idl-num_HII_regions_los_answer)/num_HII_regions_los_answer,1.e-4)
        
    def test_filling_factor_1(self):
        """
        Test filling_factor function
        Equations on page 592 of Anantharamaiah+1993
        """
        num_HII_regions = 100.
        HII_region_size = 1. # pc
        region_size = 1000. # pc
        filling_factor_answer = 1e-07
        filling_factor_idl = 1.0000002e-07
        filling_factor_check = utils.filling_factor(num_HII_regions,
                                                         HII_region_size,region_size)
        self.assertLess(np.abs(filling_factor_check-filling_factor_answer)/filling_factor_answer,1.e-4)
        self.assertLess(np.abs(filling_factor_idl-filling_factor_answer)/filling_factor_answer,1.e-4)

    def test_thermal_fluxden_1(self):
        """
        Test thermal_fluxden function
        Equations 2-4 in Anantharamaiah+1993
        """
        contfreqs = np.array([10,40]) # GHz
        electron_temp = 5000. # K
        electron_dens = 100. # cm-3
        HII_region_size = 1 # pc
        region_size = 1000. # pc
        distance = 3. # Mpc
        omega_HII_region = np.pi * (0.5*HII_region_size/(distance*1.e6))**2.
        omega_region = np.pi * (0.5*region_size/(distance*1.e6))**2.
        heII_abundance = 0.
        heIII_abundance = 0.
        tau_c_cont = utils.tau_c(electron_temp,electron_dens,HII_region_size,contfreqs,
                                       heII_abundance=heII_abundance,heIII_abundance=heIII_abundance)
        num_HII_regions = 100.
        num_HII_regions_los = utils.num_HII_regions_los(num_HII_regions,HII_region_size,region_size)
        thermal_fluxden_answer = np.array([8.53670801e-32,6.95232971e-32]) # W/m2/Hz
        thermal_fluxden_idl = np.array([8.5338193031806409e-29,6.9037639306629900e-29])*1.e-3 # W/m2/Hz        
        thermal_fluxden_check = utils.thermal_fluxden(contfreqs,electron_temp,
                                                           omega_HII_region,omega_region,
                                                           tau_c_cont,num_HII_regions,num_HII_regions_los)
        for answer,check in zip(thermal_fluxden_answer,thermal_fluxden_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(thermal_fluxden_idl,thermal_fluxden_answer):
            # difference is 0.000338 and 0.006985. Physical constants?
            self.assertLess(np.abs(idl-answer)/answer,1.e-2)

    def test_thermal_fluxden_2(self):
        contfreqs = np.array([10,40]) # GHz
        electron_temp = 5000. # K
        electron_dens = 1.e4 # cm-3
        HII_region_size = 10. # pc
        region_size = 1000. # pc
        distance = 3. # Mpc
        omega_HII_region = np.pi * (0.5*HII_region_size/(distance*1.e6))**2.
        omega_region = np.pi * (0.5*region_size/(distance*1.e6))**2.
        heII_abundance = 0.
        heIII_abundance = 0.
        tau_c_cont = utils.tau_c(electron_temp,electron_dens,HII_region_size,contfreqs,
                                       heII_abundance=heII_abundance,heIII_abundance=heIII_abundance)
        num_HII_regions = 1000.
        num_HII_regions_los = utils.num_HII_regions_los(num_HII_regions,HII_region_size,region_size)
        thermal_fluxden_answer = np.array([6.31444022e-24,5.9380819514121456e-24]) # W/m2/Hz
        thermal_fluxden_idl = np.array([6.3144821638610799e-21,5.9381215189383294e-21])*1.e-3 # W/m2/Hz
        thermal_fluxden_check = utils.thermal_fluxden(contfreqs,electron_temp,
                                                           omega_HII_region,omega_region,
                                                           tau_c_cont,num_HII_regions,num_HII_regions_los)
        for answer,check in zip(thermal_fluxden_answer,thermal_fluxden_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(thermal_fluxden_idl,thermal_fluxden_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_nonthermal_fluxden_1(self):
        """
        Test nonthermal_fluxden function
        Equations 5 and 6 in Anantharamaiah+1993
        """
        contfreqs = np.array([5.,35.]) # GHz
        contfluxden = np.array([70.,30.]) # mJy
        electron_temp = 5000. # K
        electron_dens = 1.e3 # cm-3
        HII_region_size = 1. # pc
        region_size = 500. # pc
        distance = 3. # Mpc
        omega_HII_region = np.pi * (0.5*HII_region_size/(distance*1.e6))**2.
        omega_region = np.pi * (0.5*region_size/(distance*1.e6))**2.
        heII_abundance = 0.
        heIII_abundance = 0.
        tau_c_cont = utils.tau_c(electron_temp,electron_dens,HII_region_size,contfreqs,
                                       heII_abundance=heII_abundance,heIII_abundance=heIII_abundance)
        num_HII_regions = 1000.
        num_HII_regions_los = utils.num_HII_regions_los(num_HII_regions,HII_region_size,region_size)
        thermal_fluxden = utils.thermal_fluxden(contfreqs,electron_temp,
                                                omega_HII_region,omega_region,
                                                tau_c_cont,num_HII_regions,num_HII_regions_los)
        nonthermal_fluxden_answer = np.array([6.07993533e-28,2.28965594e-28]) # W/m2/Hz
        nonthermal_fluxden_idl = np.array([6.0799355e-28,2.2896558e-28]) # W/m2/Hz
        unatten_nonthermal_fluxden_answer = np.array([6.08026917e-28,2.28965792e-28]) # W/m2/Hz
        unatten_nonthermal_fluxden_idl = np.array([6.0802692e-28,2.2896579e-28]) # W/m2/Hz
        nonthermal_fluxden_check,unatten_nonthermal_fluxden_check = \
          utils.nonthermal_fluxden(contfluxden,thermal_fluxden,tau_c_cont,
                                   num_HII_regions_los)
        for answer,check in zip(nonthermal_fluxden_answer,nonthermal_fluxden_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for answer,check in zip(unatten_nonthermal_fluxden_answer,unatten_nonthermal_fluxden_check):
            self.assertLess(np.abs(check-answer)/answer,1.e-4)
        for idl,answer in zip(nonthermal_fluxden_idl,nonthermal_fluxden_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)
        for idl,answer in zip(unatten_nonthermal_fluxden_idl,unatten_nonthermal_fluxden_answer):
            self.assertLess(np.abs(idl-answer)/answer,1.e-4)

    def test_spectral_info_1(self):
        """
        test spectral info function
        """
        freq = np.array([5.,10.,35.]) # GHz
        fluxden = np.array([6.e-28,5.e-28,2.e-28]) # W/m2/Hz
        spectral_index_answer = -0.56457503405357967
        spectral_index_idl = -0.56457502
        scaling_factor_answer = 1.4885802903369804e-27
        scaling_factor_idl = 1.4885801e-27
        spectral_index_check,scaling_factor_check = utils.spectral_info(freq,fluxden)
        self.assertLess(np.abs(spectral_index_check-spectral_index_answer)/spectral_index_answer,1.e-4)
        self.assertLess(np.abs(scaling_factor_check-scaling_factor_answer)/scaling_factor_answer,1.e-4)
        self.assertLess(np.abs(spectral_index_idl-spectral_index_answer)/spectral_index_answer,1.e-4)
        self.assertLess(np.abs(scaling_factor_idl-scaling_factor_answer)/scaling_factor_answer,1.e-4)
        
    def test_spectral_info_2(self):
        """
        test spectral_info function
        """
        freq = np.array([35.,10.,5.]) # GHz
        fluxden = np.array([2.e-28,5.e-28,6.e-28]) # W/m2/Hz
        spectral_index_answer = -0.56457503405357967
        spectral_index_idl = -0.56457502
        scaling_factor_answer = 1.4885802903369804e-27
        scaling_factor_idl = 1.4885801e-27
        spectral_index_check,scaling_factor_check = utils.spectral_info(freq,fluxden)
        self.assertLess(np.abs(spectral_index_check-spectral_index_answer)/spectral_index_answer,1.e-4)
        self.assertLess(np.abs(scaling_factor_check-scaling_factor_answer)/scaling_factor_answer,1.e-4)
        self.assertLess(np.abs(spectral_index_idl-spectral_index_answer)/spectral_index_answer,1.e-4)
        self.assertLess(np.abs(scaling_factor_idl-scaling_factor_answer)/scaling_factor_answer,1.e-4)

if __name__ == "__main__":
    unittest.main()
