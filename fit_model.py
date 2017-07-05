"""
fit_model.py
Code to fit an HII region model
Based on IDL code by A.A. Kepley

Trey Wenger May 2016

dsb 21Jun2016 - Add plot_rss function.
dsb 23Jun2016 - Execute plot_model.py in fit_model procedure.
dsb 29Jun2016 - Add plotRRS and plotModel logic to fit_model.
dsb 14Jul2016 - Modify best model print statement.
dsb 30Sep2016 - Add mass_ion,num_ion,sfr to print statement.
"""
import math
import numpy as np
from numpy.lib.recfunctions import append_fields
import matplotlib
matplotlib.rc_file('matplotlibrc')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import itertools
import pickle
import time
import utils
from hii_region_model import HIIRegionModel
import argparse
import os

import pdb
#pdb.set_trace()


def read_line_data(linefile):
    """
    Read line data from linefile

    Inputs:
      linefile = filename of file containing line data
      Format: (space-delimited and ; indicate comment lines)
        ; freq nn S_L_peak error FWHM error V_0  error S_int    error    S_int    error    Bmaj   Bmin
        ; GHz     mJy      mJy   km/s km/s  km/s km/s  mJy.km/s mJy.km/s W m^{-2} W m^{-2} arcsec arcsec
      freq     = observed frequency
      nn       = line lower energy level
      S_L_peak = peak flux density
      FWHM     = line full width at half maximum
      V0       = line center velocity
      S_int    = line flux
      Bmaj     = beam major axis
      Bmin     = beam minor axis

    Returns:
      line_data = numpy array containing data
    """
    #
    # Read line data
    #
    line_data = np.genfromtxt(linefile,comments=';',dtype=None,
                              names='freq_GHz,n,peak_mJy,peak_mJy_err,fwhm_kms,fwhm_kms_err,v0_kms,v0_kms_err,intflux_mJykms,intflux_mJykms_err,intflux_Wm2,intflux_Wm2_err,Bmaj_arcsec,Bmin_arcsec')
    #
    # Calculate beam solid angle in sterradians
    #
    omegaB = np.pi * (line_data['Bmaj_arcsec']*line_data['Bmin_arcsec'])/4.
    omegaB *= (np.pi/(60.*60.*180.))**2.
    #
    # add omegaB to data structure
    #
    line_data = append_fields(line_data,'omegaB_sr',omegaB,usemask=False)
    return line_data    

def read_cont_data(contfile):
    """
    Read continuum data from contfile

    Inputs:
      contfile = filename of file containing continuum data
      Format: (space delimited and ; indicate comment lines)
        ; freq S_C error
        ; GHz  Jy  Jy   
      freq = observed frequency
      S_C  = continuum flux density

    N.B. Currently no error in contfile? Error is definied by
    AAK as 0.1 * flux density

    Returns:
      cont_data = numpy array containing data
    """
    #
    # Read cont data
    #
    cont_data = np.genfromtxt(contfile,comments=';',dtype=None,
                              usecols=(0,1),names='freq_GHz,S_C_mJy')
    #
    # Convert flux density to mJy
    #
    cont_data['S_C_mJy'] *= 1000.
    #
    # Calculate uncertainty (AAK indicates "temporary fix")
    #
    S_C_mJy_err = cont_data['S_C_mJy']*0.1
    #
    # add fluxdensity error to data structure
    #
    cont_data = append_fields(cont_data,'S_C_mJy_err',S_C_mJy_err,usemask=False)
    return cont_data

def fit_model(linefile,contfile,omega_region,distance,
              electron_temps=[5000.,7000.,10000.],
              electron_dens=[10.,100.,1000.],
              HII_region_sizes=[0.1,1.,10.],
              num_HII_regions=[None],
              v_turbulent=0.,heII_abundance=0.,heIII_abundance=0.,
              n_for_num=58,fit_lines=[],background_fluxden=0.,
              use_voigt=False,line_flux_override=None,
              outfile='best_model.pkl',verbose=False,
              plotRSS=False,plotModel=False):
    """
    Fit the best HII region model to the data in linefile and contfile
    using the above parameters and grid values.

    Inputs:
      linefile  = filename where line data is saved (see read_line_data)
      contfile  = filename where continuum data is saved (see read_cont_data)
      outfile   = file in which to save best fit model

    The following input parameters define properties that are the
    same for all models
      omega_region    = solid angle of CASA region used to extract (arcsec2)
      distance        = distance to object (Mpc)
      v_turbulent     = turbulent velocity (km/s)
      heII_abundance  = abundance of singly ionized helium relative
                        to singly ionized hydrogen.
      heIII_abundance = abundance of doubly ionized helium relative
                        to singly ionized hydrogen.
      background_fluxden = background continuum flux density (mJy)
      use_voigt          = if True, use a Voigt profile, otherwise Gaussian
      fit_lines          = list of lines in linedata to use in model fitting.
                           if empty: use all lines
                           i.e. to fit H56 and H102 use fit_lines=[56,102]
      n_for_num          = line to use for estimating the number of HII regions
      line_flux_override = use this line flux to estimate number of
                           HII regions if set (W/m2)

      verbose            = if True, print out model check info
      plotRSS            = if True, plot RRS color image
      plotModel          = if True, plot cont/line best model with data

    The following four input parameters define the 4-D grid used to
    fit the model to the data:
      electron_temp   = list of electron temperatures (K)
      electron_dens   = list of electron densities (cm-3)
      HII_region_size = list of linear HII region sizes (pc)
      num_HII_regions = list of number of HII regions
                        if [None], estimate the number of HII regions
                        from the data and the model

    Returns:
      best_model = best fit model
    """
    start_time = time.time()
    #
    # read line and continuum data
    #
    linedata = read_line_data(linefile)
    contdata = read_cont_data(contfile)
    #
    # if fit_lines is empty, use all RRLs to calculate residuals
    #
    if len(fit_lines) == 0:
        fit_lines = linedata['n']
    else:
        fit_lines = np.array(fit_lines)
    #
    # Store frequencies used in the fit (GHz)
    #
    fit_freqs = np.zeros(len(fit_lines))
    for i,n in enumerate(fit_lines):
        myrow = np.where(linedata['n'] == n)[0]
        if len(myrow) != 1:
            raise ValueError("Problem with level {0} in linedata.".format(n))
        fit_freqs[i] = linedata['freq_GHz'][myrow[0]]
    #
    # Get all rss values and models
    #
    fit_rss = []
    fit_models = []
    #
    # Loop over grid values, create models, and calculate chisq
    #
    num_models = len(electron_temps)*len(electron_dens)*len(HII_region_sizes)*len(num_HII_regions)
    num_failed = 0
    num_passed = 0
    for temp,dens,size,num in \
      itertools.product(electron_temps,electron_dens,HII_region_sizes,num_HII_regions):
        #
        # Set up model
        #
        model = HIIRegionModel(linedata,contdata,omega_region,distance,
                               electron_temp=temp,electron_dens=dens,HII_region_size=size,num_HII_regions=num,
                               v_turbulent=v_turbulent,heII_abundance=heII_abundance,heIII_abundance=heIII_abundance,
                               fit_lines=fit_lines,n_for_num=n_for_num,
                               background_fluxden=background_fluxden,
                               line_flux_override=line_flux_override,use_voigt=use_voigt)
        #
        # Calculate physical properties for this model
        #
        model.calc_HII_region_properties()

        # print info for testing
        print ("=============================================================")
        foo = "(Te,ne,l) = (%10.4e, %10.4e, %10.4e)" % (temp,dens,size)
        print (foo)
        foo = "(Rsun,Omega,Vturb,y) = (%10.4e, %10.4e, %10.4e, %10.4e)" % (distance,omega_region,v_turbulent,heII_abundance)
        print (foo)
        print ("NHII = ", model.num_HII_regions)
        print ("n = ", fit_lines)
        print ("PeakLineFluxHII = ", model.peak_line_fluxden_HII_region)
        print ("FWHMLineFluxHII = ", model.deltaV_doppler)
        print ("LineFluxHII = ", model.line_flux)

        #
        # Check that physical properties are sane
        #
        if verbose: print("Checking this model:")
        if not model.check(verbose=verbose):
            if verbose: print("At least one failed test, skipping this model.")
            num_failed += 1
            continue
        if verbose: print("This model passed all tests!")
        num_passed += 1
        #
        # Calculate the residuals
        #
        rss = model.rss()
        fit_rss.append(rss)
        fit_models.append(model)


    #
    # Plot RSS image
    #
    if plotRSS:
        plotRSS = plot_rss(electron_temps, electron_dens, HII_region_sizes, fit_models, fit_rss)

    #
    # Calculate runtime
    #
    time_diff = time.time() - start_time
    total_hours = int(time_diff / 3600.)
    total_mins = int((time_diff-total_hours*3600.) / 60.)
    total_secs = time_diff-total_hours*3600.-total_mins*60.
    print(" ")
    print("Ran {0} models in {1:02d}h {2:02d}m {3:02.2f}s".format(num_models,total_hours,total_mins,total_secs))
    per_model = time_diff / num_models
    hours_per = int(per_model / 3600.)
    mins_per = int((per_model-hours_per*3600.) / 60.)
    secs_per = per_model-hours_per*3600.-mins_per*60.
    print("{0:02d}h {1:02d}m {2:02.2f}s per model".format(hours_per,mins_per,secs_per))
    print("{0} models passed, {1} models failed".format(num_passed,num_failed))
    # if no models pass then write last model and exit
    if num_passed == 0:
        failed_model = model
        with open('failed_model.pkl','wb') as f:
            pickle.dump(failed_model,f)
        if plotModel:
            plots = os.system("python plot_model.py failed_model.pkl")
        return None

    #
    # Find best model (smallest rss)
    #
    best_rss = np.min(fit_rss)
    ind = np.argmin(fit_rss)
    best_model = fit_models[ind]
    # print results
    print("================================")
    print("Multile Compact HII Region Model")
    print("================================")
    print("Best model RSS: {:g}".format(best_rss))
    print("Electron temperature = {:g} K".format(best_model.electron_temp))
    print("Electron density = {:g} cm-3".format(best_model.electron_dens))
    print("HII region size = {:g} pc".format(best_model.HII_region_size))
    print("Number of HII_regions = {:g}".format(best_model.num_HII_regions))
    print("Non-thermal spectral index = {:g}".format(best_model.nonthermal_spectral_index))
    print("Mass of ionized gas = {:g} Msun".format(best_model.mass_ion))
    print("Number of H-ionizing photons = {:g} s-1".format(best_model.num_ion))
    print("Star Formation Rate = {:g} Msun/yr".format(best_model.sfr))
    print("================================")
    with open(outfile,'wb') as f:
        pickle.dump(best_model,f)

    #
    # Make model plots
    #
    if plotModel:
        plots = os.system("python plot_model.py best_model.pkl")


    return best_model


def plot_rss(electron_temps, electron_dens, HII_region_sizes, fit_models, fit_rss):
    '''
    Plot the fit results 

    Inputs:
      electron_temps: range of model Te values
      electron_dens: range of model ne values
      HII_region_sizes: range of HII region sizes
      fit_models: model parameters for each fit

    Returns:
      none
    '''
    # calculate mean rss
    rss_mean = np.array(fit_rss).mean()
    # create the grid for the colormap
    x = np.log10(electron_dens)
    y = np.log10(HII_region_sizes)
    X, Y = np.meshgrid(x, y)
    # create the rss array for the color map
    Z = np.zeros(len(x)*len(y)).reshape(len(x),len(y))
    # Set the rss values for when the model converged; else set to NaN
    for i in range(len(x)):
        for j in range(len(y)):
            Z[i,j] = float('NaN')
            for k in range(len(fit_rss)):
                #pdb.set_trace()
                if (X[i,j] == math.log10(fit_models[k].electron_dens)) and (Y[i,j] == math.log10(fit_models[k].HII_region_size)):
                    # normalize by the mean
                    Z[i,j] = fit_rss[k]/rss_mean
    # generate the image (color maps: copper, gist_ncar, YlOrRd)
    #im = plt.imshow(Z, interpolation='None', cmap=plt.get_cmap('copper'),
    #                origin='lower', extent=[x[0], x[len(x)-1], y[0], y[len(y)-1]])
    # find limits
    deltax = (x[1]-x[0])/2.0
    deltay = (y[1]-y[0])/2.0
    xmin = x[0] - deltax
    xmax = x[len(x)-1] + deltax
    ymin = y[0] - deltay
    ymax = y[len(y)-1] + deltay
    # plot
    im = plt.imshow(Z, interpolation='none', cmap=plt.get_cmap('copper'),
                    origin='lower', extent=[xmin,xmax,ymin,ymax])
    plt.xlabel(r'Log($n_{e}$) $cm^{-3}$', fontsize=20)
    plt.ylabel(r'Log($\ell$) pc', fontsize=20)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ctit = "Normalized RSS (Te = %5.1f K)" % (electron_temps[0])
    plt.title(ctit)
    plt.colorbar(im)
    plt.savefig('rss.eps')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fit HII region model to observed line and continuum data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--linefile',required=True,help='file containing line data')
    parser.add_argument('--contfile',required=True,help='file containing continuum data')
    parser.add_argument('--omega_region',required=True,type=float,help='solid angle of region used to extract line and continuum data (sq arcsec)')
    parser.add_argument('--distance',required=True,type=float,help='distance to object (Mpc)')
    parser.add_argument('--electron_temps',nargs='+',default=[5000.,7000.,10000.],type=float,help='model electron temperatures to try in fit (K)')
    parser.add_argument('--electron_dens',nargs='+',default=[10.,100.,1000.],type=float,help='model electron densities to try in fit (cm-3)')
    parser.add_argument('--HII_region_sizes',nargs='+',default=[0.1,1.,10.],type=float,help='model HII region sizes to try in fit (pc)')
    parser.add_argument('--num_HII_regions',nargs='+',default=[None],type=float,help='model number of HII regions to try in fit. If not supplied, number of HII regions is estimated from the data')
    parser.add_argument('--v_turbulent',default=0.,type=float,help='model turbulent velocity (km/s)')
    parser.add_argument('--heII_abundance',default=0.,type=float,help='abundance of singly-ionized helium')
    parser.add_argument('--heIII_abundance',default=0.,type=float,help='abundance of doubly-ionized helium')
    parser.add_argument('--n_for_num',default=58,type=int,help='the observed RRL transition to use in estimating the model number of HII regions (i.e. if num_HII_regions is [None])')
    parser.add_argument('--fit_lines',nargs='+',default=[],type=int,help='the observed RRL transitions to use in fitting the model. If empty, use all observed RRL transitions')
    parser.add_argument('--background_fluxden',default=0.,type=float,help='model background continuum flux density (mJy)')
    parser.add_argument('--use_voigt',action='store_true',help='use a voigt profile model instead of Gaussian profile (default)')
    parser.add_argument('--line_flux_override',default=None,type=float,help='if set, use this RRL flux in estimating the model number of HII regions (W/m2) (i.e. if num_HII_regions is [None]')
    parser.add_argument('--outfile',default='best_model.pkl',help='output pickle file to save best fit model')
    parser.add_argument('--verbose',action='store_true',help='print info on why a model is rejected during model checks')
    args = parser.parse_args()
    fit_model(args.linefile,args.contfile,args.omega_region,args.distance,
              electron_temps=args.electron_temps,
              electron_dens=args.electron_dens,
              HII_region_sizes=args.HII_region_sizes,
              num_HII_regions=args.num_HII_regions,
              v_turbulent=args.v_turbulent,
              heII_abundance=args.heII_abundance,
              heIII_abundance=args.heIII_abundance,
              n_for_num=args.n_for_num,fit_lines=args.fit_lines,
              background_fluxden=args.background_fluxden,
              use_voigt=args.use_voigt,
              line_flux_override=args.line_flux_override,
              outfile=args.outfile,verbose=args.verbose)
