"""
plot_model.py
Code to plot an HII region model fit and data
Based on IDL code by A.A. Kepley

Trey Wenger June 2016

dsb 21Jun2016 - Use linecolor to distinguish models instead of linestyle
dsb 02Sep2016 - Modify units labels in plots
dsb 14Sep2016 - Add spontaneous emission model; GBT data point
dsb 26Sep2016 - Create separate plots for different models
dsb 16Nov2016 - Need divide HII region size by 2 to get radius
dsb 07Dec2016 - Add obsDV and smDens to inputs.
"""
import numpy as np
import math
import matplotlib
matplotlib.rc_file('matplotlibrc')
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u
from scipy.optimize import curve_fit
import utils
import argparse

import pdb
#pdb.set_trace()


# Define power law function
def powFunction(x, a, b):
    return a*np.array(x)**b


def plot_model(filenames,labels=None,
               freqrange=[1.,120.],
               linerange=[1.e-25,1.e-21],
               contrange=[1.e-29,1.e-26],
               obsDV = 75.0,
               smDen=[10.0, 25.0, 100.0]):
    """
    Plot model fits on top of the data

    Inputs:
      filenames = list of HII region model pickle files to plot
      labels    = list of labels for legend
      linerange = range for y-axis of line plot (W/m2)
      contrange = range for y-axis of cont plot (W/m2/Hz)
      obsDV     = observed RRL FWHM line width (km/s)
      smDen     = densities to use for simple model (cm-3)
    """

    # --------------------------------
    # Model: Collection of HII regions
    # --------------------------------
    #
    if len(filenames) > 4:
        raise ValueError("Can plot up to 4 models at once")
    if labels is None:
        model_labels = [None] * len(filenames)
    linecolors = ['k','r','g','b']
    #
    # Calculate best fit line and continuum flux models for a range
    # of RRL frequencies
    #
    models = [np.load(filename) for filename in filenames]
    # range of n to plot
    fit_lines = np.arange(110)+40.
    for model in models:
        model.fit_lines = fit_lines
        model.calc_HII_region_properties()
    #
    # Plot RRL
    #
    fig = plt.figure()
    # Models
    for model,linecolor,label in zip(models,linecolors,model_labels):
        plt.plot(model.rrl_freq,model.line_flux*model.num_HII_regions,
                 linecolor,linestyle='-',label=label)
    # VLA data
    plt.errorbar(model.linedata['freq_GHz'],model.linedata['intflux_Wm2'],
                 yerr=model.linedata['intflux_Wm2_err'],linestyle='',
                 color='k',marker='o')
    # GBT data
    # (H52+H53) integrated intensity = 147.1 (5.7) mJy km/s
    # H52alpha freq = 45.453 GHz
    freqGBT = 45.45373
    sInt_mks = 147.1*1.e-3*1.e-26*1.e3*freqGBT*1.e9/c.c.value
    sIntErr_mks = 5.7*1.e-3*1.e-26*1.e3*freqGBT*1.e9/c.c.value
    #plt.errorbar([freqGBT], [sInt_mks], yerr=[sIntErr_mks], color='g', marker='o')
    # labels, etc.
    plt.xlabel("Frequency (GHz)", fontsize=20)
    plt.ylabel(r"Integrated Line Flux (W$\,$m$^{-2}$)", fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(freqrange)
    plt.ylim(linerange)
    plt.text(freqrange[0]*(1.2), linerange[1]*(1-0.5), 'Multiple Compact HII Region Model', fontsize=20)
    fooNe = "%.2g" % (model.electron_dens)
    fooTe = "%4.0f" % (model.electron_temp)
    plt.text(freqrange[0]*(1.2), linerange[1]*(1-0.75), r'(n$_{\rm e}$, T$_{\rm e}$)  = (' + fooNe + ' cm$^{-3}$, '  + fooTe + r' K)', fontsize=20)

    if label is not None:
        plt.legend(fontsize=12,loc='best')
    plt.tight_layout()
    plt.savefig('linefit.eps')
    #
    # Plot radio continuum
    #
    fig = plt.figure()
    for model,linecolor,label in zip(models,linecolors,model_labels):
        if label is None:
            mylabel = ''
        else:
            mylabel = label+', '
        plt.plot(model.rrl_freq,model.thermal_fluxden_rrl+model.nonthermal_fluxden_rrl,
                 linecolor,linestyle='-',label=mylabel+'Total')
        plt.plot(model.rrl_freq,model.thermal_fluxden_rrl,
                 linecolor,linestyle=':',label=mylabel+'Thermal')
        plt.plot(model.rrl_freq,model.nonthermal_fluxden_rrl,
                 linecolor,linestyle='--',label=mylabel+'Non-thermal')
    plt.errorbar(model.contdata['freq_GHz'],model.contdata['S_C_mJy']*1.e-29,
                 yerr=model.contdata['S_C_mJy_err']*1.e-29,linestyle='',
                 color='k',marker='o')
    plt.xlabel("Frequency (GHz)", fontsize=20)
    plt.ylabel(r"Continuum Flux Density (W$\,$m$^{-2}\,$Hz$^{-1}$)", fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(freqrange)
    plt.ylim(contrange)
    plt.text(freqrange[0]*(1.2), contrange[1]*(1-0.5), 'Multiple Compact HII Region Model', fontsize=20)
    plt.legend(fontsize=16,loc='best')
    plt.tight_layout()
    plt.savefig('contfit.eps')    

    # --------------------------------
    # Model: Spontaneous emission
    # --------------------------------
    #
    # Calculate simple model (sm): spontaneous, optically thin emission in LTE
    #  Use equations 4-9 from Puxley et al. (1991, MNRAS, 248, 585)
    #  
    smTe = model.electron_temp              # K
    smDist = model.distance                 # Mpc
    smFreq = model.linedata['freq_GHz'][0]  # GHz use Ka-band freq
    smOmega = model.omega_region            # steradians
    smAngSize = math.sqrt(smOmega/math.pi)  # radians
    smLinSize = smDist*1.e6*smAngSize       # pc
    # line width (Hz)  
    smDV = (obsDV*1.e3)*(smFreq*1.e9)/c.c.value
    # Planck Function (mks)
    smB = (2.0*c.h.value*(smFreq*1.e9)**3/c.c.value**2)*(1.0/(math.exp((c.h.value*smFreq*1.e9)/(c.k_B.value*smTe)) - 1.0)) 
    # set density list to an np array 
    smDen = np.array(smDen)
    # EM cm-6 pc
    smEM = smDen**2*smLinSize
    # Continuum optical depth (Eq. 8)
    smTauC = 8.2e-2*smTe**(-1.35)*smFreq**(-2.1)*smEM
    # LTE line optical depth (Eq. 7)
    smTauL = 1.7e3*(smDV*1.e-3)**(-1)*smTe**(-2.5)*smEM
    # thermal continuum flux density (Eq. 5) 
    smContFlux = smOmega*smB*(1.0 - np.exp(-smTauC))
    # spontaneous line emission (Eq. 4)
    smLineEmission = smOmega*smB*np.exp(-smTauC)*(1.0 - np.exp(-smTauL))
    # integrated line emission (mks)
    smIntLine = 1.064*smLineEmission*smDV
    # Line flux (optically thin) = a * freq^2
    aLine = smIntLine/smFreq**2
    # Cont flux (optically thin) = a * freq^(-0.1)
    aCont = smContFlux/smFreq**(-0.1)
    # number of H-ionizing photons (assume a sphere)
    #  Assume case B recombination rate (cm^3 s-1)
    #  Use Hui & Gnedin (1997) approximation
    lambda_hi = 315614.0/smTe
    alpha_B = (2.753e-14*(lambda_hi**1.5/(1.0 + (lambda_hi/2.740)**0.407)**2.242))
    smRadius = (smLinSize*u.pc).to(u.cm).value/2.0
    smNL = (4.0/3.0)*math.pi*smRadius**3*smDen**2*alpha_B
    # SFR from Murphy et al. (2011)
    smSFR = 7.29e-54*smNL

    # Plot RRL
    #
    fig = plt.figure()
    # Simple model
    xx = np.array([model.rrl_freq[len(model.rrl_freq)-1], model.rrl_freq[0]])
    for i in range(len(aLine)):
        yy = aLine[i]*xx**2
        if smTauC[i] < 1.0 and smTauL[i] < 1.0:
            fooNe = "%3.0f" % (smDen[i])
            fooTe = "%4.0f" % (smTe)
            plt.plot(xx, yy, '--', label=r'(n$_{\rm e}$, T$_{\rm e}$)  = (' + fooNe + ' cm$^{-3}$, '  + fooTe + r' K)')
    # VLA data
    plt.errorbar(model.linedata['freq_GHz'],model.linedata['intflux_Wm2'],
                 yerr=model.linedata['intflux_Wm2_err'],linestyle='',
                 color='k',marker='o')
    # GBT data
    # (H52+H53) integrated intensity = 147.1 (5.7) mJy km/s
    # H52alpha freq = 45.453 GHz
    freqGBT = 45.45373
    sInt_mks = 147.1*1.e-3*1.e-26*1.e3*freqGBT*1.e9/c.c.value
    sIntErr_mks = 5.7*1.e-3*1.e-26*1.e3*freqGBT*1.e9/c.c.value
    #plt.errorbar([freqGBT], [sInt_mks], yerr=[sIntErr_mks], color='g', marker='o')
    # labels, etc.
    plt.xlabel("Frequency (GHz)", fontsize=20)
    plt.ylabel(r"Integrated Line Flux (W$\,$m$^{-2}$)", fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(freqrange)
    plt.ylim(linerange)
    plt.text(freqrange[0]*(1.2), linerange[1]*(1-0.5), 'Spontaneous Emission Model', fontsize=20)
    plt.legend(fontsize=12, loc=4)
    plt.tight_layout()
    plt.savefig('linefitSimple.eps')
    #
    # Plot radio continuum (use second density model only)
    #
    fig = plt.figure()
    # non-thermal = total observed - thermal
    xObs = model.contdata['freq_GHz']
    yObs = model.contdata['S_C_mJy']*1.e-29
    yThermal = aCont[1]*xObs**(-0.1)
    yNonthermal = yObs - yThermal
    # Fit Power law
    p0 = [yObs.mean(), -0.7]
    popt, pcov = curve_fit(powFunction, xObs, yNonthermal, p0)
    # Calculate various contributions over extended frequency range
    xx = np.array([model.rrl_freq[len(model.rrl_freq)-1], model.rrl_freq[0]])
    # Thermal
    yyThermal = aCont[1]*xx**(-0.1)
    # Non-thermal
    yyNonthermal = popt[0]*xx**(popt[1])
    # Total
    yyTotal = yyThermal + yyNonthermal
    # plot model
    plt.plot(xx, yyTotal, '-k', label='Total')
    plt.plot(xx, yyThermal, ':k', label='Thermal')
    plt.plot(xx, yyNonthermal, '--k', label='Non-thermal')
    # plot data
    plt.errorbar(model.contdata['freq_GHz'],model.contdata['S_C_mJy']*1.e-29,
                 yerr=model.contdata['S_C_mJy_err']*1.e-29,linestyle='',
                 color='k',marker='o')
    plt.xlabel("Frequency (GHz)", fontsize=20)
    plt.ylabel(r"Continuum Flux Density (W$\,$m$^{-2}\,$Hz$^{-1}$)", fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(freqrange)
    plt.ylim(contrange)
    plt.text(freqrange[0]*(1.2), contrange[1]*(1-0.5), 'Spontaneous Emission Model', fontsize=20)
    plt.legend(fontsize=16,loc='best')
    plt.tight_layout()
    plt.savefig('contfitSimple.eps')    


    # print results (second density model only)
    print("=============================")
    print("Spontaneous Emissionn Model")
    print("=============================")
    print("Electron temperature = {:g} K".format(smTe))
    print("Electron density = {:g} cm-3".format(smDen[1]))
    print("HII region size = {:g} pc".format(smLinSize))
    print("Non-thermal spectral index = {:g}".format(popt[1]))
    #print("Mass of ionized gas = {:g} Msun".format(best_model.mass_ion))
    print("Number of H-ionizing photons = {:g} s-1".format(smNL[1]))
    print("Star Formation Rate = {:g} Msun/yr".format(smSFR[1]))
    print("=============================")

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot HII region model(s) and data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('filenames',nargs="+",help='filenames of model pickle files')
    parser.add_argument('--labels',nargs="+",default=None,help='labels for each model')
    parser.add_argument('--freqrange',nargs=2,default=[1.,120.],type=float,help='x-axis range (min max) (GHz)')
    parser.add_argument('--linerange',nargs=2,default=[1.e-25,1.e-21],type=float,help='y-axis range (min max) for line plots (W/m2)')
    parser.add_argument('--contrange',nargs=2,default=[1.e-29,1.e-26],type=float,help='y-axis range (min max) for cont plots (W/m2/Hz)')
    parser.add_argument('--obsDV',nargs=1,default=75.0,type=float,help='observed RRL FWHM line width (km/s)')
    parser.add_argument('--smDen',nargs=3,default=[10., 25.0, 100.0],type=float,help='densities to use for simple model (cm-3)')
    args = parser.parse_args()
    plot_model(args.filenames,labels=args.labels,
               linerange=args.linerange,contrange=args.contrange,
               obsDV=args.obsDV,smDen=args.smDen)
