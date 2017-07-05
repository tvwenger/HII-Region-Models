#
# IC342 (only two RRLs: stacked C-band and stacked Ka-band)
#

#--------------------
# wholeRegion
#--------------------
#
import fit_model
import numpy as np
linefile = 'wholeRegion_lineKaC.txt'
contfile = 'wholeRegion_contKaC.txt'
omega_region = 864.0*0.5*0.5 # sq arcsec (864 pixels; 0.5x0.5 arcsec)
distance = 3.280 # Mpc
v_turbulent = 12.0 # km/s
heII_abundance = 0.07 # He+/H+
electron_temps = np.array([9000.0])
electron_dens = np.logspace(2.,6.,10.) # 100 to 10^6 cm-3
HII_region_sizes = np.logspace(-2.,2.,10.) # 0.01 to 10^2 pc

# use H57 line flux to calculate number of HII regions
best_model = fit_model.fit_model(linefile,contfile,omega_region,distance,n_for_num=57,electron_temps=electron_temps,electron_dens=electron_dens,HII_region_sizes=HII_region_sizes,v_turbulent=v_turbulent,heII_abundance=heII_abundance,verbose=True,plotRSS=True,plotModel=False)

from plot_model import *
plt = plot_model(['best_model.pkl'], labels=None, freqrange=[1.,120.], linerange=[1.e-24,1.e-20], contrange=[1.e-29,1.e-26], obsDV = 75.0, smDen=[10.0, 25.0, 100.0])


Ran 100 models in 00h 00m 9.34s
00h 00m 0.09s per model
14 models passed, 86 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 1.55979e-43
Electron temperature = 7500 K
Electron density = 46415.9 cm-3
HII region size = 0.0774264 pc
Number of HII_regions = 3494.22
Non-thermal spectral index = -0.658973
Mass of ionized gas = 974.334 Msun
Number of H-ionizing photons = 1.76561e+52 s-1
Star Formation Rate = 0.128713 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 7500 K
Electron density = 25 cm-3
HII region size = 131.856 pc
Non-thermal spectral index = -2.1369
Number of H-ionizing photons = 7.23981e+51 s-1
Star Formation Rate = 0.0527782 Msun/yr
=============================



Ran 100 models in 00h 00m 8.90s
00h 00m 0.09s per model
18 models passed, 82 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 1.55978e-43
Electron temperature = 6000 K
Electron density = 46415.9 cm-3
HII region size = 0.0278256 pc
Number of HII_regions = 66881.2
Non-thermal spectral index = -0.634816
Mass of ionized gas = 865.623 Msun
Number of H-ionizing photons = 1.8776e+52 s-1
Star Formation Rate = 0.136877 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 6000 K
Electron density = 25 cm-3
HII region size = 131.856 pc
Non-thermal spectral index = -2.5226
Number of H-ionizing photons = 8.66593e+51 s-1
Star Formation Rate = 0.0631747 Msun/yr
=============================


Ran 100 models in 00h 00m 7.91s
00h 00m 0.08s per model
9 models passed, 91 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 1.55979e-43
Electron temperature = 9000 K
Electron density = 46415.9 cm-3
HII region size = 0.0774264 pc
Number of HII_regions = 4950.43
Non-thermal spectral index = -1.02981
Mass of ionized gas = 1380.39 Msun
Number of H-ionizing photons = 2.15416e+52 s-1
Star Formation Rate = 0.157038 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 9000 K
Electron density = 25 cm-3
HII region size = 131.856 pc
Non-thermal spectral index = -1.85291
Number of H-ionizing photons = 6.23472e+51 s-1
Star Formation Rate = 0.0454511 Msun/yr
=============================




#--------------------
# eastComp
#--------------------
#
import fit_model
import numpy as np
linefile = 'eastComp_lineKaC.txt'
contfile = 'eastComp_contKaC.txt'
omega_region = 90.0*0.5*0.5 # sq arcsec (90 pixels; 0.5x0.5 arcsec)
distance = 3.280 # Mpc
v_turbulent = 12.0 # km/s
heII_abundance = 0.07 # He+/H+
electron_temps = np.array([9000.0])
electron_dens = np.logspace(2.,6.,10.) # 100 to 10^6 cm-3
HII_region_sizes = np.logspace(-2.,2.,10.) # 0.01 to 10^2 pc

# use H57 line flux to calculate number of HII regions
best_model = fit_model.fit_model(linefile,contfile,omega_region,distance,n_for_num=57,electron_temps=electron_temps,electron_dens=electron_dens,HII_region_sizes=HII_region_sizes,v_turbulent=v_turbulent,heII_abundance=heII_abundance,verbose=True,plotRSS=True,plotModel=False)

from plot_model import *
plt = plot_model(['best_model.pkl'], labels=None, freqrange=[1.,120.], linerange=[1.e-25,1.e-21], contrange=[1.e-30,1.e-27], obsDV = 75.0, smDen=[25.0, 50.0, 150.0])



Ran 100 models in 00h 00m 7.35s
00h 00m 0.07s per model
6 models passed, 94 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 5.41078e-45
Electron temperature = 7500 K
Electron density = 46415.9 cm-3
HII region size = 0.0774264 pc
Number of HII_regions = 645.781
Non-thermal spectral index = -1.20996
Mass of ionized gas = 180.071 Msun
Number of H-ionizing photons = 3.2631e+51 s-1
Star Formation Rate = 0.023788 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 7500 K
Electron density = 50 cm-3
HII region size = 42.5564 pc
Non-thermal spectral index = -1.79574
Number of H-ionizing photons = 9.73599e+50 s-1
Star Formation Rate = 0.00709754 Msun/yr
=============================


Ran 100 models in 00h 00m 7.20s
00h 00m 0.07s per model
8 models passed, 92 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 5.41066e-45
Electron temperature = 6000 K
Electron density = 46415.9 cm-3
HII region size = 0.0774264 pc
Number of HII_regions = 408.346
Non-thermal spectral index = -0.745791
Mass of ionized gas = 113.864 Msun
Number of H-ionizing photons = 2.4698e+51 s-1
Star Formation Rate = 0.0180048 Msun/yr
======================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 6000 K
Electron density = 50 cm-3
HII region size = 42.5564 pc
Non-thermal spectral index = -2.02949
Number of H-ionizing photons = 1.16538e+51 s-1
Star Formation Rate = 0.00849564 Msun/yr
=============================


Ran 100 models in 00h 00m 7.64s
00h 00m 0.08s per model
2 models passed, 98 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 5.4126e-45
Electron temperature = 9000 K
Electron density = 129155 cm-3
HII region size = 0.0774264 pc
Number of HII_regions = 56.6202
Non-thermal spectral index = -0.749711
Mass of ionized gas = 43.9313 Msun
Number of H-ionizing photons = 1.90763e+51 s-1
Star Formation Rate = 0.0139067 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 9000 K
Electron density = 50 cm-3
HII region size = 42.5564 pc
Non-thermal spectral index = -1.6124
Number of H-ionizing photons = 8.38436e+50 s-1
Star Formation Rate = 0.0061122 Msun/yr
=============================




#--------------------
# westComp
#--------------------
#
import fit_model
import numpy as np
linefile = 'westComp_lineKaC.txt'
contfile = 'westComp_contKaC.txt'
omega_region = 90.0*0.5*0.5 # sq arcsec (90 pixels; 0.5x0.5 arcsec)
distance = 3.280 # Mpc
v_turbulent = 12.0 # km/s
heII_abundance = 0.07 # He+/H+
electron_temps = np.array([9000.0])
electron_dens = np.logspace(2.,6.,10.) # 100 to 10^6 cm-3
HII_region_sizes = np.logspace(-2.,2.,10.) # 0.01 to 10^2 pc

# use H57 line flux to calculate number of HII regions
best_model = fit_model.fit_model(linefile,contfile,omega_region,distance,n_for_num=57,electron_temps=electron_temps,electron_dens=electron_dens,HII_region_sizes=HII_region_sizes,v_turbulent=v_turbulent,heII_abundance=heII_abundance,verbose=True,plotRSS=True,plotModel=False)

from plot_model import *
plt = plot_model(['best_model.pkl'], labels=None, freqrange=[1.,120.], linerange=[1.e-25,1.e-21], contrange=[1.e-30,1.e-27], obsDV = 75.0, smDen=[25.0, 70.0, 150.0])


Ran 100 models in 00h 00m 7.34s
00h 00m 0.07s per model
11 models passed, 89 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 9.82054e-45
Electron temperature = 7500 K
Electron density = 46415.9 cm-3
HII region size = 0.0278256 pc
Number of HII_regions = 24636.2
Non-thermal spectral index = -0.952361
Mass of ionized gas = 318.858 Msun
Number of H-ionizing photons = 5.7781e+51 s-1
Star Formation Rate = 0.0421223 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 7500 K
Electron density = 70 cm-3
HII region size = 42.5564 pc
Non-thermal spectral index = -2.86678
Number of H-ionizing photons = 1.90825e+51 s-1
Star Formation Rate = 0.0139112 Msun/yr
=============================



Ran 100 models in 00h 00m 10.66s
00h 00m 0.11s per model
17 models passed, 83 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 9.82041e-45
Electron temperature = 6000 K
Electron density = 46415.9 cm-3
HII region size = 0.0278256 pc
Number of HII_regions = 16823.2
Non-thermal spectral index = -0.571386
Mass of ionized gas = 217.737 Msun
Number of H-ionizing photons = 4.72289e+51 s-1
Star Formation Rate = 0.0344299 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 6000 K
Electron density = 70 cm-3
HII region size = 42.5564 pc
Non-thermal spectral index = -4.0287
Number of H-ionizing photons = 2.28415e+51 s-1
Star Formation Rate = 0.0166515 Msun/yr
=============================


Ran 100 models in 00h 00m 7.59s
00h 00m 0.08s per model
8 models passed, 92 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 9.82042e-45
Electron temperature = 9000 K
Electron density = 46415.9 cm-3
HII region size = 0.0774264 pc
Number of HII_regions = 1245.22
Non-thermal spectral index = -1.02525
Mass of ionized gas = 347.22 Msun
Number of H-ionizing photons = 5.41853e+51 s-1
Star Formation Rate = 0.0395011 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 9000 K
Electron density = 70 cm-3
HII region size = 42.5564 pc
Non-thermal spectral index = -2.40772
Number of H-ionizing photons = 1.64333e+51 s-1
Star Formation Rate = 0.0119799 Msun/yr
=============================




















#--------------------
# innerNucleus Region
#--------------------
#
import fit_model
import numpy as np
linefile = 'innerNucleus_lineKaC_dsb.txt'
contfile = 'innerNucleus_contKaC_dsb.txt'
omega_region = 737.0*0.5*0.5 # sq arcsec (737 pixels; 0.5x0.5 arcsec)
distance = 3.280 # Mpc
v_turbulent = 12.0 # km/s
heII_abundance = 0.07 # He+/H+
electron_temps = np.array([7500.0])
electron_dens = np.logspace(1.,5.,10.) # 10 to 10^5 cm-3
HII_region_sizes = np.logspace(-1.,3.,10.) # 0.1 to 10^3 pc

# use H57 line flux to calculate number of HII regions
best_model = fit_model.fit_model(linefile,contfile,omega_region,distance,n_for_num=57,electron_temps=electron_temps,electron_dens=electron_dens,HII_region_sizes=HII_region_sizes,v_turbulent=v_turbulent,heII_abundance=heII_abundance,verbose=True,plotRSS=True,plotModel=True)

from plot_model import *
plt = plot_model(['best_model.pkl'], labels=None, freqrange=[1.,120.], linerange=[1.e-25,1.e-21], contrange=[1.e-29,1.e-26], obsDV = 75.0, smDen=[1.0, 10.0, 100.0])



================================
Multile Compact HII Region Model
================================
Best model RSS: 3.3892e-45
Electron temperature = 7500 K
Electron density = 35938.1 cm-3
HII region size = 0.1 pc
Number of HII_regions = 414.902
Non-thermal spectral index = -0.41649
Mass of ionized gas = 192.986 Msun
Number of H-ionizing photons = 3.37749e+50 s-1
Star Formation Rate = 0.00246219 Msun/yr
================================

=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 7500 K
Electron density = 10 cm-3
HII region size = 121.78 pc
Non-thermal spectral index = -0.429406
Number of H-ionizing photons = 1.13834e+50 s-1
Star Formation Rate = 0.000829846 Msun/yr
=============================



================================
Multile Compact HII Region Model
================================
Best model RSS: 3.38872e-45
Electron temperature = 6000 K
Electron density = 35938.1 cm-3
HII region size = 0.1 pc
Number of HII_regions = 265.331
Non-thermal spectral index = -0.408279
Mass of ionized gas = 123.415 Msun
Number of H-ionizing photons = 2.69612e+50 s-1
Star Formation Rate = 0.00196547 Msun/yr
================================


================================
Multile Compact HII Region Model
================================
Best model RSS: 3.38959e-45
Electron temperature = 9000 K
Electron density = 4641.59 cm-3
HII region size = 0.1 pc
Number of HII_regions = 52598.9
Non-thermal spectral index = -0.448047
Mass of ionized gas = 3159.86 Msun
Number of H-ionizing photons = 5.93861e+50 s-1
Star Formation Rate = 0.00432925 Msun/yr
================================




# use H104 line flux to calculate number of HII regions
best_model = fit_model.fit_model(linefile,contfile,omega_region,distance,n_for_num=104,electron_temps=electron_temps,electron_dens=electron_dens,HII_region_sizes=HII_region_sizes,v_turbulent=v_turbulent,heII_abundance=heII_abundance,verbose=True,plotRSS=True,plotModel=True)


================================
Multile Compact HII Region Model
================================
Best model RSS: 3.77219e-46
Electron temperature = 7500 K
Electron density = 12915.5 cm-3
HII region size = 0.278256 pc
Number of HII_regions = 62.8382
Non-thermal spectral index = -0.403087
Mass of ionized gas = 226.304 Msun
Number of H-ionizing photons = 1.42337e+50 s-1
Star Formation Rate = 0.00103763 Msun/yr
================================

=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 7500 K
Electron density = 10 cm-3
HII region size = 121.78 pc
Non-thermal spectral index = -0.429406
Number of H-ionizing photons = 1.13834e+50 s-1
Star Formation Rate = 0.000829846 Msun/yr
=============================



================================
Multile Compact HII Region Model
================================
Best model RSS: 5.28561e-46
Electron temperature = 6000 K
Electron density = 1668.1 cm-3
HII region size = 0.774264 pc
Number of HII_regions = 215.901
Non-thermal spectral index = -0.406586
Mass of ionized gas = 2163.56 Msun
Number of H-ionizing photons = 2.19385e+50 s-1
Star Formation Rate = 0.00159932 Msun/yr
================================


================================
Multile Compact HII Region Model
================================
Best model RSS: 3.72684e-46
Electron temperature = 9000 K
Electron density = 12915.5 cm-3
HII region size = 0.278256 pc
Number of HII_regions = 85.4097
Non-thermal spectral index = -0.406136
Mass of ionized gas = 307.593 Msun
Number of H-ionizing photons = 1.60856e+50 s-1
Star Formation Rate = 0.00117264 Msun/yr
================================










=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  7089084853.56
n =  [ 57 104]
PeakLineFluxHII =  [  1.81236730e-39   3.37384951e-40]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.49526198e-33   1.70402813e-34]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  329042772.138
n =  [ 57 104]
PeakLineFluxHII =  [  3.90462730e-38   7.26881351e-39]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.18391840e-31   3.67125524e-33]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  15272373.6869
n =  [ 57 104]
PeakLineFluxHII =  [  8.41226501e-37   1.56606340e-37]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.55067503e-30   7.90970693e-32]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  708824.113679
n =  [ 57 104]
PeakLineFluxHII =  [  1.81236801e-35   3.37425111e-36]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.49526414e-29   1.70423096e-30]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  32893.3828481
n =  [ 57 104]
PeakLineFluxHII =  [  3.90463119e-34   7.27122096e-35]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.18391958e-27   3.67247116e-29]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1525.83120125
n =  [ 57 104]
PeakLineFluxHII =  [  8.41228889e-33   1.56750665e-33]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.55068227e-26   7.91699632e-28]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  70.7011183932
n =  [ 57 104]
PeakLineFluxHII =  [  1.81238230e-31   3.38290342e-32]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.49530746e-25   1.70860098e-26]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  3.26603784804
n =  [ 57 104]
PeakLineFluxHII =  [  3.90471684e-30   7.32309502e-31]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.18394555e-23   3.69867116e-25]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.14961416549
n =  [ 57 104]
PeakLineFluxHII =  [  8.41280235e-29   1.59861236e-29]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.55083795e-22   8.07410178e-24]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+01, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.00670049811636
n =  [ 57 104]
PeakLineFluxHII =  [  1.81269012e-27   3.56951156e-28]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.49624079e-21   1.80285104e-22]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  867047686.709
n =  [ 57 104]
PeakLineFluxHII =  [  1.41008860e-38   2.75849943e-39]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  4.27551649e-32   1.39323363e-33]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  40240553.8541
n =  [ 57 104]
PeakLineFluxHII =  [  3.03794464e-37   5.94363228e-38]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  9.21132361e-31   3.00194675e-32]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1867254.32875
n =  [ 57 104]
PeakLineFluxHII =  [  6.54505819e-36   1.28089169e-36]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.98452099e-29   6.46939188e-31]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  86599.7323388
n =  [ 57 104]
PeakLineFluxHII =  [  1.41009296e-34   2.76184520e-35]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  4.27552971e-28   1.39492348e-29]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  4010.52074834
n =  [ 57 104]
PeakLineFluxHII =  [  3.03797067e-33   5.96369075e-34]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  9.21140254e-27   3.01207767e-28]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  184.988531468
n =  [ 57 104]
PeakLineFluxHII =  [  6.54521424e-32   1.29291829e-32]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.98456830e-25   6.53013455e-27]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  8.43956475989
n =  [ 57 104]
PeakLineFluxHII =  [  1.41018651e-30   2.83397381e-31]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  4.27581336e-24   1.43135343e-25]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.373909124342
n =  [ 57 104]
PeakLineFluxHII =  [  3.03853150e-29   6.39660921e-30]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  9.21310301e-23   3.23073154e-24]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.0153977168496
n =  [ 57 104]
PeakLineFluxHII =  [  6.54857640e-28   1.55331506e-28]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.98558774e-21   7.84531896e-23]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 2.7826e+01, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.000542387137895
n =  [ 57 104]
PeakLineFluxHII =  [  1.41220227e-26   4.40967416e-27]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  4.28192534e-20   2.22719146e-21]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  105159339.001
n =  [ 57 104]
PeakLineFluxHII =  [  1.10036386e-37   2.27440622e-38]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.33640301e-31   1.14873297e-32]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  4877383.78549
n =  [ 57 104]
PeakLineFluxHII =  [  2.37066788e-36   4.90375713e-37]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  7.18807996e-30   2.47673764e-31]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  225914.084976
n =  [ 57 104]
PeakLineFluxHII =  [  5.10748385e-35   1.05869917e-35]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.54863542e-28   5.34716549e-30]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  10425.2516073
n =  [ 57 104]
PeakLineFluxHII =  [  1.10039487e-33   2.29418976e-34]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.33649701e-27   1.15872503e-28]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  476.216068457
n =  [ 57 104]
PeakLineFluxHII =  [  2.37085372e-32   5.02240623e-33]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  7.18864346e-26   2.53666367e-27]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  21.1676048113
n =  [ 57 104]
PeakLineFluxHII =  [  5.10859800e-31   1.12991081e-31]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.54897323e-24   5.70683368e-26]
Checking this model:
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.87851627968
n =  [ 57 104]
PeakLineFluxHII =  [  1.10106281e-29   2.72248859e-30]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.33852227e-23   1.37504566e-24]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.0314144114335
n =  [ 57 104]
PeakLineFluxHII =  [  2.37485847e-28   7.61354563e-29]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  7.20078620e-22   3.84536888e-23]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.000878090569615
n =  [ 57 104]
PeakLineFluxHII =  [  5.13261464e-27   2.72380849e-27]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.55625530e-20   1.37571230e-21]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 7.7426e+01, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1.8374137433e-05
n =  [ 57 104]
PeakLineFluxHII =  [  1.11547518e-25   1.30169406e-25]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.38222190e-19   6.57445828e-20]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  12865829.655
n =  [ 57 104]
PeakLineFluxHII =  [  8.65091367e-37   1.85899442e-37]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.62303546e-30   9.38921183e-32]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  594666.450446
n =  [ 57 104]
PeakLineFluxHII =  [  1.86383419e-35   4.02200351e-36]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.65131425e-29   2.03139087e-30]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  27282.4655236
n =  [ 57 104]
PeakLineFluxHII =  [  4.01581678e-34   8.76662172e-35]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.21763206e-27   4.42775232e-29]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1226.78623454
n =  [ 57 104]
PeakLineFluxHII =  [  8.65365999e-33   1.94960661e-33]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.62386818e-26   9.84686628e-28]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  52.3774379355
n =  [ 57 104]
PeakLineFluxHII =  [  1.86548070e-31   4.56637561e-32]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.65630663e-25   2.30633656e-26]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1.9849101429
n =  [ 57 104]
PeakLineFluxHII =  [  4.02568964e-30   1.20496666e-30]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.22062560e-23   6.08591782e-25]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.0605315060999
n =  [ 57 104]
PeakLineFluxHII =  [  8.71288393e-29   3.95124903e-29]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.64182541e-22   1.99565495e-23]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.00139354642868
n =  [ 57 104]
PeakLineFluxHII =  [  1.90104756e-27   1.71630489e-27]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.76414857e-21   8.66853070e-22]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  2.40173879597e-05
n =  [ 57 104]
PeakLineFluxHII =  [  4.23996494e-26   9.95841243e-26]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.28559581e-19   5.02968933e-20]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 2.1544e+02, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  2.61426639713e-07
n =  [ 57 104]
PeakLineFluxHII =  [  1.00153557e-24   9.14884019e-24]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.03674665e-18   4.62079917e-18]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1589063.51581
n =  [ 57 104]
PeakLineFluxHII =  [  6.93423327e-36   1.50513213e-36]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.10252240e-29   7.60196171e-31]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  72302.1782092
n =  [ 57 104]
PeakLineFluxHII =  [  1.49448216e-34   3.30799238e-35]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  4.53140543e-28   1.67076571e-29]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  3181.08732784
n =  [ 57 104]
PeakLineFluxHII =  [  3.22304345e-33   7.51865731e-34]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  9.77255997e-27   3.79744369e-28]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  128.903131351
n =  [ 57 104]
PeakLineFluxHII =  [  6.96350675e-32   1.85546350e-32]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.11139839e-25   9.37137824e-27]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  4.41079103149
n =  [ 57 104]
PeakLineFluxHII =  [  1.51205349e-30   5.42249799e-31]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  4.58468330e-24   2.73873777e-25]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.11698475654
n =  [ 57 104]
PeakLineFluxHII =  [  3.32875486e-29   2.04449761e-29]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.00930866e-22   1.03261317e-23]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.00237321706986
n =  [ 57 104]
PeakLineFluxHII =  [  7.60353501e-28   1.00780943e-27]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.30546076e-21   5.09013699e-22]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  3.63186389707e-05
n =  [ 57 104]
PeakLineFluxHII =  [  1.90656941e-26   6.58546305e-26]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.78089132e-20   3.32611583e-20]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  2.72953887978e-07
n =  [ 57 104]
PeakLineFluxHII =  [  5.88999251e-25   8.76247108e-24]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.78589914e-18   4.42565596e-18]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 5.9948e+02, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  9.13316637941e-11
n =  [ 57 104]
PeakLineFluxHII =  [  2.71160193e-23   2.61875285e-20]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  8.22182293e-17   1.32265191e-14]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  193561.221545
n =  [ 57 104]
PeakLineFluxHII =  [  5.72725486e-35   1.23565585e-35]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.73655561e-28   6.24091949e-30]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  8288.39284628
n =  [ 57 104]
PeakLineFluxHII =  [  1.23874498e-33   2.88566263e-34]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.75598711e-27   1.45745987e-28]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  316.565547796
n =  [ 57 104]
PeakLineFluxHII =  [  2.69791526e-32   7.55530905e-33]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  8.18032370e-26   3.81595536e-27]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  9.8525348892
n =  [ 57 104]
PeakLineFluxHII =  [  5.98827982e-31   2.42754842e-31]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.81570074e-24   1.22608041e-25]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.239319696537
n =  [ 57 104]
PeakLineFluxHII =  [  1.39760552e-29   9.99395613e-30]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  4.23766667e-23   5.04764137e-24]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.00477620290495
n =  [ 57 104]
PeakLineFluxHII =  [  3.69172154e-28   5.00764016e-28]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.11936344e-21   2.52920578e-22]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  7.92993991482e-05
n =  [ 57 104]
PeakLineFluxHII =  [  1.27284415e-26   3.01610173e-26]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.85937889e-20   1.52334067e-20]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  6.33848008637e-07
n =  [ 57 104]
PeakLineFluxHII =  [  7.32021987e-25   3.77338181e-24]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.21955704e-18   1.90581967e-18]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1.58089980762e-10
n =  [ 57 104]
PeakLineFluxHII =  [  1.44632439e-22   1.51290457e-20]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  4.38538670e-16   7.64121796e-15]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.6681e+03, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  6.67642862509e-19
n =  [ 57 104]
PeakLineFluxHII =  [  2.15863370e-18   3.58238017e-12]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  6.54517310e-12   1.80935058e-06]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  21692.0348372
n =  [ 57 104]
PeakLineFluxHII =  [  4.85160388e-34   1.10259391e-34]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.47105029e-27   5.56886437e-29]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  784.457318919
n =  [ 57 104]
PeakLineFluxHII =  [  1.07594629e-32   3.04892375e-33]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.26236670e-26   1.53991807e-27]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  23.0617849548
n =  [ 57 104]
PeakLineFluxHII =  [  2.50583800e-31   1.03710556e-31]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  7.59792804e-25   5.23810279e-26]
Checking this model:
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.573615324589
n =  [ 57 104]
PeakLineFluxHII =  [  6.58988861e-30   4.16960713e-30]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.99811398e-23   2.10594095e-24]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.0144886895678
n =  [ 57 104]
PeakLineFluxHII =  [  2.25923943e-28   1.65077079e-28]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  6.85021881e-22   8.33753801e-23]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.000468252259673
n =  [ 57 104]
PeakLineFluxHII =  [  1.30321844e-26   5.10782490e-27]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.95147648e-20   2.57980602e-21]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  2.10552694407e-05
n =  [ 57 104]
PeakLineFluxHII =  [  2.69901274e-24   1.13593918e-25]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  8.18365137e-18   5.73728113e-20]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1.29104005226e-06
n =  [ 57 104]
PeakLineFluxHII =  [  4.87813871e-20   1.85257657e-24]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.47909589e-13   9.35679724e-19]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1.18677353473e-07
n =  [ 57 104]
PeakLineFluxHII =  [  7.74076520e-10   2.01533863e-23]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.34707020e-03   1.01788586e-17]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 4.6416e+03, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  1.47482299424e-08
n =  [ 57 104]
PeakLineFluxHII =  [  4.87321960e+17   1.62172041e-22]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.47760437e+24   8.19081344e-17]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  2039.76803394
n =  [ 57 104]
PeakLineFluxHII =  [  4.31985004e-33   1.17256007e-33]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.30981770e-26   5.92224204e-28]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  62.838194976
n =  [ 57 104]
PeakLineFluxHII =  [  1.08276901e-31   3.80620505e-32]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.28305383e-25   1.92239768e-26]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  2.06378346038
n =  [ 57 104]
PeakLineFluxHII =  [  3.36067451e-30   1.15891546e-30]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.01898698e-23   5.85332727e-25]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.109413886519
n =  [ 57 104]
PeakLineFluxHII =  [  1.60202837e-28   2.18596617e-29]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  4.85749527e-22   1.10406461e-23]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.0114513978166
n =  [ 57 104]
PeakLineFluxHII =  [  1.98624597e-26   2.08861013e-28]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  6.02247788e-20   1.05489305e-22]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.00147356883338
n =  [ 57 104]
PeakLineFluxHII =  [  5.97598835e-23   1.62310066e-27]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.81197385e-16   8.19778468e-22]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.000190318723617
n =  [ 57 104]
PeakLineFluxHII =  [  5.33151904e-15   1.25670796e-26]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.61656491e-08   6.34724728e-21]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  2.45806083734e-05
n =  [ 57 104]
PeakLineFluxHII =  [  1.83740823e+06   9.73023333e-26]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.57118837e+12   4.91444305e-20]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  3.17470765105e-06
n =  [ 57 104]
PeakLineFluxHII =  [  6.70066354e+61   7.53376629e-25]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.03170194e+68   3.80507477e-19]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
/users/dbalser/errl/models/utils.py:299: RuntimeWarning: overflow encountered in exp
  peak_line_fluxden_background = peak_line_fluxden_background * np.exp(-1.*tau_c) * (np.exp(-1.*tau_l) - 1.)
/users/dbalser/software/python/anaconda3/lib/python3.5/site-packages/astropy/units/quantity.py:795: RuntimeWarning: invalid value encountered in multiply
  return super(Quantity, self).__mul__(other)
=============================================================
(Te,ne,l) = (7.5000e+03, 1.2915e+04, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  nan
n =  [ 57 104]
PeakLineFluxHII =  [  7.01361029e+214   5.83312163e-024]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.12659023e+221               nan]
Checking this model:
/users/dbalser/errl/models/hii_region_model.py:282: RuntimeWarning: invalid value encountered in greater
  if np.any(peak_line_fluxden.to('mJy').value > self.linedata['peak_mJy'][self.fit_lines_ind]):
HII region peak line flux > observed peak line flux
/users/dbalser/errl/models/hii_region_model.py:296: RuntimeWarning: invalid value encountered in less
  if np.any(self.nonthermal_fluxden_cont < 0.):
line flux is nan
number of HII regions is nan
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  256.173216749
n =  [ 57 104]
PeakLineFluxHII =  [  4.72328833e-32   9.33645827e-33]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.43214385e-25   4.71555932e-27]
Checking this model:
This model passed all tests!
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  17.0883873956
n =  [ 57 104]
PeakLineFluxHII =  [  1.75732321e-30   1.39963502e-31]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  5.32836331e-24   7.06912813e-26]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  2.00880983251
n =  [ 57 104]
PeakLineFluxHII =  [  1.20544937e-28   1.19063065e-30]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.65503179e-22   6.01351099e-25]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.259375406902
n =  [ 57 104]
PeakLineFluxHII =  [  4.56955220e-26   9.22119247e-30]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.38552966e-19   4.65734209e-24]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.0334996219893
n =  [ 57 104]
PeakLineFluxHII =  [  6.98301466e-21   7.13963444e-29]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.11731336e-14   3.60601084e-23]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.00432664255585
n =  [ 57 104]
PeakLineFluxHII =  [  4.64758232e-08   5.52795966e-28]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.40918910e-01   2.79200323e-22]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.000558807374364
n =  [ 57 104]
PeakLineFluxHII =  [  5.80764548e+26   4.28009840e-27]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  1.76093077e+33   2.16174670e-21]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  7.21727477168e-05
n =  [ 57 104]
PeakLineFluxHII =  [  1.13710132e+120   3.31392475e-026]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.44779431e+126   1.67376196e-020]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
/users/dbalser/errl/models/utils.py:294: RuntimeWarning: overflow encountered in exp
  peak_line_fluxden_HII_region *= (1. - np.exp(-1.*(tau_l + tau_c)))
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  nan
n =  [ 57 104]
PeakLineFluxHII =  [             inf   2.56585158e-25]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [ nan  nan]
Checking this model:
line flux is nan
number of HII regions is nan
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 3.5938e+04, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  nan
n =  [ 57 104]
PeakLineFluxHII =  [             inf   1.98664569e-24]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [ nan  nan]
Checking this model:
line flux is nan
number of HII regions is nan
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 1.0000e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  281.182118988
n =  [ 57 104]
PeakLineFluxHII =  [  8.79050643e-31   8.50605493e-33]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.66536126e-24   4.29614801e-27]
Checking this model:
HII region peak line flux > observed peak line flux
At least one failed test, skipping this model.
/users/dbalser/errl/models/utils.py:509: RuntimeWarning: invalid value encountered in log10
  spectral_index = np.log10(fluxden[max_ind]/fluxden[min_ind])/np.log10(freq[max_ind]/freq[min_ind])
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 2.7826e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  36.3134504508
n =  [ 57 104]
PeakLineFluxHII =  [  8.30968726e-29   6.58640399e-32]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.51957252e-22   3.32659107e-26]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 7.7426e-01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  4.69006247652
n =  [ 57 104]
PeakLineFluxHII =  [  1.03851145e-25   5.09961341e-31]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  3.14886089e-19   2.57565865e-25]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 2.1544e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.605744862045
n =  [ 57 104]
PeakLineFluxHII =  [  7.24269062e-19   3.94844546e-30]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  2.19604947e-12   1.99423895e-24]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 5.9948e+00)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.0782349573659
n =  [ 57 104]
PeakLineFluxHII =  [  2.07033440e-01   3.05713792e-29]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  6.27744162e+05   1.54406680e-23]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 1.6681e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.0101044332978
n =  [ 57 104]
PeakLineFluxHII =  [  2.01507541e+46   2.36703087e-28]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  6.10989134e+52   1.19551484e-22]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 4.6416e+01)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  0.0013050377441
n =  [ 57 104]
PeakLineFluxHII =  [  2.93354750e+175   1.83270603e-027]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [  8.89478203e+181   9.25643726e-022]
Checking this model:
HII region peak line flux > observed peak line flux
minimum number of HII regions is greater than number of HII regions
overproduces thermal continuum emission (attenuated)
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 1.2915e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  nan
n =  [ 57 104]
PeakLineFluxHII =  [             inf   1.41899772e-26]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [ nan  nan]
Checking this model:
line flux is nan
number of HII regions is nan
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 3.5938e+02)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  nan
n =  [ 57 104]
PeakLineFluxHII =  [             inf   1.09867840e-25]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [ nan  nan]
Checking this model:
line flux is nan
number of HII regions is nan
At least one failed test, skipping this model.
=============================================================
(Te,ne,l) = (7.5000e+03, 1.0000e+05, 1.0000e+03)
(Rsun,Omega,Vturb,y) = (3.2800e+00, 1.8425e+02, 1.2000e+01, 7.0000e-02)
NHII =  nan
n =  [ 57 104]
PeakLineFluxHII =  [             inf   8.50666787e-25]
FWHMLineFluxHII =  24.6831146009
LineFluxHII =  [ nan  nan]
Checking this model:
line flux is nan
number of HII regions is nan
At least one failed test, skipping this model.
 
Ran 100 models in 00h 00m 8.86s
00h 00m 0.09s per model
35 models passed, 65 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 3.77219e-46
Electron temperature = 7500 K
Electron density = 12915.5 cm-3
HII region size = 0.278256 pc
Number of HII_regions = 62.8382
Non-thermal spectral index = -0.403087
Mass of ionized gas = 226.304 Msun
Number of H-ionizing photons = 1.42337e+50 s-1
Star Formation Rate = 0.00103763 Msun/yr
================================
/users/dbalser/software/python/anaconda3/lib/python3.5/site-packages/matplotlib/__init__.py:872: UserWarning: axes.color_cycle is deprecated and replaced with axes.prop_cycle; please use the latter.
  warnings.warn(self.msg_depr % (key, alt_key))
/users/dbalser/software/python/anaconda3/lib/python3.5/site-packages/matplotlib/__init__.py:872: UserWarning: axes.color_cycle is deprecated and replaced with axes.prop_cycle; please use the latter.
  warnings.warn(self.msg_depr % (key, alt_key))












#---------------------
# westComponent Region
#---------------------
#
import fit_model
import numpy as np
linefile = 'westComponent_lineKaC_dsb.txt'
contfile = 'westComponent_contKaC_dsb.txt'
omega_region = 33.0*0.5*0.5 # sq arcsec (33 pixels; 0.5x0.5 arcsec)
distance = 3.280 # Mpc
v_turbulent = 12.0 # km/s
heII_abundance = 0.07 # He+/H+
electron_temps = np.array([7500.0])
electron_dens = np.logspace(1.,5.,10.) # 10 to 10^5 cm-3
HII_region_sizes = np.logspace(-1.,3.,10.) # 0.1 to 10^3 pc

# use H104 line flux to calculate number of HII regions
best_model = fit_model.fit_model(linefile,contfile,omega_region,distance,n_for_num=104,electron_temps=electron_temps,electron_dens=electron_dens,HII_region_sizes=HII_region_sizes,v_turbulent=v_turbulent,heII_abundance=heII_abundance,verbose=True,plotRSS=True,plotModel=True)


#---------------------
# eastComponent Region
#---------------------
#
import fit_model
import numpy as np
linefile = 'eastComponent_lineKaC_dsb.txt'
contfile = 'eastComponent_contKaC_dsb.txt'
omega_region = 34.0*0.5*0.5 # sq arcsec (34 pixels; 0.5x0.5 arcsec)
distance = 3.280 # Mpc
v_turbulent = 12.0 # km/s
heII_abundance = 0.07 # He+/H+
electron_temps = np.array([7500.0])
electron_dens = np.logspace(1.,5.,10.) # 10 to 10^5 cm-3
HII_region_sizes = np.logspace(-1.,3.,10.) # 0.1 to 10^3 pc

# use H104 line flux to calculate number of HII regions
best_model = fit_model.fit_model(linefile,contfile,omega_region,distance,n_for_num=104,electron_temps=electron_temps,electron_dens=electron_dens,HII_region_sizes=HII_region_sizes,v_turbulent=v_turbulent,heII_abundance=heII_abundance,verbose=True,plotRSS=True,plotModel=True)






--> Whole Region (l > 0.1)


Ran 100 models in 00h 00m 7.23s
00h 00m 0.07s per model
5 models passed, 95 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 1.55985e-43
Electron temperature = 7500 K
Electron density = 46415.9 cm-3
HII region size = 0.1 pc
Number of HII_regions = 1455.58
Non-thermal spectral index = -0.602045
Mass of ionized gas = 874.435 Msun
Number of H-ionizing photons = 1.97655e+51 s-1
Star Formation Rate = 0.014409 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 7500 K
Electron density = 25 cm-3
HII region size = 131.856 pc
Non-thermal spectral index = -2.1369
Number of H-ionizing photons = 9.03066e+50 s-1
Star Formation Rate = 0.00658335 Msun/yr
=============================




Ran 100 models in 00h 00m 8.35s
00h 00m 0.08s per model
8 models passed, 92 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 1.56004e-43
Electron temperature = 6000 K
Electron density = 5994.84 cm-3
HII region size = 0.1 pc
Number of HII_regions = 115853
Non-thermal spectral index = -0.951147
Mass of ionized gas = 8988.98 Msun
Number of H-ionizing photons = 3.2757e+51 s-1
Star Formation Rate = 0.0238798 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 6000 K
Electron density = 25 cm-3
HII region size = 131.856 pc
Non-thermal spectral index = -2.5226
Number of H-ionizing photons = 1.12725e+51 s-1
Star Formation Rate = 0.00821768 Msun/yr
=============================

Ran 100 models in 00h 00m 11.75s
00h 00m 0.12s per model
3 models passed, 97 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 1.55978e-43
Electron temperature = 9000 K
Electron density = 46415.9 cm-3
HII region size = 0.1 pc
Number of HII_regions = 2101.88
Non-thermal spectral index = -0.847165
Mass of ionized gas = 1262.7 Msun
Number of H-ionizing photons = 2.3731e+51 s-1
Star Formation Rate = 0.0172999 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 9000 K
Electron density = 25 cm-3
HII region size = 131.856 pc
Non-thermal spectral index = -1.85291
Number of H-ionizing photons = 7.50856e+50 s-1
Star Formation Rate = 0.00547374 Msun/yr
=============================




--> West Component (l > 0.1)

Ran 100 models in 00h 00m 7.21s
00h 00m 0.07s per model
4 models passed, 96 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 9.82129e-45
Electron temperature = 7500 K
Electron density = 46415.9 cm-3
HII region size = 0.1 pc
Number of HII_regions = 366.134
Non-thermal spectral index = -0.538078
Mass of ionized gas = 219.954 Msun
Number of H-ionizing photons = 4.97177e+50 s-1
Star Formation Rate = 0.00362442 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 7500 K
Electron density = 70 cm-3
HII region size = 42.5564 pc
Non-thermal spectral index = -2.86678
Number of H-ionizing photons = 2.38028e+50 s-1
Star Formation Rate = 0.00173523 Msun/yr
=============================


Ran 100 models in 00h 00m 7.23s
00h 00m 0.07s per model
7 models passed, 93 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 9.82302e-45
Electron temperature = 6000 K
Electron density = 16681 cm-3
HII region size = 0.1 pc
Number of HII_regions = 3055.72
Non-thermal spectral index = -0.627798
Mass of ionized gas = 659.722 Msun
Number of H-ionizing photons = 6.68957e+50 s-1
Star Formation Rate = 0.0048767 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 6000 K
Electron density = 70 cm-3
HII region size = 42.5564 pc
Non-thermal spectral index = -4.0287
Number of H-ionizing photons = 2.97119e+50 s-1
Star Formation Rate = 0.002166 Msun/yr
=============================


Ran 100 models in 00h 00m 7.24s
00h 00m 0.07s per model
2 models passed, 98 models failed
================================
Multile Compact HII Region Model
================================
Best model RSS: 9.8205e-45
Electron temperature = 9000 K
Electron density = 46415.9 cm-3
HII region size = 0.1 pc
Number of HII_regions = 528.704
Non-thermal spectral index = -0.800156
Mass of ionized gas = 317.617 Msun
Number of H-ionizing photons = 5.96926e+50 s-1
Star Formation Rate = 0.00435159 Msun/yr
================================
=============================
Spontaneous Emissionn Model
=============================
Electron temperature = 9000 K
Electron density = 70 cm-3
HII region size = 42.5564 pc
Non-thermal spectral index = -2.40772
Number of H-ionizing photons = 1.97909e+50 s-1
Star Formation Rate = 0.00144276 Msun/yr
=============================





