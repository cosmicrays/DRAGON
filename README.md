# DRAGON
Diffusion Reacceleration and Advection of Galactic cosmic rays: an Open New code

## CREDITS

We acknowledge here the use of external routines/table:
* **dmspec.F** Routines for calculating the annihilation spectrum from [DarkSUSY](http://www.darksusy.org) package, to be cited as [Gondolo et al., 2004](http://arxiv.org/abs/astro-ph/0406204)
* **MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz** The ISRF model used for the energy losses, to be downloaded from the [GALPROP](http://galprop.stanford.edu) package (v54)
* **webber_xsec.dat** Tabulated spallation cross sections, to be cited as [Webber et al., 2003](http://adsabs.harvard.edu/abs/2003ApJS..144..153W)
* **webber_xsec_total.dat** Tabulated inelastic cross sections, to be cited as [Webber et al., 2003](http://adsabs.harvard.edu/abs/2003ApJS..144..153W)
* **cparamlib** Library for secondary production in pp interactions from the [cparamlib](https://github.com/niklask/cparamlib) repository, to be cited as [Kamae, et al., 2007](https://arxiv.org/abs/astro-ph/0605581)
* **tinyxml** A C++ XML parser, from [here](http://www.grinninglizard.com/tinyxml)

## Cross-section network

We included the routines and data tables taken from the public version of [GALPROP](http://galprop.stanford.edu). In more detail, the material included in our code contains: 1) the nuclear reaction network, built using the Nuclear Data Sheets; 2) the isotopic cross section database built using the T16 Los Alamos compilation [Mashnik et al., 1998](http://adsabs.harvard.edu/abs/1998nucl.th..12071M) and the CEM2k and LAQGSM codes [Mashnik et al., 2004](http://adsabs.harvard.edu/abs/2004AdSpR..34.1288M); 3) fits to some particular channels of isotopic production cross section [Moskalenko et al., 2001](http://adsabs.harvard.edu/abs/2001ICRC....5.1836M), [Moskalenko et al., 2003](http://adsabs.harvard.edu/abs/2003ICRC....4.1969M), [Moskalenko et al., 2003b](http://adsabs.harvard.edu/abs/2003ApJ...586.1050M); 4) phenomenological approximations adapted from [Webber et al., 1990](http://adsabs.harvard.edu/abs/1990PhRvC..41..566W) and [Silberberg et al., 1998](http://adsabs.harvard.edu/abs/1998ApJ...501..911S); 5) inelastic cross section database adapted from [Barashenkov and Polanski, 1994](http://lt-jds.jinr.ru/record/5725?ln=en)

* **crn6.f**, **WNEWTR_FUNC_aws.f**, **YIELDX_011000_imos.f**, **xsec.cc** Routines to compute the nuclear cross-sections, from [GALPROP](http://galprop.stanford.edu) package (v54)
* **galprop_*.dat** Tables with nuclear cross-sections, to be downloaded from the [GALPROP](http://galprop.stanford.edu) package (v54)