# DRAGON
Diffusion Reacceleration and Advection of Galactic cosmic rays: an Open New code

## CREDITS

We acknowledge here the use of external routines/table:
* **dmspec.F** Routines for calculating the annihilation spectrum from [DarkSUSY](http://www.darksusy.org) package, to be cited as [Gondolo et al., 2004](http://arxiv.org/abs/astro-ph/0406204)
* **crn6.f**, **WNEWTR_FUNC_aws.f**, **YIELDX_011000_imos.f**, **xsec.cc** Routines to compute the nuclear cross-sections, from [GALPROP](http://galprop.stanford.edu) package (v54), to be cited as
** [Strong and Moskalenko, 2001](http://adsabs.harvard.edu/abs/2001AdSpR..27..717S)
** [Mashnik et al., 1998](http://adsabs.harvard.edu/abs/1998nucl.th..12071M)
** [Mashnik et al., 2004](http://adsabs.harvard.edu/abs/2004AdSpR..34.1288M)
** [Moskalenko et al., 2003](http://adsabs.harvard.edu/abs/2003ICRC....4.1969M)
** [Moskalenko et al., 2001](http://adsabs.harvard.edu/abs/2001ICRC....5.1836M)
** [Moskalenko et al., 2003b](http://adsabs.harvard.edu/abs/2003ApJ...586.1050M)
** [Barashenkov and Polanski, 1994](http://lt-jds.jinr.ru/record/5725?ln=en)
* galprop_*.dat Tables with nuclear cross-sections, from GALPROP package (v54)
* **MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz** The ISRF model used for energy losses, from [GALPROP](http://galprop.stanford.edu) package (v54)
* webber_xsec.dat Tabulated spallation cross sections, from Webber et al., ApJS, 144 (2003)
* webber_xsec_total.dat Tabulated inelastic cross sections, from Webber et al., ApJS, 144 (2003)
* **cparamlib** Library for secondary production in pp interactions, from [cparamlib](https://github.com/niklask/cparamlib) repository, to be cited as [Kamae, et al., 2007](https://arxiv.org/abs/astro-ph/0605581)
* **tinyxml** A C++ XML parser, from [here](http://www.grinninglizard.com/tinyxml)