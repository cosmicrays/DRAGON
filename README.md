# DRAGON
Diffusion Reacceleration and Advection of Galactic cosmic rays: an Open New code

## CREDITS
We acknowledge here the use of external routines/table:
* **dmspec.F** Routines for calculating the annihilation spectrum from [DarkSUSY](http://www.darksusy.org) package [ref.](http://arxiv.org/abs/astro-ph/0406204)
* crn6.f, WNEWTR_FUNC_aws.f, YIELDX_011000_imos.f Routines for evaluating the nuclear cross-sections, from GALPROP package (v54)
* galprop_*.dat Tables with nuclear cross-sections, from GALPROP package (v54)
* MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz The ISRF model used for energy losses, from GALPROP package (v54)
* webber_xsec.dat Tabulated spallation cross sections, from Webber et al., ApJS, 144 (2003)
* webber_xsec_total.dat Tabulated inelastic cross sections, from Webber et al., ApJS, 144 (2003)
* cparamlib Library for secondary production in pp interactions, from cparamlib repository
* tinyxml A C++ XML parser, from here