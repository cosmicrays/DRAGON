# DRAGON
Diffusion Reacceleration and Advection of Galactic cosmic rays: an Open New code

## Introduction
 
DRAGON adopts a second-order Cranck-Nicholson scheme with Operator Splitting and time overrelaxation to solve the diffusion equation. This provides fast a solution that is enough accurate for the average user. Occasionally, users may want to have very accurate solutions to their problem. To enable this feature, users may get close to the accurate solution by using the fast method, and then switch to a more accurate solution scheme, featuring the Alternating-Direction-Implicit (ADI) Cranck-Nicholson scheme.
 
Some parts of DRAGON are built following [GALPROP](http://galprop.stanford.edu/), v50p. The first reason is that it is a waste of time to reimplement standard parts, like energy losses, in which nothing new has to be found. The second reason is that it is essential to be able to compare our predictions with that of the Galprop code, and this  can be done only by following the details of its implementation. Therefore, we kept in the code some features and models used in Galprop, like nuclear cross-sections, the gas distribution, the convergence technique. However, each of these models is accompanied by other models, which can be selected by setting the appropriate switch. This is done very easily using the well known C++ structure of abstract/derived classes. The code is then very flexible and easy to manage and to modify or update. 
 
 The code was built having in mind a few motivations:
 - in order to find good propagation models one needs to run the code thousand times. Therefore we wrote the code aiming at performances and with an efficient memory management.
 - we wanted to propagate DM originated cosmic rays. Therefore we wrote DRAGON as a library that can be coupled to, e.g., DarkSUSY.
 - from the physics point of view, we wanted to have a position dependent diffusion model, which requires a substantial modification of the discretization scheme.
 
## Installation

Get the latest version of the code here and untar it.

To install DRAGON you need to download and install the [GSL](http://www.gnu.org/software/gsl/) libraries and [CFITSIO](http://heasarc.gsfc.nasa.gov/fitsio/) library first.

It also optionally requires [DarkSUSY](http://www.physto.se/~edsjo/darksusy/), if user wants to propagate DM originated cosmic rays.

Before installing the code launch the script to initialize installation tools:

`./start.sh`

Configure the code, a typical command line is:

`./configure --with-cfitsio=$CFITSIO_DIR --with-numcpu=2`
 
where `$CFITSIO_DIR` is the path of your cfitsio library and `NUMCPU` is the machine core number.

The default installation path is in the same folder as the source code is (the program automatically creates the `bin/` and `lib/` subfolders). It can be set via `--prefix=<NEW_INSTALLATION_PATH>`

Finally create the executable:

`make`

Run the example models in the examples/ directory:

`./DRAGON examples/run_2D.xml` 

or

`./DRAGON examples/run_3D.xml` 

A detailed description of the input XML file and some examples can be found [here](https://indico.desy.de/indico/event/10336/contribution/2/material/slides/0.pdf).
 
Please let us know for any problem with installation!

## CREDITS

We acknowledge here the use of external routines/table:
* **dmspec.F** Routines for calculating the annihilation spectrum from [DarkSUSY](http://www.darksusy.org) package, to be cited as [Gondolo et al., 2004](http://arxiv.org/abs/astro-ph/0406204)
* **MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz** The ISRF model used for the energy losses, to be downloaded from the [GALPROP](http://galprop.stanford.edu) package (v54) and to be cited as [Porter and Strong, 2008](http://adsabs.harvard.edu/abs/2008AAS...212.1810P)
* **webber_xsec.dat** Tabulated spallation cross sections, to be cited as [Webber et al., 2003](http://adsabs.harvard.edu/abs/2003ApJS..144..153W)
* **webber_xsec_total.dat** Tabulated inelastic cross sections, to be cited as [Webber et al., 2003](http://adsabs.harvard.edu/abs/2003ApJS..144..153W)
* **cparamlib** Library for secondary production in pp interactions from the [cparamlib](https://github.com/niklask/cparamlib) repository, to be cited as [Kamae, et al., 2007](https://arxiv.org/abs/astro-ph/0605581)
* **tinyxml** A C++ XML parser, from [here](http://www.grinninglizard.com/tinyxml)

## Cross-section network

We include the routines and data tables taken from the public version of [GALPROP](http://galprop.stanford.edu). In more detail, the material included in our code contains: 1) the nuclear reaction network, built using the Nuclear Data Sheets; 2) the isotopic cross section database built using the T16 Los Alamos compilation [Mashnik et al., 1998](http://adsabs.harvard.edu/abs/1998nucl.th..12071M) and the CEM2k and LAQGSM codes [Mashnik et al., 2004](http://adsabs.harvard.edu/abs/2004AdSpR..34.1288M); 3) fits to some particular channels of isotopic production cross section [Moskalenko et al., 2001](http://adsabs.harvard.edu/abs/2001ICRC....5.1836M), [Moskalenko et al., 2003](http://adsabs.harvard.edu/abs/2003ICRC....4.1969M), [Moskalenko et al., 2003b](http://adsabs.harvard.edu/abs/2003ApJ...586.1050M); 4) phenomenological approximations adapted from [Webber et al., 1990](http://adsabs.harvard.edu/abs/1990PhRvC..41..566W) and [Silberberg et al., 1998](http://adsabs.harvard.edu/abs/1998ApJ...501..911S); 5) inelastic cross section database adapted from [Barashenkov and Polanski, 1994](http://lt-jds.jinr.ru/record/5725?ln=en)

Specific files from [GALPROP](http://galprop.stanford.edu) are listed below.

* **crn6.f**, **WNEWTR_FUNC_aws.f**, **YIELDX_011000_imos.f**, **xsec.cc** Routines to compute the nuclear cross-sections, from [GALPROP](http://galprop.stanford.edu) package (v54)
* **galprop_*.dat** Tables with nuclear cross-sections, to be downloaded from the [GALPROP](http://galprop.stanford.edu) package (v54)
