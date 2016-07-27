/**
 * @file    dragon.cc
 * @author  Luca Maccione
 * @email   luca.maccione@lmu.de
 * @brief   Implementation of the DRAGON class. See the .h file
 */

#include "dragon.h"
#include "grid.h"
#include "nucleilist.h"
#include "galaxy.h"
#include "particle.h"
#include "xsec.h"
#include "crevolutor.h"
#include "input.h"
#include "utilities.h"
#include "sources.h"
#include "spectrum.h"
#include <iomanip>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <sys/stat.h> //sk for mkdir
#include <sys/types.h>

using namespace std;

DRAGON::DRAGON(Input* inputStructure_ /**< Wrapper to user input. */) {
  
  inputStructure = inputStructure_;
  
  const int feedback = inputStructure->feedback;
  
  status = 0;
  
  norm = 1.0;
  normel = 0.0;
  normextra = 0.0;
  normDM = 0.0;
  
  if (feedback>0) cout<<"Getting the nuclei list ..."<<endl;
  nucleiList = new TNucleiList(inputStructure_);
  if (feedback>0) cout<<"... done!"<<endl<<endl;
  
  if (feedback>0) cout<<"Initializing the Galaxy ..."<<endl;
  gal = new Galaxy(inputStructure_, nucleiList);
  if (feedback >0) cout<<"... done!"<<endl<<endl;
  
  if (feedback>0) cout<<"Initializing the Algorithm ..."<<endl;
  if (inputStructure->gridtype == "2D") {
    if (inputStructure->OpSplit) alg.push_back( new TCREvolutor(gal) );
    if (inputStructure->ADI) alg.push_back( new TCREvolutorADI(gal) );
  }
  else {
    if (inputStructure->OpSplit) alg.push_back( new TCREvolutor3D(gal) );
    if (inputStructure->ADI) alg.push_back( new TCREvolutor3DADI(gal) );
  }
  if (feedback>0) cout<<"... done!"<<endl<<endl;
  
  if (inputStructure->fullstore) {
    string fullstoreFilename = make_filename(inputStructure->filename.c_str(), ".fits.gz");
    if (fits_create_file(&output_ptr, fullstoreFilename.c_str(), &status)) fits_report_error(stderr, status);
  } else output_ptr = NULL;
  
  if (inputStructure->partialstore) {
    string partialstoreFilename = make_filename(inputStructure->filename.c_str(), "_spectrum.fits.gz");
    if (fits_create_file(&output_ptr_sp, partialstoreFilename.c_str(), &status)) fits_report_error(stderr, status);
  } else output_ptr_sp = NULL;
  
  if (inputStructure->asciistore) {
    string asciiFilename = make_filename(inputStructure->filename.c_str(), ".txt");
  }
}


DRAGON::~DRAGON() {
   
  if (gal)
    delete gal;
  if (nucleiList)
    delete nucleiList;
   
  for (vector<TCREvolutorBasis*>::iterator i = alg.begin(); i != alg.end(); ++i)
    delete *i;
  alg.clear();
  
  for (vector<TParticle*>::iterator i = particles.begin(); i != particles.end(); ++i)
    delete *i;
  particles.clear();
  
  if (output_ptr) {
    if (fits_close_file(output_ptr, &status))
      fits_report_error(stderr, status);
  }
  if (output_ptr_sp) {
    if (fits_close_file(output_ptr_sp, &status))
      fits_report_error(stderr, status);
  }
}


TParticle* DRAGON::FindParticle(int uid, bool sec, int isDM, int isextra) {
  for (vector<TParticle*>::reverse_iterator ripart = particles.rbegin(); ripart != particles.rend(); ++ripart) {
    if ( (*ripart)->GetUid() == uid && (*ripart)->GetIsSec() == sec && (*ripart)->IsDM() == isDM && (*ripart)->IsExtra() == isextra )
      return *ripart;
  }
   
  return NULL;
}


vector<TParticle*> DRAGON::FindSpecies(int Z) {
  vector<TParticle*> result;
  for (vector<TParticle*>::reverse_iterator ripart = particles.rbegin(); ripart != particles.rend(); ++ripart) {
    if ( (*ripart)->GetZ() == Z )
      result.push_back(*ripart);
  }
   
  return result;
}


TParticle* DRAGON::CreateParticle(const string& grtype, int A_ /**< Mass number */, int Z_ /**< Charge */, Galaxy* gal /**< A model for the galaxy */, Input* in /**< User input */, bool issec, vector<TXSecBase*> xsecmodel, TNucleiList* l, int K_electron, bool isDM, bool isextra, bool isTPP_) {
   
  if (grtype == "2D") return new TParticle2D(A_, Z_, gal, in, issec, xsecmodel, l, K_electron, isDM, isextra, isTPP_);
  if (grtype == "3D") return new TParticle3D(A_, Z_, gal, in, issec, xsecmodel, l, K_electron, isDM, isextra, isTPP_);
   
  return NULL;
   
}

//*****************************************************************************************************************************************
//******************************************** The main method of DRAGON ******************************************************************
//*****************************************************************************************************************************************

void DRAGON::Run() {
  
  const int feedback = inputStructure->feedback;
  
  if (feedback >0) cout<<"Loading cross sections ..."<<endl;
  
  vector<TXSecBase*> xsecmodel;
  xsecmodel.push_back(new TGalpropXSec(gal->GetCoordinates()));
  
  if (feedback >1) cout<<"Galprop cross sections loaded"<<endl;
  
  if (inputStructure->spallationxsec == Webber03) 
    xsecmodel.push_back(new TWebber03());
  if (inputStructure->spallationxsec == Fluka) {
    if (feedback >1) cout << "Fluka cross sections " << endl;
    xsecmodel.push_back(new FlukaXSec(gal->GetCoordinates()));
  }
  
  if (feedback >0) cout<<"... done!"<<endl<<endl;
  
  vector<int> list_nuc = nucleiList->GetList();
  
  if (feedback >0) cout<<"Building the spallation network ..."<<endl;
  
  TSpallationNetwork* spallnet = new TSpallationNetwork(gal->GetCoordinates(), inputStructure, xsecmodel, list_nuc);
  
  if (feedback>0) cout<<"... done!"<<endl<<endl;
  
  const string grty = inputStructure->gridtype;
  
  //-- Loop over particle list!--------------------------------------------------------------------------------

  for (vector<int>::iterator inuc = list_nuc.begin(); inuc != list_nuc.end(); ++inuc) {
    
    int A = -1000;
    int Z = -1000;
    Utility::id_nuc(*inuc, A, Z); /**< Compute the ID of the nucleus. */
    
    if (feedback >0) cout << "Propagating nucleus with A = " << A << " Z = " << Z << endl;
    
    if (*inuc > 1000) { // nuclei
      
      DECMODE decay_mode = nucleiList->GetDecayMode(*inuc);
      double life = nucleiList->GetLifeTime(*inuc);
      
      if (feedback>1){
	if (decay_mode != STABLE) cout<<"Lifetime = "<<life*1e6<<" yr"<<endl;
	else cout<<"Stable nucleus"<<endl;
	if (decay_mode == EC || decay_mode == ECBM || decay_mode == ECBP) cout<<"Nucleus may attach an electron and decay via EC"<<endl;
      }
         
      if ((decay_mode != EC && decay_mode != ECBM && decay_mode != ECBP) || ((life*1e6 < t_half_limit) && !(inputStructure->Kcapture))) { //modified
	
	particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, -1, false, false)); /**< Creates the nucleus and add it to the particle vector. */
        //cout << "Evolve" << endl; 
	particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */

	if (*inuc == 1001) {
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, true, xsecmodel, nucleiList, -1, false, false)); /**< Creates the nucleus and add it to the particle vector. */
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
	}
	
      }
      else { /* if the decay mode is  K_capture and the lifetime is long enough */
            
	if (feedback>0) cout<<"Since the lifetime is long enough, the nucleus is propagated twice"<<endl;
	if (feedback>0) cout<<"Propagating naked nucleus"<<endl;
	particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, 0, false, false)); // naked nucleus. K_electron = 0
	particles.back()->Evolve(particles, alg, spallnet, xsecmodel);
            
	//	prev_uid = particles.back()->GetUid();
            
	if (feedback>0) cout<<"Propagating nucleus with an attached electron"<<endl;
	particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, 1, false, false)); // K_electron = 1
	particles.back()->Evolve(particles, alg, spallnet, xsecmodel);
      }
    }
    else {
         
      if (*inuc == -999) {
	// Antiprotons
            
	if (inputStructure->prop_secap) {
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, -1, false, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
               
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, true, xsecmodel, nucleiList, -1, false, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
               
	}
	if (inputStructure->prop_DMap) {
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, -1, true, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
               
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, true, xsecmodel, nucleiList, -1, true, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
               
	}
            
      }
      else if (*inuc == -998) {
	// Antideuterons
            
	if (inputStructure->prop_secdeuteron) {
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, -1, false, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
               
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, true, xsecmodel, nucleiList, -1, false, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
               
	}
	if (inputStructure->prop_DMdeuteron) {
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, -1, true, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
               
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, true, xsecmodel, nucleiList, -1, true, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
               
	}
      }
      else {
	// Leptons
            
	if (inputStructure->prop_extracomp && *inuc == -1000) {
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, -1, false, true));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
	}
	if (inputStructure->prop_lep) {
               
	  // Secondary electron or positron
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, true, xsecmodel, nucleiList, -1, false, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
               
	  //TPP: positrons from TPP
	  if (*inuc == 1000 && inputStructure->prop_TPP) {
	    if (inputStructure->feedback >0) cout << "TPP positrons" << endl;
	    particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, true, xsecmodel, nucleiList, -1, false,false,true));
	    if (inputStructure->feedback >0) cout << "Evolving TPP positrons " << endl;
	    particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
	  }
               
	  // Primary electrons
	  if (*inuc == -1000) {
	    particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, -1, false, false));
	    particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
	  }
	}
            
	if (inputStructure->prop_DMel && *inuc == -1000 ) {
	  if (inputStructure->feedback >0) cout << "DM electrons" << endl;
	  particles.push_back(CreateParticle(grty, A, Z, gal, inputStructure, false, xsecmodel, nucleiList, -1, true, false));
	  particles.back()->Evolve(particles, alg, spallnet, xsecmodel); /**< Evolve the nucleus. */
	}
      }
    }
    if (inputStructure->feedback >0) cout<<"... propagation done!"<<endl<<endl;
  }
  //-- End of loop over particle list!-----------------------------------------------------------------

  //-- Second Run!--------DG24.10.2013-----------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  //-- The second run is useful to refine the computation of very heavy nuclei. In this way the beta- decays
  //-- of nuclei with lower Z into nuclei with higher Z are taken into account properly!
  //-- The second run is useless for B/C protons and leptons; 
  //-- since it makes the code much slower so it is DISABLED in the default xml
  //-- In order to enable it, put the <DoubleRun> flag after the <NuclearChain> block!
  //---------------------------------------------------------------------------------------------------

  if (inputStructure->DoubleRun == true) {

    if (inputStructure->feedback >1) cout<<endl<<endl<<" ******* STARTING SECOND ITERATION ******* "<<endl<<endl;

    //int count_nuc = 0;

    for (vector<TParticle*>::iterator i_current_part = particles.begin(); i_current_part != particles.end()-1; ++i_current_part) { 

      int A = (*i_current_part)->GetA();
      int Z = (*i_current_part)->GetZ();
      int uid = (*i_current_part)->GetUid();	
      if (feedback>1){
	cout << endl << "--SECOND ITERATION--SECOND ITERATION-----" << endl;
	cout         << "-----------------------------------------" << endl;
	cout         << "- Starting with nucleus A = " << A << " Z = " << Z << endl;
	cout         << "- Starting propagation..." << endl;
	cout         << "-----------------------------------------" << endl;
	cout         << "--SECOND ITERATION--SECOND ITERATION-----" << endl;
      }
      (*i_current_part)->Evolve(particles, alg, spallnet, xsecmodel, true);

    } 
  }  
  //-- End of second loop over particle vector!-----------------------------------------
  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------

  /**
   * Normalization part. If it is DM then just convert units, otherwise find protons and electrons and normalize to observed flux.
   */
   
  //DG17.10.2013
  //	
  // BEFORE normalization DM particles are number densities per energy unit: 1/GeV/kpc^3
  // ------
  // applying kpc^-3 ---> 1/cm^3/GeV
  // applying C/4pi  ---> 1/cm^2/s/sr/GeV
  // applying 1.e4   ---> 1/m^2/s/sr/GeV
  // -----
  // AFTER normalization DM particles have DRAGON units for the flux: 1/m^2/s/sr/GeV

  normDM = C * pow(kpc,-3.) /4./Pi*1.e4;    //C is in cm/s
  if (inputStructure->propDM && feedback >0) cout<<"Calculate the DM normalization = "<<normDM<<endl;
   
  TParticle* protons = FindParticle(1001, false) ;
  if (protons) norm = protons->FindNormalization(inputStructure->sp_ref_rig_norm, inputStructure->spect_norm);
  if (feedback>0) cout<<"Calculate the proton normalization = "<<norm<<endl;
  
  TParticle* electrons = FindParticle( -1000, false);
  if (electrons) normel = electrons->FindNormalization(inputStructure->sp_ref_rig_el, inputStructure->spect_norm_el);   
  if (feedback>0) cout<<"Calculate the electron normalization = "<<normel<<endl;
  
  TParticle* electronsExtra = FindParticle( -1000, false, false, true);
  if (electronsExtra) 
    normextra = electronsExtra->FindNormalization(inputStructure->sp_ref_rig_el_extra, inputStructure->spect_norm_el_extra);
  if (feedback>0) cout<<"Calculate the extra-component normalization = "<<electronsExtra<<endl;

  if (inputStructure->DontNormalize) {
	norm = 1.;
	normel = 1.;
	normextra = 1.;
  }

  //cout << "Finding normalization for the extra component" << endl;
  //cout << in->sp_ref_rig_el_extra << " <- energy || norm -> " << in->spect_norm_el_extra << endl	;
   
  //DG14.01.2014	--- injected energy rate --- --- --- --- --- --- --- --- --- 
  //updated 16.01.2014
  
  vector<double> Ek;
  vector<double> spectrum_;
  vector<double> source_distribution;

  Ek.clear();
  spectrum_.clear();
  source_distribution.clear();

  /**
   * Print nuclear densities and/or spectra in output file.
   * Lighter nuclei are printed first.
   */
  if (feedback>0) cout<<endl<<endl;
  
  for (vector<TParticle*>::reverse_iterator ripart = particles.rbegin(); ripart != particles.rend(); ++ripart) {
    if ((*ripart)->IsDM()) {
      if (output_ptr) (*ripart)->Print(output_ptr, normDM);
      if (output_ptr_sp) (*ripart)->PrintSpectrum(output_ptr_sp, normDM);
      
    }
    else {
      if ((*ripart)->IsExtra()) {
	if (output_ptr) (*ripart)->Print(output_ptr, normextra);
	if (output_ptr_sp) (*ripart)->PrintSpectrum(output_ptr_sp, normextra);
	
      }
      else {
	if ( ( (*ripart)->GetUid() == -1000 && !(*ripart)->GetIsSec()) || (*ripart)->IsTPP() ) {
	  // primary electrons and secondary positrons from TPP are normalized to electrons
	  
	  if (output_ptr) (*ripart)->Print(output_ptr, normel);
	  if (output_ptr_sp) (*ripart)->PrintSpectrum(output_ptr_sp, normel);
	  
	}
	else {
	  if (output_ptr) (*ripart)->Print(output_ptr, norm);
	  if (output_ptr_sp) (*ripart)->PrintSpectrum(output_ptr_sp, norm);
	}
      }
    }
  }
  
  /* ASCII output */
  vector<double> e_GeV = gal->GetCoordinates()->GetEk();
  char ascii_filename[1000];
  sprintf(ascii_filename,"output/%s.txt",inputStructure->filename.c_str());
  cout<<"Writing ASCII output file: "<<ascii_filename<<endl;
  ofstream outfile(ascii_filename, ios::out);
  outfile << "Energy [GeV] ";
  for (vector<TParticle*>::reverse_iterator current_part = particles.rbegin(); current_part != particles.rend(); ++current_part) {
    if ( (*current_part)->IsDM() && (*current_part)->GetUid() == -1000 )
      outfile << "DM_e-        " ;
    if ( (*current_part)->IsDM() && (*current_part)->GetUid() == 1000 )
      outfile << "DM_e+        " ;
    if ( (*current_part)->IsDM() && (*current_part)->GetUid() == -999 )
      outfile << "DM_pbar      " ;
    if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -1000 && !(*current_part)->GetIsSec() &&  !(*current_part)->IsExtra() )
      outfile << "Pri_e-       " ;
    if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -1000 &&  (*current_part)->GetIsSec() &&  !(*current_part)->IsExtra() )
      outfile << "Sec_e-       " ;
    if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() ==  1000 &&  (*current_part)->GetIsSec() &&  !(*current_part)->IsExtra() )
      outfile << "Sec_e+       " ;
    if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -1000 &&  !(*current_part)->GetIsSec() &&   (*current_part)->IsExtra() )
      outfile << "Extra_e+     " ;
    if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -999 && !(*current_part)->GetIsSec()  &&  !(*current_part)->IsExtra() )
      outfile << "Sec_pbar     " ;
    if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == -999 &&  (*current_part)->GetIsSec()  &&  !(*current_part)->IsExtra() )
      outfile << "Ter_pbar     " ;
    if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == 1001 && !(*current_part)->GetIsSec()  &&  !(*current_part)->IsExtra() )
      outfile << "Pri_p        " ;
    if ( !((*current_part)->IsDM()) && (*current_part)->GetUid() == 1001 &&  (*current_part)->GetIsSec()  &&  !(*current_part)->IsExtra() )
      outfile << "Sec_p        " ;
    if ( (*current_part)->GetUid() > 1001)
      outfile << "NUC_" << (*current_part)->GetUid() <<"     ";
  }
  outfile << endl;
  for (int ie=0; ie<e_GeV.size(); ie++) {
    //
    outfile << scientific << setprecision(6) << e_GeV[ie] << " ";
    //
    for (vector<TParticle*>::reverse_iterator current_part = particles.rbegin(); current_part != particles.rend(); ++current_part) {
      //
      double flux = (*current_part)->GetFluxAtSunPosition(ie);
      //
      if ((*current_part)->IsDM())
	flux *= normDM;
      else {
	if ((*current_part)->IsExtra())			
	  flux *= normextra;
	else {
	  if ( ( (*current_part)->GetUid() == -1000 && !(*current_part)->GetIsSec()) || (*current_part)->IsTPP() )
	    flux *= normel;
	  else
	    flux *= norm;
	}
      }	 		
      outfile << flux << " ";
    }	
    outfile << endl;
  }
  outfile.close(); 

  delete spallnet;
  spallnet = NULL;
  while (xsecmodel.size()) {
    delete xsecmodel.back();
    xsecmodel.pop_back();
  }
   
  return ;
}

void DRAGON::Print() {
   
  double Rmax = inputStructure->Rmax;

  unsigned int cr_irsun;
  if (inputStructure->gridtype == "2D")
    cr_irsun = (unsigned int) (inputStructure->robs/Rmax*(double)(inputStructure->numr-1)); /**< Radial Sun position. */
  else
    cr_irsun = (unsigned int) ((inputStructure->robs+Rmax)/(2.0*Rmax)*(double)(inputStructure->numx-1)); /**< Radial Sun position. */
  unsigned int cr_ixsun;
  unsigned int cr_iysun;
  if (inputStructure->gridtype == "3D") { //DG28.11.2013 added ixsun and iyzun in the 3D case
    cr_ixsun = (unsigned int) ((inputStructure->xobs+Rmax)/(2.0*Rmax)*(double)(inputStructure->numx-1));
    cr_iysun = (unsigned int) ((inputStructure->yobs+Rmax)/(2.0*Rmax)*(double)(inputStructure->numy-1));
  }   
  unsigned int cr_izsun = (unsigned int) ((inputStructure->zobs+inputStructure->zmax)/(2.0*inputStructure->zmax)*(double)(inputStructure->numz-1)), /**< Vertical Sun position. */

    cr_gas_model = (unsigned int)inputStructure->gas_model,
    cr_SNR_model = (unsigned int)inputStructure->SNR_model;
   
  double cr_rmin = inputStructure->Rmin,
    cr_rmax = Rmax,
    cr_zmin = -inputStructure->zmax,
    cr_zmax =  inputStructure->zmax,
    cr_rB = rB,
    cr_zt = inputStructure->zt,
    cr_zr = zr,
    cr_Bh = Bh,
    cr_D0 = inputStructure->D0*kpc*kpc/Myr,
    cr_D_ref_rig = inputStructure->D_ref_rig,
    cr_ind_diff = inputStructure->delta,
    cr_He_abundance = He_abundance,
    cr_sp_ref_rig = inputStructure->sp_ref_rig,
    cr_sp_ref_rig_norm = inputStructure->sp_ref_rig_norm,
    cr_spect_norm = inputStructure->spect_norm,
    cr_Ekmax = inputStructure->Ekmax,
    cr_Ekmin = inputStructure->Ekmin,
    cr_Ekfac = inputStructure->Ekfact,
    cr_u = u,
    //cr_tol = tolerance,
    //cr_alpha = alpha,
    cr_p = p,
    cr_vA = inputStructure->vAlfven*kpc/Myr/km,
    cr_v0 = inputStructure->v0*kpc/Myr/km,
    cr_dvdz = inputStructure->dvdz*kpc/Myr/km,
    cr_robs = inputStructure->robs,
    cr_zobs = inputStructure->zobs;

  double cr_xmin, cr_xmax, cr_ymin, cr_ymax, cr_xobs, cr_yobs; //DG28.11.2013
  if (inputStructure->gridtype == "3D") {
    cr_xmin = -inputStructure->Rmax;
    cr_xmax =  inputStructure->Rmax;
    cr_ymin = -inputStructure->Rmax;
    cr_ymax =  inputStructure->Rmax;
    cr_xobs =  inputStructure->xobs;
    cr_yobs =  inputStructure->yobs;
  }   
   
  double crmx = inputStructure->mx;
  double crtaudec = inputStructure->taudec;
  double crsigmav = inputStructure->sigmav;
  int crdmmode = inputStructure->dmmode;
  int crdmprof = inputStructure->dmprof;
   
  const long naxis = (inputStructure->gridtype == "2D") ? 3 : 4;
  int dimE = int(log(inputStructure->Ekmax/inputStructure->Ekmin)/log(inputStructure->Ekfact) + 1.9);
  long size_axes[naxis];
  size_axes[0] =dimE;
  if (inputStructure->gridtype == "2D") {
    size_axes[1] = inputStructure->numr;
    size_axes[2] = inputStructure->numz;
  }
  else {
    size_axes[1] = inputStructure->numx;
    size_axes[2] = inputStructure->numy;
    size_axes[3] = inputStructure->numz;
  }
  long nelements = 1;
  for (int i = 0; i < naxis; ++i) nelements *= size_axes[i];
   
  long fpixel = 1;
  int bitpix = FLOAT_IMG;
   
  /**
   * Create the first HDU, which only contains header data.
   */
   
  if (output_ptr) {
    if (fits_create_img(output_ptr, bitpix, naxis, size_axes, &status))
      fits_report_error(stderr, status);
    //
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "rmin",      &cr_rmin,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*)  "rmax",      &cr_rmax,            NULL, &status))
      fits_report_error(stderr, status);
    if (inputStructure->gridtype == "3D") {//DG28.11.2013
      if (fits_write_key(output_ptr, TDOUBLE, (char*) "xmin",      &cr_xmin,            NULL, &status))
	fits_report_error(stderr, status);
      if (fits_write_key(output_ptr, TDOUBLE, (char*) "xmax",      &cr_xmax,            NULL, &status))
	fits_report_error(stderr, status);
      if (fits_write_key(output_ptr, TDOUBLE, (char*) "ymin",      &cr_ymin,            NULL, &status))
	fits_report_error(stderr, status);
      if (fits_write_key(output_ptr, TDOUBLE, (char*) "ymax",      &cr_ymax,            NULL, &status))
	fits_report_error(stderr, status);
    }
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "zmin",      &cr_zmin,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "zmax",      &cr_zmax,            NULL, &status))
      fits_report_error(stderr, status);
    // 
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "robs",      &cr_robs,            NULL, &status))
      fits_report_error(stderr, status);
    if (inputStructure->gridtype == "3D") {//DG28.11.2013
      if (fits_write_key(output_ptr, TDOUBLE, (char*) "xobs",      &cr_xobs,            NULL, &status))
	fits_report_error(stderr, status);
      if (fits_write_key(output_ptr, TDOUBLE, (char*) "yobs",      &cr_yobs,            NULL, &status))
	fits_report_error(stderr, status);
    }
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "zobs",      &cr_zobs,            NULL, &status))
      fits_report_error(stderr, status);
    //
    if (fits_write_key(output_ptr, TUINT,   (char*) "irsun",     &cr_irsun,           NULL, &status))
      fits_report_error(stderr, status);
    if (inputStructure->gridtype == "3D") {//DG28.11.2013
      if (fits_write_key(output_ptr, TINT, (char*) "ixsun",      &cr_ixsun,            NULL, &status))
	fits_report_error(stderr, status);
      if (fits_write_key(output_ptr, TINT, (char*) "iysun",      &cr_iysun,            NULL, &status))
	fits_report_error(stderr, status);
    }
    if (fits_write_key(output_ptr, TUINT,   (char*) "izsun",     &cr_izsun,           NULL, &status))
      fits_report_error(stderr, status);
    //
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "u",         &cr_u,               NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "rB",        &cr_rB,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "zt",        &cr_zt,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "zr",        &cr_zr,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "Bh",        &cr_Bh,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "D0",        &cr_D0,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "Drefrig",   &cr_D_ref_rig,       NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "ind_diff",  &cr_ind_diff,        NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "He_ab",     &cr_He_abundance,    NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "sprefrig",  &cr_sp_ref_rig,      NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "sprignor",  &cr_sp_ref_rig_norm, NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "specnorm",  &cr_spect_norm,      NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "Ekmax",     &cr_Ekmax,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "Ekmin",     &cr_Ekmin,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "Ekin_fac",  &cr_Ekfac,           NULL, &status))
      fits_report_error(stderr, status);
    //    if (fits_write_key(output_ptr, TDOUBLE, (char*) "tol",       &cr_tol,             NULL, &status)) fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TUINT,   (char*) "gasmod",    &cr_gas_model,       NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TUINT,   (char*) "SNRmod",    &cr_SNR_model,       NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TINT,    (char*) "dimE",      &dimE,               NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TINT,    (char*) "dimr",      &(inputStructure->numr),         NULL, &status))
      fits_report_error(stderr, status);
    if (inputStructure->gridtype == "3D") {//DG28.11.2013
      if (fits_write_key(output_ptr, TINT, (char*) "dimx",      &(inputStructure->numx),            NULL, &status))
	fits_report_error(stderr, status);
      if (fits_write_key(output_ptr, TINT, (char*) "dimy",      &(inputStructure->numy),            NULL, &status))
	fits_report_error(stderr, status);
    }
    if (fits_write_key(output_ptr, TINT,    (char*) "dimz",      &(inputStructure->numz),         NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "vAlfven",   &cr_vA,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "v0_conv",   &cr_v0,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "dvdz_conv", &cr_dvdz,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "mx",        &crmx,               NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "taudec",    &crtaudec,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TDOUBLE, (char*) "sigmav",    &crsigmav,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TINT,    (char*) "dmmode",    &crdmmode,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr, TINT,    (char*) "dmprof",    &crdmprof,           NULL, &status))
      fits_report_error(stderr, status);
  }
  if (output_ptr_sp) {
    if (fits_create_img(output_ptr_sp, bitpix, 1, size_axes, &status))
      fits_report_error(stderr, status);
      
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "rmin",      &cr_rmin,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "rmax",      &cr_rmax,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zmin",      &cr_zmin,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zmax",      &cr_zmax,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "robs",      &cr_robs,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zobs",      &cr_zobs,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TUINT,   (char*) "irsun",     &cr_irsun,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TUINT,   (char*) "izsun",     &cr_izsun,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "u",         &cr_u,               NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "rB",        &cr_rB,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zt",        &cr_zt,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "zr",        &cr_zr,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Bh",        &cr_Bh,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "D0",        &cr_D0,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Drefrig",   &cr_D_ref_rig,       NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "ind_diff",  &cr_ind_diff,        NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "He_ab",     &cr_He_abundance,    NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "sprefrig",  &cr_sp_ref_rig,      NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "sprignor",  &cr_sp_ref_rig_norm, NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "specnorm",  &cr_spect_norm,      NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Ekmax",     &cr_Ekmax,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Ekmin",     &cr_Ekmin,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "Ekin_fac",  &cr_Ekfac,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TUINT,   (char*) "gasmod",    &cr_gas_model,       NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TUINT,   (char*) "SNRmod",    &cr_SNR_model,       NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TINT,    (char*) "dimE",      &dimE,               NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TINT,    (char*) "dimr",      &(inputStructure->numr),         NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TINT,    (char*) "dimz",      &(inputStructure->numz),         NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "vAlfven",   &cr_vA,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "v0_conv",   &cr_v0,              NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "dvdz_conv", &cr_dvdz,            NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "mx",        &crmx,               NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "taudec",    &crtaudec,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TDOUBLE, (char*) "sigmav",    &crsigmav,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TINT,    (char*) "dmmode",    &crdmmode,           NULL, &status))
      fits_report_error(stderr, status);
    if (fits_write_key(output_ptr_sp, TINT,    (char*) "dmprof",    &crdmprof,           NULL, &status))
      fits_report_error(stderr, status);
  }
   
   
  if (output_ptr) {
    float* array = new float[nelements]();
    if (fits_write_img(output_ptr, TFLOAT, fpixel, nelements, array, &status))
      fits_report_error(stderr, status);
    delete [] array;
  }
  if (output_ptr_sp) {
    float* array = new float[size_axes[0]]();
    if (fits_write_img(output_ptr_sp, TFLOAT, fpixel, size_axes[0], array, &status))
      fits_report_error(stderr, status);
    delete [] array;
  }
}


