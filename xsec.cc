/*
 * @file xsec.cc
 * @author Luca Maccione, Daniele Gaggero
 * @email luca.maccione@desy.de
 * @email daniele.gaggero@sissa.it
 * @brief All the classes related to the cross sections are implemented.
 */

#include "xsec.h"
#include "grid.h"
#include "nucleilist.h"
#include "utilities.h"
#include "kamae.h"
#include "input.h"
#include <algorithm>
#include <fstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h> 
#include <gsl/gsl_vector.h>

using namespace std;
//ofstream outf("spallation_xsec.dat", ios::out);

const double TSpallationNetwork::ethr = 5.e1;

TSpallationNetwork::TSpallationNetwork(TGrid* co, Input* in, vector<TXSecBase*> xsecmodel, vector<int>& nuclei) {

  /* the relevant structures of TSpallationNetwork to be filled are

     1) for nucleus --> nucleus
     
     -- given that pair<int,int> couple(1000*iz+ia,1000*jz+ja); is the parent-daugther uid pair --
     map<pair<int,int>, vector<double> > spall;   //the spallation cross section vector for all energies for each pair (primary, secondary)
     
     the spall cross section is obtained, depending on the xsec model (Galprop, Webber, Fluka) via 

     spall[couple] = xsecmodel[?]->GetXSec(iz,ia,jz,ja) (parent Z,A, daughter Z,A)
     or
     spall[couple][ip] = xsecmodel[?]->GetXSec(couple, energy[ip])
    
     GetXSec is a method of Galprop, Webber or Fluka classes -- that inherit from XSecBase class

     in the xsec class the xsections are stored in the private object  map<pair<int,int>, vector<double> > xsec;

     2) for nucleus --> antiprotons

     3) for nucleus --> leptons

     -- pair<int,int> coupleel(1001, -1000); // Electrons from protons
     -- pair<int,int> couplepos(1001, 1000); // Positrons from protons
     spall_apel[coupleel] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
     spall_apel[couplepos] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
            
     -- pair<int,int> coupleelHe(2004, -1000); // Electrons from He
     -- pair<int,int> coupleposHe(2004, 1000); // Positrons from He
     spall_apel[coupleelHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
     spall_apel[coupleposHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));

     In Galprop and Fluka modes, the vectors are filled in 

     TSpallationNetwork::InitXSecGalprop(factorelpos);
     TSpallationNetwork::InitXSecFluka(factorelpos);
     respectively

  */
    
  energy = co->GetEk();
  const int dimEn = energy.size();
  const double factor = Clight*1.e-27;
    
  //  vector<int> nuclei = TNucleiList::GetInstance()->GetList();
  //ofstream outf("spallation_xsec.dat", ios::in);
    
  for (int iloop = 0; iloop < nuclei.size()-1; ++iloop) {
        
    // If antiprotons or leptons, do not compute spallation. If protons, do not compute, because their spallation products will be computed elsewhere
    if (nuclei[iloop] <= 1001) continue;
        
    int iz = -1000;  // parent nucleus charge
    int ia = -1000;  // parent nucleus mass
    Utility::id_nuc(nuclei[iloop], ia, iz);
        
    for (int idaught = iloop+1; idaught < nuclei.size(); ++idaught) {
            
      int jz = -1000; // daughter nucleus charge
      int ja = -1000; // daughter nucleus mass
      Utility::id_nuc(nuclei[idaught], ja, jz);
            
      // If nucleus is antiproton of leptons skip it. Skip also if mass(daughter) > mass(parent)
      if (nuclei[idaught] < 1001 || ia < ja) continue;
            
      pair<int,int> couple(1000*iz+ia,1000*jz+ja);
            
      if (in->spallationxsec == GalpropXSec) spall[couple] = xsecmodel[0]->GetXSec(iz,ia,jz,ja);

      else if (in->spallationxsec == Webber03) {

	if (!xsecmodel[1]->IsPresent(couple)) spall[couple] = xsecmodel[0]->GetXSec(iz,ia,jz,ja);
	vector<double> beta = co->GetBeta();
	//	TWebber03::GetInstance()->Print(co);
                
	for(unsigned int ip = 0; ip < dimEn; ip++) {
	  spall[couple].push_back(xsecmodel[1]->GetXSec(couple, energy[ip])*factor*beta[ip]*(1.0*He_abundance*xsecmodel[1]->GetHefactor()));
                    
	}
      }
      else if (in->spallationxsec == Fluka) {

        vector<double> beta = co->GetBeta();

	if (!xsecmodel[1]->IsPresent(couple)) 
	  spall[couple] = xsecmodel[0]->GetXSec(iz,ia,jz,ja); //if the cross section is not present, use Galprop database
	else { //do nothing at the moment --> only gaprop files are actually read...

	  //if the cross section is present, then use Fluka model
         
	  /*
	    for(unsigned int ip = 0; ip < dimEn; ip++) {
	    
	    spall[couple].push_back(xsecmodel[1]->GetXSec(couple, energy[ip])*factor*beta[ip]*(1.0*He_abundance*xsecmodel[1]->GetHefactor()));
	    
            //fills spall[couple] with Fluka cross section. 
	    //Fluka->GetXSec(couple,energy) gives xsec in mbarn/GeV --> conversion to cm^3/s via "factor"

	   
 	    }*/
                    
	}
      }	
      else {
	cerr << "Wrong SpallationXSec option" << endl;
      }
            
      if(!in->SPALL) for(int i=0;i<spall[couple].size();i++) spall[couple][i]=0.; //fk130701 sets Xsecs to Zero i.e. no Spallation            //  for(unsigned int ip = 0; ip < dimEn; ip++) cout << couple.first << " " << couple.second << " " << spall[couple].back() << endl;

    }
  }
  // outf.close();
  //exit(-1);
  const double factorprot = factor*co->GetDeltaE(); 
  const double factorelpos = 1.e3*factorprot;  // barn/GeV --- *10^3 ---> mbarn/GeV --- *factorprot ---> cm^3/s
    
  double PP_inel = 0.0;
  double PA_inel = 0.0;
  double aPP_non = 0.0;
  double aPA_non = 0.0;
  double aPP_ann = 0.0;
  double aPA_ann = 0.0;
    
  if (in->spallationxsec != Fluka) {

    pair<int,int> coupleprpr(1001,1001);  // Secondary protons
    spall[coupleprpr] = vector<double>(dimEn, 0);

    for(unsigned int ip = 0; ip < dimEn; ++ip) {
      xsecmodel[0]->nucleon_cs(2, energy[ip], 1, 2, 4, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);  // Galprop CS
      spall[coupleprpr][ip] = in->SPALL*factorprot*(PP_inel + He_abundance*PA_inel);
    }
  } 
  else { //fluka

    pair<int,int> coupleppr(1001, 1001);  // Secondary protons, from protons
    pair<int,int> couplepHe(2004, 1001);  // Secondary protons, from Helium
            
    spall_apel[coupleppr] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
    spall_apel[couplepHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));

    TSpallationNetwork::InitXSecFluka(factorelpos,2); //2-->p
    //cout << "Secondary protons ok" << endl;

    /*vector<double> beta = co->GetBeta();
      for(unsigned int ip2 = 0; ip2 < dimEn; ip2++) {
      spall[coupleprpr].push_back(xsecmodel[1]->GetXSec(coupleprpr, energy[ip2])*factor*beta[ip2]*(1.0*He_abundance*xsecmodel[1]->GetHefactor()));

	
      }*/
	
  }
    
  // If antiprotons and/or leptons are wanted in output, add them, and add also secondary protons
    
  if (in->prop_ap || in->prop_lep || in->prop_deuteron) {
        
        
    if (in->prop_ap) {

      if (in->spallationxsec != Fluka) {

	pair<int,int> coupletert(-999,-999);  // Tertiary antiprotons
	spall[coupletert] = vector<double>(dimEn, 0.0);	
          
	for(unsigned int ip = 0; ip < dimEn; ++ip) {
	  xsecmodel[0]->nucleon_cs(2, energy[ip], -1, 2, 4, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);
                
	  spall[coupletert][ip] = factorprot*( aPP_non + He_abundance*aPA_non );
	  // Exploits approximation d\sigma/dEkin \propto 1/Ekin' [Tan & Ng '83]
	}
        
      }
      else { //fluka

	pair<int,int> coupleapap(-999, -999);  // Tertiary antiprotons, from sec. antiprotons
    
	spall_apel[coupleapap] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));

	TSpallationNetwork::InitXSecFluka(factorelpos,3); //3-->tertiary ap
	//cout << "Tertiary antiprotons ok" << endl;	

      }
      
      spall_apel[pair<int,int>(2003,-999)] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
      
      pair<int,int> coupleappr(1001,-999);  // Secondary antiprotons, from protons
      pair<int,int> coupleapHe(2004,-999);  // Secondary antiprotons, from Helium
            
      vector<double> pp = co->GetMomentum();
            
      spall_apel[coupleappr] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
      spall_apel[coupleapHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
      
      size_t limit = 100000;
      gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
            
      if(in->SPALL) //fk 130701
	{

	  if ( in->apy == FlukaAp) 
	    TSpallationNetwork::InitXSecFluka(factorelpos,1); //1-->ap
	  //cout << "Secondary antiprotons ok" << endl;	


	  for (unsigned int i = 0; i < dimEn; i++) {

	    if ( in->apy == GalpropFunction ){  
	      for (unsigned int ii=i+1; ii < dimEn; ii++) {
		spall_apel[coupleappr][i][ii] = factorelpos*energy[ii]*(xsecmodel[0]->antiproton_cc1(w,limit,in->antiproton_cs, pp[i], pp[ii], 1, 1, 1, 1) * ( (!in->scaling) + (in->scaling)*(0.12 * pow(energy[i], -1.67) + 1.78)) + (!in->scaling)*He_abundance*xsecmodel[0]->antiproton_cc1(w,limit,in->antiproton_cs, pp[i], pp[ii], 1, 1, 2, 4)); 
		spall_apel[coupleapHe][i][ii] = factorelpos*energy[ii]*4.0*(xsecmodel[0]->antiproton_cc1(w,limit,in->antiproton_cs, pp[i], 4.0*pp[ii], 2, 4, 1, 1) * ( (!in->scaling) + (in->scaling)*(0.12 * pow(energy[i], -1.67) + 1.78)) + (!in->scaling)*He_abundance*xsecmodel[0]->antiproton_cc1(w,limit,in->antiproton_cs, pp[i], 4.0*pp[ii], 2, 4, 2, 4)); 
	      }
	    }
	    else if ( in->apy == QGSJET ){
	     
	      // === TEST ===
	      spec_ini();
                        
	      double x,es,fff;
                        
	      for ( int i=1;i<=100;i++ ){
		x = (double)i/100.-.005;
		es = x*100.;
		fff=spec_int(100.,es,1,1);
		cout<<x<<"\t"<<fff<<endl;
	      }
	      exit(3);
	      // === END TEST ===

	    }
	    
	  }
	}
     
      gsl_integration_workspace_free(w);

    } // ap
        
        
    if (in->prop_lep) {
      pair<int,int> coupleel(1001, -1000); // Electrons from protons
      pair<int,int> couplepos(1001, 1000); // Positrons from protons
      spall_apel[coupleel] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
      spall_apel[couplepos] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
            
      pair<int,int> coupleelHe(2004, -1000); // Electrons from He
      pair<int,int> coupleposHe(2004, 1000); // Positrons from He
      spall_apel[coupleelHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
      spall_apel[coupleposHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
            
      spall_apel[pair<int,int>(2003,-1000)] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
      spall_apel[pair<int,int>(2003,1000)]  = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
            
      // Read in the tabulated database
            
            
      if(in->SPALL){ //fk 130628
	if (in->ly == GalpropTable) TSpallationNetwork::InitXSecGalprop(factorelpos);
	else if (in->ly == Pohl) TSpallationNetwork::InitXSecPohl(factorelpos);
	else if (in->ly == Kamae) TSpallationNetwork::InitXSecKamae(factorelpos);
	else if (in->ly == FlukaLep) {TSpallationNetwork::InitXSecFluka(factorelpos,0); cout << "Sec. leptons ok " << endl; } //0-->lep 
      }            

            
    } // el
        
    if (in->prop_deuteron) {
      pair<int,int> coupledeutdeut(-998,-998);  // Tertiary antideuterons
      pair<int,int> coupledeutapr(-999,-998);  // Secondary antideuterons, from antiprotons
      pair<int,int> coupledeutpr(1001,-998);  // Secondary antideuterons, from protons
      pair<int,int> coupledeutHe(2004,-998);  // Secondary antideuterons, from Helium
            
    }
        
  } // either
    
  return ;
}

vector<double> TSpallationNetwork::GetXSec(int i,int j) {
  pair<int,int> input(i,j);
  //  return spall[pair<int,int>(i,j)];
  map< pair<int,int>, vector<double> >::iterator it = spall.find(input);
  if (it != spall.end()) return (*it).second;
  else {
    //cerr << "No channel in TSpallationNetwork. Parent = " << i << " daughter = " << j << endl;
    return vector<double>();
  }
}

double TSpallationNetwork::GetXSec(int i, int j, double en) {
  if (en > energy.back()) {
    //cerr << "Out of range!" << endl;
    return 0;
  }
    
  int l = 0;
  for (l = 0; l < energy.size(); ++l) {
    if (en >= energy[l]) break;
  }
  //return spall[pair<int,int>(i,j)][l];
    
  pair<int,int> input(i,j);
  map< pair<int,int>, vector<double> >::iterator it = spall.find(input);
  if (it != spall.end()) {
    return (*it).second[l];
  }
  else {
    //    cerr << "Out of range in TSpallationNetwork" << endl;
    return 0;
  }
}

double TSpallationNetwork::GetXSec(int i, int j, int k) {
  if (k >= energy.size()) {
    //cerr << "Out of range!" << endl;
    return 0;
  }
  //  return spall[pair<int,int>(i,j)][k];
  pair<int,int> input(i,j);
  map< pair<int,int>, vector<double> >::iterator it = spall.find(input);
  if (it != spall.end()) return (*it).second[k];
  else {
    //cerr << "Out of range in TSpallationNetwork" << endl;
    return 0;
  }
}

vector<double> TSpallationNetwork::GetXSecApEl(int i, int j, int k) {
    
  if (k >= energy.size()) {
    //cerr << "Out of range!" << endl;
    return vector<double>();
  }
  //return spall_apel[pair<int,int>(i,j)][k];
    
  pair<int,int> input(i,j);
  map< pair<int,int>, vector< vector<double> > >::iterator it = spall_apel.find(input);
  if (it != spall_apel.end()) return (*it).second[k];
  else {
    //cerr << "No channel in TSpallationNetwork. Parent = " << i << " daughter = " << j << endl;
    return vector<double>();
  }
    
}

//TPP
vector<double> TSpallationNetwork::GetXSecTPP(vector<double> nu_vector) { // nu_vector -> vector of frequencies of ISRF
	
  //double photon_energy = 10.*1.e-9; //GeV
  double m_e           = 0.511e-3 ; //GeV
  double hPlanck = 4.135667e-15;
    
  vector<double> result(energy.size()*energy.size()*nu_vector.size());
    
  for (int k =0; k<energy.size(); ++k) {
    for (int inu=0; inu < nu_vector.size(); ++inu) {
      for (int l = 0; l < energy.size(); ++l) { // l -> energy index of parent electron
	result[(k*nu_vector.size() + inu)*energy.size() + l] = 0.;
      }
    }
  }
    
  for (int k =0; k<energy.size(); ++k) { // k -> energy index of secondary particle
        
    for (int inu=0; inu < nu_vector.size(); ++inu) {
            
      double photon_energy = hPlanck*nu_vector[inu]*1.e-9; //GeV
            
      if (k == 0)
	//cout << "***" << endl;
            
	for (int l = 0; l < energy.size(); ++l) { // l -> energy index of parent electron
                
	  int index = (k*nu_vector.size() + inu)*energy.size() + l;
                
	  double Average_TPP_positron_energy = 0.5 * sqrt(energy[l] / photon_energy ) * m_e;
                
	  double s = energy[l] * hPlanck*nu_vector[inu]*1.e-9 / pow(m_e, 2.); //tutto in GeV; s: numero puro
                
	  double correction = 0.;
	  if (s>4.) { 
	    correction = 0.1193662073189215 * exp(-3.839009766244388 - 0.1223186374016053*pow(log(-4. + s), 2) ) * pow((-4. + s),1.8194204144025334);
                    
	    if (s>79.)
	      correction = 0.13 * (-8.07407 + 3.11111 * log(0.0153186 * energy[l] * 3.));
                    
	  }
                
	  if (k == 0)
	    cout << " photon energy [eV] " << photon_energy*1.e9 << " energy of electron [GeV] -> " << energy[l] << " s-> " << s << " correction-> " << correction << endl;
                
                
	  double sigma_tot =  9.46 * 6.65 * 1.e-2 * (1/137.) * correction;// * temp; // c sigma_Thompson alpha_F; units: cm/Myr * cm^2
	  //double std_dev = 30;//Average_TPP_positron_energy;
	  double std_dev = sqrt(Average_TPP_positron_energy); // GeV
	  //double x = energy[l]*photon_energy/(m_e*m_e);
	  //if (x>4)
	  result[index] = sigma_tot * energy[l] * exp( - pow(energy[k] - Average_TPP_positron_energy,2.0) / (2.*(pow(std_dev,2.0))) ) / (sqrt(2*3.14) * std_dev);
	  //			   cm^3/Myr   GeV          1/GeV
	}
    }
  }
    
  return result;
    
}

//*********************************************************************************************************************************************

void TSpallationNetwork::InitDataTablesGalprop() {
    
  cout << "Please provide Galprop leptonic cross sections: data/Electron_production.dat and data/Positron_production.dat" << endl;

  Matrix_El_pp = vector<double>(401*801, 0.0);
  Matrix_El_pHe = vector<double>(401*801, 0.0);
  Matrix_El_Hep = vector<double>(401*801, 0.0);
  Matrix_El_HeHe = vector<double>(401*801, 0.0);
    
  Matrix_Pos_pp = vector<double>(401*801, 0.0);
  Matrix_Pos_pHe = vector<double>(401*801, 0.0);
  Matrix_Pos_Hep = vector<double>(401*801, 0.0);
  Matrix_Pos_HeHe = vector<double>(401*801, 0.0);
    
  ifstream infile(ElTablefile.c_str(), ios::in);
  double a,b,c,d;
  Nelectrons = 401;
  Nprotons = 801;

  for (int i = 0; i < 401; i++) {
    for (int j = 0; j < 801; j++) {
            
      infile >> a >> b >> c >> d;
            
      int ind = index_matrix(i,j);
            
      Matrix_El_pp[ind] = a;
      Matrix_El_pHe[ind] = b;
      Matrix_El_Hep[ind] = c;
      Matrix_El_HeHe[ind] = d;
    }
  }
    
  infile.close();
    
  infile.open(PosTablefile.c_str(), ios::in);
    
  for (int i = 0; i < 401; i++) {
    for (int j = 0; j < 801; j++) {
            
      infile >> a >> b >> c >> d;
            
      int ind = index_matrix(i,j);
            
      Matrix_Pos_pp[ind] = a;
      Matrix_Pos_pHe[ind] = b;
      Matrix_Pos_Hep[ind] = c;
      Matrix_Pos_HeHe[ind] = d;
    }
  }
    
  infile.close();
    
  return;
}
//ofstream outf("total_xsec.dat", ios::out);

void TSpallationNetwork::InitXSecKamae(double factorelpos) { 
    
  pair<int,int> coupleel(1001, -1000); // Electrons from protons
  pair<int,int> couplepos(1001, 1000); // Positrons from protons
  pair<int,int> coupleelHe(2004, -1000); // Electrons from He
  pair<int,int> coupleposHe(2004, 1000); // Positrons from He
    
  const int dimEn = energy.size();
    
  for (unsigned int j=0; j<dimEn; j++){
        
    double Epr = energy[j];
        
    for (unsigned int i = 0; i < dimEn; i++) {
            
      double Eel = min(1e5,energy[i]);
            
      double cs_pp = KamaeYields::GetSigma(Eel, Epr, ID_ELECTRON);
            
      double cs_pHe = pow(4.0,2.0/3.0)*KamaeYields::GetSigma(Eel, Epr, ID_ELECTRON);
            
      double cs_Hep = pow(4.0,2.0/3.0)*KamaeYields::GetSigma(Eel, Epr, ID_ELECTRON);
            
      double cs_HeHe = pow(4.0*4.0,2.0/3.0)*KamaeYields::GetSigma(Eel, Epr, ID_ELECTRON);  
            
      spall_apel[coupleel][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;
      spall_apel[coupleelHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;	
            
      cs_pp = KamaeYields::GetSigma(Eel, Epr, ID_POSITRON);
            
      cs_pHe = pow(4.0,2.0/3.0)*KamaeYields::GetSigma(Eel, Epr, ID_POSITRON);
            
      cs_Hep = pow(4.0,2.0/3.0)*KamaeYields::GetSigma(Eel, Epr, ID_POSITRON);
            
      cs_HeHe = pow(4.0*4.0,2.0/3.0)*KamaeYields::GetSigma(Eel, Epr, ID_POSITRON);
            
      spall_apel[couplepos][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;   // H
      spall_apel[coupleposHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;   // He            
            
    }
        
  }
    
  return ;
}

void TSpallationNetwork::InitXSecGalprop(double factorelpos) { 
    
  pair<int,int> coupleel(1001, -1000); // Electrons from protons
  pair<int,int> couplepos(1001, 1000); // Positrons from protons
  pair<int,int> coupleelHe(2004, -1000); // Electrons from He
  pair<int,int> coupleposHe(2004, 1000); // Positrons from He
    
  const int dimEn = energy.size();
    
  const int dimlept = 400;
  const double DBlog = (log10(1e5)-log10(1e-3))/(double)dimlept;
  double Elept[dimlept+1];
  for (int i = 0; i <= dimlept; i++) Elept[i] = pow(10, log10(1e-3)+double(i)*DBlog);
    
  const int dimpr = 800;
  const double DBprlog = (log10(1e5)-log10(1e-3))/double(dimpr);
  double ET[dimpr+1];
  for (int j = 0; j <= dimpr; j++) ET[j] = pow(10, log10(1e-3) + double(j)*DBprlog);
    
  TSpallationNetwork::InitDataTablesGalprop();
    
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline_El_pp = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_pHe = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_Hep = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_HeHe = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pp = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pHe = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_Hep = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_HeHe = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    
  //gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline_El_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_pHe_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_HeHe_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pHe_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_HeHe_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    
  double Vectorlept_El_pp[dimlept+1];
  double Vectorlept_El_pHe[dimlept+1];
  double Vectorlept_El_Hep[dimlept+1];
  double Vectorlept_El_HeHe[dimlept+1];
  double Vectorlept_Pos_pp[dimlept+1];
  double Vectorlept_Pos_pHe[dimlept+1];
  double Vectorlept_Pos_Hep[dimlept+1];
  double Vectorlept_Pos_HeHe[dimlept+1];
  double Vectorlept_El_pp_up[dimlept+1];
  double Vectorlept_El_pHe_up[dimlept+1];
  double Vectorlept_El_Hep_up[dimlept+1];
  double Vectorlept_El_HeHe_up[dimlept+1];
  double Vectorlept_Pos_pp_up[dimlept+1];
  double Vectorlept_Pos_pHe_up[dimlept+1];
  double Vectorlept_Pos_Hep_up[dimlept+1];
  double Vectorlept_Pos_HeHe_up[dimlept+1];
    
  for (unsigned int j=0; j<dimEn; j++){
        
    double Epr = energy[j];
    double momentum = sqrt(Epr*Epr + 2.0*mp*Epr);
    int j_pr = int(floor(log10(Epr/ET[0])/DBprlog));
    if (j_pr > dimpr-1) j_pr = dimpr-1;
    double u = (Epr-ET[j_pr])/(ET[j_pr+1]-ET[j_pr]);
        
    for (int i = 0; i <= dimlept; i++) {
      int index = index_matrix(i,j_pr);
      Vectorlept_El_pp[i] = Matrix_El_pp[index];
      Vectorlept_El_pHe[i] = Matrix_El_pHe[index];
      Vectorlept_El_Hep[i] = Matrix_El_Hep[index];
      Vectorlept_El_HeHe[i] = Matrix_El_HeHe[index];
      Vectorlept_Pos_pp[i] = Matrix_Pos_pp[index];
      Vectorlept_Pos_pHe[i] = Matrix_Pos_pHe[index];
      Vectorlept_Pos_Hep[i] = Matrix_Pos_Hep[index];
      Vectorlept_Pos_HeHe[i] = Matrix_Pos_HeHe[index];    
      //
      if (j_pr <= dimpr)
	index = index_matrix(i,j_pr+1);
      else
	index = dimpr;
      Vectorlept_El_pp_up[i] = Matrix_El_pp[index];
      Vectorlept_El_pHe_up[i] = Matrix_El_pHe[index];
      Vectorlept_El_Hep_up[i] = Matrix_El_Hep[index];
      Vectorlept_El_HeHe_up[i] = Matrix_El_HeHe[index];
      Vectorlept_Pos_pp_up[i] = Matrix_Pos_pp[index];
      Vectorlept_Pos_pHe_up[i] = Matrix_Pos_pHe[index];
      Vectorlept_Pos_Hep_up[i] = Matrix_Pos_Hep[index];
      Vectorlept_Pos_HeHe_up[i] = Matrix_Pos_HeHe[index];    
    }
        
    gsl_spline_init(spline_El_pp, Elept, Vectorlept_El_pp, dimlept+1);
    gsl_spline_init(spline_El_pHe, Elept, Vectorlept_El_pHe, dimlept+1);
    gsl_spline_init(spline_El_Hep, Elept, Vectorlept_El_Hep, dimlept+1);
    gsl_spline_init(spline_El_HeHe, Elept, Vectorlept_El_HeHe, dimlept+1);
    gsl_spline_init(spline_Pos_pp, Elept, Vectorlept_Pos_pp, dimlept+1);
    gsl_spline_init(spline_Pos_pHe, Elept, Vectorlept_Pos_pHe, dimlept+1);
    gsl_spline_init(spline_Pos_Hep, Elept, Vectorlept_Pos_Hep, dimlept+1);
    gsl_spline_init(spline_Pos_HeHe, Elept, Vectorlept_Pos_HeHe, dimlept+1);
    //
    gsl_spline_init(spline_El_pp_up, Elept, Vectorlept_El_pp_up, dimlept+1);
    gsl_spline_init(spline_El_pHe_up, Elept, Vectorlept_El_pHe_up, dimlept+1);
    gsl_spline_init(spline_El_Hep_up, Elept, Vectorlept_El_Hep_up, dimlept+1);
    gsl_spline_init(spline_El_HeHe_up, Elept, Vectorlept_El_HeHe_up, dimlept+1);
    gsl_spline_init(spline_Pos_pp_up, Elept, Vectorlept_Pos_pp_up, dimlept+1);
    gsl_spline_init(spline_Pos_pHe_up, Elept, Vectorlept_Pos_pHe_up, dimlept+1);
    gsl_spline_init(spline_Pos_Hep_up, Elept, Vectorlept_Pos_Hep_up, dimlept+1);
    gsl_spline_init(spline_Pos_HeHe_up, Elept, Vectorlept_Pos_HeHe_up, dimlept+1);
        
    for (unsigned int i = 0; i < dimEn; i++) {
            
      double Eel = min(1e5,energy[i]);
      double Eel_tot = energy[i] + MeleGeV;
            
      int i_lept = int(floor(log10(Eel/Elept[0])/DBlog));
      if (i_lept > dimlept-1) i_lept = dimlept-1;
      //double t = (Eel-Elept[i_lept])/(Elept[i_lept+1]-Elept[i_lept]);
            
      double valuefix = gsl_spline_eval(spline_El_pp, Eel, acc);
      double valueup = gsl_spline_eval(spline_El_pp_up, Eel, acc);
      double cs_pp = valuefix*(1-u) + valueup*u;
            
      valuefix = gsl_spline_eval(spline_El_pHe, Eel, acc);
      valueup = gsl_spline_eval(spline_El_pHe_up, Eel, acc);
      double cs_pHe = valuefix*(1-u) + valueup*u;
            
      valuefix = gsl_spline_eval(spline_El_Hep, Eel, acc);
      valueup = gsl_spline_eval(spline_El_Hep_up, Eel, acc);
      double cs_Hep = valuefix*(1-u) + valueup*u;
            
      valuefix = gsl_spline_eval(spline_El_HeHe, Eel, acc);
      valueup = gsl_spline_eval(spline_El_HeHe_up, Eel, acc);
      double cs_HeHe = valuefix*(1-u) + valueup*u;  
            
      spall_apel[coupleel][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;
      spall_apel[coupleelHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;	
            
      //
            
      valuefix = gsl_spline_eval(spline_Pos_pp, Eel, acc);
      valueup = gsl_spline_eval(spline_Pos_pp_up, Eel, acc);
      cs_pp = valuefix*(1-u) + valueup*u;
            
      valuefix = gsl_spline_eval(spline_Pos_pHe, Eel, acc);
      valueup = gsl_spline_eval(spline_Pos_pHe_up, Eel, acc);
      cs_pHe = valuefix*(1-u) + valueup*u;
            
      valuefix = gsl_spline_eval(spline_Pos_Hep, Eel, acc);
      valueup = gsl_spline_eval(spline_Pos_Hep_up, Eel, acc);
      cs_Hep = valuefix*(1-u) + valueup*u;
            
      valuefix = gsl_spline_eval(spline_Pos_HeHe, Eel, acc);
      valueup = gsl_spline_eval(spline_Pos_HeHe_up, Eel, acc);
      cs_HeHe = valuefix*(1-u) + valueup*u; 
            
      spall_apel[couplepos][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;   // H
      spall_apel[coupleposHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;   // He
            
    }
        
  }
    
    
  return ;
}

//*********************************************************************************************************************************************
// Fluka model
//*********************************************************************************************************************************************

void TSpallationNetwork::InitDataTablesFluka(int which_particle) {

  cout << "Entering in TSpallationNetwork InitDataTablesFluka routine ..." <<  endl;
  cout << "which particle? "<< which_particle << endl;

  //reads from the tables

  if (which_particle == 0) { //Leptonic tables

    cout << "leptonic tables" << endl;

    Matrix_El_pp = vector<double>(408*380, 0.0);
    Matrix_El_pHe = vector<double>(408*380, 0.0);
    Matrix_El_Hep = vector<double>(408*380, 0.0);
    Matrix_El_HeHe = vector<double>(408*380, 0.0);
    
    Matrix_Pos_pp = vector<double>(408*380, 0.0);
    Matrix_Pos_pHe = vector<double>(408*380, 0.0);
    Matrix_Pos_Hep = vector<double>(408*380, 0.0);
    Matrix_Pos_HeHe = vector<double>(408*380, 0.0);
    
    ifstream infile(FlukaElTablefile.c_str(), ios::in);
    double a,b,c,d;
    Nelectrons = 408;
    Nprotons = 380;

    for (int i = 0; i < Nelectrons; i++) {
      for (int j = 0; j < Nprotons; j++) {
            
        infile >> a >> b >> c >> d;
            
        int ind = index_matrix(i,j);
            
        Matrix_El_pp[ind] = a;
        Matrix_El_pHe[ind] = b;
        Matrix_El_Hep[ind] = c;
        Matrix_El_HeHe[ind] = d;
      }
    }
    
    infile.close();
    
    infile.open(FlukaPosTablefile.c_str(), ios::in);
    
    for (int i = 0; i < Nelectrons; i++) {
      for (int j = 0; j < Nprotons; j++) {
            
        infile >> a >> b >> c >> d;
            
        int ind = index_matrix(i,j);
            
        Matrix_Pos_pp[ind] = a;
        Matrix_Pos_pHe[ind] = b;
        Matrix_Pos_Hep[ind] = c;
        Matrix_Pos_HeHe[ind] = d;
      }
    }
    
    infile.close();
    
  }

  if (which_particle == 1) { //sec Antiproton tables

    cout << "ap tables " << endl;

    Matrix_Ap_pp   = vector<double>(408*380, 0.0);
    Matrix_Ap_pHe  = vector<double>(408*380, 0.0);
    Matrix_Ap_Hep  = vector<double>(408*380, 0.0);
    Matrix_Ap_HeHe = vector<double>(408*380, 0.0);
    ifstream infile(FlukaApTablefile.c_str(), ios::in);
    double a,b,c,d;
    Nap = 408;
    Nprotons = 380;
    for (int i = 0; i < Nap; i++) {
      for (int j = 0; j < Nprotons; j++) {
        infile >> a >> b >> c >> d;
        int ind = index_matrix(i,j);
        Matrix_Ap_pp[ind] = a;
        Matrix_Ap_pHe[ind] = b;
        Matrix_Ap_Hep[ind] = c;
        Matrix_Ap_HeHe[ind] = d;
      }
    }
    infile.close();
  }

  if (which_particle == 2) { //sec Proton tables

    cout << "proton tables " << endl;

    Matrix_p_pp   = vector<double>(408*380, 0.0);
    Matrix_p_pHe  = vector<double>(408*380, 0.0);
    Matrix_p_Hep  = vector<double>(408*380, 0.0);
    Matrix_p_HeHe = vector<double>(408*380, 0.0);
    ifstream infile(FlukaProtTablefile.c_str(), ios::in);
    double a,b,c,d;
    Nap = 408;
    Nprotons = 380;
    for (int i = 0; i < Nap; i++) {
      for (int j = 0; j < Nprotons; j++) {
        infile >> a >> b >> c >> d;
        int ind = index_matrix(i,j);
        Matrix_p_pp[ind] = a;
        Matrix_p_pHe[ind] = b;
        Matrix_p_Hep[ind] = c;
        Matrix_p_HeHe[ind] = d;
        if (i%100==0 && j%100==0) 
	  cout << i << " <-i,j-> " << j << ";  abcd= " << a << " " << b << " " << c << "  " << d << endl;
      }
    }
    infile.close();
  }

  if (which_particle == 2) { //tertiary antiProton tables

    cout << "tertiary antiproton tables " << endl;

    Matrix_3Ap_app   = vector<double>(408*380, 0.0);
    Matrix_3Ap_apHe  = vector<double>(408*380, 0.0);
    ifstream infile(FlukaTertiaryApTablefile.c_str(), ios::in);
    double a,b,c,d;
    Nap = 408;
    Nprotons = 380;
    for (int i = 0; i < Nap; i++) {
      for (int j = 0; j < Nprotons; j++) {
        infile >> a >> b;
        int ind = index_matrix(i,j);
        Matrix_3Ap_app[ind] = a;
        Matrix_3Ap_apHe[ind] = b;
        if (i%100==0 && j%100==0) 
	  cout << i << " <-i,j-> " << j << ";  abcd= " << a << " " << b << endl;
      }
    }
    infile.close();
  }

  return;

}

void TSpallationNetwork::InitXSecFluka(double factorelpos, int which_particle) { 

  //which_particle: 0-->leptons
  //1-->sec antiprotons
  //2-->sec protons
  //3-->tertiary antiprotons

  cout << "Entering in TSpallationNetwork InitXSecFluka routine ..." <<  endl;

  pair<int,int> coupleel(1001, -1000); // Electrons from protons
  pair<int,int> couplepos(1001, 1000); // Positrons from protons
  pair<int,int> coupleelHe(2004, -1000); // Electrons from He
  pair<int,int> coupleposHe(2004, 1000); // Positrons from He

  pair<int,int> coupleappr(1001,-999);  // Secondary antiprotons, from CR protons
  pair<int,int> coupleapHe(2004,-999);  // Secondary antiprotons, from CR Helium
    
  pair<int,int> coupleppr(1001,1001);  // Secondary antiprotons, from CR protons
  pair<int,int> couplepHe(2004,1001);  // Secondary antiprotons, from CR Helium

  pair<int,int> coupleapap(-999,-999);  // Tertiary antiprotons

  const int dimEn = energy.size();
    
  const int dimlept = 407;
  const int dimap   = 407;
  const int dimp   = 407;
  const double DBlog = log10(1.05);
  double Elept[dimlept+1];
  double Eap[dimap+1];
  for (int i = 0; i <= dimlept; i++) {
    Elept[i] = 0.001*pow(1.05, double(i)+0.5);
    Elept[i] = pow(10., log10(1e-3)+double(i)*DBlog+0.5*DBlog);
    Eap[i] = Elept[i];	
    //cout << i << " " << Elept[i] << endl;
  }
    
  const int dimpr = 379;
  const double DBprlog = log10(1.05);
  double ET[dimpr+1];
  for (int j = 0; j <= 94; j++) ET[j] = pow(10, log10(1e-3) + double(j)*DBprlog);
  for (int j = 0; j <= 284; j++) {
    ET[j+95] = 0.1*pow(1.05, double(j));
    ET[j+95] = pow(10, log10(1e-1) + double(j)*DBprlog);
  }
  //   for (int j = 0; j <= dimpr; j++) {
  //     cout << j << " " << ET[j] << endl;
  //   }
    
  TSpallationNetwork::InitDataTablesFluka(which_particle);

  cout << "test" << endl;
    
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline_El_pp = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_pHe = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_Hep = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_HeHe = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pp = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pHe = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_Hep = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_HeHe = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  //
  gsl_spline *spline_Ap_pp = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
  gsl_spline *spline_Ap_pHe = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
  gsl_spline *spline_Ap_Hep = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
  gsl_spline *spline_Ap_HeHe = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
  //
  gsl_spline *spline_p_pp = gsl_spline_alloc(gsl_interp_cspline, dimp+1);
  gsl_spline *spline_p_pHe = gsl_spline_alloc(gsl_interp_cspline, dimp+1);
  gsl_spline *spline_p_Hep = gsl_spline_alloc(gsl_interp_cspline, dimp+1);
  gsl_spline *spline_p_HeHe = gsl_spline_alloc(gsl_interp_cspline, dimp+1);
  //
  gsl_spline *spline_3Ap_app  = gsl_spline_alloc(gsl_interp_cspline, dimp+1);
  gsl_spline *spline_3Ap_apHe = gsl_spline_alloc(gsl_interp_cspline, dimp+1);   
 
  //gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline_El_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_pHe_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_HeHe_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pHe_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_HeHe_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  //
  gsl_spline *spline_Ap_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
  gsl_spline *spline_Ap_pHe_up = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
  gsl_spline *spline_Ap_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
  gsl_spline *spline_Ap_HeHe_up = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
  //
  gsl_spline *spline_p_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimp+1);
  gsl_spline *spline_p_pHe_up = gsl_spline_alloc(gsl_interp_cspline, dimp+1);
  gsl_spline *spline_p_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimp+1);
  gsl_spline *spline_p_HeHe_up = gsl_spline_alloc(gsl_interp_cspline, dimp+1);
  //
  gsl_spline *spline_3Ap_app_up = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
  gsl_spline *spline_3Ap_apHe_up = gsl_spline_alloc(gsl_interp_cspline, dimap+1);
    
  double Vectorlept_El_pp[dimlept+1];
  double Vectorlept_El_pHe[dimlept+1];
  double Vectorlept_El_Hep[dimlept+1];
  double Vectorlept_El_HeHe[dimlept+1];
  double Vectorlept_Pos_pp[dimlept+1];
  double Vectorlept_Pos_pHe[dimlept+1];
  double Vectorlept_Pos_Hep[dimlept+1];
  double Vectorlept_Pos_HeHe[dimlept+1];
  //
  double Vector_Ap_pp[dimap+1];
  double Vector_Ap_pHe[dimap+1];
  double Vector_Ap_Hep[dimap+1];
  double Vector_Ap_HeHe[dimap+1];
  //
  double Vector_p_pp[dimp+1];
  double Vector_p_pHe[dimp+1];
  double Vector_p_Hep[dimp+1];
  double Vector_p_HeHe[dimp+1];
  //
  double Vector_3Ap_app[dimap+1];
  double Vector_3Ap_apHe[dimap+1];
  //
  double Vectorlept_El_pp_up[dimlept+1];
  double Vectorlept_El_pHe_up[dimlept+1];
  double Vectorlept_El_Hep_up[dimlept+1];
  double Vectorlept_El_HeHe_up[dimlept+1];
  double Vectorlept_Pos_pp_up[dimlept+1];
  double Vectorlept_Pos_pHe_up[dimlept+1];
  double Vectorlept_Pos_Hep_up[dimlept+1];
  double Vectorlept_Pos_HeHe_up[dimlept+1];
  //
  double Vector_Ap_pp_up[dimap+1];
  double Vector_Ap_pHe_up[dimap+1];
  double Vector_Ap_Hep_up[dimap+1];
  double Vector_Ap_HeHe_up[dimap+1];
  //
  double Vector_p_pp_up[dimp+1];
  double Vector_p_pHe_up[dimp+1];
  double Vector_p_Hep_up[dimp+1];
  double Vector_p_HeHe_up[dimp+1];
  //
  double Vector_3Ap_app_up[dimp+1];
  double Vector_3Ap_apHe_up[dimp+1];

  if (which_particle == 0) { //leptons
    
    for (unsigned int j=0; j<dimEn; j++){
        
      double Epr = energy[j];
      double momentum = sqrt(Epr*Epr + 2.0*mp*Epr);
      int j_pr = int(floor(log10(Epr/ET[0])/DBprlog));
      if (j_pr > dimpr-1) j_pr = dimpr-1;
      double u = (Epr-ET[j_pr])/(ET[j_pr+1]-ET[j_pr]);
        
      for (int i = 0; i <= dimlept; i++) {
	int index = index_matrix(i,j_pr);
	Vectorlept_El_pp[i] = Matrix_El_pp[index];
	Vectorlept_El_pHe[i] = Matrix_El_pHe[index];
	Vectorlept_El_Hep[i] = Matrix_El_Hep[index];
	Vectorlept_El_HeHe[i] = Matrix_El_HeHe[index];
	Vectorlept_Pos_pp[i] = Matrix_Pos_pp[index];
	Vectorlept_Pos_pHe[i] = Matrix_Pos_pHe[index];
	Vectorlept_Pos_Hep[i] = Matrix_Pos_Hep[index];
	Vectorlept_Pos_HeHe[i] = Matrix_Pos_HeHe[index];    
	//
	if (j_pr <= dimpr)
	  index = index_matrix(i,j_pr+1);
	else
	  index = dimpr;
	Vectorlept_El_pp_up[i] = Matrix_El_pp[index];
	Vectorlept_El_pHe_up[i] = Matrix_El_pHe[index];
	Vectorlept_El_Hep_up[i] = Matrix_El_Hep[index];
	Vectorlept_El_HeHe_up[i] = Matrix_El_HeHe[index];
	Vectorlept_Pos_pp_up[i] = Matrix_Pos_pp[index];
	Vectorlept_Pos_pHe_up[i] = Matrix_Pos_pHe[index];
	Vectorlept_Pos_Hep_up[i] = Matrix_Pos_Hep[index];
	Vectorlept_Pos_HeHe_up[i] = Matrix_Pos_HeHe[index];    
      }
        
      gsl_spline_init(spline_El_pp, Elept, Vectorlept_El_pp, dimlept+1);
      gsl_spline_init(spline_El_pHe, Elept, Vectorlept_El_pHe, dimlept+1);
      gsl_spline_init(spline_El_Hep, Elept, Vectorlept_El_Hep, dimlept+1);
      gsl_spline_init(spline_El_HeHe, Elept, Vectorlept_El_HeHe, dimlept+1);
      gsl_spline_init(spline_Pos_pp, Elept, Vectorlept_Pos_pp, dimlept+1);
      gsl_spline_init(spline_Pos_pHe, Elept, Vectorlept_Pos_pHe, dimlept+1);
      gsl_spline_init(spline_Pos_Hep, Elept, Vectorlept_Pos_Hep, dimlept+1);
      gsl_spline_init(spline_Pos_HeHe, Elept, Vectorlept_Pos_HeHe, dimlept+1);
      //
      gsl_spline_init(spline_El_pp_up, Elept, Vectorlept_El_pp_up, dimlept+1);
      gsl_spline_init(spline_El_pHe_up, Elept, Vectorlept_El_pHe_up, dimlept+1);
      gsl_spline_init(spline_El_Hep_up, Elept, Vectorlept_El_Hep_up, dimlept+1);
      gsl_spline_init(spline_El_HeHe_up, Elept, Vectorlept_El_HeHe_up, dimlept+1);
      gsl_spline_init(spline_Pos_pp_up, Elept, Vectorlept_Pos_pp_up, dimlept+1);
      gsl_spline_init(spline_Pos_pHe_up, Elept, Vectorlept_Pos_pHe_up, dimlept+1);
      gsl_spline_init(spline_Pos_Hep_up, Elept, Vectorlept_Pos_Hep_up, dimlept+1);
      gsl_spline_init(spline_Pos_HeHe_up, Elept, Vectorlept_Pos_HeHe_up, dimlept+1);
        
      for (unsigned int i = 0; i < dimEn; i++) {
            
	double Eel = min(1e5,energy[i]);
	double Eel_tot = energy[i] + MeleGeV;
            
	int i_lept = int(floor(log10(Eel/Elept[0])/DBlog));
	if (i_lept > dimlept-1) i_lept = dimlept-1;
	//double t = (Eel-Elept[i_lept])/(Elept[i_lept+1]-Elept[i_lept]);
            
	double valuefix = gsl_spline_eval(spline_El_pp, Eel, acc);
	double valueup = gsl_spline_eval(spline_El_pp_up, Eel, acc);
	double cs_pp = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_El_pHe, Eel, acc);
	valueup = gsl_spline_eval(spline_El_pHe_up, Eel, acc);
	double cs_pHe = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_El_Hep, Eel, acc);
	valueup = gsl_spline_eval(spline_El_Hep_up, Eel, acc);
	double cs_Hep = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_El_HeHe, Eel, acc);
	valueup = gsl_spline_eval(spline_El_HeHe_up, Eel, acc);
	double cs_HeHe = valuefix*(1-u) + valueup*u;  
            
	spall_apel[coupleel][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;
	spall_apel[coupleelHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;	
            
	//
            
	valuefix = gsl_spline_eval(spline_Pos_pp, Eel, acc);
	valueup = gsl_spline_eval(spline_Pos_pp_up, Eel, acc);
	cs_pp = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_Pos_pHe, Eel, acc);
	valueup = gsl_spline_eval(spline_Pos_pHe_up, Eel, acc);
	cs_pHe = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_Pos_Hep, Eel, acc);
	valueup = gsl_spline_eval(spline_Pos_Hep_up, Eel, acc);
	cs_Hep = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_Pos_HeHe, Eel, acc);
	valueup = gsl_spline_eval(spline_Pos_HeHe_up, Eel, acc);
	cs_HeHe = valuefix*(1-u) + valueup*u; 
            
	spall_apel[couplepos][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;   // H
	spall_apel[coupleposHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;   // He
            
            
      }
        
    }

  } else if (which_particle == 1) { //sec antiprotons

    //
    for (unsigned int j=0; j<dimEn; j++){
        
      double Epr = energy[j];
      double momentum = sqrt(Epr*Epr + 2.0*mp*Epr);
      int j_pr = int(floor(log10(Epr/ET[0])/DBprlog));
      if (j_pr > dimpr-1) j_pr = dimpr-1;
      double u = (Epr-ET[j_pr])/(ET[j_pr+1]-ET[j_pr]);
        
      for (int i = 0; i <= dimlept; i++) {
	int index = index_matrix(i,j_pr);
	Vector_Ap_pp[i] = Matrix_Ap_pp[index];
	Vector_Ap_pHe[i] = Matrix_Ap_pHe[index];
	Vector_Ap_Hep[i] = Matrix_Ap_Hep[index];
	Vector_Ap_HeHe[i] = Matrix_Ap_HeHe[index];
	//
	if (j_pr <= dimpr)
	  index = index_matrix(i,j_pr+1);
	else
	  index = dimpr;
	Vector_Ap_pp_up[i] = Matrix_Ap_pp[index];
	Vector_Ap_pHe_up[i] = Matrix_Ap_pHe[index];
	Vector_Ap_Hep_up[i] = Matrix_Ap_Hep[index];
	Vector_Ap_HeHe_up[i] = Matrix_Ap_HeHe[index];
      }
        
      gsl_spline_init(spline_Ap_pp, Eap, Vector_Ap_pp, dimap+1);
      gsl_spline_init(spline_Ap_pHe, Eap, Vector_Ap_pHe, dimap+1);
      gsl_spline_init(spline_Ap_Hep, Eap, Vector_Ap_Hep, dimap+1);
      gsl_spline_init(spline_Ap_HeHe, Eap, Vector_Ap_HeHe, dimap+1);
      //
      gsl_spline_init(spline_Ap_pp_up, Eap, Vector_Ap_pp_up, dimap+1);
      gsl_spline_init(spline_Ap_pHe_up, Eap, Vector_Ap_pHe_up, dimap+1);
      gsl_spline_init(spline_Ap_Hep_up, Eap, Vector_Ap_Hep_up, dimap+1);
      gsl_spline_init(spline_Ap_HeHe_up, Eap, Vector_Ap_HeHe_up, dimap+1);
        
      for (unsigned int i = 0; i < dimEn; i++) {
            
	double Eap_i = min(1e5,energy[i]);
	double Eap_i_tot = energy[i] + MpGeV;
            
	int i_ap = int(floor(log10(Eap_i/Eap[0])/DBlog));
	if (i_ap > dimap-1) i_ap = dimap-1;
	//double t = (Eel-Elept[i_lept])/(Elept[i_lept+1]-Elept[i_lept]);
            
	double valuefix = gsl_spline_eval(spline_Ap_pp, Eap_i, acc);
	double valueup = gsl_spline_eval(spline_Ap_pp_up, Eap_i, acc);
	double cs_pp = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_Ap_pHe, Eap_i, acc);
	valueup = gsl_spline_eval(spline_Ap_pHe_up, Eap_i, acc);
	double cs_pHe = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_Ap_Hep, Eap_i, acc);
	valueup = gsl_spline_eval(spline_Ap_Hep_up, Eap_i, acc);
	double cs_Hep = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_Ap_HeHe, Eap_i, acc);
	valueup = gsl_spline_eval(spline_Ap_HeHe_up, Eap_i, acc);
	double cs_HeHe = valuefix*(1-u) + valueup*u;  
            
	spall_apel[coupleappr][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;
	spall_apel[coupleapHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;	
            
      }
        
    }


  } else if (which_particle == 2) {  //sec protons

    //
    cout << "test" << endl; 
    for (unsigned int j=0; j<dimEn; j++){

      //cout << j << endl;
        
      double Epr = energy[j];
      double momentum = sqrt(Epr*Epr + 2.0*mp*Epr);
      int j_pr = int(floor(log10(Epr/ET[0])/DBprlog));
      if (j_pr > dimpr-1) j_pr = dimpr-1;
      double u = (Epr-ET[j_pr])/(ET[j_pr+1]-ET[j_pr]);
        
      for (int i = 0; i <= dimlept; i++) {
	int index = index_matrix(i,j_pr);
	Vector_p_pp[i] = Matrix_p_pp[index];
	Vector_p_pHe[i] = Matrix_p_pHe[index];
	Vector_p_Hep[i] = Matrix_p_Hep[index];
	Vector_p_HeHe[i] = Matrix_p_HeHe[index];
	//
	if (j_pr <= dimpr)
	  index = index_matrix(i,j_pr+1);
	else
	  index = dimpr;
	Vector_p_pp_up[i] = Matrix_p_pp[index];
	Vector_p_pHe_up[i] = Matrix_p_pHe[index];
	Vector_p_Hep_up[i] = Matrix_p_Hep[index];
	Vector_p_HeHe_up[i] = Matrix_p_HeHe[index];
      }
    
        
      gsl_spline_init(spline_p_pp, Eap, Vector_p_pp, dimp+1);
      gsl_spline_init(spline_p_pHe, Eap, Vector_p_pHe, dimp+1);
      gsl_spline_init(spline_p_Hep, Eap, Vector_p_Hep, dimp+1);
      gsl_spline_init(spline_p_HeHe, Eap, Vector_p_HeHe, dimp+1);
      //
      gsl_spline_init(spline_p_pp_up, Eap, Vector_p_pp_up, dimp+1);
      gsl_spline_init(spline_p_pHe_up, Eap, Vector_p_pHe_up, dimp+1);
      gsl_spline_init(spline_p_Hep_up, Eap, Vector_p_Hep_up, dimp+1);
      gsl_spline_init(spline_p_HeHe_up, Eap, Vector_p_HeHe_up, dimp+1);
        
      for (unsigned int i = 0; i < dimEn; i++) {
	//cout << i << endl;
            
	double Eap_i = min(1e5,energy[i]);
	double Eap_i_tot = energy[i] + MpGeV;
            
	int i_ap = int(floor(log10(Eap_i/Eap[0])/DBlog));
	if (i_ap > dimap-1) i_ap = dimp-1;
	//double t = (Eel-Elept[i_lept])/(Elept[i_lept+1]-Elept[i_lept]);
            
	double valuefix = gsl_spline_eval(spline_p_pp, Eap_i, acc);
	double valueup = gsl_spline_eval(spline_p_pp_up, Eap_i, acc);
	double cs_pp = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_p_pHe, Eap_i, acc);
	valueup = gsl_spline_eval(spline_p_pHe_up, Eap_i, acc);
	double cs_pHe = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_p_Hep, Eap_i, acc);
	valueup = gsl_spline_eval(spline_p_Hep_up, Eap_i, acc);
	double cs_Hep = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_p_HeHe, Eap_i, acc);
	valueup = gsl_spline_eval(spline_p_HeHe_up, Eap_i, acc);
	double cs_HeHe = valuefix*(1-u) + valueup*u;  
            
	//cout << "test" << endl;
	spall_apel[coupleppr][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;
	//cout << "test" << endl;
	spall_apel[couplepHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;	
	//cout << "test" << endl;            
      }
      //cout << "end of i loop" << endl;	
        
    }


  }  else if (which_particle == 3) { //tertiary antiprotons

    //
    for (unsigned int j=0; j<dimEn; j++){
        
      double Epr = energy[j];
      double momentum = sqrt(Epr*Epr + 2.0*mp*Epr);
      int j_pr = int(floor(log10(Epr/ET[0])/DBprlog));
      if (j_pr > dimpr-1) j_pr = dimpr-1;
      double u = (Epr-ET[j_pr])/(ET[j_pr+1]-ET[j_pr]);
        
      for (int i = 0; i <= dimlept; i++) {
	int index = index_matrix(i,j_pr);
	Vector_3Ap_app[i] = Matrix_3Ap_app[index];
	Vector_3Ap_apHe[i] = Matrix_3Ap_apHe[index];
	//
	if (j_pr <= dimpr)
	  index = index_matrix(i,j_pr+1);
	else
	  index = dimpr;
	Vector_3Ap_app_up[i] = Matrix_3Ap_app[index];
	Vector_3Ap_apHe_up[i] = Matrix_3Ap_apHe[index];
      }
        
      gsl_spline_init(spline_3Ap_app, Eap, Vector_3Ap_app, dimap+1);
      gsl_spline_init(spline_3Ap_apHe, Eap, Vector_3Ap_apHe, dimap+1);
      //
      gsl_spline_init(spline_3Ap_app_up, Eap, Vector_3Ap_app_up, dimap+1);
      gsl_spline_init(spline_3Ap_apHe_up, Eap, Vector_3Ap_apHe_up, dimap+1);
        
      for (unsigned int i = 0; i < dimEn; i++) {
            
	double Eap_i = min(1e5,energy[i]);
	double Eap_i_tot = energy[i] + MpGeV;
            
	int i_ap = int(floor(log10(Eap_i/Eap[0])/DBlog));
	if (i_ap > dimap-1) i_ap = dimap-1;
	//double t = (Eel-Elept[i_lept])/(Elept[i_lept+1]-Elept[i_lept]);
            
	double valuefix = gsl_spline_eval(spline_3Ap_app, Eap_i, acc);
	double valueup = gsl_spline_eval(spline_3Ap_app_up, Eap_i, acc);
	double cs_pp = valuefix*(1-u) + valueup*u;
            
	valuefix = gsl_spline_eval(spline_3Ap_apHe, Eap_i, acc);
	valueup = gsl_spline_eval(spline_3Ap_apHe_up, Eap_i, acc);
	double cs_pHe = valuefix*(1-u) + valueup*u;
            
	spall_apel[coupleapap][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;
            
      }
        
    }


  } 
  cout << "end of routine" << endl;  
  return ;

}


//*********************************************************************************************************************************************

void TSpallationNetwork::InitXSecPohl(double factorelpos) {
    
  pair<int,int> coupleel(1001, -1000); // Electrons from protons
  pair<int,int> couplepos(1001, 1000); // Positrons from protons
  pair<int,int> coupleelHe(2004, -1000); // Electrons from He
  pair<int,int> coupleposHe(2004, 1000); // Positrons from He
    
  const int dimEn = energy.size();
    
  InitDataTablesPohl();
    
  // Huang & Pohl tables, arXiv:0711.2528,  http://cherenkov.physics.iastate.edu/gamma-prod/
    
  //ofstream outfile("elposspectr_P.dat", ios::out);
  double Emin_table_lepton = 0.01;
  double Emax_table_lepton = 1.e8;
    
  /* Formula for secondary lepton KINETIC energy bin: 
     
     int i = int((log(Ek)-log(Emin_table_lepton))/(log(Emax_table_lepton)-log(Emin_table_lepton))*201.0 + 1.0);
     // but be careful with the names of indexes. Also below.
     */
    
  const int dimleptHP = 201;
  const int dimlept = dimleptHP-1;
  const double DBlogHP = log(Emax_table_lepton/Emin_table_lepton)/(double)dimleptHP;
  double EleptHP[dimleptHP];
  for (int i=0; i < dimleptHP; i++) EleptHP[i] = exp(log(0.01)+double(i)*DBlogHP);
    
  /* Formula for primary proton (He) TOTAL energy:   */
  const int dimETHP = 374;
  double ET[374];
  for (int j = 0; j < 374; j++) {
    ET[j] = 1.24*pow(1.05,j)-mp;
    //cout << ET[j] << " " << ProdXsec[j] <<  endl;    
  }
  double DBlogpr = (log(ET[Nelectrons-1])-log(ET[0]))/(double)Nelectrons;
  //cout << "ProdXsec.size() = " << ProdXsec.size() << " " << ProdXsec.size()/2 << endl;
  //cout << "Matrix.size() = " << Matrix_El_pp.size() << " " << Matrix_El_pp.size()/Nelectrons << endl;
    
  gsl_spline *spline_ProdXsec_pp = gsl_spline_alloc(gsl_interp_cspline, Nelectrons);
  gsl_spline *spline_ProdXsec_He = gsl_spline_alloc(gsl_interp_cspline, Nelectrons);
    
  gsl_spline_init(spline_ProdXsec_pp, ET, &(ProdXsec[0]), Nelectrons);
  gsl_spline_init(spline_ProdXsec_He, ET, &(ProdXsec[Nelectrons]), Nelectrons);
    
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline_El_pp = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_Hep = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pp = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_Hep = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    
  gsl_spline *spline_El_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_El_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
  gsl_spline *spline_Pos_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    
  double Vectorlept_El_pp[dimlept+1];
  double Vectorlept_El_Hep[dimlept+1];
  double Vectorlept_Pos_pp[dimlept+1];
  double Vectorlept_Pos_Hep[dimlept+1];
  double Vectorlept_El_pp_up[dimlept+1];
  double Vectorlept_El_Hep_up[dimlept+1];
  double Vectorlept_Pos_pp_up[dimlept+1];
  double Vectorlept_Pos_Hep_up[dimlept+1];
    
  for (unsigned int j=0; j<dimEn; j++){
        
    double Epr = energy[j];
    double momentum = sqrt(Epr*Epr + 2.0*mp*Epr);
    int j_pr = int(log(Epr/ET[0])/DBlogpr)-1;
    if (j_pr > Nelectrons-2) j_pr = Nelectrons-2;
    double u = (Epr-ET[j_pr])/(ET[j_pr+1]-ET[j_pr]);
    double interp_prodxsec;
    double interp_prodxsecHe;
        
    if (Epr < ET[0]) {
      interp_prodxsec = ProdXsec[0];
      interp_prodxsecHe = ProdXsec[Nelectrons];
    }
    else if (Epr > ET[Nelectrons-1]) {
      interp_prodxsec = ProdXsec[Nelectrons-1];
      interp_prodxsecHe = ProdXsec[2*Nelectrons-1];
    }
    else {
      interp_prodxsec = gsl_spline_eval(spline_ProdXsec_pp, Epr, acc); 
      interp_prodxsecHe = gsl_spline_eval(spline_ProdXsec_He, Epr, acc);
    }
    //cout << Epr << " " << interp_prodxsec << " " << interp_prodxsecHe << endl;
        
    for (int i = 0; i <= dimlept; i++) {
      //int index = index_matrix(i,j_pr);
      Vectorlept_El_pp[i] = ElppPohl[i][j_pr];//[index];
      Vectorlept_El_Hep[i] = ElHepPohl[i][j_pr];//Matrix_El_Hep[index];
      Vectorlept_Pos_pp[i] = PosppPohl[i][j_pr];//Matrix_Pos_pp[index];
      Vectorlept_Pos_Hep[i] = PosHepPohl[i][j_pr];//Matrix_Pos_Hep[index];
      //
      //if (j_pr <= Nelectrons)
      //index = index_matrix(i,j_pr+1);
      //else
      //    index = Nelectrons;
      Vectorlept_El_pp_up[i] = ElppPohl[i][j_pr+1];//[index];
      Vectorlept_El_Hep_up[i] = ElHepPohl[i][j_pr+1];//Matrix_El_Hep[index];
      Vectorlept_Pos_pp_up[i] = PosppPohl[i][j_pr+1];//Matrix_Pos_pp[index];
      Vectorlept_Pos_Hep_up[i] = PosHepPohl[i][j_pr+1];//Matrix_Pos_Hep[index];
    }
        
    gsl_spline_init(spline_El_pp, EleptHP, Vectorlept_El_pp, dimlept+1);
    gsl_spline_init(spline_El_Hep, EleptHP, Vectorlept_El_Hep, dimlept+1);
    gsl_spline_init(spline_Pos_pp, EleptHP, Vectorlept_Pos_pp, dimlept+1);
    gsl_spline_init(spline_Pos_Hep, EleptHP, Vectorlept_Pos_Hep, dimlept+1);
        
    gsl_spline_init(spline_El_pp_up, EleptHP, Vectorlept_El_pp_up, dimlept+1);
    gsl_spline_init(spline_El_Hep_up, EleptHP, Vectorlept_El_Hep_up, dimlept+1);
    gsl_spline_init(spline_Pos_pp_up, EleptHP, Vectorlept_Pos_pp_up, dimlept+1);
    gsl_spline_init(spline_Pos_Hep_up, EleptHP, Vectorlept_Pos_Hep_up, dimlept+1);
        
    for (unsigned int i = 0; i < dimEn; i++) {
            
      double Eel = min(1e8,energy[i]);
      //double Eel_tot = energy[i] + MeleGeV;
            
      bool stopp = false;
      bool stopHe = false;
      if (Eel < EleptHP[0] || Eel > EleptHP[200] || Epr > ET[Nelectrons-1] || Epr < ET[0]) {
	spall_apel[coupleel][i][j] = 0.0;
	spall_apel[coupleelHe][i][j] = 0.0;
	spall_apel[couplepos][i][j] = 0.0;
	spall_apel[coupleposHe][i][j] = 0.0;
	//   continue;
	stopp = true;
	stopHe = true;
      }
            
      if (Eel+MeleGeV > Epr + mp) {
	spall_apel[coupleel][i][j] = 0.0;
	spall_apel[couplepos][i][j] = 0.0;
	stopp = true;
	//                continue;
      }
      if (Eel+MeleGeV > 4.0*(Epr + mp)) {
	spall_apel[coupleel][i][j] = 0.0;
	spall_apel[coupleelHe][i][j] = 0.0;
	spall_apel[couplepos][i][j] = 0.0;
	spall_apel[coupleposHe][i][j] = 0.0;
	stopp = true;
	stopHe = true;
      }
            
      int i_lept = int(floor(log(Eel/EleptHP[0])/DBlogHP));
      if (i_lept > dimlept-1) i_lept = dimlept-1;
            
            
            
      //double valuefix = gsl_spline_eval(spline_El_pp, Eel, acc);
      //double valueup = gsl_spline_eval(spline_El_pp_up, Eel, acc);
            
      double t;
      double valuefix;
      double valueup;
      double cs_pp;
      double cs_Hep;
      if (!stopp) {
                
	//t = (Eel-EleptHP[i_lept])/(EleptHP[i_lept+1]-EleptHP[i_lept]);
	valuefix = gsl_spline_eval(spline_El_pp, Eel, acc);//Vectorlept_El_pp[i_lept]*(1.0-t) + Vectorlept_El_pp[i_lept+1]*t;
	valueup = gsl_spline_eval(spline_El_pp_up, Eel, acc);//Vectorlept_El_pp_up[i_lept]*(1.0-t) + Vectorlept_El_pp_up[i_lept+1]*t;
                
	cs_pp = max(0.0,valuefix*(1-u) + valueup*u);
                
	//cout << Epr << " " << Eel << " " << cs_pp << endl;
                
	spall_apel[coupleel][i][j] = (Epr)*sqrt(1.0-1.0/pow(1.0+Epr/mp,2))*(cs_pp)*interp_prodxsec*factorelpos/1000.0;
                
	//valuefix = gsl_spline_eval(spline_Pos_pp, Eel, acc);
	//valueup = gsl_spline_eval(spline_Pos_pp_up, Eel, acc);
	valuefix = gsl_spline_eval(spline_Pos_pp, Eel, acc);//Vectorlept_Pos_pp[i_lept]*(1.0-t) + Vectorlept_Pos_pp[i_lept+1]*t;
	valueup = gsl_spline_eval(spline_Pos_pp_up, Eel, acc);//Vectorlept_Pos_pp_up[i_lept]*(1.0-t) + Vectorlept_Pos_pp_up[i_lept+1]*t;
	cs_pp = max(0.0,valuefix*(1-u) + valueup*u);
	//cout << cs_pp << " ";
                
	spall_apel[couplepos][i][j] = 4.0*(Epr)*sqrt(1.0-1.0/pow(1.0+Epr/mp,2))*(cs_pp)*interp_prodxsec*factorelpos/1000.0;   // H
                
      }
      if (!stopHe) {
	//valuefix = gsl_spline_eval(spline_El_Hep, Eel, acc);
	//valueup = gsl_spline_eval(spline_El_Hep_up, Eel, acc);
	valuefix = gsl_spline_eval(spline_El_Hep, Eel, acc);//Vectorlept_El_Hep[i_lept]*(1.0-t) + Vectorlept_El_Hep[i_lept+1]*t;
	valueup = gsl_spline_eval(spline_El_Hep_up, Eel, acc);//Vectorlept_El_Hep_up[i_lept]*(1.0-t) + Vectorlept_El_Hep_up[i_lept+1]*t;
	cs_Hep = max(0.0,valuefix*(1-u) + valueup*u);
                
	//cout << cs_Hep << " ";
                
	spall_apel[coupleelHe][i][j] = (Epr)*sqrt(1.0-1.0/pow(1.0+Epr/mp,2))*(cs_Hep)*interp_prodxsecHe*factorelpos/1000.0;	
                
                
                
	//valuefix = gsl_spline_eval(spline_Pos_Hep, Eel, acc);
	//valueup = gsl_spline_eval(spline_Pos_Hep_up, Eel, acc);
	valuefix = gsl_spline_eval(spline_Pos_Hep, Eel, acc);//Vectorlept_Pos_Hep[i_lept]*(1.0-t) + Vectorlept_Pos_Hep[i_lept+1]*t;
	valueup = gsl_spline_eval(spline_Pos_Hep_up, Eel, acc);//Vectorlept_Pos_Hep_up[i_lept]*(1.0-t) + Vectorlept_Pos_Hep_up[i_lept+1]*t;
	cs_Hep = max(0.0,valuefix*(1-u) + valueup*u);
	//cout << cs_Hep << endl;
                
	spall_apel[coupleposHe][i][j] = 4.0*(Epr)*sqrt(1.0-1.0/pow(1.0+Epr/mp,2))*(cs_Hep)*interp_prodxsecHe*factorelpos/1000.0;   // He
      }
            
    }
    //double sum = 0;
    //for (int i = 0; i < dimEn; i++) sum += spall_apel[coupleel][i][j];
    //cout << sum << endl;
  }
    
  //outfile.close();
    
}

void TSpallationNetwork::InitDataTablesPohl() {
  Nelectrons = 374;

  ifstream infile("data/espectra_eminus.decay.p.matrix.final.data", ios::in);
  if (!infile.is_open()) cerr<< "Problem opening eminus pr file" << endl;
    
    
  char* s = new char[20];
  for (int i = 0; i < 201; i++) {
    vector<double> dnde(374,0.0);
    for (int j = 0; j < Nelectrons; j++) {
      infile.get(s,16);
      s[11] = 'E';
      dnde[j] = atof(s);
    }
    infile.get(*s);
    ElppPohl.push_back(dnde);
  }
  infile.close();
    
  infile.open("data/espectra_eminus.decay.he.matrix.final.data", ios::in);
  if (!infile.is_open()) cerr<< "Problem opening eminus He file" << endl;
    
  for (int i = 0; i < 201; i++) {
    vector<double> dnde(374,0.0);
    for (int j = 0; j < 374; j++) {
      infile.get(s,16);
      s[11] = 'E';
      dnde[j] = atof(s);
    }
    infile.get(*s);
    ElHepPohl.push_back(dnde);
  }
  infile.close();
    
  infile.open("data/espectra_eplus.decay.p.matrix.final.data", ios::in);
  if (!infile.is_open()) cerr<< "Problem opening eplus pr file" << endl;
    
  for (int i = 0; i < 201; i++) {
    vector<double> dnde(374,0.0);
    for (int j = 0; j < 374; j++) {
      infile.get(s,16);
      s[11] = 'E';
      dnde[j] = atof(s);
    }
    infile.get(*s);
    PosppPohl.push_back(dnde);
  }
  infile.close();
    
  infile.open("data/espectra_eplus.decay.he.matrix.final.data", ios::in);
  if (!infile.is_open()) cerr<< "Problem opening eplus He file" << endl;
    
  for (int i = 0; i < 201; i++) {
    vector<double> dnde(374,0.0);
    for (int j = 0; j < 374; j++) {
      infile.get(s,16);
      s[11] = 'E';
      dnde[j] = atof(s);
    }
    infile.get(*s);
    PosHepPohl.push_back(dnde);
  }
  infile.close();
    
  infile.open("data/prodxsection.p.matrix.final.data", ios::in);
  if (!infile.is_open()) cerr<< "Problem opening eplus He file" << endl;
  ProdXsec = vector<double>(2*Nelectrons, 0);
  for (int j = 0; j < Nelectrons; j++) {
    infile.get(s,16);
    s[11] = 'E';
    ProdXsec[j] = atof(s);
  }
  infile.close();
    
  infile.open("data/prodxsection.he.matrix.final.data", ios::in);
  if (!infile.is_open()) cerr<< "Problem opening eplus He file" << endl;
    
  for (int j = 0; j < Nelectrons; j++) {
    infile.get(s,16);
    s[11] = 'E';
    ProdXsec[Nelectrons+j] = atof(s);
  }
  infile.close();
    
  delete [] s;
  return ;
}

//*********************************************************************************************************************************************
//*********************************************************************************************************************************************

// Inelastic cross sections
TInelasticCrossSection::TInelasticCrossSection(TGrid* Coord, Input* in, int uid, int K_electron, vector<TXSecBase*> xsecmodel) {

  if (in->feedback >1){ cout << "Inelastic xsec constructor " << endl;
    cout << in->spallationxsec << endl;
  }
  const double factor = 1e-27*Clight;

  vector<double> beta_vec = Coord->GetBetaEl();
  vector<double> gamma_vec = Coord->GetGammaEl();
  vector<double> e_vec = Coord->GetEk();

  if (in->feedback > 1) cout << "[debug] xsec test " << endl;
  
  if (in->spallationxsec == GalpropXSec) 
    xsec = xsecmodel[0]->GetXSec(uid, K_electron);
  else if (in->spallationxsec == Fluka)  {
    
    cout << "Fluka mode for inelastic xsec! " << endl;
    cout << "uid = " << uid << endl;
    /*if (!xsecmodel[1]->IsPresent(uid))
      xsec = xsecmodel[0]->GetXSec(uid, K_electron);
      else {

      xsecmodel[1]->InitTableInelastic(); 
      vector<double> energy1 = Coord->GetEk();
      for(unsigned int ip = 0; ip < energy1.size(); ip++) 
      xsec.push_back(xsecmodel[1]->GetTotalXSec(uid, energy1[ip])*factor*beta_vec[ip]*(1.0+He_abundance*xsecmodel[1]->GetHefactor()));
      }*/

    if (uid!=-999)
      xsec = xsecmodel[0]->GetXSec(uid, K_electron);
    else {		//pbar
      cout << "Fluka Antiproton total xsec VS galprop" << endl;
      vector<double> xsec_galprop = xsecmodel[0]->GetXSec(uid, K_electron);
      xsecmodel[1]->InitTableInelastic();
      vector<double> energy1 = Coord->GetEk();
      for(unsigned int ip = 0; ip < energy1.size(); ip++) {
	xsec.push_back(xsecmodel[1]->GetTotalApHXSec(energy1[ip])*factor*beta_vec[ip] + He_abundance*xsecmodel[1]->GetTotalApHeXSec(energy1[ip])*factor*beta_vec[ip]);
	cout << "E= " << energy1[ip] << " fluka on H -> " << xsecmodel[1]->GetTotalApHXSec(energy1[ip])*factor*beta_vec[ip] << " fluka on He -> " << xsecmodel[1]->GetTotalApHeXSec(energy1[ip])*factor*beta_vec[ip] << " galprop -> " << xsec_galprop[ip] << endl;
      }
 

    }


  }
  else if (in->spallationxsec == Webber03) {
    if (xsecmodel[1]->IsPresent(uid)) {
      vector<double> beta = Coord->GetBeta();
      vector<double> energy1 = Coord->GetEk();
      //	TWebber03::GetInstance()->Print(co);
      const double factor = 1e-27*Clight;
      for(unsigned int ip = 0; ip < energy1.size(); ip++) xsec.push_back(xsecmodel[1]->GetTotalXSec(uid, energy1[ip])*factor*beta[ip]*(1.0+He_abundance*xsecmodel[1]->GetHefactor()));
    }
  }

  if (in->feedback > 1) cout << "[debug] xsec test " << endl;

  //vector<double> energy1 = Coord->GetEk();
  //if (uid == 1000){
  //cout << "*** Galprop spallation cross section for positrons" << endl;
  //for(unsigned int ip = 0; ip < energy1.size(); ip++)
  //  cout << xsec[ip] << endl;
  //}
    
  //modified
  //Positron Annihilation  -- Still under developement!!

  if (uid == 1000) {
	 
    if (in->feedback >1) cout << endl << "*** Positron annihilation cross section" << endl;
    
    for(unsigned int ip = 0; ip < e_vec.size(); ip++) {
      
      double gamma = gamma_vec[ip];
      // cross section from P. A. M. Dirac, Proc. Camb. Phil. Soc. 26, 361 (1930). sigma = pi r_e^2 /(gamma+1) * [�some function of gamma��] 
      double annihilation_cross_section = 249.46/(gamma+1) * ( (gamma*gamma + 4*gamma + 1)/(gamma*gamma -1) * log(gamma + sqrt(gamma*gamma-1)) - (gamma+3)/(sqrt(gamma*gamma-1)) );  //mbarn
      if (in->feedback >1) cout << factor*beta_vec[ip]*annihilation_cross_section << endl;
	 
      xsec[ip] += (factor*beta_vec[ip])*annihilation_cross_section;
      //cout << "xsec test " << endl;
    }
  }
    
    
}

//*********************************************************************************************************************************************
//   Webber model for nucleus-nucleus xsec
//*********************************************************************************************************************************************

TWebber03::TWebber03() : TXSecBase() {
    
  ifstream file_to_read(Webber03Data.c_str(), ios::in);
  const int max_num_of_char_in_a_line = 512;
    
  int channel_temp;
  double energy_temp, xsec_double;
    
  while (file_to_read.good()) {
    file_to_read >> channel_temp;
        
    if (channel_temp == 0) {
      for (int i = 0; i < 15; i++) {
	file_to_read >> energy_temp;
	energy.push_back(energy_temp/1000.0);
      }
    }
    else {
            
      int uid1 = -1;   // Parent
      int uid2 = -1;   // Daughter
      TWebber03::convert(channel_temp, uid1, uid2);
            
      pair<int,int> couple = pair<int,int>(uid1,uid2);
      map< pair<int,int>, vector<double> >::iterator it = xsec.find(couple);
            
      if (it == xsec.end()) {
	xsec[couple] = vector<double>(energy.size(), 0.0);
	for (int i = 0; i < 15; i++) {
	  file_to_read >> xsec_double;
	  xsec[couple][i] = xsec_double;
	}
      }
    }
        
    file_to_read.ignore(max_num_of_char_in_a_line, '\n');
        
  }
    
  file_to_read.close();
    
    
  // Total XSec
  file_to_read.open(Webber03DataTotal.c_str(), ios::in);
    
  while (file_to_read.good()) {
    file_to_read >> channel_temp;
        
    if (channel_temp == 0) file_to_read.ignore(max_num_of_char_in_a_line, '\n');
    else {
            
      int uid1 = -1;   // Nucleus id
      TWebber03::convert_single(channel_temp, uid1);
      map< int, vector<double> >::iterator it = totalxsec.find(uid1);
            
      if (it == totalxsec.end()) {
	totalxsec[uid1] = vector<double>(energy.size(), 0.0);
	for (int i = 0; i < 15; i++) {
	  file_to_read >> xsec_double;
	  totalxsec[uid1][i] = xsec_double;
	}
      }
    }
        
    file_to_read.ignore(max_num_of_char_in_a_line, '\n');
        
  }
  file_to_read.close();
    
}

void TWebber03::Print(TGrid* co) {
  vector<double> energia = co->GetEk();
  ofstream outfile("webber_test_6012_5011.dat", ios::out);
  for (int i = 0; i < energia.size(); ++i) outfile << energia[i] << " " << GetXSec(pair<int,int>(6012,5011), energia[i])/**(1.0*He_abundance*GetHefactor())*/ << endl;
}

double TWebber03::GetXSec(pair<int,int> input, double en) {
    
  double result = 0.0;
    
  map< pair<int,int>, vector<double> >::iterator it = xsec.find(input);
    
  if (it == xsec.end()) {
    // cerr << "Channel not found" << endl;
    return -1;
  }
    
  int N = energy.size();
  if (en > energy.back()) return (*it).second.back();
    
  gsl_vector_view x = gsl_vector_view_array (&energy[0], N);
  gsl_vector_view y = gsl_vector_view_array (&((*it).second)[0], N);
    
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
    
  gsl_spline_init (spline, (&x.vector)->data, (&y.vector)->data, N);
    
  result = gsl_spline_eval (spline, en, acc);
    
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
    
  int l = 0;
  while (energy[l] < en) l++;
  double e1 = (energy[l]-en)/(energy[l]-energy[l-1]);
  double e2 = (en-energy[l-1])/(energy[l]-energy[l-1]);
    
  double result1 = (*it).second[l-1]*e1 + (*it).second[l]*e2;
    
  //  cout << "Comparison interpolation in Webber: " << result << " " << result1 << endl;
  return max(result,0.0);
}

double TWebber03::GetTotalXSec(int input, double en) {
    
  double result = 0.0;
    
  int N = energy.size();
  if (en > energy.back()) return totalxsec[input].back();
    
  gsl_vector_view x = gsl_vector_view_array (&energy[0], N);
  gsl_vector_view y = gsl_vector_view_array (&(totalxsec[input][0]), N);
    
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
    
  gsl_spline_init (spline, (&x.vector)->data, (&y.vector)->data, N);
    
  result = gsl_spline_eval (spline, en, acc);
    
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  int l = 0;
  while (energy[l] < en) l++;
  double e1 = (energy[l]-en)/(energy[l]-energy[l-1]);
  double e2 = (en-energy[l-1])/(energy[l]-energy[l-1]);
    
  double result1 = totalxsec[input][l-1]*e1 + totalxsec[input][l]*e2;
    
  //  cout << "Comparison interpolation in Webber total: " << result << " " << result1 << endl;
  return max(result,0.0);
}

//*********************************************************************************************************************************************
//   Galprop model for nucleus-nucleus xsec
//*********************************************************************************************************************************************


TGalpropXSec::TGalpropXSec(TGrid* co) : TXSecBase() {
    
  energy = co->GetEk();
  beta = co->GetBeta();
    
  data_filename.push_back("data/galprop_isotope_cs_updated.dat");
  data_filename.push_back("data/galprop_nucdata.dat");
  data_filename.push_back("data/galprop_p_cs_fits.dat");
  data_filename.push_back("data/galprop_eval_iso_cs_updated.dat");    // data-files to read
    
  file_no[0] = 0;
  file_no[1] = 1;
  file_no[2] = 2;
  file_no[3] = 3;
    
  int i,j,k,size;
  const int BufferSize=200;
  char readBuffer[BufferSize];
  ifstream data;
    
  for(j=0; j<N_DATA_FILES; j++)
    {
      data.open(data_filename[j].c_str());                    // open file if exists
      if(data.fail())
        {
	  cerr<<"read_nucdata: Error opening file "<<data_filename[j]<<endl;
	  exit(1);
        }
        
      while(!isspace(data.get()) && !data.eof())      // skip comments:
	data.getline(readBuffer,BufferSize,'\n');    // any symbol in 1st col. 
        
      for(i=0; i<3; data >> n_data[i++][j]);          // read array's sizes
      data.getline(readBuffer,BufferSize,'\n');       // skip the rest of line
        
      for(size=1, i=0; i<3; size*=n_data[i++][j]);    // allocate space
      //      data_file[j] = (float*) calloc( (size_t) size, (size_t) sizeof(float));
      data_file[j] = new float[size];
        
      for(k = 0; k < size && !data.eof();)            // read data loop
        {
	  while(!isspace(data.get()) && !data.eof())   // skip comments:
	    data.getline(readBuffer,BufferSize,'\n'); // any symbol in 1st col. 
	  for(i=0; i < n_data[0][j]; i++) data >> *(data_file[j]+k++);
	  data.getline(readBuffer,BufferSize,'\n');    // skip the rest of line
        }
      data.close();
    }
    
  set_sigma_cc();
    
  int ISS = -1;
  sigtap_cc(ISS);
  //  Print();
  return;
}

void TGalpropXSec::set_sigma_cc() { 
  int  cdr=99;
  set_sigma_(&cdr);
}

double TGalpropXSec::wsigma_cc(int IZ, int IA, int IZF, int IAF, double E) {
  return( wsigma_(&IZ,&IA,&IZF,&IAF,&E) );
}

double TGalpropXSec::yieldx_cc(int IZ, int IA, int IZF, int IAF, float E) {
  float CSmb;
  yieldx_(&IZ,&IA,&IZF,&IAF,&E,&CSmb);
  return( 1.*CSmb );
}

double TGalpropXSec::sighad_cc(int IS, double PA, double PZ, double TA, double TZ, double E) { 
  return( sighad_(&IS, &PA, &PZ, &TA, &TZ, &E) );
}

//double GalpropYields::pp_meson_cc(double Esec, double Pp1, int NA1, int NA2, int key1)
//{ 
//   return ( pp_meson_(&Esec,&Pp1,&NA1,&NA2,&key1) );
//}

void TGalpropXSec::sigtap_cc(int ISS) { 
  sigtap2_(&ISS);
}

//double TGalpropXSec::antiproton_cc(int key,double Pap1,double Pp1,int NZ1,int NA1,int NZ2,int NA2) { 
//  return ( antiproton_(&key, &Pap1, &Pp1, &NZ1, &NA1, &NZ2, &NA2) );
//}
/*
  void TGalpropXSec::Cleanup() {
  for(int j=0; j<N_DATA_FILES; j++) delete data_file[j];
  delete [] data_file;
  return;
  }
*/

vector<double> TGalpropXSec::GetXSec(int iz, int ia, int jz, int ja) {
  int IZ1,IA1,IZ3,IA3,kopt,info, K_electron =0;
  int galdef_network_par=0; // temporary solution; value 1 caused problems in nuc_package
  double branching_ratio,t_half;
    
  const int diz=3;                // maximum delta Z considered
  const double factor = Clight*1.e-27;
    
  kopt = cross_section_option;
    
  vector<double> spall = vector<double>(energy.size(), 0.0);
    
  // a loop over an intermediate state; final state must be as requested
  for (IZ1=(jz-diz>1) ? jz-diz: 1; IZ1<=iz && IZ1<=jz+diz; IZ1++) {
    for (IA1=(2*IZ1-4>ja) ? 2*IZ1-4: ja; IA1<ia && IA1<=2.5*IZ1+4.2; IA1++) {
            
      // channel selection procedure
      if (IA1 < IZ1 || ia-IA1 < iz-IZ1) continue;
            
      // IMOS20010816 line below
      branching_ratio = TGalpropXSec::nucdata(galdef_network_par,IZ1,IA1,K_electron,jz,ja, &IZ3,&IA3,&t_half);
            
      if (branching_ratio == 0) continue;
            
      t_half = t_half / year;
            
      // skip if long-lived intermediate state
      if(t_half>=t_half_limit                            // IMOS20010816
	 && 100*IZ1+IA1!=100*jz+ja && 100*IZ3+IA3!=100*jz+ja) continue;
            
      for(unsigned int ip = 0; ip < energy.size(); ip++) spall[ip] += (TGalpropXSec::isotope_cs(energy[ip]*1000.0,iz,ia,IZ1,IA1,kopt,&info)*branching_ratio*(1.0+ He_abundance*TGalpropXSec::He_to_H_CS_ratio(energy[ip],iz,ia,IZ1,IA1))*beta[ip]*factor);
            
    } //ja1 //jz1
  }
    
  // if (1000*iz+ia == 26056) {  
  //for(unsigned int ip = 0; ip < energy.size(); ip++) outf << 1000*iz+ia << " " << 1000*jz+ja << " " << energy[ip] << " " << spall[ip]/beta[ip]/factor << endl;
  //  }
     
  return spall;
}

//modified
void TGalpropXSec::Kcapture_cs(double Ek_GeV, int Zp, int Zt, double *attach, double *strip) {
    
  double gam = 1.+Ek_GeV/amu, sT =1.e27 * 8./3.*Pi*pow(Rele,2), //AWS20010829
    beta,Mbet,Nbet,a,fcor;
    
  beta = sqrt(1.-pow(gam,-2));
  Mbet = 4./3.+gam*(gam-2.)/(gam+1.)*(1.-log((1.+beta)/(1.-beta))/2./beta/pow(gam,2));
  Nbet = pow(beta,-3) *((-4.*gam +34. -63./gam +25./pow(gam,2) +8./pow(gam,3))/15.
			-(gam-2.)*(gam-1.)/2./beta/pow(gam,3) *log((1.+beta)/(1.-beta)));
  a = Zp*ALPHAf;
  fcor = pow(a,2.*(sqrt(1.-a*a)-1.)) *exp(-2.*a/beta*acos(a)) *(1.+Pi*a*Nbet/Mbet);
  //    factor 1.202 accounts for contribution of states higher than 1s
  *attach = 1.202 *fcor *3./2.*sT*pow(1.*Zp,5)*Zt*pow(ALPHAf,4) *beta*gam *pow(gam-1.,-3)*Mbet;
  *strip = 3./2.*sT*pow(Zp*beta*ALPHAf,-2) 
    *0.285*(2.*log(2./0.219*beta*gam/Zp/ALPHAf)-pow(beta,2))*Zt*(Zt+1.);
  //return;
    
}


//modified
vector<double> TGalpropXSec::GetXSec(int uid, int K_electron) { //modified
  int Z = -1000;
  int A = -1000;
  Utility::id_nuc(uid, A, Z);
  /*
    if (uid == 28058) {
    for (unsigned int i = 0; i < energy.size(); i++) outf << energy[i] << "\t";
    outf << endl;
    }
    outf << uid << endl;
  */
  vector<double> result = vector<double>(energy.size(), 0.0);
    
  if (A == 0) return result;
    
  int target = 1;
  const bool protons = (A == 1 && Z == 1);
  const bool antiprotons = (A == 1 && Z == -1);
  if (protons) {
    A = 4;
    Z = 2;
  }
  if (antiprotons) {
    A = 4;
    Z = 2;
    target = -1;
  }
    
  double PP_inel = 0.0;
  double PA_inel = 0.0;
  double aPP_non = 0.0;
  double aPA_non = 0.0;
  double aPP_ann = 0.0;
  double aPA_ann = 0.0;
    
  const double factor = 1.e-27*Clight;
    
  const double Hecontribution = He_abundance*TGalpropXSec::He_to_H_CS_totratio(A);
    
  for (unsigned int k = 0; k < energy.size(); ++k) {
    TGalpropXSec::nucleon_cs(2, energy[k], target, Z, A, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);  // Galprop CS
        
    if (protons) result[k] = (PP_inel + He_abundance*PA_inel)*factor*beta[k];
    else if (antiprotons) result[k] = factor*(aPP_non + aPP_ann + He_abundance*(aPA_non + aPA_ann))*beta[k];
    else result[k] = PA_inel*factor*beta[k]*(1.0 + Hecontribution);
        
    //modified
    if (K_electron >= 0) {
            
      double attach_H=0., attach_He=0., strip_H=0., strip_He=0.;
            
      //modified
      TGalpropXSec::Kcapture_cs(energy[k],Z,1,&attach_H,&strip_H);
      TGalpropXSec::Kcapture_cs(energy[k],Z,2,&attach_He,&strip_He);
            
      //modified
      if (K_electron == 0) 
	result[k] += ( factor*beta[k]*(attach_H + He_abundance*attach_He) ); //The naked nucleus unstable for EC may attach an electron
      if (K_electron == 1)
	result[k] += ( factor*beta[k]*(strip_H + He_abundance*strip_He) );   //The "dressed" nucleus unstable for EC may lose an electron and become naked again 
    }
  }
  //  for(unsigned int ip = 0; ip < energy.size(); ip++) outf << result[ip]/factor/beta[ip] << "\t";
  //outf << endl;
  //if (uid == 1001) outf.close();
  return result;
    
}

double TGalpropXSec::He_to_H_CS_ratio(double E1,int IZI,int IAI,int IZF,int IAF) {
  double X,Y,a,b;
  double AMU,DELTA,FZI;
  double  E[]= {-999,  0.43,0.73,1.51};  // zeroth element dummy to conform to F77 original
  double d[]= {-999,  2.45,3.45,4.40};
  double am[]= {-999,  0.082,0.047,0.031};
  double Z[]= {-999,  6.,8.,26.};
  double F[]= {-999, -0.76,-0.41,+1.00};
    
  double CSratio    = 0.;
  if(IAI <= IAF) return 0;
    
  if(E1 < E[2])          // interpolation/extrapolation
    {
      AMU  = FI(E1,E[1],E[2],am[1],am[2]);
      DELTA= FI( log(E1), log(E[1]), log(E[2]),d[1],d[2]);  //log scale
    } else
    {
      AMU  = FI(E1,E[2],E[3],am[2],am[3]);
      DELTA= FI( log(E1), log(E[2]), log(E[3]),d[2],d[3]);  //log scale
    }    
  if(E1 < E[1])          // asymptotics E1->0
    {
      AMU  = am[1];
      DELTA=  d[1];
    }             
  if(E1 > E[3])          // asymptotics E1->inf
    {
      AMU  = am[3];
      DELTA=  d[3];
    }             
    
  if(IZI < Z[2]) FZI = FI( log(IZI*1.), log(Z[1]), log(Z[2]),F[1],F[2]); // inter-/extra-polation on log scale
  else           FZI = FI( log(IZI*1.), log(Z[2]), log(Z[3]),F[2],F[3]);   
    
  if(IZI-IZF <= IZI/2) CSratio = exp(AMU*pow(fabs(IZI-IZF-FZI*DELTA),1.43));
  else                 // if (IZI-IZF) > IZI/2, then another approximation (requires continuity in Z)
    {
      X =   exp(AMU*pow(fabs(IZI/2-  FZI*DELTA), 1.43));      //value
      Y = X-exp(AMU*pow(fabs(IZI/2-1-FZI*DELTA), 1.43));      // derivative in dZ
      b = IZI/2*Y/X;
      a = X/pow(IZI/2,b);
      CSratio = a*pow(IZI-IZF,b);
    }
  return CSratio;
}

double TGalpropXSec::nucdata(int ksp,int iz,int ia,int K_electron,int izf,int iaf,int* izl,int* ial,double* To) {
  int i,j,k,l,m,n,iy,iz0,ia0,iz4,ia4,iz5,ia5,iw[121],
    nksp=ksp*n_data[0][file_no[1]]*n_data[1][file_no[1]];
  float w[2][121], *decay = data_file[file_no[1]]+nksp;
  double b, xxx;
    
  // STABLE & LONG-LIVED ISOTOPES (numbers in the table are the proton numbers)
  // The long-lived isotopes listed in "longliv" table are included as stable;
  int stable[64][3] = {      // second index changes faster
    1,  0,  1,     1,  0,  1,     1,  0,  2,     2,  0,  2, //  A = 1- 4
    0,  0,  0,     3,  0,  3,     3,  0,  4,     0,  0,  0, //  A = 5- 8
    4,  0,  4,     4,  0,  5,     5,  0,  5,     6,  0,  6, //  A = 9-12
    6,  0,  6,     6,  0,  7,     7,  0,  7,     8,  0,  8, //  A =13-16
    8,  0,  8,     8,  0,  8,     9,  0,  9,    10,  0, 10, //  A =17-20
    10,  0, 10,    10,  0, 11,    11,  0, 11,    12,  0, 12, //  A =21-24
    12,  0, 12,    12,  0, 13,    13,  0, 13,    14,  0, 14, //  A =25-28
    14,  0, 14,    14,  0, 14,    15,  0, 15,    14,  0, 16, //  A =29-32
    16,  0, 16,    16,  0, 16,    17,  0, 17,    16, 17, 18, //  A =33-36
    17,  0, 18,    18,  0, 18,    18,  0, 19,    18, 19, 20, //  A =37-40
    19,  0, 20,    18,  0, 20,    20,  0, 20,    20,  0, 22, //  A =41-44
    21,  0, 21,    20,  0, 22,    22,  0, 22,    20,  0, 22, //  A =45-48
    22,  0, 23,    22, 23, 24,    23,  0, 24,    24,  0, 24, //  A =49-52
    24,  0, 25,    24, 25, 26,    25,  0, 26,    26,  0, 28, //  A =53-56
    26,  0, 27,    26,  0, 28,    27,  0, 28,    26, 27, 28, //  A =57-60
    28,  0, 28,    28,  0, 28,    28,  0, 29,    28,  0, 28  //  A =61-64
  };
    
  // LONG-LIVED ISOTOPES (>~1y): Zi.Ai T_1/2(y) Zf.Af - [][][0] half-life shown for naked nucleus
  int nll = 25;                                 // - [][][1] half-life shown for H2-like atoms
  float longliv[25][2][3] = {   // third index changes faster
    1.03,  12.33,    2.03,     //  3H (b-) 3He   100% [ToI]
    1.03,  12.33,    2.03,     // no EC
        
    4.07,  0.,       4.07,     // stable
    4.07,  0.1459,   3.07,     //  7Be(EC) 7Li   100% [ToI]
        
    4.10,  1.60e6,   5.10,     // 10Be(b-)10B    100% [ToI]
    4.10,  1.60e6,   5.10,     // no EC
        
    6.14,  5.73e3,   7.14,     // 14C (b-)14N    100% [ToI]
    6.14,  5.73e3,   7.14,     // no EC
        
    11.22,  4.80e3,  10.22,     // 22Na(b+)22Ne        [M98]
    11.22,  2.60e0,  10.22,     // 22Na(EC?)22Ne       [ToI] T1/2(Lab)=2.60e0 y
        
    13.26,  9.10e5,  12.26,     // 26Al(b+)26Mg        [M98]
    13.26,  4.075e6, 12.26,     // 26Al(EC)26Mg        [M98] T1/2(Lab)=7.4e5 y [ToI]
        
    14.32,  172.,    16.32,     // 32Si(2b-)32S   100% [ToI] Si-P -S 
    14.32,  172.,    16.32,     // no EC
        
    17.36,  3.07e5,  18.36,     // 36Cl(b-)36Ar        [ToI]
    17.36,  1.58e7,  16.36,     // 36Cl(EC)36S         [ToI] T1/2(Lab)=3.01e5 y
        
    18.37,  0.,      18.37,     // stable
    18.37,  0.1,     17.37,     // 37Ar(EC)37Cl   100% [ToI] T1/2(Lab)=35.04 d
        
    18.39,  2.69e2,  19.39,     // 39Ar(b-)39K    100% [ToI]
    18.39,  2.69e2,  19.39,     // no EC
        
    19.40,  1.43e9,  20.40,     // 40K (b-)40Ca   89.3%[ToI] T1/2(Lab)=1.277e9 y incl 10.7% ECb+
    19.40,  1.43e9,  20.40,     // no EC
        
    20.41,  0.,      20.41,     // stable
    20.41,  1.03e5,  19.41,     // 41Ca(EC)41K    100% [ToI]
        
    18.42,  32.9,    20.42,     // 42Ar(2b-)42Ca  100% [ToI] Ar-K -Ca
    18.42,  32.9,    20.42,     // no EC
        
    22.44,  0.,      22.44,     // stable
    22.44,  49.,     20.44,     // 44Ti(ECb+)44Ca 100% [ToI] Ti(EC)Sc(b+)Ca
        
    23.49,  0.,      23.49,     // stable
    23.49,  0.903,   22.49,     // 49V (EC)49Ti   100% [ToI] 
        
    24.51,  0.,      24.51,     // stable
    24.51,  0.076,   23.51,     // 51Cr(EC)51V   <100% [ToI] 
        
    25.53,  0.,      25.53,     // stable
    25.53,  3.74e6,  24.53,     // 53Mn(EC)53Cr   100% [ToI]
        
    25.54,  6.30e5,  26.54,     // 54Mn(b-)54Fe        [W98]
    25.54,  0.855,   24.54,     // 54Mn(EC)54Cr        [ToI] T1/2(Lab)=312.3 d
        
    26.55,  0.,      26.55,     // stable
    26.55,  2.73e0,  25.55,     // 55Fe(EC)55Mn   100% [ToI]
        
    28.56,  4.00e4,  26.56,     // 56Ni(2b+)56Fe <100% [F99] Ni-Co-Fe
    28.56,  0.1,     26.56,     // 56Ni(ECb+)56Fe      [ToI] T1/2(Lab)=~30 d Ni(EC)Co(b+)Fe
        
    27.57,  0.,      27.57,     // stable
    27.57,  0.744,   26.57,     // 57Co(EC)57Fe   100% [ToI]
        
    28.59,  0.,      28.59,     // stable
    28.59,  7.60e4,  27.59,     // 59Ni(EC)59Co  <100% [ToI] [B76]
        
    27.60,  5.27e0,  28.60,     // 60Co(b-)60Ni   100% [ToI]
    27.60,  5.27e0,  28.60,     // no EC
        
    26.60,  1.50e6,  27.60,     // 60Fe(b-)60Co   100% [ToI]
    26.60,  1.50e6,  27.60,     // no EC
        
    28.63,  1.00e2,  29.63,     // 63Ni(b-)63Cu   100% [ToI]
    28.63,  1.00e2,  29.63      // no EC
  };
  // K-capture nuclei - factor of 2 because only 1 electron
  for(i=0;i<nll;i++) if(longliv[i][0][1]!=longliv[i][1][1]) longliv[i][1][1] *=2.; 
    
  // BOUNDARY NUCLEI 
  // on the left side from the left boundary, the proton-emission is assumed;
  // on the right side from the right boundary, the neutron-emission is assumed.
  // ZB.AB; left boundary[][0]: Nn=0(1)28; right boundary[][1]: Np=1(1)28 
  int nb = 29;
  float boundary[29][2] = {  // second index changes faster
    1.01,  1.04,    3.04,  2.08,
    3.05,  3.11,    6.09,  4.14,
    6.10,  5.15,    8.13,  6.16,
    8.14,  7.21,   10.17,  8.22,
    12.20,  9.24,   12.21, 10.26,
    14.24, 11.30,   14.25, 12.30,
    15.27, 13.34,   16.29, 14.35,
    16.30, 15.38,   18.33, 16.40,
    20.36, 17.43,   20.37, 18.46,
    20.38, 19.48,   22.41, 20.51,
    22.42, 21.52,   24.45, 22.53,
    24.46, 23.55,   24.47, 24.59,
    25.49, 25.63,   26.51, 26.65,
    27.53, 27.69,   27.54, 28.69,
    28.56, 00.99
  };
    
  b = *To = *izl = *ial = 0;
  if(iz <= 0 || ia <= 0) return(0.);    // check against negative numbers,
  if(iz*ia > 1 && iz >= ia) return(0.); // non-existed nuclei,
  if(2864 < fnuc(iz,ia)) return(0.);    // Ni64 is the heaviest nucleus
  if(64 < ia) return(0.);               // A=64 is the maximal atomic number
    
  // CHECK FOR NUCLEI OUTSIDE THE BOUNDARIES (p/n decay)
  iz0 = iz;
  ia0 = ia;
  if(ia>inuc(boundary[iz-1][1]-iz)) ia0=inuc(boundary[iz-1][1]-iz); // n -decay
  if(29>ia-iz) if(ia>inuc(modf(boundary[ia-iz][0], &xxx)))          // p -decay
		 { 
		   iz0=(int)boundary[ia-iz][0];
		   ia0=inuc(boundary[ia-iz][0]-iz0);
		 }
    
  for(i=0; i<121; iw[i++]=-1)  for(j=0; j<2; w[j++][i]=0.);
    
  // SEARCH FOR A SPECIAL CASE (non beta decay)
  for(i=0; i<n_data[1][file_no[1]]; i++)
    if(fnuc(iz0,ia0) == inuc(*(decay +i*n_data[0][file_no[1]] +0)))
      {
	iw[0]  = i;            // if found, save the line number
	w[1][0]= 1.00;         // assign 1 to the branching ratio
	break;
      }
    
  // STANDARD CASE (beta decay & long-lived isotopes)
  if(iw[0] < 0)
    {
      iz5 = iz0;
      ia5 = ia0;
      // *** BETA DECAY ***
      if(iz0 > stable[ia0-1][2]) iz5 = stable[ia0-1][2];   // b+ decay
      if(iz0 < stable[ia0-1][0]) iz5 = stable[ia0-1][0];   // b- decay
      // *** LONG-LIVED ISOTOPES (>~1 y) ***
      for(i=0; i<nll; i++)
	if(fnuc(iz5,ia5) == inuc(longliv[i][K_electron][0]))
	  {
	    *izl = iz5;
	    *ial = ia5;
	    *To = longliv[i][K_electron][1]*year;
	    if(!*To) *izl=*ial=0;
	    iz5 = (int) longliv[i][K_electron][2];
	    ia5 = inuc(longliv[i][K_electron][2]-iz5);
	    break;
	  }
      if(fnuc(izf,iaf)==fnuc(iz5,ia5) || fnuc(izf,iaf)==fnuc(*izl,*ial)) b = 1.;
      if(fnuc(iz0,ia0) == fnuc(*izl,*ial)) *izl = *ial = 0;
      return(b);
    }
    
  // DEVELOPING A NETWORK OF DECAYS
  for(l=-1, m=0, ia4=1, i=0; i<4; ia4 =(int) pow(3.,++i))
    {
      for(l+=ia4, iy=0, n=0; n<ia4; n++, m++)
        {                                                      // check if there is
	  if(iw[m] < 0) continue;                             // a required channel
	  for(w[0][m]=0., k=2; k<8; k+=2)
            {
	      w[0][l+3*n+k/2]                                  // store sec.nuclei
                =*(decay +iw[m]*n_data[0][file_no[1]] +k-1);
	      w[1][l+3*n+k/2]                                  // store branchings
                =*(decay +iw[m]*n_data[0][file_no[1]] +k)*w[1][m];
	      for(j=0; j<iw[m]; j++)                           // check if sec.nucleus
		if(*(decay +iw[m]*n_data[0][file_no[1]] +k-1) // also develops a
		   == *(decay +j*n_data[0][file_no[1]] +0))   // network of decays
		  {
		    iw[l+3*n+k/2] = j;                         // store such a nucleus
		    iy = l+3*n+k/2;
		  }  ///printf("%d %d %d %d %d %d\n",l,n,m,k,j,iy);
            }
        }
      if(iy == 0) break;
    }
    
  // CHECK FOR STABILITY OF THE FINAL NUCLEI
  for(k=0; k<=l+3*n; k++)
    {
      *To = *izl = *ial = 0;
      if(w[0][k] == 0.) continue;
      iz4 = (int) w[0][k];
      ia4 = inuc(w[0][k]-iz4);
      iz5 = iz4;
      ia5 = ia4;
      // *** BETA DECAY ***
      if(iz4 > stable[ia4-1][2]) iz5 = stable[ia4-1][2];   // b+ decay
      if(iz4 < stable[ia4-1][0]) iz5 = stable[ia4-1][0];   // b- decay
      // *** LONG-LIVED ISOTOPES (>~1 y) ***
      for(i=0; i<nll; i++)
        {
	  if(fnuc(iz5,ia5) != inuc(longliv[i][K_electron][0])) continue;
	  *izl = iz5;
	  *ial = ia5;
	  *To = longliv[i][K_electron][1]*year;
	  if(!*To) *izl=*ial=0;
	  iz5 = (int) longliv[i][K_electron][2];
	  ia5 = inuc(longliv[i][K_electron][2]-iz5);
	  break;
        }
      if(fnuc(izf,iaf) == fnuc(*izl,*ial) || fnuc(izf,iaf) == fnuc(iz5,ia5))
	return(w[1][k]);
    }
  return(b);
}

double TGalpropXSec::isotope_cs(double emev,int iz,int ia,int izf,int iaf,int kopt,int* info) {
  int a1,a2,i,j,size, itable=0, info1;
  double e1,y,err2,xi1,xi2, f1=0., f2=0., T[11], a[3]={1.,1.,0.}, b[6];
  float *cs_data = data_file[file_no[0]], *p_cs = data_file[file_no[2]];
  double *tp=T;
  double ej;
  double CSmb = 0.0;
    
  e1 = emev;
  *info = kopt;
    
  // CHECK if user wants to use specific program (the value of "kopt")
  if(kopt == 1) CSmb = wsigma_cc(iz,ia,izf,iaf,emev);       // Webber's code IMOS20020502
  if(kopt == 2) CSmb = yieldx_cc(iz,ia,izf,iaf,e1);         // TS code       IMOS20020502
  CSmb = max(0.,CSmb);
  if(kopt == 1 || kopt == 2) return(CSmb);
    
  a1 = fnuc(iz, ia);
  a2 = fnuc(izf,iaf);
    
  // EVALUATED CROSS SECTIONS
    
  if(kopt == 12 || kopt == 22)
    {      
      CSmb = eval_cs(emev,a1,a2,&info1);
      if (info1 > 0) return(max(0.,CSmb));
      kopt--;                                          // if evaluation doesn't exist,  
    }                                                   // try other options
    
  // if user wants, use THE CROSS SECTION FITS
    
  if(kopt == 11 || kopt == 21)
    {      
      // special cases: Be, B - recursion calls
      if(izf != 0)
        {
	  // A = 10
	  if(10 == iaf)  // B10 = B10 + C10 = a10 - Be10
            {
	      b[0] = isotope_cs(emev,iz,ia,0,iaf,21,&j);
	      if(j == -21)
                {                                          // B10
		  if(510 == a2) 
                    {
		      b[0]-=isotope_cs(emev,iz,ia,4,iaf,21,&j);
		      return(max(0.,b[0]));
                    }
		  if(5 < izf) return(0.);                // C10, =0
                }
            }
	  // A = 11
	  if(11 == iaf)  // B11 = a11 = Be11 + B11 + C11
            {
	      b[0] = isotope_cs(emev,iz,ia,0,iaf,21,&j);
	      if(j == -21) 
                {
		  if(511 == a2) return(max(0.,b[0]));     // B11
		  return(0.);                             // =0 for the rest
                }
            }
        }
      // straight search in the table
        
      for(i=0; i<n_data[1][file_no[2]]-1; i++, p_cs+=n_data[0][file_no[2]]) // -1 fixes the reading error in the line below
	if(a1 == inuc(*p_cs) && a2 == inuc(*(p_cs+1)))
	  {
	    for(p_cs+=2, j=0; j<6; b[j++]=*p_cs++);    // take the parameters
	    if(b[0] >= 0.)                             // if positive use fit
	      {
		*info=-kopt;
		if(emev < b[5]) return(0);              // fitting function
		b[0]*=(1.+sin(b[1]*pow(log10(emev),1.*b[2]))*exp(-b[3]*(emev-b[4])));
		return(max(0.,b[0]));
	      }
	    kopt = (int)(-b[0]+0.1);                   // negative b[0] gives kopt
	  }
      if(izf == 0) return(0.);
    }
    
  // CHECK if user wants to use specific program (the value of "kopt")
  if(kopt == 1) CSmb = wsigma_cc(iz,ia,izf,iaf,emev);       // Webber's code IMOS20020502
  if(kopt == 2) CSmb = yieldx_cc(iz,ia,izf,iaf,e1);         // TS code       IMOS20020502
  CSmb = max(0.,CSmb);
  if(kopt == 1 || kopt == 2) return(CSmb);
    
  // STARTING THE ALGHORITHM
    
  for(i=0; i<11; T[i++] = 0.);
    
  // CHECK the array: cs_data (is there a channel we are looking for ?)
    
  for(size=1, i=0; i<3; size*=n_data[i++][file_no[0]]);
  for(tp = T, i=0; i<size; i+=n_data[0][file_no[0]], tp = T, f1=0., f2=0.)
    {
      if(a1 != inuc(*(cs_data+i)))   continue;
      if(a2 != inuc(*(cs_data+i+1))) continue;
        
      // if there is such a channel then the LEAST-SQUARE FIT
        
      itable++;
      if(*(cs_data+i+4) < 0.) *(cs_data+i+4) *= -*(cs_data+i+3);  // calc.abs.err.
      err2 = pow(*(cs_data+i+4),2);                               // err^2
        
      y = *(cs_data+i+3);                                         // cs measured
      ej = *(cs_data+i+2);                                        // @ energy
        
      if(kopt/10 != 2) f1=wsigma_cc(iz,ia,izf,iaf,ej);             // Webber IMOS20020502
      if(kopt/10 != 1) f2=yieldx_cc(iz,ia,izf,iaf,*(cs_data+i+2)); // TS     IMOS20020502
        
      // calculations of the separate terms:
      *tp++ += f1*y /err2;       // Webber
      *tp++ += f1*f1/err2;
      *tp++ += f2*y /err2;       // TS
      *tp++ += f2*f2/err2;
      *tp++ += y    /err2;       // const cs
      *tp++ += 1.   /err2;
        
      // calculation of terms for the Xi2 estimates
      *tp++ += y*y    /err2;
      *tp++ += 2.*f1*y/err2;
      *tp++ += f1*f1/err2;
      *tp++ += 2.*f2*y/err2;
      *tp   += f2*f2/err2;
        
      // calculation of renormalization coefficients 
      for(j=0; j<3; j++) {
	a[j]= (T[2*j+1] != 0.) ? T[2*j]/T[2*j+1]: a[j];
      }
    }
    
  if(kopt == 3 && a[2] != 0.) return(a[2]);                  // const cr.sect.
  if(kopt/10 == 1) CSmb = wsigma_cc(iz,ia,izf,iaf,emev);     // Webber code IMOS20020502
  if(kopt/10 == 2) CSmb = yieldx_cc(iz,ia,izf,iaf,e1);       // TS code     IMOS20020502
  if(kopt/10 == 1 || kopt/10 == 2) return(max(0.,CSmb*a[kopt/10-1]));
    
  // CHOOSE THE BEST APPROXIMATION (kopt = 0)
  if(itable < 2)                                         // no data or 1 pt.
    {  
      *info = itable;
      CSmb = a[0]*wsigma_cc(iz,ia,izf,iaf,emev);          // use Webber code     IMOS20020502
      if(CSmb <= 0.) 
        {                                                   // if W-code give 0,
	  CSmb = a[1]*yieldx_cc(iz,ia,izf,iaf,e1);         // take the TS approx. IMOS20020502
	  if(CSmb != 0. && itable == 1) *info = 2;
        }
    } else                                                 // data exists
    {
      xi1= T[6] -a[0]*T[7] +a[0]*a[0]*T[8];               // Xi2 evaluation 1
      xi2= T[6] -a[1]*T[9] +a[1]*a[1]*T[10];              // Xi2 evaluation 2
      if(xi1 < xi2)
        {
	  *info = 1;
	  CSmb = a[0]*wsigma_cc(iz,ia,izf,iaf,emev);       // renorm. Webber approx. IMOS20020502
        } else
        {
	  *info = 2;
	  CSmb = a[1]*yieldx_cc(iz,ia,izf,iaf,e1);         // renorm. TS approx.     IMOS20020502
        }
    }
    
  return(max(0.,CSmb));
}

double TGalpropXSec::eval_cs(double emev,int za1,int za2,int* info) {
  int i,size;
  float *eval = data_file[file_no[3]];
  double x[2]={-1.e10,1.e10},y[2]={0.,0.};
    
  // CHECK the array: eval (is there a channel we are looking for ?)
    
  for(size=1, i=0; i<3; size*=n_data[i++][file_no[3]]);
  for(*info=0, i=0; i<size; i+=n_data[0][file_no[3]])
    {
        
      if(za1 != inuc(*(eval+i)))   continue;
      if(za2 != inuc(*(eval+i+1))) continue;
        
      if(x[0] < *(eval+i+2) && *(eval+i+2) <= emev)   // find lower energy pt
        { 
	  x[0] = *(eval+i+2); 
	  y[0] = *(eval+i+3);
        }
      if(emev <= *(eval+i+2) && *(eval+i+2) < x[1])   // find higher energy pt 
        { 
	  x[1] = *(eval+i+2); 
	  y[1] = *(eval+i+3);
        }
    }
    
  if(x[0]*x[1] < -1.e19) { *info = -1; return(0.); } // no evaluation found, return 0
    
  if(x[0] <   0.) { *info = 1; return(0.); }         // no lower grid pt, return 0
  if(x[1] > 9.e9) { *info = 2; return(y[0]); }       // no higher grid pt, extrapolate
    
  if(x[1]-x[0] == 0.) { *info = 3; return(y[1]); }   // emev falls exactly on the grid
    
  for(*info = 4, i=0; i<2; i++) x[i] = log10(x[i]);
  return(y[0]+(log10(emev)-x[0])*(y[1]-y[0])/(x[1]-x[0]));// interpolate
}

void TGalpropXSec::nucleon_cs(int option, double Ek, int Zp, int Zt, int At, double *PP_inel,double *PA_inel,double *aPP_non,double *aPA_non,double *aPP_ann,double *aPA_ann) {
    
  *PP_inel= *PA_inel= *aPP_non= *aPA_non= *aPP_ann= *aPA_ann=0.;
  if(Ek <= 0.) return;
    
  double U,Cp,C1,s2,Z,A, aPP_inel=0.,aPA_inel=0., Emev=Ek, b0,Fcorr,rN,s0,p1,p2,p3,p4,p5,x,f1,f2;
  double PZ=fabs(1.*Zp),Em=1000.*Ek, TZ=Zt, TA=At;
  int ISS=2;
    
  //## Proton-Proton INELASTIC cross section, mb [TN83,p.234]
  if(Ek > 0.3)
    {
      U = log((Ek+mp)/200.);
      *PP_inel = 32.2*(1.+0.0273*U);
      if(U >= 0.) *PP_inel += 32.2*0.01*U*U;
      if(Ek < 3.) *PP_inel /= 1.+0.00262/pow(Ek,17.9+13.8*log(Ek)+4.41*pow(log(Ek),2));
    }
  if(Zp*At == 1) return;
    
  //## Proton-Nucleus INELASTIC cross section, mb
  switch(option) {
  case 0: // [L83]
    C1 = (At == 4) ? 0.8 : 1.;                  // a correction for He
    if(At == 9) C1 = 1.+0.75*exp(-Ek/0.075);    // a correction for Be
    *PA_inel = C1 *45. *pow(TA,0.7) *(1.+0.016*sin(5.3-2.63*log(TA)));
    if(Ek < 5.) *PA_inel *= 1.-0.62 *exp(-Ek/0.2) *sin(10.9/pow(Ek*1.e3,0.28));
    if(At == 4) *PA_inel = (Ek > 0.01) ?        // pHe, my fit
		  111.*(1.-(1.-sin(9.72*pow(log10(Ek*1000.),0.319)-4.14))*exp(-3.84*(Ek-0.1))) : 0.;
    break;
            
  case 1: // Zt>5 [WA96], Zt<=5 [BP01]
    if(Zt>5) {
      b0 = 2.247-0.915*(1.-pow(TA,-1./3.));
      Fcorr = (1.-0.15*exp(-Emev))/(1.-0.0007*At); // high-energy correction
      rN = (At-Zt>1.5) ? log(TA-Zt) : 1.;
      s0 = Pi*10.*pow(1.36,2.)*Fcorr*rN*(1.+pow(TA,1./3.)-b0*(1.-pow(TA,-1./3.)));
      p1 = 8.-8./At-0.008*At;
      p2 = 2.*(1.17-2.7/At-0.0014*At);
      p3 = 0.8+18./At-0.002*At;
      p4 = 5.6-0.016*At;
      p5 = 1.37*(1.+1./At);
      x = log10(Emev);
      f1 = 1./(1.+exp(-p1*(x+p2))); // low-energy return to zero
      f2 = 1. +p3 *( 1. -1./(1.+exp(-p4*(x+p5))) ); // low-energy rise
      *PA_inel = f1*f2*s0;
    }
    break;
            
  case 2: // [BP01]
  default:
    if (Em<14.) Em=14.;
    if (Em>1.e6) Em=1.e6;
    *PA_inel = sighad_cc(ISS,PZ,PZ,TA,TZ,Em); // IMOS20020502
  }
  if(Zp*At >= 1) return;
    
  //## AntiProton-Proton ANNIHILATION cross section [TN83]
  if(Ek < 10.) *aPP_ann = 661.*(1.+0.0115/pow(Ek,0.774)-0.948*pow(Ek,0.0151)); // 0.1GeV<Ek<12GeV
  else
    {
      // assuming aPP_ann = aPP_tot -PP_tot (i.e., aPP_elast = PP_elast); (aPP_tot-PP_tot) from [PDG00]
      s2 = 2.*mp*(Ek+2*mp);                   // square of the total CMS energy
      *aPP_ann = 2*35.43/pow(s2,0.560);
    }
    
  //## AntiProton-Proton TOTAL INELASTIC cross section
  aPP_inel = *PP_inel + *aPP_ann;
  if(Ek <= 14.)
    { 
      aPP_inel = 24.7*(1.+0.584/pow(Ek,0.115)+0.856/pow(Ek,0.566));
      if(*aPP_ann > aPP_inel) *aPP_ann = aPP_inel;
    }
    
  //## AntiProton-Proton TOTAL INELASTIC NON-ANNIHILATION cross section
  *aPP_non = aPP_inel - *aPP_ann;
  if(*aPP_non < 0.) *aPP_non = 0.;
    
  //## AntiProton-NUCLEUS cross sections
  if(At > 1)
    {
      //## AntiProton-NUCLEUS TOTAL INELASTIC NON-ANNIHILATION cross section
      *aPA_non = *PA_inel;
        
      //## AntiProton-NUCLEUS ANNIHILATION cross section on 12C-nucleus [mb] (0.4<Pp<300) [MO97]
      A = At;
      Z = Zt;                               // Z = 0.59*pow(A,.927);  for Z > 2 nuclei
      if(At == 4) { Z = 2.; A = 3.30; }     // Modified to agree with HE p-He cs / imos
      *aPA_ann = pow(A,2./3.)               // Scaling to other nuclei
        //         *(48.2 +19./pow(Ek-0.02,0.55) 
        *(48.2 +19./pow(Ek,0.55)           // modified to agree w. He@<100 MeV / imos
          +(0.1-0.18/pow(Ek,1.2))*Z +0.0012/pow(Ek,1.5)*Z*Z)  - *aPA_non;
      if(*aPA_ann < 0.) *aPA_ann = 0.;
      if(*aPA_ann < *aPP_ann) *aPA_ann = *aPP_ann;
      if(At == 4 && Ek > 5.)  *aPA_ann = *aPP_ann;
      /*
      //## my fit to AntiProton-NUCLEUS total cross section on 12C-nucleus [mb] (0.4<Pp<300)
      double Pp =sqrt(pow(Ek+mp,2)-mp*mp); // GeV,kin. momentum per nucleon
      *aPA_ann = (Pp > 40.) ? 236.*(1.+6.9e-2*exp(-Pp/100.)) :
      209.7*pow(.542/Pp,.565)+29.6*cos(log10(Pp/9.29)*5.11)+257.9;
      *aPA_ann *= pow(At/12.,2./3.);                             // scaling to other nuclei
      */
    }
  return;
}

/*
  c***********************************************************************
  c             ###   I.Moskalenko   ###    version of 22.06.1998     ###
  c Antiproton (+antineitron) production spectrum vs. momentum [barn c/GeV] for 
  c pp-, pA-, Ap-, and AA-collisions (Pp1, Pap1 fixed) per 1 target nucleus/cm^3. 
  c Refs: Moskalenko I.V. et al. 2002, ApJ 565, 280
  c (Tan & Ng 1983, J.Phys.G:Nucl.Phys.9,227; ibid.,p.1289;
  c Letaw et al.1983,ApJS,51,271; Gaisser & Schaeffer 1992,ApJ,394,174;
  c Westfall et al.1979,PRC,19,1309)
  c
  c Pap1 [GeV/c] - secondary anti-p momentum; Pp1 [GeV/c] -beam momentum/nucleus
  c NA1 & NA2 are the atomic numbers of beam and target nuclei, correspondingly
  c (NA1=NA2=1 for pp-collisions)
  c***********************************************************************
*/

double TGalpropXSec::antiproton_cc1(gsl_integration_workspace* w, size_t limit, int key, double Pap, double Pp1, int NZ1, int NA1, int NZ2, int NA2) {
    
  double Pp = Pp1/double(NA1);
    
  double s2 = 2.*mp*( sqrt( pow(Pp,2.) + pow(mp,2.) ) + mp );
  double s1 = sqrt(s2); // Center of Mass total energy
    
  if (s1 <= 4.*mp) 
    return 0.;
    
  double MX = 3.*mp;
    
  double gamma_CM = s1/(2.*mp);  // Center of Mass Lorentz factor
  double betagamma_CM = sqrt( pow(gamma_CM, 2.) - 1.  ); // Center of Mass beta*gamma
  double Eap = sqrt ( pow(Pap,2.) + mp*mp );
    
  if (NA2 > 1.) {
    Eap = Eap + 0.06;
    Pap = sqrt( Eap*Eap - mp*mp);
  }
    
  double Eap_MAX = (s2 - MX*MX + mp*mp)/(2.*s1);
  double Ek = sqrt(Pp*Pp + mp*mp)-mp; 
    
  double result = 0.;
  double cosX_MAX = -1.;
    
  if (betagamma_CM*Pap > 0)
    cosX_MAX = (gamma_CM*Eap - Eap_MAX)/(betagamma_CM*Pap); 
    
  if (cosX_MAX < -1.)
    cosX_MAX = -1.;
    
  double param[2] = {Pap, Pp};
    
  if (cosX_MAX <= 1.0) {
        
    //    size_t limit = 1000;
        
    double AI = 0;
    double error = 0;
    gsl_function F;
    F.function = &(tan_ng);
    F.params = param;
    gsl_integration_qag(&F, 1.0, cosX_MAX, 0, 1e-4, limit, 6, w, &AI, &error);
    // dF/dEap - production spectrum vs. energy; factor 2 accounts for antineutrons
    result = -AI * 2.0 * Pap  * (2. * Pi) * 1.e-3 * Pap / Eap  ;   // conversion from mbarn c^3/GeV^2 to mbarn/GeV
    //   gsl_integration_workspace_free(w);
        
  }
  if (NA1 * NA2 == 1 || Ek <= 0.3) return result;
    
  // end of p p -> ap X computation
    
  double PP_inel = 0.;
  double PA_inel = 0.;
  double aPP_non = 0.;
  double aPA_non = 0.;
  double aPP_ann = 0.;
  double apA_ann = 0.;
  nucleon_cs(key, Ek, 1, NZ2, NA2, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &apA_ann);
  double CS2 = PP_inel;
  if (NA2 > 1) CS2 = PA_inel;
  double CS1 = PP_inel;
  if (NA1 > 1) {
    nucleon_cs(key, Ek, 1, NZ1, NA1, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &apA_ann);
    CS1 = PA_inel;
  }
    
  //  multiplicity of anti-p in pA-,Ap-,AA-coll.; Gaisser&Schaeffer 1992,ApJ,394,174
  double correction = 1.2 * (double(NA1)*CS2 + double(NA2)*CS1) / (2.*PP_inel);
  return result*correction;
}

/*
  c***********************************************************************
  c The invar. cross section of the inclusive anti-proton production 
  c [mbarn c^3/GeV^2] from p+p->(anti-p)+X reaction (for Pp,(Pap,COSX) fixed);
  c (Tan & Ng 1983,J.Phys.G:Nucl.Phys.9,227; ibid., p.1289)
  c COSX = cos(theta) is the LS polar angle.
  c***********************************************************************
*/
double tan_ng(double cosX, void* param) {
    
  // This routine computes the inclusive antiproton production cross section
  // result is E d^3sigma/dp^3 expressed in mbarn c^3/GeV^2
  // antiprotons come from the interaction between Cosmic Rays protons and
  // nuclei and interstellar gas
    
  double* param1 = (double*)param;
  double Pap = param1[0];
  double Pp = param1[1];
    
    
  double s2 = 2.*mp*( sqrt(Pp*Pp+mp*mp) + mp );
  double s1 = sqrt(s2); // Center of Mass total energy
    
  if (s1 < 4.*mp) // if CM energy < 4 GeV, no antiprotons are produced (reaction threshold)
    return 0;
    
  double gamma_CM = s1/(2.*mp);
  double betagamma_CM = sqrt( pow(gamma_CM, 2.) - 1.  ); // Center of Mass beta*gamma
  double Eap = sqrt(Pap*Pap + mp*mp);
  double sinX = 0.;
  if (1.-pow(cosX,2) > 0)
    sinX = sqrt(1.-pow(cosX,2));
  double pt = Pap*sinX; // anti-proton transverse momentum
  double MX = 3.*mp;
    
  double Xr = 2.*s1*( gamma_CM*Eap - betagamma_CM*Pap*cosX )/( s2 - MX*MX + mp*mp  )   ;  
    
  if (Xr > 1)
    return 0;
    
  //Cross section calculation. See Tan&Ng 1983
    
  double A = 0.465 * exp(-0.037*Xr) + 2.31 * exp(0.014*Xr); // c/GeV
  double B = 0.0302 * exp(-3.19*(Xr + 0.399))*pow(Xr+0.399, 8.39); 
  //(c/GeV)^2
  double f = (3.15-1.05e-4)*pow((1.-Xr), 7.9);
  if (Xr <= 0.5)
    f += 1.05e-4*exp(-10.1*Xr);
  double result = f*exp(-A*pt + B*pt*pt);
    
  if (s1 < 10.) {
        
    double Xt = ((gamma_CM*Eap-betagamma_CM*Pap*cosX)-mp)/( (s2 - MX*MX +
							     mp*mp) / (2.*s1)-mp); 
    double a = 0.306*exp(-0.12*Xt);
    double b = 0.0552*exp(2.72*Xt);
    double c = 0.758 - 0.68*Xt + 1.54*pow(Xt,2.);
    double d = 0.594*exp(2.87*Xt);
    double Q1 = s1 - 4.*mp;
    double correction = (1. - exp(-exp(c*Q1-d)*(1.-exp(-a*pow(Q1,b))) ));
    result /= correction;
        
  }
    
  return result;
    
}
/*
  void TGalpropXSec::Kcapture_cs(double Ek, int Zp, int Zt, double *attach, double *strip)
  {
  double gam = 1.+Ek*1.e-3/amu, sT =1.e27 * 8./3.*Pi*pow(Rele,2), //AWS20010829
  beta,Mbet,Nbet,a,fcor;
 
  beta = sqrt(1.-pow(gam,-2));
  Mbet = 4./3.+gam*(gam-2.)/(gam+1.)*(1.-log((1.+beta)/(1.-beta))/2./beta/pow(gam,2));
  Nbet = pow(beta,-3) *((-4.*gam +34. -63./gam +25./pow(gam,2) +8./pow(gam,3))/15.
  -(gam-2.)*(gam-1.)/2./beta/pow(gam,3) *log((1.+beta)/(1.-beta)));
  a = Zp*ALPHAf;
  fcor = pow(a,2.*(sqrt(1.-a*a)-1.)) *exp(-2.*a/beta*acos(a)) *(1.+Pi*a*Nbet/Mbet);
  //    factor 1.202 accounts for contribution of states higher than 1s
  *attach = 1.202 *fcor *3./2.*sT*pow(1.*Zp,5)*Zt*pow(ALPHAf,4) *beta*gam *pow(gam-1.,-3)*Mbet;
  *strip = 3./2.*sT*pow(Zp*beta*ALPHAf,-2) 
  *0.285*(2.*log(2./0.219*beta*gam/Zp/ALPHAf)-pow(beta,2))*Zt*(Zt+1.);
  return;
  }
*/

double TSpallationNetwork::spec_int(double ep, double es, int id, int reac){
  /*-----------------------------------------------------------------------------
    interpolation of lab. energy spectra of photons and antiprotons (+antineutrons)
    input:
    ep  - incident energy per nucleon (GeV);
    es  - secondary particle energy (GeV);
    id  - secondary particle type (0 - photon, 1 - antiproton)
    reac - production process (1 - p+p; 2 - p+He; 3 - He+p)
    reac=0: photon spectra for p+p collisions based on the non-diffractive part
    of the Kamae parametrization for ep<ethr and present tables for ep>ethr 
    -----------------------------------------------------------------------------*/
    
  /*implicit double precision (a-h,o-z)
    double precision kamae_nd
    integer reac
    common/data/data_gam(ngamdat,nebin,nreac)
    *,data_ap(napdat,nebin,nreac)
    common/ethr/ethr
    data xlim/.55d0,.45d0/
    data ilim/28,23/
  */  
  double spec_int_temp;
    
  double xlim[3] = {-999,.55e0,.45e0};
  int ilim[3] = {-999,28,23};
    
  if (es > ep || es < 0.e0){
    cerr<<"respect kinematic limits!"<<endl;
    exit(1);
  }
    
  if (id*(id-1) != 0){
    cerr<<"no tables for this particle!"<<endl;
    exit(2);
  }
    
  if (reac == 0 && ep < ethr){ //for reac=0 and ep<ethr use ND Kamae param.
    if (id != 0){
      cerr<<"reac=0 - for p+p->gamma only!"<<endl;
      exit(3);
    }
    //spec_int_temp=kamae_nd(ep,es);
  }
    
  if (reac*(reac-1)*(reac-2)*(reac-3) != 0){
    cerr<<"no tables for this reaction!"<<endl;
    exit(4);
  }
    
  if (ep < 10.e0 && reac>1){
    cerr<<"results are questionable below 10 GeV lab.!"<<endl;
    exit(5);
  }
    
  if (ep > 1.e8){
    cerr<<"extrapolation beyond the tabulated range!"<<endl;
    exit(6);
  }
    
  double yy, x, xx, dspec;
  int ix, jy, ir;
    
  double wi[4],
    wj[4];
    
  if (ep < 1.e3)
    yy=log10(ep)*4.e0-3.e0;
  else
    yy=log10(ep)*2.e0+3.e0;
    
  //jy=max_int(1,(int)yy);
  jy= (1>(int)yy) ? 1 : (int)yy;
  //jy=min_int(jy,nebin-2);
  jy= (jy<nebin-2) ? jy : nebin-2;
  if(jy == 8) jy=7;
  wj[2]=yy-(double)jy;
  wj[3]=wj[2]*(wj[2]-1.e0)*.5e0;
  wj[1]=1.e0-wj[2]+wj[3];
  wj[2]=wj[2]-2.e0*wj[3];
    
  x=es/ep;
    
  if (x < xlim[id+1])
    xx=50.e0*x+.5e0;
  else
    xx=ilim[id+1]+(x-xlim[id+1])*10.e0;
    
  //ix=max_int(1,(int)xx);
  ix= (1>(int)xx) ? 1 : (int)xx;
  if (ix == ilim[id+1]-1)
    ix=ix-1;
    
  spec_int_temp=0.e0;
  //ir=max_int(1,reac);
  ir= (1>reac) ? 1 : reac;
    
  if(id == 0){
    /*ix=min_int(ix,ngamdat-3);
      if( data_gam[ix][jy][ir] * data_gam[ix][jy+1][ir] * data_gam[ix][jy+2][ir] == 0.e0)
      ix=ix-1;*/
  }
  else if(id == 1){
    ix=min(ix,napdat-3);
    if( data_ap[ix][jy][ir]*data_ap[ix][jy+1][ir]*data_ap[ix][jy+2][ir] == 0.e0){
      if(data_ap[ix+1][jy][ir] * data_ap[ix+1][jy+1][ir] * data_ap[ix+1][jy+2][ir] == 0.e0)
	return -1;
      ix=ix+1;
    }
    if(data_ap[ix+2][jy][ir]*data_ap[ix+2][jy+1][ir]*data_ap[ix+2][jy+2][ir] == 0.e0)
      ix=ix-1;
    if(data_ap[ix+2][jy][ir] == 0.e0){
      if(data_ap[ix+1][jy][ir] == 0.e0)
	return -1;
      ix=ix-1;
    }
  }
    
  wi[2]=xx-(double)ix;
  wi[3]=wi[2]*(wi[2]-1.e0)*.5e0;
  wi[1]=1.e0-wi[2]+wi[3];
  wi[2]=wi[2]-2.e0*wi[3];
    
  if (id == 0){
    /*for ( int i=1;i<=3;i++)
      for ( int j=1;j<=3;j++)
      spec_int_temp=spec_int_temp+log(data_gam[ix+i-1][jy+j-1][ir])*wi[i]*wj[j];*/
  }
  else if (id == 1){
    for ( int i=1;i<=3;i++)
      for ( int j=1;j<=3;j++)
	{
	  if (ix == ilim[id+1]-1 && i == 1)
	    dspec=log(data_ap[ix+i-5][jy+j-1][ir]);
	  else
	    dspec=log(data_ap[ix+i-1][jy+j-1][ir]);
                
	  spec_int_temp=spec_int_temp+dspec*wi[i]*wj[j]; 
	}
  }
    
  return exp(spec_int_temp);
}

void TSpallationNetwork::spec_ini(){
  /*-----------------------------------------------------------------------------
    read tables
    -----------------------------------------------------------------------------*/
    
  /*parameter(ngamdat=32,napdat=28,nebin=19,nreac=3)
    common/data/data_gam(ngamdat,nebin,nreac)
    *,data_ap(napdat,nebin,nreac)
    common/ethr/ethr
     
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) '!                                                   !'
    write(*,*) '! ppfrag -- for more infor see                      !'
    write(*,*) '! http://sourceforge.net/projects/ppfrag            !'
    write(*,*) '!                                                   !'
    write(*,*) '! if you use data from reac=1,2,3, please cite      !'
    write(*,*) '! M.Kachelriess, S.Ostapchenko                      !'
    write(*,*) '! arxiv[1206.xxxx], submitted to PRD                !'
    write(*,*) '!                                                   !'
    write(*,*) '! if you use reac=0, please cite additionally       !'
    write(*,*) '! T. Kamae et al., ApJ 647 (2006) 692.              !'
    write(*,*) '!                                                   !'
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     
    open(1,file='gamfrag.dat',status='old')
    read(1,*)data_gam
    close (1)
    open(1,file='apfrag.dat',status='old')
    read(1,*)data_ap
    close (1)
    write(*,*)'cross section tables read'
    return
    end
  */
    
  double temp;
    
  ifstream file_to_read("data/apfrag.dat");
    
  if (file_to_read.good())
    {
      for ( int k=1;k<=nreac;k++)
	for ( int j=1;j<=nebin;j++)
	  for ( int i=1;i<=napdat;i++){
	    file_to_read >> temp;
	    data_ap[i][j][k] = temp;
	  }
      file_to_read.close();
    }
  else{
    cout<<"File apfrag.dat not properly read!"<<endl;
    exit(11);
  }
    
  return;
}

//*********************************************************************************************************************************************
//   Fluka model for nucleus-nucleus xsec
//*********************************************************************************************************************************************


FlukaXSec::FlukaXSec(TGrid* co) : TXSecBase() {

  initialized = false;
    
  /*ifstream file_to_read(FlukaNucleiTablefile.c_str(), ios::in);

    int parent_uid;
    int daughter_uid;
    int gas_uid;
    double current_energy = 0.;
    double xsec_;

    while (file_to_read.good()) {
            
    file_to_read >> parent_uid >> daughter_uid >> gas_uid >> current_energy >> xsec_;
    if (parent_uid == 0 && daughter_uid == 0) {
    energy.push_back(current_energy);
    }
    else { 
    pair<int,int> couple1(parent_uid,daughter_uid);
    pair<pair<int,int>,int> couple2(couple1,gas_uid);
    xsec_extended[couple2].push_back(xsec_);  
    if (gas_uid == 1001) xsec[couple1].push_back(xsec_);            
    }
    }

    file_to_read.close();
  */

  /*ifstream file_to_read(FlukaLightNucleiTablefile.c_str(), ios::in);

    int parent_uid;
    int daughter_uid;
    int gas_uid;
    double parent_energy = 0.;
    double daughter_energy = 0.;
    double xsec_;

    while (file_to_read.good()) {
            
    file_to_read >> parent_uid >> daughter_uid >> gas_uid >> current_energy >> daughter_energy >> xsec_;
    if (parent_uid == 0 && daughter_uid == 0 && daughter_energy == 0) {
    energy.push_back(current_energy);
    }
    else { 
    pair<int,int> couple1(parent_uid,daughter_uid);
    pair<pair<int,int>,int> couple2(couple1,gas_uid);
    xsec_extended[couple2].push_back(xsec_);  
    if (gas_uid == 1001) xsec[couple1].push_back(xsec_);            
    }
    }

    file_to_read.close();*/

  /*
    ifstream file_to_read(FlukaData.c_str(), ios::in);
    const int max_num_of_char_in_a_line = 512;
    
    int channel_temp;
    double energy_temp, xsec_double;
    
    while (file_to_read.good()) {
 
    file_to_read >> channel_temp;
        
    if (channel_temp == 0) {
    for (int i = 0; i < 15; i++) {
    file_to_read >> energy_temp;
    energy.push_back(energy_temp/1000.0);
    }
    }
    else {
            
    int uid1 = -1;   // Parent
    int uid2 = -1;   // Daughter
    Fluka::convert(channel_temp, uid1, uid2);
            
    pair<int,int> couple = pair<int,int>(uid1,uid2);
    map< pair<int,int>, vector<double> >::iterator it = xsec.find(couple);
            
    if (it == xsec.end()) {
    xsec[couple] = vector<double>(energy.size(), 0.0);
    for (int i = 0; i < 15; i++) {
    file_to_read >> xsec_double;
    xsec[couple][i] = xsec_double;
    }
    }
    }
        
    file_to_read.ignore(max_num_of_char_in_a_line, '\n');
        
    }
    
    file_to_read.close();*/
    

  return;
}


void FlukaXSec::InitTableInelastic() {

  // reads the Inelastic cross section table and fills the structure 
  // map< pair<int, int>, vector<double> > xsec_inelastic;

  if (initialized) return;

  cout << "=====>>>>>>> Entering FlukaXsec InitTableInelastic routine ..." <<  endl;

  ifstream infile(Fluka_ap_inelastic_H_datafile.c_str(), ios::in); // Read databasis
  double en, xsec_H;
  while (infile >> en >> xsec_H) {
    energy_vec_ap_H.push_back(en);
    inelastic_xsec_ap_H_vec.push_back(xsec_H);
  }
  infile.close();

  ifstream infile2(Fluka_ap_inelastic_He_datafile.c_str(), ios::in); // Read databasis
  double xsec_He;
  while (infile2 >> en >> xsec_He) {
    energy_vec_ap_He.push_back(en);
    inelastic_xsec_ap_He_vec.push_back(xsec_He);
  }
  infile2.close();

  initialized = true;

  /*ifstream infile(Fluka_inelastic_datafile.c_str(), ios::in); // Read databasis

    int uid, uid2;
    vector<double> energy_vec;
    double energy_, xsec_value;
  
    // reads the datafile
    while (infile >> uid  >> energy_ >> xsec_value) {

    //cout << uid0 << " " << uid1 << " " << energy << " " << xsec_value << endl;

    if (uid == 0) 
    energy_vec.push_back(energy_);
    else {	
    totalxsec[uid].push_back(xsec_value);
    }
    }*/

}

double FlukaXSec::GetTotalApHXSec(double en) {

  int N = energy_vec_ap_H.size(); 
  double e_vec[N];
  double xsec_vec[N];
  for (int i=0; i<N; i++){
    e_vec[i] = energy_vec_ap_H[i];
    xsec_vec[i] = inelastic_xsec_ap_H_vec[i];
  }
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
  gsl_spline_init (spline, e_vec, xsec_vec, N);
  double result = gsl_spline_eval (spline, en, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  return result;
}

double FlukaXSec::GetTotalApHeXSec(double en) {

  int N = energy_vec_ap_He.size(); 
  double e_vec[N];
  double xsec_vec[N];
  for (int i=0; i<N; i++){
    e_vec[i] = energy_vec_ap_He[i];
    xsec_vec[i] = inelastic_xsec_ap_He_vec[i];
  }
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
  gsl_spline_init (spline, e_vec, xsec_vec, N);
  double result = gsl_spline_eval (spline, en, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return result;
}

double FlukaXSec::GetXSec(int input, double en) { //modified

  int Z = -1000;
  int A = -1000;
  Utility::id_nuc(input, A, Z);
  /*
    if (uid == 28058) {
    for (unsigned int i = 0; i < energy.size(); i++) outf << energy[i] << "\t";
    outf << endl;
    }
    outf << uid << endl;
  */
    
  if (A == 0) return 0.0;
    
  const double factor = 1.e-27*Clight;
    
  int N = energy_vec.size();
  if (en > energy_vec.back()) return totalxsec[input].back();
    
  gsl_vector_view x = gsl_vector_view_array (&energy_vec[0], N);
  gsl_vector_view y = gsl_vector_view_array (&(totalxsec[input][0]), N);
    
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
    
  gsl_spline_init (spline, (&x.vector)->data, (&y.vector)->data, N);
    
  double result = gsl_spline_eval (spline, en, acc);
    
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  /*int l = 0;
    while (energy_vec[l] < en) l++;
    double e1 = (energy_vec[l]-en)/(energy_vec[l]-energy_vec[l-1]);
    double e2 = (en-energy_vec[l-1])/(energy_vec[l]-energy_vec[l-1]);
  */

  //double result1 = totalxsec[input][l-1]*e1 + totalxsec[input][l]*e2;
    
  //
  return max(result,0.0);  
  
  // if (K_electron >= 0) {
            
  //not yet implemented
  // }
    
}



//not used now
double FlukaXSec::GetXSec(pair<int,int> input, double en) {
    
  double result = 0.0;

  map< pair<int,int>, vector<double> >::iterator it = xsec.find(input);
    
  if (it == xsec.end()) {
    // cerr << "Channel not found" << endl;
    return -1;
  }
    
  int N = energy.size();
  if (en > energy.back()) return (*it).second.back();
    
  gsl_vector_view x = gsl_vector_view_array (&energy[0], N);
  gsl_vector_view y = gsl_vector_view_array (&((*it).second)[0], N);
    
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
    
  gsl_spline_init (spline, (&x.vector)->data, (&y.vector)->data, N);
    
  result = gsl_spline_eval (spline, en, acc);
    
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
    
  /*int l = 0;
    while (energy[l] < en) l++;
    double e1 = (energy[l]-en)/(energy[l]-energy[l-1]);
    double e2 = (en-energy[l-1])/(energy[l]-energy[l-1]);
    
    double result1 = (*it).second[l-1]*e1 + (*it).second[l]*e2;*/
  
  return max(result,0.0);


}

//not used now
double FlukaXSec::GetXSecExtended(pair<int,int> input, int gas_uid, double en) {
    
  double result = 0.0;
  
  pair<pair<int,int>, int> combined_input(input,gas_uid);

  map<pair<pair<int,int>, int>, vector<double> > ::iterator it = xsec_extended.find(combined_input);
    
  if (it == xsec_extended.end()) {
    // cerr << "Channel not found" << endl;
    return -1;
  }
    
  int N = energy.size();
  if (en > energy.back()) return (*it).second.back();
    
  gsl_vector_view x = gsl_vector_view_array (&energy[0], N);
  gsl_vector_view y = gsl_vector_view_array (&((*it).second)[0], N);
    
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
    
  gsl_spline_init (spline, (&x.vector)->data, (&y.vector)->data, N);
    
  result = gsl_spline_eval (spline, en, acc);
    
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
    
  /*int l = 0;
    while (energy[l] < en) l++;
    double e1 = (energy[l]-en)/(energy[l]-energy[l-1]);
    double e2 = (en-energy[l-1])/(energy[l]-energy[l-1]);
    
    double result1 = (*it).second[l-1]*e1 + (*it).second[l]*e2;
  */
  
  return max(result,0.0);


}

//not used now
vector<double> FlukaXSec::GetTotalXSec_vec(int uid /**< Nucleus uid. */, int gas_uid /**< Gas uid. */) {

  // this method gets the total inelastic Xsec vector for a given nucleus on a given gas target (H, He, C, ...) from the structure "xsec_inelastic"

  pair<int,int> bullet_gas(uid, gas_uid);
  map< pair<int,int>, vector<double> >::iterator it = totalxsec_extended.find(bullet_gas);
  if (it != totalxsec_extended.end()) return (*it).second;
  else {
    cerr << "No Inelastic cross section found for nucleus with uid = " << uid << endl;
    return vector<double>();
  }

}




// completare...
