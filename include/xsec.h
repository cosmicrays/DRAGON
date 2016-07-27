/*
 * @file xsec.h
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief All the classes related to the cross sections are defined.
 */


#ifndef _XSEC_H
#define _XSEC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <math.h>

#include <gsl/gsl_integration.h>

#include "constants.h"

class TGrid;
class TNucleiList;
class Input;


extern "C" void yieldx_(int*,int*,int*,int*,float*,float*);
/**
 * @fn extern "C" void yieldx_(int*,int*,int*,int*,float*,float*)
 * @brief Silberberg & Tsao isotopic production cross section
 */
extern "C" double wsigma_(int*,int*,int*,int*,double*);
/**
 * @fn extern "C" double wsigma_(int*,int*,int*,int*,double*)
 * @brief Wrapper to Webber's (1993) code, from Galprop
 */
extern "C" void set_sigma_(int*); /**< Initialization of Webber's code. */
extern "C" void sigtap2_(int*); /**< initialization of the Barashenkov & Polanski cross section code. */
extern "C" double sighad_(int*,double*,double*,double*,double*,double*); /**< Barashenkov & Polanski pA total cross section. */

using namespace std;

/**
 * @class TXSecBase
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Abstract class for cross sections.
 *
 */
class TXSecBase {
   
   public :
   TXSecBase() { }
   virtual ~TXSecBase() { }
   virtual void Print(TGrid*)=0;
   virtual vector<double> GetXSec(int /**< Charge of primary. */, int /**< Mass of primary. */, int /**< Charge of secondary. */, int /**< Mass of secondary. */)=0;
   virtual vector<double> GetXSec(int /**< uid of nucleus. */)=0;
   virtual void InitTableInelastic() { throw -1; }
   //modified
   virtual vector<double> GetXSec(int, int)=0;
   virtual void nucleon_cs(int, double /**< kinetic energy per nucleon of the beam momentum particles (GeV) */, int /**< =+/-1 is the charge of beam and At is the atomic number of target nuclei */, int /**< =1 for proton, = -1 for antiproton collisions. */, int /**< =1 for proton collisions. */, double* /**< [mbarn] is the proton-proton inelastic cross sect. */, double* /**< [mbarn] is the proton-nucleus inelastic cross sect. (p+A2) */, double* /**< [mbarn] is the antiproton-proton inelastic cross sect.(nonannihil.) */, double* /**< [mbarn] put =PA_inel, is the antiproton-nucleus inelastic cross sect.(nonan.) */, double* /**< [mbarn] is the antiproton-proton annihilation cross sect. */, double* /**< [mbarn] is the antiproton-nucleus annihilation cross sect. */)=0; // { /*...*/ }
   
   //modified
   virtual void Kcapture_cs(double, int, int, double*, double*)=0;// { /* */ } /* Cross section of electron attachment/stripping process */
   
   virtual double GetXSec(pair<int,int> /**< Pair of parent-daughter nuclei uid. */, double /**< Energy. */)=0;// { /*...*/ }
   virtual double GetTotalXSec(int /**< Nucleus uid. */, double /**< Energy. */)=0;// { /*...*/ }
   virtual double GetTotalApHXSec(double)=0;
   virtual double GetTotalApHeXSec(double)=0;
   virtual vector<double> GetTotalXSec_vec(int /**< Nucleus uid. */)=0; 
   virtual map<pair<int,int>, vector<double> > GetTotalXSec_for_each_gas_type()=0;
   virtual bool IsPresent(int uid /**< Nucleus uid. */)=0;
   virtual bool IsPresent(pair<int,int> uid /**< Nucleus parent-daughter uid. */)=0;
   inline double GetHefactor() { return 1.26/He_abundance; }
   /**
    * @fn inline double GetHefactor()
    * @brief Factor to account for ISM Helium.
    * @return 1.26/He_abundance
    *
    * Still experimental and to be understood with Webber.
    */
   
   virtual double antiproton_cc1(gsl_integration_workspace*,size_t,int,double,double,int,int,int,int)=0;// { /*...*/ }
   
};

/**
 * @class TSpallationNetwork
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class containing the description of the network of spallation used to compute secondary CR source terms.
 *
 * Two models can be used up to now. One is the Galprop model for the spallation cross section. The second model uses Webber (2003) cross sections. For electron and positron production also Kamae (2007) model can be used. \sa cparammodel.h
 */
class TSpallationNetwork {
   
public:
   
   TSpallationNetwork(TGrid* co /**< Kinematics. */, Input*, vector<TXSecBase*> xsecmodel /**< Model of the cross section */, vector<int>& nuclei /**< List of propagating nuclei. */);
   /**
    * @fn TSpallationNetwork(TGrid* co, TXSecBase* xsecmodel, vector<int>& nuclei)
    * @brief Construct a TSpallationNetwork object.
    */
   
   ~TSpallationNetwork() {
	spall.clear();
	spall_apel.clear(); 
	energy.clear();
    }
   
   vector<double> GetXSec(int /**< Parent nucleus uid */,int /**< Daughter nucleus uid. */);
   /**
    * @fn vector<double> GetXSec(int,int)
    * @brief Obtain the spallation cross section.
    */
   double GetXSec(int /**< Parent nucleus uid */, int /**< Daughter nucleus uid. */, double /** Energy. */);
   /**
    * @fn double GetXSec(int, int, double)
    * @brief Obtain spallation cross section at given energy.
    */
   double GetXSec(int /**< Parent nucleus uid */, int /**< Daughter nucleus uid. */, int /** Energy index. */);
   /**
    * @fn double GetXSec(int, int, int)
    * @brief Obtain spallation cross section at given energy.
    */
   vector<double> GetXSecApEl(int /**< Parent nucleus uid */, int /**< Daughter particle uid. */, int /** Energy index of parent particle. */);
   /**
    * @fn vector<double> GetXSecApEl(int, int, int)
    * @brief Obtain energy spectrum of secondary antiprotons or electrons/positrons, given the energy of parent particle.
    */
   
   //TPP
   vector<double> GetXSecTPP(vector<double> /** ISRF frequency vector*/);
   /**
    * @fn vector<double> GetXSecApEl(vector<double>)
    * @brief Obtain the vetor of production cross section of TPP positrons for each energy of parent and daughter particle, given the ISRF frequency vector
    */
   
   void spec_ini();
   double spec_int(double ep, double es, int id, int reac);
   
private:
   map<pair<int,int>, vector<double> > spall;
   /**
    * @var map<pair<int,int>, vector<double> > spall
    * @brief Description of the spallation cross section.
    *
    * The spallation cross section is associated (through std::map container) to the pair of parent-daughter nucleus (through std::pair). The first element of the std::map is an std::pair containing as a first element the uid of the parent nucleus, and as second element the uid of the daughter nucleus.
    */
   map<pair<int,int>, vector< vector<double> > > spall_spectrum;
   /**
    * @var map<pair<int,int>, vector< vector<double> > > spall_apel
    * @brief Same as for the spallation cross section. But in this case we store the spectrum of secondary particles, given the energy of the parent particle. \sa spall
    */
   vector<double> energy; /**< Energy array. */
   vector<double> energy_daughter; /**< Energy array. */
   map<pair<int,int>, vector< vector<double> > > spall_apel;
   /**
    * @var map<pair<int,int>, vector< vector<double> > > spall_apel
    * @brief Same as for the spallation cross section. But in this case we store the spectrum of secondary antiprotons and/or electrons/positrons, given the energy of the parent particle. \sa spall
    */
   //void InitDataTables();
   /**
    * @fn void InitDataTables()
    * @brief Initialize tables for the production of secondary electrons and positrons.
    */
   int index_matrix(int i /**< Energy index of daughter particle. */, int j /**< Energy index of parent particle. */) { return i*Nprotons+j; }
   /**
    * @fn int index_matrix(int i, int j)
    * @brief Convert matrix to linearized index, for data tables. \sa InitDataTables.
    * @return i*801+j
    */
   void InitXSecKamae(double);
   void InitXSecGalprop(double);
   void InitXSecPohl(double);
   void InitXSecFluka(double, int);
   void InitDataTablesPohl();
   void InitDataTablesGalprop();
   void InitDataTablesFluka(int);
   int Nelectrons;
   int Nprotons;
   int Nap;
   vector<double> ProdXsec;
   vector< vector<double> > ElppPohl;
   vector< vector<double> > ElHepPohl;
   vector< vector<double> > PosppPohl;
   vector< vector<double> > PosHepPohl;
   vector<double> Matrix_El_pp;    /**< Table for the electron spectrum from pp interactions. */
   vector<double> Matrix_El_pHe;   /**< Table for the electron spectrum from pHe interactions. */
   vector<double> Matrix_El_Hep;   /**< Table for the electron spectrum from Hep interactions. */
   vector<double> Matrix_El_HeHe;  /**< Table for the electron spectrum from HeHe interactions. */
   vector<double> Matrix_Pos_pp;   /**< Table for the positron spectrum from pp interactions. */
   vector<double> Matrix_Pos_pHe;  /**< Table for the positron spectrum from pHe interactions. */
   vector<double> Matrix_Pos_Hep;  /**< Table for the positron spectrum from Hep interactions. */
   vector<double> Matrix_Pos_HeHe; /**< Table for the positron spectrum from HeHe interactions. */
   vector<double> Matrix_Ap_pp;    /**< Table for the ap spectrum from pp interactions. */
   vector<double> Matrix_Ap_pHe;   /**< Table for the ap spectrum from pHe interactions. */
   vector<double> Matrix_Ap_Hep;   /**< Table for the ap spectrum from Hep interactions. */
   vector<double> Matrix_Ap_HeHe;  /**< Table for the ap spectrum from HeHe interactions. */
   vector<double> Matrix_p_pp;    /**< Table for the p spectrum from pp interactions. */
   vector<double> Matrix_p_pHe;   /**< Table for the p spectrum from pHe interactions. */
   vector<double> Matrix_p_Hep;   /**< Table for the p spectrum from Hep interactions. */
   vector<double> Matrix_p_HeHe;  /**< Table for the p spectrum from HeHe interactions. */
   vector<double> Matrix_3Ap_app;    /**< Table for the p spectrum from app interactions. */
   vector<double> Matrix_3Ap_apHe;   /**< Table for the p spectrum from apHe interactions. */
   
   static const int napdat = 28;
   static const int nebin = 19;
   static const int nreac = 3;
   static const double ethr;
   double data_ap[napdat+1][nebin+1][nreac+1];
};

/**
 * @class TInelasticCrossSection
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class containing a description of the total inelastic cross section, used to compute the disappearance rate of CRs.
 *
 * Two models can be used up to now. One is the Galprop model for the cross section. The second model uses Webber (2003) cross sections.
 */
class TInelasticCrossSection {
   
public:
   TInelasticCrossSection() {} /**< Default constructor. */
   
   //modified
   TInelasticCrossSection(TGrid* /**< Kinematics. */, Input*, int /**< Nucleus uid. */, int /**< K_electron */,  vector<TXSecBase*>);
   /**
    * @fn TInelasticCrossSection(TGrid*, int)
    * @brief Constructor of total inelastic cross section given the kinematics and the specified particle.
    */
   virtual ~TInelasticCrossSection() { 
	xsec.clear();
	xsec_for_each_gas_type.clear(); 
   } /**< Destructor. */

   inline vector<double>& GetXSec() { return xsec; } /**< Obtain the cross section. */

   inline map<pair<int, int>, vector<double> >& GetXSec_extended() { 
	cout << "returning extended inelastic cross section" << endl;
	return xsec_for_each_gas_type; } /**only for Fluka model */

   inline double GetXSec(int i /**< Energy index. */) { return xsec[i]; } /**< Obtain cross section at given energy. */
   
protected:

   vector<double> xsec; /**< Array of cross section. */

   //not used now
   map<pair<int, int>, vector<double> > xsec_for_each_gas_type; /**< Array of cross section for each nucleus against each gas type, only for Fluka model. */
   
};

/**
 * @class TWebber03
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class describing Webber 03 cross sections.
 */
class TWebber03 : public TXSecBase {
   
public:
   
   TWebber03();
   /**< Construct a TWebber03 object. */
   
   ~TWebber03() { }
   
   double GetXSec(pair<int,int> /**< Pair of parent-daughter nuclei uid. */, double /**< Energy. */);
   /**
    * @fn double GetXSec(pair<int,int>, double)
    * @brief Obtain spallation cross section given the parent-daughter nucleus pair and energy.
    * @warning Uses GSL interpolation.
    */
   double GetTotalXSec(int /**< Nucleus uid. */, double /**< Energy. */);
   /**
    * @fn double GetTotalXSec(int, double)
    * @brief Obtain total inelastic cross section given the nucleus and energy.
    * @warning Uses GSL interpolation.
    */
   virtual vector<double> GetTotalXSec_vec(int /**< Nucleus uid. */) { throw -1; } 
   virtual void InitTableInelastic() { throw -1; }
   virtual double GetTotalApHXSec(double) { throw -1; };
   virtual double GetTotalApHeXSec(double) { throw -1; };
   virtual map<pair<int,int>, vector<double> > GetTotalXSec_for_each_gas_type() {throw -1;}
   virtual vector<double> GetXSec(int /**< uid of nucleus. */) { throw -2; }
   virtual vector<double> GetXSec(int, int) { throw -2; }
   virtual vector<double> GetXSec(int /**< Charge of primary. */, int /**< Mass of primary. */, int /**< Charge of secondary. */, int /**< Mass of secondary. */) { throw -2; }
   virtual void nucleon_cs(int, double /**< kinetic energy per nucleon of the beam momentum particles (GeV) */, int /**< =+/-1 is the charge of beam and At is the atomic number of target nuclei */, int /**< =1 for proton, = -1 for antiproton collisions. */, int /**< =1 for proton collisions. */, double* /**< [mbarn] is the proton-proton inelastic cross sect. */, double* /**< [mbarn] is the proton-nucleus inelastic cross sect. (p+A2) */, double* /**< [mbarn] is the antiproton-proton inelastic cross sect.(nonannihil.) */, double* /**< [mbarn] put =PA_inel, is the antiproton-nucleus inelastic cross sect.(nonan.) */, double* /**< [mbarn] is the antiproton-proton annihilation cross sect. */, double* /**< [mbarn] is the antiproton-nucleus annihilation cross sect. */) { throw -2; }
   virtual void Kcapture_cs(double, int, int, double*, double*) { throw -2; }
   
   virtual bool IsPresent(int uid /**< Nucleus uid. */) {
      map<int, vector<double> >::iterator it = totalxsec.find(uid);
      if (it != totalxsec.end()) return true;
      else return false;
   }
   /**
    * @fn bool IsPresent(int uid)
    * @brief Check whether nucleus is present in Webber03 database. Otherwise use Galprop cross section.
    */
   
   virtual bool IsPresent(pair<int,int> uid /**< Nucleus parent-daughter uid. */) {
      map<pair<int,int>, vector<double> >::iterator it = xsec.find(uid);
      if (it != xsec.end()) return true;
      else return false;
   }
   /**
    * @fn bool IsPresent(pair<int,int> uid)
    * @brief Check whether nucleus is present in Webber03 database. Otherwise use Galprop cross section.
    */
   
   //  inline double GetHefactor() { return 1.26/He_abundance; }
   /**
    * @fn inline double GetHefactor()
    * @brief Factor to account for ISM Helium.
    * @return 1.26/He_abundance
    *
    * Still experimental and to be understood with Webber.
    */
   void Print(TGrid*);
   /**
    * @fn void Print(TGrid* co)
    * @brief Print some specified cross section, only for testing purposes.
    */
   virtual double antiproton_cc1(gsl_integration_workspace*,size_t,int,double,double,int,int,int,int) { return -1; }
   
protected:
   map<pair<int,int>, vector<double> > xsec;
   /**
    * @var map<pair<int,int>, vector<double> > xsec
    * @brief Description of spallation cross section. \sa TSpallationNetwork::spall
    */
   map<int, vector<double> > totalxsec;
   /**
    * @var map<int, vector<double> > totalxsec
    * @brief Description of total inelastic cross section.
    */
   vector<double> energy; /**< Energy array. */
   
   void convert(const int channel /**< Database ID. */, int& uid1 /**< uid of parent nucleus (it will be modified). */, int& uid2 /**< uid of daughter nucleus (it will be modified). */) {
      int A1 = channel/1000000;
      int Z2 = channel%100;
      
      int Z1 = (channel-A1*1000000)/10000;
      uid1 = A1+Z1*1000;
      int A2 = (channel-A1*1000000-Z1*10000)/100;
      uid2 = A2 + Z2*1000;
      return ;
   }
   /**
    * @fn void convert(const int channel, int& uid1, int& uid2)
    * @brief Convert from database to uid for spallation (2 nuclei in one number).
    */
   
   void convert_single(const int channel /**< Database ID. */, int& uid1 /**< uid of nucleus (it will be modified). */) {
      int A1 = channel/100;
      int Z1 = channel%100;
      uid1 = A1+Z1*1000;
      return ;
   }
   /**
    * @fn void convert_single(const int channel, int& uid1)
    * @brief Convert from database to uid for inelastic cross section (1 nucleus in one number).
    */
   
};

/**
 * @class TGalpropXSec
 * @author Luca Maccione
 * @email luca.maccione@desy.de
 * @brief Class describing Galprop cross sections, adapted from Galprop. Only selected parts will be documented.
 *
 * The content of this class has been adapted from Galprop.
 */
class TGalpropXSec : public TXSecBase {
   
public:
   
   TGalpropXSec(TGrid* co);
   /**< Construct an object of this class, given kinematics. */
   
   vector<double> GetXSec(int /**< Charge of primary. */, int /**< Mass of primary. */, int /**< Charge of secondary. */, int /**< Mass of secondary. */);
   /**
    * @fn vector<double> GetXSec(int, int, int, int)
    * @brief Obtain spallation cross section for the specified parent-daughter nuclei.
    */

   virtual vector<double> GetTotalXSec_vec(int /**< Nucleus uid. */) { throw -1; } 
   virtual void InitTableInelastic() { throw -1; }   
   virtual bool IsPresent(int uid /**< Nucleus uid. */) { return true;  }
   virtual bool IsPresent(pair<int,int> uid /**< Nucleus parent-daughter uid. */) { return true; }
   virtual void Print(TGrid*) { }
   virtual vector<double> GetXSec(int /**< uid of nucleus. */) { throw -1; }
   virtual double GetTotalApHXSec(double) { throw -1; };
   virtual double GetTotalApHeXSec(double) { throw -1; };

      virtual double GetXSec(pair<int,int> /**< Pair of parent-daughter nuclei uid. */, double /**< Energy. */) { throw -1; }
      virtual double GetTotalXSec(int /**< Nucleus uid. */, double /**< Energy. */) { throw -1; }
   //modified

   vector<double> GetXSec(int /**< uid of nucleus. */, int /**< K electron */);
   /**
    * @fn vector<double> GetXSec(int)
    * @brief Obtain total inelastic cross section for specified nucleus.
    */
   virtual map<pair<int,int>, vector<double> > GetTotalXSec_for_each_gas_type() {throw -1;}   

   ~TGalpropXSec() {
      for(int j=0; j<N_DATA_FILES; j++) delete [] data_file[j];
   }
   void nucleon_cs(int, double /**< kinetic energy per nucleon of the beam momentum particles (GeV) */, int /**< =+/-1 is the charge of beam and At is the atomic number of target nuclei */, int /**< =1 for proton, = -1 for antiproton collisions. */, int /**< =1 for proton collisions. */, double* /**< [mbarn] is the proton-proton inelastic cross sect. */, double* /**< [mbarn] is the proton-nucleus inelastic cross sect. (p+A2) */, double* /**< [mbarn] is the antiproton-proton inelastic cross sect.(nonannihil.) */, double* /**< [mbarn] put =PA_inel, is the antiproton-nucleus inelastic cross sect.(nonan.) */, double* /**< [mbarn] is the antiproton-proton annihilation cross sect. */, double* /**< [mbarn] is the antiproton-nucleus annihilation cross sect. */);
   /**
    * @fn void nucleon_cs(int, double, int, int, int, double*, double*, double*, double*, double*, double*)
    * @brief parametrization of the pp-, pA-, AA-, (anti_p)p-, and (anti_p)A total inelastic cross sections, as well as (anti_p) annihilation cross sect.
    
    * @param option pA total inelastic cross section - =0 [L83]; - =1 [WA96] for Zt>5 and [BP01] for Zt<=5; - =2 -[BP01];
    *
    * from galprop package * 2001/05/11
    *
    * REFERENCES:
    * - [W79] Westfall et al. 1979, PRC, 19, 1309
    * - [L83] Letaw et al. 1983, ApJS, 51, 271
    * - [TN83] Tan & Ng 1983, J.Phys.G:Nucl.Phys., 9, 227
    * - [PDG00] D.E.Groom et al. 2000, Europ. Phys. J. C15, 1
    * - [MO97] Moiseev & Ormes 1997, Astropart. Phys.6, 379
    * - [WA96] Wellisch H.P., Axen D. 1996, PRC 54, 1329; Wellish 1999, private comm.
    *        (typo corrections & code)
    * - [BP01] V.S.Barashenkov, A.Polanski code used for calc. of the pA total inelastic
    * cross sections
    */
   
   double antiproton_cc1(gsl_integration_workspace*,size_t,int,double,double,int,int,int,int);
   void Kcapture_cs(double Ek, int Zp, int Zt, double *attach, double *strip);
   
protected:
   vector<double> energy; /**< Energy array. */
   vector<double> beta; /**< Velocity array. */
   //modified
   //  void Kcapture_cs(double, int, int, double*, double*); /* Cross section of electron attachment/stripping process */
   void set_sigma_cc(); /**< initialization of Webber's code */
   double wsigma_cc(int,int,int,int,double); /**< Webber's isotopic production cross section */
   double yieldx_cc(int,int,int,int,float); /**< Silberberg & Tsao isotopic production cross section */
   void sigtap_cc(int); /**< initialization of the Barashenkov & Polanski cross section code */
   double sighad_cc(int,double,double,double,double,double); /**< Barashenkov & Polanski pA total cross section */
   inline int fnuc(int z,int a) { return (100 * (z) + (a)); }
   inline int inuc(float b) { return (int)(100 * (b) + 0.1); }
   double He_to_H_CS_ratio(double /**< energy of the primary (GeV/nucleon). */, int /**< primary charge. */, int /**< primary mass. */, int /**< secondary charge. */, int /**< secondary mass. */);
   /**
    * @fn double He_to_H_CS_ratio(double, int, int, int, int)
    * @brief For any given primary (IZI,IAI) and secondary (IZF,IAF) nuclei calculates
    * ratios of the total cross sections (He+AI)/(H+AI) and (He+AI->AF)/(H+AI->AF).
    * @return ratio of the cross sections (He+AI->AF)/(H+AI->AF)
    *
    * REFERENCES:
    * Ferrando P. et al. 1988, PRC 37, 1490
    *    MODIFICATIONS OF THE ORIGINAL SCHEME:
    * - AMU,DELTA are found from linear extrapolation for E1<E(1);
    * - for linear inter-/extra-polation of DELTA the log scale in E is applied;
    * - FZI is found from linear inter-/extra-polation on log scale in Z.
    */
   double nucdata(int, int, int, int, int, int, int*, int*, double*);
   double isotope_cs(double, int, int, int, int, int, int*);
   double eval_cs(double, int, int, int*);
   inline double FI(double X,double X1,double X2,double F1,double F2) { return ((F1-F2)*X+X1*F2-X2*F1)/(X1-X2); } /**< linear interpolation. */
   inline double He_to_H_CS_totratio(int IAI /** Mass of primary */) { return  2.10/pow((double)IAI,0.055); }
   /**
    * @fn inline double He_to_H_CS_totratio(int IAI)
    * @brief ratio of the total cross sections (He+AI)/(H+AI)
    */
   static const int N_DATA_FILES = 4; /**< total number of data-files to read. */
   
   vector<string> data_filename; /**< data files to read. */
   int n_data[3][N_DATA_FILES]; /**< their dimensions n1,n2,n3. */
   int file_no[N_DATA_FILES]; /**< 0=isotope_cs.dat, 1=nucdata.dat etc. */
   float *data_file[N_DATA_FILES]; /**< pointers to the data arrays */
   void Print() {
      ofstream outfile("Galprop_test_6012_5010.dat", ios::out);
      vector<double> result(GetXSec(6, 12, 5, 10));
      for (int i = 0; i < result.size(); ++i) outfile << energy[i] << " " << result[i]/(1e-27*beta[i]*Clight) << endl;
   }
   /**
    * @fn void Print()
    * @brief Print some specified cross section, only for testing purposes.
    */
   
};

double tan_ng(double, void*);

//********************************************************************************************************************************************
//********************************************************************************************************************************************
//********************************************************************************************************************************************

/**
 * @class FlukaXSec
 * @author Daniele Gaggero
 * @author M.Nicola Mazziotta
 * @email daniele.gaggero@sissa.it
 * @brief ...put some description here...
 */
class FlukaXSec : public TXSecBase {
   
public:
   
   //FlukaXSec();
   FlukaXSec(TGrid* co);
   
   ~FlukaXSec() { };
   
   //this is used
   double GetXSec(pair<int,int> /**< Pair of parent-daughter nuclei uid. */, double /**< Energy. */);

   virtual void InitTableInelastic();
  
   double GetXSec(int /**< Nucleus uid. */, double /**< Energy. */);

   inline double GetHefactor() { return 1.26/He_abundance; }
 
   vector<double> GetXSec(int /**< uid of nucleus. */, int /**< K electron */) { throw -1; }

   virtual double GetTotalXSec(int /**< Nucleus uid. */, double /**< Energy. */) {throw -1; }// { /*...*/ }
	
   vector<double> GetTotalXSec_vec(int /**< Nucleus uid. */) { throw -1; }

   map<pair<int,int>, vector<double> > GetTotalXSec_for_each_gas_type() {return totalxsec_extended; }
   
   vector<double> GetTotalXSec_vec(int /**< Nucleus uid. */, int /*target gas nucleus uid */);

   virtual void Print(TGrid*) { }
   
   virtual vector<double> GetXSec(int /**< uid of nucleus. */) { throw -2; }

   virtual vector<double> GetXSec(int /**< Charge of primary. */, int /**< Mass of primary. */, int /**< Charge of secondary. */, int /**< Mass of secondary. */) { throw -2; }

   double GetXSecExtended(pair<int,int>, int, double);

   vector<double> GetXSecExtendedSpectrum(pair<int,int>, int, double);

   virtual void nucleon_cs(int, double /**< kinetic energy per nucleon of the beam momentum particles (GeV) */, int /**< =+/-1 is the charge of beam and At is the atomic number of target nuclei */, int /**< =1 for proton, = -1 for antiproton collisions. */, int /**< =1 for proton collisions. */, double* /**< [mbarn] is the proton-proton inelastic cross sect. */, double* /**< [mbarn] is the proton-nucleus inelastic cross sect. (p+A2) */, double* /**< [mbarn] is the antiproton-proton inelastic cross sect.(nonannihil.) */, double* /**< [mbarn] put =PA_inel, is the antiproton-nucleus inelastic cross sect.(nonan.) */, double* /**< [mbarn] is the antiproton-proton annihilation cross sect. */, double* /**< [mbarn] is the antiproton-nucleus annihilation cross sect. */) { throw -2; }
   virtual void Kcapture_cs(double, int, int, double*, double*) { throw -2; }
   
   virtual bool IsPresent(int uid /**< Nucleus uid. */) {
      map<int, vector<double> >::iterator it = totalxsec.find(uid);
      if (it != totalxsec.end()) return true;
      else return false;
   }
   
   virtual bool IsPresent(pair<int,int> uid /**< Nucleus parent-daughter uid. */) {
      /*map<pair<int,int>, vector<double> >::iterator it = xsec.find(uid);
      if (it != xsec.end()) return true;
      else return false;*/
      return false; // only for now	
   }
  
   virtual double antiproton_cc1(gsl_integration_workspace*,size_t,int,double,double,int,int,int,int) { return -1; }
   
protected:

   //this is actually used!
   map<pair<int,int>, vector<double> > xsec; /**< xsection vector from each <parent, daughter> pair */

   //not used now   
   map<pair<pair<int,int>, int>, vector<double> > xsec_extended; /**< xsection vector for each  <parent, daughter> pair and gas species. */
   map<pair<pair<int,int>, int>, vector< vector<double> > > xsec_extended_spectrum; /**< xsection matrix (parent, daughter energy) for each  <parent, daughter> pair and gas species. */
     
   //to be implemented better
   map<int, vector<double> > totalxsec; /**< inelastic total cross section of a given nucleus. */

   //this is actually used!
   vector<double> energy_vec_ap_H;
   vector<double> inelastic_xsec_ap_H_vec;
   vector<double> energy_vec_ap_He;
   vector<double> inelastic_xsec_ap_He_vec;
   double GetTotalApHXSec(double);
   double GetTotalApHeXSec(double);
   bool initialized ;

   //not used now   
   map<pair<int,int>, vector<double> > totalxsec_extended; /**< inelastic total cross section of a given nucleus on each gas target. */

   vector<double> energy; /**< Energy array for Fluka xsec. */
   vector<double> energy_vec; /**< Energy array for Fluka xsec. */

   vector<double> beta; 
   
   void convert(const int channel /**< Database ID. */, int& uid1 /**< uid of parent nucleus (it will be modified). */, int& uid2 /**< uid of daughter nucleus (it will be modified). */) {
      int A1 = channel/1000000;
      int Z2 = channel%100;
      
      int Z1 = (channel-A1*1000000)/10000;
      uid1 = A1+Z1*1000;
      int A2 = (channel-A1*1000000-Z1*10000)/100;
      uid2 = A2 + Z2*1000;
      return ;
   }

   /**
    * @fn void convert(const int channel, int& uid1, int& uid2)
    * @brief Convert from database to uid for spallation (2 nuclei in one number).
    */
   
   void convert_single(const int channel /**< Database ID. */, int& uid1 /**< uid of nucleus (it will be modified). */) {
      int A1 = channel/100;
      int Z1 = channel%100;
      uid1 = A1+Z1*1000;
      return ;
   }
   /**
    * @fn void convert_single(const int channel, int& uid1)
    * @brief Convert from database to uid for inelastic cross section (1 nucleus in one number).
    */
   
};

#endif
