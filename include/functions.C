
void modulate(double *Ekin, double  *spectrum, int np, int z, int a, double phi){

  double T,y,y2,M;
  double Mp=939.e-3; // Mp = (mp+mn)/2 ~939 MeV
  double Me=0.51e-3;
  int j;

  for(int i=0; i<np; i++){
    if(a>0) {T = Ekin[i]+abs(z)*phi/a; M=Mp;}  // ekin is per nucleon
    else if(a==1) {T = Ekin[i]+abs(z)*phi; M=Mp;}
    else  {T = Ekin[i]+abs(z)*phi; M=Me; }

for(j=0;j<np;j++){
if(T<Ekin[j]) break;
}

    if(j==np-1) {  // last bin is not touched 			
      spectrum[i] = spectrum[i]; 
     break; 
    }

    y = spectrum[j-1] +(T-Ekin[j-1]) *(spectrum[j]-spectrum[j-1]) /(Ekin[j]-Ekin[j-1]);

    spectrum[i]= y * Ekin[i] *(Ekin[i]+2*M) /T /(T+2*M);
  }
}



double Calc(string name,const int Edim,double *energy , int n,double *eng_exp,double *err_exp,double *val_exp, double *yvalue, double Cut_off){



double chi2=0.0;

double dummy;
double tval;
double erg1=0.0,erg2;



for (int i=0;i<n;i++){

dummy=0.0;
tval=eng_exp[i];
if(eng_exp[i]>=Cut_off){
spline_quadratic_val( Edim, energy, yvalue, tval, &erg1, &erg2 );

dummy=pow(((erg1-val_exp[i])/(err_exp[i])),2.0);
chi2+=dummy;
}
// cout << name << " Mod " << erg1 <<  " Exp " << val_exp[i] << " err " << err_exp[i] << " chi2 " << dummy << endl;
}
//cout << "chi2-" << name << "\t\t"  <<chi2 << endl;


return chi2;

}

double CalcAsym(string name,const int Edim, double *energy , int n,double *eng_exp,double *err_exp,double *err_exp_down, double *val_exp, double *yvalue ,double Cut_off){



double chi2=0.0;

double dummy;
double tval;
double erg1=0.0,erg2;



for (int i=0;i<n;i++){

dummy=0.0;
tval=eng_exp[i];
if(eng_exp[i]>=Cut_off){
spline_quadratic_val( Edim, energy, yvalue, tval, &erg1, &erg2 );


if((erg1-val_exp[i])>=0) dummy=pow(((erg1-val_exp[i])/(err_exp[i])),2.0); 
if((erg1-val_exp[i]) <0) dummy=pow(((erg1-val_exp[i])/(err_exp_down[i])),2.0);


chi2+=dummy;
}
//cout << "check : " << double(dummy/23.0) << endl;
}

return chi2;

}





int Get_N_data(int n,double *eng_exp,double Cut_off){

int N_data=0;


for (int i=0;i<n;i++){

if(eng_exp[i]>=Cut_off) N_data++;

}

return N_data;
}




