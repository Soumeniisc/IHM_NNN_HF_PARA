//c++ HF_.cpp -std=c++1y
#include<iostream>
//#include <array>
#include <math.h>       /* cos */
#define PI 3.14159265
#include <vector>
#include <cstdlib>     /* abs */
#include<fstream>// header file for input and output
#include <sstream> 
#include <algorithm>
#include"roots_multidim.h"


using namespace std;
class HF {
  long double t,t2,delta,A,U;
  const static int dim = 2,GRID = 300, bin = 500;
  long double cos_[2*GRID] = {};
  long double eigen_energy[2*4*GRID*GRID] = {}; // there is 4*GRID*GRID no of k point in the BZ and first two is for two band (p and n)
  long double lambda_p[2*GRID][2*GRID] = {{}};
  long double lambda_n[2*GRID][2*GRID] = {{}};
  long double lambda_upp[2*GRID][2*GRID] = {{}};
  long double lambda_dnp[2*GRID][2*GRID] = {{}};
  long double lambda_upn[2*GRID][2*GRID] = {{}};
  long double lambda_dnn[2*GRID][2*GRID] = {{}};
  long double L_upp[2*GRID][2*GRID] = {{}};
  long double L_dnp[2*GRID][2*GRID] = {{}};
  long double L_upn[2*GRID][2*GRID] = {{}};
  long double L_dnn[2*GRID][2*GRID] = {{}};
  long double CA_upp[2*GRID][2*GRID] = {{}};
  long double CA_dnp[2*GRID][2*GRID] = {{}};
  long double CA_upn[2*GRID][2*GRID] = {{}};
  long double CA_dnn[2*GRID][2*GRID] = {{}};
  long double CB_upp[2*GRID][2*GRID] = {{}};
  long double CB_dnp[2*GRID][2*GRID] = {{}};
  long double CB_upn[2*GRID][2*GRID] = {{}};
  long double CB_dnn[2*GRID][2*GRID] = {{}};
  long double dos_upp[bin] = {};
  long double dos_upn[bin] = {};
  long double dos_dnp[bin] = {};
  long double dos_dnn[bin] = {};

  public:
  long double nA_upn,nA_dnn,nB_upn,nB_dnn;
  int cnt ;
  HF(long double t_, long double t2_, long double delta_, int dim_){
  	t = t_;
	t2 = t2_;
	delta = delta_;
  
  for(int i = 0; i<2*GRID; i++){
	cos_[i] = cos( (i-GRID)*PI/(GRID));
	}
	} // HF

  long double division(long double a, long double b)
	{
	   if( b == 0 )
	   {
	      throw "Division by zero condition!";
	   }
	   return (a/b);
	}  

  Doub free_energy(VecDoub_I x){
  auto mu = x[0];
  auto dn = x[1];
  auto ms = 0.0;
   int N_p, N_n;
  N_p = 0;
  N_n = 0;
  Doub free_energy_ = 0.0;
  ofstream energy_filen("energy_n.dat");
  ofstream energy_filep("energy_p.dat");
  ofstream energy_file1("energy.dat");
   for(int i = 0; i<(2*GRID); i++){
	for(int j = 0; j<(2*GRID); j++){
		auto energy =  2*t*(cos_[i] + cos_[j]) - 2*dim*t2*cos_[i]*cos_[j];
		auto energy_prime = - 2*t*(cos_[i] + cos_[j]) - 2*dim*t2*cos_[i]*cos_[j];
		auto A = energy + energy_prime + U ;
		auto B = energy - energy_prime;
		auto C = sqrt( B*B + 4*(delta-(U*dn/2))*(delta-(U*dn/2)));
		lambda_p[i][j] = 0.5*( A + C ) ;
		lambda_n[i][j] = 0.5*( A - C ) ;
		eigen_energy[(2*GRID)*i+j] = lambda_n[i][j];
		eigen_energy[4*GRID*GRID+(2*GRID)*i+j] = lambda_p[i][j];
		//energy_file1 << i+j<< "     "<<i<<"   "<<j<<"   "<<eigen_energy[i+j] <<"    "<<C<<endl;
		energy_filen <<i<<"   "<<j<<"   "<<eigen_energy[(2*GRID)*i+j] <<endl;
		energy_filep <<i<<"   "<<j<<"   "<<eigen_energy[4*GRID*GRID+(2*GRID)*i+j] <<endl;
		//energy_file1 <<i+j<<"   "<<eigen_energy[i+j] <<endl;
		//energy_file1 <<4*GRID*GRID+i+j<<"   "<<eigen_energy[4*GRID*GRID+i+j] <<endl;
		
	} // for(int i = 0; i<2*GRID; i++)
  } // for(int j = 0; j<2*GRID; j++)
 sort(eigen_energy, eigen_energy + 2*4*GRID*GRID);
  for(int i = 0; i<4*GRID*GRID; i++){
	free_energy_  = free_energy_ + eigen_energy[i] ;
  }
  
 
 for(int i = 0; i<(4*2*GRID*GRID); i++){
	energy_file1 <<i<<"   "<<eigen_energy[i] <<endl;
 }

 
 return free_energy_/(4*GRID*GRID) - U*(1+(dn*dn))/4.0 ; 
 }	//free_energy

  long double get_GRID(){ return 2*t; };
  auto get_lambda(){return lambda_upp[0][0];}
  void set_U(double U_){U=U_;}
void dos_(double dn, double min_energy, double max_energy){
	auto en_step = (max_energy - min_energy)/ bin;
	int l;		
	//ofstream energy_file1("energy.dat");
	for(int j = 0; j<2*4*GRID*GRID; j++){
		l = int((eigen_energy[j]-min_energy)/en_step);
		dos_upp[l] = dos_upp[l] + 1; 
		//energy_file1 << j<< "     "<<eigen_energy[j] <<endl;

	}
		
	 stringstream dos_file;
	 dos_file<<"dos"<<"delta"<<delta<<"U"<<U <<"t2"<<t2<<"dn"<<dn<<"h.txt";
         ofstream dos(dos_file.str());
         dos << "#energy   dos  GRID="<<GRID<<" bin="<<bin<<endl;
         
	 //std::cout<<"its is correct"<<endl;
     if(dos.is_open()){
         for(int p = 0; p<bin; p++){
	   //std::cout<<p<<endl;
	   dos_upp[p] = dos_upp[p]/(4*GRID*GRID);
	   
	   dos << (min_energy + p*en_step) <<"     "<<dos_upp[p]<<endl;
	}//for(p=)
      }	//if
  }//dos



};//class

long double abs(long double a){
	if (a>0.0){
		return a;
	}
	else{
		return -1.0*a;
	}
	}



int main(){
  double m,t,t2,delta,U,mu,ms,dn;
  double long free_energy_;
  int dim,GRID,n;
  VecDoub_IO  x(2);
  t = 0.5;
  t2 = 0.1;
  delta = 0.5;
  dim = 2;
  GRID = 200;
  U = 1.5;
  mu = 0.05;

  HF hf(t,t2,delta,dim);
  hf.set_U(U);

  stringstream energy_file;
  energy_file<<"energy_vs_dn_U"<<U<<"_.dat";
  ofstream soumen(energy_file.str());
  soumen <<"# dn, free_energy"<<endl;
  //vector<long double> dn_list ={0.1,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7};
  //vector<long double> U_list ={0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.56};
  //vector<long double> dn_list  = {0.3};
  vector<long double> dn_list  = {-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1};
  for(long double dn : dn_list) {
  	std::cout<<U <<"  mu:"<<mu <<"dn:"<< dn<<endl;
  	x[0] = mu;
  	x[1] = dn;
	free_energy_ = hf.free_energy(x);
	soumen << dn<<"  "<<free_energy_- U*(1+(dn*dn))/4.0<<"    "<<- U*(1+(dn*dn))/4.0<<"   "<<free_energy_<<std::endl;
	std::cout<<"dn:"<< dn<<"  "<<"free_energy_:"<<free_energy_<<std::endl;
	//hf.dos_(dn, -5.0, 5.0);
    
   }//dn_list
  
  return 0;

} //main
