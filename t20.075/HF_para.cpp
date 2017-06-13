//c++ HF_.cpp -std=c++1y
#include<iostream>
//#include <array>
#include <math.h>       /* cos */
#define PI 3.14159265
#include <vector>
#include <cstdlib>     /* abs */
#include<fstream>// header file for input and output
#include <sstream> 
#include"roots_multidim.h"

using namespace std;
class HF {
  long double t,t2,delta,A,U;
  const static int dim = 2,GRID = 300, bin = 500;
  long double cos_[2*GRID] = {};
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

  VecDoub solve(VecDoub_I x){
 
  auto mu = x[0];
  auto dn = x[1];
  Doub total_no_part = 0.0;
  Doub free_energy_dn = 0.0;
  Doub free_energy = 0.0;
  auto n = x.size();
  VecDoub fvec(n);

  int N_p, N_n;
  N_p = 0;
  N_n = 0;
  //cout<<"delta:"<<delta<<" mu:  "<<mu<< "dn:   "<<dn<<endl;
   for(int i = 0; i<2*GRID; i++){
	for(int j = 0; j<2*GRID; j++){
		auto energy =  -2*t*(cos_[i] + cos_[j]) + 2*dim*t2*cos_[i]*cos_[j];
		auto energy_prime = + 2*t*(cos_[i] + cos_[j]) + 2*dim*t2*cos_[i]*cos_[j];
		auto A = energy + energy_prime + U - 2*mu;
		auto B = energy - energy_prime;
		auto C = sqrt( B*B + 4*(delta-(U*dn/2))*(delta-(U*dn/2)));
		lambda_p[i][j] = 0.5*( A + C ) ;
		lambda_n[i][j] = 0.5*( A - C ) ;

		if (lambda_p[i][j]<0.0){
			N_p = N_p + 1;
			total_no_part = total_no_part + 1.0;
			free_energy_dn  = free_energy_dn - 2*(delta-(U*dn/2))/C;
			free_energy  = free_energy + lambda_p[i][j] + mu;
			//free_energy_con  = free_energy_con + 1.0/C;
			}
		if (lambda_n[i][j]<0.0){ 
			N_n = N_n + 1;
			total_no_part = total_no_part + 1.0;
			free_energy_dn  = free_energy_dn + 2*(delta-(U*dn/2))/C;
			free_energy  = free_energy + lambda_n[i][j] + mu;
			//free_energy_con  = free_energy_con + 1.0/C;
			}
		
	} // for(int i = 0; i<2*GRID; i++)
  } // for(int j = 0; j<2*GRID; j++)

 
 

 fvec[0] = total_no_part/(4*GRID*GRID) - 1.0;
 fvec[1] = -free_energy_dn/(4*GRID*GRID) + dn;
 //cout<<"N_p,N_n,free_energy:"<<N_p<<"N_n: "<<N_n<<"     "<<((free_energy/(4*GRID*GRID))- (U*(1+(dn*dn))/4.0))<<endl;
 return fvec;
 }	//solve

  long double get_GRID(){ return 2*t; };
  auto get_lambda(){return lambda_upp[0][0];}
  void set_U(double U_){U=U_;}


};//class

long double abs(long double a){
	if (a>0.0){
		return a;
	}
	else{
		return -1.0*a;
	}
	}


 VecDoub vecfunc(VecDoub_I x){
  long double t,t2,delta,U,mu,conv,ms,dn;
  int dim,GRID,loop;
  
  VecDoub fnvec;
  t = 0.5;
  t2 = 0.075;
  delta = 0.5;
  dim = 2;
  GRID = 200;
  
  ifstream infile; 
  infile.open("U_file.in"); 
  //if(infile.is_open()){std::cout<<"file occupation_in.dat is open"<<endl;}  
  std::string line;
  std::getline(infile, line);
  std::istringstream ss(line);
  ss >> U ;
  
  HF hf(t,t2,delta,dim);
 
  hf.set_U(U);
  fnvec = hf.solve(x);
  //cout << "U:"<<U<<"  "<<fnvec[0]<<"  "<<fnvec[1]<<endl;
  //cout << "U:"<<U<<"  "<<fnvec[0]<<"    "<<fnvec[1]<<endl;
return fnvec;
}
//0.196066	0.196066	0.803932
int main(){
  bool check = false;
  VecDoub_IO  x(2);
  x[0] = 0.196066;
  x[1] = 0.6;
  vector<long double> U_list ={0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
  //vector<long double> U_list ={ 0.0,    0.5,   1.0,    1.1,   1.2,   1.3};
  //auto u_list_size = U_list.size();
  //vector<long double> mu_list ={6.52500000e-04,   2.82915000e-01,   6.28838000e-01,   6.93954000e-01, 7.58288000e-01,   8.21604000e-01};
  //vector<long double> dn_list = { 0.607855,  0.529807,  0.456882,  0.4434,    0.430344,  0.417724};
  ofstream soumen("broyden_para3.dat");
  soumen <<"# U,mu,ntotal-2,ms,ms1-ms,dn,dn1-dn"<<endl;
  for(long double U : U_list) {
  //for(int i = 0; i<u_list_size; i++){
          //auto U =  U_list[i];
          //x[0] = mu_list[i];
          //x[1] = dn_list[i];
	  check = false;
	  ofstream U_file("U_file.in");	  
	  U_file <<U<<endl;
	  U_file.close();

  	broydn(x, check, *vecfunc );
 	
  	auto m = vecfunc( x);
	auto ms = 0.0;
  	auto mu = x[0];
	auto dn = x[1];
  	// write n values
         cout <<"# U,mu,ntotal-2,ms,ms-ms1,dn,dn1-dn"<<endl;
  	 cout  << U<<"  "<<mu<<"   "<<m[0]<<"    "<<ms<<"    "<<0.0<<"    "<<dn<<"    "<<m[1]<<endl;
  	 soumen << U<<"  "<<mu<<"   "<<m[0]<<"    "<<ms<<"    "<<0.0<<"    "<<dn<<"    "<<m[1]<<endl;
	//if(x[1]<0.1) x[1] = 0.1;
  }
  soumen.close();
  

  return 0;
 //TODO -free_energy_dn/(4*GRID*GRID) + dn is important instad free_energy_dn/(4*GRID*GRID) - dn.
} //main
