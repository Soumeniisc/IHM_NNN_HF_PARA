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
 
  auto mu = U/2.0;
  auto dn = x[0];
  Doub total_no_part = 0.0;
  Doub free_energy_dn = 0.0;
  Doub free_energy = 0.0;
  auto n = x.size();
  VecDoub fvec(n);

  int N_p, N_n;
  N_p = 0;
  N_n = 0;
  cout<<"delta:"<<delta<<endl;
   for(int i = 0; i<2*GRID; i++){
	for(int j = 0; j<2*GRID; j++){
		auto energy =  2*t*(cos_[i] + cos_[j]) - 2*dim*t2*cos_[i]*cos_[j];
		auto energy_prime = - 2*t*(cos_[i] + cos_[j]) - 2*dim*t2*cos_[i]*cos_[j];
		auto A = energy + energy_prime + U - 2*mu;
		auto B = energy - energy_prime;
		auto C = sqrt( B*B + 4*(delta-(U*dn/2))*(delta-(U*dn/2)));
		lambda_p[i][j] = 0.5*( A + C ) ;
		lambda_n[i][j] = 0.5*( A - C ) ;

		if (lambda_p[i][j]<0.0){
			N_p = N_p + 1;
			//total_no_part = total_no_part + 1.0;
			free_energy_dn  = free_energy_dn -2*(delta-(U*dn/2))/C;
			//free_energy_con  = free_energy_con + 1.0/C;
			}
		if (lambda_n[i][j]<0.0){ 
			N_n = N_n + 1;
			//total_no_part = total_no_part + 1.0;
			free_energy_dn  = free_energy_dn + 2*(delta-(U*dn/2))/C;
			free_energy  = free_energy + lambda_n[i][j] + mu;
			//free_energy_con  = free_energy_con + 1.0/C;
			}
		
	} // for(int i = 0; i<2*GRID; i++)
  } // for(int j = 0; j<2*GRID; j++)

 
 

 //fvec[0] = total_no_part/(4*GRID*GRID) - 1.0;
 fvec[0] = free_energy_dn/(4*GRID*GRID) - dn;
 cout<<"N_p,N_n,free_energy:"<<N_p<<" "<<N_n<<"     "<<((free_energy/(4*GRID*GRID))- (U*(1+(dn*dn))/4.0))<<endl;
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
  t2 = 0.0;
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
  cout << "U:"<<U<<"  "<<fnvec[0]<<endl;
return fnvec;
}
//1.05   0    0    0    0.332472
int main(){
  bool check = false;
  VecDoub_IO  x(1);
  x[0] = 0.332472;
  //x[1] = 0.0;
  //vector<long double> U_list ={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1};
  vector<long double> U_list ={2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9};
  ofstream soumen("broyden_para2.dat");
  soumen <<"# U,mu,ntotal-2,ms,ms1-ms,dn,dn1-dn"<<endl;
  for(long double U : U_list) {
	  check = false;
	  ofstream U_file("U_file.in");	  
	  U_file <<U<<endl;
	  U_file.close();

  	broydn(x, check, *vecfunc );
 	
  	auto m = vecfunc( x);
  	auto mu = U/2.0;
        //auto ms = x[1];
	auto ms = 0.0;
	auto dn = x[0];
  	
  	// write n values
         cout <<"# U,mu,ntotal-2,ms,ms-ms1,dn,dn1-dn"<<endl;
  	 cout  << U<<"  "<<mu<<"   "<<0.0<<"    "<<ms<<"    "<<0.0<<"    "<<dn<<"    "<<m[0]<<endl;
  	 soumen << U<<"  "<<mu<<"   "<<0.0<<"    "<<ms<<"    "<<0.0<<"    "<<dn<<"    "<<m[0]<<endl;
	//if(x[1]<0.1) x[1] = 0.1;
  }
  soumen.close();
  

  return 0;

} //main
