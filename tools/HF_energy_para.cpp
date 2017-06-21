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

  Doub free_energy(VecDoub_I x){
  auto mu = x[0];
  auto dn = x[1];
  auto ms = 0.0;
  
  Doub free_energy_ = 0.0;

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
			
			free_energy_  = free_energy_ + lambda_p[i][j] + mu ;
			}
		if (lambda_n[i][j]<0.0){ 
			
			free_energy_  = free_energy_ + lambda_n[i][j] + mu ;
			}
		
	} // for(int i = 0; i<2*GRID; i++)
  } // for(int j = 0; j<2*GRID; j++)

 
 

 
 
 return free_energy_/(4*GRID*GRID) - U*(1+(dn*dn))/4.0 ; 
 }	//free_energy

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



int main(){
  double m,t,t2,delta,U,mu,ms,dn;
  double long free_energy_;
  int dim,GRID,n;
  VecDoub_IO  x(2);
  t = 0.5;
  t2 = 0.0;
  delta = 0.5;
  dim = 2;
  GRID = 200;
  U = 0.0;
  mu = 0.0;

  HF hf(t,t2,delta,dim);


  stringstream energy_file;
  energy_file<<"brouden_energy_.dat";
  //energy_file<<"brouden_AFM_energy.dat";
  ofstream soumen(energy_file.str());
  soumen <<"# U, free_energy"<<endl;

  stringstream data_file;
  data_file<<"delta0.5t20_para.txt";
  //data_file<<"brouden_AFM.dat";
  //std::cout<<"out put file:"<<data_file<<std::endl;
  ifstream infile; 
  infile.open(data_file.str());  
  if(infile.is_open()){std::cout<<"file"<<data_file.str()<<" is open"<<endl;}
  else{std::cout<<"no file:"<<data_file.str()<<" in working directroy"<<std::endl;}
   std::string line;
   while (std::getline(infile, line)) {
	//std::cout<<line<<endl;
	if(line[0] != '#'){
		//std::cout<<line<<endl;
   		std::istringstream ss(line);
		ss >> U >> mu >> ms >> m>> m>> dn >> m ; 
		std::cout<<U <<"  mu:"<<mu <<"ms:"<< ms<<endl;
		x[0] = mu;
		x[1] = ms;
		hf.set_U(U);
	        free_energy_ = hf.free_energy(x);
		soumen << U<<"  "<<free_energy_<<std::endl;
		std::cout<<"free_energy_:"<<free_energy_<<std::endl;
	}//if
    
   }//while
  
  return 0;

} //main