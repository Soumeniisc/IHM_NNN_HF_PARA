//c++ HF_.cpp -std=c++1y
#include<iostream>
//#include <array>
#include <math.h>       /* cos */
#define PI 3.14159265
#include <vector>
#include <cstdlib>     /* abs */
#include<fstream>// header file for input and output
#include <sstream> 

using namespace std;
class HF {
  double t,t2,delta,A;
  const static int dim = 2,GRID = 1000, bin = 700;
  double cos_[2*GRID] = {};
  long double lambda_p[2*GRID][2*GRID] = {{}};
  long double lambda_n[2*GRID][2*GRID] = {{}};
  double lambda_upp[2*GRID][2*GRID] = {{}};
  double lambda_dnp[2*GRID][2*GRID] = {{}};
  double lambda_upn[2*GRID][2*GRID] = {{}};
  double lambda_dnn[2*GRID][2*GRID] = {{}};
  double L_upp[2*GRID][2*GRID] = {{}};
  double L_dnp[2*GRID][2*GRID] = {{}};
  double L_upn[2*GRID][2*GRID] = {{}};
  double L_dnn[2*GRID][2*GRID] = {{}};
  double CA_upp[2*GRID][2*GRID] = {{}};
  double CA_dnp[2*GRID][2*GRID] = {{}};
  double CA_upn[2*GRID][2*GRID] = {{}};
  double CA_dnn[2*GRID][2*GRID] = {{}};
  double CB_upp[2*GRID][2*GRID] = {{}};
  double CB_dnp[2*GRID][2*GRID] = {{}};
  double CB_upn[2*GRID][2*GRID] = {{}};
  double CB_dnn[2*GRID][2*GRID] = {{}};
  double dos_p[bin] = {};
  double dos_n[bin] = {};

  public:
  double nA_upn,nA_dnn,nB_upn,nB_dnn;
  int cnt ;
  HF(double t_, double t2_, double delta_, int dim_){
  	t = t_;
	t2 = t2_;
	delta = delta_;
  
  for(int i = 0; i<2*GRID; i++){
	cos_[i] = cos( (i-GRID)*PI/(GRID));
	}
	} // HF

  double division(double a, double b)
	{
	   if( b == 0 )
	   {
	      throw "Division by zero condition!";
	   }
	   return (a/b);
	}  

 
  double get_GRID(){ return 2*t; };
  auto get_lambda(){return lambda_upp[0][0];}

  void dos_(double U, double mu, double dn, double min_energy, double max_energy){
	auto en_step = (max_energy - min_energy)/ bin;
	int l;		
	for(int i = 0; i<2*GRID; i++){
	   for(int j = 0; j<2*GRID; j++){
		auto energy =  - 2*t*(cos_[i] + cos_[j]) + 2*dim*t2*cos_[i]*cos_[j];
		auto energy_prime = + 2*t*(cos_[i] + cos_[j]) + 2*dim*t2*cos_[i]*cos_[j];
		auto A = energy + energy_prime + U - 2*mu;
		auto B = energy - energy_prime;
		auto C = sqrt( B*B + 4*(delta-(U*dn/2))*(delta-(U*dn/2)));
		lambda_p[i][j] = 0.5*( A + C ) - min_energy ;
		lambda_n[i][j] = 0.5*( A - C ) - min_energy;

		
		l = int(lambda_p[i][j]/en_step);
		dos_p[l] = dos_p[l] + 1; 		
		l = int(lambda_n[i][j]/en_step);
		dos_n[l] = dos_n[l] + 1;
		
	} // for(int i = 0; i<2*GRID; i++)
    } // for(int j = 0; j<2*GRID; j++)
	 stringstream dos_file;
	 dos_file<<"dos"<<"delta"<<delta<<"U"<<U <<"t2"<<t2<<".txt";
         ofstream dos(dos_file.str());
         dos<< "# U,mu,nA_up, nA_dn, nB_up, nB_dn, ms, mf, ntotal,mu,cint"<<endl;
         dos << "#energy    dos_upn   dos_upp   dos_dnn  dos_dnp GRID="<<GRID<<" bin="<<bin<<endl;
	 //std::cout<<"its is correct"<<endl;
     if(dos.is_open()){
         for(int p = 0; p<bin; p++){
	   //std::cout<<p<<endl;
	   dos_p[p] = dos_p[p]/(4*GRID*GRID);
	   dos_n[p] = dos_n[p]/(4*GRID*GRID);
	   dos << (min_energy + p*en_step) <<"     "<<dos_n[p]<<"      "<<dos_p[p]<<endl;
	}//for(p=)
      }	//if
      //else{ continue;}
	//std::cout<<"its is correct"<<endl;
  }//dos

};//class


int main(){
  double m,t,t2,delta,U,mu,dn;
  int dim,GRID,n;
  t = 0.5;
  t2 = 0.05;
  delta = 0.5;
  dim = 2;
  GRID = 200;
  U = 0.0;
  mu = 0.0;

  HF hf(t,t2,delta,dim);

  //cout<<hf.get_GRID()<<endl;
  auto lambda = hf.get_lambda();
  //cout<< GRID <<hf.GRID<< hf.t<<endl;
  //hf.lambda_upp[0][0] = 100.0;
  //cout<<hf.get_lambda()<<endl;
  
  stringstream data_file;
  //data_file<<"delta"<<delta<<"t2"<<t2<<".txt";
  data_file<<"broyden_para3.dat";
  ifstream infile; 
  infile.open(data_file.str()); 
  //infile.open("delta0.5t20.5.txt"); 
  if(infile.is_open()){std::cout<<"file"<<data_file.str()<<" is open"<<endl;}
  //soumen <<"# U,mu,nA_up, nA_dn, nB_up, nB_dn, ms, mf, ntotal,mu,cint"<<endl;
  //0.471591	0.466898	0.527911	0.53287
  //0.581797   0.226212	0.226212	0.7738	0.7738
  //nA_up = 0.610637;
  //nA_dn = 0.203252;
  //nB_up = 0.693188;
  //nB_dn = 0.492885;
  //vector<double> U_list ={0.0,0.5,1.0,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,3.0};
  //vector<double> U_list ={2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,4,4.5,5};//,0.1,0.2,0.4,0.43,0.46,0.49,0.6};
  //vector<double> U_list ={0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.56};
  
   std::string line;
   while (std::getline(infile, line)) {
	//std::cout<<line<<endl;
	if(line[0] != '#'){
		//std::cout<<line<<endl;
   		std::istringstream ss(line);
		ss >> U >> mu >> m >> m >> m >> dn>>m;
		std::cout<<U <<"  "<< mu <<"  "<< dn <<endl;
		//dos_(double U, double mu, double dn, double min_energy, double max_energy)
		
 		hf.dos_(U, mu, dn, -5, 5);
	}//if
   }//while
  
  return 0;

} //main
