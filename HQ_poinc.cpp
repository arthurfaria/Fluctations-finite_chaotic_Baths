#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <time.h>
#include <sstream>
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Evol_Temp(int N_cond, int N_iter, double A, ofstream &arquivo){

	double x = (pow(2.0,1.0/3.0) + pow(2.0,-1.0/3.0) -1)/6.0;
	double c[4] = {x+1.0/2.0, -x, -x, x+1.0/2.0};
	double d[4] = {2.0*x+1.0, -4.0*x-1.0, 2.0*x+1.0, 0.0};
	
	
	// time step
	double T_step  = 0.001;
	double E_bath = 0.5;
	

	// tolerance for slicing a Poincare section
	double tol = 0.001;  
  	double q1, q2, p1, p2;

	int j;

	// Sampling initial conditions for each heat bath mode.	

	for(int i = 0; i < N_cond; i++){
	
		cout << i << " de " << N_cond << "\n";
		
		double m = rand()%4;

		if(m == 0){
			q2 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath/A, 1.0/4.0);
			p1 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath, 1.0/2.0);
			p2 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath, 1.0/2.0);
			q1 = pow(-1,rand()%2) * sqrt(-q2*q2/(A) + sqrt(q2*q2*q2*q2 - A*(2*p1*p1 + 2*p2*p2 + A*q2*q2*q2*q2 - 4*E_bath))/(A));
		}
		else if(m == 1){
			q1 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath/A, 1.0/4.0);
			p1 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath, 1.0/2.0);
			p2 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath, 1.0/2.0);
			q2 = pow(-1,rand()%2) * sqrt(-q1*q1/(A) + sqrt(q1*q1*q1*q1 -A*(2*p1*p1 + 2*p2*p2 + A*q1*q1*q1*q1 - 4*E_bath))/(A));
		}
		else if(m == 2){
			q1 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath/A, 1.0/4.0);
			q2 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath/A, 1.0/4.0);
			p2 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath, 1.0/2.0);
			p1 = pow(-1,rand()%2) * sqrt(2*E_bath - A*(q1*q1*q1*q1 + q2*q2*q2*q2)/2.0 - q1*q1*q2*q2 - p2*p2);
		}
		else{
			q1 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath/A, 1.0/4.0);
			q2 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath/A, 1.0/4.0);
			p1 = (2.0*rand()/RAND_MAX - 1) * pow(2 * E_bath, 1.0/2.0);
			p2 = pow(-1,rand()%2) * sqrt(2*E_bath - A*(q1*q1*q1*q1 + q2*q2*q2*q2)/2.0 - q1*q1*q2*q2 - p1*p1);
		}

       if(!isnan(q1) && !isnan(q2) && !isnan(p1) && !isnan(p2)){	

			j = 0;			
			while(j < N_iter){
				
				for(int l = 0; l < 4; l++){
					p1 += - c[l]*T_step * (A*q1*q1*q1 + q1*q2*q2);  
					p2 += - c[l]*T_step * (A*q2*q2*q2 + q2*q1*q1);  
					q1 += + d[l]*T_step * p1;
					q2 += + d[l]*T_step * p2;			
				}
				// writting a file with Poincare section phase-space variables

				if(abs(q2) < tol && p2 > 0){ 
					arquivo << showpos << scientific << q1 << "    " << p1 << "        " ;
					arquivo << fixed << ( (p1*p1 + p2*p2)/2 + A*(q1*q1*q1*q1 + q2*q2*q2*q2)/4.0 + (q2*q2*q1*q1)/2.0 )<<  "\n";
					j ++ ;
			  	}
			}
		}	 
	}
}

////////////////////////////////////----------Main -----//////////////////////////////////////////////////////

int main(int argn, char **arg){		

	if(argn != 4){
		cout << "\n";
		cout << "Entrada Inválida!\n";
		cout << "Entre com as seguintes variáveis: Número de Condições Iniciais, Número de Iterações e o parâmetro A.\n\n";
	}
	
	// Defining parameters	
	int N_cond = atoi(arg[1]);	// number of bath modes
	int N_iter = atoi(arg[2]);	// number of iterations between system and bath 
	double A = atof(arg[3]);	// quartic parameter 'a'

	double t_0 = clock();
	
	ofstream arquivo;
	arquivo.open (("HQ_" + string(arg[1]) + "_" + string(arg[2]) + "_" + static_cast <ostringstream*> (&(ostringstream() << A)) -> str() + ".txt").c_str());
	
		Evol_Temp(N_cond, N_iter, A, arquivo);

	arquivo.close();
	
	double T_process = (clock()-t_0)/CLOCKS_PER_SEC;
	
	cout << "\n";
    	cout << "Cálculo Encerrado! "; 
	cout << "Tempo de Processamento [em segundos]: " << T_process;
	cout << "\n\n";	

return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
