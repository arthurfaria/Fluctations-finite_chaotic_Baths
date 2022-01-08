#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <time.h>
#include <sstream>
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Evol_Temp(int N_bath, double Lambda, double W, int N_cond, ofstream &arquivo){

        double r = (pow(2.0,1.0/3.0) + pow(2.0,-1.0/3.0) -1)/6.0;
        double c[4] = {r + 1.0/2.0, -r, -r, r + 1.0/2.0};
        double d[4] = {2.0*r + 1.0, -4.0*r - 1.0, 2.0*r + 1.0, 0.0};

        double T_step  = 0.01;

        double *Q  = (double *) calloc(N_cond, sizeof(double));
        double *P  = (double *) calloc(N_cond, sizeof(double));

        double *x1  = (double *) calloc(N_cond, sizeof(double));
        double *y1  = (double *) calloc(N_cond, sizeof(double));
        double *px1 = (double *) calloc(N_cond, sizeof(double));
        double *py1 = (double *) calloc(N_cond, sizeof(double));

        double **X   = (double **) calloc(N_bath  , sizeof(double*));
        double **Y   = (double **) calloc(N_bath , sizeof(double*));
        double **Px  = (double **) calloc(N_bath, sizeof(double*));
        double **Py  = (double **) calloc(N_bath,  sizeof(double*));

        int i, j, h, l;

        for(i = 0; i < N_bath; i++){

                X[i]  = (double*) calloc(N_cond, sizeof(double));
                Y[i]  = (double*) calloc(N_cond, sizeof(double));
                Px[i]  = (double*) calloc(N_cond, sizeof(double));
                Py[i]  = (double*) calloc(N_cond, sizeof(double));               

        }

        double *E_i = (double *) calloc(N_cond, sizeof(double));
        double *E_f = (double *) calloc(N_cond, sizeof(double));

	// quartic parameter 'a'
	
	double A = 0.1; 

	// specifying the system energy after the forward protocol
    
	double C = (3.0 * 10.0)/2.0;
  	double R = 2.0/(3.0 * 10.0);
      	double U0 = 1.28689490771309;

      	double E_bath = 0.01/(1000);
      	double E_part = ((1 + R)/R) * (U0 + W/(1 + C)) - E_bath; 

	double k0 = 0.9 ;
	double delta_k = 0.6;
	double k = k0;

     	double Quant;

     	int tau_relax = 100000;	// chaotic relaxation of bath
	int tau = 200000; 		// thermalization system-bath
	int tau_prot = 40000;		// time of protocol

	// Sampling initial conditions for each heat bath mode.
	
	for(i = 0; i < N_bath; i++){

		j = 0;
		
		while(j < N_cond){

                	double m = rand()%4;

                	if(m == 0){

                        	y1[j] = (2.0*rand()/RAND_MAX - 1) / pow(A,1.0/4.0);
                        	px1[j] = (2.0*rand()/RAND_MAX - 1) * pow(2.0 * E_bath, 1.0/2.0);
                        	py1[j] = (2.0*rand()/RAND_MAX - 1) * pow(2.0 * E_bath, 1.0/2.0);
                        	x1[j] = pow(-1,rand()%2) * sqrt(-y1[j]*y1[j] + sqrt(y1[j]*y1[j]*y1[j]*y1[j] - A*(2*px1[j]*px1[j] + 2*py1[j]*py1[j] + A*y1[j]*y1[j]*y1[j]*y1[j] - 4*E_bath)))/sqrt(A);
			
                	}

                	else if(m == 1){

                        	x1[j] = (2.0*rand()/RAND_MAX - 1) * pow(4.0 * E_bath/A, 1.0/4.0);
                        	px1[j] = (2.0*rand()/RAND_MAX - 1) * pow(2.0 * E_bath, 1.0/2.0);
                        	py1[j] = (2.0*rand()/RAND_MAX - 1) * pow(2.0 * E_bath, 1.0/2.0);
                        	y1[j] = pow(-1,rand()%2) * sqrt(-x1[j]*x1[j] + sqrt(x1[j]*x1[j]*x1[j]*x1[j] -A*(2*px1[j]*px1[j] + 2*py1[j]*py1[j] + A*x1[j]*x1[j]*x1[j]*x1[j] - 4*E_bath)))/sqrt(A);
			
                	}

                	else if(m == 2){

 				x1[j] = (2.0*rand()/RAND_MAX - 1) * pow(4.0 * E_bath/A, 1.0/4.0);
                        	y1[j] = (2.0*rand()/RAND_MAX - 1) * pow(4.0 * E_bath/A, 1.0/4.0);
                        	py1[j] = (2.0*rand()/RAND_MAX - 1) * pow(2.0 * E_bath, 1.0/2.0);
                        	px1[j] = pow(-1,rand()%2) * sqrt(2*E_bath - A*(x1[j]*x1[j]*x1[j]*x1[j] + y1[j]*y1[j]*y1[j]*y1[j])/2 - x1[j]*x1[j]*y1[j]*y1[j] - py1[j]*py1[j]);
			

                	}

                	else{
				x1[j] = (2.0*rand()/RAND_MAX - 1) * pow(4.0 * E_bath/A, 1.0/4.0);
                        	y1[j] = (2.0*rand()/RAND_MAX - 1) * pow(4.0 * E_bath/A, 1.0/4.0);
                        	px1[j] = (2.0*rand()/RAND_MAX - 1) * pow(2.0 * E_bath, 1.0/2.0);
                        	py1[j] = pow(-1,rand()%2) * sqrt(2*E_bath - A*(x1[j]*x1[j]*x1[j]*x1[j] + y1[j]*y1[j]*y1[j]*y1[j])/2 - x1[j]*x1[j]*y1[j]*y1[j] - px1[j]*px1[j]);
			
                	}

               		if(!isnan(x1[j]) && !isnan(y1[j]) && !isnan(px1[j]) && !isnan(py1[j])){

                        	X[i][j]  = x1[j];
                        	Y[i][j]  = y1[j];
                        	Px[i][j] = px1[j];
                        	Py[i][j] = py1[j];

				Q[j] = sqrt(2.0 * E_part / k) * sin(j*2.0*M_PI/N_cond);
                		P[j] = sqrt(2.0 * E_part) * cos(j*2.0*M_PI/N_cond);

				j++;
			}

		}

	}

	/* Relaxing the bath modes into chaotic behavior before coupling to the system
	(see ref: "https://iopscience.iop.org/article/10.1088/1751-8121/ab9a78/meta") */

	for(i = 0; i < tau_relax ; i++){

		for(j = 0; j < N_cond; j++){

                	for(l = 0; l < 4; l++){

				for(h = 0; h < N_bath; h++){

                                	Px[h][j] += - c[l] * T_step * (A * X[h][j]*X[h][j]*X[h][j] + X[h][j]*Y[h][j]*Y[h][j]);
                                        Py[h][j] += - c[l] * T_step * (A * Y[h][j]*Y[h][j]*Y[h][j] + Y[h][j]*X[h][j]*X[h][j]);

				}

				for(h = 0; h < N_bath; h++){

                                        X[h][j]  += + d[l] * T_step * Px[h][j];
                                        Y[h][j]  += + d[l] * T_step * Py[h][j];

				}                        
			}  
           
		}
	}

	// Interaction system-bath until thermalization
	/* We use fourth-order symplectic integration to integrate Hamilton's equations. See:
                 "https://www.sciencedirect.com/science/article/abs/pii/016727899090019L?via%3Dihub"*/	

	for(i = 0; i < tau ; i++){

    		for(j = 0; j < N_cond; j++){

	                for(int l = 0; l < 4; l++){

                        	Quant = 0.0;

                                for(h= 0; h < N_bath; h++){

		                        Quant += X[h][j];
                                }

 				P[j] += - c[l] * T_step * (Lambda * (Quant) +  k * Q[j]);
	
				for(h= 0; h < N_bath; h++){

					Px[h][j] += - c[l] * T_step * (A * X[h][j]*X[h][j]*X[h][j] + X[h][j]*Y[h][j]*Y[h][j] + Lambda * Q[j]);
                                	Py[h][j] += - c[l] * T_step * (A * Y[h][j]*Y[h][j]*Y[h][j] + Y[h][j]*X[h][j]*X[h][j]);
      
					X[h][j]  += + d[l] * T_step * Px[h][j];
                                	Y[h][j]  += + d[l] * T_step * Py[h][j];
                           	}   
				// Position should be the last one to be implemented
	                        
				Q[j]    += + d[l] * T_step * P[j];               
        		}
        	}
     	}

	// Operating the protocol on the system after thermalization

	for(j = 0; j < N_cond; j++){

		E_i[j] =  P[j]*P[j]/2.0 +  k * Q[j]*Q[j]/2.0 ;

	}

	// Operating the reversed protocol on the system after thermalization

	for(i = 0; i < tau_prot; i++){
		
		k = k0 - delta_k*(i)/(tau_prot);		

		for(j = 0; j < N_cond; j++){

			for(int l = 0; l < 4; l++){

					P[j]  += - c[l]*T_step * (k * Q[j]);	
					Q[j]  += + d[l]*T_step * P[j];
 						
			}
		}
	
	}
	
	// writting a file with energies during the system interaction with the chaotic bath

	for(j = 0; j < N_cond; j++){

		E_f[j] =  P[j]*P[j]/2.0 +  k * Q[j]*Q[j]/2.0 ;

		arquivo << showpos << scientific << E_f[j] << "    " << E_i[j] << "    "<< k << "    " <<  E_f[j] - E_i[j]  << "\n" ;
	}


        free(Q);
        free(P);

        free(x1);
        free(y1);
        free(px1);
        free(py1);

        free(X);
        free(Y);
        free(Px);
        free(Py);

        free(E_i);
        free(E_f);
}

////////////////////////////////////----------Main -----//////////////////////////////////////////////////////

int main(int argn, char **arg){

        if(argn != 4){
                cout << "\n";
                cout << "Entrada Inválida!\n";
                cout << "Entre com as seguintes variáveis: Número de Condições Iniciais, Número de Iterações e o parâmetro A.\n\n";
        }

	// Picking parameters	
	int N_bath = atoi(arg[1]);		// number of bath modes
     	double Lambda = atof(arg[2]);		// coupling constant 
     	int N_cond = atoi(arg[3]);		// number of initial samples of the system

	int i;
     	double Lambda1 = Lambda/sqrt(N_bath);

     	double t_0 = clock();

	double W[35] = {1.9109e-01, 1.6063e-01, 1.3230e-01, 1.0568e-01, 8.7720e-02, 6.9160e-02, 5.6150e-02, 4.4240e-02, 3.4680e-02, 2.7460e-02, 2.1750e-02, 1.6590e-02, 1.3070e-02, 9.8700e-03, 8.0800e-03, 5.3400e-03, 4.4300e-03, 3.2200e-03, 2.2800e-03, 1.6900e-03, 1.3500e-03, 8.2000e-04, 7.9000e-04, 4.7000e-04,3.6000e-04, 2.2000e-04, 1.8000e-04, 1.3000e-04, 8.0000e-05, 6.0000e-05, 2.0000e-05, 2.0000e-05, 3.0000e-05, 2.0000e-05, 2.0000e-05};	

	for(i = 0; i < 35; i++){

		ofstream arquivo;
		arquivo.open (("Rev_"  + static_cast<ostringstream*>(&(ostringstream()<< i+1))->str() +  "_" + string(arg[1]) + "_" + string(arg[2]) + "_" + string(arg[3]) + ".txt").c_str());
	
		Evol_Temp(N_bath, Lambda1, W[i], N_cond, arquivo);

	arquivo.close();
	}

        double T_process = (clock()-t_0)/CLOCKS_PER_SEC;

        cout << "\n";
        cout << "Cálculo Encerrado! ";
        cout << "Tempo de Processamento [em segundos]: " << T_process;
        cout << "\n\n";

        return 0;
}
