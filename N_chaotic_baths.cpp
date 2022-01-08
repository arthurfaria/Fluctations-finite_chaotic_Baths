#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <time.h>
#include <sstream>
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Evol_Temp(int N_bath, double Lambda, int N_cond, int N_iter, ofstream &arquivo){

        double r = (pow(2.0,1.0/3.0) + pow(2.0,-1.0/3.0) -1)/6.0;
        double c[4] = {r + 1.0/2.0, -r, -r, r + 1.0/2.0};
        double d[4] = {2.0*r + 1.0, -4.0*r - 1.0, 2.0*r + 1.0, 0.0};

        double T_step  = 0.01;

        // system coordinates

        double *Q  = (double *) calloc(N_cond, sizeof(double));
        double *P  = (double *) calloc(N_cond, sizeof(double));

        // sample coordinates

        double *x1  = (double *) calloc(N_cond, sizeof(double));
        double *y1  = (double *) calloc(N_cond, sizeof(double));
        double *px1 = (double *) calloc(N_cond, sizeof(double));
        double *py1 = (double *) calloc(N_cond, sizeof(double));

        // heat bath coordinates for each mode

        double **X   = (double **) calloc(N_bath  , sizeof(double*));
        //double **X01   = (double **) calloc(N_bath  , sizeof(double*));
        double **Y   = (double **) calloc(N_bath , sizeof(double*));
        double **Px  = (double **) calloc(N_bath, sizeof(double*));
        double **Py  = (double **) calloc(N_bath,  sizeof(double*));

        int i, j, h, l;

        for(i = 0; i < N_bath; i++){

                X[i]  = (double*) calloc(N_cond, sizeof(double));
                //X01[i]  = (double*) calloc(N_cond, sizeof(double));
                Y[i]  = (double*) calloc(N_cond, sizeof(double));
                Px[i]  = (double*) calloc(N_cond, sizeof(double));
                Py[i]  = (double*) calloc(N_cond, sizeof(double));               

        }

	// quartic parameter 'a'
	double A = 0.1;
	
	//double C_qq0;

        double E_part = 0.1968*(100);
        double E_bath = 0.01/(1000);
        double k = 0.3;
        double Quant;
        double E_Test, E_Harm, Ep_aver, Eb_aver, Ei_aver;

        int tau_relax = 100000;


	// Sampling initial conditions for each heat bath mode.	

	for(i = 0; i < N_bath; i++){

		//cout << i << " de " << N_cond << "\n";
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

	
	// Picking up the 5th bath mode and testing whether it relaxes into chaotic behavior before coupling to the system
	/* With the correlation function for the position of the bath mode it is possible to see if the bath is chaotic
	(see ref: "https://iopscience.iop.org/article/10.1088/1751-8121/ab9a78/meta") */


	/*for(j = 0; j < N_cond; j++){
		
		X01[5][j] = X[5][j];
	}*/



	for(i = 0; i < tau_relax ; i++){

	

		/*C_qq0 = 0.0;	
			
		for(j = 0; j < N_cond; j++){
			
			C_qq0 += X01[5][j]*X[5][j]/N_cond;
		}

		arquivo << showpos << fixed << (i)*T_step << "    ";	
		arquivo << showpos << scientific << C_qq0 << "\n" ;*/

		

		// Updating bath coordinates of each mode until it gets chaotic
		/* We use fourth-order symplectic integration to integrate Hamilton's equations. See:
                 "https://www.sciencedirect.com/science/article/abs/pii/016727899090019L?via%3Dihub"*/	
       
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


 	int temp = 0;

        for(i = 0; i < N_iter ; i++){

                
               // writting a file with energies during the system interaction with the chaotic bath

        	E_Test  = 0.0;
		E_Harm  = 0.0;
               Ep_aver = 0.0;
               Eb_aver = 0.0;
               Ei_aver = 0.0;

               for(j = 0; j < N_cond; j++){

               	Quant = 0.0;

                       for(h= 0; h < N_bath; h++){

                       	Quant     +=  X[h][j];
                              Eb_aver   +=  ((Px[h][j]*Px[h][j] + Py[h][j]*Py[h][j])/2.0 + A*(X[h][j]*X[h][j]*X[h][j]*X[h][j] + Y[h][j]*Y[h][j]*Y[h][j]*Y[h][j])/4.0 + (Y[h][j]*Y[h][j]*X[h][j]*X[h][j])/2.0)/ N_cond ;

                       }

	               E_Harm    +=   (k * Q[j]*Q[j]/(2.0)) / N_cond;
	               Ep_aver   +=   (P[j]*P[j]/(2.0)) / N_cond;
                	Ei_aver   +=   (Lambda * (Quant) * Q[j]) / N_cond;
                }

		E_Test    =   (Eb_aver + Ep_aver +  E_Harm  + Ei_aver);

               if(i == temp ){
               	arquivo << showpos << fixed << (i)*T_step << "    ";
               	arquivo << showpos << scientific << E_Test << "    " << E_Harm << "    " << Ep_aver << "    " << Eb_aver << "    " <<  Ei_aver << "\n" ;
                       // ignoring some data (too heavy) 
			temp += 100;
               }

		// Interaction system-bath
               
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

	free(Q);
	free(P);

	free(x1);
  	free(y1);
      	free(px1);
      	free(py1);

      	free(X);
      	//free(X01);
      	free(Y);
      	free(Px);
      	free(Py);
}

////////////////////////////////////----------Main -----//////////////////////////////////////////////////////

int main(int argn, char **arg){

        if(argn != 5){
                cout << "\n";
                cout << "Entrada Inválida!\n";
                cout << "Entre com as seguintes variáveis: Número de Condições Iniciais, Número de Iterações e o parâmetro A.\n\n";
        }

        // Picking parameters	

    	int N_bath = atoi(arg[1]);  		// number of bath modes
        double Lambda = atof(arg[2]);	// coupling constant 
        int N_cond = atoi(arg[3]);    	// number of initial samples of the system
        int N_iter = atoi(arg[4]);		// number of iterations between system and bath

        double Lambda1 = Lambda/sqrt(N_bath);

        double t_0 = clock();

        ofstream arquivo;
        arquivo.open (("Teste_" + string(arg[1]) + "_" + string(arg[2]) + "_" + string(arg[3]) + "_" + string(arg[4]) + ".txt").c_str());

                Evol_Temp(N_bath, Lambda1, N_cond, N_iter, arquivo);

        arquivo.close();

        double T_process = (clock()-t_0)/CLOCKS_PER_SEC;

        cout << "\n";
        cout << "Cálculo Encerrado! ";
        cout << "Tempo de Processamento [em segundos]: " << T_process;
        cout << "\n\n";

        return 0;
}
