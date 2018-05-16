#include <stdio.h>
#include <math.h>
#include "matrix_alg.h"
/*
mpc test program
states: x: [phi, theta, psi, p, q, r, u, v, w, x, y, z]
input u: [thrust, taux, tauy, tauz]
output:	[phi, theta, psi, x, y, z]

linearized model refer to https://www.kth.se/polopoly_fs/1.588039!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf

A few main parts
1. getting current states
formulate augmented model A, B, C
2. formulating the cost function variables
3. solve for the minimum cost
4. get the first input
*/

#define IX 	1.0 	//moment of inertia X axis
#define IY 	1.0 	//moment of inertia X axis
#define IZ 	1.0 	//moment of inertia X axis
#define M 	1.5 	//mass of copter in kg
#define G 	9.81	//gravity constant
#define NP 	500 	//prediction horizon 
#define N 	7		//number of laguarre function in each laguarre vector
#define ALPHA 0.3	//laguarre approximation tuning parameter
#define delta_T 0.01	//discrete time interval

float (*matrix_power(float A[][18], int n))[18];


main() {

	printf("Hello mpc\n");
	//getting model states
	float x[12] ={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

	float x_prev[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

	float y[6] = {x[0], x[1], x[2], x[9], x[10], x[11]};

	float y_des[6] = {0, 0, 0, 1, 0, 1};

	float u[4] = {M*G, 0, 0, 0};

	float delta_u[4];

	float X_aug[18]=
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, y[0]-y_des[0], y[1]-y_des[1], y[2]-y_des[2], y[3]-y_des[3], y[4]-y_des[4], y[5]-y_des[5]};

	for (int i=0; i<12; i++){		//augmented states
		X_aug[i]= x[i]-x_prev[i];
	 };


	/*-----------------------------------------------define model matrices A, B, C-----------------------------------------*/
	float Am[12][12]=
	{
		{x[4]*x[1], x[4]*x[0]+x[5], 0, 1, x[0]*x[1], x[1], 0, 0, 0, 0, 0, 0 },
		{-x[5], 0, 0, 0, 1, -x[0], 0, 0, 0, 0, 0, 0},
		{x[4], 0, 0, 0, x[0], 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, (IY-IZ)/IX*x[5], (IY-IZ)/IX*x[4], 0, 0, 0, 0, 0, 0},
		{0, 0, 0, (IZ-IX)/IY*x[5], 0, (IZ-IX)/IY*x[3], 0, 0, 0, 0, 0, 0},
		{0, 0, 0, (IX-IY)/IZ*x[4], (IX-IY)/IZ*x[3], 0, 0, 0, 0, 0, 0, 0},
		{0, -G, 0, 0, -x[8], x[7], 0, x[5], -x[4], 0, 0, 0},
		{G, 0, 0, x[8], 0, -x[6], -x[5], 0, x[3], 0, 0, 0},
		{0, 0, 0, -x[7], x[6], 0, x[4], -x[3], 0, 0, 0, 0},
		{x[2]*x[8]+x[7]*x[1], x[8]+x[7]*x[0], x[8]*x[0]-x[7], 0, 0, 0, 1, x[0]*x[1]-x[2], x[0]*x[2]+x[1], 0, 0, 0},
		{x[7]*x[2]*x[1]-x[8], x[7]*x[0]*x[2]+x[8]*x[2], x[7]*x[0]*x[1]+x[8]*x[1]+x[6], 0, 0, 0, x[2], 1+x[0]*x[2]*x[1], x[2]*x[1]-x[0], 0, 0, 0},
		{x[7], -x[6], 0, 0, 0, 0, -x[1], x[0], 1, 0, 0, 0}

	};

	for (int c=0; c<12; c++){		//convert continuous A into discrete time A matrix
		for (int d=0; d<12; d++){
			if (c==d){
				Am[c][d]=1+Am[c][d]*delta_T;
			}else{
				Am[c][d]=Am[c][d]*delta_T;
			}
		}
	}
	// float Nm[18][18]=		//for testing only, remove in real usage
	// {
	// 	{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 1.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}

	// };

	// float Nm[2][2]=		//for testing only, remove in real usage
	// {
	// 	{1, 2},
	// 	{2, 4.001}

	// };

	// float testmatrix[4];		//testing code
	// int testn=12;
	// bool right=  mat_inverse(Nm, testmatrix, 2);
	//float (*testmatrix)[12]= matrix_power(Nm,20);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[0], testmatrix[1], testmatrix[2], testmatrix[3], testmatrix[4], testmatrix[5], testmatrix[6]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[12], testmatrix[13], testmatrix[14], testmatrix[15], testmatrix[16], testmatrix[17], testmatrix[18]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[24], testmatrix[25], testmatrix[26], testmatrix[27], testmatrix[28], testmatrix[29], testmatrix[30]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[36], testmatrix[3][1], testmatrix[3][2], testmatrix[3][3], testmatrix[3][4], testmatrix[3][5], testmatrix[3][6]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[48], testmatrix[4][1], testmatrix[4][2], testmatrix[4][3], testmatrix[4][4], testmatrix[4][5], testmatrix[4][6]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[60], testmatrix[5][1], testmatrix[5][2], testmatrix[5][3], testmatrix[5][4], testmatrix[5][5], testmatrix[5][6]);




	float Bm[12][4]=
	{
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 1/IX, 0, 0},
		{0, 0, 1/IY, 0},
		{0, 0, 0, 1/IZ},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{1/M, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0},
		{0, 0, 0, 0}};

	for (int c=0; c<12; c++){		//convert continuous B into discrete-time B matrix
		for (int d=0; d<4; d++){
			Bm[c][d]=Bm[c][d]*delta_T;
		}
	}

	float Cm[6][12]=
	{
		{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1} 

	};

	float CmBm[6][4];

	float sum=0;
	for(int c=0; c<6;c++){			//roll
	     for(int d=0; d<4; d++){	//column
	        for(int k=0; k<12; k++){
	            sum= sum + Cm[c][k] * Bm[k][d];
	        }
	        CmBm[c][d]= sum;
	        sum = 0;
	    }

	}

	float CmAm[6][12];
	sum=0;
	for(int c=0; c<6;c++){			//roll
	     for(int d=0; d<12; d++){	//column
	        for(int k=0; k<12; k++){
	            sum= sum + Cm[c][k] * Am[k][d];
	        }
	        CmAm[c][d]= sum;
	        sum = 0;
	    }

	}

	float A[18][18]; //augmented matrix A, 18x18 matrix

	for (int c=0; c<12; c++){
		for(int d=0; d<12; d++){
			A[c][d]=Am[c][d];
		}
	}

	for (int c=12; c<18; c++){
		for (int d=0; d<12; d++){
			A[c][d]=CmAm[c-12][d];
		}
	}

	for (int c=0; c<12; c++){
		for(int d=12; d<18; d++){
			A[c][d]=0;
		}
	}

	for (int c=12; c<18;c++){
		for (int d=12; d<18; d++){
			if(c==d){
				A[c][d]=1;
			}else{
				A[c][d]=0;
			}
		}
	}

	float B[18][4];	//augmented matrix B  18x4 matrix

	for (int c=0; c<12; c++){
		for (int d=0; d<4; d++){
			B[c][d]=Bm[c][d];
		}
	}

	for (int c=12; c<18; c++){
		for (int d=0; d<4; d++){
			B[c][d]=CmBm[c-12][d];
		}
	}

	float C[6][18];		//augmented matrix C, 6x18 matrix
	for (int c=0; c<6; c++){
		for (int d=0; d<12; d++){
			C[c][d]=0;
		}
	}

	for (int c=0; c<6; c++){
		for (int d=12; d<18; d++){

		}
	}

	// printf("%f %f %f %f %f %f %f  \n", A[0][0], A[0][1], A[0][2], A[0][3], A[0][4], A[0][5], A[0][6]);
	// printf("%f %f %f %f %f %f %f  \n", A[1][0], A[1][1], A[1][2], A[1][3], A[1][4], A[1][5], A[1][6]);
	// printf("%f %f %f %f %f %f %f  \n", A[2][0], A[2][1], A[2][2], A[2][3], A[2][4], A[2][5], A[2][6]);
	// printf("%f %f %f %f %f %f %f  \n", A[3][0], A[3][1], A[3][2], A[3][3], A[3][4], A[3][5], A[3][6]);
	// printf("%f %f %f %f %f %f %f  \n", A[4][0], A[4][1], A[4][2], A[4][3], A[4][4], A[4][5], A[4][6]);
	// printf("%f %f %f %f %f %f %f  \n", A[5][0], A[5][1], A[5][2], A[5][3], A[5][4], A[5][5], A[5][6]);

	// float (*testmatrix)[18];		//testing matrix power function
	// testmatrix= matrix_power(A,100);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[0][0], testmatrix[0][1], testmatrix[0][2], testmatrix[0][3], testmatrix[0][4], testmatrix[0][5], testmatrix[0][6]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[1][0], testmatrix[1][1], testmatrix[1][2], testmatrix[1][3], testmatrix[1][4], testmatrix[1][5], testmatrix[1][6]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[2][0], testmatrix[2][1], testmatrix[2][2], testmatrix[2][3], testmatrix[2][4], testmatrix[2][5], testmatrix[2][6]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[3][0], testmatrix[3][1], testmatrix[3][2], testmatrix[3][3], testmatrix[3][4], testmatrix[3][5], testmatrix[3][6]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[4][0], testmatrix[4][1], testmatrix[4][2], testmatrix[4][3], testmatrix[4][4], testmatrix[4][5], testmatrix[4][6]);
	// printf("%f %f %f %f %f %f %f  \n", testmatrix[5][0], testmatrix[5][1], testmatrix[5][2], testmatrix[5][3], testmatrix[5][4], testmatrix[5][5], testmatrix[5][6]);

/* --------------------------------------Laguarre FUnction formulation-----------------------------------------------------------*/

	float L[N][NP];
	L[0][0]=1;
	float beta = 1-ALPHA*ALPHA;
	float sqrt_beta = sqrt(beta);

	printf("%f %f\n",beta, sqrt_beta );
	for (int k=1; k<N; k++){		//formulate L[0]
		L[k][0]=L[k-1][0]*(-1)*ALPHA;
	}

	for (int k=0; k<N; k++){		//formulate L[0]
		L[k][0] *= sqrt_beta;
	}

	float Al[N][N];					//formulate Al
	for (int c=0; c<N; c++){
		for (int d=0; d<N; d++){
			if(c==d){
				Al[c][d]=ALPHA;
			}else if (c-d==1){
				Al[c][d]=beta;
			} else if (c-d>1){
				Al[c][d]=Al[c-1][d]*(-ALPHA);
			} else{
				Al[c][d]=0;
			}
		}
	}

	sum=0;
	for (int i = 1; i < NP; ++i)		//formulate L up to Np
	{
		for (int c=0; c<N; c++){
			for (int d=0; d<N; d++){
				sum = sum + Al[c][d]*L[d][i-1];
			}
			L[c][i] = sum;
			sum=0;
		}
	}

	// float T[18][18]=		//for testing only, remove in real usage
	// {
	// 	{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 1.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 1.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
	// 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}

	// };

/*--------------------------------objective function variables------------------------------------------------------------------------*/
	float phi[N*4][18];
	float phi_one[N*4][18];
	//float (*A_power)[18];
	float A_power[18][18];
	float Al_power[N][N];
	float inter_result2[N][N];	//for calculating Al^(m-1)
	float inter_result3[18][18];//for calculating A^m
	float Alp_phi_one[N*4][18];
	float omega[N*4][N*4];
	float psi[N*4][18];
	float Q[18][18] = 
	{
		{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}

	};

	float R[4]={1, 1, 1, 1};
	float R_L[N*4][N*4];
	for (int c=0; c<4*N; c++){
		for (int d=0; d<4*N; d++){
			R_L[c][d]=0;
		}
	}
	for (int n=0; n<4; n++){
		for (int c=0; c<N; c++){
			R_L[n*N+c][n*N+c]=R[n];
		}
	}
	//float LiBt[N][18];	//intermediate reuslt for finding phi
	// sum=0;
	// for (int c=0; c<N; c++){ 	//L[N][0]*B^T[][], Nx4 * 4x18
	// 	for (int d=0; d<18, d++){
	// 		for (int k=0; k<4; k++){
	// 			sum= sum+ L[c][0]*B[d][k];
	// 		}
	// 		LiBt[c][d]=sum;
	// 		sum=0;
	// 	}
	// }

	// get phi_one
	// for (int c= 0; c<18; c++){			//firstly get A^0, basically a unity matrix
	// 	for (int d=0; d<18; d++){
	// 		if (c==d){
	// 			A_power[c][d]=1.0;
	// 		} else{
	// 			A_power[c][d]=0;
	// 		}
	// 	}	
	// }

	for (int n=0; n<4; n++){		//number of input
		for (int c=0; c<N; c++){ 	//L[N][0]*B^T[][], Nx4 * 4x18
			for (int d=0; d<18; d++){
				phi_one[n*N+c][d]=L[c][0]*B[d][n];
				phi[n*N+c][d]=L[c][0]*B[d][n];
			}
		}	
	}
	// for (int c = 0; c < N; ++c){	//calculate phi=LiBt x A^T
	// 	for (int d = 0; d<18; ++d){
	// 		for (int k=0; k<18; k++){
	// 			sum = sum+LiBt[c][k]*A_power[d][k];
	// 		}
	// 		phi_one[c][d]= sum;	//add to the phi
	// 		sum=0;
	// 	}
	// }


	for (int m = 1; m <= NP; ++m)		//calculate phi(1) up to phi(NP). At each iteration m, we get phi(m) inside the 2D array phi[N][18]
	{		
		if (m==1){								//firstly define a unity matrix, as Al^0
				for (int c= 0; c<N; c++){	
					for (int d=0; d<N; d++){
						if (c==d){
							Al_power[c][d]=1.0;
						} else{
							Al_power[c][d]=0;
						}
					}	
				}

		}else{									//get Al^(m-1), by multiplying with Al in each iteration
			for (int c=0; c<N; c++){
				for (int d=0; d<N; d++){
					for (int k=0; k<N; k++){
						sum = sum + Al[c][k]*Al_power[k][d]; 
					}
					inter_result2[c][d]=sum;		//calculate each component in the matrix and save in intermediate result2
					sum=0;
				}
			}
			for (int c=0; c<N; c++){
				for (int d=0; d<N; d++){
					Al_power[c][d]=inter_result2[c][d];	//match  Al_power with the result2, this is for the next interation
				}
			}		
		}

		if (m==1){								//firstly define a unity matrix, as A^1
				for (int c= 0; c<18; c++){	
					for (int d=0; d<18; d++){
						A_power[c][d]=A[c][d];
					}	
				}

		}else{									//get A^(m) by multiplying with A in each iteration
			for (int c=0; c<18; c++){
				for (int d=0; d<18; d++){
					for (int k=0; k<18; k++){
						sum = sum + A[c][k]*A_power[k][d]; 
					}
					inter_result3[c][d]=sum;		//calculate each component in the matrix and save in intermediate result3
					sum=0;
				}
			}
			for (int c=0; c<18; c++){
				for (int d=0; d<18; d++){
					A_power[c][d]=inter_result3[c][d];	//match  A_power with the result3, this is for the next interation
				}
			}		
		}


		for (int n=0; n<4; n++){
			for (int c=0; c<N; c++){	//get Al^(m-1)*phi(1)
				for (int d=0; d<18; d++){
					for (int k=0; k<N; k++){
						sum = sum + Al_power[c][k]*phi_one[n*N+k][d];
					}
					Alp_phi_one[n*N+c][d]=sum;
					sum = 0;
				}
			}
		}
		float phi_At[N*4][18];
		for (int c = 0; c < 4*N; ++c){		//get phi^(m-1)*A^T
			for (int d=0; d<18; d++){
				for (int k=0; k<18; k++){
					sum = sum + phi[c][k]*A[d][k];
				}
				phi_At[c][d]=sum;
				sum=0;
			}
		}
		
		for (int c=0; c<4*N; ++c){
			for (int d=0; d<18; d++){
				phi[c][d] = phi_At[c][d]+ Alp_phi_one[c][d];	//adding to get phi(m)
			}
		}

		//calculate omega
		float phi_q[4*N][18];	//firstly calculate phi(m) x Q
		for (int c=0; c<4*N; c++){
			for (int d=0; d<18; d++){
				for (int k=0; k<18; k++){
					sum = sum+phi[c][k]*Q[k][d];
				}
				phi_q[c][d]=sum;
				sum=0;
			}
		}

		//float phi_q_phiT[4*N][4*N];
		for (int c=0; c<4*N; c++){		// then get phi(m) x Q x phi(m)^T
			for (int d=0; d<4*N; d++){
				for (int k=0; k<18; k++){
					sum = sum+phi_q[c][k]*phi[d][k];
				}
				//phi_q_phiT[c][d]=sum;
				omega[c][d]= omega[c][d] + sum;			//add to omega[c][d]
				sum=0;
			}
		}

		for (int c=0; c<4*N; c++){
			for (int d=0; d<18; d++){		//get phi(m) x Q x A^m
				for (int k=0; k<18; k++){
					sum = sum+ phi_q[c][k]*A_power[k][d];
				}
				psi[c][d]=psi[c][d]+sum;			//add to psi
				sum=0;
			}
		}

	}

	// for (int m = 1; m <= NP; ++m)		//calculate phi(1) up to phi(NP). At each iteration m, we get phi(m) inside the 2D array phi[N][18]
	// {		
	// 	// for (int i = 0; i < m; ++i)
	// 	// {
	// 	// 	A_power = matrix_power(A,m-i-1);
	// 	// }

	// 	sum=0;
	// 	for (int c=0; c<N; c++){ 	//L[N][i]*B^T[][], Nx4 * 4x18
	// 		for (int d=0; d<18; d++){
	// 			for (int k=0; k<4; k++){
	// 				sum= sum+ L[c][NP-m]*B[d][k];
	// 			}
	// 			LiBt[c][d]=sum;
	// 			sum=0;
	// 		}
	// 	}	

	// 	for (int c = 0; c < N; ++c){	//calculate phi=LiBt x A^T
	// 		for (int d = 0; d<18; ++d){
	// 			for (int k=0; k<18; k++){
	// 				sum = sum+LiBt[c][k]*A_power[d][k];
	// 			}
	// 			phi[c][d]= phi[c][d] + sum;	//add to the phi
	// 			sum=0;
	// 		}
	// 	}
	// }

	for (int c=0; c<4*N; c++){	//omega = omega +R_L
		for (int d=0; d<4*N; d++){
			omega[c][d]=omega[c][d]+R_L[c][d];
		}
	}

	float omega_one_D[16*N*N];			//convert 2d matrix to 1d
	for (int c = 0; c < 4*N; c++){
		for (int d = 0; d < 4*N; d++){
		    omega_one_D[c * 4*N + d] = omega[c][d];
		}
	}
	/*----------------------------------------------calculating change in u---------------------------------------------------------*/
	float omega_one_D_inv[16*N*N];		//testing code
	bool inverse_ok=  mat_inverse(omega_one_D, omega_one_D_inv, 4*N);
	float omega_inv[4*N][4*N];
	float minus_omega_inv_psi[4*N][18];
	float eta[4*N];

	if (inverse_ok){
		for (int c = 0; c < 4*N; c++){	//convert 1d matrix to 2d
			for (int d = 0; d < 4*N; d++){
			    omega_inv[c][d] = omega_one_D_inv[c * 4*N + d];
			}
		}

		for (int c=0; c<4*N; c++){			//get -omega^(-1) x psi
			for (int d=0; d<18; d++){
				for (int k=0; k<4*N; k++){
					sum = sum + omega_inv[c][k]*psi[k][d];
				}
				minus_omega_inv_psi[c][d] = -(sum);
				sum=0;
			}
		}

		for (int c=0; c<4*N; c++){
			for (int k=0; k<18; k++){
				sum = sum + minus_omega_inv_psi[c][k]*X_aug[k];
			}
			eta[c]=sum;
			sum=0;
		}

//		for (int i=0; i<50; i++){
			for (int n=0; n<4; n++){
				for (int c=0; c<N; c++){
					sum=sum+ L[c][0]*eta[n*N+c];
				}
				delta_u[n]=sum;
				sum=0;
			}
//			printf("%f %f %f %f \n", delta_u[0], delta_u[1], delta_u[2], delta_u[3]);
//		}

	}else{		//if inverse failed
		printf("inverse failed!!!!!\n");
	}

	/*-------------------------------------------------forward simulation---------------------------------------*/



	// if (inverse_ok){
	// 	printf("%f %f %f %f %f %f %f  \n", omega_inv[0][0], omega_inv[0][1], omega_inv[0][2], omega_inv[0][3], omega_inv[0][4], omega_inv[0][5], omega_inv[0][6]);
	// 	printf("%f %f %f %f %f %f %f  \n", omega_inv[1][0], omega_inv[1][1], omega_inv[1][2], omega_inv[1][3], omega_inv[1][4], omega_inv[1][5], omega_inv[1][6]);
	// 	printf("%f %f %f %f %f %f %f  \n", omega_inv[2][0], omega_inv[2][1], omega_inv[2][2], omega_inv[2][3], omega_inv[2][4], omega_inv[2][5], omega_inv[2][6]);
	// 	printf("%f %f %f %f %f %f %f  \n", omega_inv[3][0], omega_inv[3][1], omega_inv[3][2], omega_inv[3][3], omega_inv[3][4], omega_inv[3][5], omega_inv[3][6]);
	// 	printf("%f %f %f %f %f %f %f  \n", omega_inv[4][0], omega_inv[4][1], omega_inv[4][2], omega_inv[4][3], omega_inv[4][4], omega_inv[4][5], omega_inv[4][6]);
	// 	printf("%f %f %f %f %f %f %f  \n", omega_inv[5][0], omega_inv[5][1], omega_inv[5][2], omega_inv[5][3], omega_inv[5][4], omega_inv[5][5], omega_inv[5][6]);
	// }
	//delete[] A_power;
}

	float (*matrix_power(float Ma[][18], int n))[18]{
		// int rows =  sizeof A / sizeof A[0];
		// int cols = sizeof A[0] / sizeof(A[0][0]);
		float inter_result[18][18];
		float inter_result2[18][18];
		static float result[18][18];

		for (int c= 0; c<18; c++){	//firstly define a unity matrix
			for (int d=0; d<18; d++){
				if (c==d){
					inter_result[c][d]=1.0;
				} else{
					inter_result[c][d]=0;
				}
			}
		}

		float sum=0;
		for (int power=0; power<n; power++){
			for (int c=0; c<18; c++){
				for (int d=0; d<18; d++){
					for (int k=0; k<18; k++){
						sum = sum + Ma[c][k]*inter_result[k][d]; 
					}
					inter_result2[c][d]=sum;		//calculate each component in the matrix and save in intermediate result2
					sum=0;
				}
			}
			for (int c=0; c<18; c++){
				for (int d=0; d<18; d++){
					inter_result[c][d]=inter_result2[c][d];	//match the intermediate result 1 with the result2, this is for the next interation
				}
			}
		}

		for (int c=0; c<18; c++){					//match the result for output with intermediate result
			for (int d=0; d<18; d++){
				result[c][d]=inter_result2[c][d];
			}
		}
		return result;
	}

