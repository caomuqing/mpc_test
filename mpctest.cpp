#include <stdio.h>
#include <math.h>
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

#define IX 1.0 //moment of inertia X axis
#define IY 1.0 //moment of inertia X axis
#define IZ 1.0 //moment of inertia X axis
#define M 1.5 //mass of copter in kg
#define G 9.81//gravity constant

main() {

	printf("Hello mpc\n");
	//getting model states
	float x[12] =
	{
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
	};

	float u[4] = {M*G, 0, 0, 0};

	//define model matrices A, B, C
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

	float A[18][18]; //augmented matrix A

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

	float B[18][4];	//augmented matrix B

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

	float C[6][18];		//augmented matrix C
	for (int c=0; c<6; c++){
		for (int d=0; d<12; d++){
			C[c][d]=0;
		}
	}

	for (int c=0; c<6; c++){
		for (int d=12; d<18; d++){

		}
	}



}

	float (*matrix_power(float A[][12], int n))[12]{
		int rows =  sizeof A / sizeof A[0];
		int cols = sizeof A[0] / sizeof(A[0][0]);
		static float B[12][12];

		for (int c= 0; c<rows; c++){
			for (int d=0; d<cols; d++){
				B[c][d]=A[c][d];
			}
		}
		return B;
	}