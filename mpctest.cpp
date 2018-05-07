#include <stdio.h>
#include <math.h>
/*
mpc test program
states:
input x: [phi, theta, psi, p, q, r, u, v, w, x, y, z]
output:	[phi, theta, psi, x, y, z]

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

}