#ifndef SMOKE_H
#define SMOKE_H

#include <vector>
#include "CSCIx229.h"
#include <stdlib.h>

using std::vector;

//variables

extern int N; //size in one direction
extern int size; //total particles
//velocities (v1 = x, v2 = y, v3 = z)
extern vector<float> v1;
extern vector<float> v2;
extern vector<float> v3;
//previous velocities
extern vector<float> v1_prev;
extern vector<float> v2_prev;
extern vector<float> v3_prev;
//densities
extern vector<float> dens;
extern vector<float> dens_prev;
//other constants
extern float visc;
extern float dt;
extern float diff;
//Sets k, the iteration count for the iterative solvers
extern int K;

//TODO Source Vector

//Gives index into vectors from x, y, z coordinates
int IX(int x, int y, int z);

//Diffuse
void diffuse (int b, vector<float> &x, vector<float> &x0, float diff, float dt);

//Advect
void advect (int b, vector<float> &d, vector<float> &d0, vector<float> &v1, vector<float> &v2, vector<float> &v3, float dt);

//Density Step
void dens_step ( vector <float> &x, vector<float> &x0, vector<float> &v1, vector<float> &v2, vector<float> &v3, float diff, float dt );

//Velocity Steps
void vel_step (vector<float> &v1, vector<float> &v2, vector<float> &v3, vector<float> &v1_prev, vector<float> &v2_prev, vector<float> &v3_prev, float visc, float dt);

//Project, ensure boundary conditions
void projectNS (vector<float> &v1, vector<float> &v2, vector<float> &v3, vector<float> &p, vector<float> &div);

//Set bounds
void set_bnd (int b, vector<float> &x);

//Draw
void draw_smoke();

//Run progression and call draw_smoke()
void run_smoke(bool source);

//setup
void setup();

#endif /* SMOKE_H */
