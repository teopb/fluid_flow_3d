#include "smoke.h"
#include <stdio.h>
#include <iostream>

using std::vector;

//variables

int N; //size in one direction
int size; //total particles

//velocities (v1 = x, v2 = y, v3 = z)
vector<float> v1;
vector<float> v2;
vector<float> v3;
//previous velocities
vector<float> v1_prev;
vector<float> v2_prev;
vector<float> v3_prev;
//densities
vector<float> dens;
vector<float> dens_prev;
//other constants
float visc;
float dt;
float diff;
//Sets k, the iteration count for the iterative solvers
int K1;
int K2;

int sourceCount = 0;

//Routine for printing vectors, used for debugging
void printV(vector<float> &x){
  for (size_t i = 1; i <= N; i++) {
    for (size_t j = 1; j <= N; j++) {
      for (size_t l = 1; l <= N; l++) {
        std::cout << x.at(IX(i, j, l))<< ", ";
      }
    printf("\n");
    }
  printf("new slice\n");
  }
}
//Gives index into vectors from x, y, z coordinates
int IX(int x, int y, int z){
  //N+2 for boundaries
  return x + (N+2)*y +(N+2)*(N+2)*z;
}

//Diffuse
void diffuse (int b, vector<float> &x, vector<float> &x0, float diff, float dt){
  //This scales the difference to be relative to the number of cells and dt
  float a = dt * diff * N * N * N;

  for (size_t k = 0; k < K1; k++) {
    for (size_t i = 1; i <= N; i++) {
      for (size_t j = 1; j <= N; j++) {
        for (size_t l = 1; l <= N; l++) {
          x.at(IX(i, j, l)) = (x0.at(IX(i, j, l)) + a*(x.at(IX(i-1, j, l)) + x.at(IX(i+1, j, l)) + x.at(IX(i, j-1, l)) + x.at(IX(i, j+1, l))+ x.at(IX(i, j, l-1)) + x.at(IX(i, j, l+1))))/(1+6*a);
        }
      }
    }
    set_bnd(b, x);
  }
  return;
}

//Advect
void advect (int b, vector<float> &d, vector<float> &d0, vector<float> &v1, vector<float> &v2, vector<float> &v3, float dt){
  int i0, j0, l0, i1, j1, l1;
  float x, y, z, x0, x1, y0, y1, z0, z1;
  float dt_scaled = dt*N;

  for (size_t i = 1; i <= N; i++) {
    for (size_t j = 1; j <= N; j++) {
      for (size_t l = 1; l <= N; l++) {
        try
        {
          //backsolves for location of particle that ended up at final location
          x = i - dt_scaled * v1.at(IX(i, j, l));
          y = j - dt_scaled * v2.at(IX(i, j, l));
          z = l - dt_scaled * v3.at(IX(i, j, l));

          //enforces boundary conditions
          if (x < 0.5) x = 0.5;
          if (x > N + 0.5) x = N+ 0.5;

          //set i0, i1 to  integer boundaries around x
          i0 = static_cast<int>(x);
          i1 = i0+1;

          //Same for y and z
          if (y < 0.5) y = 0.5;
          if (y > N + 0.5) y = N+ 0.5;
          j0 = static_cast<int>(y);
          j1 = j0+1;

          if (z < 0.5) z = 0.5;
          if (z > N + 0.5) z = N+ 0.5;
          l0 = static_cast<int>(z);
          l1 = l0+1;

          //get relative positions in cell
          x1 = x - i0;
          x0 = 1 - x1;

          y1 = y - j0;
          y0 = 1 - y1;

          z1 = z - l0;
          z0 = 1 - z1;


          //Using relative positions set new density as a weighted combination of the contributing cells
          d.at(IX(i, j, l)) =
          x0 * (y0 * (z0 * d0.at(IX(i0, j0, l0)) + z1 * d0.at((IX(i0, j0, l1))))) +
          x0 * (y1 * (z0 * d0.at(IX(i0, j1, l0)) + z1 * d0.at((IX(i0, j1, l1))))) +
          x1 * (y0 * (z0 * d0.at(IX(i1, j0, l0)) + z1 * d0.at((IX(i1, j0, l1))))) +
          x1 * (y1 * (z0 * d0.at(IX(i1, j1, l0)) + z1 * d0.at((IX(i1, j1, l1)))));
        }
        catch(...)
        {
          printf("error i0 = %d, j0 = %d, l0 = %d\n", i0, j0, l0);
          printf("error i1 = %d, j1 = %d, l1 = %d\n", i1, j1, l1);
          printf("dt_scaled = %f\n", dt_scaled);
          printf("x = %f, y = %f, z= %f\n", x, y, z);
          printf("v1 = %f, v2 = %f, v3= %f\n", v1.at(IX(i, j, l)), v2.at(IX(i, j, l)), v3.at(IX(i, j, l)));

          exit(1);
        }
      }
    }
  }

  set_bnd(b, d);

  return;
}

//Density Step
void dens_step ( vector <float> &x, vector<float> &x0, vector<float> &v1, vector<float> &v2, vector<float> &v3, float diff, float dt ){

  diffuse(0, x, x0, diff, dt);
  x.swap(x0);

  advect(0, x, x0, v1, v2, v3, dt);
  x0.swap(x);

  return;
}

//Velocity Steps
void vel_step (vector<float> &v1, vector<float> &v2, vector<float> &v3, vector<float> &v1_prev, vector<float> &v2_prev, vector<float> &v3_prev, float visc, float dt){

  v1.swap(v1_prev);
  diffuse(1, v1, v1_prev, visc, dt);

  v2.swap(v2_prev);
  diffuse(2, v2, v2_prev, visc, dt);

  v3.swap(v3_prev);
  diffuse(3, v3, v3_prev, visc, dt);

  projectNS(v1, v2, v3, v3_prev, v2_prev);


  v1.swap(v1_prev);
  v2.swap(v2_prev);
  v3.swap(v3_prev);

  advect(1, v1, v1_prev, v1_prev, v2_prev, v3_prev, dt);
  advect(2, v2, v2_prev, v1_prev, v2_prev, v3_prev, dt);
  advect(3, v3, v3_prev, v1_prev, v2_prev, v3_prev, dt);

  projectNS(v1, v2, v3, v3_prev, v2_prev);

  return;
}

//Project, ensure boundary conditions, helps prevent blowup
void projectNS (vector<float> &v1, vector<float> &v2, vector<float> &v3, vector<float> &p, vector<float> &div){
  float h = 1.0/N;


  for (size_t i = 1; i <= N; i++) {
    for (size_t j = 1; j <= N; j++) {
      for (size_t l = 1; l <= N; l++) {
        div.at(IX(i, j, l)) = -0.5 * h * (v1.at(IX(i+1, j, l)) - v1.at(IX(i-1, j, l))
        + v2.at(IX(i, j+1, l)) - v2.at(IX(i, j-1, l))
        + v3.at(IX(i, j, l+1)) - v3.at(IX(i, j, l-1)));
        p.at(IX(i, j, l)) = 0;
      }
    }
  }


  set_bnd(0, div);
  set_bnd(0, p);

  for (size_t k = 0; k < K2; k++) {
    for (size_t i = 1; i <= N; i++) {
      for (size_t j = 1; j <= N; j++) {
        for (size_t l = 1; l <= N; l++) {
          p.at(IX(i, j, l)) = (div.at(IX(i, j, l)) + p.at(IX(i+1, j, l)) + p.at(IX(i-1, j, l))
          + p.at(IX(i, j+1, l)) + p.at(IX(i, j-1, l))
          + p.at(IX(i, j, l+1)) + p.at(IX(i, j, l-1)))/6;
        }
      }
    }
    set_bnd(0, p);
  }


  for (size_t i = 1; i <= N; i++) {
    for (size_t j = 1; j <= N; j++) {
      for (size_t l = 1; l <= N; l++) {

        v1.at(IX(i, j, l)) -= 0.5 * (p.at(IX(i+1, j, l)) - p.at(IX(i-1, j, l)))/h;
        v2.at(IX(i, j, l)) -= 0.5 * (p.at(IX(i, j+1, l)) - p.at(IX(i, j-1, l)))/h;
        v3.at(IX(i, j, l)) -= 0.5 * (p.at(IX(i, j, l+1)) - p.at(IX(i, j, l-1)))/h;
      }
    }
  }
  set_bnd(1, v1);
  set_bnd(2, v2);
  set_bnd(3, v3);


  return;
}

//Set bounds
void set_bnd (int b, vector<float> &x){
  for (size_t i = 1; i <= N; i++) {
    for (size_t j = 1; j <= N; j++) {
      if (b == 1) {
        x.at(IX(0, i, j)) = -x.at(IX(1, i, j));
        x.at(IX(N+1, i, j)) = -x.at(IX(N, i, j));
      }
      else{
        x.at(IX(0, i, j)) = x.at(IX(1, i, j));
        x.at(IX(N+1, i, j)) = x.at(IX(N, i, j));
      }
      if (b == 2) {
        x.at(IX(i, 0, j)) = -x.at(IX(i, 1, j));
        x.at(IX(i, N+1, j)) = -x.at(IX(i, N, j));
      }
      else{
        x.at(IX(i, 0, j)) = x.at(IX(i, 1, j));
        x.at(IX(i, N+1, j)) = x.at(IX(i, N, j));
      }
      if (b == 3) {
        x.at(IX(i, j, 0)) = -x.at(IX(i, j, 1));
        x.at(IX(i, j, N+1)) = -x.at(IX(i, j, N));
      }
      else{
        x.at(IX(i, j, 0)) = x.at(IX(i, j, 1));
        x.at(IX(i, j, N+1)) = x.at(IX(i, j, N));
      }
    }
  }
  //for corners
  x.at(IX(0, 0, 0)) = 0.3333*(x.at(IX(1, 0, 0)) + x.at(IX(0, 1, 0)) + x.at(IX(0, 0, 1)));
  x.at(IX(0, N+1, 0)) = 0.3333*(x.at(IX(1, N+1, 0)) + x.at(IX(0, N, 0)) + x.at(IX(0, N+1, 1)));
  x.at(IX(N+1, 0, 0)) = 0.3333*(x.at(IX(N, 0, 0)) + x.at(IX(N+1, 1, 0)) + x.at(IX(N+1, 0, 1)));
  x.at(IX(N+1, N+1, 0)) = 0.3333*(x.at(IX(N, N+1, 0)) + x.at(IX(N+1, N, 0)) + x.at(IX(N+1, N+1, 1)));
  x.at(IX(0, 0, N+1)) = 0.3333*(x.at(IX(1, 0, N+1)) + x.at(IX(0, 1, N+1)) + x.at(IX(0, 0, N)));
  x.at(IX(0, N+1, N+1)) = 0.3333*(x.at(IX(1, N+1, N+1)) + x.at(IX(0, N, N+1)) + x.at(IX(0, N+1, N)));
  x.at(IX(N+1, 0, N+1)) = 0.3333*(x.at(IX(N, 0, N+1)) + x.at(IX(N+1, 1, N+1)) + x.at(IX(N+1, 0, N)));
  x.at(IX(N+1, N+1, N+1)) = 0.3333*(x.at(IX(N, N+1, N+1)) + x.at(IX(N+1, N, N+1)) + x.at(IX(N+1, N+1, N)));
  return;
}

//Draw
void draw_smoke(){
  //condense the grid cells so the location = (x,y,z)*scale
  float scale = .025;
  float x, y, z;

  glPointSize(12);
  glBegin(GL_POINTS);
  float dtotal = 0;

  for (size_t i = 1; i <= N; i++) {
    for (size_t j = 1; j <= N; j++) {
      for (size_t l = 1; l <= N; l++) {

        dtotal += dens.at(IX(i, j, l));
        if (dens.at(IX(i, j, l))>.3) {

          //base color
          float color = .8;

          //set color and alpha value based on density
          glColor4f(color, color, color, dens.at(IX(i, j, l))*.003);
          x = i * scale;
          y = j * scale;
          z = l * scale;

          //Mirror vertices around y axis this could maybe
          //involve some randomness to improve the look
          glVertex3d(x, y, z);
          glVertex3d(-x, y, z);
          glVertex3d(-x, y, -z);
          glVertex3d(x, y, -z);
        }
      }
    }
  }
  glEnd();
  return;
}

//Run progression and call draw_smoke()
void run_smoke(bool source){
  int loc = 2;
  //steady source
  if(source){
    dens_prev.at(IX(loc, loc, loc)) += 1900;
    v2_prev.at(IX(loc, loc, loc)) += 1700000;
  }

  vel_step(v1, v2, v3, v1_prev, v2_prev, v3_prev, visc, dt);
  dens_step(dens, dens_prev, v1, v2, v3, diff, dt);

  return;
}

//setup
void setup(){
  N = 60;
  size = (N+2)*(N*2)*(N+2);

  v1.resize(size, 0);
  v2.resize(size, 0);
  v3.resize(size, 0);

  v1_prev.resize(size, 0);
  v2_prev.resize(size, 0);
  v3_prev.resize(size, 0);

  dens.resize(size, 0);
  dens_prev.resize(size, 0);

  visc = .01;
  dt = .01;
  diff = .01;

  K1 = 5;
  K2 = 10;

  return;
}
