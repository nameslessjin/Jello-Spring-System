/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include <iostream>
#include "jello.h"
#include "physics.h"

using namespace std;

#define gravity 9.8;
#define r_length 1.0/7.0;

double get_length(struct point vector) {
  return sqrt(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2));
}

double dot(struct point vector1, struct point vector2) {
  return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}

void normalize(struct point &vector) {
  double length = get_length(vector);
  vector.x /= length;
  vector.y /= length;
  vector.z /= length;
}

struct point computeHookForce(double kh, struct point p1, struct point p2, double resting_length) {

  struct point hook_force;
  
  // L is the vector pointing from p2 to p1
  struct point L, L_normalized;
  pDIFFERENCE(p1, p2, L);
  pCPY(L, L_normalized);
  normalize(L_normalized);

  // F = -k_hook(|L| - R) * (L/|L|)
  double pre = -1 * kh * (get_length(L) - resting_length);
  pMULTIPLY(L_normalized, pre, hook_force);

  return hook_force;
}

struct point computeDampingForce(double kd, struct point p1, struct point p2) {
  struct point damping_force;
  
  // L is the vector pointing from B to A
  struct point L, L_normalized, p1_p2;
  pDIFFERENCE(p1, p2, L);
  pCPY(L, L_normalized);
  normalize(L_normalized);

  // F = -k_damping * ((v_a - v_b) dot L_normalized) * L_normalized
  double pre = -1 * kd;
  pDIFFERENCE(p2, p1, p1_p2);
  pre *= dot(p1_p2, L_normalized);
  pMULTIPLY(L_normalized, pre, damping_force);

  return damping_force;
}

void computeStructureForce() {}
void computeShearForce() {}
void computeBendForce() {}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  /* for you to implement ... */

  // gravity
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int k = 0; k < 8; ++k) {
        a[i][j][k].x = 0;
        a[i][j][k].y = 0;
        a[i][j][k].z = 0;

        // Hook for all springs except collision springs

      }
    }
  }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F1p = v * dt
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);

         // F1v = a * dt
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);

         // buffer.p = F1p * 0.5
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);

         // buffer.v = F1v * 0.5
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);

         // buffer.p = p + buffer.p
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);

         // buffer.v = v + buffer.v
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);

         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);


         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
