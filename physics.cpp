/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include <iostream>
#include "jello.h"
#include "physics.h"

using namespace std;

#define gravity 9.8
#define r_length 1.0 / 7.0

double getLength(struct point vector)
{
  return sqrt(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2));
}

double dot(struct point vector1, struct point vector2)
{
  return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}

void normalize(struct point &vector)
{
  double length = getLength(vector);
  vector.x /= length;
  vector.y /= length;
  vector.z /= length;
}

struct point find_v(struct world *jello, int x, int y, int z)
{
  return jello->v[x][y][z];
}

bool isValidVertex(int x, int y, int z) {
  if (x < 0 || x > 7 || y < 0 || y > 7 || z < 0 || z > 7) return false;
  return true;
}

bool isDirectNeighbor(int x, int y, int z, int x1, int y1, int z1) {

  // if it is not one of the 6 direct neighbors then it is diagonal neighbors
  // to be a direct neighbors only one of the 3 scalar is different by 1
  if ((x1 == x + 1 || x1 == x - 1) && y == y1 && z == z1) return true;
  if ((y1 == y + 1 || y1 == y - 1) && x == x1 && z == z1) return true;
  if ((z1 == z + 1 || z1 == z - 1) && x == x1 && y == y1) return true;

  return false;
}

bool isEveryOtherNode(int x, int y, int z, int x1, int y1, int z1) {

  // if it is not one of the 6 direct every other neighbors then it is false
  // to be a direct neighbors only one of the 3 scalar is different by 2
  if ((x1 == x + 2 || x1 == x - 2) && y == y1 && z == z1) return true;
  if ((y1 == y + 2 || y1 == y - 2) && x == x1 && z == z1) return true;
  if ((z1 == z + 2 || z1 == z - 2) && x == x1 && y == y1) return true;

  return false;
}

bool isSelf(int x, int y, int z, int x1, int y1, int z1) {
  return x == x1 && y == y1 && z == z1;
}


struct point computeHookForce(double kh, struct point a, struct point b, double resting_length)
{

  // Find the elastic force exterted on A

  struct point hook_force;

  // L is the vector pointing from b to a
  struct point L, L_normalized;
  pDIFFERENCE(a, b, L);
  pCPY(L, L_normalized);
  normalize(L_normalized);

  // F = -k_hook(|L| - R) * (L/|L|)
  double pre = -1 * kh * (getLength(L) - resting_length);
  pMULTIPLY(L_normalized, pre, hook_force);

  return hook_force;
}

struct point computeDampingForce(double kd, struct point a, struct point b, struct point va, struct point vb)
{
  struct point damping_force;

  // L is the vector pointing from B to A
  struct point L, L_normalized, va_vb;
  pDIFFERENCE(a, b, L);
  pCPY(L, L_normalized);
  normalize(L_normalized);

  // F = -k_damping * ((v_a - v_b) dot L_normalized) * L_normalized
  double pre = -1 * kd;
  pDIFFERENCE(va, vb, va_vb);
  pre *= dot(va_vb, L_normalized);
  pMULTIPLY(L_normalized, pre, damping_force);

  return damping_force;
}

void computeStructureForce(struct world *jello, int x, int y, int z, struct point &F)
{

  // compute the structure force at a single point
  // there are 4 cases
  // first, the very corner case, there are only 8 of theses cases, p only has 3 connections
  // second, the edge case, where the point has 4 connections
  // third, the face case, where the point has 5 connections
  // fourth, the center case, where the point has 6 connections

  struct point p, hook_force, damp_force;
  pCPY(jello->p[x][y][z], p);

  // loop through all neighbors and find valid direct neighbors
  for (int i = x - 1; i <= x + 1; ++i) {
    for (int j = y - 1; j <= y + 1; ++j) {
      for (int k = z - 1; k <= z + 1; ++k) {
        
        if (isValidVertex(i, j, k) && isDirectNeighbor(x, y, z, i, j, k) && !isSelf(x, y, z, i, j, k)) {
          struct point tmp;
          pCPY(jello->p[i][j][k], tmp);
          hook_force = computeHookForce(jello->kElastic, p, tmp, r_length);
          damp_force = computeDampingForce(jello->dElastic, p, tmp, find_v(jello, x, y, z), find_v(jello, i, j, k));
          pSUM(F, hook_force, F);
          pSUM(F, damp_force, F);
        }

      }
    }
  }
}


void computeShearForce(struct world *jello, int x, int y, int z, struct point &F) {

  // compute the shear force at a single point
  // there are 4 cases
  // first, the very corner case, there are only 8 of theses cases, p only has 4 connections
  // second, the edge case, where the point has 7 connections
  // third, the face case, where the point has 12 connections
  // fourth, the center case, where the point has 20 connections

  struct point p, hook_force, damp_force;
  pCPY(jello->p[x][y][z], p);

  // loop through all neighbors and find valid indirect neighbors
  for (int i = x - 1; i <= x + 1; ++i) {
    for (int j = y - 1; j <= y + 1; ++j) {
      for (int k = z - 1; k <= z + 1; ++k) {
        
        if (isValidVertex(i, j, k) && !isSelf(x, y, z, i, j, k) && !isDirectNeighbor(x, y, z, i, j, k)) {
          struct point tmp;
          pCPY(jello->p[i][j][k], tmp);
          hook_force = computeHookForce(jello->kElastic, p, tmp, r_length);
          damp_force = computeDampingForce(jello->dElastic, p, tmp, find_v(jello, x, y, z), find_v(jello, i, j, k));
          pSUM(F, hook_force, F);
          pSUM(F, damp_force, F);
        }
      }
    }
  }
}

void computeBendForce(struct world *jello, int x, int y, int z, struct point &F) {

  // similar to struct force but this is for every other point and current point

  struct point p, hook_force, damp_force;
  pCPY(jello->p[x][y][z], p);

  // loop through one beyond neighbors and skip neighbors
  for (int i = x - 2; i <= x + 2; i += 2) {
    for (int j = y - 2; j <= y + 2; j += 2) {
      for (int k = z - 2; k <= z + 2; k += 2) {
        if (isValidVertex(i, j, k) && !isSelf(x, y, z, i, j, k) && isEveryOtherNode(x, y, z, i, j, k)) {
          struct point tmp;
          pCPY(jello->p[i][j][k], tmp);
          hook_force = computeHookForce(jello->kElastic, p, tmp, r_length);
          damp_force = computeDampingForce(jello->dElastic, p, tmp, find_v(jello, x, y, z), find_v(jello, i, j, k));
          pSUM(F, hook_force, F);
          pSUM(F, damp_force, F);
        }
      }
    }
  }
}

/* Computes acceleration to every control point of the jello cube,
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world *jello, struct point a[8][8][8])
{
  /* for you to implement ... */

  // gravity
  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j < 8; ++j)
    {
      for (int k = 0; k < 8; ++k)
      {
        a[i][j][k].x = 0;
        a[i][j][k].y = 0;
        a[i][j][k].z = 0;

        struct point force;

        // compute force for cube itself
        computeStructureForce(jello, i, j, k, force);
        computeShearForce(jello, i, j, k, force);
        computeBendForce(jello, i, j, k, force);


        // Hook for all springs except collision springs
      }
    }
  }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world *jello)
{
  int i, j, k;
  point a[8][8][8];

  computeAcceleration(jello, a);

  for (i = 0; i <= 7; i++)
    for (j = 0; j <= 7; j++)
      for (k = 0; k <= 7; k++)
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
void RK4(struct world *jello)
{

  // the jello is 8x8x8
  // each point has x, y, z

  point F1p[8][8][8], F1v[8][8][8],
      F2p[8][8][8], F2v[8][8][8],
      F3p[8][8][8], F3v[8][8][8],
      F4p[8][8][8], F4v[8][8][8];

  // "a" here stands for acceleration on each point
  point a[8][8][8];

  struct world buffer;

  int i, j, k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  // For F1p and F1v only
  for (i = 0; i <= 7; i++)
    for (j = 0; j <= 7; j++)
      for (k = 0; k <= 7; k++)
      {
        // F1p = v * dt
        pMULTIPLY(jello->v[i][j][k], jello->dt, F1p[i][j][k]);

        // F1v = a * dt
        pMULTIPLY(a[i][j][k], jello->dt, F1v[i][j][k]);

        // buffer.p = F1p * 0.5
        pMULTIPLY(F1p[i][j][k], 0.5, buffer.p[i][j][k]);

        // buffer.v = F1v * 0.5
        pMULTIPLY(F1v[i][j][k], 0.5, buffer.v[i][j][k]);

        // buffer.p = jello->p + buffer.p
        pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);

        // buffer.v = jello->v + buffer.v
        pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
      }

  // recompute acceleration with the update bufer value
  computeAcceleration(&buffer, a);

  // For F2p and F2v only
  // Buffer and acceleration are different from prev loops
  for (i = 0; i <= 7; i++)
    for (j = 0; j <= 7; j++)
      for (k = 0; k <= 7; k++)
      {
        // F2p = buffer.v * dt;
        pMULTIPLY(buffer.v[i][j][k], jello->dt, F2p[i][j][k]);

        // F2v = a * dt;
        pMULTIPLY(a[i][j][k], jello->dt, F2v[i][j][k]);

        // buffer.p = F2p * 0.5;
        pMULTIPLY(F2p[i][j][k], 0.5, buffer.p[i][j][k]);

        // buffer.v = F2v * 0.5;
        pMULTIPLY(F2v[i][j][k], 0.5, buffer.v[i][j][k]);

        // buffer.p = jello->p + buffer.p
        pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);

        // buffer.v = jello->v + buffer.v
        pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
      }

  // recompute acceleration with the update bufer value
  computeAcceleration(&buffer, a);

  // For F3p and F3v only
  // Buffer and acceleration are different from prev loops
  for (i = 0; i <= 7; i++)
    for (j = 0; j <= 7; j++)
      for (k = 0; k <= 7; k++)
      {
        // F3p = buffer.v * dt;
        pMULTIPLY(buffer.v[i][j][k], jello->dt, F3p[i][j][k]);

        // F3v = a * dt;
        pMULTIPLY(a[i][j][k], jello->dt, F3v[i][j][k]);

        // buffer.p = F3p * 1.0;
        pMULTIPLY(F3p[i][j][k], 1.0, buffer.p[i][j][k]);

        // buffer.v = F3v * 1.0;
        pMULTIPLY(F3v[i][j][k], 1.0, buffer.v[i][j][k]);

        // buffer.p = jello->p + buffer.p
        pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);

        // buffer.v = jello->v + buffer.v
        pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
      }

  // recompute acceleration with the update bufer value
  computeAcceleration(&buffer, a);

  for (i = 0; i <= 7; i++)
    for (j = 0; j <= 7; j++)
      for (k = 0; k <= 7; k++)
      {
        // F4p = buffer.v * dt;
        pMULTIPLY(buffer.v[i][j][k], jello->dt, F4p[i][j][k]);

        // F4v = a * dt;
        pMULTIPLY(a[i][j][k], jello->dt, F4v[i][j][k]);

        // buffer.p = F2p * 2.0;
        pMULTIPLY(F2p[i][j][k], 2, buffer.p[i][j][k]);

        // buffer.v = F3p * 2.0;
        pMULTIPLY(F3p[i][j][k], 2, buffer.v[i][j][k]);

        // buffer.p = buffer.p + buffer.v;
        pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);

        // buffer.p = buffer.p + F1p;
        pSUM(buffer.p[i][j][k], F1p[i][j][k], buffer.p[i][j][k]);

        // buffer.p = buffer.p + F4p;
        pSUM(buffer.p[i][j][k], F4p[i][j][k], buffer.p[i][j][k]);

        // buffer.p = buffer.p * 1/6;
        pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);

        // jello->p = buffer.p + jello->p
        pSUM(buffer.p[i][j][k], jello->p[i][j][k], jello->p[i][j][k]);

        // buffer.p = F2v * 2;
        pMULTIPLY(F2v[i][j][k], 2, buffer.p[i][j][k]);

        // buffer.v = F3v * 2;
        pMULTIPLY(F3v[i][j][k], 2, buffer.v[i][j][k]);

        // buffer.p = buffer.p + buffer.v;
        pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);

        // buffer.p = buffer.p + F1v;
        pSUM(buffer.p[i][j][k], F1v[i][j][k], buffer.p[i][j][k]);

        // buffer.p = buffer.p + F4v;
        pSUM(buffer.p[i][j][k], F4v[i][j][k], buffer.p[i][j][k]);

        // buffer.p = buffer.p * 1/6;
        pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);

        // jello->v = buffer.p + jello->v;
        pSUM(buffer.p[i][j][k], jello->v[i][j][k], jello->v[i][j][k]);
      }

  return;
}
