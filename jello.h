/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _JELLO_H_
#define _JELLO_H_

#include "openGL-headers.h"
#include "pic.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>

#define pi 3.141592653589793238462643383279 

// user defined
struct spring;
struct AABB;
//

// camera angles
extern double Theta;
extern double Phi;
extern double R;
extern bool animate;
extern float color;

// number of images saved to disk so far
extern int sprite;

// mouse control
extern int g_vMousePos[2];
extern int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

struct point 
{
   double x;
   double y;
   double z;

  point(double xx = 0.0f, double yy = 0.0f, double zz = 0.0f)
    : x(xx), y(yy), z(zz) {}

   point operator+(point& p) {
    return point(x+p.x, y+p.y, z+p.z);
   }

   void print() {
    std::cout << "Point(" << x << ", " << y << ", " << z << ")\n";
   }
};


// these variables control what is displayed on the screen
extern int shear, bend, structural, pause, viewingMode, saveScreenToFile;

struct world
{
  char integrator[10]; // "RK4" or "Euler"
  double dt; // timestep, e.g.. 0.001
  int n; // display only every nth timepoint
  double kElastic; // Hook's elasticity coefficient for all springs except collision springs
  double dElastic; // Damping coefficient for all springs except collision springs
  double kCollision; // Hook's elasticity coefficient for collision springs
  double dCollision; // Damping coefficient collision springs
  double mass; // mass of each of the 512 control points, mass assumed to be equal for every control point
  int incPlanePresent; // Is the inclined plane present? 1 = YES, 0 = NO (always NO in this assignment)
  double a,b,c,d; // inclined plane has equation a * x + b * y + c * z + d = 0; if no inclined plane, these four fields are not used
  int resolution; // resolution for the 3d grid specifying the external force field; value of 0 means that there is no force field
  struct point * forceField; // pointer to the array of values of the force field
  struct point p[8][8][8]; // position of the 512 control points
  struct point v[8][8][8]; // velocities of the 512 control points
  std::vector<spring>* structureSprings;
  std::vector<spring>* shearSprings;
  std::vector<spring>* bendSprings;
  AABB* cube;
  double invM;
  point cellSize;
};

extern struct world jello;

struct pointIndex {
  int i, j, k;

  pointIndex(int ii = 0, int jj = 0, int kk = 0)
    : i(ii), j(jj), k(kk){}
};

struct spring
{
  pointIndex p1, p2;
  double res_len;
  spring(const pointIndex& pp1, const pointIndex& pp2, double r)
    : p1(pp1), p2(pp2), res_len(r) {}
};

struct collisionSpring: spring
{
  point collidePoint;
  collisionSpring(const pointIndex& pp1, const point& p): spring(pp1, {0, 0, 0}, 0.0f), collidePoint(p) {};
};

struct plane
{
  // ax + by + cz + d = 0;
  double m_a, m_b, m_c, m_d;

  plane(): m_a(0), m_b(0), m_c(0), m_d(0){};
  plane(double a, double b, double c, double d) : m_a(a), m_b(b), m_c(c), m_d(d) {}
  plane(const point& p1, const point& p2, const point& p3);
};

struct AABB
{
  // axis aligned bounding box with min and max point
  point m_min;
  point m_max;

  // the plane used to form AABB
  plane m_plane[6];

  AABB() : m_min({0, 0, 0}), m_max({0, 0, 0}){};
  AABB(const point& min, const point& max) {
    buildAABB(min, max);
  }

  void buildAABB(const point& min, const point& max);

};

// computes crossproduct of three vectors, which are given as points
// struct point vector1, vector2, dest
// result goes into dest
#define CROSSPRODUCTp(vector1,vector2,dest)\
  CROSSPRODUCT( (vector1).x, (vector1).y, (vector1).z,\
                (vector2).x, (vector2).y, (vector2).z,\
                (dest).x, (dest).y, (dest).z )

// computes crossproduct of three vectors, which are specified by floating-point coordinates
// double x1,y1,z1,x2,y2,z2,x,y,z
// result goes into x,y,z
#define CROSSPRODUCT(x1,y1,z1,x2,y2,z2,x,y,z)\
\
  x = (y1) * (z2) - (y2) * (z1);\
  y = (x2) * (z1) - (x1) * (z2);\
  z = (x1) * (y2) - (x2) * (y1)

// computes dotproduct of three vectors, which are given as points
// struct point vector1, vector2, dest
// result goes into dest
#define DOTPRODUCTp(vector1,vector2,dest)\
  DOTPRODUCT( (vector1).x, (vector1).y, (vector1).z,\
                (vector2).x, (vector2).y, (vector2).z,\
                dest)

// computes dotproduct of three vectors, which are specified by floating-point coordinates
// double x1,y1,z1,x2,y2,z2,x,y,z
// result goes into x,y,z
#define DOTPRODUCT(x1,y1,z1,x2,y2,z2,dest)\
  dest = x1 * x2 + y1 * y2 + z1 * z2

// normalizes vector dest
// struct point dest
// result returned in dest
// must declare a double variable called 'length' somewhere inside the scope of the NORMALIZE macrp
// macro will change that variable
#define pNORMALIZE(dest)\
\
  length = sqrt((dest).x * (dest).x + (dest).y * (dest).y + (dest).z * (dest).z);\
  (dest).x /= length;\
  (dest).y /= length;\
  (dest).z /= length;

// copies vector source to vector dest
// struct point source,dest
#define pCPY(source,dest)\
\
  (dest).x = (source).x;\
  (dest).y = (source).y;\
  (dest).z = (source).z;
  
// assigns values x,y,z to point vector dest
// struct point dest
// double x,y,z
#define pMAKE(xx,yy,zz,dest)\
\
  (dest).x = (xx);\
  (dest).y = (yy);\
  (dest).z = (zz);

// sums points src1 and src2 to dest
// struct point src1,src2,dest
#define pSUM(src1,src2,dest)\
\
  (dest).x = (src1).x + (src2).x;\
  (dest).y = (src1).y + (src2).y;\
  (dest).z = (src1).z + (src2).z;

// dest = src2 - src1
// struct point src1,src2,dest
#define pDIFFERENCE(src1,src2,dest)\
\
  (dest).x = (src1).x - (src2).x;\
  (dest).y = (src1).y - (src2).y;\
  (dest).z = (src1).z - (src2).z;

// mulitplies components of point src by scalar and returns the result in dest
// struct point src,dest
// double scalar
#define pMULTIPLY(src,scalar,dest)\
\
  (dest).x = (src).x * (scalar);\
  (dest).y = (src).y * (scalar);\
  (dest).z = (src).z * (scalar);

#define pDIVIDE(src1, src2, dest)\
\
  (dest).x = (src1).x / (src2).x;\
  (dest).y = (src1).y / (src2).y;\
  (dest).z = (src1).z / (src2).z;

#define pCREATE(i,j,k,dest)\
\
  (dest).x = (i);\
  (dest).y = (j);\
  (dest).z = (k);

#endif
