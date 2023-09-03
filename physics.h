/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <iostream>
using namespace std;

void computeAcceleration(const struct world * jello, struct point a[8][8][8]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);


void computeAccelerationSpring(const struct world *jello, struct point a[8][8][8], vector<spring>& springs);
void computeDampingForce(const struct point& p1, const struct point& p2, const struct point& v1, const struct point& v2, double k, struct point& f);
void computeElasticForce(const struct point& p1, const struct point& p2, double k, double r, struct point& f);
void printPoint(struct point &p);

bool checkCollision(const struct world *jello, std::vector<pointIndex>& pointInds, std::vector<point>& collidedPoints) ;

#endif

