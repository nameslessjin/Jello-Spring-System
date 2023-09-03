/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"



void printPoint(struct point &p) {
  std::cout << "P, " << " x: " << p.x << " y: " << p.y << " z: " << p.z << '\n';
}

/**
 * Compute hook elastic force on the spring between two points
 * @param p1, p2 - The two mass points connected by a spring
 * @param k      - The hook elastic coef
 * @param r      - resting length between p1, p2
 * @param f      - We set it the elastic force
*/
void computeElasticForce(const struct point& p1, const struct point& p2, double k, double r, struct point& f)
{

  // Hook's law: F = -k_h * (|L| - R) * (L/|L|)
  point l;
  double length;
  pDIFFERENCE(p1, p2, l);
  pNORMALIZE(l);
  pMULTIPLY(l, -k * (length - r), f);
}

/**
 * Compute damping force on the spring between two points
 * @param p1, p2 - The two mass points connected by a spring
 * @param v1, v2 - The velocitties of p1 and p2
 * @param k      - The damping coef
 * @param f      - We set it the elastic force
*/
void computeDampingForce(const struct point& p1, const struct point& p2, const struct point& v1, const struct point& v2, double k, struct point& f)
{
  // Damping force: -k_d * (v1 - v2) dot (L/|L|) * (L/|L|)
  point l, v;
  double length, veclocity;
  pDIFFERENCE(p1, p2, l);
  pNORMALIZE(l);
  pDIFFERENCE(v1, v2, v);
  DOTPRODUCTp(v, l, veclocity);
  pMULTIPLY(l, -k * veclocity, f);
}

void computeAccelerationSpring(const struct world *jello, struct point a[8][8][8], vector<spring>& springs)
{
  
  double invM = 1.0f / jello->mass;

  for (const spring& s: springs) {

    point elasticForce, dampingForce;

    computeElasticForce(jello->p[s.p1.i][s.p1.j][s.p1.k], jello->p[s.p2.i][s.p2.j][s.p2.k], jello->kElastic, s.res_len, elasticForce);
    computeDampingForce(jello->p[s.p1.i][s.p1.j][s.p1.k], jello->p[s.p2.i][s.p2.j][s.p2.k], 
                        jello->v[s.p1.i][s.p1.j][s.p1.k], jello->v[s.p2.i][s.p2.j][s.p2.k], jello->dElastic, dampingForce);

    point totalForce = elasticForce + dampingForce;

    // a = F/m
    pMULTIPLY(totalForce, invM, totalForce);

    // apply to p1
    pSUM(a[s.p1.i][s.p1.j][s.p1.k], totalForce, a[s.p1.i][s.p1.j][s.p1.k]);

    // apply to p2, on equal but negative direction force
    pMULTIPLY(totalForce, -1, totalForce);
    pSUM(a[s.p2.i][s.p2.j][s.p2.k], totalForce, a[s.p2.i][s.p2.j][s.p2.k]);
  }


}

void computeAccelerationCollisions(const struct world *jello, struct point a[8][8][8])
{
  std::vector<pointIndex> pointInds; // points that collides
  std::vector<point> collidedPoints; // points on the plane that collided with jello
  std::vector<collisionSpring> cSprings;

  double invM = 1.0f / jello->mass;

  // check collision and build collision spring
  if (checkCollision(jello, pointInds, collidedPoints))
  {
    for (int i = 0; i < pointInds.size(); ++i)
      cSprings.push_back(collisionSpring(pointInds[i], collidedPoints[i]));
  }

  // for every collision spring, compute Elastic and damping force
  for (const collisionSpring& s: cSprings)
  {
    point elasticForce, dampingForce;
    computeElasticForce(jello->p[s.p1.i][s.p1.j][s.p1.k], s.collidePoint, jello->kCollision, s.res_len, elasticForce);
    computeDampingForce(jello->p[s.p1.i][s.p1.j][s.p1.k], s.collidePoint, jello->v[s.p1.i][s.p1.j][s.p1.k], {0, 0, 0}, jello->dCollision, dampingForce);

    point totalForce = elasticForce + dampingForce;

    // a = F/m
    pMULTIPLY(totalForce, invM, totalForce);

    // apply to p1 only
    pSUM(a[s.p1.i][s.p1.j][s.p1.k], totalForce, a[s.p1.i][s.p1.j][s.p1.k]);

  }
}


inline bool inCube(const point& p, const AABB& aabb)
{
  return (aabb.m_min.x <= p.x && p.x <= aabb.m_max.x) 
  && (aabb.m_min.y <= p.y && p.y <= aabb.m_max.y)
  && (aabb.m_min.z <= p.z && p.z <= aabb.m_max.z);
}

inline bool checkPlaneCollision(const point& p, const plane& pl)
{
  // check the distance between plane and point
  double distance = pl.m_a * p.x + pl.m_b * p.y + pl.m_c * p.z + pl.m_d;

  // if distance > 0, positive, distance < 0 negative, distance == 0, on the plane
  return  distance < 0;
}

point findClosestPoint(const point& p, const plane& pl)
{
  // find distance / normal
  double distance = pl.m_a * p.x + pl.m_b * p.y + pl.m_c * p.z + pl.m_d;
  double denominator = std::sqrt(pl.m_a * pl.m_a + pl.m_b * pl.m_b + pl.m_c * pl.m_c);
  distance /= denominator;


  // cloest point Q on the plane to P: Q = P - d * N_unit
  // - N_unit is the direction, d is the length
  point cloest;
  cloest.x = p.x - distance * pl.m_a;
  cloest.y = p.y - distance * pl.m_b;
  cloest.z = p.z - distance * pl.m_c;

  return cloest;
}

bool checkCollision(const struct world *jello, std::vector<pointIndex>& pointInds, std::vector<point>& collidedPoints) 
{
  // for every point in the cube check if a point has collided with a boundary
  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j < 8; ++j)
    {
      for (int k = 0; k < 8; ++k)
      {
        point pt = jello->p[i][j][k];

        // if any point is not inside the cube, then there is a collision
        if (!inCube(pt, *(jello->cube)))
        {
          // find the plane that it collides with
          // store the point index, find its cloest point in the collided plane
          for (int planeInd = 0; planeInd < 6; ++planeInd)
          {
            if (checkPlaneCollision(pt, (jello->cube)->m_plane[planeInd]))
            {
              pointInds.push_back(pointIndex(i, j, k));
              collidedPoints.push_back(findClosestPoint(pt, (jello->cube)->m_plane[planeInd]));
            }
          }
        }
        
        // if there is a inclined plane, check collision with inclined plane
        if (jello->incPlanePresent == 1)
        {
          plane incPlane{jello->a, jello->b, jello->c, jello->d};
          if (checkPlaneCollision(pt, incPlane))
          {
            pointInds.push_back({i, j, k});
            collidedPoints.push_back(findClosestPoint(pt, incPlane));
          }
        }
      }
    }
  }

  return pointInds.size() != 0;
}


/* Computes acceleration to every control point of the jello cube,
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(const struct world *jello, struct point a[8][8][8])
{
  /* for you to implement ... */
  for (int i = 0; i < 8; ++i)
  {
    for (int j = 0; j < 8; ++j)
    {
      for (int k = 0; k < 8; ++k)
      {
        pMAKE(0, 0, 0, a[i][j][k]);
      }
    }
  }

  computeAccelerationSpring(jello, a, *jello->structureSprings);
  computeAccelerationSpring(jello, a, *jello->bendSprings);
  computeAccelerationSpring(jello, a, *jello->shearSprings);

  computeAccelerationCollisions(jello, a);

  
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
