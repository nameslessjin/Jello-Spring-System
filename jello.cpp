/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

  Your name:
  <write your name here>

*/

#include <iostream>

#include "jello.h"
#include "showCube.h"
#include "input.h"
#include "physics.h"

using namespace std;

// camera parameters
double Theta = pi / 6;
double Phi = pi / 6;
double R = 6;

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

// number of images saved to disk so far
int sprite=0;

// these variables control what is displayed on screen
int shear=0, bend=0, structural=1, pause=0, viewingMode=0, saveScreenToFile=0;

struct world jello;

int windowWidth, windowHeight;

// custom global variable and forward declaration:

vector<spring> structureSprings, shearSprings, bendSprings;
AABB cube;
double length = 1.0f / 7;
float maxX = 2.0f, minX = -2.0f;
float maxY = 2.0f, minY = -2.0f;
float maxZ = 2.0f, minZ = -2.0f;
point minP = { minX, minY, minZ };
point maxP = { maxX, maxY, maxZ };
bool animate = false;
bool increasing = true;
float color = 255.0f;

void performAnimation();
void generateSprings();

// end of custom global variable nad forward declaration

void myinit()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(90.0,1.0,0.01,1000.0);

  // set background color to blue
  glClearColor(156.0f / color, 197.0f / color, 225.0f / color, 0.0);

  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  return; 
}

void reshape(int w, int h) 
{
  // Prevent a divide by zero, when h is zero.
  // You can't make a window of zero height.
  if(h == 0)
    h = 1;

  glViewport(0, 0, w, h);

  // Reset the coordinate system before modifying
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // Set the perspective
  double aspectRatio = 1.0 * w / h;
  gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); 

  windowWidth = w;
  windowHeight = h;

  glutPostRedisplay();
}

void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // camera parameters are Phi, Theta, R
  gluLookAt(R * cos(Phi) * cos (Theta), R * sin(Phi) * cos (Theta), R * sin (Theta),
	        0.0,0.0,0.0, 0.0,0.0,1.0);


  /* Lighting */
  /* You are encouraged to change lighting parameters or make improvements/modifications
     to the lighting model . 
     This way, you will personalize your assignment and your assignment will stick out. 
  */

  // global ambient light
  GLfloat aGa[] = { 0.0, 0.0, 0.0, 0.0 };
  
  // light 's ambient, diffuse, specular
  GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd0[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lKs0[] = { 1.0, 1.0, 1.0, 1.0 };

  GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd1[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs1[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd2[] = { 1.0, 1.0, 0.0, 1.0 };
  GLfloat lKs2[] = { 1.0, 1.0, 0.0, 1.0 };

  GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd3[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat lKs3[] = { 0.0, 1.0, 1.0, 1.0 };

  GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd4[] = { 0.0, 0.0, 1.0, 1.0 };
  GLfloat lKs4[] = { 0.0, 0.0, 1.0, 1.0 };

  GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd5[] = { 1.0, 0.0, 1.0, 1.0 };
  GLfloat lKs5[] = { 1.0, 0.0, 1.0, 1.0 };

  GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd6[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lKs6[] = { 1.0, 1.0, 1.0, 1.0 };

  GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd7[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat lKs7[] = { 0.0, 1.0, 1.0, 1.0 };

  // light positions and directions
  GLfloat lP0[] = { -1.999, -1.999, -1.999, 1.0 };
  GLfloat lP1[] = { 1.999, -1.999, -1.999, 1.0 };
  GLfloat lP2[] = { 1.999, 1.999, -1.999, 1.0 };
  GLfloat lP3[] = { -1.999, 1.999, -1.999, 1.0 };
  GLfloat lP4[] = { -1.999, -1.999, 1.999, 1.0 };
  GLfloat lP5[] = { 1.999, -1.999, 1.999, 1.0 };
  GLfloat lP6[] = { 1.999, 1.999, 1.999, 1.0 };
  GLfloat lP7[] = { -1.999, 1.999, 1.999, 1.0 };
  
  // jelly material color

  GLfloat mKa[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat mKd[] = { 0.3, 0.3, 0.3, 1.0 };
  GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mKe[] = { 0.0, 0.0, 0.0, 1.0 };

  /* set up lighting */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  // set up cube color
  glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
  glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
  glMaterialf(GL_FRONT, GL_SHININESS, 120);
    
  // macro to set up light i
  #define LIGHTSETUP(i)\
  glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
  glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
  glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
  glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
  glEnable(GL_LIGHT##i)
  
  LIGHTSETUP (0);
  LIGHTSETUP (1);
  LIGHTSETUP (2);
  LIGHTSETUP (3);
  LIGHTSETUP (4);
  LIGHTSETUP (5);
  LIGHTSETUP (6);
  LIGHTSETUP (7);

  // enable lighting
  glEnable(GL_LIGHTING);    
  glEnable(GL_DEPTH_TEST);

  // show the cube
  showCube(&jello);

  glDisable(GL_LIGHTING);

  // show the bounding box
  showBoundingBox(minP, maxP);
 
  glutSwapBuffers();
}

void doIdle()
{

  if (animate)
  {
    if (increasing)
    {
      if (Phi < 0.75) Phi += 0.003;
      else increasing = false;
    } else {
      if (Phi > -0.75) Phi -= 0.003;
      else increasing = true;
    }
  };

  char s[20]="picxxxx.ppm";
  int i;
  
  // save screen to file
  s[3] = 48 + (sprite / 1000);
  s[4] = 48 + (sprite % 1000) / 100;
  s[5] = 48 + (sprite % 100 ) / 10;
  s[6] = 48 + sprite % 10;

  if (saveScreenToFile==1)
  {
    saveScreenshot(windowWidth, windowHeight, s);
    saveScreenToFile=0; // save only once, change this if you want continuos image generation (i.e. animation)
    sprite++;
  }

  if (sprite >= 300) // allow only 300 snapshots
  {
    exit(0);	
  }

  if (pause == 0)
  {
    // insert code which appropriately performs one step of the cube simulation:
    performAnimation();
  }


  glutPostRedisplay();
}

int main (int argc, char ** argv)
{
  if (argc<2)
  {  
    printf ("Oops! You didn't say the jello world file!\n");
    printf ("Usage: %s [worldfile]\n", argv[0]);
    exit(0);
  }

  readWorld(argv[1],&jello);
  
  generateSprings();
  jello.structureSprings = &structureSprings;
  jello.shearSprings = &shearSprings;
  jello.bendSprings = &bendSprings;
  jello.invM = 1.0f / jello.mass;
  cube.buildAABB(minP, maxP);
  jello.cube = &cube;

  // compute cell size for external force field
  if (jello.resolution >= 2)
  {
    double invRes = 1.0f / (jello.resolution - 1);
    point p;
    pDIFFERENCE(jello.cube->m_max, jello.cube->m_min, p);
    pMULTIPLY(p, invRes, p);
    jello.cellSize = p;
  }

  glutInit(&argc,argv);
  
  /* double buffered window, use depth testing, 640x480 */
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  
  windowWidth = 640;
  windowHeight = 480;
  glutInitWindowSize (windowWidth, windowHeight);
  glutInitWindowPosition (0,0);
  glutCreateWindow ("Jello cube");

  /* tells glut to use a particular display function to redraw */
  glutDisplayFunc(display);

  /* replace with any animate code */
  glutIdleFunc(doIdle);

  /* callback for mouse drags */
  glutMotionFunc(mouseMotionDrag);

  /* callback for window size changes */
  glutReshapeFunc(reshape);

  /* callback for mouse movement */
  glutPassiveMotionFunc(mouseMotion);

  /* callback for mouse button changes */
  glutMouseFunc(mouseButton);

  /* register for keyboard events */
  glutKeyboardFunc(keyboardFunc);

  /* do initialization */
  myinit();

  /* forever sink in the black hole */
  glutMainLoop();

  return(0);
}

void performAnimation() {
  string integrator = jello.integrator;
  if (integrator == "RK4")
    RK4(&jello);
  else
    Euler(&jello);
}

bool isValidVertex(int x, int y, int z) {
  if (x < 0 || x > 7 || y < 0 || y > 7 || z < 0 || z > 7) return false;
  return true;
}

void processNeighborSpring(int i, int j, int k, int di, int dj, int dk, vector<spring>& springVec, double scale) {
  int ip = i + di, jp = j + dj, kp = k + dk;

  if (isValidVertex(ip, jp, kp))
    springVec.push_back(spring(pointIndex(i, j, k), pointIndex(ip, jp, kp), length * scale));
}

/**
 * create springs for structure, sheader and bending
*/
void generateSprings() {

  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j)
      for (int k = 0; k < 8; ++k)
        {
          // structure springs
          double scale = 1.0f;
          processNeighborSpring(i, j, k, 1, 0, 0, structureSprings, scale);
          processNeighborSpring(i, j, k, 0, 1, 0, structureSprings, scale);
          processNeighborSpring(i, j, k, 0, 0, 1, structureSprings, scale);

          // shear springs
          scale = std::sqrt(2.0f);
          processNeighborSpring(i, j, k, 1, 1, 0, shearSprings, scale);
          processNeighborSpring(i, j, k, -1, 1, 0, shearSprings, scale);
          processNeighborSpring(i, j, k, 1, 0, 1, shearSprings, scale);
          processNeighborSpring(i, j, k, -1, 0, 1, shearSprings, scale);
          processNeighborSpring(i, j, k, 0, 1, 1, shearSprings, scale);
          processNeighborSpring(i, j, k, 0, -1, 1, shearSprings, scale);

          scale = std::sqrt(3.0f);
          processNeighborSpring(i, j, k, 1, 1, 1, shearSprings, scale);
          processNeighborSpring(i, j, k, -1, 1, 1, shearSprings, scale);
          processNeighborSpring(i, j, k, -1, -1, 1, shearSprings, scale);
          processNeighborSpring(i, j, k, 1, -1, 1, shearSprings, scale);

          // bending springs
          scale = 2.0f;
          processNeighborSpring(i, j, k, 2, 0, 0, bendSprings, scale);
          processNeighborSpring(i, j, k, 0, 2, 0, bendSprings, scale);
          processNeighborSpring(i, j, k, 0, 0, 2, bendSprings, scale);
        }

  // std::cout<< structureSprings.size() << ' ' << shearSprings.size() << ' ' << bendSprings.size() << '\n';
}


plane::plane(const point& p1, const point& p2, const point& p3)
{
    // calculate plane coefficient
    // find the plane normal
    point v12, v23, n;
    pDIFFERENCE(p1, p2, v12);
    pDIFFERENCE(p2, p3, v23);
    CROSSPRODUCTp(v12, v23, n);
    pNORMALIZE(n);

    m_a = n.x, m_b = n.y, m_c = n.z;

    // d = -p dot n
    DOTPRODUCTp(p1, n, m_d);
    m_d *= -1;    
}

void AABB::buildAABB(const point& min, const point& max)
{
  m_min = min, m_max = max;

  // create corners
  point topLeftBack = {m_min.x, m_max.y, m_min.z};
  point topLeftFront = {m_min.x, m_max.y, m_max.z};
  point topRightBack = {m_max.x, m_max.y, m_min.z};
  point topRightFront = {m_max.x, m_max.y, m_max.z};

  point bottomLeftBack = {m_min.x, m_min.y, m_min.z};
  point bottomLeftFront = {m_min.x, m_min.y, m_max.z};
  point bottomRightBack = {m_max.x, m_min.y, m_min.z};
  point bottomRightFront = {m_max.x, m_min.y, m_max.z};

  // create planes, we want the normal point to the center of the box
  m_plane[0] = plane{ topLeftBack, topLeftFront, topRightFront }; // top
  m_plane[1] = plane{ bottomLeftBack, bottomRightBack, bottomRightFront }; // bottom
  m_plane[2] = plane{ topLeftBack, bottomLeftBack, bottomLeftFront }; // left;
  m_plane[3] = plane{ topRightBack, topRightFront, bottomRightFront }; // right;
  m_plane[4] = plane{ topRightFront, topLeftFront, bottomLeftFront }; // front;
  m_plane[5] = plane{ topLeftBack, topRightBack, bottomRightBack }; // back;

}
