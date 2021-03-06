/*  Teo Price-Broncucia
 *  Final Project
 *  12/11/18
 *
 *  Based on work of Ex9 by Willem A. (Vlakkies) Schreuder. And previous homework 5
 *  by me.
 *
 *  Take a walk among the trees! Since the trees are randomly generated every
 *  walk will different. Now the trees and ground have image textures.
 *
 *  Key bindings:
 *  1, 2, 3    Toggle between orthogonal, perspective, and first person
 *  w, a, s, d Move position in first person
 *  arrows     Change view angle
 // u/j, i/k, o/l change relative red/blue/green
 // m pause light movement
 *  ESC        Exit

 */
#include "CSCIx229.h"
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "smoke.h"

int axes=0;       //  Display axes
int th=-30;         //  Horizontal view angle
int ph=30;         //  Elevation of view angle
double Px = 2.0;   // Position of viewer
double Py = 1.5;
double Pz = 2.0;
int fov=40;       //  Field of view (for perspective)
double asp=1;     //  Aspect ratio
double dim=2.5;   //  Size of world
bool source = true; //source





//  Macro for sin & cos in degrees
#define Cos(th) cos(3.1415926/180*(th))
#define Sin(th) sin(3.1415926/180*(th))


/*
 *  Convenience routine to output raster text
 *  Use VARARGS to make this more flexible
 */
#define LEN 8192  //  Maximum length of text string
void Print(const char* format , ...)
{
   char    buf[LEN];
   char*   ch=buf;
   va_list args;
   //  Turn the parameters into a character string
   va_start(args,format);
   vsnprintf(buf,LEN,format,args);
   va_end(args);
   //  Display the characters one at a time at the current raster position
   while (*ch)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,*ch++);
}

/*
 *  Set projection
 */
static void Project()
{
   //  Tell OpenGL we want to manipulate the projection matrix
   glMatrixMode(GL_PROJECTION);
   //  Undo previous transformations
   glLoadIdentity();
   //  Perspective transformation
   gluPerspective(fov,asp,dim/16,6*dim);

   // //  Switch to manipulating the model matrix
   glMatrixMode(GL_MODELVIEW);
   //  Undo previous transformations
   glLoadIdentity();
}



/*
 *  OpenGL (GLUT) calls this routine to display the scene
 */
void display()
{
   const double len=0.5;  //  Length of box
   //  Erase the window and the depth buffer
   glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
   // We want transarency so we disable the depth test
   glDisable(GL_DEPTH_TEST);

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE);
   //  Undo previous transformations
   glLoadIdentity();
   //  Perspective - set eye position
    double Ex = -2*dim*Sin(th)*Cos(ph);
    double Ey = +2*dim        *Sin(ph);
    double Ez = +2*dim*Cos(th)*Cos(ph);
    gluLookAt(Ex,Ey,Ez , 0,1.0,0 , 0,Cos(ph),0);

   draw_smoke();

   //  Draw box
   glDisable(GL_LIGHTING);
   glDisable(GL_TEXTURE_2D);
   glColor3f(1,1,1);
   // if (mode == 1 || mode == 2)
   // {
      glBegin(GL_LINES);
      glVertex3d(len,0.0,len);
      glVertex3d(len,0.0,-len);

      glVertex3d(len,0.0,-len);
      glVertex3d(-len,0.0,-len);

      glVertex3d(-len,0.0,-len);
      glVertex3d(-len,0.0,len);

      glVertex3d(-len,0.0,len);
      glVertex3d(len,0.0,len);

      glVertex3d(len,3*len,len);
      glVertex3d(len,3*len,-len);

      glVertex3d(len,3*len,-len);
      glVertex3d(-len,3*len,-len);

      glVertex3d(-len,3*len,-len);
      glVertex3d(-len,3*len,len);

      glVertex3d(-len,3*len,len);
      glVertex3d(len,3*len,len);


      glEnd();
   // }
   //  Display parameters
   glWindowPos2i(5,5);
   Print("Angle=%d,%d ",th,ph);

   //  Render the scene and make it visible
   glFlush();
   glutSwapBuffers();
}

/*
 *  GLUT calls this routine when the window is resized
 */
void idle()
{
   run_smoke(source);
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  GLUT calls this routine when an arrow key is pressed
 */
void special(int key,int x,int y)
{
   //  Right arrow key - increase angle by 5 degrees
   if (key == GLUT_KEY_RIGHT)
      th += 5;
   //  Left arrow key - decrease angle by 5 degrees
   else if (key == GLUT_KEY_LEFT)
      th -= 5;
   //  Up arrow key - increase elevation by 5 degrees
   else if (key == GLUT_KEY_UP)
      ph += 5;
   //  Down arrow key - decrease elevation by 5 degrees
   else if (key == GLUT_KEY_DOWN)
      ph -= 5;
   //  Keep angles to +/-360 degrees
   th %= 360;
   ph %= 360;
   //  Update projection
   Project();
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  GLUT calls this routine when a key is pressed
 */
void key(unsigned char ch,int x,int y)
{
   //  Exit on ESC
   if (ch == 27)
      exit(0);
   //  Reset view angle
   else if (ch == '0')
      th = ph = 0;
   //  Change field of view angle
   else if (ch == '-' && ch>1)
      fov--;
   else if (ch == '+' && ch<179)
      fov++;
  //turn source on or off
  else if (ch == 's')
     source = !source;
   //  Reproject
   Project();
   //  Animate if requested
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  GLUT calls this routine when the window is resized
 */
void reshape(int width,int height)
{
   //  Ratio of the width to the height of the window
   asp = (height>0) ? (double)width/height : 1;
   //  Set the viewport to the entire window
   glViewport(0,0, width,height);
   //  Set projection
   Project();
}

/*
 *  Start up GLUT and tell it what to do
 */
int main(int argc,char* argv[])
{


  setup();
   //  Initialize GLUT
   glutInit(&argc,argv);
   //  Request double buffered, true color window with Z buffering at 600x600
   glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
   glutInitWindowSize(700,700);
   glutCreateWindow("Teo Price-Broncucia");
   //  Set callbacks
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutSpecialFunc(special);
   glutKeyboardFunc(key);
   glutIdleFunc(idle);
   //  Pass control to GLUT so it can interact with the user
   glutMainLoop();
   return 0;
}
