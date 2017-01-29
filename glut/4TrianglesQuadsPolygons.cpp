#include "screencasts.h"

#define PI 3.1415926535898
#define Length(x) (sizeof (x)/sizeof *(x))

/*  Globals */
double dim = 2; /* dimension of orthogonal box */
char *windowName = "Hello OpenGL";
int windowWidth = 500;
int windowHeight = 500;

/*
*  Display the scene
*/
void display()
{
	//GL_TRIANGLES
	GLfloat verts[6][2] = { {1,1},{1,0},{0,1},
						{-1,-1},{-1,0},{0,-1} };
	//GL_TRIANGLE_STRIP
	GLfloat verts2[6][2] = { {-1.5,0},{-1,0.5},{-0.5,0},
							{0,0.5},{0.5,0},{1,0.5} };
	//GL_TRIANGLE_FAN
	GLfloat verts3[8][2] = { {-1.5,-1.5},{-1,0.8},{-0.5,1.2},
							{0,1},{0.5,0.6},{1,0} };
	//GL_QAUDS
	GLfloat verts4[8][2] = { {-1.5,-1.5},{-1.5,0},{0.0},{0,-1},
						{0.3,0},{1,0},{1,1},{0.4,0.8} };
	//GL_POLYGON
	GLfloat verts5[6][2] = { { 0.0,-1.0 },{ 0.5,-0.8 },{ 0.5,0.5 },
	{ 0.0, 1 },{ -0.5,0.8 },{ -0.5,-0.5 } };

	// White, Red, Green, Blue
	GLfloat colors[4][3] = { { 1.0,1.0,1.0 },{ 1.0,0.0,0.0 },
	{ 0.0,1.0,0.0 },{ 0.0,0.0,1.0 } };

	/* Clear the image */
	glClear(GL_COLOR_BUFFER_BIT);
	/* Reset previous transforms */
	glLoadIdentity();

	/* draw something */

	// GL_TRIANGLES
	// 4 vertices - last one ignored
	/*
	glBegin(GL_TRIANGLES);
	glColor3fv( colors[0] );
	glVertex2fv( verts[0] );
	glColor3fv( colors[1] );
	glVertex2fv( verts[1] );
	glColor3fv( colors[2] );
	glVertex2fv( verts[2] );
	glColor3fv( colors[3] );
	glVertex2fv( verts[3] );
	glEnd();
	*/
	// 6 vertices
	/*
	glBegin(GL_TRIANGLES);
	glColor3fv( colors[0] );
	glVertex2fv( verts[0] );
	glColor3fv( colors[1] );
	glVertex2fv( verts[1] );
	glColor3fv( colors[2] );
	glVertex2fv( verts[2] );
	glColor3fv( colors[0] );
	glVertex2fv( verts[3] );
	glColor3fv( colors[2] );
	glVertex2fv( verts[4] );
	glColor3fv( colors[3] );
	glVertex2fv( verts[5] );
	glEnd();
	*/

	// GL_TRIANGLE_STRIP
	/*
	glBegin(GL_TRIANGLE_STRIP);
	glColor3fv( colors[0] );
	glVertex2fv( verts2[0] );
	glColor3fv( colors[1] );
	glVertex2fv( verts2[1] );
	glColor3fv( colors[2] );
	glVertex2fv( verts2[2] );
	glColor3fv( colors[0] );
	glVertex2fv( verts2[3] );
	glColor3fv( colors[1] );
	glVertex2fv( verts2[4] );
	glColor3fv( colors[3] );
	glVertex2fv( verts2[5] );
	glEnd();
	*/
	// GL_TRIANGLE_FAN
	/*
	glBegin(GL_TRIANGLE_FAN);
	glColor3fv( colors[0] );
	glVertex2fv( verts3[0] );
	glColor3fv( colors[1] );
	glVertex2fv( verts3[1] );
	glColor3fv( colors[2] );
	glVertex2fv( verts3[2] );
	glColor3fv( colors[0] );
	glVertex2fv( verts3[3] );
	glColor3fv( colors[1] );
	glVertex2fv( verts3[4] );
	glColor3fv( colors[3] );
	glVertex2fv( verts3[5] );
	glEnd();
	*/
	// GL_QUADS
	// 6 vertices - last two ignored
	/*
	glBegin(GL_QUADS);
	glColor3fv( colors[0] );
	glVertex2fv( verts3[0] );
	glColor3fv( colors[1] );
	glVertex2fv( verts3[1] );
	glColor3fv( colors[2] );
	glVertex2fv( verts3[2] );
	glColor3fv( colors[3] );
	glVertex2fv( verts3[3] );
	glColor3fv( colors[1] );
	glVertex2fv( verts3[4] );
	glColor3fv( colors[3] );
	glVertex2fv( verts3[5] );
	glEnd();
	*/
	// 8 vertices
	/*
	glBegin(GL_QUADS);
	glColor3fv( colors[0] );
	glVertex2fv( verts4[0] );
	glColor3fv( colors[1] );
	glVertex2fv( verts4[1] );
	glColor3fv( colors[2] );
	glVertex2fv( verts4[2] );
	glColor3fv( colors[3] );
	glVertex2fv( verts4[3] );
	glColor3fv( colors[0] );
	glVertex2fv( verts4[4] );
	glColor3fv( colors[1] );
	glVertex2fv( verts4[5] );
	glColor3fv( colors[2] );
	glVertex2fv( verts4[6] );
	glColor3fv( colors[3] );
	glVertex2fv( verts4[7] );
	glEnd();
	*/
	// GL_QUAD_STRIP
	/*
	glBegin(GL_QUAD_STRIP);
	glColor3fv( colors[0] );
	glVertex2fv( verts2[0] );
	glColor3fv( colors[1] );
	glVertex2fv( verts2[1] );
	glColor3fv( colors[2] );
	glVertex2fv( verts2[2] );
	glColor3fv( colors[0] );
	glVertex2fv( verts2[3] );
	glColor3fv( colors[1] );
	glVertex2fv( verts2[4] );
	glColor3fv( colors[3] );
	glVertex2fv( verts2[5] );
	glEnd();
	*/
	// GL_POLYGON

	glBegin(GL_POLYGON);
	glColor3fv(colors[0]);
	glVertex2fv(verts5[0]);
	glColor3fv(colors[1]);
	glVertex2fv(verts5[1]);
	glColor3fv(colors[2]);
	glVertex2fv(verts5[2]);
	glColor3fv(colors[3]);
	glVertex2fv(verts5[3]);
	glColor3fv(colors[0]);
	glVertex2fv(verts5[4]);
	glColor3fv(colors[1]);
	glVertex2fv(verts5[5]);
	glEnd();

	




	/* Flush and swap */
	glFlush();
	glutSwapBuffers();

}

/*
*  GLUT calls this routine when the window is resized
*/
void reshape(int width, int height)
{
	double w2h = (height > 0) ? (double)width / height : 1;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-dim*w2h, +dim*w2h, -dim, +dim, -dim, +dim);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

}

void windowKey(unsigned char key, int x, int y)
{
	/* Exit on ESC */
	if (key == 27) exit(0);
}

/*
*  Start up GLUT and tell it what to do
*/
int main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(windowWidth, windowHeight);
	glutCreateWindow(windowName);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(windowKey);

	glutMainLoop();
	return 0;
}