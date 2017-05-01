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
	int i;
	float angle;
	GLfloat vert1[2] = { 0,0.5 };
	GLfloat vert2[2] = { 0.5,-0.5 };
	GLfloat vert3[2] = { -0.5,-0.5 };
	GLfloat vert4[4][2] = { { 0.5,0.5 },{ 0.5,-0.5 },{ -0.5,-0.5 }, {-0.5,0.5 } };
	GLfloat color4[4][3] = { {1,1,1},{0,1,0},{0,0,1},{0,1,1} };
	GLfloat color1[3] = { 1,1,1 };
	GLfloat color23[2][3] = { {0,1,0},{0,0,1} };

	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	/* draw something */

	//GL_POINTS
	/*glPointSize(4);

	glBegin(GL_POINTS);
	glColor3fv(color1);
	glVertex2fv(vert1);
	glColor3fv(color23[0]);
	glVertex2fv(vert2);
	glColor3fv(color23[1]);
	glVertex2fv(vert3);
	glEnd();*/

	//GL_LINES
	//4 vertexes -> 2 line
	/*glBegin(GL_LINES);
	glColor3f(1, 1, 1);
	glVertex2f(0.5, 0.5);
	glColor3f(0, 1, 0);
	glVertex2f(0.5, -0.5);
	glColor3f(0, 0, 1);
	glVertex2f(-0.5, -0.5);
	glColor3f(0, 1, 1);
	glVertex2f(-0.5, 0.5);
	glEnd();*/

	//GL_LINE_STRIP
	/*glBegin(GL_LINE_STRIP);
	glColor3f(1, 1, 1);
	glVertex2f(0.5, 0.5);
	glColor3f(0, 1, 0);
	glVertex2f(0.5, -0.5);
	glColor3f(0, 0, 1);
	glVertex2f(-0.5, -0.5);
	glColor3f(0, 1, 1);
	glVertex2f(-0.5, 0.5);
	glEnd();*/

	//GL_LINE_LOOP
	/*glBegin(GL_LINE_LOOP);
	glColor3f(1, 1, 1);
	glVertex2f(0.5, 0.5);
	glColor3f(0, 1, 0);
	glVertex2f(0.5, -0.5);
	glColor3f(0, 0, 1);
	glVertex2f(-0.5, -0.5);
	glColor3f(0, 1, 1);
	glVertex2f(-0.5, 0.5);
	glEnd();*/

	//easy mode
	/*glBegin(GL_LINE_LOOP);
	for (i = 0; i < Length(color4);i++)
	{
		glColor3fv(color4[i]);
		glVertex2fv(vert4[i]);
	}
	glEnd();*/


	//circle+stipples
	glLineWidth(3);
	glEnable(GL_LINE_STIPPLE);
	glLineStipple(1, 0x00FF);
	glBegin(GL_LINE_LOOP);
	for (i = 0; i < 180; i++)
	{
	angle = 2 * PI*i / 180;
	glVertex2f(cos(angle), sin(angle));
	}
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