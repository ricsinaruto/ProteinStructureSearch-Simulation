
#include "screencasts.h"

/*
*  main()
*  ----
*  Start up GLUT and tell it what to do
*/
int main(int argc, char* argv[])
{
	
	
	initializeGlobals();
	/* screencast specific variables */
	windowName = "TDK";
	screencastID = 17;
	toggleAxes = 1;			//start with showing the grid
	

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);	//where to put the window

	main_window = glutCreateWindow(windowName);

	//glut initialization functions
	glutDisplayFunc(display);
	glutReshapeFunc(displayReshape);
	glutKeyboardFunc(windowKey);
	glutSpecialFunc(windowSpecial);
	glutMouseFunc(windowMouse);
	glutPassiveMotionFunc(windowPmotion);
	glutMouseWheelFunc(mouseWheel);

	//my own intialization functions
	initializeTextures();
	initializeObjs();
	lul();
	



	//creating the molecule grid, with 0 starting dipole values
	int i, j, k;
	for (i = 0; i < 38; i++)
	{
		for (j = 0; j < 38; j++)
		{
			for (k = 0; k < 38; k++)
			{
				dronpa[i][j][k].van = false;
				dronpa[i][j][k].ter = 0;
				dronpa[i][j][k].dipA = 0;
				dronpa[i][j][k].dipB = 0;
				dronpa[i][j][k].dip = 0;
				dronpa[i][j][k].qeA = 0;
				dronpa[i][j][k].qeB = 0;
				dronpa[i][j][k].qp1A = 0;
				dronpa[i][j][k].qp1B = 0;
				dronpa[i][j][k].qp2A = 0;
				dronpa[i][j][k].qp2B = 0;
			}
		}
	}
	
	//glut loop functions
	redisplayAll();
	glutMainLoop();
	return EXIT_SUCCESS;
}