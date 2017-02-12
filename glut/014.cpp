
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
	toggleAxes = 1;			//alapból a rács megjelenítésével induljon
	

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);	//where to put the window

	main_window = glutCreateWindow(windowName);

	//glut inicializáló függvények
	glutDisplayFunc(display);
	glutReshapeFunc(displayReshape);
	glutKeyboardFunc(windowKey);
	glutSpecialFunc(windowSpecial);
	glutMouseFunc(windowMouse);
	glutPassiveMotionFunc(windowPmotion);
	glutMouseWheelFunc(mouseWheel);

	//saját inicializáló függvények
	initializeTextures();
	initializeObjs();
	lul();
	
	srand(time(NULL));


	//molekula rács létrehozása, alapból 0 minden molekula dipólja, és false az értéke
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
	

	redisplayAll();
	glutMainLoop();
	return EXIT_SUCCESS;
}