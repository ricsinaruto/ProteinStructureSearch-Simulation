#include "screencasts.h"


//Start up GLUT and tell it what to do
int main(int argc, char* argv[]) {
	windowName = "Dronpa Simulation Program";
	toggleAxes = 1;			//start with showing the grid
	
							//initialize the glut window
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
	initializeProteins();

	//random seed
	srand(time(NULL));

	//glut loop functions
	redisplayAll();
	glutMainLoop();
	return EXIT_SUCCESS;
}