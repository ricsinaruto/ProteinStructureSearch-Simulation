#include "screencasts.h"


//Start up GLUT and tell it what to do
int main(int argc, char* argv[]) {
	windowName = "Dronpa Simulation Program";
	toggleAxes = 1;			//alapb�l a r�cs megjelen�t�s�vel induljon
	
	//glut ablack inicializ�l�sa
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);	//where to put the window

	main_window = glutCreateWindow(windowName);

	//glut inicializ�l� f�ggv�nyek
	glutDisplayFunc(display);
	glutReshapeFunc(displayReshape);
	glutKeyboardFunc(windowKey);
	glutSpecialFunc(windowSpecial);
	glutMouseFunc(windowMouse);
	glutPassiveMotionFunc(windowPmotion);
	glutMouseWheelFunc(mouseWheel);

	//saj�t inicializ�l� f�ggv�nyek
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