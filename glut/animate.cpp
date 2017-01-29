#include "screencasts.h"

/*
* moveCube()
* ------
* Rotates a cube if animations are on
*/
int moveCube(int i)
{
	
	cubes[i].tsr.r.y = fmod(cubes[i].tsr.r.y + 1.0, 360.0);
	
	return 0;
}

/*
* moveLight()
* ------
* Moves the light if it is on
*/
void moveLight(void)
{
	if (toggleLight) lightPh = (lightPh + 10) % 360;
}

/*
* timer()
* ------
* Timer called periodically
*/
void timer(int value)
{	//ide kell írni a cuccokat animáláshoz
	
	
	if (toggleAnimation != DEF_ANIMATE_PAUSED) 
	{
		moveLight();
	}
	
	/*
	* glutTimerFunc(millisecs, timerCallback, value)
	*   millisecs = how long until you want the callback to be called
	*   timerCallback = the function to be called when time is up
	*   value = can be used for whatever. I'm using it toggle the animation
	*/
	glutTimerFunc(20, timer, toggleAnimation);

	redisplayAll();
}