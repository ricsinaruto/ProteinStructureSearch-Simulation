#include "main.h"
#define LEN 8192

/* these are needed for visualization of text on window */
void printv(va_list args, const char* format) {
	char buf[LEN];
	char* ch = buf;
	vsnprintf(buf, LEN, format, args);
	while (*ch)	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *ch++);
}

void print(const char* format, ...) {
	va_list args;
	va_start(args, format);
	printv(args, format);
	va_end(args);
}

// this function is used for printing characters
void printAt(int x, int y, const char* format, ...) {
	va_list args;
	glWindowPos2i(x, y);
	va_start(args, format);
	printv(args, format);
	va_end(args);
}