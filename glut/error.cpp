#include "screencasts.h"

//valami error cucc
void errCheck(char* where) {
	int err = glGetError();
	if (err) fprintf(stderr, "ERROR: %s [%s]\n", gluErrorString(err), where);
}

//vmi error kezelõ cucc
void fatal(const char* format, ...) {
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	exit(1);
}