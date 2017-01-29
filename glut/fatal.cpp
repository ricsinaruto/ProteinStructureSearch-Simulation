#include "screencasts.h"

//vmi error kezelõ cucc
void fatal(const char* format, ...)
{
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	exit(1);
}