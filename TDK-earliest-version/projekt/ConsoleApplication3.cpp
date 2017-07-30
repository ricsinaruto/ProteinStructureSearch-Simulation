// WindowsFormApplication.cpp : main project file.

#include "stdafx.h"
#include "Form1.h"

using namespace System;
using namespace WindowsFormApplication;

[STAThreadAttribute]

int main(array<System::String ^> ^args)
{
	Console::WriteLine(L"The windows is starting");
	Application::Run(gcnew Form1());

	return 0;
}