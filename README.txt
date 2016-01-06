====================================================
Yihao Wang - 999015805 - ECS 175 Homework Project 5
====================================================

===============================
           Summary
===============================

Window Size: 1260 * 840
Subwindow Size: 800 * 800

Things are done:

	I made a GUI control panel based on AnttweakBar (An Open Source GUI Tool for OpenGL). It can be moved by dragging the title bar. (The files in lib directory and AntTweakBar.h are from this tool.)

	X, Y axises are shown by GREY lines.

	All the Curves are shown in colors to identify different types of curves.

	Done: 	Camera Viewing Model.

			Simple Ray Tracing for Spheres.

			Ray Tracing Spheres With Reflections.

			Every Change needs to press the 'r' key to take effect (After pressing please wait 5 sec.).
			

	The main driver function is actually controlDisplay(float* buffer) in line # 324, main.cpp

	The steps to render pixels are in calcColorSphere(float* buffer, Vector* baseSet, Vertex front) in line # 512, main.cpp.
	inside this function, the CVM will be used to calculate direction vector.
	Then call rayTracing() Recursively.
	int rayTracing(Vertex the_origin, Vector direction, int level, ColorRGB colorType) is in line # 627, main.cpp.
	

=============================================
                Instruction
=============================================

Build:

1. Use CMake to build my source code. (CMakeList.txt is included)
2. Use Makefile (Makefile included)
	$make
Then type:
	./hw5.out
Remove all output files:
	$make clean

--------------------------------------

No Input Files. The Objects are set using GUI control panel.

======================================
		        Manual
======================================

In the subwindow:
The GREY lines will indicate the coordinate axises: (origin is at the center of the viewport)

The Blue Frame is the ViewPort.

This is a GUI based UI Design. 
It is easy to follow by the selections on the control panel.
It supports direct inputs from keyboard for each parameter, +/- buttons by mouse, and mouse dragging reel at the right of +/- button on each entry.
Scrolling down to select different parameters.

Every change will not appear immediately on the viewports this time, because the program needs time to render.

Everything is operated through the Control Panel Displayed except for some functionality here:

	At anytime when running, if mouse pointer is stopped in the sub-windows area, 

	Press 'r' key -> Refresh the Viewport.


-----------------------
    Control Panel:
-----------------------
You can set:
	1. Anything related to CVM: viewing angle, eye point coordinates, CVM base vectors (these three have to be under some ralationship in order to work).

	2. Single Light Source Point with color and position.

	3. Phong lighting direct component coefficients.

	4. The Position and the Ambient Color of the Colored Sphere.

	5. The color of the Bright Sphere (default Yellow).

	6. The Background Color (default light Red).

	7. Enable or Disable the Reflection Recursion. (MAX = 5)

	8. The Reflection Parameter Kr.

