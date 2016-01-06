// Yihao Wang ECS 175 - Project 5 - Driver File
#include <iostream>
#include <fstream>
#include <GLUT/glut.h>
#include <string.h>
#include <math.h>
#include "AntTweakBar.h"
#include "elements.h"

using namespace std;

#define NUM_DIMENSION       3
#define PI                  3.14159265
#define GAP                 20
// Square Viewports and Windows
#define SUBWINDOW_SIZE      800
#define WINDOW_WIDTH        GAP * 3 + SUBWINDOW_SIZE + 400
#define WINDOW_HEIGHT       GAP * 2 + SUBWINDOW_SIZE
#define NUM_BITS_PIXEL      3
#define BUFF_SIZE           SUBWINDOW_SIZE * SUBWINDOW_SIZE * NUM_BITS_PIXEL * 3
#define MAX_NUM_OBJECTS     5
#define MAX_LEVEL           5


typedef enum { Red, Green, Blue, Yellow, Grey, White } Color;
typedef enum { R, G, B } ColorRGB;

/* ================================================== */
/* ================ Global Variables ================ */
/* ================================================== */
static GLint main_window, sub_window;
static TwBar *tweak_bar;
static float* pixels;
static int g_bound = SUBWINDOW_SIZE * 3 / 8;

static int g_numOfSpheres = 0;
static Sphere g_spheres[MAX_NUM_OBJECTS];
static int g_nearestIndex = 0;
static float g_viewingAngle = 0;

static float eye_x = 0;
static float eye_y = 0;
static float eye_z = 4;
static Vector baseE_x;
static Vector baseE_y;
static Vector baseE_z;
static float g_light_x = -12;
static float g_light_y = 0;
static float g_light_z = 0;
static int g_light_R = 100;
static int g_light_G = 100;
static int g_light_B = 100;
static int g_bright_R = 255;
static int g_bright_G = 255;
static int g_bright_B = 0;
static int g_ambient_R[MAX_NUM_OBJECTS];
static int g_ambient_G[MAX_NUM_OBJECTS];
static int g_ambient_B[MAX_NUM_OBJECTS];
static bool isReflecting = true;
static float g_Ka = 0.3;
static float g_Kd = 0.3;
static float g_Ks = 0.4;
static float g_Kr = 0.5;
static float g_Kt = 0.1;
static int intensity_R[2 * SUBWINDOW_SIZE * 3 / 8][2 * SUBWINDOW_SIZE * 3 / 8];
static int intensity_G[2 * SUBWINDOW_SIZE * 3 / 8][2 * SUBWINDOW_SIZE * 3 / 8];
static int intensity_B[2 * SUBWINDOW_SIZE * 3 / 8][2 * SUBWINDOW_SIZE * 3 / 8];

static Color g_color = Grey;
static int g_octant = -1;
static bool special_case = false;


/* ================================================== */
/* ================ Display Functions =============== */
/* ================================================== */
void reshape(int, int); //reshape the windows and subwindows
void keyboard(unsigned char, int, int); //add keyboard events to read from or write to files at runtime
void displayMain(void);
void displaySub(void);


/* ================================================== */
/* ================ Member Functions ================ */
/* ================================================== */
void controlDisplay(float*);
Vertex convertCVMToWorld(Vector*, Vertex, Vertex);

// Line Drawing:
void setPixelValue(float*, int, int);   // Handle the position of a point
void setLineValueDDA(float*, int, int, int, int);   // Use DDA here, From Project 1
void drawPixelAt(float*, int);  // Handle the color of a point
void getCurrWorld();
void confineSpaceByPt(float, float, float, float, float, float);
void stretchToNDC(Vertex*, int);    // Using NDC
void stretchSphereToNDC(Sphere*);

//Phong Lighting
int directPhongLt(Vertex, Vertex, Vector, int, int);
void renderPixel(float*, int, int, int*);

// Ray Tracing:
void calcColorSphere(float*, Vector*, Vertex);
Vertex intersection(Vertex, Vector);
int rayTracing(Vertex, Vector, int, ColorRGB);

// Coordinate Switch for Line Drawing - From Project 1
int getOctant(int, int, int, int);
void switchOctantInput(int, int*, int*);  //Changing coordinate.
void switchOctantOutput(int, int*, int*);   //Changing coordinate.

/* ================================================== */
/* =============== Main Program Starts ============== */
/* ================================================== */
int main(int argc, char** argv) {

    // Default initialization:
    g_numOfSpheres = 3;
    g_spheres[0] = Sphere(Vertex(1, 0, 0), 1);
    g_spheres[1] = Sphere(Vertex(-(float)(1.2), 0, 0), 1);
    g_spheres[2] = Sphere(Vertex((0), (0), (-50)), 45);
    g_ambient_R[0] = 30;
    g_ambient_G[0] = 200;
    g_ambient_B[0] = 200;
    g_ambient_R[1] = 255;
    g_ambient_G[1] = 255;
    g_ambient_B[1] = 255;
    g_ambient_R[2] = 200;
    g_ambient_G[2] = 55;
    g_ambient_B[2] = 55;
    g_viewingAngle = 30;

    baseE_x = Vector(1, 0, 0);
    baseE_y = Vector(0, 1, 0);
    baseE_z = Vector(0, 0, -1);

    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    glutInitWindowPosition(1500, 0);

    main_window = glutCreateWindow("Project 5 - Ray Tracing in CVM");
    glutReshapeFunc(reshape);
    glutDisplayFunc(displayMain);

    // Add the TweakBar
    TwInit(TW_OPENGL, NULL);
    // Set GLUT event callbacks
    // - Directly redirect GLUT mouse button events to AntTweakBar
    glutMouseFunc((GLUTmousebuttonfun)TwEventMouseButtonGLUT);
    // - Directly redirect GLUT mouse motion events to AntTweakBar
    glutMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
    // - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
    glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
    // - Directly redirect GLUT key events to AntTweakBar
    glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
    // - Directly redirect GLUT special key events to AntTweakBar
    glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
    // - Send 'glutGetModifers' function pointer to AntTweakBar;
    //   required because the GLUT key event functions do not report key modifiers states.
    TwGLUTModifiersFunc(glutGetModifiers);

    TwWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

    tweak_bar = TwNewBar("==Controller==");
    TwDefine(" GLOBAL help='This is the project to show the 3-D Ray Tracing for ECS 175.' ");
    // change default tweak bar size and color
    TwDefine(" ==Controller== position='840 20' size='400 800' color='100 110 110' ");

    TwAddVarRW(tweak_bar, "Set the Viewing Angle ", TW_TYPE_FLOAT, &g_viewingAngle,
               "step = 0.01 group='Eye point' help='Set viewing angle - x' ");

    TwAddVarRW(tweak_bar, "Set the Eye Point, X = ", TW_TYPE_FLOAT, &eye_x,
               "step = 0.01 group='Eye point' help='Set eye point - x' ");
    TwAddVarRW(tweak_bar, "Set the Eye Point, Y = ", TW_TYPE_FLOAT, &eye_y,
               "step = 0.01 group='Eye point' help='Set eye point - y' ");
    TwAddVarRW(tweak_bar, "Set the Eye Point, Z = ", TW_TYPE_FLOAT, &eye_z,
               "step = 0.01 group='Eye point' help='Set eye point - z' ");
    TwAddVarRW(tweak_bar, "Set the CVM base 1 vector: X= ", TW_TYPE_FLOAT, &baseE_x.x,
               "step = 0.01 group='Eye point' help='Set CVM - x' ");
    TwAddVarRW(tweak_bar, "Set the CVM base 1 vector: Y= ", TW_TYPE_FLOAT, &baseE_x.y,
               "step = 0.01 group='Eye point' help='Set CVM - x' ");
    TwAddVarRW(tweak_bar, "Set the CVM base 1 vector: Z= ", TW_TYPE_FLOAT, &baseE_x.z,
               "step = 0.01 group='Eye point' help='Set CVM - x' ");
    TwAddVarRW(tweak_bar, "Set the CVM base 2 vector: X= ", TW_TYPE_FLOAT, &baseE_y.x,
               "step = 0.01 group='Eye point' help='Set CVM - y' ");
    TwAddVarRW(tweak_bar, "Set the CVM base 2 vector: Y= ", TW_TYPE_FLOAT, &baseE_y.y,
               "step = 0.01 group='Eye point' help='Set CVM - y' ");
    TwAddVarRW(tweak_bar, "Set the CVM base 2 vector: Z= ", TW_TYPE_FLOAT, &baseE_y.z,
               "step = 0.01 group='Eye point' help='Set CVM - y' ");
    TwAddVarRW(tweak_bar, "Set the CVM base 3 vector: X= ", TW_TYPE_FLOAT, &baseE_z.x,
               "step = 0.01 group='Eye point' help='Set CVM - z' ");
    TwAddVarRW(tweak_bar, "Set the CVM base 3 vector: Y= ", TW_TYPE_FLOAT, &baseE_z.y,
               "step = 0.01 group='Eye point' help='Set CVM - z' ");
    TwAddVarRW(tweak_bar, "Set the CVM base 3 vector: Z= ", TW_TYPE_FLOAT, &baseE_z.z,
               "step = 0.01 group='Eye point' help='Set CVM - z' ");

    TwAddVarRW(tweak_bar, "Set the Light Source Point, X = ", TW_TYPE_FLOAT, &g_light_x,
               "step = 0.01 group='Light Source Point' help='Set light point - x' ");
    TwAddVarRW(tweak_bar, "Set the Light Source Point, Y = ", TW_TYPE_FLOAT, &g_light_y,
               "step = 0.01 group='Light Source Point' help='Set light point - y' ");
    TwAddVarRW(tweak_bar, "Set the Light Source Point, Z = ", TW_TYPE_FLOAT, &g_light_z,
               "step = 0.01 group='Light Source Point' help='Set light point - z' ");
    TwAddVarRW(tweak_bar, "Set the Light Source Color, Red = ", TW_TYPE_INT32, &g_light_R,
               "step = 1 group='Light Source Color' help='Set light point - x' ");
    TwAddVarRW(tweak_bar, "Set the Light Source Color, Green = ", TW_TYPE_INT32, &g_light_G,
               "step = 1 group='Light Source Color' help='Set light point - y' ");
    TwAddVarRW(tweak_bar, "Set the Light Source Color, Blue = ", TW_TYPE_INT32, &g_light_B,
               "step = 1 group='Light Source Color' help='Set light point - z' ");

    TwAddVarRW(tweak_bar, "Set the Phong Lightning Coefficent: K(ambient)", TW_TYPE_FLOAT, &g_Ka,
               " min=0 max=1 step = 0.01 group='Phong Lighting' help='Set K(amb)' ");
    TwAddVarRW(tweak_bar, "Set the Phong Lightning Coefficent: K(diffuse)", TW_TYPE_FLOAT, &g_Kd,
               " min=0 max=1 step = 0.01 group='Phong Lighting' help='Set K(diff)' ");
    TwAddVarRW(tweak_bar, "Set the Phong Lightning Coefficent: K(spectrum)", TW_TYPE_FLOAT, &g_Ks,
               " min=0 max=1 step = 0.01 group='Phong Lighting' help='Set K(spec)' ");

    TwAddVarRW(tweak_bar, "Set the Center of Color Sphere, X = ", TW_TYPE_FLOAT, &(g_spheres[0].center.x),
               "step=0.01 group='Set Color Sphere' help='Set sphere - x' ");
    TwAddVarRW(tweak_bar, "Set the Center of Color Sphere, Y = ", TW_TYPE_FLOAT, &(g_spheres[0].center.y),
               "step=0.01 group='Set Color Sphere' help='Set sphere - y' ");
    TwAddVarRW(tweak_bar, "Set the Center of Color Sphere, Z = ", TW_TYPE_FLOAT, &(g_spheres[0].center.z),
               "step=0.01 group='Set Color Sphere' help='Set sphere - z' ");
    TwAddVarRW(tweak_bar, "Set the Radius of Color Sphere:", TW_TYPE_FLOAT, &(g_spheres[0].radius),
               "step=0.01 group='Set Color Sphere' help='Set sphere - x' ");
    TwAddVarRW(tweak_bar, "Set the Color Sphere Color, Red = ", TW_TYPE_INT32, &g_ambient_R[0],
               "step = 1 group='Color Sphere Color' help='Set sphere color - x' ");
    TwAddVarRW(tweak_bar, "Set the Color Sphere Color, Green = ", TW_TYPE_INT32, &g_ambient_G[0],
               "step = 1 group='Color Sphere Color' help='Set sphere color - x' ");
    TwAddVarRW(tweak_bar, "Set the Color Sphere Color, Blue = ", TW_TYPE_INT32, &g_ambient_B[0],
               "step = 1 group='Color Sphere Color' help='Set sphere color - x' ");

    TwAddVarRW(tweak_bar, "Set the Center of Bright Sphere, X = ", TW_TYPE_FLOAT, &(g_spheres[1].center.x),
               "step=0.01 group='Set Bright Sphere' help='Set sphere - x' ");
    TwAddVarRW(tweak_bar, "Set the Center of Bright Sphere, Y = ", TW_TYPE_FLOAT, &(g_spheres[1].center.y),
               "step=0.01 group='Set Bright Sphere' help='Set sphere - y' ");
    TwAddVarRW(tweak_bar, "Set the Center of Bright Sphere, Z = ", TW_TYPE_FLOAT, &(g_spheres[1].center.z),
               "step=0.01 group='Set Bright Sphere' help='Set sphere - z' ");
    TwAddVarRW(tweak_bar, "Set the Radius of Bright Sphere:", TW_TYPE_FLOAT, &(g_spheres[1].radius),
               "step=0.01 group='Set Bright Sphere' help='Set sphere - x' ");

    TwAddVarRW(tweak_bar, "Set the Background Color, Red = ", TW_TYPE_INT32, &g_ambient_R[2],
               "step = 1 group='Background Color' help='Set Background color - x' ");
    TwAddVarRW(tweak_bar, "Set the Background Color, Green = ", TW_TYPE_INT32, &g_ambient_G[2],
               "step = 1 group='Background Color' help='Set Background color - x' ");
    TwAddVarRW(tweak_bar, "Set the Background Color, Blue = ", TW_TYPE_INT32, &g_ambient_B[2],
               "step = 1 group='Background Color' help='Set Background color - x' ");

    TwAddVarRW(tweak_bar, "Enable/Disable Reflection Mode", TW_TYPE_BOOL8, &isReflecting,
               " group='Reflection' help='Check this to enable Reflection' ");
    TwAddVarRW(tweak_bar, "Reflection Coefficient Kr = ", TW_TYPE_FLOAT, &g_Kr,
               " group='Reflection' help='Check this to change Kr' ");

    sub_window = glutCreateSubWindow(main_window, GAP, GAP, SUBWINDOW_SIZE, SUBWINDOW_SIZE);
    glutDisplayFunc(displaySub);
    glutKeyboardFunc(keyboard);

    glutMainLoop();

    return 0;
}

void displayMain(void) {
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    TwDraw();
    // Present frame buffer
    glutSwapBuffers();

    // Recall Display at next frame
    glutPostRedisplay();
}

void displaySub(void) {
    pixels = new float[BUFF_SIZE];
    memset(pixels, 0, sizeof(float) * BUFF_SIZE);

    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);

    controlDisplay(pixels);


    glDrawPixels(SUBWINDOW_SIZE, SUBWINDOW_SIZE, GL_RGB, GL_FLOAT, pixels);

    // Present frame buffer
    glutSwapBuffers();

    delete [] pixels;
    //glutPostRedisplay();
}

void reshape(int width, int height) {
    int sub_width = width - GAP * 3 - 400;
    int sub_height = height - GAP * 2;

    glutSetWindow(sub_window);
    glutPositionWindow(GAP, GAP);
    glutReshapeWindow(sub_width, sub_height);

    TwWindowSize(width, height);
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 'r':
            cout << "Refreshing..." << endl;
            for (int i = 0; i < 2 * g_bound; i++) {
                for (int j = 0; j < 2 * g_bound; j++) {
                    intensity_R[i][j] = 0;
                    intensity_G[i][j] = 0;
                    intensity_B[i][j] = 0;
                }
            }
            glutPostRedisplay();
            break;
        default:
            break;
    }
}

void controlDisplay(float* buffer) {

    // Draw the Viewport:
    g_color = Blue;
    setLineValueDDA(buffer, -g_bound, g_bound, g_bound, g_bound);
    setLineValueDDA(buffer, g_bound, g_bound, g_bound, -g_bound);
    setLineValueDDA(buffer, g_bound, -g_bound, -g_bound, -g_bound);
    setLineValueDDA(buffer, -g_bound, -g_bound, -g_bound, g_bound);

    // Draw Objects:
    Vector baseVSet[3];
    baseVSet[0] = baseE_x;
    baseVSet[1] = baseE_y;
    baseVSet[2] = baseE_z;

    calcColorSphere(buffer, baseVSet, Vertex(eye_x, eye_y, eye_z));

    // Draw the Axises:
    g_color = Grey;
    setLineValueDDA(buffer, 0, SUBWINDOW_SIZE / 2, 0, -SUBWINDOW_SIZE / 2);
    setLineValueDDA(buffer, SUBWINDOW_SIZE / 2, 0, -SUBWINDOW_SIZE / 2, 0);
}

Vertex convertCVMToWorld(Vector* baseV, Vertex eye, Vertex point) {
    Vector x_CVM, y_CVM, z_CVM;
    x_CVM = baseV[0].normalize(1);
    y_CVM = baseV[1].normalize(1);
    z_CVM = baseV[2].normalize(1);
    SQMatrix trans_matrix = SQMatrix(NUM_DIMENSION);

    trans_matrix.matrix[0][0] = x_CVM.x;
    trans_matrix.matrix[0][1] = y_CVM.x;
    trans_matrix.matrix[0][2] = z_CVM.x;
    trans_matrix.matrix[1][0] = x_CVM.y;
    trans_matrix.matrix[1][1] = y_CVM.y;
    trans_matrix.matrix[1][2] = z_CVM.y;
    trans_matrix.matrix[2][0] = x_CVM.z;
    trans_matrix.matrix[2][1] = y_CVM.z;
    trans_matrix.matrix[2][2] = z_CVM.z;

    Vertex new_point = Vertex((trans_matrix.matrix[0][0]*point.x
                               + trans_matrix.matrix[0][1]*point.y
                               + trans_matrix.matrix[0][2]*point.z),
                              (trans_matrix.matrix[1][0]*point.x
                               + trans_matrix.matrix[1][1]*point.y
                               + trans_matrix.matrix[1][2]*point.z ),
                              (trans_matrix.matrix[2][0]*point.x
                               + trans_matrix.matrix[2][1]*point.y
                               + trans_matrix.matrix[2][2]*point.z));
    Vertex translatedPoint = Vertex(new_point.x + eye.x, new_point.y + eye.y, new_point.z + eye.z);
    return translatedPoint;
}

void setLineValueDDA(float* buffer, int x1, int y1, int x2, int y2) {
    int delta_x, delta_y, index;
    int max_x, max_y, min_x, min_y;
    float val;

    g_octant = getOctant(x1, y1, x2, y2);
    if (g_octant == -1) {
        min_x = min(x1, x2);
        min_y = min(y1, y2);
        max_x = max(x1, x2);
        max_y = max(y1, y2);
        delta_x = max_x - min_x;
        delta_y = max_y - min_y;
        special_case = true;
        setPixelValue(buffer, min_x, min_y);
        if(delta_x == 0 && delta_y == 0)
            return;
        else if (delta_x == 0) {
            while (min_y != max_y) {
                min_y++;
                setPixelValue(buffer, x1, min_y);
            }
        }
        else {
            while (min_x != max_x) {
                min_x++;
                setPixelValue(buffer, min_x, y1);
            }
        }
    }
    else {
        index = 0;
        special_case = false;
        switchOctantInput(g_octant, &x1, &y1);
        switchOctantInput(g_octant, &x2, &y2);

        delta_x = x2 - x1;
        delta_y = y2 - y1;
        float slope_m = (float)delta_y / (float)delta_x;
        setPixelValue(buffer, x1, y1);
        while (x1 != x2) {
            x1++;
            val = (float)y1 + (float)index * slope_m;
            val += (val < 0 ? -0.5 : 0.5);  //Round val to the nearest integer.
            setPixelValue(buffer, x1, (int)val);
            index++;
        }
    }
}

void setPixelValue(float* buffer, int x, int y) {
    //set the pixel color to yellow:
    int pos, absolute_pos, absolute_x, absolute_y;
    if (!special_case)
        switchOctantOutput(g_octant, &x, &y);

    absolute_x = x + SUBWINDOW_SIZE / 2;
    absolute_y = y + SUBWINDOW_SIZE / 2;
    if (absolute_x > SUBWINDOW_SIZE * 7 / 8)
        return;
    if (absolute_x < SUBWINDOW_SIZE / 8)
        return;
    if (absolute_y > SUBWINDOW_SIZE * 7 / 8)
        return;
    if (absolute_y < SUBWINDOW_SIZE / 8)
        return;
    absolute_pos = SUBWINDOW_SIZE * absolute_y + absolute_x;
    pos = NUM_BITS_PIXEL * absolute_pos;
    drawPixelAt(buffer, pos);
}

void drawPixelAt(float* buff, int position) {
    switch (g_color) {
        case Red:
            buff[position] = 1.0;
            buff[position + 1] = 0.5;
            buff[position + 2] = 0.5;
            break;
        case Green:
            buff[position] = 0.2;
            buff[position + 1] = 0.9;
            buff[position + 2] = 0.5;
            break;
        case Blue:
            buff[position] = 0.1;
            buff[position + 1] = 0.7;
            buff[position + 2] = 1.0;
            break;
        case Yellow:
            buff[position] = 0.9;
            buff[position + 1] = 0.9;
            buff[position + 2] = 0.2;
            break;
        case Grey:
            buff[position] = 0.4;
            buff[position + 1] = 0.4;
            buff[position + 2] = 0.4;
            break;
        case White:
            buff[position] = 0.8;
            buff[position + 1] = 0.8;
            buff[position + 2] = 0.8;
            break;
    }
}

int directPhongLt(Vertex point, Vertex front, Vector normal, int i_a, int i_l) {
    float i_direct;
    float i_ambient, i_diffuse, i_spectrum;

    Vertex f = Vertex(front.x, front.y, front.z);
    Vertex x = Vertex(g_light_x, g_light_y, g_light_z);

    //Normalization:
    Vector f_p = Vector(point, f);
    Vector v = f_p.normalize(1);
    Vector l_raw = Vector(point, x);
    Vector l = l_raw.normalize(1);

    Vector n = normal.normalize(1);

    Vector r_raw = Vector(-l.x + 2 * n.dotProduct(l) * n.x,
                          -l.y + 2 * n.dotProduct(l) * n.y,
                          -l.z + 2 * n.dotProduct(l) * n.z);
    Vector r = r_raw.normalize(1);

    //Apply the Fumula:
    i_ambient = g_Ka * i_a;
    i_diffuse = i_l * g_Kd * l.dotProduct(n);
    i_spectrum = i_l * g_Ks * (float) pow(r.dotProduct(v), 1.0);
    i_direct = i_ambient + i_diffuse + i_spectrum;

    return (int) (i_direct + 0.5);
}

void calcColorSphere(float* buffer, Vector* baseSet, Vertex front) {

    // Apply Phong Model:
    for (int j = 0; j < 2 * g_bound; j++) {
        for (int i = 0; i < 2 * g_bound; i++) {

            float ex = (float) i / (float)(2 * g_bound) - (float)0.5;
            float ey = (float) j / (float)(2 * g_bound) - (float)0.5;
            float ez = 1 / (2 * (float)tan(g_viewingAngle * PI / 180));
            //float py = (float)sqrt(the_py) * SUBWINDOW_SIZE / (float)FACTOR * 3 / 4 / 2;
            Vertex p_world = convertCVMToWorld(baseSet, front, Vertex(ex, ey, ez));
            Vector the_d = Vector(front.x - p_world.x, front.y - p_world.y, front.z - p_world.z).normalize(1);
            //Vertex originalLight = Vertex(g_light_x, g_light_y, g_light_z);
            intensity_R[i][j] = rayTracing(front, the_d, 0, R);
            intensity_G[i][j] = rayTracing(front, the_d, 0, G);
            intensity_B[i][j] = rayTracing(front, the_d, 0, B);

        }
    }

    // Normalize:
    int maxIntensity = 0;
    for (int j = 0; j < 2 * g_bound; j++) {
        for (int i = 0; i < 2 * g_bound; i++) {
            if (intensity_R[i][j] > maxIntensity)
                maxIntensity = intensity_R[i][j];
            if (intensity_G[i][j] > maxIntensity)
                maxIntensity = intensity_G[i][j];
            if (intensity_B[i][j] > maxIntensity)
                maxIntensity = intensity_B[i][j];
        }
    }

    for (int j = 1; j < 2 * g_bound; j++) {
        for (int i = 1; i < 2 * g_bound; i++) {
            int intensity[3];
            intensity[0] = (int) (255 * ((float)intensity_R[i][j] / (float)maxIntensity) + 0.5);
            intensity[1] = (int) (255 * ((float)intensity_G[i][j] / (float)maxIntensity) + 0.5);
            intensity[2] = (int) (255 * ((float)intensity_B[i][j] / (float)maxIntensity) + 0.5);
            renderPixel(buffer, i - g_bound, j - g_bound, intensity);
        }
    }
}

Vertex intersection(Vertex origin, Vector the_d) {
    Vertex* the_point = new Vertex[g_numOfSpheres];
    int indices[MAX_NUM_OBJECTS];
    int index = 0;
    for (int i = 0; i < g_numOfSpheres; i++) {
        Sphere the_sphere = g_spheres[i];
        float coeff_a = (float) (pow(the_d.x, 2.0) + pow(the_d.y, 2.0) + pow(the_d.z, 2.0));
        float coeff_b = 2 * (origin.x * the_d.x + origin.y * the_d.y + origin.z * the_d.z)
                        - (2 * the_sphere.center.x * the_d.x + 2 * the_sphere.center.y * the_d.y +
                           2 * the_sphere.center.z * the_d.z);
        float coeff_c = (float) (pow(origin.x, 2.0) + pow(origin.y, 2.0) + pow(origin.z, 2.0) -
                                 pow(the_sphere.radius, 2.0))
                        + (float) pow(the_sphere.center.x, 2.0) - 2 * the_sphere.center.x * the_d.x
                        + (float) pow(the_sphere.center.y, 2.0) - 2 * the_sphere.center.y * the_d.y
                        + (float) pow(the_sphere.center.z, 2.0) - 2 * the_sphere.center.z * the_d.z;
        float the_t1 = 0, the_t2 = 0;
        float determiner = (float) pow(coeff_b, 2.0) - 4 * coeff_a * coeff_c;
        if (determiner < 0) {
            continue;
        }
        else {
            the_t1 = (-coeff_b - (float) sqrt(determiner)) / (2 * coeff_a);
            the_t2 = (-coeff_b + (float) sqrt(determiner)) / (2 * coeff_a);
        }
        if (the_t1 > 0 && the_t2 > 0) {
            continue;
        }
        else if (the_t1 <= 0 && the_t2 <= 0) {
            Vertex p1 = Vertex(origin.x + the_t1 * the_d.x, origin.y + the_t1 * the_d.y, origin.z + the_t1 * the_d.z);
            Vertex p2 = Vertex(origin.x + the_t2 * the_d.x, origin.y + the_t2 * the_d.y, origin.z + the_t2 * the_d.z);
            float distance1 = Vector(p1.x - origin.x, p1.y - origin.y, p1.z - origin.z).magnitude();
            float distance2 = Vector(p2.x - origin.x, p2.y - origin.y, p2.z - origin.z).magnitude();
            the_point[index] = distance1 < distance2 ? p1 : p2;
            indices[index] = i;
            index++;

        }
        else {
            float the_t = the_t2;
            Vertex p = Vertex(origin.x + the_t * the_d.x, origin.y + the_t * the_d.y, origin.z + the_t * the_d.z);
            float distance = Vector(p.x - origin.x, p.y - origin.y, p.z - origin.z).magnitude();
            the_point[index] = p;
            indices[index] = i;
            index++;
        }

    }

    Vertex nearestP;
    if (index != 0) {
        nearestP = the_point[0];
        float nearestDist = Vector(the_point[0].x - origin.x, the_point[0].y - origin.y,
                                   the_point[0].z - origin.z).magnitude();
        g_nearestIndex = indices[0];
        for (int i = 0; i < index; i++) {
            float distance = Vector(the_point[i].x - origin.x, the_point[i].y - origin.y,
                                    the_point[i].z - origin.z).magnitude();
            if (distance < nearestDist) {
                nearestDist = distance;
                nearestP = the_point[i];
                g_nearestIndex = indices[i];
            }
        }
    }
    else {
        nearestP = Vertex(5000, 5000, 5000);
    }
    delete [] the_point;
    return nearestP;
}

int rayTracing(Vertex the_origin, Vector direction, int level, ColorRGB colorType) {
    if (level > MAX_LEVEL)
        return 0;

    Vertex the_point = intersection(the_origin, direction);
    if (the_point.x == 5000 && the_point.y == 5000 && the_point.z == 5000) {
        return 0;
    }

    Vector raw_normal = Vector(the_point.x - g_spheres[g_nearestIndex].center.x,
                               the_point.y - g_spheres[g_nearestIndex].center.y,
                               the_point.z - g_spheres[g_nearestIndex].center.z);
    Vector normal = raw_normal.normalize(1);
    int intensity = 0;
    // Direct Phong:
    if (g_nearestIndex == 1) {
        switch (colorType) {
            case R:
                intensity = g_bright_R;
                break;
            case G:
                intensity = g_bright_G;
                break;
            case B:
                intensity = g_bright_B;
                break;
        }
    }
    else {
        switch (colorType) {
            case R:
                intensity = directPhongLt(the_point, the_origin, normal, g_ambient_R[g_nearestIndex], g_light_R);
                break;
            case G:
                intensity = directPhongLt(the_point, the_origin, normal, g_ambient_G[g_nearestIndex], g_light_G);
                break;
            case B:
                intensity = directPhongLt(the_point, the_origin, normal, g_ambient_B[g_nearestIndex], g_light_B);
                break;
        }
        // Global Phong:
        // reflection:
        Vector vision = Vector(direction.x, direction.y, direction.z).normalize(1);
        float the_product = normal.dotProduct(vision);
        Vector reflectRay = Vector(-vision.x + 2 * the_product * normal.x,
                                   -vision.y + 2 * the_product * normal.y,
                                   -vision.z + 2 * the_product * normal.z);
        reflectRay = Vector(-reflectRay.x, -reflectRay.y, -reflectRay.z).normalize(1);
        if (isReflecting)
            intensity += (int) (g_Kr * (float) rayTracing(the_point, reflectRay, level + 1, colorType) + 0.5);
    }

    // refraction:
    // TODO: refraction.

    if (intensity > 255) {
        intensity = 255;
    }
    return intensity;

}

void renderPixel(float* buffer, int the_x, int the_y, int* intensity) {
    int pos, absolute_pos, absolute_x, absolute_y;

    absolute_x = the_x + SUBWINDOW_SIZE / 2;
    absolute_y = the_y + SUBWINDOW_SIZE / 2;
    if (absolute_x > SUBWINDOW_SIZE * 7 / 8)
        return;
    if (absolute_x < SUBWINDOW_SIZE / 8)
        return;
    if (absolute_y > SUBWINDOW_SIZE * 7 / 8)
        return;
    if (absolute_y < SUBWINDOW_SIZE / 8)
        return;
    absolute_pos = SUBWINDOW_SIZE * absolute_y + absolute_x;
    pos = NUM_BITS_PIXEL * absolute_pos;
    for (int which_color = 0; which_color < 3; which_color++) {
        float curr_color = (float) intensity[which_color] / (float) 255;
        buffer[pos + which_color] = curr_color;
        //cout << "drew pixel at: " << pos  << " color is " << curr_color << "RGB = " << which_color << endl;
    }
}

int getOctant(int x1, int y1, int x2, int y2) {
    int delta_x = x2 - x1;
    int delta_y = y2 - y1;
    if (delta_x == 0 || delta_y == 0)
        return -1;

    float slope_m = (float)delta_y / (float)delta_x;
    if (delta_x > 0) {
        if(slope_m < -1)
            return 6;
        else if(slope_m >= -1 && slope_m <0)
            return 7;
        else if(slope_m >= 0 && slope_m < 1)
            return 0;
        else
            return 1;
    }
    else {
        if(slope_m < -1)
            return 2;
        else if(slope_m >= -1 && slope_m <0)
            return 3;
        else if(slope_m >= 0 && slope_m < 1)
            return 4;
        else
            return 5;
    }
}
void switchOctantInput(int oct, int* x, int* y) {
    int temp_x, temp_y;
    switch(oct) {
        case 0:
            break;
        case 1:
            //return (y, x)
            temp_x = *x;
            *x = *y;
            *y = temp_x;
            break;
        case 2:
            //return (y, -x)
            temp_x = *x;
            *x = *y;
            *y = -temp_x;
            break;
        case 3:
            //return (-x, y)
            temp_x = *x;
            *x = -temp_x;
            break;
        case 4:
            //return (-x, -y)
            temp_x = *x;
            temp_y = *y;
            *x = -temp_x;
            *y = -temp_y;
            break;
        case 5:
            //return (-y, -x)
            temp_x = *x;
            temp_y = *y;
            *x = -temp_y;
            *y = -temp_x;
            break;
        case 6:
            //return (-y, x)
            temp_x = *x;
            temp_y = *y;
            *x = -temp_y;
            *y = temp_x;
            break;
        case 7:
            //return (x, -y)
            temp_y = *y;
            *y = -temp_y;
            break;
        default:
            break;
    }
}

void switchOctantOutput(int oct, int* x, int* y) {
    int temp_x, temp_y;
    switch(oct) {
        case 0:
            break;
        case 1:
            //return (y, x)
            temp_x = *x;
            *x = *y;
            *y = temp_x;
            break;
        case 2:
            //return (-y, x)
            temp_x = *x;
            temp_y = *y;
            *x = -temp_y;
            *y = temp_x;
            break;
        case 3:
            //return (-x, y)
            temp_x = *x;
            *x = -temp_x;
            break;
        case 4:
            //return (-x, -y)
            temp_x = *x;
            temp_y = *y;
            *x = -temp_x;
            *y = -temp_y;
            break;
        case 5:
            //return (-y, -x)
            temp_x = *x;
            temp_y = *y;
            *x = -temp_y;
            *y = -temp_x;
            break;
        case 6:
            //return (y, -x)
            temp_x = *x;
            *x = *y;
            *y = -temp_x;
            break;
        case 7:
            //return (x, -y)
            temp_y = *y;
            *y = -temp_y;
            break;
        default:
            break;
    }
}