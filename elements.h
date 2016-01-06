//
// Created by charlesfoxw on 11/29/15.
//

#ifndef HOMEWORK5_ELEMENTS_H
#define HOMEWORK5_ELEMENTS_H

#define MAX_SIZE            30

typedef enum { XY, YZ, XZ } Plane;

class Vertex {
public:
    float x, y, z;
    int intensity[3];
    Vertex();
    Vertex(float, float, float);
    void setIntensity(int, int, int);
};

class Vertex2D {
public:
    float x, y;
    int intensity[3];
    Vertex2D();
    Vertex2D(float, float);
    void setIntensity(int, int, int);
};

class DirectL {
public:
    Vertex start;
    Vertex end;
    DirectL();
    DirectL(Vertex, Vertex);
    float magnitude();
    float getMinX();
    float getMaxX();
    float getMinY();
    float getMaxY();
    float getMinZ();
    float getMaxZ();
};

class DirectL2D {
public:
    Vertex2D start;
    Vertex2D end;
    DirectL2D();
    DirectL2D(Vertex2D, Vertex2D);
    float getMinY();
    float getMaxY();
    float getBottomX();
    float getTopX();
};

class Vector {
public:
    float x;
    float y;
    float z;
    Vector();
    Vector(float, float, float);
    Vector(Vertex, Vertex);
    Vector normalize(int);
    float magnitude();
    float dotProduct(Vector other);
    Vector crossProduct(Vector other);
};

class Surface {
public:
    int edgeCount;
    DirectL* edgeList;
    Surface();
    Surface(Vertex*, int);
    void addEdge(DirectL);
    float getMaxAllX();
    float getMaxAllY();
    float getMaxAllZ();
    float getAvgX();
    float getAvgY();
    float getAvgZ();
};

class Surface2D {
public:
    int edgeCount;
    DirectL2D* edgeList;
    Surface2D();
    Surface2D(Vertex2D*, int);
    void addEdge(DirectL2D);
    void sortByYmin();
    float getMinAllY();
    float getMaxAllY();
};

class Polyhedron {
public:
    int vertexCount;
    int edgeCount;
    int surfaceCount;
    Vertex* vertexList;
    DirectL* edgeList;
    Surface* surfaceList;
    Vector* normalList;

    Polyhedron();
    Polyhedron(int, int, int, Vertex*, DirectL*, Surface*, Vector*);
    void addVertex(Vertex);
    void addEdge(DirectL);
    void addSurface(Surface);
    void sortByDistance(Plane);
};

class Circle {
public:
    Vertex2D center;
    float radius;
    Circle();
    Circle(Vertex2D, float);
};

class Sphere {
public:
    Vertex center;
    float radius;
    Sphere();
    Sphere(Vertex, float);
};

class SQMatrix {
public:
    float** matrix;
    int order;
    SQMatrix(int);
    float determinant();
    SQMatrix getMinor(int, int);
    SQMatrix product(SQMatrix);
    SQMatrix inverse();
};

#endif //HOMEWORK5_ELEMENTS_H
