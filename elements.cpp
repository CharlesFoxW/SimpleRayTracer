//
// Created by charlesfoxw on 11/29/15.
//

#include "elements.h"
#include <iostream>
#include <math.h>
using namespace std;

#define UNREACHABLE_VALUE    5000

Vertex::Vertex() {
    x = 0;
    y = 0;
    z = 0;
    for (int i = 0; i < 3; i++)
        intensity[i] = 0;
}

Vertex::Vertex(float the_x, float the_y, float the_z) {
    x = the_x;
    y = the_y;
    z = the_z;
    for (int i = 0; i < 3; i++)
        intensity[i] = 0;
}

void Vertex::setIntensity(int r, int g, int b) {
    intensity[0] = r;
    intensity[1] = g;
    intensity[2] = b;
}

Vertex2D::Vertex2D() {
    x = 0;
    y = 0;
    for (int i = 0; i < 3; i++)
        intensity[i] = 0;
}

Vertex2D::Vertex2D(float the_x, float the_y) {
    x = the_x;
    y = the_y;
    for (int i = 0; i < 3; i++)
        intensity[i] = 0;
}

void Vertex2D::setIntensity(int r, int g, int b) {
    intensity[0] = r;
    intensity[1] = g;
    intensity[2] = b;
}

DirectL::DirectL():start(), end() {
}

DirectL::DirectL(Vertex v1, Vertex v2) {
    start.x = v1.x;
    start.y = v1.y;
    start.z = v1.z;
    end.x = v2.x;
    end.y = v2.y;
    end.z = v2.z;
    for (int i = 0; i < 3; i++) {
        start.intensity[i] = v1.intensity[i];
        end.intensity[i] = v2.intensity[i];
    }
}

float DirectL::magnitude() {
    float delta_x = (end.x - start.x);
    float delta_y = (end.y - start.y);
    float delta_z = (end.z - start.z);
    return (float) sqrt(pow(delta_x, 2.0) + pow(delta_y, 2.0) + pow(delta_z, 2.0));
}

float DirectL::getMinX() {
    return (start.x < end.x ? start.x : end.x);
}
float DirectL::getMaxX() {
    return (start.x > end.x ? start.x : end.x);
}
float DirectL::getMinY() {
    return (start.y < end.y ? start.y : end.y);
}
float DirectL::getMaxY() {
    return (start.y > end.y ? start.y : end.y);
}
float DirectL::getMinZ() {
    return (start.z < end.z ? start.z : end.z);
}
float DirectL::getMaxZ() {
    return (start.z > end.z ? start.z : end.z);
}

DirectL2D::DirectL2D():start(), end() {
}

DirectL2D::DirectL2D(Vertex2D v1, Vertex2D v2) {
    start.x = v1.x;
    start.y = v1.y;
    end.x = v2.x;
    end.y = v2.y;
    for (int i = 0; i < 3; i++) {
        start.intensity[i] = v1.intensity[i];
        end.intensity[i] = v2.intensity[i];
    }
}

float DirectL2D::getMinY() {
    return (start.y < end.y ? start.y : end.y);
}
float DirectL2D::getMaxY() {
    return (start.y > end.y ? start.y : end.y);
}
float DirectL2D::getBottomX() {
    return (start.y < end.y ? start.x : end.x);
}
float DirectL2D::getTopX() {
    return (start.y > end.y ? start.x : end.x);
}

Vector::Vector() {
    x = 0;
    y = 0;
    z = 0;
}

Vector::Vector(float the_x, float the_y, float the_z) {
    x = the_x;
    y = the_y;
    z = the_z;
}

Vector::Vector(Vertex the_start, Vertex the_end) {
    x = the_end.x - the_start.x;
    y = the_end.y - the_start.y;
    z = the_end.z - the_start.z;
}

Vector Vector::normalize(int factor) {  //parameter factor -> the scaling factor of the coordinate if scaled.
    float new_x, new_y, new_z;
    new_x = x / magnitude();
    new_y = y / magnitude();
    new_z = z / magnitude();
    return Vector(new_x * factor, new_y * factor, new_z * factor);
}

float Vector::magnitude() {
    return (float)sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
}

float Vector::dotProduct(Vector other) {
    float t = x * other.x + y * other.y + z * other.z;
    return t;
}

Vector Vector::crossProduct(Vector other) {
    Vector product;
    product.x = (y * other.z) - (z * other.y);
    product.y = -((x * other.z) - (z * other.x));
    product.z = (x * other.y) - (y * other.x);
    return product;
}

Surface::Surface() {
    edgeList = new DirectL[MAX_SIZE];
    for (int i = 0; i < MAX_SIZE; i++) {
        edgeList[i].start = Vertex(UNREACHABLE_VALUE,UNREACHABLE_VALUE,UNREACHABLE_VALUE);
        edgeList[i].end = Vertex(UNREACHABLE_VALUE,UNREACHABLE_VALUE,UNREACHABLE_VALUE);
    }
    edgeCount = 0;
}

Surface::Surface(Vertex* listV, int numV) {
    edgeCount = 0;
    edgeList = new DirectL[numV];
    for (int i = 0; i < numV - 1; i++) {
        addEdge(DirectL(listV[i], listV[i+1]));
    }
    addEdge(DirectL(listV[numV-1], listV[0]));
}

void Surface::addEdge(DirectL newEdge) {
    edgeList[edgeCount] = newEdge;
    edgeCount++;
}

float Surface::getMaxAllX() {
    float max_x = -UNREACHABLE_VALUE;
    for (int k = 0; k < edgeCount; k++) {
        if (edgeList[k].getMaxX() > max_x)
            max_x = edgeList[k].getMaxX();
    }
    return max_x;
}
float Surface::getMaxAllY() {
    float max_y = -UNREACHABLE_VALUE;
    for (int k = 0; k < edgeCount; k++) {
        if (edgeList[k].getMaxY() > max_y)
            max_y = edgeList[k].getMaxY();
    }
    return max_y;
}
float Surface::getMaxAllZ() {
    float max_z = -UNREACHABLE_VALUE;
    for (int k = 0; k < edgeCount; k++) {
        if (edgeList[k].getMaxZ() > max_z)
            max_z = edgeList[k].getMaxZ();
    }
    return max_z;
}

float Surface::getAvgX() {
    float sum = 0.0;
    for (int i = 0; i < edgeCount; i++) {
        sum += edgeList[i].start.x;
    }
    return sum / (float) edgeCount;
}
float Surface::getAvgY() {
    float sum = 0.0;
    for (int i = 0; i < edgeCount; i++) {
        sum += edgeList[i].start.y;
    }
    return sum / (float) edgeCount;
}
float Surface::getAvgZ() {
    float sum = 0.0;
    for (int i = 0; i < edgeCount; i++) {
        sum += edgeList[i].start.z;
    }
    return sum / (float) edgeCount;
}


Surface2D::Surface2D() {
    edgeList = new DirectL2D[MAX_SIZE];
    for (int i = 0; i < MAX_SIZE; i++) {
        edgeList[i].start = Vertex2D(UNREACHABLE_VALUE,UNREACHABLE_VALUE);
        edgeList[i].end = Vertex2D(UNREACHABLE_VALUE,UNREACHABLE_VALUE);
    }
    edgeCount = 0;
}

Surface2D::Surface2D(Vertex2D* listV, int numV) {
    edgeCount = 0;
    edgeList = new DirectL2D[numV];
    for (int i = 0; i < numV - 1; i++) {
        addEdge(DirectL2D(listV[i], listV[i+1]));
    }
    addEdge(DirectL2D(listV[numV-1], listV[0]));
}

void Surface2D::addEdge(DirectL2D newEdge) {
    edgeList[edgeCount] = newEdge;
    edgeCount++;
}

void Surface2D::sortByYmin() {
    int i, j, min_index;
    for (j = 0; j < edgeCount - 1; j++) {
        min_index = j;
        for (i = j + 1; i < edgeCount; i++) {
            if (edgeList[i].getMinY() < edgeList[min_index].getMinY()) {
                min_index = i;
            }
        }
        if(min_index != j) {
            swap(edgeList[j], edgeList[min_index]);
        }
    }
}

float Surface2D::getMaxAllY() {
    float max_y = -UNREACHABLE_VALUE;
    for (int k = 0; k < edgeCount; k++) {
        if (edgeList[k].start.y > max_y)
            max_y = edgeList[k].start.y;
    }
    return max_y;
}
float Surface2D::getMinAllY() {
    float min_y = UNREACHABLE_VALUE;
    for (int k = 0; k < edgeCount; k++) {
        if (edgeList[k].start.y < min_y)
            min_y = edgeList[k].start.y;
    }
    return min_y;
}

Polyhedron::Polyhedron() {
    edgeList = new DirectL[MAX_SIZE];
    for (int i = 0; i < MAX_SIZE; i++) {
        edgeList[i].start = Vertex(UNREACHABLE_VALUE, UNREACHABLE_VALUE, UNREACHABLE_VALUE);
        edgeList[i].end = Vertex(UNREACHABLE_VALUE, UNREACHABLE_VALUE, UNREACHABLE_VALUE);
    }
    edgeCount = 0;
}

Polyhedron::Polyhedron(int numV, int numE, int numS, Vertex* listV, DirectL* listE,
                       Surface* listS, Vector* listN) {
    int i;
    vertexCount = 0;
    edgeCount = 0;
    surfaceCount = 0;
    vertexList = new Vertex[numV];
    edgeList = new DirectL[numE];
    surfaceList = new Surface[numS];
    normalList = new Vector[numS];

    for (i = 0; i < numV; i++) {
        addVertex(listV[i]);
    }
    for (i = 0; i < numE; i++) {
        addEdge(listE[i]);
    }
    for (i = 0; i < numS; i++) {
        normalList[surfaceCount] = listN[i];
        addSurface(listS[i]);
    }
}


void Polyhedron::addVertex(Vertex newVertex) {
    vertexList[vertexCount] = newVertex;
    vertexCount++;
}

void Polyhedron::addEdge(DirectL newEdge) {
    edgeList[edgeCount] = newEdge;
    edgeCount++;
}

void Polyhedron::addSurface(Surface newSurface) {
    surfaceList[surfaceCount] = newSurface;
    surfaceCount++;
}

/* For Painter's algorithm */
void Polyhedron::sortByDistance(Plane the_plane) {  //selection sort
    int i, j, max_index;
    for (j = 0; j < surfaceCount - 1; j++) {
        max_index = j;
        for (i = j + 1; i < surfaceCount; i++) {
            switch (the_plane) {
                case XY:
                    if ((UNREACHABLE_VALUE - surfaceList[i].getMaxAllZ())
                        > (UNREACHABLE_VALUE - surfaceList[max_index].getMaxAllZ())) {
                        max_index = i;
                    }
                    else if ((UNREACHABLE_VALUE - surfaceList[i].getMaxAllZ())
                             == (UNREACHABLE_VALUE - surfaceList[max_index].getMaxAllZ())) {
                        if ((UNREACHABLE_VALUE - surfaceList[i].getAvgZ())
                            > (UNREACHABLE_VALUE - surfaceList[max_index].getAvgZ())) {
                            max_index = i;
                        }
                    }
                    break;
                case YZ:
                    if ((UNREACHABLE_VALUE - surfaceList[i].getMaxAllX())
                        > (UNREACHABLE_VALUE - surfaceList[max_index].getMaxAllX())) {
                        max_index = i;
                    }
                    else if ((UNREACHABLE_VALUE - surfaceList[i].getMaxAllX())
                             == (UNREACHABLE_VALUE - surfaceList[max_index].getMaxAllX())) {
                        if ((UNREACHABLE_VALUE - surfaceList[i].getAvgX())
                            > (UNREACHABLE_VALUE - surfaceList[max_index].getAvgX())) {
                            max_index = i;
                        }
                    }
                    break;
                case XZ:
                    if ((UNREACHABLE_VALUE - surfaceList[i].getMaxAllY())
                        > (UNREACHABLE_VALUE - surfaceList[max_index].getMaxAllY())) {
                        max_index = i;
                    }
                    else if ((UNREACHABLE_VALUE - surfaceList[i].getMaxAllY())
                             == (UNREACHABLE_VALUE - surfaceList[max_index].getMaxAllY())) {
                        if ((UNREACHABLE_VALUE - surfaceList[i].getAvgY())
                            > (UNREACHABLE_VALUE - surfaceList[max_index].getAvgY())) {
                            max_index = i;
                        }
                    }
                    break;
            }
        }

        if(max_index != j) {
            swap(surfaceList[j], surfaceList[max_index]); // swap the changed max. value
        }
    }

}

Circle::Circle() {
    center = Vertex2D();
    radius = 1;
}

Circle::Circle(Vertex2D the_center, float the_radius) {
    center = the_center;
    radius = the_radius;
}

Sphere::Sphere() {
    center = Vertex();
    radius = 1;
}

Sphere::Sphere(Vertex the_center, float the_radius) {
    center = the_center;
    radius = the_radius;
}

SQMatrix::SQMatrix(int the_order) {
    order = the_order;
    matrix = new float*[order];
    for (int i = 0; i < order; i++) {
        matrix[i] = new float[order];
    }
}

/* recursive function calculating the determinant of the matrix */
float SQMatrix::determinant() {
    if (order < 1) {
        perror("cannot get determinant because the order is not a positive integer");
        exit(1);
    }
    // stop the recursion when matrix is a single element
    if( order == 1 )
        return matrix[0][0];

    // the determinant value
    float det = 0;
    for (int i = 0; i < order; i++) {
        // get minor at element (0,i)
        SQMatrix minorMatrix = getMinor(0, i);
        // the recusion
        det += (i % 2 == 1? -1.0 : 1.0) * matrix[0][i] * minorMatrix.determinant();
    }
    return det;
}

SQMatrix SQMatrix::getMinor(int row, int col) { // get the minor matrix at element[row][col]
    // indicate which col and row is being copied to dest
    SQMatrix minorMatrix = SQMatrix(order - 1); // reduce the order by 1.
    int colCount = 0, rowCount = 0;

    for(int i = 0; i < order; i++ ) {
        if( i != row ) {    // when i is not the row of the element
            colCount = 0;
            for(int j = 0; j < order; j++ ) {
                if( j != col ) {    // when j is not the col of the element
                    minorMatrix.matrix[rowCount][colCount] = matrix[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }
    return minorMatrix;
}

SQMatrix SQMatrix::product(SQMatrix other) {
    float sum = 0.0;
    SQMatrix prodMatrix = SQMatrix(order);
    if (order != other.order) {
        perror("order needs to be match for matrix mutiplication");
        exit(1);
    }
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            int index = 0;
            sum = 0.0;
            while (index < order) {
                sum += matrix[i][index] * other.matrix[index][j];
                index++;
            }
            prodMatrix.matrix[i][j] = sum;
        }
    }
    return prodMatrix;
}

SQMatrix SQMatrix::inverse() {
    SQMatrix invMatrix = SQMatrix(order);
    float inv_det = (float) (1.0 / determinant());
    // Use the determinant - minor algorithm:
    for(int j = 0; j < order; j++) {
        for(int i = 0; i < order; i++) {
            SQMatrix minorMatrix = getMinor(j,i);   // get the minor
            invMatrix.matrix[i][j] = inv_det * minorMatrix.determinant();
            if((i + j) % 2 == 1)
                invMatrix.matrix[i][j] = - invMatrix.matrix[i][j];
        }
    }
    return invMatrix;
}