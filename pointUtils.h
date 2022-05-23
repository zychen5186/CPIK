#ifndef POINTUTILS_H
#define POINTUTILS_H

#include "header.h"

struct Point
{
public:
    vector<double> values;
    int dimension;
    int cluster;    // no default cluster
    double minDist; // default infinite dist to nearest cluster
    int index;

    Point();
    Point(vector<double> line, int num);
    Point(int size);
    ~Point();

    Point operator=(vector<double> p);
    Point operator-(Point p);

    int getDimension();
    void clearDist();
    void clearCluster();
    double distance(Point p);
    double distance(Point p, vector<double> miss_weight, vector<int> miss_index);
    double partialDistance(Point p, vector<int> exe_seq, int it_target);
};

void printPoints(vector<Point> p);
void normalizePoints(vector<Point> &points);
vector<Point> pointsMean(vector<Point> points);

/* Matrix operator */
vector<Point> Multiplication(vector<Point> M, vector<Point> N);
vector<Point> Transpose(vector<Point> M);
vector<Point> Inverse(vector<Point> M);
vector<Point> matrixMinus(vector<Point> M, vector<Point> N);
vector<Point> matrixAdd(vector<Point> M, vector<Point> N);
double minValue(vector<Point> points);
double maxValue(vector<Point> points);
vector<double> dimMin(vector<Point> points);
vector<double> dimMax(vector<Point> points);

#endif