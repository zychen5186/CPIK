#ifndef utils_H
#define utils_H

#include "header.h"
#include "pointUtils.h"
#include "methods.h"



/* input data handling */
tuple<int, int> getDataSize(string GT_file_name);
vector<Point> readMissData(vector<vector<int>> &miss_index, string file_name, vector<int> nan_index);
vector<Point> readGTData(string file_name, vector<int> nan_index);
vector<Point> readCompleteData(string file_name);
vector<int> getGT(string file_name, vector<string> &GT);
set<string> vecToSet(vector<string> v);
int readGT(vector<int> &GT_vec, string dataset, vector<int> &nan_index);




/* CPIK functions */
tuple<vector<int>, vector<int>> makePriorityAndMissingVec(vector<double> miss_count_vec);
vector<double> missCount(vector<vector<int>> &missIndex, int dim);
double computeDimMean(vector<Point> points, vector<int> priority_col_vector, int TargetCol);
double getObj(vector<Point> points, vector<Point> centroids);
double partialGetObj(vector<Point> points, vector<Point> centroids, vector<int> exeSeq, int itTarget);
void updateCentroids(vector<Point> points, vector<Point>& centroids, int k, vector<int> priority_col_vector, int last_col);

/* EM functions */
void updatePoints_1(vector<vector<int>> &miss_index, vector<Point>& points, vector<Point> covariance, vector<Point> means);
void updatePoints_2(vector<vector<int>> &miss_index, vector<Point>& points, vector<Point> covariance, vector<Point> means);
double changeRatio(vector<Point> pre_points, vector<Point> points);

/* metric functions */
vector<vector<double>> getCostMatrix(vector<int> &v1, vector<int> &v2);
vector<int> getAssignment(vector<int> &v1, vector<int> &v2, vector<vector<double>> &maxCostMatrix);
vector<int> getAdjustClusters(vector<int> &v1, vector<int> &assignment);

#endif