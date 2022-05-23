#ifndef METHODS_H
#define METHODS_H

#include "header.h"
#include "pointUtils.h"
#include "utils.h"
#include "methods.h"

#define MINIMUM_POINTS 3     // minimum number of cluster
#define EPSILON 0.5  // distance for clustering, metre^2

int inputeMinMissColByGCS(vector<vector<int>> &miss_index, vector<Point> &points, vector<double> miss_count_vec, int cluster_count);
void imputeByGCS(vector<Point> &points, vector<vector<int>> &miss_index, int tar_col, vector<Point> all_means);
vector<int> CPIK(vector<Point> &points, int k, vector<vector<int>> &miss_index, vector<double> miss_count_vec, int start_idx, vector<Point> &centroids);
void imputeByEM(vector<Point> &points, vector<vector<int>> &miss_index, vector<double> all_means, int empty_start);

#endif
