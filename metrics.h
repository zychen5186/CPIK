#include "header.h"
#include "pointUtils.h"

double getAcc(vector<vector<double>> &maxCostMatrix, vector<int> &assignment, int size);
double getNMI(vector<int> &v1, vector<int> &v2);
double getNMI2(vector<int> &v1, vector<int> &v2);
double getFScore(vector<vector<double>> &maxCostMatrix, vector<int> &assignment);
double getFScore2(vector<int> &v1, vector<int> &v2);
double RMSE(vector<Point> &v1, vector<Point> &v2);