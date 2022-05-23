#include "utils.h"
#include "Hungarian.h"

const double FILLNUM = 1e+11;



/* input data handling */
tuple<int, int> getDataSize(string GT_file_name)
{
    int size = 0;
    int dimension = 0;
    string line;
    ifstream file;
    file.open(GT_file_name);
    getline(file, line);
    stringstream input(line);
    string i;
    while (input >> i)
    {
        dimension++;
    }
    file.close();
    file.open(GT_file_name);
    while (getline(file, line))
    {
        if (!(line == "NaN"))
            size++;
    }
    return {size, dimension};
}

vector<Point> readMissData(vector<vector<int>> &miss_index, string file_name, vector<int> nan_index)
/*=========== Read dataset with missing data into points ===========*/
{
    vector<Point> points;
    string line;
    ifstream file;
    file.open(file_name);
    int x = 0;
    int y = 0;
    int count = 0;

    while (getline(file, line))
    {
        /* if a data point is labeled as NaN, then ignore this point*/
        if (find(nan_index.begin(), nan_index.end(), count) != nan_index.end())
        {
            count++;
            continue;
        }
        else
        {
            string tmp = "";
            vector<double> values;
            for (int i = 0; i < (int)line.length(); i++) // one point's data
            {
                if ((48 <= int(line[i]) && int(line[i]) <= 57) || line[i] == '.' || line[i] == '+' || line[i] == '-' || line[i] == 'e')
                {
                    tmp += line[i];
                }
                else if (line[i] == 'N')
                    continue;
                else if (line[i] == 'a')
                {
                    tmp += to_string(FILLNUM);
                    miss_index[x][y] = 1;
                }
                else if (tmp.length() > 0) // execute when read space
                {
                    values.push_back(stod(tmp));
                    tmp = "";
                    y++;
                }
            }
            if (tmp.length() > 0)
            {
                values.push_back(stod(tmp));
                tmp = "";
            }
            x++;
            y = 0;
            points.push_back(Point(values, count));

            count++;
        }
    }
    file.close();

    return points;
}

vector<Point> readGTData(string file_name, vector<int> nan_index)
/*=========== Read dataset with missing data into points ===========*/
{
    vector<Point> points;
    string line;
    ifstream file;
    file.open(file_name);
    int count = 0;

    while (getline(file, line))
    {
        /* if a data point is labeled as NaN, then ignore this point*/
        if (find(nan_index.begin(), nan_index.end(), count) != nan_index.end())
        {
            count++;
            continue;
        }
        else
        {
            string tmp = "";
            vector<double> values;
            for (int i = 0; i < (int)line.length(); i++) // one point's data
            {
                if ((48 <= int(line[i]) && int(line[i]) <= 57) || line[i] == '.' || line[i] == '+' || line[i] == '-' || line[i] == 'e')
                {
                    tmp += line[i];
                }
                else if (tmp.length() > 0) // execute when read space
                {
                    values.push_back(stod(tmp));
                    tmp = "";
                }
            }
            if (tmp.length() > 0)
            {
                values.push_back(stod(tmp));
                tmp = "";
            }
            points.push_back(Point(values, count));

            count++;
        }
    }
    file.close();

    return points;
}

vector<int> getGT(string file_name, vector<string> &GT)
{
    fstream f;
    f.open(file_name, ios::in);
    string tok;
    vector<int> nan_index = {};
    int count = 0;
    while (getline(f, tok, '\n'))
    {
        if (tok == "")
            break;
        if (tok == "NaN")
            nan_index.push_back(count); // if label is NaN
        else
            GT.push_back(tok);
        count++;
    }
    return nan_index;
}

set<string> vecToSet(vector<string> v)
{
    set<string> s(v.begin(), v.end());
    return s;
}

int readGT(vector<int> &GT_vec, string dataset, vector<int> &nan_index)
{

    vector<string> GT;
    nan_index = getGT(dataset, GT); // nan_index is used to record the data points with NaN label

    set<string> GT_set;
    GT_set = vecToSet(GT);

    map<string, int> label_map;
    int cluster_count = 0;
    for (const auto &n : GT_set)
    {
        label_map[n] = cluster_count;
        cluster_count++;
    }

    for (vector<string>::iterator it = GT.begin(); it != GT.end(); it++)
    {
        GT_vec.push_back(label_map[*it]);
    }
    return cluster_count;
}

tuple<vector<int>, vector<int>> makePriorityAndMissingVec(vector<double> miss_count_vec)
{
    vector<pair<int, int>> M;
    int count = 0;
    for (auto const &n : miss_count_vec)
    {
        M.push_back(make_pair(n, count));
        count++;
    }

    sort(M.begin(), M.end());
    vector<int> missing_num = {};
    vector<int> priority_index_vector = {};
    for (auto const &n : M)
    {
        priority_index_vector.push_back(n.second);
        missing_num.push_back(n.first);
    }

    return {priority_index_vector, missing_num};
}

double getObj(vector<Point> points, vector<Point> centroids)
/*=========== calculate total distance between centroid and points ===========*/
{
    double sum = 0.0;
    for (vector<Point>::iterator it = points.begin(); it != points.end(); ++it)
    {
        int clusterId = it->cluster;
        sum += it->distance(centroids[clusterId]);
    }
    // sum = sqrt(sum);
    return sum;
}

double partialGetObj(vector<Point> points, vector<Point> centroids, vector<int> exeSeq, int itTarget)
/*=========== calculate total distance of the executed dimensions between centroid and points ===========*/
{
    double sum = 0.0;
    int point_count = points.size();
    // for (vector<Point>::iterator it = points.begin(); it != points.end(); ++it)
    for (int i = 0; i < point_count; i++)
    {
        int clusterId = points[i].cluster;
        sum += points[i].partialDistance(centroids[clusterId], exeSeq, itTarget);
    }
    sum = sqrt(sum);
    // cout << "Total distance: " << sum << endl;
    return sum;
}

void updateCentroids(vector<Point> points, vector<Point> &centroids, int k, vector<int> priority_col_vector, int last_col)
{
    /* Computing new centroids (with attributes before and including next priority) */
    int point_count = points.size();
    int dimension = points[0].getDimension();

    vector<vector<int>> nPoints(k, vector<int>(dimension, 0));
    vector<vector<double>> sum(k, vector<double>(dimension, 0));

    for (vector<Point>::iterator it = points.begin(); it != points.end(); ++it)
    {
        int clusterId = it->cluster;
        for (int i = 0; i <= last_col; i++)
        {
            if (it->values[priority_col_vector[i]] != FILLNUM)
            {
                nPoints[clusterId][priority_col_vector[i]] += 1;
                sum[clusterId][priority_col_vector[i]] += it->values[priority_col_vector[i]];
            }
        }
    }
    vector<double> dim_mean_vec(dimension);
    for (int i = 0; i <= last_col; i++)
    {
        dim_mean_vec[priority_col_vector[i]] = computeDimMean(points, priority_col_vector, priority_col_vector[i]);
    }
    for (vector<Point>::iterator c = begin(centroids); c != end(centroids); ++c)
    {
        int clusterId = c - begin(centroids);
        for (int i = 0; i <= last_col; i++)
        {
            if (nPoints[clusterId][priority_col_vector[i]] == 0)
                c->values[priority_col_vector[i]] = dim_mean_vec[priority_col_vector[i]];
            else
                c->values[priority_col_vector[i]] = sum[clusterId][priority_col_vector[i]] / nPoints[clusterId][priority_col_vector[i]];
        }
    }
    return;
}

vector<double> missCount(vector<vector<int>> &missIndex, int dim)
/* count missing numbers of each attricutes */
{
    int n = missIndex.size();
    vector<double> missCountVec(dim, 0);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (missIndex[i][j] == 1)
                missCountVec[j]++;
        }
    }

    return missCountVec;
}

double computeDimMean(vector<Point> points, vector<int> priority_col_vector, int TargetCol)
{
    double mean;
    double sum = 0;
    int count = 0;
    for (int j = 0; j < points.size(); j++)
    {
        if (points[j].values[priority_col_vector[TargetCol]] != FILLNUM)
        {
            sum += points[j].values[priority_col_vector[TargetCol]];
            count++;
        }
    }
    mean = sum / count;
    return mean;
}

void updatePoints_1(vector<vector<int>> &miss_index, vector<Point> &points, vector<Point> covariance, vector<Point> means)
// Impute missing values
{
    int point_count = points.size();
    int dimension = points[0].getDimension();
    for (int point_idx = 0; point_idx < point_count; point_idx++)
    {
        vector<Point> data_a(1);
        int miss_count = 0;
        int avail_count;

        for (int j = 0; j < dimension; j++)
        {
            if (miss_index[point_idx][j] == 0)
            {
                data_a[0].values.push_back(points[point_idx].values[j]);
            }
            else
            {
                miss_count++;
            }
        }
        avail_count = dimension - miss_count;

        if (miss_count == 0)
        {
            continue;
        }

        // Calculate partitioned means
        vector<Point> miu_a(1), miu_m(1);
        for (int j = 0; j < dimension; j++)
        {
            if (miss_index[point_idx][j] == 0)
            {
                miu_a[0].values.push_back(means[0].values[j]);
            }
            else
            {
                miu_m[0].values.push_back(means[0].values[j]);
            }
        }

        // Calculate partition covariance
        vector<Point> cov_mm(0);
        vector<Point> cov_am(0);
        vector<Point> cov_ma(0);
        vector<Point> cov_aa(0);

        for (int k = 0; k < dimension; k++)
        {
            if (miss_index[point_idx][k] == 1)
            {
                Point MM(0);
                Point MA(0);
                for (int h = 0; h < dimension; h++)
                {
                    if (miss_index[point_idx][h] == 1)
                    {
                        MM.values.push_back(covariance[k].values[h]);
                    }
                    else
                    {
                        MA.values.push_back(covariance[k].values[h]);
                    }
                }
                cov_mm.push_back(MM);
                cov_ma.push_back(MA);
            }
            else
            {
                Point AA(0);
                Point AM(0);
                for (int h = 0; h < dimension; h++)
                {
                    if (miss_index[point_idx][h] == 1)
                    {
                        AM.values.push_back(covariance[k].values[h]);
                    }

                    else
                    {
                        AA.values.push_back(covariance[k].values[h]);
                    }
                }
                cov_aa.push_back(AA);
                cov_am.push_back(AM);
            }
        }

        vector<Point> B = Multiplication(Inverse(cov_aa), cov_am);
        vector<Point> C = matrixMinus(cov_mm, Multiplication(Multiplication(cov_ma, Inverse(cov_aa)), cov_am));
        vector<Point> data_m = matrixAdd(miu_m, Multiplication(matrixMinus(data_a, miu_a), B));

        // Filling missing data
        int count = 0;
        for (int i = 0; i < dimension; i++)
        {
            if (miss_index[point_idx][i] == 1)
            {
                points[point_idx].values[i] = data_m[0].values[count];
                count++;
            }
        }

        /* ====== Recalculate covariance ====== */
        Point temp(dimension);
        vector<Point> compensate(dimension, temp);
        vector<Point> centered_points;

        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                if (i == j)
                    compensate[i].values[j] = 1e-8;
            }
        }   

        means = pointsMean(points);
        centered_points.assign(points.begin(), points.end());

        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < point_count; j++)
            {
                centered_points[j].values[i] -= means[0].values[i];
            }
        }

        covariance = Multiplication(Transpose(centered_points), centered_points);
        for (auto &i : covariance)
        {
            for (auto &j : i.values)
            {
                j /= (point_count - 1);
            }
        }
        covariance = matrixAdd(covariance, compensate);
    }
}

void updatePoints_2(vector<vector<int>> &miss_index, vector<Point> &points, vector<Point> covariance, vector<Point> means)
{
    int point_count = points.size();
    int dimension = points[0].getDimension();
    for (int point_idx = 0; point_idx < point_count; point_idx++)
    {
        vector<Point> data_a(1);
        int miss_count = 0;
        int avail_count;

        for (int j = 0; j < dimension; j++)
        {
            if (miss_index[point_idx][j] == 0)
            {
                data_a[0].values.push_back(points[point_idx].values[j]);
            }
            else
            {
                miss_count++;
            }
        }
        avail_count = dimension - miss_count;

        if (miss_count == 0)
        {
            continue;
        }

        // Calculate partitioned means
        vector<Point> miu_a(1), miu_m(1);
        for (int j = 0; j < dimension; j++)
        {
            if (miss_index[point_idx][j] == 0)
            {
                miu_a[0].values.push_back(means[0].values[j]);
            }
            else
            {
                miu_m[0].values.push_back(means[0].values[j]);
            }
        }

        // Calculate partition covariance
        vector<Point> cov_mm(0);
        vector<Point> cov_am(0);
        vector<Point> cov_ma(0);
        vector<Point> cov_aa(0);

        for (int k = 0; k < dimension; k++)
        {
            if (miss_index[point_idx][k] == 1)
            {
                Point MM(0);
                Point MA(0);
                for (int h = 0; h < dimension; h++)
                {
                    if (miss_index[point_idx][h] == 1)
                    {
                        MM.values.push_back(covariance[k].values[h]);
                    }
                    else
                    {
                        MA.values.push_back(covariance[k].values[h]);
                    }
                }
                cov_mm.push_back(MM);
                cov_ma.push_back(MA);
            }
            else
            {
                Point AA(0);
                Point AM(0);
                for (int h = 0; h < dimension; h++)
                {
                    if (miss_index[point_idx][h] == 1)
                    {
                        AM.values.push_back(covariance[k].values[h]);
                    }

                    else
                    {
                        AA.values.push_back(covariance[k].values[h]);
                    }
                }
                cov_aa.push_back(AA);
                cov_am.push_back(AM);
            }
        }

        vector<Point> B = Multiplication(Inverse(cov_aa), cov_am);
        vector<Point> C = matrixMinus(cov_mm, Multiplication(Multiplication(cov_ma, Inverse(cov_aa)), cov_am));
        // Filling missing data
        vector<Point> data_m = matrixAdd(miu_m, Multiplication(matrixMinus(data_a, miu_a), B));
        int count = 0;
        for (int i = 0; i < dimension; i++)
        {
            if (miss_index[point_idx][i] == 1)
            {
                points[point_idx].values[i] = data_m[0].values[count];
                count++;
            }
        }

        means = pointsMean(points);
        vector<Point> miu_2 = Multiplication(Transpose(means), means);
        // Calculate covariance
        Point temp(dimension);
        vector<Point> compensate(dimension, temp);
        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                if (i == j)
                    compensate[i].values[j] = 1e-8;
            }
        }   
        vector<Point> empty_covariance(dimension, temp);
        covariance = empty_covariance;
        for (int point_idx = 0; point_idx < point_count; point_idx++)
        {
            vector<Point> data_a(1), data_m(1);
            int miss_count = 0;
            int avail_count = 0;
            for (int j = 0; j < dimension; j++)
            {
                if (miss_index[point_idx][j] == 0)
                {
                    data_a[0].values.push_back(points[point_idx].values[j]);
                    avail_count++;
                }
                else
                {
                    data_m[0].values.push_back(points[point_idx].values[j]);
                    miss_count++;
                }
            }
            vector<Point> miu_a(1), miu_m(1);
            for (int j = 0; j < dimension; j++)
            {
                if (miss_index[point_idx][j] == 0)
                {
                    miu_a[0].values.push_back(means[0].values[j]);
                }
                else
                {
                    miu_m[0].values.push_back(means[0].values[j]);
                }
            }

            vector<Point> E_aa = Multiplication(Transpose(data_a), data_a);
            vector<Point> E_am = Multiplication(Transpose(data_a), data_m);
            vector<Point> E_ma = Multiplication(Transpose(data_m), data_a);
            vector<Point> C_hat = Multiplication(Transpose(matrixMinus(data_m, miu_m)), matrixMinus(data_m, miu_m));
            vector<Point> E_mm = matrixAdd(Multiplication(Transpose(data_m), data_m), C_hat);
            
            // assign the Es into the S matrix
            vector<Point> S_i(dimension, temp);
            int cnt_aa = 0;
            int cnt_am = 0;
            int cnt_ma = 0;
            int cnt_mm = 0;
            for (int k = 0; k < dimension; k++)
            {
                for (int h = 0; h < dimension; h++)
                {
                    if (miss_index[point_idx][k] == 0 && miss_index[point_idx][h] == 0)
                    {
                        S_i[k].values[h] = E_aa[cnt_aa / avail_count].values[cnt_aa % avail_count];
                        cnt_aa++;
                    }
                    else if (miss_index[point_idx][k] == 0 && miss_index[point_idx][h] == 1)
                    {
                        S_i[k].values[h] = E_am[cnt_am / miss_count].values[cnt_am % miss_count];
                        cnt_am++;
                    }
                    else if (miss_index[point_idx][k] == 1 && miss_index[point_idx][h] == 0)
                    {
                        S_i[k].values[h] = E_ma[cnt_ma / avail_count].values[cnt_ma % avail_count];
                        cnt_ma++;
                    }
                    else if (miss_index[point_idx][k] == 1 && miss_index[point_idx][h] == 1)
                    {
                        S_i[k].values[h] = E_mm[cnt_mm / miss_count].values[cnt_mm % miss_count];
                        cnt_mm++;
                    }
                }
            }

            S_i = matrixMinus(S_i, miu_2);
            covariance = matrixAdd(covariance, S_i);
        }
        for (auto &i : covariance)
        {
            for (auto &j : i.values)
            {
                j /= (point_count - 1);
            }
        }
        covariance = matrixAdd(covariance, compensate);
    }
}

double changeRatio(vector<Point> pre_points, vector<Point> points)
{
    int point_count = points.size();
    int dimension = points[0].getDimension();

    double error = 0;
    for (int i = 0; i < point_count; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            error += pow((points[i].values[j] - pre_points[i].values[j]), 2);
        }
    }
    error = sqrt(error);

    double sum = 0;
    for (int i = 0; i < point_count; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            sum += pow(points[i].values[j], 2);
        }
    }
    sum = sqrt(sum);
    if (isnan(error))
    {
        getchar();
    }
    return error / sum;
}

vector<vector<double>> getCostMatrix(vector<int> &v1, vector<int> &v2)
{
    int size = v1.size();
    // cout << v1.size() << ' ' << v2.size() << endl;
    int costMatrixSize = max(*max_element(v1.begin(), v1.end()), *max_element(v2.begin(), v2.end())) + 1; // 因為max會比size小一所以要加回來
    vector<vector<double>> maxCostMatrix(costMatrixSize, vector<double>(costMatrixSize, 0.0));

    for (int i = 0; i < size; i++)
    {
        maxCostMatrix[v1[i]][v2[i]] += 1;
    }

    vector<int> count(costMatrixSize, 0);

    for (int i = 0; i < costMatrixSize; i++)
    {
        for (int j = 0; j < costMatrixSize; j++)
        {
            // cout << maxCostMatrix[i][j] << " ";
            count[j] += maxCostMatrix[i][j];
        }
        // cout << endl;
    }

    // for(auto const &m: count)
    //     cout << m << " ";
    // cout << endl;

    return maxCostMatrix;
}

vector<int> getAssignment(vector<int> &v1, vector<int> &v2, vector<vector<double>> &maxCostMatrix)
{
    int size = v1.size();
    int costMatrixSize = max(*max_element(v1.begin(), v1.end()), *max_element(v2.begin(), v2.end())) + 1; // 因為max會比size小一所以要加回來
    vector<vector<double>> minCostMatrix(costMatrixSize, vector<double>(costMatrixSize, 0.0));

    vector<vector<double>>::iterator it;
    vector<double>::iterator it2;
    vector<int> tmp;
    for (it = maxCostMatrix.begin(); it != maxCostMatrix.end(); it++)
    {
        for (it2 = it->begin(); it2 != it->end(); it2++)
        {
            tmp.push_back(*it2);
            // cout << *it2 << " ";
        }
        // cout << endl;
    }

    int max_aggregate = *max_element(tmp.begin(), tmp.end());

    for (int i = 0; i < costMatrixSize; i++)
    {
        for (int j = 0; j < costMatrixSize; j++)
        {
            minCostMatrix[i][j] = max_aggregate - maxCostMatrix[i][j];
        }
    }

    HungarianAlgorithm HungAlgo;

    vector<int> assignment;
    HungAlgo.Solve(minCostMatrix, assignment);

    return assignment;
}

vector<int> getAdjustClusters(vector<int> &v1, vector<int> &assignment)
{
    vector<int> adjust_v1;
    for (int n : v1)
    {
        adjust_v1.push_back(assignment[n]);
    }

    return adjust_v1;
}