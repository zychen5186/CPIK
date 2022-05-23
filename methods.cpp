#include "methods.h"

const double FILLNUM = 1e+11;

int inputeMinMissColByGCS(vector<vector<int>> &miss_index, vector<Point> &points, vector<double> miss_count_vec, int cluster_count)
/* Impute missing data by mean in the min miss dimension */
{
    /* Make priority index vector */
    int point_count = points.size();
    int dimension = points[0].getDimension();

    vector<int> missing_num;
    vector<int> priority_col_vector;
    tie(priority_col_vector, missing_num) = makePriorityAndMissingVec(miss_count_vec);

    /* Impute missing data by GCS */
    bool flag = true;
    int start_idx = 0;
    vector<Point> means = pointsMean(points);
    while (flag)
    {
        imputeByGCS(points, miss_index, priority_col_vector[start_idx], means);
        cout << "column " << priority_col_vector[start_idx] << " is imputed by ICKNNI." << endl;

        // make sure after imputation, the number of different values in the first priority column is larger then cluster_count
        set<double> column_values;
        for (int i = 0; i < point_count; i++)
        {
            column_values.insert(points[i].values[priority_col_vector[start_idx]]);
        }
        column_values.erase(FILLNUM);

        if (column_values.size() >= cluster_count)
        {
            flag = false;
        }
        else
        {
            start_idx++;
        }
    }
    return start_idx; // number of dimensions that are used to start the PIK process.
}


void imputeByGCS(vector<Point> &points, vector<vector<int>> &miss_index, int tar_col, vector<Point> all_means)
{
    int point_count = points.size();
    if (point_count == 0)
    {
        cout << "This cluster has no points, can't process GCS" << endl;
        return;
    }

    int dim = points[0].getDimension();
    if (point_count == 1)
    {
        cout << "This cluster has only one points, impute by mean" << endl;
        for (int j = 0; j < dim; j++)
        {
            if (miss_index[0][j] == 1)
            {
                points[0].values[j] = all_means[0].values[j];
            }
        }
        return;
    }

    int k = 5;
    vector<Point> means = pointsMean(points);

    for (int row = 0; row < point_count; row++)
    {
        // Check which points also have available attributes as same as the target row
        vector<int> Candidate_pt_temp;
        for (int i = 0; i < point_count; i++)
        {
            if (i != row)
            {
                vector<int> temp = miss_index[i];
                for (int j = 0; j < dim; j++)
                {
                    temp[j] -= miss_index[row][j];
                }
                if (find(temp.begin(), temp.end(), 1) == temp.end())
                {
                    Candidate_pt_temp.push_back(i);
                }
            }
        }

        // Available entries in target row
        vector<int> available_entry;
        for (int i = 0; i < dim; i++)
        {
            if (miss_index[row][i] == 0)
                available_entry.push_back(i);
        }

        // Impute entry if miss_index[row][tar_col] == 1
        if (miss_index[row][tar_col] == 1)
        {
            vector<int> Candidate_pt;
            for (int i : Candidate_pt_temp)
            {
                if (miss_index[i][tar_col] == 0)
                    Candidate_pt.push_back(i);
            }

            // Calculate difference between target row and the candidate points
            vector<pair<double, int>> diff_vec;
            for (int i : Candidate_pt)
            {
                double sum = 0;
                double norm_a = 0;
                double norm_b = 0;
                // Only use the available attributes
                for (int j : available_entry)
                {
                    sum += points[i].values[j] * points[row].values[j];
                    norm_a += pow(points[i].values[j], 2);
                    norm_b += pow(points[row].values[j], 2);
                }
                norm_a = sqrt(norm_a);
                norm_b = sqrt(norm_b);
                double cosine = sum / (norm_a * norm_b);
                double angle = acos(cosine) * 180 / 3.141592;
                double length = max(norm_a, norm_b) / min(norm_a, norm_b) - 1;
                double dist = length + 10 * angle;

                diff_vec.push_back(make_pair(dist, i));
            }

            sort(diff_vec.begin(), diff_vec.end());

            // Impute
            // If candidate less than k, use all
            if (diff_vec.size() < k)
            {
                if (diff_vec.size() == 0)
                {
                    points[row].values[tar_col] = means[0].values[tar_col];
                }
                else
                {
                    int count = 0;
                    double sum = 0;
                    for (int i = 0; i < diff_vec.size(); i++)
                    {
                        sum += points[diff_vec[i].second].values[tar_col];
                        count++;
                    }
                    points[row].values[tar_col] = sum / count;
                }
            }
            else
            {
                double sum = 0;
                for (int i = 0; i < k; i++)
                {
                    sum += points[diff_vec[i].second].values[tar_col];
                }
                points[row].values[tar_col] = sum / k;
            }
        }
    }
}


void imputeByEM(vector<Point> &points, vector<vector<int>> &miss_index, vector<double> all_means, int empty_start)
{
    // Checking exceptation conditions
    int point_count = points.size();
    if (point_count == 0)
    {
        cout << "This cluster has no points, can't process EM" << endl;
        return;
    }

    int dimension = points[0].getDimension();
    if (point_count == 1)
    {
        cout << "This cluster has only one point, impute by mean" << endl;
        for (int j = 0; j < dimension; j++)
        {
            if (miss_index[0][j] == 1)
            {
                points[0].values[j] = all_means[j];
            }
        }
        return;
    }

    Point temp(dimension);
    vector<Point> means;
    vector<Point> centered_points;
    vector<Point> covariance;
    vector<Point> compensate(dimension, temp);

    means = pointsMean(points);

    // Center data to mean zero
    centered_points.assign(points.begin(), points.end());
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < point_count; j++)
        {
            centered_points[j].values[i] -= means[0].values[i];
        }
    }

    // Filling missing entries with zeros
    for (int i = 0; i < point_count; i++)
    {
        for (int j = empty_start; j < dimension; j++)
        {
            if (miss_index[i][j] == 1)
            {
                centered_points[i].values[j] = 0;
                points[i].values[j] = 0;
            }
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

    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            if (i == j)
                compensate[i].values[j] = 1e-8;
        }
    }
    covariance = matrixAdd(covariance, compensate);

    updatePoints_1(miss_index, points, covariance, means);

    // Parameter setting
    int iter = 1; // The first iter was finished above
    int maxit = 10;
    double obj = numeric_limits<double>::max();
    double threshold = 0.05;
    vector<Point> pre_points = points;

    while (iter < maxit && obj > threshold)
    {

        Point temp(dimension);
        vector<Point> means = pointsMean(points);
        vector<Point> miu_2 = Multiplication(Transpose(means), means);
        vector<Point> covariance(dimension, temp);


        for (int point_idx = 0; point_idx < point_count; point_idx++)
        {
            int miss_count = 0;
            int avail_count = 0;
            vector<Point> data_a(1), data_m(1);
            vector<Point> miu_a(1), miu_m(1);

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
            Point temp(dimension);
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

        updatePoints_2(miss_index, points, covariance, means);

        obj = changeRatio(pre_points, points);
        pre_points = points;
        iter += 1;
    }
    return;
}

/* Priority-based Incremental K-Means Clustering Algorithm */
vector<int> CPIK(vector<Point> &points, int k, vector<vector<int>> &miss_index, vector<double> miss_count_vec, int start_idx, vector<Point> &centroids)
{
    int point_count;
    int dimension;
    vector<int> centroidId;
    vector<int> missing_num;
    vector<int> priority_col_vector;

    point_count = points.size();
    dimension = points[0].getDimension();
    tie(priority_col_vector, missing_num) = makePriorityAndMissingVec(miss_count_vec);

    /*========== Randomly pick k points as centroid ==========*/
    centroids.clear();
    centroids.reserve(k);

    for (int i = 0; i < k; i++)
    {
        int rand_num;
        while (true)
        {
            rand_num = rand() % point_count;
            vector<int>::iterator it = std::find(centroidId.begin(), centroidId.end(), rand_num); // check if rand_num has already been picked.
            if (it == centroidId.end())
            {
                bool flag;
                for (int j = 0; j < i; j++)
                {
                    flag = false;
                    for (int h = 0; h <= start_idx; h++)
                    {
                        if (centroids[j].values[priority_col_vector[h]] != points.at(rand_num).values[priority_col_vector[h]])
                        {   // flag = true, if any attribute of the new picked centroid is different with the other centroids.
                            // flag = false, one of the existing centroid has identical attributes.
                            flag = true;
                        }
                    }
                    if (flag == false)
                    {
                        break;
                    }
                }
                if (i == 0 || flag == true)
                {
                    centroidId.push_back(rand_num);
                    centroids.push_back(points.at(rand_num));
                    break;
                }
            }
        }
    }
    /*========== Start CPIK process ==========*/
    int max_iter = 10;
    int this_priority;
    int next_priority;
    
    vector<int> clusters_vec;
    vector<int> processed_col = {};

    for (int i = 0; i <= start_idx; i++)
    { // Add the indexes before and including start_idx into the processed_col array.
        processed_col.push_back(priority_col_vector[i]);
    }

    this_priority = start_idx;         // priority_col_vector[this_priority] is the last inputed component.
    next_priority = this_priority + 1; // priority_col_vector[next_priority] is the last component needed to be imputed in this iteration.
    while (next_priority <= priority_col_vector.size())
    {
        // If multiple coponent has the same amount of missing attributes, impute together in this iteration.
        bool flag = true;
        while (flag)
        {
            if (next_priority + 1 < priority_col_vector.size())
            {
                if (missing_num[next_priority] != missing_num[next_priority + 1]) flag = false;
                else next_priority++;
            }
            else
            {
                flag = false;
            }
        }

        /*============ Clustering process ============*/
        int iter = 1;
        double pre_obj = 0;
        double cur_obj;
        while (true)
        {
            /*============ Assigning points to clusters ============*/
            for (vector<Point>::iterator c = begin(centroids); c != end(centroids); ++c)
            {
                int clusterId;
                clusterId = c - begin(centroids);
                for (int i = 0; i < point_count; i++)
                {
                    double dist = c->partialDistance(points[i], priority_col_vector, this_priority);
                    if (dist < points[i].minDist)
                    {
                        points[i].minDist = dist;
                        points[i].cluster = clusterId;
                    }
                }
            }
            updateCentroids(points, centroids, k, priority_col_vector, this_priority);

            /*============ Calculate the objective value ============*/
            bool done = false;
            cur_obj = partialGetObj(points, centroids, priority_col_vector, this_priority);

            if ((iter >= 2 && ((pre_obj - cur_obj) / cur_obj <= 0.0001)) || cur_obj == 0)
            {
                done = true;
            }

            if (done || max_iter <= iter)
            {
                cout << "Dimension " << priority_col_vector[this_priority] << ", iteration: " << iter << " Done!" << endl;
                break;
            }
            else
            {
                pre_obj = cur_obj;
                for (vector<Point>::iterator it = points.begin(); it != points.end(); ++it)
                {
                    it->clearCluster();
                    it->clearDist();
                }
                iter++;
            }
        }

        /*============ Stop if all components are imputed and partitioned ============*/
        if (next_priority == priority_col_vector.size())
        {
            double sse;

            for (vector<Point>::iterator it = points.begin(); it != points.end(); ++it)
            {
                clusters_vec.push_back(it->cluster);
            }

            for (vector<Point>::iterator it = points.begin(); it != points.end(); ++it)
            {
                sse += it->distance(centroids[it->cluster]);
            }

            cout << "SSE= " << sse << endl;

            break;
        }

        /*============== impute priority_col_vector[this_priority+1 ~ next_priority] missing data according to last clustering result ===============*/
        vector<int> impute_cols;
        vector<int>::iterator begin = priority_col_vector.begin();
        impute_cols.assign(begin, begin + next_priority + 1);
        vector<Point> means = pointsMean(points);

        // EM is performed on the components that are before priority_col_vector[next_priority] based on each cluster.
        for (int cluster_idx = 0; cluster_idx < k; cluster_idx++)
        {
            vector<double> em_means;
            vector<Point> em_points;
            vector<vector<int>> em_miss_index;

            for (int j : impute_cols)
            {
                em_means.push_back(means[0].values[j]);
            }

            for (int i = 0; i < point_count; i++)
            {
                if (points[i].cluster == cluster_idx) // only points within cluster_idx need to be considered
                {
                    vector<double> temp_line;
                    vector<int> temp_miss;
                    for (int j : impute_cols)
                    {
                        temp_line.push_back(points[i].values[j]);
                        temp_miss.push_back(miss_index[i][j]);
                    }
                    Point temp_p(temp_line, i);
                    em_points.push_back(temp_p);
                    em_miss_index.push_back(temp_miss);
                }
            }

            imputeByEM(em_points, em_miss_index, em_means, this_priority + 1);

            for (Point p : em_points)
            {
                int point_index = p.index;
                for (int col_index = 0; col_index < impute_cols.size(); col_index++)
                {
                    points[point_index].values[impute_cols[col_index]] = p.values[col_index];
                }
            }
        }

        updateCentroids(points, centroids, k, priority_col_vector, next_priority); // update centroids after next_priority components are imputed

        for (int i = this_priority + 1; i <= next_priority; i++)
        {
            cout << "Dimension " << priority_col_vector[i] << " ";
            processed_col.push_back(priority_col_vector[i]);
        }
        cout << "is imputed by EM" << endl;

        this_priority = next_priority;
        next_priority = this_priority + 1;
    }

    cout << "CPIK Finished!!!" << endl;
    return clusters_vec;
}