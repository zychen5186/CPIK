#include "pointUtils.h"

const double FillNum = 1e+11;

Point::Point() : cluster(-1), minDist(numeric_limits<double>::max()) {}

Point::Point(vector<double> line, int num)
{
    values.assign(line.begin(), line.end());
    cluster = UNCLASSIFIED;
    minDist = numeric_limits<double>::max();
    dimension = values.size();
    index = num;
}

Point::Point(int size)
{
    values.resize(size);
    fill(values.begin(), values.end(), 0);
    cluster = -1;
    minDist = numeric_limits<double>::max();
    dimension = values.size();
}

Point::~Point() {}

Point Point::operator=(vector<double> p)
{
    values.clear();
    for (int i = 0; i < p.size(); i++)
    {
        values.push_back(p[i]);
    }
}

Point Point::operator-(Point p)
{
    Point temp;
    for (int i = 0; i < p.getDimension(); i++)
    {
        temp.values.push_back(this->values[i] - p.values[i]);
    }
    return temp;
}

int Point::getDimension() { return values.size(); }
void Point::clearDist() { this->minDist = numeric_limits<double>::max(); }
void Point::clearCluster() { this->cluster = -1; }

double Point::partialDistance(Point p, vector<int> exe_seq, int this_priority)
/*=========== calculate dist between specific attributes of two points ===========*/
{
    double dist = 0.0;
    for (int it = 0; it <= this_priority; it++)
    {
        dist += pow(values[exe_seq[it]] - p.values[exe_seq[it]], 2.0);
    }
    return dist;
}

double Point::distance(Point p)
/*=========== calculate the distance between two points ===========*/
{
    double dist = 0.0;
    for (int i = 0; i < this->values.size(); i++)
    {
        dist += pow(this->values[i] - p.values[i], 2.0);
    }
    return dist;
}

double Point::distance(Point p, vector<double> miss_weight, vector<int> miss_index)
/*=========== calculate the distance between two points ===========*/
{
    double dist = 0.0;
    for (int i = 0; i < this->values.size(); i++)
    {
        if (miss_index[i] == 1)
            dist += pow(this->values[i] - p.values[i], 2.0) * miss_weight[i];
        else
        {
            dist += pow(this->values[i] - p.values[i], 2.0);
        }
    }
    return dist;
}

void printPoints(vector<Point> p)
{
    int count = 0;
    for (auto &i : p)
    {
        cout << setw(5) << i.cluster << " ";
        cout << setw(5) << count << " ";
        for (auto &j : i.values)
        {
            cout << setw(5) << j << " ";
        }
        cout << endl;
        count++;
    }
    cout << endl;
}

void normalizePoints(vector<Point> &points)
{
    int point_count = points.size();
    if (point_count == 0)
    {
        cout << "there are no points !" << endl;
        return;
    }
    int dim = points[0].getDimension();

    // Calculate mean of each dimension
    vector<double> mean(dim, 0.0);
    vector<int> count(dim, 0);
    for (int i = 0; i < point_count; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (points[i].values[j] != FillNum)
            {
                count[j]++;
                mean[j] += points[i].values[j];
            }
        }
    }

    for (int i = 0; i < dim; i++)
    {
        mean[i] /= count[i];
    }

    // Calculate max and min of each dimension
    vector<vector<double>> minMax(dim, vector<double>(2, 0.0));
    vector<double> dimRange(dim, 0.0);
    for (int i = 0; i < dim; i++)
    {
        if (points[0].values[i] == FillNum)
        {
            minMax[i][0] = numeric_limits<double>::max();
            minMax[i][1] = -numeric_limits<double>::max();
        }
        else
        {
            minMax[i][0] = points[0].values[i];
            minMax[i][1] = points[0].values[i];
        }
    }

    for (int i = 0; i < point_count; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (points[i].values[j] != FillNum && points[i].values[j] > minMax[j][1])
                minMax[j][1] = points[i].values[j];
            else if (points[i].values[j] != FillNum && points[i].values[j] < minMax[j][0])
                minMax[j][0] = points[i].values[j];
        }
    }

    for (int i = 0; i < dim; i++)
    {
        dimRange[i] = minMax[i][1] - minMax[i][0];
    }

    // Normalize
    for (int i = 0; i < point_count; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (points[i].values[j] != FillNum)
            {
                if (dimRange[j] == 0)
                    points[i].values[j] = (points[i].values[j] - mean[j]);
                else
                    points[i].values[j] = (points[i].values[j] - minMax[j][0]) / dimRange[j]; // *
            }
        }
    }
}

vector<Point> pointsMean(vector<Point> points)
{
    int point_count = points.size();
    int dimension = points[0].getDimension();
    vector<Point> means(1);

    for (int i = 0; i < dimension; i++)
    {
        double sum = 0;
        int count = 0;
        for (int j = 0; j < point_count; j++)
        {
            if (points[j].values[i] != FillNum)
            {
                sum += points[j].values[i];
                count++;
            }
        }
        if (count == 0)
        {
            means[0].values.push_back(0);
        }
        else
        {
            means[0].values.push_back(sum / count);
        }
    }
    return means;
}

/* Matrix operator */
vector<Point> Multiplication(vector<Point> M, vector<Point> N)
{
    if (M.size() == 0 || N.size() == 0)
    {
        vector<Point> temp(0);
        return temp;
    }

    else if (M[0].getDimension() != N.size())
    {
        cerr << "These two matrix can't multiply." << endl;
        abort();
    }
    else
    {
        int m = M.size();
        int o = M[0].getDimension();
        int n = N[0].getDimension();

        vector<vector<int>> index_M(m, vector<int>(2, 0));
        vector<vector<int>> index_N(n, vector<int>(2, 0));
        Point temp(n);
        vector<Point> O(m, temp);

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < o; j++)
            {
                if (M[i].values[j] != 0)
                {
                    index_M[i][0] = j;
                    break;
                }
            }
            for (int j = o - 1; j >= 0; j--)
            {
                if (M[i].values[j] != 0)
                {
                    index_M[i][1] = j + 1;
                    break;
                }
            }
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < o; j++)
            {
                if (N[j].values[i] != 0)
                {
                    index_N[i][0] = j;
                    break;
                }
            }
            for (int j = o - 1; j >= 0; j--)
            {
                if (N[j].values[i] != 0)
                {
                    index_N[i][1] = j + 1;
                    break;
                }
            }
        }
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                int index = (index_M[i][1] < index_N[j][1]) ? (index_M[i][1]) : (index_N[j][1]);

                double sum = 0;

                for (int k = (index_M[i][0] > index_N[j][0]) ? (index_M[i][0]) : (index_N[j][0]); k < index; k++)
                {
                    sum += M[i].values[k] * N[k].values[j];
                }
                O[i].values[j] = sum;
            }
        }
        return O;
    }
}

vector<Point> Transpose(vector<Point> M)
{
    if (M.size() == 0)
    {
        vector<Point> temp(0);
        //cerr << "Matrix is empty." << endl;
        return temp;
    }
    int m = M.size();
    int n = M[0].getDimension();
    Point temp(m);
    vector<Point> O(n, temp);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            O[j].values[i] = M[i].values[j];
        }
    }

    return O;
}

vector<Point> Inverse(vector<Point> M)
{
    int m = M.size();

    if (m == 0)
    {
        vector<Point> temp(0);
        //cerr << "Matrix is empty." << endl;
        return temp;
    }

    Point temp(m);
    vector<Point> N(m, temp);

    if (m == 1)
    {
        N[0].values[0] = 1 / M[0].values[0];
    }
    else if (m >= 2)
    {

        vector<Point> T(m, temp);

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                T[i].values[j] = M[i].values[j];
            }
        }

        // make identity matrix
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                N[i].values[j] = (i == j) ? (1) : (0);
            }
        }

        // lower triangle elimination
        for (int k = 0; k < m - 1; k++)
        {
            for (int i = k + 1; i < m; i++)
            {
                double ratio = T[i].values[k] / T[k].values[k];

                for (int j = k; j < m; j++)
                {
                    T[i].values[j] -= T[k].values[j] * ratio;
                }
                for (int j = 0; j < m; j++)
                {
                    N[i].values[j] -= N[k].values[j] * ratio;
                }
            }
        }

        // make diagonal to 1.0
        for (int i = 0; i < m; i++)
        {
            double ratio = T[i].values[i];

            T[i].values[i] = 1.0;
            for (int j = i + 1; j < m; j++)
            {
                T[i].values[j] /= ratio;
            }
            for (int j = 0; j <= i; j++)
            {
                N[i].values[j] /= ratio;
            }
        }

        // upper triangle elimination
        for (int k = m - 1; k > 0; k--)
        {
            for (int i = k - 1; i >= 0; i--)
            {
                double ratio = T[i].values[k];

                T[i].values[k] = 0;
                for (int j = 0; j < m; j++)
                {
                    N[i].values[j] -= N[k].values[j] * ratio;
                }
            }
        }
    }

    return N;
}

vector<Point> matrixMinus(vector<Point> M, vector<Point> N)
{
    if (M.size() == 0 && N.size() == 0)
    {
        vector<Point> temp(0);
        return temp;
    }
    else if (M.size() == 0 && N.size() != 0)
    {
        int row = N.size();
        int col = N[0].getDimension();

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                N[i].values[j] *= -1;
            }
        }
        return N;
    }
    else if (M.size() != 0 && N.size() == 0)
    {
        return M;
    }
    else
    {
        int row = M.size();
        int col = M[0].getDimension();

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                M[i].values[j] -= N[i].values[j];
            }
        }
        return M;
    }
}

vector<Point> matrixAdd(vector<Point> M, vector<Point> N)
{
    if (M.size() == 0 && N.size() == 0)
    {
        vector<Point> temp(0);
        return temp;
    }
    else if (M.size() == 0 && N.size() != 0)
    {
        return N;
    }
    else if (M.size() != 0 && N.size() == 0)
    {
        return M;
    }
    else
    {
        int row = M.size();
        int col = M[0].getDimension();

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                M[i].values[j] += N[i].values[j];
            }
        }
        return M;
    }
}

double minValue(vector<Point> points){
    vector<double> min_vec;
    for(int i = 0; i < points.size(); i++){
        min_vec.push_back(*min_element(points[i].values.begin(), points[i].values.end()));
    }
    double min;
    min = *min_element(min_vec.begin(), min_vec.end());

    return min;
}

double maxValue(vector<Point> points){
    vector<double> max_vec;
    for(int i = 0; i < points.size(); i++){
        max_vec.push_back(*max_element(points[i].values.begin(), points[i].values.end()));
    }
    double max;
    max = *max_element(max_vec.begin(), max_vec.end());

    return max;
}

vector<double> dimMin(vector<Point> points){
    vector<double> min_vec;
    for(int j = 0; j < points[0].getDimension(); j++){
        vector<double> temp;
        for(int i = 0; i < points.size(); i++)
        {
            temp.push_back(points[i].values[j]);
        }
        min_vec.push_back(*min_element(temp.begin(), temp.end()));
    }
    return min_vec;
}

vector<double> dimMax(vector<Point> points){
    vector<double> max_vec;
    for(int j = 0; j < points[0].getDimension(); j++){
        vector<double> temp;
        for(int i = 0; i < points.size(); i++)
        {
            temp.push_back(points[i].values[j]);
        }
        max_vec.push_back(*max_element(temp.begin(), temp.end()));
    }
    return max_vec;
}
