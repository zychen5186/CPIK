#include "metrics.h"


double getAcc(vector<vector<double>> &maxCostMatrix, vector<int> &assignment, int size)
{
    double cost = 0;
    for (unsigned int x = 0; x < maxCostMatrix.size(); x++)
    {
        // std::cout << x << "," << assignment[x] << "\t";
        // cout << maxCostMatrix[x][assignment[x]] << endl;
        cost += maxCostMatrix[x][assignment[x]];
    }

    double accuracy = cost / size;

    return accuracy;
}

double getNMI(vector<int> &v1, vector<int> &v2)
{
    int total = v1.size();
    set<int> A_ids(v1.begin(), v1.end());
    set<int> B_ids(v2.begin(), v2.end());

    double MI = 0;
    double eps = 1e-45;
    int count;
    int idA, idB;
    vector<int> idApos, idBpos, idABpos;
    double px, py, pxy;
    // compare the position of classfied class between groud truth and result cluster vector
    for (set<int>::iterator it1 = A_ids.begin(); it1 != A_ids.end(); it1++)
    {
        for (set<int>::iterator it2 = B_ids.begin(); it2 != B_ids.end(); it2++)
        {
            idA = *it1;
            idB = *it2;
            idApos.clear();
            idBpos.clear();
            idABpos.clear();
            count = 0;
            for (vector<int>::iterator it = v1.begin(); it != v1.end(); it++)
            {
                if (*it == idA)
                    idApos.push_back(count);
                count++;
            }
            count = 0;
            for (vector<int>::iterator it = v2.begin(); it != v2.end(); it++)
            {
                if (*it == idB)
                    idBpos.push_back(count);
                count++;
            }
            set_intersection(idApos.begin(), idApos.end(), idBpos.begin(), idBpos.end(), back_inserter(idABpos));
            px = 1.0 * idApos.size() / total;
            py = 1.0 * idBpos.size() / total;
            pxy = 1.0 * idABpos.size() / total;
            MI += pxy * log2(pxy / (px * py) + eps);
        }
    }

    double Hcount = 0.0;
    double Hx = 0.0;
    for (int n : A_ids)
    {
        Hcount = 0.0;
        for (int m : v1)
        {
            if (n == m)
                Hcount++;
        }
        Hx -= (Hcount / total) * log2(Hcount / total + eps);
    }

    double Hy = 0.0;
    for (int n : B_ids)
    {
        Hcount = 0.0;
        for (int m : v2)
        {
            if (n == m)
                Hcount++;
        }
        Hy -= (Hcount / total) * log2(Hcount / total + eps);
    }
    // cout << MI << " " << Hx << " " << Hy <<endl;
    double MIhat = 0.0;
    MIhat = (2.0 * MI / (Hx + Hy));

    return MIhat;
}

double getNMI2(vector<int> &v1, vector<int> &v2)
{
    int total = v1.size();
    set<int> A_ids(v1.begin(), v1.end());
    set<int> B_ids(v2.begin(), v2.end());
    int nClass = B_ids.size();
    vector<vector<double>> G(nClass, vector<double>(nClass, 0.0));
    int x = 0;
    int y = 0;
    for (set<int>::iterator it1 = A_ids.begin(); it1 != A_ids.end(); it1++)
    {
        for (set<int>::iterator it2 = B_ids.begin(); it2 != B_ids.end(); it2++)
        {
            for (int i = 0; i < total; i++)
            {
                if (v1[i] == *it1 && v2[i] == *it2)
                {
                    G[x][y]++;
                }
            }
            y++;
        }
        x++;
        y = 0;
    }
    double sumG = 0;
    for (auto &i : G)
    {
        for (auto &j : i)
        {
            sumG += j;
        }
    }
    vector<double> P1(nClass, 0.0);
    vector<double> P2(nClass, 0.0);
    for (int i = 0; i < nClass; i++)
    {
        for (int j = 0; j < nClass; j++)
        {
            P1[i] += G[i][j] / sumG;
            P2[j] += G[i][j] / sumG;
        }
    }
    double H1 = 0;
    double H2 = 0;
    for (int i = 0; i < nClass; i++)
    {
        H1 += (-P1[i] * log2(P1[i]));
        H2 += (-P2[i] * log2(P2[i]));
    }
    vector<vector<double>> P12(nClass, vector<double>(nClass, 0.0));
    vector<vector<double>> PPP(nClass, vector<double>(nClass, 0.0));
    double MI = 0;
    for (int i = 0; i < nClass; i++)
    {
        for (int j = 0; j < nClass; j++)
        {
            P12[i][j] = G[i][j] / sumG;
            PPP[i][j] = P12[i][j] / P2[j] / P1[i];
            if (abs(PPP[i][j]) < 1e-12)
                PPP[i][j] = 1;
            MI += P12[i][j] * log2(PPP[i][j]);
        }
    }
    double MIhat = MI / max(H1, H2);

    return MIhat;
}

double getFScore(vector<vector<double>> &maxCostMatrix, vector<int> &assignment)
{

    double eps = 1e-45;
    int costMatrixSize = maxCostMatrix[0].size();
    double precision, recall, TP, FP, FN;
    double f_score = 0;
    // for(int i = 0; i< costMatrixSize; i++){
    //     cout << i << " " << assignment[i] <<endl;}
    for (int i = 0; i < costMatrixSize; i++) // i 是預測後編號, adjust_i是做mapping後調整的編號
    {
        int adjust_i = assignment[i];
        FP = 0;
        FN = 0;
        TP = maxCostMatrix[i][adjust_i];
        for (int j = 0; j < costMatrixSize; j++)
        {
            if (j != adjust_i)
            {
                FP += maxCostMatrix[i][j];
            }
            if (j != i)
            {
                FN += maxCostMatrix[j][adjust_i];
            }
        }
        // cout << TP << " " << FP << " " << FN << endl;
        precision = TP / (TP + FP + eps);
        recall = TP / (TP + FN + eps);
        f_score += 2 * precision * recall / (precision + recall + eps);
        // cout << precision << " " << recall<< " " << f_score <<endl;
    }

    f_score /= costMatrixSize;

    return f_score;
}

double getFScore2(vector<int> &v1, vector<int> &v2)
{
    int total = v1.size();
    set<int> p(v1.begin(), v1.end());
    set<int> c(v2.begin(), v2.end());
    // int C_size = p.size();
    int C_size = c.size();
    Point temp(total);
    vector<Point> Pid(C_size, temp);
    vector<Point> Cid(C_size, temp);
    for (int i = 0; i < total; i++)
    {
        int Pidx = v1[i];
        Pid[Pidx].values[i] = 1;
        int Cidx = v2[i];
        Cid[Cidx].values[i] = 1;
    }
    vector<Point> CP = Multiplication(Cid, Transpose(Pid));
    vector<int> Pj(C_size, 0);
    vector<int> Ci(C_size, 0);
    double sumPj = 0;
    for (int i = 0; i < C_size; i++)
    {
        for (int j = 0; j < C_size; j++)
        {
            Pj[j] += CP[i].values[j];
            sumPj += CP[i].values[j];
            Ci[i] += CP[i].values[j];
        }
    }
    Point temp2(C_size);
    vector<vector<double>> precision(C_size, vector<double>(C_size));
    vector<vector<double>> recall(C_size, vector<double>(C_size));
    vector<Point> F(C_size, temp2);
    vector<double> maxF(C_size, numeric_limits<double>::min());
    for (int i = 0; i < C_size; i++)
    {
        for (int j = 0; j < C_size; j++)
        {
            precision[i][j] = CP[i].values[j] / Ci[i];
            recall[i][j] = CP[i].values[j] / Pj[j];
            F[i].values[j] = 2 * precision[i][j] * recall[i][j] / (precision[i][j] + recall[i][j]);
        }
    }
    F = Transpose(F);
    for (int i = 0; i < C_size; i++)
    {
        for (int j = 0; j < C_size; j++)
        {
            if (F[i].values[j] > maxF[j] && !isnan(F[i].values[j]))
            {
                maxF[j] = F[i].values[j];
            }
        }
    }
    double FMeasure = 0;

    for (int i = 0; i < C_size; i++)
    {
        FMeasure += Pj[i] / sumPj * maxF[i];
    }

    return FMeasure;
}

double RMSE(vector<Point> &v1, vector<Point> &v2)
{
    int point_count = v1.size();
    int dimension = v1[0].getDimension();
    if(v1[0].getDimension() != v2[0].getDimension() || v1.size() != v2.size()){
        cout << v1.size() << " " << v2.size() <<endl;
        cout << "dimension or size is not the same!" <<endl;
    }
    double rmse = 0;
    for (int i = 0; i < point_count; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            rmse += pow(v1[i].values[j] - v2[i].values[j], 2);
        }
    }
    rmse /= point_count;
    //rmse /= dimension;
    rmse = sqrt(rmse);
    return rmse;
}
