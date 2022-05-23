#include "utils.h"
#include "pointUtils.h"
#include "Hungarian.h"
#include "metrics.h"
#include "methods.h"

// argv[1]->dataset, argv[2]->missing_ratio, argv[3]->executing times

int main(int argc, char *argv[])
{
	/* Setting parameters */
	string dataset = argv[1];												 
	string missing_ratio = argv[2];											 
	string times = argv[3];
	
	int max_miss = stoi(missing_ratio);
	int max_times = stoi(times);

	ofstream f("result/"+ dataset +"_Result.csv", ios::app);
	ofstream ACC_file("result/"+ dataset +"_acc_Result.csv", ios::app);
	ofstream NMI_file("result/"+ dataset +"_nmi_Result.csv", ios::app);
	ofstream FSCORE_file("result/"+ dataset +"_fscore_Result.csv", ios::app);

	ACC_file << ",";
	NMI_file << ",";
	FSCORE_file << ",";

	for (int i = 5; i <= max_miss; i += 5)
	{
		ACC_file << i << "%,";
		NMI_file << i << "%,";
		FSCORE_file << i << "%,";
	}

	ACC_file << "\n" << "CPIK" << ",";
	NMI_file << "\n" << "CPIK" << ",";
	FSCORE_file << "\n" << "CPIK" << ",";

	for (int missrate = 5; missrate <= max_miss; missrate += 5)
	{
		vector<double> ACCvec;
		vector<double> NMIvec;
		vector<double> FSCOREvec;

		for (int round = 1; round <= max_times; round++)
		{
			srand(round);
			/*======== Read data ========*/
			int dimension;
			int point_count;
			int cluster_count;
			vector<int> nan_index;
			vector<int> GT_vec;
			vector<Point> points;
			vector<Point> GTpoints;
			string label_filename = "data_set\\" + dataset + "\\" + dataset + "_label.txt";
			string miss_filename = "data_set\\" + dataset + "\\" + dataset + "_missing_" + to_string(missrate) + ".txt";
			string GT_filename = "data_set\\" + dataset + "\\" + dataset + "_dataset.txt";

			tie(point_count, dimension) = getDataSize(GT_filename);

			vector<vector<int>> miss_index(point_count, vector<int>(dimension, 0));

			cluster_count = readGT(GT_vec, label_filename, nan_index);
			points = readMissData(miss_index, miss_filename, nan_index);
			GTpoints = readGTData(GT_filename, nan_index);

			normalizePoints(points);

			/* Informations */
			cout << "=======================" << endl;
			cout << "Informations" << endl;
			cout << "-----------------------" << endl;
			cout << "Round: " << round << endl;
			cout << "Dataset: " << dataset << endl;
			cout << "Missing ratio: " << missrate << "%" << endl;
			cout << "Point#: " << point_count << endl;
			cout << "Dimension#: " << dimension << endl;
			cout << "Cluster#: " << cluster_count << endl;
			cout << "=======================" << endl;
			cout << "CPIK processing" << endl;
			cout << "-----------------------" << endl;

			/*======== CPIK process ========*/
			int start_idx;
			vector<int> clusters;
			vector<double> miss_count_vec;
			vector<Point> centroids;

			miss_count_vec = missCount(miss_index, dimension);
			start_idx = inputeMinMissColByGCS(miss_index, points, miss_count_vec, cluster_count);
			clusters = CPIK(points, cluster_count, miss_index, miss_count_vec, start_idx, centroids);

			if (clusters.size() == 0)
			{
				continue;
			}

			/*======== Calculate metrics ========*/
			double Acc, NMI, f_score, rmse;
			vector<int> adj_clusters;
			vector<int> assignment;
			vector<vector<double>> cost_matrix;

			cost_matrix = getCostMatrix(clusters, GT_vec);
			assignment = getAssignment(clusters, GT_vec, cost_matrix);
			adj_clusters = getAdjustClusters(clusters, assignment);

			Acc = getAcc(cost_matrix, assignment, point_count) * 100;
			NMI = getNMI(adj_clusters, GT_vec) * 100;
			f_score = getFScore2(adj_clusters, GT_vec) * 100;

			cout << "=======================" << endl;
			cout << "Evaluating metrics" << endl;
			cout << "-----------------------" << endl;
			cout << "Acc: " << Acc << "%" << endl
					<< "NMI: " << NMI << endl
					<< "F_score: " << f_score << endl;
			cout << "=======================" << endl;
			cout << endl;

			/*======== Write into CSV files ========*/
			f << "CPIK" << "," << missrate << "," << Acc << "," << NMI << "," << f_score << endl;

			ACCvec.push_back(Acc);
			NMIvec.push_back(NMI);
			FSCOREvec.push_back(f_score);

			
		}

		double ACC_avg = (!isnan(accumulate(ACCvec.begin(), ACCvec.end(), 0.0) / ACCvec.size())) ? accumulate(ACCvec.begin(), ACCvec.end(), 0.0) / ACCvec.size() : 0;
		double NMI_avg = (!isnan(accumulate(NMIvec.begin(), NMIvec.end(), 0.0) / NMIvec.size())) ? accumulate(NMIvec.begin(), NMIvec.end(), 0.0) / NMIvec.size() : 0;
		double Fscore_avg = (!isnan(accumulate(FSCOREvec.begin(), FSCOREvec.end(), 0.0) / FSCOREvec.size())) ? accumulate(FSCOREvec.begin(), FSCOREvec.end(), 0.0) / FSCOREvec.size() : 0;

		if (ACC_avg == 0)
			ACC_file << ",";
		else
			ACC_file << ACC_avg << ",";

		if (NMI_avg == 0)
			NMI_file << ",";
		else
			NMI_file << NMI_avg << ",";

		if (Fscore_avg == 0)
			FSCORE_file << ",";
		else
			FSCORE_file << Fscore_avg << ",";
	}

	ACC_file << endl;
	NMI_file << endl;
	FSCORE_file << endl;

	f.close();
	ACC_file.close();
	NMI_file.close();
	FSCORE_file.close();

	return 0;
}