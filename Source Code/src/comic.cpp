#include <unordered_set>
//#include <windows.h>
#include "mc.h"
#include "memoryusage.h"
#include "anyoption.h"
#include "timgraph.h"

//choose dataset
string database = "/lastfm_wc.txt";
AnyOption* read_options(int argc, char* argv[])
{
	// read the command line options
	AnyOption* opt = new AnyOption();

	// ignore POSIX style options
	opt->noPOSIX();
	opt->setVerbose(); /* print warnings about unknown options */
	opt->autoUsagePrint(true); /* print usage for bad options */

	/* 3. SET THE USAGE/HELP   */
	opt->addUsage("");
	opt->addUsage("Usage: ");
	opt->addUsage("");
	opt->addUsage(" -help: print this help ");
	opt->addUsage(" -c <config_file>: specify config file ");
	opt->addUsage(" -dataset <path_to_dataset>: specify the dataset name");
	opt->addUsage(" -aseeds <path_to_A_seeds>: specify the path to A seed set file");
	opt->addUsage(" -bseeds <path_to_B_seeds>: specify the path to B seed set file");
	opt->addUsage("");

	/* 4. SET THE OPTION STRINGS/CHARACTERS */
	opt->setOption("phase");
	opt->setCommandOption("dataset");
	opt->setOption("aseeds");
	opt->setOption("qao");
	opt->setOption("qab");
	opt->setOption("qbo");
	opt->setOption("qba");
	opt->setOption("kA");
	opt->setOption("kB");
	opt->setOption("selectBSeeds");
	opt->setOption("ignoreBSeeds");
	opt->setOption("bseeds");
	opt->setOption("algo");
	opt->setOption("epsilon");
	opt->setOption("pr");
	opt->setOption("extent");
	opt->setOption("theta");
	opt->setOption("runs");
	/* for options that will be checked only on the command line not in option/resource file */
	opt->setCommandFlag("help");
	opt->setCommandOption("c");

	/* go through the command line and get the options  */
	opt->processCommandArgs(argc, argv);

	/* 6. GET THE VALUES */
	if (opt->getFlag("help")) {
		opt->printUsage();
		delete opt;
		exit(0);
	}

	//const char *configFile = opt->getValue("c");
	//绝对路径
	const char* configFile = "config.txt";

	if (configFile == NULL) {
		cout << "[error] config file not specified!" << endl;
		opt->printUsage();
		delete opt;
		exit(1);
	}

	opt->processFile(configFile);
	opt->processCommandArgs(argc, argv);

	cout << "[info] config file processed: " << configFile << endl;
	return opt;
}


void select_B_seeds(MonteCarlo& mc, string folder, int kB, int n, bool isOverlap)
{
	srand(static_cast<unsigned int>(time(nullptr)));
	unordered_set<int> bset;
	bset.clear();
	while (bset.size() < kB) {
		int node = rand() % n;
		//if (mc.graph->group[node] == BG || ((mc.graph->group[node] == AB) && isOverlap)) {
			bset.insert(node);
		//}

	}
	ofstream my_file;
	string filename = folder + "/B-SEED/BS1/bseeds_" + to_string(kB) + ".txt";
	my_file.open(filename.c_str());
	for (auto it = bset.begin(); it != bset.end(); ++it) {
		my_file << *it << endl;
	}
	my_file.close();
}


void select_A_seeds(MonteCarlo& mc, string folder, int kA, int n, bool isOverlap)
{


	unordered_set<int> aset;
	aset.clear();
	srand(static_cast<unsigned int>(time(nullptr)));
	while (aset.size() < kA) {
		
		int node = rand() % n;
		if (node < 0 || node >= mc.n) {
			cerr << "Error: node index out of bounds: " << node << endl;
			continue;
		}
		//if (mc.graph->group[node] == AG || (mc.graph->group[node] == AB && isOverlap)) {
			aset.insert(node);
		//}
	}

	ofstream my_file;
	string filename = folder + "/A-SEED/AS1/aseeds_" + to_string(kA) + ".txt";
	my_file.open(filename.c_str());
	for (auto it = aset.begin(); it != aset.end(); ++it) {
		my_file << *it << endl;
	}
	my_file.close();
}




void computeTrueSpreadByMC(int k, int* seedSet, string dataset, vector<double> qq, string bseeds, bool ignore_B)
{
	cout << "[info] compute true spread using Monte Carlo (SelfInfMax)" << endl;
	MonteCarlo mc(dataset, database);
	mc.setParametersMC(k, qq, ignore_B, dataset, bseeds);
	mc.readBSeedsMC();
	printf("\n");
	for (int i = 0; i < k; i++) {
		if (i == 0 || (i + 1) % 10 == 0) {
			double cov = mc.compute_coverage(seedSet, i + 1);
			cout << i << "\t" << cov << endl;
		}
	}
}
void computeTrueNodeCnt(int k, string aseeds_file_name, string dataset, vector<double> qq, string bseeds, bool ignore_B) {
	cout << "[info] compute true nodeCount using Monte Carlo (Seq)" << endl;
	MonteCarlo mc(dataset, database);
	mc.setParametersMC(k, qq, ignore_B, dataset, bseeds);
	mc.readASeedsMC(aseeds_file_name);
	mc.readBSeedsMC();
	srand(time(NULL));

	double cov1 = mc.compute_coverage_seq_Afirst_nodeCount(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size(),false);
	double cov2 = mc.compute_coverage_seq_Bfirst_nodeCount(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size(),false);

	double cov3 = mc.compute_coverage_seq_Afirst_nodeCount(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size(), true);
	double cov4 = mc.compute_coverage_seq_Bfirst_nodeCount(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size(), true);

	double high_cov = max(cov1, cov2);
	double low_cov = min(cov1, cov2);
	cout << endl;
	cout << "A before B nodeCnt: " << cov1 << endl;
	cout << "B before A nodeCnt: " << cov2 << endl;
	cout << "Influence gap : [" << cov1 - cov2 << "]" << endl;
	double ratio = (high_cov - low_cov) / (low_cov);
	cout << "Gap ratio :" << ratio * 100 << "%" << endl;
	cout << "A before B nodeCnt q1: " << cov3 << endl;
	cout << "B before A nodeCnt q1: " << cov4 << endl;
	cout << "SA factor 1:" << (cov1 / cov3)*100 << "%" <<endl;
	cout << "SA factor 2:" << (cov2 / cov4)*100 << "%" << endl;
	cout << endl;
	
}
void computeTrueAdoption_upper(int k, string aseeds_file_name, string dataset, vector<double> qq, string bseeds, bool ignore_B) {
	cout << "[info] compute true adoption using Monte Carlo (Seq), upper bound" << endl;
	MonteCarlo mc(dataset, database);
	mc.setParametersMC(k, qq, ignore_B, dataset, bseeds);
	double cov1 = 0, cov2 = 0, cov3 = 0, cov4 = 0;
	double Acover1 = 0, Acover2 = 0, Acover3 = 0, Acover4 = 0;
	double Bcover1 = 0, Bcover2 = 0, Bcover3 = 0, Bcover4 = 0;
	double cnt_qab1 = 0, cnt_qba2 = 0, cnt_qab3 = 0, cnt_qba4 = 0;
	double cnt_qba1 = 0, cnt_qab2 = 0, cnt_qba3 = 0, cnt_qab4 = 0;
	int runs = 1;
	//srand(static_cast<unsigned int>(time(nullptr)));
	for (int i = 0; i < runs; i++) {
		cout << endl;
		cout << "[runs:" << i << "]" << endl;

		mc.readASeedsMC(aseeds_file_name);

		mc.readBSeedsMC();
		srand(time(NULL));
		cov1 += mc.compute_coverage_seq_Afirst(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size(), Acover1, Bcover1, cnt_qab1, cnt_qba1);
		cov2 += mc.compute_coverage_seq_Bfirst(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size(), Acover2, Bcover2, cnt_qab2, cnt_qba2);

		cov3 += mc.compute_coverage_seq_Afirst_upper(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size());
		cov4 += mc.compute_coverage_seq_Bfirst_upper(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size());

	}

	cov1 /= runs;
	cov2 /= runs;
	cov3 /= runs;
	cov4 /= runs;


	double high_cov = max(cov1, cov2);
	double low_cov = min(cov1, cov2);
	cout << endl;
	cout << "A before B AdopCnt: " << cov1 << endl;
	cout << "B before A AdopCnt: " << cov2 << endl;
	cout << "Influence gap : [" << cov1 - cov2 << "]" << endl;
	double ratio = (high_cov - low_cov) / (low_cov);
	cout << "Gap ratio :" << ratio * 100 << "%" << endl;
	cout << "A before B AdopCnt q+: " << cov3 << endl;
	cout << "B before A AdopCnt q+: " << cov4 << endl;
	cout << "SA factor 1:" << (cov1 / cov3) * 100 << "%" << endl;
	cout << "SA factor 2:" << (cov2 / cov4) * 100 << "%" << endl;
	cout << endl;

}
void computeTrueAdoption(int k, string aseeds_file_name, string dataset, vector<double> qq, string bseeds, bool ignore_B) {
	cout << "[info] compute true adoption using Monte Carlo (Seq)" << endl;
	MonteCarlo mc(dataset, database);
	mc.setParametersMC(k, qq, ignore_B, dataset, bseeds);
	double cov1 = 0, cov2 = 0, cov3 = 0, cov4 = 0;
	double Acover1 = 0, Acover2 = 0, Acover3 = 0, Acover4 = 0;
	double Bcover1 = 0, Bcover2 = 0, Bcover3 = 0, Bcover4 = 0;
	double cnt_qab1 = 0, cnt_qba2 = 0, cnt_qab3 = 0, cnt_qba4 = 0;
	double cnt_qba1 = 0, cnt_qab2 = 0, cnt_qba3 = 0, cnt_qab4 = 0;
	int runs = 1;
	//srand(static_cast<unsigned int>(time(nullptr)));
	for (int i = 0; i < runs; i++) {
		cout << endl;
		cout << "[runs:" << i << "]" << endl;
		mc.readASeedsMC(aseeds_file_name);
		mc.readBSeedsMC();
		srand(time(NULL));
		cov1 += mc.compute_coverage_seq_Afirst(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size(), Acover1, Bcover1, cnt_qab1, cnt_qba1);
		cov2 += mc.compute_coverage_seq_Bfirst(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size(), Acover2, Bcover2, cnt_qab2, cnt_qba2);
		//cov3 += mc.compute_coverage_simu_v1(mc.aSeeds, k, mc.bSeeds, mc.bSeeds.size(), Acover3, Bcover3, cnt_qab3, cnt_qba3);
		//cov4 += mc.compute_coverage_simu_v2(mc.aSeeds, k, mc.bSeeds, k, Acover4, Bcover4, cnt_qab4, cnt_qba4);

	}

	cov1 /= runs;
	cov2 /= runs;
	//cov3 /= runs;
	//cov4 /= runs;

	double high_cov = max(cov1, cov2);
	double low_cov = min(cov1, cov2);
	cout << endl;
	cout << "A before B coverage: " << cov1 << endl;
	cout << "B before A coverage: " << cov2 << endl;
	cout << "Influence gap : [" << cov1 - cov2 << "]" << endl;
	double ratio = (high_cov - low_cov) / (low_cov);
	cout << "Gap ratio :" << ratio * 100 << "%" << endl;
	cout << endl;
	/*
	cout << "A simutaneous with B coverage v1: " << cov3 << endl;
	cout << " Influence gap : [" << high_cov - cov3 << "]" << endl;
	ratio = (high_cov - cov3) / (cov3);
	cout << "Gap ratio :" << ratio * 100 << "%" << endl;
	cout << endl;
	
	cout << "A simutaneous with B coverage v2: " << cov4 << endl;
	cout << " Influence gap : [" << high_cov - cov4 << "]" << endl;
	ratio = (high_cov - cov4) / (cov4);
	cout << "Gap ratio :" << ratio * 100 << "%" << endl;
	cout << endl;*/

	cout << "A first, Item A coverage: " << Acover1 / (double)runs << "  Item B coverage: " << Bcover1 / (double)runs << endl;
	cout << "[Couting] node calculate A-adoption with probability qab for  " << cnt_qab1 / (double)runs << " times" << endl;
	cout << "[Couting] node calculate B-adoption with probability qba for  " << cnt_qba1 / (double)runs << " times" << endl;
	cout << "_____________________________________________________________________________________________________________" << endl;
	cout << "B first, Item A coverage: " << Acover2 / (double)runs << "  Item B coverage: " << Bcover2 / (double)runs << endl;
	cout << "[Couting] node calculate A-adoption with probability qab for  " << cnt_qab2 / (double)runs << " times" << endl;
	cout << "[Couting] node calculate B-adoption with probability qba for  " << cnt_qba2 / (double)runs << " times" << endl;
	cout << "_____________________________________________________________________________________________________________" << endl;
	cout << "Simultaneous, Item A coverage: " << Acover3 / (double)runs << "  Item B coverage: " << Bcover3 / (double)runs << endl;
	cout << "[Couting] node calculate A-adoption with probability qab for  " << cnt_qab3 / (double)runs << " times" << endl;
	cout << "[Couting] node calculate B-adoption with probability qba for  " << cnt_qba3 / (double)runs << " times" << endl;
	cout << "_____________________________________________________________________________________________________________" << endl;
	
}
void computeTrueSpread(vector<vector<int>> vec_aseed, vector<vector<int>> vec_bseed, int runs, int k, string dataset, vector<double> qq) {
	MonteCarlo mc(dataset, database);
	mc.setParametersMC(k, qq, 0, dataset, "null");
	double cov1 = 0, cov2 = 0, cov3 = 0;
	double Acover1 = 0, Acover2 = 0, Acover3 = 0;
	double Bcover1 = 0, Bcover2 = 0, Bcover3 = 0;
	double cnt_qab1 = 0, cnt_qba2 = 0, cnt_qab3 = 0;
	double cnt_qba1 = 0, cnt_qab2 = 0, cnt_qba3 = 0;
	//srand(static_cast<unsigned int>(time(nullptr)));
	for (int i = 0; i < runs; i++) {
		cout << endl;
		cout << "[runs:" << i << "]" << endl;
		vector<int> aseed = vec_aseed[i];
		vector<int> bseed = vec_bseed[i];
		cov1 += mc.compute_coverage_seq_Afirst(aseed, aseed.size(), bseed, bseed.size(), Acover1, Bcover1, cnt_qab1, cnt_qba1);
		cov2 += mc.compute_coverage_seq_Bfirst(aseed, aseed.size(), bseed, bseed.size(), Acover2, Bcover2, cnt_qab2, cnt_qba2);
		cov3 += mc.compute_coverage_simu_v1(aseed, aseed.size(), bseed, bseed.size(), Acover3, Bcover3, cnt_qab3, cnt_qba3);

	}

	cov1 /= runs;
	cov2 /= runs;
	cov3 /= runs;

	double high_cov = max(cov1, cov2);
	double low_cov = min(cov1, cov2);
	cout << endl;
	cout << "A before B coverage: " << cov1 << endl;
	cout << "B before A coverage: " << cov2 << endl;
	cout << "Influence gap : [" << cov1 - cov2 << "]" << endl;
	double ratio = (high_cov - low_cov) / (low_cov);
	cout << "Gap ratio :" << ratio * 100 << "%" << endl;
	cout << endl;
	cout << "A simutaneous with B coverage: " << cov3 << endl;
	cout << " Influence gap : [" << high_cov - cov3 << "]" << endl;
	ratio = (high_cov - cov3) / (cov3);
	cout << "Gap ratio :" << ratio * 100 << "%" << endl;
	cout << endl;

	cout << "A ahead, Item A coverage: " << Acover1 / (double)runs << ",  Item B coverage: " << Bcover1 / (double)runs << endl;
	cout << "[Couting] node calculate A-adoption with probability qab for  " << cnt_qab1 / (double)runs << " times" << endl;
	cout << "[Couting] node calculate B-adoption with probability qba for  " << cnt_qba1 / (double)runs << " times" << endl;
	cout << "_____________________________________________________________________________________________________________" << endl;
	cout << "B ahead, Item A coverage: " << Acover2 / (double)runs << ",  Item B coverage: " << Bcover2 / (double)runs << endl;
	cout << "[Couting] node calculate A-adoption with probability qab for  " << cnt_qab2 / (double)runs << " times" << endl;
	cout << "[Couting] node calculate B-adoption with probability qba for  " << cnt_qba2 / (double)runs << " times" << endl;
	cout << "_____________________________________________________________________________________________________________" << endl;
	cout << "Simultaneously Item A coverage: " << Acover3 / (double)runs << "  Item B coverage: " << Bcover3 / (double)runs << endl;
	cout << "[Couting] node calculate A-adoption with probability qab for  " << cnt_qab3 / (double)runs << " times" << endl;
	cout << "[Couting] node calculate B-adoption with probability qba for  " << cnt_qba3 / (double)runs << " times" << endl;
	cout << "_____________________________________________________________________________________________________________" << endl;
}

int main(int argc, char** argv)
{
	AnyOption* opt = read_options(argc, argv);

	//string dataset = string(opt->getValue("dataset"));
	string dataset = "dataset";
	
	int kA = atoi(opt->getValue("kA"));
	int kB = atoi(opt->getValue("kB"));
	//GAPs
	double qao, qab, qbo, qba;
	qao = atof(opt->getValue("qao"));
	qab = atof(opt->getValue("qab"));
	qbo = atof(opt->getValue("qbo"));
	qba = atof(opt->getValue("qba"));
	vector<double> qq;
	qq.push_back(qao);
	qq.push_back(qab);
	qq.push_back(qbo);
	qq.push_back(qba);

	bool ignore_B = false, select_B = false;
	ignore_B = (atoi(opt->getValue("ignoreBSeeds"))) == 1 ? true : false;
	select_B = (atoi(opt->getValue("selectBSeeds"))) == 1 ? true : false;
	//string bseeds = "dataset/B-SEED/bseeds_"+to_string(kB) + ".txt";
	string bseeds = "dataset/B-SEED/BS2/bseeds_" + to_string(kB) + ".txt";
	int algo = atoi(opt->getValue("algo"));
	int phase = atoi(opt->getValue("phase"));
	int64 theta = atoi(opt->getValue("theta"));
	int runs = atoi(opt->getValue("runs"));

	cout << "[param] *** dataset = " << database << " ***" << endl;
	cout << "[param] kA = " << kA << ", kB = " << kB << endl;
	cout << "[param] Q-probability : " << qq[0] << ", " << qq[1] << ", " << qq[2] << ", " << qq[3] << endl;
	cout << "[param] ignore_B = " << ignore_B << endl;
	cout << "[param] select_B = " << select_B << endl;
	cout << "[param] Algorithm = " << algo << "| phase = " << phase << endl;
	cout << "[param] Theta = " << theta  << endl;

	clock_t begin;
	clock_t end;


	//compute adoption count in order
	if (algo == 1) {
		cout << "[info] algorithm: (105) reading seeds A and B for adoption cnt" << endl;
		string aseeds_file_name = "dataset/A-SEED/AS2/aseeds_" + to_string(kA) + ".txt";
		begin = clock();
		computeTrueAdoption(kA, aseeds_file_name, dataset, qq, bseeds, ignore_B);
		end = clock();
	}
	//compute node count in order
	if (algo == 2) {
		cout << "[info] algorithm: (106) reading seeds A and B for node cnt" << endl;
		string aseeds_file_name = "dataset/A-SEED/AS2/aseeds_" + to_string(kA) + ".txt";
		begin = clock();
		computeTrueNodeCnt(kA, aseeds_file_name, dataset, qq, bseeds, ignore_B);
		end = clock();
	}
	//mine seeds using IMM-q1
	if (algo == 3) {
		cout << "[info] algorithm: (1) IMM-Order Node Count q=1" << endl;
		string aseeds_file_name = "dataset/A-SEED/ASRR_test/aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/B-SEED/BSRR_test/bseeds_" + to_string(kB) + ".txt";
		TimGraph tim(dataset, database, theta);
		tim.setParametersIG(kA, kB, qq, ignore_B, dataset, bseeds);
		double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;

		tim.minseeds_Order_q1(epsilon, aseeds_file_name, bseeds_file_name);

		cout << "_________________________________________________________________" << endl;
		cout << "[info] Estimate true spread......" << endl;
		computeTrueAdoption(kA, aseeds_file_name, dataset, qq, bseeds_file_name, ignore_B);
	}
	//mine seeds using IMM-NoCompe
	if (algo == 4) {
		cout << "[info] algorithm: (2) IMM-independent Adoption Count" << endl;
		string aseeds_file_name = "dataset/A-SEED/ASRR_test/aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/B-SEED/BSRR_test/bseeds_" + to_string(kB) + ".txt";
		TimGraph tim(dataset, database, theta);
		tim.setParametersIG(kA, kB, qq, ignore_B, dataset, bseeds);
		double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;

		tim.minseeds_idp(epsilon, aseeds_file_name, bseeds_file_name);

		cout << "_________________________________________________________________" << endl;
		cout << "[info] Estimate true spread......" << endl;
		computeTrueAdoption(kA, aseeds_file_name, dataset, qq, bseeds_file_name, ignore_B);
	}
	//Estimate SA factor，node count
	if (algo == 5) {
		cout << "[info] algorithm: (3)Estimate SA factor Node Count q=1" << endl;
		string aseeds_file_name = "dataset/A-SEED/ASRR_test/aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/B-SEED/BSRR_test/bseeds_" + to_string(kB) + ".txt";
		TimGraph tim(dataset, database, theta);
		tim.setParametersIG(kA, kB, qq, ignore_B, dataset, bseeds);
		double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;

		tim.minseeds_Order_q1(epsilon, aseeds_file_name, bseeds_file_name);

		cout << "_________________________________________________________________" << endl;
		cout << "[info] Estimate SA factor......" << endl;
		computeTrueNodeCnt(kA, aseeds_file_name, dataset, qq, bseeds_file_name, ignore_B);
	}
	//Estimate SA factor， adoption count
	if (algo == 6) {
		cout << "[info] algorithm: (3)Estimate SA factor Adoption Count q+" << endl;
		string aseeds_file_name = "dataset/A-SEED/ASRR_test/aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/B-SEED/BSRR_test/bseeds_" + to_string(kB) + ".txt";
		TimGraph tim(dataset, database, theta);
		tim.setParametersIG(kA, kB, qq, ignore_B, dataset, bseeds);
		double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;

		tim.minseeds_idp(epsilon, aseeds_file_name, bseeds_file_name);

		cout << "_________________________________________________________________" << endl;
		cout << "[info] Estimate SA factor......" << endl;
		computeTrueAdoption_upper(kA, aseeds_file_name, dataset, qq, bseeds_file_name, ignore_B);
	}
	//mine seeds using OPPRT based on A frist
	if (algo == 7) {
		cout << "[info] algorithm: (4) RR-Order Afirst" << endl;
		string aseeds_file_name = "dataset/A-SEED/ASRR_test/aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/B-SEED/BSRR_test/bseeds_" + to_string(kB) + ".txt";
		TimGraph tim(dataset, database,theta);

		tim.setParametersIG(kA,kB,qq, ignore_B, dataset, bseeds);
		double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;
		
		vector<vector<vector<int>>> results = tim.mineSeeds_Order_Afirst(runs, epsilon, aseeds_file_name, bseeds_file_name);
		
		cout << "_________________________________________________________________" << endl;
		cout << "[info] Estimate true spread......" << endl;
		computeTrueSpread(results[0],results[1],runs,kA, dataset, qq);
		
	}
	//mine seeds using OPPRT based on B frist
	if (algo == 8) {
		cout << "[info] algorithm: (8) RR-Order Bfirst" << endl;
		string aseeds_file_name = "dataset/A-SEED/ASRR_test/aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/B-SEED/BSRR_test/bseeds_" + to_string(kB) + ".txt";
		TimGraph tim(dataset, database, theta);

		tim.setParametersIG(kA, kB, qq, ignore_B, dataset, bseeds);
		double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;
		int runs = 1;

		vector<vector<vector<int>>> results = tim.mineSeeds_Order_Bfirst(runs, epsilon, aseeds_file_name, bseeds_file_name);

		cout << "_________________________________________________________________" << endl;
		cout << "[info] Estimate true spread......" << endl;
		computeTrueSpread(results[0], results[1], runs, kA, dataset, qq);
	}
	//mine seeds using Greedy,A first
	if (algo == 9) {
		cout << "[info] algorithm: (5) Greedy mining SA and SB for influence with order.A first!" << endl;
		string aseeds_file_name = "dataset/A-SEED/ASGreedy/Afirst_aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/B-SEED/BSGreedy/Afirst_bseeds_" + to_string(kB) + ".txt";

		MonteCarlo mc(dataset, database);
		mc.setParametersMC(kA, qq, ignore_B, dataset, bseeds);

		
		cout << endl << "[info] phase: mining seeds for A and B ..." << endl;
		begin = clock();
		mc.mineABseed_Greedy_Afirst(kA, kB, aseeds_file_name, bseeds_file_name);
		end = clock();
		
		cout << "_________________________________________________________________" << endl;

		cout << "[info] Estimate true node Count......" << endl;
		double cov1 = mc.compute_coverage_seq_Afirst_nodeCount(mc.aSeeds, mc.aSeeds.size(), mc.bSeeds, mc.bSeeds.size(),false);
		double cov2 = mc.compute_coverage_seq_Bfirst_nodeCount(mc.aSeeds, mc.aSeeds.size(), mc.bSeeds, mc.bSeeds.size(),false);
		double high_cov = max(cov1, cov2);
		double low_cov = min(cov1, cov2);
		cout << endl;
		cout << "A before B nodeCnt: " << cov1 << endl;
		cout << "B before A nodeCnt: " << cov2 << endl;
		cout << "Influence gap : [" << cov1 - cov2 << "]" << endl;
		double ratio = (high_cov - low_cov) / (low_cov);
		cout << "Gap ratio :" << ratio * 100 << "%" << endl;
		cout << endl;

		cout << "[info] Estimate true adoption Count......" << endl;
		computeTrueAdoption(kA, aseeds_file_name, dataset, qq, bseeds_file_name, ignore_B);

	}
	//mine seeds using OPPRT based on B frist
	if (algo == 10) {
		cout << "[info] algorithm: (6) Greedy mining SA and SB for influence with order. B first!" << endl;
		string aseeds_file_name = "dataset/A-SEED/ASGreedy/Bfirst_aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/B-SEED/BSGreedy/Bfirst_bseeds_" + to_string(kB) + ".txt";

		MonteCarlo mc(dataset, database);
		mc.setParametersMC(kA, qq, ignore_B, dataset, bseeds);

		cout << endl << "[info] phase: mining seeds for A and B ..." << endl;
		begin = clock();
		mc.mineABseed_Greedy_Bfirst(kA, kB, aseeds_file_name, bseeds_file_name);
		end = clock();


		cout << "_________________________________________________________________" << endl;
		cout << "[info] Estimate true spread......" << endl;
		
		double cov1 = mc.compute_coverage_seq_Afirst_nodeCount(mc.aSeeds, mc.aSeeds.size(), mc.bSeeds, mc.bSeeds.size(),false);
		double cov2 = mc.compute_coverage_seq_Bfirst_nodeCount(mc.aSeeds, mc.aSeeds.size(), mc.bSeeds, mc.bSeeds.size(),false);
		double high_cov = max(cov1, cov2);
		double low_cov = min(cov1, cov2);
		cout << endl;
		cout << "A before B nodeCnt: " << cov1 << endl;
		cout << "B before A nodeCnt: " << cov2 << endl;
		cout << "Influence gap : [" << cov1 - cov2 << "]" << endl;
		double ratio = (high_cov - low_cov) / (low_cov);
		cout << "Gap ratio :" << ratio * 100 << "%" << endl;
		cout << endl;

	}

	//mine seeds using high degree
	if (algo == 11) {
		cout << "[info] algorithm: Highest Degree for both SIM and CIM!" << endl;
		string aseeds_file_name = "dataset/A-SEED/ASHighD/aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/B-SEED/BSHighD/bseeds_" + to_string(kB) + ".txt";
		MonteCarlo mc(dataset, database);
		mc.setParametersMC(kA, qq, ignore_B, dataset, bseeds_file_name);
		// find high-degree seeds!	    
		begin = clock();
		vector<int> vec = mc.graph->findHighDegree(kA+kB);
		end = clock();
		// spread for SIM
        vector<int> aSeeds;
        vector<int> bSeeds;
		int ca = 0, cb = 0;
        for (int i = 0; i < kA+kB; i++) {
            if ((i % 2 == 0 && ca<kA)||cb==kB) {
                aSeeds.push_back(vec[i]);
				ca++;
            } else {
                bSeeds.push_back(vec[i]);
				cb++;;
            }
        }
		mc.write_seeds_to_file(aSeeds, aseeds_file_name);
		mc.write_seeds_to_file(bSeeds, bseeds_file_name);
		cout << "[info] high degree seeds for A and B are written to files." << endl;
		cout << "_________________________________________________________________" << endl;
		cout << "[info] Estimate true spread......" << endl;
		computeTrueAdoption(kA, aseeds_file_name, dataset, qq, bseeds_file_name, ignore_B);

	}

	// mine seeds using page rank
	if (algo == 12)
	{
		cout << "[info] algorithm: PageRank for both CIM and SIM!" << endl;
		string aseeds_file_name = "dataset/dm/ASEED/PR/aseeds_" + to_string(kA) + ".txt";
		string bseeds_file_name = "dataset/dm/BSEED/PR/bseeds_" + to_string(kB) + ".txt";
		MonteCarlo mc(dataset, database);
		mc.setParametersMC(kA, qq, ignore_B, dataset, bseeds_file_name);
		begin = clock();
		vector<int> prSeeds = mc.graph->findPageRank(kA + kB);
		end = clock();

		vector<int> aSeeds;
		vector<int> bSeeds;
		int ca = 0, cb = 0;
		for (int i = 0; i < kA + kB; i++) {
			if ((i % 2 == 0 && ca < kA) || cb == kB) {
				aSeeds.push_back(prSeeds[i]);
				ca++;
			}
			else {
				bSeeds.push_back(prSeeds[i]);
				cb++;;
			}
		}
		mc.write_seeds_to_file(aSeeds, aseeds_file_name);
		mc.write_seeds_to_file(bSeeds, bseeds_file_name);
		cout << "_________________________________________________________________" << endl;
		cout << "[info] Estimate true spread......" << endl;
		computeTrueAdoption(kA, aseeds_file_name, dataset, qq, bseeds_file_name, ignore_B);

	}


	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//cout << "[info] *** running time for seed mining = " << elapsed_secs << " sec ***" << endl << endl;

	delete opt;
	cout << endl;
	return 0;
}



