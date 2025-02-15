#include "timgraph.h"
using namespace chrono;
void write_seeds_to_file(vector<int> seeds, string file_name) {
	ofstream myfile;
	myfile.open(file_name.c_str());
	if (myfile.is_open()) {
		for (int i = 0; i < seeds.size(); ++i) {
			myfile << seeds.at(i) << endl;

		}
		myfile.close();
	}
	else {
		cerr << "[error] unable to open output file; no output was written " << file_name << endl;
		exit(1);
	}
}
void TimGraph::EstimateTheta(double eps, double ell = 0.1) {
    reset_hyper();
    double epsprime = eps * sqrt(2.0);
    double LB = 1.0, spread;
    double maxRounds = max(max(log2((double)n), 1.0) - 1.0, 1.0);
    int64 theta0 = 0, theta1;
    ell = ell + log(2) / log(n);
    for (int r = 1; r < maxRounds; r++) {
        double x = max(double(n) / pow(2, r), 1.0);
        theta1 = (int64)(LambdaPrime(epsprime, ell, n) / x);

        if (theta0 < theta1) {

            if (theta1 > hyperGT_A.size()) {
                hyperGT_A.resize(theta1);
                hyperGT_B.resize(theta1);
            }
////////////////////////////////////////////Parallel////////////////////////////////////////////
omp_set_num_threads(4);
#pragma omp parallel
            {
                    graph->init_EdgeStatus();
                    init_visit();

#pragma omp for
            for (int i = theta0; i < theta1; i++) {
                 int r = randomInt(n);
                 //Afirst
                 //Bfirst
                 BuildSingleRRSet_q1(r, i, true, true, hyperGT_A, hyperG_A);
                 //BuildSingleRRSet_q1(r, i, true, false, hyperGT_A, hyperG_A);
                 time_stamp++;
                 BuildSingleRRSet_q1(r, i, false, false, hyperGT_B, hyperG_B);
                 //BuildSingleRRSet_q1(r, i, false, true, hyperGT_B, hyperG_B);
                 time_stamp++;

                 graph->reset_thEdgeStatus();

 
            }
            cout << "Thread ID: " << omp_get_thread_num() << " finished!" << endl;
///////////////////////////////////////////////////////////////////////////////////////////

            }
            BuildSeedSetIG_q1(false);
            spread = computeInfHyperGraph_Order(theta1);
            cout << "Estimated spread = " << spread << "\t round = " << theta1 << endl;
            theta0 = theta1;
        }
        if (spread >= (1.0 + epsprime) * x) {
            LB = spread / (1.0 + epsprime);
            break;
        }

    }
    theta1 = (int64)LambdaStar(eps, ell) / LB;
    cout << "[info] Estimated theta = " << theta1 << endl;
    this->theta = theta1;
}


double TimGraph::LambdaPrime(double epsprime, double ell, int n) {
    static double SMALL_DOUBLE = 1e-16;
    double cst = (2.0 + 2.0 / 3.0 * epsprime) / max(epsprime * epsprime, SMALL_DOUBLE);
    double part2 = logcnk(n, kA+kB) + logcnk(n, kB);
    part2 += ell * log(max((double)n, 1.0)) + log(2);
    part2 += log(max(log2(max((double)n, 1.0)), 1.0));
    // calc \lambda'
    double lambda = cst * part2 * n;
    return lambda;
}
double TimGraph::LambdaStar(double eps, double ell)
{
    static double SMALL_DOUBLE = 1e-16;
    static double APPRO_RATIO = (1.0 - 1.0 / 2.0);
    double logsum = ell * log(n) + log(4);
    // calc \alpha and \beta
    double alpha = sqrt(max(logsum, SMALL_DOUBLE));
    double beta = sqrt(APPRO_RATIO * (logcnk(n, kA) + logcnk(n, kB) + logsum));
    // calc \lambda*
    double lambda = 2.0 * n / max(pow(eps, 2), SMALL_DOUBLE);
    lambda = lambda * pow(APPRO_RATIO * alpha + beta, 2);
    return lambda;
}

void TimGraph::minseeds_Order_q1(double epsilon, string file_SA, string file_SB) {
    auto total_start = high_resolution_clock::now();
    double total_build_time = 0.0;
    double total_select_time = 0.0;
    int runs = 5;
    cout << "[IMM-q1] Starting " << runs << " runs..." << endl;
    for (int i = 0; i < runs; i++) {
        cout << "--------------------------------------------------------" << endl;
        cout << "Run " << i + 1 << "/" << runs << " started." << endl;
        auto run_start = high_resolution_clock::now();
        auto build_start = high_resolution_clock::now();

        EstimateTheta(epsilon);
        BuildHyperGraphR_q1();
        auto build_end = high_resolution_clock::now();
        double build_duration = duration_cast<milliseconds>(build_end - build_start).count() / 1000.0;
        total_build_time += build_duration;
        cout << "Step 1 [BuildHypergraphR]: Completed in " << build_duration << " sec." << endl;

        auto select_start = high_resolution_clock::now();
        BuildSeedSetIG_q1(true);
        auto select_end = high_resolution_clock::now();
        double select_duration = duration_cast<milliseconds>(select_end - select_start).count() / 1000.0;
        total_select_time += select_duration;
        cout << "Step 2 [BuildSeedSetIG_Order]: Completed in " << select_duration << " sec." << endl;

        auto run_end = high_resolution_clock::now();
        double run_duration = duration_cast<milliseconds>(run_end - run_start).count() / 1000.0;
        cout << "Run " << i + 1 << "/" << runs << " completed in " << run_duration << " sec." << endl;
    }
    // total run time
    auto total_end = high_resolution_clock::now();
    double total_duration = duration_cast<milliseconds>(total_end - total_start).count() / 1000.0;
    cout << "========================================================" << endl;
    cout << "[IMM-Order Afirst] All runs completed." << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Average time per run (total): " << total_duration / runs << " sec." << endl;
    cout << "Average BuildHypergraphR_Order time: " << total_build_time / runs << " sec." << endl;
    cout << "Average BuildSeedSetIG_Order time: " << total_select_time / runs << " sec." << endl;
    cout << "Total time for all runs: " << total_duration << " sec." << endl;
    write_seeds_to_file(aSeeds, file_SA);
    write_seeds_to_file(bSeeds, file_SB);
    // seed repeat rate
    double repeat_rate = calculateRepeatRate(aSeeds, bSeeds) * 100;
    cout << "Seed set repeat rate: " << repeat_rate << '%' << endl;
}


void TimGraph::minseeds_idp(double epsilon, string file_SA, string file_SB) {
    auto total_start = high_resolution_clock::now();
    double total_build_time = 0.0;
    double total_select_time = 0.0;
    int runs = 5;
    cout << "[IMM-NoCompete] Starting " << runs << " runs..." << endl;
    for (int i = 0; i < runs; i++) {
        cout << "--------------------------------------------------------" << endl;
        cout << "Run " << i + 1 << "/" << runs << " started." << endl;
        auto run_start = high_resolution_clock::now();
        auto build_start = high_resolution_clock::now();
        EstimateTheta(epsilon);
        BuildHyperGraphR_idp();
        auto build_end = high_resolution_clock::now();
        double build_duration = duration_cast<milliseconds>(build_end - build_start).count() / 1000.0;
        total_build_time += build_duration;
        cout << "Step 1 [BuildHypergraphR]: Completed in " << build_duration << " sec." << endl;
        
        auto select_start = high_resolution_clock::now();
        BuildSeedSetIG_idp();
        auto select_end = high_resolution_clock::now();
        double select_duration = duration_cast<milliseconds>(select_end - select_start).count() / 1000.0;
        total_select_time += select_duration;
        cout << "Step 2 [BuildSeedSetIG_Order]: Completed in " << select_duration << " sec." << endl;

        auto run_end = high_resolution_clock::now();
        double run_duration = duration_cast<milliseconds>(run_end - run_start).count() / 1000.0;
        cout << "Run " << i + 1 << "/" << runs << " completed in " << run_duration << " sec." << endl;
    }
  
    auto total_end = high_resolution_clock::now();
    double total_duration = duration_cast<milliseconds>(total_end - total_start).count() / 1000.0;
    cout << "========================================================" << endl;
    cout << "[IMM-Order Afirst] All runs completed." << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Average time per run (total): " << total_duration / runs << " sec." << endl;
    cout << "Average BuildHypergraphR_Order time: " << total_build_time / runs << " sec." << endl;
    cout << "Average BuildSeedSetIG_Order time: " << total_select_time / runs << " sec." << endl;
    cout << "Total time for all runs: " << total_duration << " sec." << endl;
    write_seeds_to_file(aSeeds, file_SA);
    write_seeds_to_file(bSeeds, file_SB);
   
    double repeat_rate = calculateRepeatRate(aSeeds, bSeeds) * 100;
    cout << "Seed set repeat rate: " << repeat_rate << '%' << endl;
}
vector<vector<vector<int>>> TimGraph::mineSeeds_Order_Afirst(int runs, double epsilon, string file_SA, string file_SB) {
    // Record the start time of the overall execution
    auto total_start = high_resolution_clock::now();

    // Store seed sets for multiple runs
    vector<vector<vector<int>>> results(2, vector<vector<int>>(runs));

    // Variables to accumulate time statistics
    double total_build_time = 0.0;
    double total_select_time = 0.0;
    double total_estimate_time = 0.0; // New: To accumulate time for EstimateTheta

    cout << "[IMM-Order Afirst] Starting " << runs << " runs..." << endl;

    // Start multiple runs
    for (int i = 0; i < runs; i++) {
        cout << "--------------------------------------------------------" << endl;
        cout << "Run " << i + 1 << "/" << runs << " started." << endl;

        // Step 0: Time statistics for EstimateTheta
        auto estimate_start = high_resolution_clock::now();  // Record start time for EstimateTheta
        //EstimateTheta(epsilon);
        auto estimate_end = high_resolution_clock::now();    // Record end time for EstimateTheta
        double estimate_duration = duration_cast<milliseconds>(estimate_end - estimate_start).count() / 1000.0;
        total_estimate_time += estimate_duration; // Accumulate time

        cout << "Step 0 [EstimateTheta]: Completed in " << estimate_duration << " sec." << endl;

        // Record the start time of the current run
        auto run_start = high_resolution_clock::now();

        // Step 1: Build Hypergraph
        auto build_start = high_resolution_clock::now();
        BuildHypergraphR_Order(i, true); // Construct RR sets
        auto build_end = high_resolution_clock::now();
        double build_duration = duration_cast<milliseconds>(build_end - build_start).count() / 1000.0;
        total_build_time += build_duration;

        cout << "Step 1 [BuildHypergraphR_Order]: Completed in " << build_duration << " sec." << endl;

        // Step 2: Build Seed Set
        auto select_start = high_resolution_clock::now();
        BuildSeedSetIG_Order_Afirst(); // Construct seed sets
        auto select_end = high_resolution_clock::now();
        double select_duration = duration_cast<milliseconds>(select_end - select_start).count() / 1000.0;
        total_select_time += select_duration;

        cout << "Step 2 [BuildSeedSetIG_Order]: Completed in " << select_duration << " sec." << endl;

        // Save seed sets to the results for multiple runs
        results[0][i] = aSeeds;
        results[1][i] = bSeeds;
        // Record the end time of the current run
        auto run_end = high_resolution_clock::now();
        double run_duration = duration_cast<milliseconds>(run_end - run_start).count() / 1000.0;

        cout << "Run " << i + 1 << "/" << runs << " completed in " << run_duration << " sec." << endl;
    }

    // Record the end time of the overall execution
    auto total_end = high_resolution_clock::now();
    double total_duration = duration_cast<milliseconds>(total_end - total_start).count() / 1000.0;

    // Output statistics
    cout << "========================================================" << endl;
    cout << "[IMM-Order Afirst] All runs completed." << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Average time per run (total): " << total_duration / runs << " sec." << endl;
    cout << "Average BuildHypergraphR_Order time: " << total_build_time / runs << " sec." << endl;
    cout << "Average BuildSeedSetIG_Order time: " << total_select_time / runs << " sec." << endl;
    cout << "Average EstimateTheta time: " << total_estimate_time / runs << " sec." << endl;  // Output average time for EstimateTheta
    cout << "Total time for all runs: " << total_duration << " sec." << endl;

    // Calculate the repeat rate of seed sets
    double repeat_rate = calculateRepeatRate(aSeeds, bSeeds) * 100;
    cout << "Seed set repeat rate: " << repeat_rate << '%' << endl;

    // Write the final seed sets to files
    write_seeds_to_file(aSeeds, file_SA);
    write_seeds_to_file(bSeeds, file_SB);
    cout << "Results saved to files: " << file_SA << " and " << file_SB << endl;
    cout << "========================================================" << endl;
    return results;
}

vector<vector<vector<int>>> TimGraph::mineSeeds_Order_Bfirst(int runs, double epsilon, string file_SA, string file_SB)
{

    auto total_start = high_resolution_clock::now();
    vector<vector<vector<int>>> results(2, vector<vector<int>>(runs));

    double total_build_time = 0.0;
    double total_select_time = 0.0;

    cout << "[IMM-Order Bfirst] Starting " << runs << " runs..." << endl;

    for (int i = 0; i < runs; i++) {
        cout << "--------------------------------------------------------" << endl;
        cout << "Run " << i + 1 << "/" << runs << " started." << endl;

        auto run_start = high_resolution_clock::now();

        auto build_start = high_resolution_clock::now();
        BuildHypergraphR_Order(i, false); 
        auto build_end = high_resolution_clock::now();
        double build_duration = duration_cast<milliseconds>(build_end - build_start).count() / 1000.0;
        total_build_time += build_duration;

        cout << "Step 1 [BuildHypergraphR_Order]: Completed in " << build_duration << " sec." << endl;

        auto select_start = high_resolution_clock::now();
        BuildSeedSetIG_Order_Bfirst(); 
        auto select_end = high_resolution_clock::now();
        double select_duration = duration_cast<milliseconds>(select_end - select_start).count() / 1000.0;
        total_select_time += select_duration;

        cout << "Step 2 [BuildSeedSetIG_Order]: Completed in " << select_duration << " sec." << endl;

        results[0][i] = aSeeds;
        results[1][i] = bSeeds;

        auto run_end = high_resolution_clock::now();
        double run_duration = duration_cast<milliseconds>(run_end - run_start).count() / 1000.0;

        cout << "Run " << i + 1 << "/" << runs << " completed in " << run_duration << " sec." << endl;
    }

    auto total_end = high_resolution_clock::now();
    double total_duration = duration_cast<milliseconds>(total_end - total_start).count() / 1000.0;

    cout << "========================================================" << endl;
    cout << "[IMM-Order Bfirst] All runs completed." << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Average time per run (total): " << total_duration / runs << " sec." << endl;
    cout << "Average BuildHypergraphR_Order time: " << total_build_time / runs << " sec." << endl;
    cout << "Average BuildSeedSetIG_Order time: " << total_select_time / runs << " sec." << endl;
    cout << "Total time for all runs: " << total_duration << " sec." << endl;

    double repeat_rate = calculateRepeatRate(aSeeds, bSeeds) * 100;
    cout << "Seed set repeat rate: " << repeat_rate << '%' << endl;

    write_seeds_to_file(aSeeds, file_SA);
    write_seeds_to_file(bSeeds, file_SB);
    cout << "Results saved to files: " << file_SA << " and " << file_SB << endl;
    cout << "========================================================" << endl;
    return results;
}

