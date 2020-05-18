/*
 * A simple interval-based bnb solver
 */

/* 
 * File:   bnb_test_prelim.cpp fork from tutorialbnb.cpp
 * Author: Gorchakov A.Y.
 *
 * Created on Apr 18, 2018, 4:12 PM
 */

#include <iostream>
#include <limits>
#include <random>
#include <algorithm>
#include <vector>
#include <iterator>
#include <chrono>
#include <omp.h>
#include <sys/time.h>
#include <algorithm>

#include <common/parbench.hpp>
#include <common/bnbiutils.hpp>
#include <common/bnbstat.hpp>


static int gMaxStepsTotal;

static int gProcs;

static std::string gKnrec;

constexpr char gKnownRecord[] = "knrec";

static double gEps;

std::vector<BnBStat> gStat;

std::ostream* gOutStream = &std::cout;

std::ostream* gStatStream = &std::cout;

double findMin(const BM& bm, const double eps, const long long int maxstep, const int n_thr) {
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    std::vector<std::vector < Box >> pool(n_thr);
    std::vector<std::vector < Box >> pool_new(n_thr);
    pool[0].push_back(ibox);
    std::vector<double> c(dim);
    std::vector<double> recordVec(dim);
    double recordVal;
    if (gKnrec == std::string(gKnownRecord)) {
        recordVal = bm.getGlobMinY();
    } else {
        recordVal = std::numeric_limits<double>::max();
    }

    long long int steps = 0;
    long long pool_size = 0;
    const auto start = std::chrono::system_clock::now();
    while (true) {
#pragma omp parallel shared(pool, pool_new, recordVec, recordVal, steps, pool_size) firstprivate(c) num_threads(n_thr)
        {
#pragma omp for 
            for (int i = 0; i < n_thr; i++) pool_new[i].clear();

            for (int i = 0; i < n_thr; i++)
#pragma omp for nowait schedule(dynamic)
                for (long long int j = 0; j < pool[i].size(); j++) {
                    const int thread_num = omp_get_thread_num();
                    Box b = pool[i][j];
                    mid(b, c);
                    double v = bm.calcFunc(c);
                    double rv;
                    //                    double rv = gRecv.load(morder);
                    auto lb = bm.calcInterval(b).lb();
#pragma omp atomic read
                    rv = recordVal;
                    if (v < rv) {
#pragma omp critical
                        {
                            if (v < recordVal) {
#pragma omp atomic write
                                recordVal = v;
                                rv = v;
                                recordVec = c;
                            }
                        }
                    }
                    if (lb <= rv - eps && pool_size + steps < maxstep) {
                        split(b, pool_new[thread_num]);
                    }
                }
        }
        for (int i = 0; i < n_thr; i++) steps += pool[i].size();
        if (steps >= maxstep) break;
        pool_size = 0;
        for (int i = 0; i < n_thr; i++) pool_size += pool_new[i].size();
        if (pool_size == 0) break;
        if (pool_size + steps >= maxstep) {
            pool_size = 0;
            for (int i = 0; i < n_thr; i++) {
                if (pool_size + pool_new[i].size() + steps > maxstep) {
                    pool_new[i].erase(pool_new[i].begin() + maxstep - steps - pool_size, pool_new[i].end());
                }
                pool_size += pool_new[i].size();
            }
        }
        pool.swap(pool_new);
    }
    auto end = std::chrono::system_clock::now();
    unsigned long int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end - start)).count();
    *gOutStream << "Time: " << mseconds << " microsecond\n";
    *gOutStream << "Time per subproblem: " << (double) mseconds / (double) steps << " miscroseconds." << std::endl;
    if (steps >= gMaxStepsTotal) {
        *gOutStream << "Failed to converge in " << gMaxStepsTotal << " steps\n";
    } else {
        *gOutStream << "Converged in " << steps << " steps\n";
    }
    *gOutStream << "BnB found = " << recordVal << std::endl;
    *gOutStream << " at x [ ";
    std::copy(recordVec.begin(), recordVec.end(), std::ostream_iterator<double>(*gOutStream, " "));
    *gOutStream << "]\n";

    *gOutStream << bm.getDesc() << ":" << steps << "\n";

    gStat.emplace_back(mseconds, steps);
    return recordVal;
}

bool testBench(const BM& bm) {
    bool rv = true;
    *gOutStream << "*************Testing benchmark**********" << std::endl;
    *gOutStream << bm;
    double v = findMin(bm, gEps, gMaxStepsTotal, gProcs);
    double diff = v - bm.getGlobMinY();
    if (diff > gEps) {
        *gOutStream << "BnB failed for " << bm.getDesc() << " benchmark " << std::endl;
        rv = false;
    }
    *gOutStream << "the difference is " << v - bm.getGlobMinY() << std::endl;
    *gOutStream << "****************************************" << std::endl << std::endl;
    return rv;
}

int main(int argc, char* argv[]) {
    std::string bench;
    int nruns;
#if 0   
    Benchmarks<double> tests;
#else 
    ParBenchmarks<double> tests;
#endif    


    if ((argc == 2) && (std::string(argv[1]) == std::string("list"))) {
        for (auto b : tests) {
            *gOutStream << b->getDesc() << "\n";
        }
        return 0;
    } else if (argc == 7 || argc == 8) {
        nruns = atoi(argv[1]);
        bench = argv[2];
        gKnrec = argv[3];
        gEps = atof(argv[4]);
        gMaxStepsTotal = atoi(argv[5]);
        gProcs = atoi(argv[6]);
        if (argc == 8) {
            gOutStream = new BnbStream(nullptr);
        }
    } else {
        std::cerr << "Usage: " << argv[0] << " num_of_runs name_of_bench knrec|unknrec eps max_steps omp_thread_number [statonly]\n";
        std::cerr << "or to list benchmarks run:\n";
        std::cerr << argv[0] << " list\n";
        return -1;
    }

    *gOutStream << "Openmp PBnB solver with " << std::endl;
    for (int z = 0; z < nruns; z++)
        for (auto bm : tests) {
            if (bench == bm->getDesc())
                testBench(*bm);
        }


    *gStatStream << "Statistics for " << bench << ":\n" << gStat;

}
