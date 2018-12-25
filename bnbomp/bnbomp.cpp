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
#include <common/parbench.hpp>
#include <algorithm>

static int gMaxStepsTotal;

static int gProcs;

static std::string gKnrec;

constexpr char knownRecord[] = "knrec";

static double gEps;

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;

std::atomic<double> gRecv;



//const std::memory_order morder = std::memory_order_seq_cst;
const std::memory_order morder = std::memory_order_relaxed;

#define EXCHNAGE_OPER compare_exchange_strong
//#define EXCHNAGE_OPER compare_exchange_weak


double len(const Interval<double>& I) {
    return I.rb() - I.lb();
}

void split(const Box& ibox, std::vector<Box>& v) {
    auto result = std::max_element(ibox.begin(), ibox.end(),
            [](const Interval<double>& f, const Interval<double>& s) {
                return len(f) < len(s);
            });
    const int i = result - ibox.begin();
    const double maxlen = len(ibox[i]);
    Box b1(ibox);
    Interval<double> ilow(ibox[i].lb(), ibox[i].lb() + 0.5 * maxlen);
    b1[i] = ilow;
    Box b2(ibox);
    Interval<double> iupper(ibox[i].lb() + 0.5 * maxlen, ibox[i].rb());
    b2[i] = iupper;
    v.push_back(std::move(b1));
    v.push_back(std::move(b2));
}

void getCenter(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}


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
    if (gKnrec == std::string(knownRecord)) {
        recordVal = bm.getGlobMinY();
    } else {
        recordVal = std::numeric_limits<double>::max();
    }

    long long int step = 0;
    long long pool_size = 0;
    const auto start = std::chrono::system_clock::now();
    while (true) {
#pragma omp parallel shared(pool, pool_new, recordVec, recordVal, step, pool_size) firstprivate(c) num_threads(n_thr)
        {
#pragma omp for 
            for (int i = 0; i < n_thr; i++) pool_new[i].clear();

            for (int i = 0; i < n_thr; i++)
#pragma omp for nowait schedule(dynamic)
                for (long long int j = 0; j < pool[i].size(); j++) {
                    const int thread_num = omp_get_thread_num();
                    Box b = pool[i][j];
                    getCenter(b, c);
                    double v = bm.calcFunc(c);
                    double rv;
//                    double rv = gRecv.load(morder);
                    auto lb = bm.calcInterval(b).lb();
#pragma omp atomic read
                    rv = recordVal;
                    if(v < rv) {
#pragma omp critical
{
                    if(v < recordVal) {
#pragma omp atomic write
                       recordVal = v;
                       rv = v;
                       recordVec = c;
                    }
}
}
                    if (lb <= rv - eps && pool_size + step < maxstep) {
                       split(b, pool_new[thread_num]);
                    }
                }
        }
        for (int i = 0; i < n_thr; i++) step += pool[i].size();
        if (step >= maxstep) break;
        pool_size = 0;
        for (int i = 0; i < n_thr; i++) pool_size += pool_new[i].size();
        if (pool_size == 0) break;
        if (pool_size + step >= maxstep) {
            pool_size = 0;
            for (int i = 0; i < n_thr; i++) {
                if (pool_size + pool_new[i].size() + step > maxstep) {
                    pool_new[i].erase(pool_new[i].begin() + maxstep - step - pool_size, pool_new[i].end());
                }
                pool_size += pool_new[i].size();
            }
        }
        pool.swap(pool_new);
    }
    auto end = std::chrono::system_clock::now();
    unsigned long int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end - start)).count();
    std::cout << "Time: " << mseconds << " microsecond\n";
    std::cout << "Time per subproblem: " << (double) mseconds / (double) step << " miscroseconds." << std::endl;
    if (step >= gMaxStepsTotal) {
        std::cout << "Failed to converge in " << gMaxStepsTotal << " steps\n";
    } else {
        std::cout << "Converged in " << step << " steps\n";
    }
    std::cout << "BnB found = " << recordVal << std::endl;
    std::cout << " at x [ ";
    std::copy(recordVec.begin(), recordVec.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n";

    std::cout << bm.getDesc() << ":" <<step <<"\n";

    return recordVal;
}

bool testBench(const BM& bm) {
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm;
    double v = findMin(bm, gEps, gMaxStepsTotal, gProcs);
    double diff = v - bm.getGlobMinY();
    if (diff > gEps) {
        std::cout << "BnB failed for " << bm.getDesc() << " benchmark " << std::endl;
    }
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;
}

main(int argc, char* argv[]) {
    std::string bench;
    Benchmarks<double> tests;

    if ((argc == 2) && (std::string(argv[1]) == std::string("list"))) {
        for (auto b : tests) {
            std::cout << b->getDesc() << "\n";
        }
        return 0;
    } else if (argc == 6) {
        bench = argv[1];
        gKnrec = argv[2];
        gEps = atof(argv[3]);
        gMaxStepsTotal = atoi(argv[4]);
        gProcs = atoi(argv[5]);
    } else {
        std::cerr << "Usage: " << argv[0] << " name_of_bench knrec|unknrec eps max_steps omp_thread_number\n";
        std::cerr << "or to list benchmarks run:\n";
        std::cerr << argv[0] << " list\n";
        return -1;
    }   
        
  std::cout << "Openmp PBnB solver with " << std::endl;
  for(int z = 0; z < 1; z++)
    for (auto bm : tests) {
        if (bench == bm->getDesc())
            testBench(*bm);
    }
}
