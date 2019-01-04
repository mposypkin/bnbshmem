/*
 * A  multithreaded interval-based bnb solver Global-PAMIGO motivated by 
 * Branch-and-Bound interval global optimization on shared memory multiprocessors by  L. G. Casado J. A. Martínez I. García E. M. T. Hendrix 
 */

/* 
 * File:   tutorialbnb.cpp
 * Author: mposypkin
 *
 * Created on December 13, 2017, 3:22 PM
 */

#include <iostream>
#include <limits>
#include <random>
#include <algorithm>
#include <vector>
#include <iterator>
#include <functional>
#include <memory>
#include <thread>
#include <mutex>
#include <chrono>
#include <condition_variable>
#include <common/parbench.hpp>
#include <common/bnbiutils.hpp>
#include <common/bnbstat.hpp>

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;

static int gProcs = 8;

static double gEps;

static std::string gKnrec;

constexpr char knownRecord[] = "knrec";

std::atomic<double> gRecv;

std::vector<double> gRecord;

std::atomic<int> gNumWaitThreads(0);

std::vector<BnBStat> stat;

std::mutex gRecordLock;

std::mutex gPoolLock;

std::condition_variable gPoolCV;

std::atomic<long long int> gSteps;

const std::memory_order morder = std::memory_order_seq_cst;
//const std::memory_order morder = std::memory_order_relaxed;

#define EXCHNAGE_OPER compare_exchange_strong
//#define EXCHNAGE_OPER compare_exchange_weak

struct State {
    std::vector<Box> mPool;
};

void updateRecord(double rv, const std::vector<double> & record) {
    std::lock_guard<std::mutex> lock(gRecordLock);
    if (rv <= gRecv) {
        gRecv = rv;
        gRecord = record;
    }
}

void solve(State& s, const BM& bm) {
    const int dim = bm.getDim();
    std::vector<double> c(dim);
    double locRecv = std::numeric_limits<double>::max();
    std::vector<double> locRecord;
    int steps = 0;
    while (true) {
        steps++;
        Box b;
        {
            std::unique_lock<std::mutex> lock(gPoolLock);
            gNumWaitThreads++;
            while (s.mPool.empty() && (gNumWaitThreads < gProcs))
                gPoolCV.wait(lock);
            if (s.mPool.empty()) {
                gPoolCV.notify_one();
                break;
            }
            b = s.mPool.back();
            s.mPool.pop_back();
            gNumWaitThreads--;
        }
        mid(b, c);
        double v = bm.calcFunc(c);
        if (v < locRecv) {
            locRecv = v;
            locRecord = c;
        }
        double rv = gRecv.load(morder);
        while (v < rv) {
            gRecv.EXCHNAGE_OPER(rv, v, morder);
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= gRecv.load(morder) - gEps) {
            std::unique_lock<std::mutex> lock(gPoolLock);
            split(b, s.mPool);
            gPoolCV.notify_all();
        }
    }
    updateRecord(locRecv, locRecord);
    gSteps += steps;
}

double findMin(const BM& bm) {
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    State s;
    s.mPool.push_back(ibox);
    if (gKnrec == std::string(knownRecord)) {
        gRecv = bm.getGlobMinY();
    } else {
        gRecv = std::numeric_limits<double>::max();
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    gNumWaitThreads = 0;
    std::vector<std::thread> threads;
    for(int i = 0; i < gProcs; i ++) {
        threads.emplace_back(solve, std::ref(s), std::cref(bm));
    }
    for(auto & t : threads) {
        t.join();
    }

    end = std::chrono::system_clock::now();
    int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end - start)).count();
    std::cout << "Time: " << mseconds << " microsecond\n";
    std::cout << "Time per subproblem: " << (double) mseconds / (double) gSteps << " miscroseconds." << std::endl;
    std::cout << "Converged in " << gSteps << " steps\n";

    std::cout << "BnB found = " << gRecv << std::endl;
    std::cout << " at x [ ";
    std::copy(gRecord.begin(), gRecord.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n";
    stat.emplace_back((double) mseconds, gSteps);
    return gRecv;
}

bool testBench(const BM& bm) {
    bool rv = true;
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm;
    double v = findMin(bm);
    double diff = v - bm.getGlobMinY();
    if (diff > gEps) {
        std::cout << "BnB failed for " << bm.getDesc() << " benchmark " << std::endl;
        rv = false;
    }
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;
    char c;
    return rv;
}

int main(int argc, char* argv[]) {
    std::string bench;
#if 0    
    Benchmarks<double> tests;
#else 
    ParBenchmarks<double> tests;
#endif    

    int nruns = 0;

    if ((argc == 2) && (std::string(argv[1]) == std::string("list"))) {
        for (auto b : tests) {
            std::cout << b->getDesc() << "\n";
        }
        return 0;
    } else if (argc == 6) {
        nruns = atoi(argv[1]);
        bench = argv[2];
        gKnrec = argv[3];
        gEps = atof(argv[4]);
        gProcs = atoi(argv[5]);
    } else {
        std::cerr << "Usage: " << argv[0] << " number_of_runs name_of_bench knrec|unknrec eps virtual_procs_number\n";
        std::cerr << "or to list benchmarks run:\n";
        std::cerr << argv[0] << " list\n";
        return -1;
    }
    std::cout << "Local PAMIGO BnB solver with np = " << gProcs << "\n";
    std::cout << "record is " << (gRecv.is_lock_free() ? "lock free" : "not lock free") << std::endl;
#if 0    
    PowellSingular2Benchmark<double> pb(8);
    testBench(pb);
#else        
    for (int z = 0; z < nruns; z++)
        for (auto bm : tests) {
            if (bench == bm->getDesc())
                testBench(*bm);
        }
#endif   
    std::cout << "Statistics:\n" << stat;
}
