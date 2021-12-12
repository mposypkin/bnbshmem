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

static int gMaxStepsTotal = 0;

static double gEps;

static std::string gKnrec;

constexpr int DFS_POLICY = 1;

constexpr int WFS_POLICY = 2;

std::atomic<int> gBranchPolicy(DFS_POLICY);

constexpr char gKnownRecord[] = "knrec";

std::atomic<double> gRecv;

std::vector<double> gRecord;

std::atomic<int> gNumWaitThreads(0);

std::vector<BnBStat> gStat;

std::mutex gRecordLock;

std::mutex gPoolLock;

std::condition_variable gPoolCV;

std::atomic<long long int> gSteps;

std::ostream* gOutStream = &std::cout;

std::ostream* gStatStream = &std::cout;

const std::memory_order gMorder = std::memory_order_seq_cst;
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
    while (true) {
        Box b;
        {
            std::unique_lock<std::mutex> lock(gPoolLock);
            if (gSteps >= gMaxStepsTotal)
                break;
            gNumWaitThreads++;
            while (s.mPool.empty() && (gNumWaitThreads < gProcs))
                gPoolCV.wait(lock);
            if (s.mPool.empty()) {
                gPoolCV.notify_one();
                break;
            }
            if(gBranchPolicy == DFS_POLICY) {
                b = s.mPool.back();
                s.mPool.pop_back();
            } else {
                b = s.mPool.front();
                s.mPool.erase(s.mPool.begin());
            }
            gSteps++;
            gNumWaitThreads--;
        }
        mid(b, c);
        double v = bm.calcFunc(c);
        if (v < locRecv) {
            locRecv = v;
            locRecord = c;
        }
        double rv = gRecv.load(gMorder);
        while (v < rv) {
            gRecv.EXCHNAGE_OPER(rv, v, gMorder);
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= gRecv.load(gMorder) - gEps) {
            std::unique_lock<std::mutex> lock(gPoolLock);
            split(b, s.mPool);
            gPoolCV.notify_all();
        }
    }
    updateRecord(locRecv, locRecord);
}

double findMin(const BM& bm) {
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    State s;
    s.mPool.push_back(ibox);
    if (gKnrec == std::string(gKnownRecord)) {
        gRecv = bm.getGlobMinY();
    } else {
        gRecv = std::numeric_limits<double>::max();
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    gNumWaitThreads = 0;
    gSteps = 0;
    std::vector<std::thread> threads;
    start = std::chrono::system_clock::now();
    for (int i = 0; i < gProcs; i++) {
        threads.emplace_back(solve, std::ref(s), std::cref(bm));
    }
    for (auto & t : threads) {
        t.join();
    }

    end = std::chrono::system_clock::now();
    int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end - start)).count();
    *gOutStream << "Time: " << mseconds << " microsecond\n";
    *gOutStream << "Time per subproblem: " << (double) mseconds / (double) gSteps << " miscroseconds." << std::endl;
    if (gSteps >= gMaxStepsTotal) {
        *gOutStream << "Failed to converge in " << gMaxStepsTotal << " steps\n";
    } else {
        *gOutStream << "Converged in " << gSteps << " steps\n";
    }

    *gOutStream << "BnB found = " << gRecv << std::endl;
    *gOutStream << " at x [ ";
    std::copy(gRecord.begin(), gRecord.end(), std::ostream_iterator<double>(*gOutStream, " "));
    *gOutStream << "]\n";
    gStat.emplace_back(mseconds, gSteps);
    return gRecv;
}

bool testBench(const BM& bm) {
    bool rv = true;
    *gOutStream << "*************Testing benchmark**********" << std::endl;
    *gOutStream << "Branching policy = " << ((gBranchPolicy == DFS_POLICY) ? "dfs" : "wfs") << "\n";
    *gOutStream << bm;
    double v = findMin(bm);
    double diff = v - bm.getGlobMinY();
    if (diff > gEps) {
        *gOutStream << "BnB failed for " << bm.getDesc() << " benchmark " << std::endl;
        rv = false;
    }
    *gOutStream << "the difference is " << v - bm.getGlobMinY() << std::endl;
    *gOutStream << "****************************************" << std::endl << std::endl;
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

    gBranchPolicy = DFS_POLICY;

    if ((argc == 2) && (std::string(argv[1]) == std::string("list"))) {
        for (auto b : tests) {
            *gOutStream << b->getDesc() << "\n";
        }
        return 0;
    } else if (argc >= 7 && argc <=  9) {
        nruns = atoi(argv[1]);
        bench = argv[2];
        gKnrec = argv[3];
        gEps = atof(argv[4]);
        gMaxStepsTotal = atoi(argv[5]);
        gProcs = atoi(argv[6]);
        if (argc >= 8) {
          if (std::string(argv[7]) == std::string("dfs"))
              gBranchPolicy = DFS_POLICY;
          else if (std::string(argv[7]) == std::string("wfs"))
              gBranchPolicy = WFS_POLICY;
          else if (std::string(argv[7]) == std::string("statonly"))
              gOutStream = new BnbStream(nullptr);
        }
        if (argc == 9) {
            if (std::string(argv[8]) == std::string("statonly"))
                gOutStream = new BnbStream(nullptr);
        }
    } else {
        std::cerr << "Usage: " << argv[0] << " number_of_runs name_of_bench knrec|unknrec eps max_steps virtual_procs_number [dfs| wfs statonly]\n";
        std::cerr << "or to list benchmarks run:\n";
        std::cerr << argv[0] << " list\n";
        return -1;
    }
    *gOutStream << "Global PAMIGO BnB solver with np = " << gProcs << "\n";
    *gOutStream << "record is " << (gRecv.is_lock_free() ? "lock free" : "not lock free") << std::endl;
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
    *gStatStream << "Statistics:\n" << gStat;
}
