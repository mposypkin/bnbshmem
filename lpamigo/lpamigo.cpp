/*
 * A  multithreaded interval-based bnb solver L-PAMIGO motivated by 
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

constexpr char gKnownRecord[] = "knrec";

std::atomic<double> gRecv;

std::vector<double> gRecord;

std::atomic<int> gNumWaitThreads(0);

std::vector<BnBStat> gStat;

std::mutex gRecordLock;

std::mutex gPoolLock;

std::condition_variable gPoolCV;

std::atomic<long long int> gSteps;

const std::memory_order gMorder = std::memory_order_seq_cst;
//const std::memory_order morder = std::memory_order_relaxed;

std::ostream* gOutStream = &std::cout;

std::ostream* gStatStream = &std::cout;

#define EXCHNAGE_OPER compare_exchange_strong
//#define EXCHNAGE_OPER compare_exchange_weak

struct State {

    void merge(const State& s) {
        mPool.insert(mPool.end(), s.mPool.begin(), s.mPool.end());
    }

    void split(State& s) {
        int k = 0;
        for (auto i = mPool.begin(); i != mPool.end();) {
            if (k % 2) {
                s.mPool.push_back(*i);
                i = mPool.erase(i);
            } else
                i++;
            k++;
        }

    }

    std::vector<Box> mPool;
};

void updateRecord(double rv, const std::vector<double> & record) {
    std::lock_guard<std::mutex> lock(gRecordLock);
    if (rv <= gRecv) {
        gRecv = rv;
        gRecord = record;
    }
}

void solve(std::shared_ptr<State> s, const BM& bm) {
    const int dim = bm.getDim();
    std::vector<double> c(dim);
    double locRecv = std::numeric_limits<double>::max();
    std::vector<double> locRecord;
    while (!s->mPool.empty()) {
        if(gSteps >= gMaxStepsTotal)
            break;
        gSteps ++;
        if (s->mPool.size() >= 2) {
            if (gNumWaitThreads < gProcs) {
                std::lock_guard<std::mutex> lock(gPoolLock);
                auto sn = std::make_shared<State>();
                s->split(*sn);
                std::thread th(solve, sn, std::cref<BM>(bm));
                th.detach();
                gNumWaitThreads++;
            }
        }
        Box b = s->mPool.back();
        s->mPool.pop_back();
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
            split(b, s->mPool);
        }
    }
    updateRecord(locRecv, locRecord);
    std::unique_lock<std::mutex> lock(gPoolLock);
    gNumWaitThreads--;
    gPoolCV.notify_one();
}

double findMin(const BM& bm) {
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    auto s = std::make_shared<State>();
    s->mPool.push_back(ibox);
    if (gKnrec == std::string(gKnownRecord)) {
        gRecv = bm.getGlobMinY();
    } else {
        gRecv = std::numeric_limits<double>::max();
    }
    gSteps = 0;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
#if 0   
    solveSerial(s, bm, eps);
#else    
    gNumWaitThreads = 1;
    solve(s, bm);
    {
        std::unique_lock<std::mutex> lock(gPoolLock);
        gPoolCV.wait(lock, []() {
            return gNumWaitThreads == 0;
        });
        *gOutStream << "==" << gNumWaitThreads << "==\n";
    }
#endif
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
    gStat.emplace_back((double) mseconds, gSteps);
    return gRecv;
}

bool testBench(const BM& bm) {
    bool rv = true;
    *gOutStream << "*************Testing benchmark**********" << std::endl;
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
        std::cerr << "Usage: " << argv[0] << " number_of_runs name_of_bench knrec|unknrec eps max_steps virtual_procs_number [statonly]\n";
        std::cerr << "or to list benchmarks run:\n";
        std::cerr << argv[0] << " list\n";
        return -1;
    }
    *gOutStream << "Local PAMIGO BnB solver with np = " << gProcs << "\n";
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
    *gStatStream << "Statistics for " << bench << ":\n" << gStat;
}
