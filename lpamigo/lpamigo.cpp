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
    int steps = 0;
    while (!s->mPool.empty()) {
        steps++;
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
        double rv = gRecv.load(morder);
        while (v < rv) {
            gRecv.EXCHNAGE_OPER(rv, v, morder);
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= gRecv.load(morder) - gEps) {
            split(b, s->mPool);
        }
    }
    updateRecord(locRecv, locRecord);
    gSteps += steps;
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
    if (gKnrec == std::string(knownRecord)) {
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
        std::cout << "==" << gNumWaitThreads << "==\n";
    }
#endif
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
