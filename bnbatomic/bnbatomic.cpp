/*
 * A simple multithreaded interval-based bnb solver
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
#include <thread>
#include <mutex>
#include <chrono>
#include <common/parbench.hpp>
#include <common/bnbiutils.hpp>
#include <common/bnbstat.hpp>


static int gProcs = 8;

static int gMtStepsLimit = 1000;

static int gMaxStepsTotal = 1000000;

static double gEps;

static std::string gKnrec;

constexpr char knownRecord[] = "knrec";

std::atomic<double> gRecv;

std::vector<double> gRecord;

std::atomic_int gNumRecUpdates;

std::mutex gMutex;

std::vector<BnBStat> gStat;

std::ostream* gOutStream = &std::cout;

std::ostream* gStatStream = &std::cout;

//const std::memory_order morder = std::memory_order_seq_cst;
const std::memory_order morder = std::memory_order_relaxed;

#define EXCHNAGE_OPER compare_exchange_strong
//#define EXCHNAGE_OPER compare_exchange_weak

struct State {

    void merge(const State& s) {
        mSteps += s.mSteps;
        mPool.insert(mPool.end(), s.mPool.begin(), s.mPool.end());
    }

    void split(State& s1, State& s2) {
        const int remMaxSteps = mMaxSteps - mSteps;
        s1.mMaxSteps = remMaxSteps / 2;
        s2.mMaxSteps = remMaxSteps - s1.mMaxSteps;
        s1.mProcs = mProcs / 2;
        s2.mProcs = mProcs - s1.mProcs;

        while (true) {
            if (mPool.empty())
                break;
            s1.mPool.push_back(mPool.back());
            mPool.pop_back();
            if (mPool.empty())
                break;
            s2.mPool.push_back(mPool.back());
            mPool.pop_back();
        }
    }

    std::vector<Box> mPool;

    int mMaxSteps;

    int mSteps = 0;

    int mProcs;
};

std::ostream& operator<<(std::ostream & out, const State s) {
    out << "\"steps\" :" << s.mSteps << "\n";
    out << "\"max steps\" :" << s.mMaxSteps << "\n";
    return out;
}

void updateRecord(double newrv, const std::vector<double> &record) {
    double rv = gRecv;
    if (newrv < rv) {
        std::lock_guard<std::mutex> lock(gMutex);
        if (newrv < gRecv) {
            gRecv = newrv;
            gRecord = record;
            gNumRecUpdates++;
        }
    }
}

void solveSerial(State& s, const BM& bm) {
    const int dim = bm.getDim();
    std::vector<double> c(dim);
    while (!s.mPool.empty()) {
        s.mSteps++;
        Box b = s.mPool.back();
        s.mPool.pop_back();
        mid(b, c);
        double v = bm.calcFunc(c);
        updateRecord(v, c);
        auto lb = bm.calcInterval(b).lb();
        if (lb <= gRecv - gEps) {
            split(b, s.mPool);
        }
        if (s.mSteps >= s.mMaxSteps)
            break;
    }
}

void solve(State& s, const BM& bm) {
    if ((s.mProcs == 1) || (s.mMaxSteps <= gMtStepsLimit)) {
        solveSerial(s, bm);
    } else {
        auto presolve = [&](State & s) {
            const int dim = bm.getDim();
            std::vector<double> c(dim);
            while ((!s.mPool.empty()) && (s.mPool.size() < 2)) {
                s.mSteps++;
                Box b = s.mPool.back();
                s.mPool.pop_back();
                mid(b, c);
                double v = bm.calcFunc(c);
                updateRecord(v, c);
                auto lb = bm.calcInterval(b).lb();
                if (lb <= gRecv - gEps) {
                    split(b, s.mPool);
                }
                if (s.mSteps >= s.mMaxSteps)
                    break;
            }
        };

        while (true) {
            if (s.mPool.empty())
                break;
            //*outStream << s << "\n";
            presolve(s);
            const int remMaxSteps = s.mMaxSteps - s.mSteps;
            if (remMaxSteps == 0)
                break;
            if (remMaxSteps <= gMtStepsLimit) {
                solveSerial(s, bm);
                break;
            } else {
                State s1, s2;
                s.split(s1, s2);

#if 1            
                std::thread t1(solve, std::ref(s1), std::ref(bm));
                std::thread t2(solve, std::ref(s2), std::ref(bm));
                t1.join();
                t2.join();
#else
                solve(s1, bm, eps);
                solve(s2, bm, eps);
#endif
                s.merge(s1);
                s.merge(s2);
            }
        }
    }
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
    gNumRecUpdates = 0;
    s.mMaxSteps = gMaxStepsTotal;
    s.mProcs = gProcs;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
#if 0   
    solveSerial(s, bm, eps);
#else    
    solve(s, bm);
#endif
    end = std::chrono::system_clock::now();
    long long int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end - start)).count();
    *gOutStream << "Time: " << mseconds << " microsecond\n";
    *gOutStream << "Time per subproblem: " << (double) mseconds / (double) s.mSteps << " miscroseconds." << std::endl;
    if (s.mSteps >= gMaxStepsTotal) {
        *gOutStream << "Failed to converge in " << gMaxStepsTotal << " steps\n";
    } else {
        *gOutStream << "Converged in " << s.mSteps << " steps\n";
    }
    *gOutStream << "Number of record updates: " << gNumRecUpdates << "\n";
    *gOutStream << "BnB found = " << gRecv << std::endl;
    *gOutStream << " at x [ ";
    std::copy(gRecord.begin(), gRecord.end(), std::ostream_iterator<double>(*gOutStream, " "));
    *gOutStream << "]\n";
    gStat.emplace_back(mseconds, s.mSteps);
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
            std::cout << b->getDesc() << "\n";
        }
        return 0;
    } else if (argc == 8 || argc == 9) {
        nruns = atoi(argv[1]);
        bench = argv[2];
        gKnrec = argv[3];
        gEps = atof(argv[4]);
        gMaxStepsTotal = atoi(argv[5]);
        gProcs = atoi(argv[6]);
        gMtStepsLimit = atoi(argv[7]);
        if (argc == 9) {
            gOutStream = new BnbStream(nullptr);
        }
    } else {
        std::cerr << "Usage: " << argv[0] << " number_of_runs name_of_bench knrec|unknrec eps max_steps virtual_procs_number parallel_steps_limit [statonly]\n";
        std::cerr << "or to list benchmarks run:\n";
        std::cerr << argv[0] << " list\n";
        return -1;
    }
    *gOutStream << "Simple PBnB solver with np = " << gProcs << ", mtStepsLimit =  " << gMtStepsLimit << ", maxStepsTotal = " << gMaxStepsTotal << std::endl;
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
