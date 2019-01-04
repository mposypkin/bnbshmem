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

std::vector<BnBStat> stat;

//const std::memory_order morder = std::memory_order_seq_cst;
const std::memory_order morder = std::memory_order_relaxed;

#define EXCHNAGE_OPER compare_exchange_strong
//#define EXCHNAGE_OPER compare_exchange_weak

struct State {

    void merge(const State& s) {
        mSteps += s.mSteps;
        if (s.mRecordVal < mRecordVal) {
            mRecordVal = s.mRecordVal;
            mRecord = s.mRecord;
        }
        mPool.insert(mPool.end(), s.mPool.begin(), s.mPool.end());
    }

    void split(State& s1, State& s2) {
        const int remMaxSteps = mMaxSteps - mSteps;
        s1.mMaxSteps = remMaxSteps / 2;
        s2.mMaxSteps = remMaxSteps - s1.mMaxSteps;
        s1.mProcs = mProcs / 2;
        s2.mProcs = mProcs - s1.mProcs;

        s1.mRecord = mRecord;
        s2.mRecord = mRecord;

        s1.mRecordVal = mRecordVal;
        s2.mRecordVal = mRecordVal;
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

    double mRecordVal;

    std::vector<double> mRecord;

    std::vector<Box> mPool;

    int mMaxSteps;

    int mSteps = 0;

    int mProcs;
};

std::ostream& operator<<(std::ostream & out, const State s) {
    out << "\"recval\" : " << s.mRecordVal << "\n";
    out << "\"record\" : [";
    for (int i = 0; i < s.mRecord.size(); i++) {
        out << s.mRecord[i];
        if (i != s.mRecord.size() - 1)
            out << ", ";
    }
    out << "]\n";
    out << "\"steps\" :" << s.mSteps << "\n";
    out << "\"max steps\" :" << s.mMaxSteps << "\n";
    return out;
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
        double rv = gRecv.load(morder);
        while (v < rv) {
            gRecv.EXCHNAGE_OPER(rv, v, morder);
            s.mRecordVal = v;
            s.mRecord = c;
        }
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
                double rv = gRecv.load(morder);
                while (v < rv) {
                    gRecv.EXCHNAGE_OPER(rv, v, morder);
                    s.mRecordVal = v;
                    s.mRecord = c;
                }
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
            //std::cout << s << "\n";
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
        s.mRecordVal = bm.getGlobMinY();
    } else {
        s.mRecordVal = std::numeric_limits<double>::max();
    }
    gRecv = s.mRecordVal;
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
    int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end - start)).count();
    std::cout << "Time: " << mseconds << " microsecond\n";
    std::cout << "Time per subproblem: " << (double) mseconds / (double) s.mSteps << " miscroseconds." << std::endl;
    if (s.mSteps >= gMaxStepsTotal) {
        std::cout << "Failed to converge in " << gMaxStepsTotal << " steps\n";
    } else {
        std::cout << "Converged in " << s.mSteps << " steps\n";
    }

    std::cout << "BnB found = " << s.mRecordVal << std::endl;
    std::cout << " at x [ ";
    std::copy(s.mRecord.begin(), s.mRecord.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n";
    stat.emplace_back((double) mseconds, s.mSteps);
    return s.mRecordVal;
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
    } else if (argc == 8) {
        nruns = atoi(argv[1]);
        bench = argv[2];
        gKnrec = argv[3];
        gEps = atof(argv[4]);
        gMaxStepsTotal = atoi(argv[5]);
        gProcs = atoi(argv[6]);
        gMtStepsLimit = atoi(argv[7]);
    } else {
        std::cerr << "Usage: " << argv[0] << " number_of_runs name_of_bench knrec|unknrec eps max_steps virtual_procs_number parallel_steps_limit\n";
        std::cerr << "or to list benchmarks run:\n";
        std::cerr << argv[0] << " list\n";
        return -1;
    }
    std::cout << "Simple PBnB solver with np = " << gProcs << ", mtStepsLimit =  " << gMtStepsLimit << ", maxStepsTotal = " << gMaxStepsTotal << std::endl;
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
