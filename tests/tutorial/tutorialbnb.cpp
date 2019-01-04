/*
 * A simple interval-based bnb solver
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
#include "testfuncs/benchmarks.hpp"

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;

double wid(const Interval<double>& I) {
    return I.rb() - I.lb();
}

void split(const Box& ibox, std::vector<Box>& v) {
    auto result = std::max_element(ibox.begin(), ibox.end(),
            [](const Interval<double>& f, const Interval<double>& s) {
                return wid(f) < wid(s);
            });
    const int i = result - ibox.begin();
    const double maxlen = wid(ibox[i]);
    Box b1(ibox);
    Interval<double> ilow(ibox[i].lb(), ibox[i].lb() + 0.5 * maxlen);
    b1[i] = ilow;
    Box b2(ibox);
    Interval<double> iupper(ibox[i].lb() + 0.5 * maxlen, ibox[i].rb());
    b2[i] = iupper;
    v.push_back(std::move(b1));
    v.push_back(std::move(b2));
}

void mid(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}

double findMin(const BM& bm, double eps, int maxstep) {
    
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    std::vector<Box> pool;
    pool.push_back(ibox);
    std::vector<double> c(dim);
    std::vector<double> recordVec(dim);
    double recordVal = std::numeric_limits<double>::max();

    
    int step = 0;
    while (!pool.empty()) {
        if(step > maxstep)
            break;
        step++;
        Box b = pool.back();
        pool.pop_back();        
        mid(b, c);
        double v = bm.calcFunc(c);
        if (v < recordVal) {
            recordVal = v;
            recordVec = c;
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= recordVal - eps) {
            split(b, pool);
        }
    }
    if(step > maxstep) {
        std::cout << "Failed to converge in " << maxstep << " steps\n";
    } else {
        std::cout << "Converged in " << step << " steps\n";
    }
    
    std::cout << "BnB found = " << recordVal << std::endl;
    std::cout << " at x [ " ;
    std::copy(recordVec.begin(), recordVec.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n" ;
    return recordVal;
}

bool testBench(const BM& bm) {
    constexpr double eps = 0.1;
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm;
    double v = findMin(bm, eps, 1000000);
    double diff = v - bm.getGlobMinY();
    if(diff > eps) {
        std::cout << "BnB failed for " << bm.getDesc()  << " benchmark " << std::endl;
    }
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;
}

main() {
    ZakharovBenchmark<double> zb(3);
    testBench(zb);
    Benchmarks<double> tests;
    for (auto bm : tests) {
        testBench(*bm);
    }
}
