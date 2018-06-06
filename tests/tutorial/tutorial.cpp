/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   tutorial.cpp
 *
 * Created on November 20, 2017, 9:13 AM
 */
#include <iostream>
#include <limits>
#include <random>
#include "testfuncs/benchmarks.hpp"


std::mt19937 gen;

using BM = Benchmark<double>;

double findMin(const BM& bm, int cnt) {
    const int dim = bm.getDim();
    double v = std::numeric_limits<double>::max();
    for (int k = 0; k < cnt; k++) {
        std::vector<double> x(dim);
        for (int i = 0; i < dim; i++) {
            double a = bm.getBounds()[i].first;
            double b = bm.getBounds()[i].second;
            std::uniform_real_distribution<> urd(a, b);
            x[i] = urd(gen);
        }
        v = std::min(v, bm.calcFunc(x));
    }
    std::cout << "best value = " << v << std::endl;
    return v;
}

bool testBench(const BM& bm) {
    const int ntries = 1000;
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm;
    double v = findMin(bm, ntries);
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
