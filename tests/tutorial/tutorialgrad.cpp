/*
 * Demonstrates the use of gradient
 */

#include <iostream>
#include <limits>
#include <random>
#include "testfuncs/benchmarks.hpp"


std::mt19937 gen;

using BM = Benchmark<double>;

void localSearch(const BM& bm, std::vector<double>& x) {
    constexpr double amin = .05;    
    constexpr double amax = 1.;    
    double a = amin;    
    constexpr int maxit = 10000;
    double vold = bm.calcFunc(x);
    const int n = x.size();
    std::vector<double> y(n);
    for (int j = 0; j < maxit; j ++) {
        ValDer<double> valgrad = bm.calcGrad(x);
        Grad<double> g = valgrad.grad();
        for (int i = 0; i < n; i++) {
            y[i] = x[i] - a * g[i];
            y[i] = std::max(y[i], bm.getBounds()[i].first);
            y[i] = std::min(y[i], bm.getBounds()[i].second);
        }
        double v = bm.calcFunc(y);
        //std::cout << v << ":" << vold << std::endl;
        if (v < vold) {
            std::swap(x, y);
            std::swap(vold, v);
            if(a < amax)
                a *= 1.1;
        } else if(a > amin) {
            a *= 0.5;
        } else
            break;
    }
}

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
        localSearch(bm, x);
        v = std::min(v, bm.calcFunc(x));
    }
    std::cout << "best value = " << v << std::endl;
    return v;
}

bool testBench(const BM& bm) {
    constexpr int ntries = 10;
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
        try{
            testBench(*bm);
        } catch(std::invalid_argument exc) {
            std::cout << "Benchmark " << bm->getDesc() << exc.what() << std::endl;
        }
    }
}
