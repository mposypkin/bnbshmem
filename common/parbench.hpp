/*
 * To change this license header, choose License Headers in Project Properties.
 * Parallel benchmarks
 * and open the template in the editor.
 */

/* 
 * File:   parbench.hpp
 * Author: mposypkin
 *
 * Created on June 8, 2018, 10:30 AM
 */

#ifndef PARBENCH_HPP
#define PARBENCH_HPP

#include <testfuncs/benchmarks.hpp>

template <class T>
class ParBenchmarks {
private:
    std::vector<PtrBench<T>> vbm;
public:

    ParBenchmarks() {
        fill();
    }

    void fill() {
        clear();
        add(std::make_shared<Ackley1Benchmark<double>>(6));
        add(std::make_shared<Ackley2Benchmark<double>>(6));
        add(std::make_shared<BrownBenchmark<double>>(6));
        add(std::make_shared<ChungReynoldsBenchmark<double>>(6));
        add(std::make_shared<ExponentialBenchmark<double>>(6));
        add(std::make_shared<GriewankBenchmark<double>>(6));
        add(std::make_shared<PowellSingular2Benchmark<double>>(6));
        add(std::make_shared<QingBenchmark<double>>());
        add(std::make_shared<QuinticBenchmark<double>>(6));
        add(std::make_shared<RosenbrockBenchmark<double>>(8));
        add(std::make_shared<SchumerSteiglitzBenchmark<double>>(3));
        add(std::make_shared<SchafferF6Benchmark<double>>(6));
        add(std::make_shared<Schwefel2_22Benchmark<double>>(3));
        add(std::make_shared<StrechedVSineWaveBenchmark<double>>(6));
        add(std::make_shared<Trid6Benchmark<double>>());
        add(std::make_shared<Trid10Benchmark<double>>());
        add(std::make_shared<Trigonometric1Benchmark<double>>(6));
        add(std::make_shared<Trigonometric2Benchmark<double>>(6));
        add(std::make_shared<WWavyBenchmark<double>>(6));
        add(std::make_shared<WhitleyBenchmark<double>>(6));
        add(std::make_shared<XinSheYang2Benchmark<double>>(6));
        add(std::make_shared<XinSheYang3Benchmark<double>>(6));
        add(std::make_shared<XinSheYang4Benchmark<double>>(6));
        add(std::make_shared<ZakharovBenchmark<double>>(6));
        add(std::make_shared<Cluster2D1Benchmark<double>>());        
        add(std::make_shared<Cluster2D2Benchmark<double>>());        
    }

    void clear() {
        vbm.clear();
    }

    void add(const PtrBench<T> &ptrBenchmark) {
        vbm.push_back(ptrBenchmark);
    }

    void print() const {
        for (PtrBench<T> ptrBench : * this)
            std::cout << *ptrBench;
    }

    typename std::vector<PtrBench<T>>::iterator begin() {
        return vbm.begin();
    }

    typename std::vector<PtrBench<T>>::iterator end() {
        return vbm.end();
    }

    typename std::vector<PtrBench<T>>::const_iterator begin() const {
        return vbm.begin();
    }

    typename std::vector<PtrBench<T>>::const_iterator end() const {
        return vbm.end();
    }
};


#endif /* PARBENCH_HPP */

