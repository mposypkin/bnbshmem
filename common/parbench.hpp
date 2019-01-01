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
        add(std::make_shared<Cluster2D2Benchmark<double>>());
        add(std::make_shared<Hartman6Benchmark<double>>());
        add(std::make_shared<BiggsEXP6Benchmark<double>>());
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

