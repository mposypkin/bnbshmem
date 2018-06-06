/*
 * A simple interval-based bnb solver
 */

/* 
 * File:   bnb_test_prelim.cpp fork from tutorialbnb.cpp
 * Author: Gorchakov A.Y.
 *
 * Created on Apr 18, 2018, 4:12 PM
 */

#include <iostream>
#include <limits>
#include <random>
#include <algorithm>
#include <vector>
#include <iterator>
#include "testfuncs/benchmarks.hpp"
#include <chrono>
#include <omp.h>
#include <sys/time.h>
double MyGetTime(){
struct timeval t0;
gettimeofday(&t0, NULL);
return t0.tv_sec + (t0.tv_usec / 1000000.0);
}
#include <algorithm>

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;

double len(const Interval<double>& I) {
    return I.rb() - I.lb();
}

void split(const Box& ibox, std::vector<Box>& v) {
    auto result = std::max_element(ibox.begin(), ibox.end(),
            [](const Interval<double>& f, const Interval<double>& s) {
                return len(f) < len(s);
            });
    const int i = result - ibox.begin();
    const double maxlen = len(ibox[i]);
    Box b1(ibox);
    Interval<double> ilow(ibox[i].lb(), ibox[i].lb() + 0.5 * maxlen);
    b1[i] = ilow;
    Box b2(ibox);
    Interval<double> iupper(ibox[i].lb() + 0.5 * maxlen, ibox[i].rb());
    b2[i] = iupper;
    v.push_back(std::move(b1));
    v.push_back(std::move(b2));
}

void getCenter(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}

double findMin(const BM& bm, double eps, long long int maxstep, int n_thr) {
    
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    std::vector<std::vector<Box>> pool(n_thr);
    std::vector<std::vector<Box>> pool_new(n_thr);
    pool[0].push_back(ibox);
    std::vector<double> c(dim);
    std::vector<double> recordVec(dim);
    double recordVal = std::numeric_limits<double>::max();
    
    long long int step = 0;
    long long pool_size = 0;
    auto start = std::chrono::system_clock::now();
    while (1) {
#pragma omp parallel shared(pool, pool_new, recordVec, recordVal, step, pool_size) firstprivate(c) num_threads(n_thr)
{
#pragma omp for 
        for(int i = 0; i < n_thr; i++) pool_new[i].clear();

        for(int i = 0; i < n_thr; i++) 
#pragma omp for nowait schedule(dynamic)
         for(long long int j = 0; j < pool[i].size(); j++) {
            int thread_num = omp_get_thread_num();
            Box b = pool[i][j];
            getCenter(b, c);
            double v = bm.calcFunc(c);
#pragma omp critical
{
            if (v < recordVal) {
                recordVal = v;
                recordVec = c;
            }
}
            auto lb = bm.calcInterval(b).lb();
            double rv;
#pragma omp atomic read
            rv = recordVal;
            if (lb <= rv - eps && pool_size + step < maxstep) 
               split(b, pool_new[thread_num]);
            
        }
}
        for(int i = 0; i < n_thr; i++) step += pool[i].size();
        if(step >= maxstep) break;
        pool_size = 0;
        for(int i = 0; i < n_thr; i++) pool_size += pool_new[i].size();
        if(pool_size == 0) break;
        if(pool_size + step >= maxstep) {
           pool_size = 0;
           for(int i = 0; i < n_thr; i++) {
              if(pool_size + pool_new[i].size() + step > maxstep) {
                 pool_new[i].erase(pool_new[i].begin() + maxstep - step - pool_size, pool_new[i].end());
              }
              pool_size += pool_new[i].size();
           }
        }
        
        pool.swap(pool_new);
    }
    auto end = std::chrono::system_clock::now();
    unsigned long int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end-start)).count();
    std::cout << bm.getDesc() << "," << bm.getDim() << "," << mseconds/1000000.0 << "," << step << "," << ((double) mseconds / (double) step) << "," << ((step >= maxstep)?"Fail":"OK") << "\n";
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
    double v = findMin(bm, eps, 10000000, 64);
    double diff = v - bm.getGlobMinY();
    if(diff > eps) {
        std::cout << "BnB failed for " << bm.getDesc()  << " benchmark " << std::endl;
    }
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;
}

main() {
    Benchmarks<double> tests;
    char *arrayDesc[] = {
"Mishra 8 function",
"Hartman 6 function",
"Langerman-5 function",	
"Biggs EXP5 Function", 
"Biggs EXP5 Function",
"Goldstein Price function",
"Powell Singular 2",
//"Mishra 8 function",
"Colville Function",	
"Dolan function",	
"Helical Valley function",
"Shubert 2 function",
"Shubert function",	
"Hosaki function",	
"Quintic function",	
"Biggs EXP4 Function",
"Deckkers-Aarts function",	
"Trid 10 function",
"Egg Holder function",
"Whitley function",	
"Trid 6 function",
"Trigonometric 1 function",
"Hansen function",	
"Schwefel 2.36 function",
"Chichinadze Function"};	

    Benchmarks<double> tests25;
    tests25.clear();
    for(auto bm : tests) {
       for(int i = 0; i < 24; i++) 
          if(bm->getDesc() == arrayDesc[i]) { 
             tests25.add(bm);
             break;
          }
    }

double t1 = MyGetTime();
    for (auto bm : tests25) {
        testBench(*bm);
    }
double t2 = MyGetTime();
printf("Time testBench=%lf\n", t2-t1);
}
