/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bnbiutils.hpp
 * Author: mikhail
 *
 * Created on January 3, 2019, 2:00 PM
 * Useful utilities for BnB interval solver
 */

#ifndef BNBIUTILS_HPP
#define BNBIUTILS_HPP

#include "parbench.hpp"

using BM = Benchmark<double>;

using Box = std::vector<Interval<double>>;



/**
 * Width of the interval
 * @param I
 * @return width of the interval
 */
double wid(const Interval<double>& I) {
    return I.rb() - I.lb();
}

/**
 * Computes the center of a box
 * @param ibox the box
 * @param c the computed center
 */
void mid(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}

/**
 * Splits a box into two new boxes
 * @param ibox input box to split
 * @param v the vector to push results to
 */
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

#endif /* BNBIUTILS_HPP */

