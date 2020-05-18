/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   newfile.hpp
 * Author: posypkin
 *
 * Created on December 31, 2018, 12:32 PM
 */

#ifndef BNBSTAT_HPP
#define BNBSTAT_HPP

#include <iostream>
#include <vector>

/**
 * Statistics for a parallel Branch-and-Bound run
 */
class BnBStat {
public:
    /**
     * Constructor
     * @param totalTime the total running time
     * @param steps the number of steps performed
     */
    BnBStat(long long int totalTime, long long int steps) :
    mTotalTime(totalTime), mSteps(steps), mTPS(totalTime / (double) steps) {
    }

    /**
     * The total wall running time 
     */
    const long long int mTotalTime;
    /**
     * The total steps
     */
    const long long int mSteps;
    /**
     * Time for one subproblem (step)
     */
    const double mTPS;
};

std::ostream& operator<<(std::ostream& os, const BnBStat& st)
{
    os << st.mTotalTime << ", " << st.mSteps << ", " << st.mTPS;
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<BnBStat>& vst)
{
    for(auto &s : vst)
    os << s << "\n";
    return os;
}
#endif /* BNBSTAT_HPP */

