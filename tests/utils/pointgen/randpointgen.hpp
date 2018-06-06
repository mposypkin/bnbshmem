/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   randpointgen.hpp
 * Author: mposypkin
 *
 * Created on February 22, 2017, 1:32 PM
 */

#ifndef RANDPOINTGEN_HPP
#define RANDPOINTGEN_HPP


#include <random>
#include <ctime>
#include <box/box.hpp>

#include "pointgenerator.hpp"

namespace snowgoose {

    template <class FT> class RandomPointGenerator : public PointGenerator <FT> {
    public:

        /**
         * Constructor
         * @param box input box
         * @param maxnum maximal number of points to generate (negative means infinite number of points)
         * @param seed seed (0 - not use seed)
         */
        RandomPointGenerator(const snowgoose::Box<FT>& box, int maxnum = -1, int seed = 0) :
        mBox(box), mSeed(seed), mMaxNumber(maxnum) {
            mGen.seed(time(0));
            reset();
        }

        bool getPoint(FT* x) {
            mCurNumber++;
            if ((mMaxNumber >= 0) && (mCurNumber > mMaxNumber)) {
                return false;
            } else {
                const int n = mBox.mDim;
                for (int i = 0; i < n; i++) {
                    std::uniform_real_distribution<> urd(mBox.mA[i], mBox.mB[i]);
                    x[i] = urd(mGen);
                }
                return true;
            }
        }

        void reset() {
            if (mSeed != 0) {
                mGen.seed(mSeed);
            }
            mCurNumber = 0;
        }

        virtual std::string about() const {
            return "Random uniformly distribute point generator";
        }



    private:
        const snowgoose::Box<FT>& mBox;
        const int mSeed;
        std::mt19937 mGen;
        const int mMaxNumber;
        int mCurNumber;
    };
}

#endif /* RANDPOINTGEN_HPP */

