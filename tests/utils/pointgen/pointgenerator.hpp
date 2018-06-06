/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   pointgenerator.hpp
 * Author: mposypkin
 *
 * Created on February 22, 2017, 12:41 PM
 */

#ifndef POINTGENERATOR_HPP
#define POINTGENERATOR_HPP

#include <box/box.hpp>
#include <common/sgerrcheck.hpp>

namespace snowgoose {

    /**
     * Abstract class for generating sequence of points 
     */
    template <class FT> class PointGenerator {
    public:
        /**
         * Get next point
         * @param x new point to fill in
         * @return true if the point was successfully generated, false - otherwise (no more points)
         */
        virtual bool getPoint(FT* x) = 0;

        /**
         * Resets the generator
         */
        virtual void reset() {
            SG_ERROR_REPORT("Resetting the generator is not implemented\n");
        }
        
        /**
         * Description of the point generator
         * @return point generator description
         */
        virtual std::string about() const {
            return "Pointe generator with empty description";
        }
    };

}

#endif /* POINTGENERATOR_HPP */

