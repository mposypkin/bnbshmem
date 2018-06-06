#ifndef BENCHMARK__HPP
#define BENCHMARK__HPP

#include <iostream>
#include <vector>

#include "interval/interval_air.hpp"
#include "derivatives/valder.hpp"
#include "derivatives/intervalder.hpp"
#include "expression/expr.hpp"


using namespace snowgoose::expression;
using namespace snowgoose::derivative;

template <class T>
using Bound = std::pair<T, T>;

/**
* Benchmark to describe global optimization problem
*/
template <class T>
class Benchmark
{
protected:
	const std::string mDesciption;
	const T mGlobMinY;
	const bool mMultiModal;
	std::vector<T> mGlobMinX;
	std::vector<Bound<T>> mBounds;

public:
	/**
	* description of global optimization problem
	* @return string
	*/
	std::string getDesc() const { return mDesciption; }
	/**
	* get size of box
	* @return size of box
	*/
	int getDim() const { return mBounds.size(); }
	/**
	* get multimodal
 	* @return multimodal true/false
	*/
	bool getMulMod() const { return mMultiModal; }	
	/**
	* get expected global minimum
 	* @return expected global minimum
	*/
	T getGlobMinY() const { return mGlobMinY; }
	/**
	* get where is global minimum
 	* @return where is expected global minimum
	*/
	const std::vector<T> getGlobMinX() const { return mGlobMinX; }	
	/**
	* get bounds of box
 	* @return bounds of box
	*/
	const std::vector<Bound<T>> getBounds() const { return mBounds; }
	/**
	* calculate value of function
	* @param point coordinate where value of function will be calculated
 	* @return value of function
	*/
	virtual T calcFunc(const std::vector<T>& point) const = 0;
	/**
	* calculate interval estimation of function
	* @param box box where interval estimation of function will be calculated
 	* @return interval estimation of function
	*/
        virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const = 0;
	/**
	* calculate gradient of function
	* @param point coordinate where gradient of function will be calculated
 	* @return gradient of function
	*/
	virtual ValDer<T> calcGrad(const std::vector<T>& point) const = 0; 
	/**
	* calculate interval estimation of function gradient
	* @param box box where interval estimation of function gradient will be calculated
 	* @return interval estimation of function gradient
	*/
	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const = 0;
        /**
         * Constructor
         * @param desc description of global optimization problem
         * @param globmin expected global minimum
         * @param mMultiModal is multimodal function
         */
	Benchmark(const std::string& desc, T globmin, bool mMultiModal = true);
	/**
	* print benchmark
	* @param out out stream
	* @param bm benchmark
 	* @return out stream
	*/
        template <class T2> friend std::ostream& operator<<(std::ostream &out, const Benchmark<T2> &bm);
};

template <class T> Benchmark<T>::Benchmark(const std::string& desc, T globmin, bool mMultiModal) :  mDesciption(desc), mGlobMinY(globmin), mMultiModal(mMultiModal) {}

template <class T2> std::ostream& operator<<(std::ostream &out, const Benchmark<T2> &bm)
{
	const char* boollut[2] = {"true", "false"};
	out << "desciption: " << bm.mDesciption << std::endl;
	out << "dimension: " << bm.getDim() << std::endl;
	out << "globMinY: " << bm.mGlobMinY << std::endl;
	out << "multiModal: " << boollut[bm.mMultiModal] << std::endl;
	out << "globMinX:"<< std::endl;
	for(int i=0; i < bm.getDim(); ++i)
		out << "x[" << i << "]:" <<  bm.mGlobMinX[i] << std::endl;
	out << "bounds:"<< std::endl;
	for(int i=0; i < bm.getDim(); ++i)
		out << "[ " << bm.mBounds[i].first << " , " << bm.mBounds[i].second << " ]:" << std::endl;
	return out;	
};

#endif


