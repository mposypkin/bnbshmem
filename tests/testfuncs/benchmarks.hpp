#ifndef BENCHMARKS__HPP
#define BENCHMARKS__HPP

#include <iostream>
#include <vector>
#include <memory>

#include "benchmark.hpp"
#include "expression/expr.hpp"
#include "testfuncs.hpp"


using namespace snowgoose::expression;
using namespace snowgoose::derivative;
using namespace OPTITEST;

template <class T>
class Ackley1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Ackley1<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Ackley1<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Ackley1<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Ackley1<IntervalDer<T>>(this->getDim()), box);
	}

	Ackley1Benchmark(int size) : Benchmark<T>("Ackley 1 function",0) {
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-32,32});
		}
	}
};


template <class T>
class Ackley2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Ackley2<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Ackley2<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Ackley2<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Ackley2<IntervalDer<T>>(this->getDim()), box);
	}

	Ackley2Benchmark(int size) : Benchmark<T>("Ackley 2 function",-200, false) {
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-32,32});
		}
	}
};

template <class T>
class Ackley3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Ackley3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Ackley3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Ackley3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Ackley3<IntervalDer<T>>(), box);
	}

	Ackley3Benchmark() : Benchmark<T>("Ackley 3 function",82.4617, false)	{
		this->mGlobMinX.push_back(-32);
		this->mGlobMinX.push_back(-32);
		this->mBounds.push_back({-32,32});
		this->mBounds.push_back({-32,32});
	}
};
template <class T>
class Ackley4Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Ackley4<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Ackley4<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Ackley4<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Ackley4<IntervalDer<T>>(), box);
	}

	Ackley4Benchmark() : Benchmark<T>("Ackley 4 function",-4.99995)	{
		this->mGlobMinX.push_back(-23.5625);
		this->mGlobMinX.push_back(-29.0625);
		this->mBounds.push_back({-32,32});
		this->mBounds.push_back({-32,32});
	}
};
template <class T>
class AdjimanBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Adjiman<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Adjiman<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Adjiman<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Adjiman<IntervalDer<T>>(), box);
	}

	AdjimanBenchmark() : Benchmark<T>("Adjiman function",-2.02181)	{
		this->mGlobMinX.push_back(2);
		this->mGlobMinX.push_back(0.10578);
		this->mBounds.push_back({-1,2});
		this->mBounds.push_back({-1,1});
	}
};
template <class T>
class Alpine1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Alpine1<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Alpine1<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Alpine1<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Alpine1<IntervalDer<T>>(this->getDim()), box);
	}

	Alpine1Benchmark(int size) : Benchmark<T>("Alpine 1 function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-10,10});
		}
	}
};

template <class T>
class Alpine2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Alpine2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Alpine2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Alpine2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Alpine2<IntervalDer<T>>(), box);
	}

	Alpine2Benchmark() : Benchmark<T>("Alpine 2 function",-17.2124)	{
		this->mGlobMinX.push_back(7.91705);
		this->mGlobMinX.push_back(7.91705);
		this->mGlobMinX.push_back(4.81584);
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
	}
};
template <class T>
class BartelsConnBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(BartelsConn<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(BartelsConn<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(BartelsConn<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(BartelsConn<IntervalDer<T>>(), box);
	}

	BartelsConnBenchmark() : Benchmark<T>("Bartels Conn Function",1) {
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class BealeBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Beale<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Beale<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Beale<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Beale<IntervalDer<T>>(), box);
	}

	BealeBenchmark() : Benchmark<T>("Beale function",0, false)	{
		this->mGlobMinX.push_back(3);
		this->mGlobMinX.push_back(0.5);
		this->mBounds.push_back({-4.5,4.5});
		this->mBounds.push_back({-4.5,4.5});
	}
};

template <class T>
class BiggsEXP2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(BiggsExpr2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(BiggsExpr2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(BiggsExpr2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(BiggsExpr2<IntervalDer<T>>(), box);
	}

	BiggsEXP2Benchmark() : Benchmark<T>("Biggs EXP2 Function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(10);
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
	}
};

template <class T>
class BiggsEXP3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(BiggsExpr3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(BiggsExpr3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(BiggsExpr3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(BiggsExpr3<IntervalDer<T>>(), box);
	}

	BiggsEXP3Benchmark() : Benchmark<T>("Biggs EXP3 Function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(10);
		this->mGlobMinX.push_back(5);
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
	}
};
template <class T>
class BiggsEXP4Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(BiggsExpr4<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(BiggsExpr4<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(BiggsExpr4<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(BiggsExpr4<IntervalDer<T>>(), box);
	}

	BiggsEXP4Benchmark() : Benchmark<T>("Biggs EXP4 Function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(10);
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(5);
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
	}
};
template <class T>
class BiggsEXP5Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(BiggsExpr5<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(BiggsExpr5<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(BiggsExpr5<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(BiggsExpr5<IntervalDer<T>>(), box);
	}

	BiggsEXP5Benchmark() : Benchmark<T>("Biggs EXP5 Function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(10);
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(5);
		this->mGlobMinX.push_back(4);
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
	}
};
template <class T>
class BiggsEXP6Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(BiggsExpr6<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(BiggsExpr6<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(BiggsExpr6<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(BiggsExpr6<IntervalDer<T>>(), box);
	}

    BiggsEXP6Benchmark() : Benchmark<T>("Biggs EXP6 Function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(10);
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(5);
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(3);
		this->mBounds.push_back({-20,20});
		this->mBounds.push_back({-20,20});
		this->mBounds.push_back({-20,20});
		this->mBounds.push_back({-20,20});
		this->mBounds.push_back({-20,20});
		this->mBounds.push_back({-20,20});
	}
};

template <class T>
class BirdBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Bird<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Bird<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Bird<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Bird<IntervalDer<T>>(), box);
	}

	BirdBenchmark() : Benchmark<T>("Bird Function",-106.765)	{
		this->mGlobMinX.push_back(4.70104);
		this->mGlobMinX.push_back(3.15294);
		this->mBounds.push_back({-6.28319,6.28319});
		this->mBounds.push_back({-6.28319,6.28319});
	}
};
template <class T>
class Bohachevsky1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Bohachevsky1<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Bohachevsky1<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Bohachevsky1<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Bohachevsky1<IntervalDer<T>>(), box);
	}

	Bohachevsky1Benchmark() : Benchmark<T>("Bohachevsky 1 Function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};

template <class T>
class Bohachevsky2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Bohachevsky2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Bohachevsky2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Bohachevsky2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Bohachevsky2<IntervalDer<T>>(), box);
	}

	Bohachevsky2Benchmark() : Benchmark<T>("Bohachevsky 2 Function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class Bohachevsky3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Bohachevsky3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Bohachevsky3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Bohachevsky3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Bohachevsky3<IntervalDer<T>>(), box);
	}

	Bohachevsky3Benchmark() : Benchmark<T>("Bohachevsky 3 Function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class BoothBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Booth<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Booth<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Booth<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Booth<IntervalDer<T>>(), box);
	}

	BoothBenchmark() : Benchmark<T>("Booth Function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(3);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class BoxBettsQuadraticSumBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(BoxBettsQuadraticSum<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(BoxBettsQuadraticSum<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(BoxBettsQuadraticSum<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(BoxBettsQuadraticSum<IntervalDer<T>>(), box);
	}

	BoxBettsQuadraticSumBenchmark() : Benchmark<T>("Box-Betts Quadratic Sum Function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(10);
		this->mGlobMinX.push_back(1);
		this->mBounds.push_back({0.9,1.2});
		this->mBounds.push_back({9,11.2});
		this->mBounds.push_back({0.9,1.2});
	}
};

template <class T>
class BradBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Brad<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Brad<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Brad<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Brad<IntervalDer<T>>(), box);
	}

	BradBenchmark() : Benchmark<T>("Brad function",6.93602)	{
		this->mGlobMinX.push_back(-0.249878);
		this->mGlobMinX.push_back(2.49992);
		this->mGlobMinX.push_back(2.49985);
		this->mBounds.push_back({-0.25,0.25});
		this->mBounds.push_back({0.01,2.5});
		this->mBounds.push_back({0.01,2.5});
	}
};


template <class T>
class BraninRCOSBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(BraninRCOS<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(BraninRCOS<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(BraninRCOS<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(BraninRCOS<IntervalDer<T>>(), box);
	}

	BraninRCOSBenchmark() : Benchmark<T>("Branin RCOS Function",0.397887)	{
		this->mGlobMinX.push_back(-3.14159);
		this->mGlobMinX.push_back(12.275);
		this->mBounds.push_back({-5,10});
		this->mBounds.push_back({0,15});
	}
};
template <class T>
class BraninRCOS2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(BraninRCOS2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(BraninRCOS2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(BraninRCOS2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(BraninRCOS2<IntervalDer<T>>(), box);
	}

	BraninRCOS2Benchmark() : Benchmark<T>("Branin RCOS 2 Function",-0.179891)	{
		this->mGlobMinX.push_back(-3.19699);
		this->mGlobMinX.push_back(12.5263);
		this->mBounds.push_back({-5,10});
		this->mBounds.push_back({0,15});
	}
};
template <class T>
class BrentBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Brent<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Brent<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Brent<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Brent<IntervalDer<T>>(), box);
	}

	BrentBenchmark() : Benchmark<T>("Brent Function",0)	{
		this->mGlobMinX.push_back(-10);
		this->mGlobMinX.push_back(-10);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class BrownBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Brown<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Brown<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Brown<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Brown<IntervalDer<T>>(this->getDim()), box);
	}

    BrownBenchmark(int size) : Benchmark<T>("Brown function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-1,4});
		}
	}
};

template <class T>
class Bukin2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Bukin2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Bukin2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Bukin2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Bukin2<IntervalDer<T>>(), box);
	}

	Bukin2Benchmark() : Benchmark<T>("Bukin 2 Function",0)	{
		this->mGlobMinX.push_back(-10);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-15,-5});
		this->mBounds.push_back({-3,3});
	}
};
template <class T>
class Bukin4Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Bukin4<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Bukin4<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Bukin4<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Bukin4<IntervalDer<T>>(), box);
	}

	Bukin4Benchmark() : Benchmark<T>("Bukin 4 Function",0)	{
		this->mGlobMinX.push_back(-10);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-15,-5});
		this->mBounds.push_back({-3,3});
	}
};
template <class T>
class Bukin6Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Bukin6<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Bukin6<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Bukin6<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Bukin6<IntervalDer<T>>(), box);
	}

	Bukin6Benchmark() : Benchmark<T>("Bukin 6 Function",0)	{
		this->mGlobMinX.push_back(-10);
		this->mGlobMinX.push_back(1);
		this->mBounds.push_back({-15,-5});
		this->mBounds.push_back({-3,3});
	}
};
template <class T>
class CamelSixHumpBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(CamelSixHump<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(CamelSixHump<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(CamelSixHump<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(CamelSixHump<IntervalDer<T>>(), box);
	}

	CamelSixHumpBenchmark() : Benchmark<T>("Camel Function – Six Hump Function",-1.0316)	{
		this->mGlobMinX.push_back(-0.0898);
		this->mGlobMinX.push_back(0.7126);
		this->mBounds.push_back({-5,5});
		this->mBounds.push_back({-5,5});
	}
};
template <class T>
class CamelThreeHumpBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(CamelThreeHump<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(CamelThreeHump<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(CamelThreeHump<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(CamelThreeHump<IntervalDer<T>>(), box);
	}

	CamelThreeHumpBenchmark() : Benchmark<T>("Camel Function – Three Hump Function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-5,5});
		this->mBounds.push_back({-5,5});
	}
};

template <class T>
class ChichinadzeBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Chichinadze<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Chichinadze<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Chichinadze<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Chichinadze<IntervalDer<T>>(), box);
	}

	ChichinadzeBenchmark() : Benchmark<T>("Chichinadze Function",-42.9444)	{
		this->mGlobMinX.push_back(6.18987);
		this->mGlobMinX.push_back(0.5);
		this->mBounds.push_back({-30,30});
		this->mBounds.push_back({-30,30});
	}
};
template <class T>
class ChungReynoldsBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(ChungReynolds<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(ChungReynolds<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(ChungReynolds<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(ChungReynolds<IntervalDer<T>>(this->getDim()), box);
	}

	ChungReynoldsBenchmark(int size) : Benchmark<T>("Chung Reynolds function",0, false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-100,100});
		}
	}
};
template <class T>
class ColvilleBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Colville<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Colville<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Colville<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Colville<IntervalDer<T>>(), box);
	}

	ColvilleBenchmark() : Benchmark<T>("Colville Function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(1);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class ComplexBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Complex<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Complex<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Complex<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Complex<IntervalDer<T>>(), box);
	}

	ComplexBenchmark() : Benchmark<T>("Complex function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-2,2});
		this->mBounds.push_back({-2,2});
	}
};
template <class T>
class CosineMixtureBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(CosineMixture<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(CosineMixture<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(CosineMixture<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(CosineMixture<IntervalDer<T>>(), box);
	}

	CosineMixtureBenchmark() : Benchmark<T>("Cosine Mixture function",-0.2){
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-1,1});
		this->mBounds.push_back({-1,1});
	}
};
template <class T>
class CrossInTrayBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(CrossInTray<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(CrossInTray<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(CrossInTray<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(CrossInTray<IntervalDer<T>>(), box);
	}

	CrossInTrayBenchmark() : Benchmark<T>("Cross In Tray function",-2.06261){
		this->mGlobMinX.push_back(1.34941);
		this->mGlobMinX.push_back(1.34941);
		this->mBounds.push_back({-15,15});
		this->mBounds.push_back({-15,15});
	}
};
template <class T>
class CrossLegBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(CrossLeg<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(CrossLeg<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(CrossLeg<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(CrossLeg<IntervalDer<T>>(), box);
	}

	CrossLegBenchmark() : Benchmark<T>("Cross-Leg function",-1)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class CubeBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Cube<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Cube<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Cube<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Cube<IntervalDer<T>>(), box);
	}

	CubeBenchmark() : Benchmark<T>("Cube function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(1);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class DavisBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Davis<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Davis<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Davis<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Davis<IntervalDer<T>>(), box);
	}

	DavisBenchmark() : Benchmark<T>("Davis's function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class Deb1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Deb1<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Deb1<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Deb1<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Deb1<IntervalDer<T>>(this->getDim()), box);
	}

	Deb1Benchmark(int size) : Benchmark<T>("Deb 1 function",-1)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(-0.9);
			this->mBounds.push_back({-1,1});
		}
	}
};
template <class T>
class DeckkersAartsBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(DeckkersAarts<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(DeckkersAarts<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(DeckkersAarts<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(DeckkersAarts<IntervalDer<T>>(), box);
	}

	DeckkersAartsBenchmark() : Benchmark<T>("Deckkers-Aarts function",-24776.5)	{
		this->mGlobMinX.push_back(-0.000532866);
		this->mGlobMinX.push_back(14.9483);
		this->mBounds.push_back({-20,20});
		this->mBounds.push_back({-20,20});
	}
};
template <class T>
class DixonPriceBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(DixonPrice<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(DixonPrice<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(DixonPrice<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(DixonPrice<IntervalDer<T>>(), box);
	}

	DixonPriceBenchmark() : Benchmark<T>("Dixon & Price function",0,false)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(0.707107);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class DolanBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Dolan<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Dolan<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Dolan<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Dolan<IntervalDer<T>>(), box);
	}

	DolanBenchmark() : Benchmark<T>("Dolan function",-529.871)	{
		this->mGlobMinX.push_back(98.9643);
		this->mGlobMinX.push_back(100);
		this->mGlobMinX.push_back(100);
		this->mGlobMinX.push_back(99.2243);
		this->mGlobMinX.push_back(-0.249988);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class DropWaveBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(DropWave<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(DropWave<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(DropWave<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(DropWave<IntervalDer<T>>(), box);
	}

	DropWaveBenchmark() : Benchmark<T>("Drop-Wave function",-1)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-5.12,5.12});
		this->mBounds.push_back({-5.12,5.12});
	}
};
template <class T>
class EasomBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Easom<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Easom<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Easom<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Easom<IntervalDer<T>>(), box);
	}

	EasomBenchmark() : Benchmark<T>("Easom function",-1)	{
		this->mGlobMinX.push_back(3.14159);
		this->mGlobMinX.push_back(3.14159);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class EggCrateBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(EggCrate<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(EggCrate<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(EggCrate<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(EggCrate<IntervalDer<T>>(), box);
	}

	EggCrateBenchmark() : Benchmark<T>("Egg Crate function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-5,5});
		this->mBounds.push_back({-5,5});
	}
};
template <class T>
class EggHolderBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(EggHolder<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(EggHolder<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(EggHolder<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(EggHolder<IntervalDer<T>>(), box);
	}

	EggHolderBenchmark() : Benchmark<T>("Egg Holder function",-959.641)	{
		this->mGlobMinX.push_back(512);
		this->mGlobMinX.push_back(404.232);
		this->mBounds.push_back({-512,512});
		this->mBounds.push_back({-512,512});
	}
};
template <class T>
class ElAttarVidyasagarDuttBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(ElAttarVidyasagarDutt<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(ElAttarVidyasagarDutt<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(ElAttarVidyasagarDutt<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(ElAttarVidyasagarDutt<IntervalDer<T>>(), box);
	}

	ElAttarVidyasagarDuttBenchmark() : Benchmark<T>("El-Attar-Vidyasagar-Dutta function",1.71278,false)	{
		this->mGlobMinX.push_back(3.40919);
		this->mGlobMinX.push_back(-2.17143);
		this->mBounds.push_back({-500,500});
		this->mBounds.push_back({-500,500});
	}
};
template <class T>
class EngvallBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Engvall<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Engvall<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Engvall<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Engvall<IntervalDer<T>>(), box);
	}

	EngvallBenchmark() : Benchmark<T>("Engvall's function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-2000,2000});
		this->mBounds.push_back({-2000,2000});
	}
};
template <class T>
class Exp2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Exp2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Exp2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Exp2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Exp2<IntervalDer<T>>(), box);
	}

	Exp2Benchmark() : Benchmark<T>("Exp 2  function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(10);
		this->mBounds.push_back({0,20});
		this->mBounds.push_back({0,20});
	}
};
template <class T>
class ExponentialBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Exponential<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Exponential<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Exponential<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Exponential<IntervalDer<T>>(this->getDim()), box);
	}

	ExponentialBenchmark(int size) : Benchmark<T>("Exponential function",-1){
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-1,1});
		}
	}
};
template <class T>
class FreudensteinRothBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(FreudensteinRoth<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(FreudensteinRoth<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(FreudensteinRoth<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(FreudensteinRoth<IntervalDer<T>>(), box);
	}

	FreudensteinRothBenchmark() : Benchmark<T>("Freudenstein Roth  function",0)	{
		this->mGlobMinX.push_back(5);
		this->mGlobMinX.push_back(4);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class GoldsteinPriceBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(GoldsteinPrice<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(GoldsteinPrice<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(GoldsteinPrice<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(GoldsteinPrice<IntervalDer<T>>(), box);
	}

	GoldsteinPriceBenchmark() : Benchmark<T>("Goldstein Price function",3)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(-1);
		this->mBounds.push_back({-2,2});
		this->mBounds.push_back({-2,2});
	}
};
template <class T>
class GramacyLee2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(GramacyLee2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(GramacyLee2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(GramacyLee2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(GramacyLee2<IntervalDer<T>>(), box);
	}

	GramacyLee2Benchmark() : Benchmark<T>("GramacyLee's 2  function",-0.428882)	{
		this->mGlobMinX.push_back(-0.707107);
		this->mGlobMinX.push_back(-3.32453e-09);
		this->mBounds.push_back({-1.5,1.5});
		this->mBounds.push_back({-1.5,1.5});
	}
};
template <class T>
class GramacyLee3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(GramacyLee3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(GramacyLee3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(GramacyLee3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(GramacyLee3<IntervalDer<T>>(), box);
	}

	GramacyLee3Benchmark() : Benchmark<T>("GramacyLee's 3  function",-1.12687)	{
		this->mGlobMinX.push_back(-1.04083);
		this->mGlobMinX.push_back(-1.04083);
		this->mBounds.push_back({-1.5,1.5});
		this->mBounds.push_back({-1.5,1.5});
	}
};
template <class T>
class GriewankBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Griewank<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Griewank<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Griewank<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Griewank<IntervalDer<T>>(this->getDim()), box);
	}

	GriewankBenchmark(int size) : Benchmark<T>("Griewank function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-100,100});
		}
	}
};
template <class T>
class HansenBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Hansen<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Hansen<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Hansen<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Hansen<IntervalDer<T>>(), box);
	}

	HansenBenchmark() : Benchmark<T>("Hansen function",-176.542)	{
		this->mGlobMinX.push_back(-7.58989);
		this->mGlobMinX.push_back(-7.70831);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Hartman3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Hartman3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Hartman3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Hartman3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Hartman3<IntervalDer<T>>(), box);
	}

	Hartman3Benchmark() : Benchmark<T>("Hartman 3 function",-3.86278)	{
		this->mGlobMinX.push_back(0.114);
		this->mGlobMinX.push_back(0.556);
		this->mGlobMinX.push_back(0.852);
		this->mBounds.push_back({0,1});
		this->mBounds.push_back({0,1});
		this->mBounds.push_back({0,1});
	}
};
template <class T>
class Hartman6Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Hartman6<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Hartman6<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Hartman6<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Hartman6<IntervalDer<T>>(), box);
	}

	Hartman6Benchmark() : Benchmark<T>("Hartman 6 function",-3.32236)	{
		this->mGlobMinX.push_back(0.20169);
		this->mGlobMinX.push_back(0.150011);
		this->mGlobMinX.push_back(0.476874);
		this->mGlobMinX.push_back(0.275332);
		this->mGlobMinX.push_back(0.311652);
		this->mGlobMinX.push_back(0.657301);
		this->mBounds.push_back({0,1});
		this->mBounds.push_back({0,1});
		this->mBounds.push_back({0,1});
		this->mBounds.push_back({0,1});
		this->mBounds.push_back({0,1});
		this->mBounds.push_back({0,1});
	}
};
template <class T>
class HelicalValleyBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(HelicalValley<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(HelicalValley<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(HelicalValley<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(HelicalValley<IntervalDer<T>>(), box);
	}

	HelicalValleyBenchmark() : Benchmark<T>("Helical Valley function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({0.1,10});
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class HimmelblauBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Himmelblau<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Himmelblau<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Himmelblau<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Himmelblau<IntervalDer<T>>(), box);
	}

	HimmelblauBenchmark() : Benchmark<T>("Himmelblau function",0)	{
		this->mGlobMinX.push_back(3);
		this->mGlobMinX.push_back(2);
		this->mBounds.push_back({-5,5});
		this->mBounds.push_back({-5,5});
	}
};
template <class T>
class HosakiBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Hosaki<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Hosaki<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Hosaki<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Hosaki<IntervalDer<T>>(), box);
	}

	HosakiBenchmark() : Benchmark<T>("Hosaki function",-2.3458)	{
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(2);
		this->mBounds.push_back({0,5});
		this->mBounds.push_back({0,6});
	}
};
template <class T>
class JennrichSampsonBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(JennrichSampson<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(JennrichSampson<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(JennrichSampson<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(JennrichSampson<IntervalDer<T>>(), box);
	}

	JennrichSampsonBenchmark() : Benchmark<T>("Jennrich-Sampson function",124.361)	{
		this->mGlobMinX.push_back(0.257825);
		this->mGlobMinX.push_back(0.257825);
		this->mBounds.push_back({-1,1});
		this->mBounds.push_back({-1,1});
	}
};
template <class T>
class KeaneBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Keane<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Keane<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Keane<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Keane<IntervalDer<T>>(), box);
	}

	KeaneBenchmark() : Benchmark<T>("Keane function",-0.673668)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(1.39325);
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0.1,10});
	}
};
template <class T>
class Langerman5Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Langerman5<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Langerman5<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Langerman5<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Langerman5<IntervalDer<T>>(), box);
	}

	Langerman5Benchmark() : Benchmark<T>("Langerman-5 function",-0.965)	{
		this->mGlobMinX.push_back(8.074);
		this->mGlobMinX.push_back(8.777);
		this->mGlobMinX.push_back(3.467);
		this->mGlobMinX.push_back(1.86301);
		this->mGlobMinX.push_back(6.708);
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
	}
};
template <class T>
class LeonBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Leon<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Leon<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Leon<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Leon<IntervalDer<T>>(), box);
	}

	LeonBenchmark() : Benchmark<T>("Leon function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(1);
		this->mBounds.push_back({-1.2,1.2});
		this->mBounds.push_back({-1.2,1.2});
	}
};
template <class T>
class MatyasBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Matyas<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Matyas<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Matyas<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Matyas<IntervalDer<T>>(), box);
	}

	MatyasBenchmark() : Benchmark<T>("Matyas function",0,false)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class McCormickBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(McCormick<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(McCormick<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(McCormick<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(McCormick<IntervalDer<T>>(), box);
	}

	McCormickBenchmark() : Benchmark<T>("McCormick function",-1.9133)	{
		this->mGlobMinX.push_back(-0.547);
		this->mGlobMinX.push_back(-1.547);
		this->mBounds.push_back({-1.5,4});
		this->mBounds.push_back({-3,3});
	}
};
template <class T>
class MieleCantrellBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(MieleCantrell<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(MieleCantrell<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(MieleCantrell<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(MieleCantrell<IntervalDer<T>>(), box);
	}

	MieleCantrellBenchmark() : Benchmark<T>("Miele Cantrell function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(1);
		this->mBounds.push_back({-1,1});
		this->mBounds.push_back({-1,1});
		this->mBounds.push_back({-1,1});
		this->mBounds.push_back({-1,1});
	}
};
template <class T>
class Mishra3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Mishra3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Mishra3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Mishra3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Mishra3<IntervalDer<T>>(), box);
	}

	Mishra3Benchmark() : Benchmark<T>("Mishra 3 function",-0.184651)	{
		this->mGlobMinX.push_back(-8.46661);
		this->mGlobMinX.push_back(-9.99852);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Mishra4Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Mishra4<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Mishra4<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Mishra4<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Mishra4<IntervalDer<T>>(), box);
	}

	Mishra4Benchmark() : Benchmark<T>("Mishra 4 function",-0.199407)	{
		this->mGlobMinX.push_back(-9.94113);
		this->mGlobMinX.push_back(-9.99957);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Mishra5Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Mishra5<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Mishra5<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Mishra5<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Mishra5<IntervalDer<T>>(), box);
	}

	Mishra5Benchmark() : Benchmark<T>("Mishra 5 function",-1.01983)	{
		this->mGlobMinX.push_back(-1.98682);
		this->mGlobMinX.push_back(-10);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Mishra6Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Mishra6<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Mishra6<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Mishra6<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Mishra6<IntervalDer<T>>(), box);
	}

	Mishra6Benchmark() : Benchmark<T>("Mishra 6 function",-2.28395)	{
		this->mGlobMinX.push_back(2.88631);
		this->mGlobMinX.push_back(1.82326);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Mishra7Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Mishra7<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Mishra7<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Mishra7<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Mishra7<IntervalDer<T>>(), box);
	}

	Mishra7Benchmark() : Benchmark<T>("Mishra 7 function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(2);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Mishra8Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Mishra8<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Mishra8<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Mishra8<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Mishra8<IntervalDer<T>>(), box);
	}

	Mishra8Benchmark() : Benchmark<T>("Mishra 8 function",0)	{
		this->mGlobMinX.push_back(2);
		this->mGlobMinX.push_back(-3);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Mishra9Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Mishra9<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Mishra9<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Mishra9<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Mishra9<IntervalDer<T>>(), box);
	}

    Mishra9Benchmark() : Benchmark<T>("Mishra 9 function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(2);
		this->mGlobMinX.push_back(3);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class ParsopoulosBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Parsopoulos<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Parsopoulos<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Parsopoulos<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Parsopoulos<IntervalDer<T>>(), box);
	}

	ParsopoulosBenchmark() : Benchmark<T>("Parsopoulos function",0)	{
		this->mGlobMinX.push_back(1.5708);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-5,5});
		this->mBounds.push_back({-5,5});
	}
};
template <class T>
class PathologicalBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Pathological<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Pathological<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Pathological<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Pathological<IntervalDer<T>>(this->getDim()), box);
	}

	PathologicalBenchmark(int size) : Benchmark<T>("Pathological function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-100,100});
		}
	}
};
template <class T>
class PeriodicBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Periodic<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Periodic<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Periodic<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Periodic<IntervalDer<T>>(), box);
	}

	PeriodicBenchmark() : Benchmark<T>("Periodic function",0.9)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class PinterBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Pinter<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Pinter<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Pinter<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Pinter<IntervalDer<T>>(this->getDim()), box);
	}

	PinterBenchmark(int size) : Benchmark<T>("Pinter function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-10,10});
		}
	}
};
template <class T>
class PowellSingular2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(PowellSingular2<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(PowellSingular2<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(PowellSingular2<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(PowellSingular2<IntervalDer<T>>(this->getDim()), box);
	}

    PowellSingular2Benchmark(int size) : Benchmark<T>("Powell Singular 2_" + std::to_string(size) + "s function",0,false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-4,5});
		}
	}
};
template <class T>
class PowellSumBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(PowellSum<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(PowellSum<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(PowellSum<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(PowellSum<IntervalDer<T>>(this->getDim()), box);
	}

	PowellSumBenchmark(int size) : Benchmark<T>("Powell Sum function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-1,1});
		}
	}
};
template <class T>
class Price1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Price1<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Price1<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Price1<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Price1<IntervalDer<T>>(), box);
	}

	Price1Benchmark() : Benchmark<T>("Price 1 function",0)	{
		this->mGlobMinX.push_back(-5);
		this->mGlobMinX.push_back(-5);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Price2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Price2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Price2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Price2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Price2<IntervalDer<T>>(), box);
	}

	Price2Benchmark() : Benchmark<T>("Price 2 function",0.9)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Price3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Price3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Price3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Price3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Price3<IntervalDer<T>>(), box);
	}

	Price3Benchmark() : Benchmark<T>("Price 3 function",0)	{
		this->mGlobMinX.push_back(0.341308);
		this->mGlobMinX.push_back(0.116491);
		this->mBounds.push_back({-5,5});
		this->mBounds.push_back({-5,5});
	}
};
template <class T>
class Price4Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Price4<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Price4<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Price4<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Price4<IntervalDer<T>>(), box);
	}

    Price4Benchmark() : Benchmark<T>("Price 4 function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-500,500});
		this->mBounds.push_back({-500,500});
	}
};
template <class T>
class Problem02Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Problem02<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Problem02<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Problem02<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Problem02<IntervalDer<T>>(), box);
	}

	Problem02Benchmark() : Benchmark<T>("Problem 02 function",-1.8996)	{
		this->mGlobMinX.push_back(5.14574);
		this->mBounds.push_back({2.7,7.5});
	}
};
template <class T>
class Problem04Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Problem04<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Problem04<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Problem04<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Problem04<IntervalDer<T>>(), box);
	}

	Problem04Benchmark() : Benchmark<T>("Problem 04 function",-3.85045)	{
		this->mGlobMinX.push_back(2.86803);
		this->mBounds.push_back({1.9,3.9});
	}
};
template <class T>
class Problem05Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Problem05<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Problem05<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Problem05<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Problem05<IntervalDer<T>>(), box);
	}

	Problem05Benchmark() : Benchmark<T>("Problem 05 function",-1.48907)	{
		this->mGlobMinX.push_back(0.96609);
		this->mBounds.push_back({0,1.2});
	}
};
template <class T>
class Problem06Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Problem06<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Problem06<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Problem06<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Problem06<IntervalDer<T>>(), box);
	}

	Problem06Benchmark() : Benchmark<T>("Problem 06 function",-0.824239)	{
		this->mGlobMinX.push_back(0.67956);
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class QingBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Qing<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Qing<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Qing<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Qing<IntervalDer<T>>(), box);
	}

	QingBenchmark() : Benchmark<T>("Qing function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(1.41421);
		this->mBounds.push_back({-500,500});
		this->mBounds.push_back({-500,500});
	}
};
template <class T>
class QuadraticBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Quadratic<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Quadratic<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Quadratic<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Quadratic<IntervalDer<T>>(), box);
	}

	QuadraticBenchmark() : Benchmark<T>("Quadratic function",-3873.72)	{
		this->mGlobMinX.push_back(0.19388);
		this->mGlobMinX.push_back(0.48513);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class QuinticBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Quintic<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Quintic<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Quintic<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Quintic<IntervalDer<T>>(this->getDim()), box);
	}

	QuinticBenchmark(int size) : Benchmark<T>("Quintic function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(-1);
			this->mBounds.push_back({-10,10});
		}
	}
};
template <class T>
class RosenbrockBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Rosenbrock<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Rosenbrock<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Rosenbrock<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Rosenbrock<IntervalDer<T>>(this->getDim()), box);
	}

	RosenbrockBenchmark(int size) : Benchmark<T>("Rosenbrock function",0,false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(1);
			this->mBounds.push_back({-30,30});
		}
	}
};
template <class T>
class RosenbrockModifiedBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(RosenbrockModified<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(RosenbrockModified<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(RosenbrockModified<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(RosenbrockModified<IntervalDer<T>>(), box);
	}

	RosenbrockModifiedBenchmark() : Benchmark<T>("Rosenbrock Modified function",34.3712)	{
		this->mGlobMinX.push_back(-0.9);
		this->mGlobMinX.push_back(-0.95);
		this->mBounds.push_back({-2,2});
		this->mBounds.push_back({-2,2});
	}
};
template <class T>
class RotatedEllipseBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(RotatedEllipse<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(RotatedEllipse<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(RotatedEllipse<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(RotatedEllipse<IntervalDer<T>>(), box);
	}

	RotatedEllipseBenchmark() : Benchmark<T>("Rotated Ellipse function",0,false)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-500,500});
		this->mBounds.push_back({-500,500});
	}
};
template <class T>
class RotatedEllipse2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(RotatedEllipse2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(RotatedEllipse2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(RotatedEllipse2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(RotatedEllipse2<IntervalDer<T>>(), box);
	}

	RotatedEllipse2Benchmark() : Benchmark<T>("Rotated Ellipse 2 function",0,false)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-500,500});
		this->mBounds.push_back({-500,500});
	}
};
template <class T>
class Scahffer1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Scahffer1<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Scahffer1<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Scahffer1<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Scahffer1<IntervalDer<T>>(), box);
	}

	Scahffer1Benchmark() : Benchmark<T>("Scahffer 1 function",0,false)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class Scahffer3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Scahffer3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Scahffer3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Scahffer3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Scahffer3<IntervalDer<T>>(), box);
	}

	Scahffer3Benchmark() : Benchmark<T>("Scahffer 3 function",0.00156685,false)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(1.25311);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class Scahffer4Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Scahffer4<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Scahffer4<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Scahffer4<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Scahffer4<IntervalDer<T>>(), box);
	}

	Scahffer4Benchmark() : Benchmark<T>("Scahffer 4 function",0.292579,false)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(1.25311);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class Scahffer2_6Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Scahffer2_6<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Scahffer2_6<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Scahffer2_6<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Scahffer2_6<IntervalDer<T>>(), box);
	}

	Scahffer2_6Benchmark() : Benchmark<T>("Scahffer 2.6 function",0)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(3);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class SchafferF6Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(SchafferF6<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(SchafferF6<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(SchafferF6<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(SchafferF6<IntervalDer<T>>(this->getDim()), box);
	}

	SchafferF6Benchmark(int size) : Benchmark<T>("Schaffer F6 function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-100,100});
		}
	}
};
template <class T>
class SchmidtVettersBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(SchmidtVetters<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(SchmidtVetters<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(SchmidtVetters<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(SchmidtVetters<IntervalDer<T>>(), box);
	}

	SchmidtVettersBenchmark() : Benchmark<T>("Schmidt Vetters function",0.194018)	{
		this->mGlobMinX.push_back(7.05852);
		this->mGlobMinX.push_back(9.9855);
		this->mGlobMinX.push_back(3.10432);
		this->mBounds.push_back({0.1,10});
		this->mBounds.push_back({0.1,10});
		this->mBounds.push_back({0.1,10});
	}
};
template <class T>
class SchumerSteiglitzBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(SchumerSteiglitz<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(SchumerSteiglitz<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(SchumerSteiglitz<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(SchumerSteiglitz<IntervalDer<T>>(this->getDim()), box);
	}

	SchumerSteiglitzBenchmark(int size) : Benchmark<T>("Schumer Steiglitz function",0,false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-100,100});
		}
	}
};
template <class T>
class SchwefelBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Schwefel<T>(this->getDim(), 1.0), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Schwefel<Interval<T>>(this->getDim(), 1.0), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Schwefel<ValDer<T>>(this->getDim(), 1.0), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Schwefel<IntervalDer<T>>(this->getDim(), 1.0), box);
	}

	SchwefelBenchmark(int size) : Benchmark<T>("Schwefel function",0,false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-100,100});
		}
	}
};
template <class T>
class Schwefel1_2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Schwefel1_2<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Schwefel1_2<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Schwefel1_2<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Schwefel1_2<IntervalDer<T>>(this->getDim()), box);
	}

	Schwefel1_2Benchmark(int size) : Benchmark<T>("Schwefel 1.2 function",0,false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-100,100});
		}
	}
};
template <class T>
class Schwefel2_20Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Schwefel2_20<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Schwefel2_20<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Schwefel2_20<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Schwefel2_20<IntervalDer<T>>(this->getDim()), box);
	}

	Schwefel2_20Benchmark(int size) : Benchmark<T>("Schwefel 2.20 function",0,false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-100,100});
		}
	}
};
template <class T>
class Schwefel2_22Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Schwefel2_22<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Schwefel2_22<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Schwefel2_22<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Schwefel2_22<IntervalDer<T>>(this->getDim()), box);
	}

	Schwefel2_22Benchmark(int size) : Benchmark<T>("Schwefel 2.22 function",0,false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-100,100});
		}
	}
};
template <class T>
class Schwefel2_23Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Schwefel2_23<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Schwefel2_23<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Schwefel2_23<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Schwefel2_23<IntervalDer<T>>(this->getDim()), box);
	}

	Schwefel2_23Benchmark(int size) : Benchmark<T>("Schwefel 2.23 function",0,false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-10,10});
		}
	}
};
template <class T>
class Schwefel2_26Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Schwefel2_26<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Schwefel2_26<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Schwefel2_26<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Schwefel2_26<IntervalDer<T>>(), box);
	}

	Schwefel2_26Benchmark() : Benchmark<T>("Schwefel 2.26 function",-837.966)	{
		this->mGlobMinX.push_back(420.969);
		this->mGlobMinX.push_back(420.969);
		this->mBounds.push_back({-500,500});
		this->mBounds.push_back({-500,500});
	}
};
template <class T>
class Schwefel2_36Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Schwefel2_36<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Schwefel2_36<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Schwefel2_36<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Schwefel2_36<IntervalDer<T>>(), box);
	}

	Schwefel2_36Benchmark() : Benchmark<T>("Schwefel 2.36 function",-3456)	{
		this->mGlobMinX.push_back(12);
		this->mGlobMinX.push_back(12);
		this->mBounds.push_back({0,500});
		this->mBounds.push_back({0,500});
	}
};
template <class T>
class Schwefel2_4Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Schwefel2_4<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Schwefel2_4<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Schwefel2_4<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Schwefel2_4<IntervalDer<T>>(this->getDim()), box);
	}

	Schwefel2_4Benchmark(int size) : Benchmark<T>("Schwefel 2.4 function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(1);
			this->mBounds.push_back({0,10});
		}
	}
};
template <class T>
class Shekel10Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Shekel10<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Shekel10<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Shekel10<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Shekel10<IntervalDer<T>>(), box);
	}

	Shekel10Benchmark() : Benchmark<T>("Shekel 10 function",-10.5341)	{
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(4);
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
	}
};
template <class T>
class Shekel5Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Shekel5<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Shekel5<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Shekel5<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Shekel5<IntervalDer<T>>(), box);
	}

	Shekel5Benchmark() : Benchmark<T>("Shekel 5 function",-10.1532)	{
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(4);
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
	}
};
template <class T>
class Shekel7Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Shekel7<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Shekel7<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Shekel7<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Shekel7<IntervalDer<T>>(), box);
	}

	Shekel7Benchmark() : Benchmark<T>("Shekel 7 function",-10.4028)	{
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(4);
		this->mGlobMinX.push_back(4);
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
		this->mBounds.push_back({0,10});
	}
};
template <class T>
class ShubertBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Shubert<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Shubert<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Shubert<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Shubert<IntervalDer<T>>(), box);
	}

	ShubertBenchmark() : Benchmark<T>("Shubert function",-24.0625)	{
		this->mGlobMinX.push_back(-6.77458);
		this->mGlobMinX.push_back(-6.77458);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Shubert2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Shubert2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Shubert2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Shubert2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Shubert2<IntervalDer<T>>(), box);
	}

	Shubert2Benchmark() : Benchmark<T>("Shubert 2 function",-186.731)	{
		this->mGlobMinX.push_back(-7.0835);
		this->mGlobMinX.push_back(4.858);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Shubert3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Shubert3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Shubert3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Shubert3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Shubert3<IntervalDer<T>>(), box);
	}

	Shubert3Benchmark() : Benchmark<T>("Shubert 3 function",-24.0625)	{
		this->mGlobMinX.push_back(-6.77458);
		this->mGlobMinX.push_back(-6.77458);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class SolomonBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Solomon<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Solomon<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Solomon<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Solomon<IntervalDer<T>>(), box);
	}

	SolomonBenchmark() : Benchmark<T>("Solomon function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class SphereBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Sphere<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Sphere<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Sphere<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Sphere<IntervalDer<T>>(this->getDim()), box);
	}

	SphereBenchmark(int size) : Benchmark<T>("Sphere function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({0,10});
		}
	}
};
template <class T>
class StrechedVSineWaveBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(StrechedVSineWave<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(StrechedVSineWave<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(StrechedVSineWave<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(StrechedVSineWave<IntervalDer<T>>(this->getDim()), box);
	}

	StrechedVSineWaveBenchmark(int size) : Benchmark<T>("Streched V Sine Wave function",0,false)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-10,10});
		}
	}
};
template <class T>
class StyblinskiTangBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(StyblinskiTang<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(StyblinskiTang<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(StyblinskiTang<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(StyblinskiTang<IntervalDer<T>>(), box);
	}

	StyblinskiTangBenchmark() : Benchmark<T>("Styblinski-Tang function",-78.332)	{
		this->mGlobMinX.push_back(-2.90353);
		this->mGlobMinX.push_back(-2.90353);
		this->mBounds.push_back({-5,5});
		this->mBounds.push_back({-5,5});
	}
};
template <class T>
class SumSquaresBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(SumSquares<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(SumSquares<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(SumSquares<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(SumSquares<IntervalDer<T>>(this->getDim()), box);
	}

	SumSquaresBenchmark(int size) : Benchmark<T>("Sum Squares function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-10,10});
		}
	}
};
template <class T>
class Table1HolderTable1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Table1HolderTable1<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Table1HolderTable1<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Table1HolderTable1<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Table1HolderTable1<IntervalDer<T>>(), box);
	}

	Table1HolderTable1Benchmark() : Benchmark<T>("Table 1 / Holder Table 1 function",-26.9203)	{
		this->mGlobMinX.push_back(9.64617);
		this->mGlobMinX.push_back(9.64617);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Table2HolderTable2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Table2HolderTable2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Table2HolderTable2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Table2HolderTable2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Table2HolderTable2<IntervalDer<T>>(), box);
	}

	Table2HolderTable2Benchmark() : Benchmark<T>("Table 2 / Holder Table 2 function",-19.2085)	{
		this->mGlobMinX.push_back(8.05502);
		this->mGlobMinX.push_back(9.66459);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Table3CarromBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Table3Carrom<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Table3Carrom<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Table3Carrom<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Table3Carrom<IntervalDer<T>>(), box);
	}

	Table3CarromBenchmark() : Benchmark<T>("Table 3 Carrom Table function",-24.1568)	{
		this->mGlobMinX.push_back(9.64616);
		this->mGlobMinX.push_back(9.64613);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class TesttubeHolderBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(TesttubeHolder<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(TesttubeHolder<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(TesttubeHolder<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(TesttubeHolder<IntervalDer<T>>(), box);
	}

	TesttubeHolderBenchmark() : Benchmark<T>("Testtube Holder function",-10.8723)	{
		this->mGlobMinX.push_back(-1.5708);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class TrecanniBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Trecanni<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Trecanni<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Trecanni<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Trecanni<IntervalDer<T>>(), box);
	}

	TrecanniBenchmark() : Benchmark<T>("Trecanni function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-5,5});
		this->mBounds.push_back({-5,5});
	}
};
template <class T>
class TrefethenBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Trefethen<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Trefethen<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Trefethen<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Trefethen<IntervalDer<T>>(), box);
	}

	TrefethenBenchmark() : Benchmark<T>("Trefethen function",-3.30687)	{
		this->mGlobMinX.push_back(-0.024403);
		this->mGlobMinX.push_back(0.210612);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class Trid10Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Trid10<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Trid10<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Trid10<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Trid10<IntervalDer<T>>(), box);
	}

	Trid10Benchmark() : Benchmark<T>("Trid 10 function",-210)	{
		this->mGlobMinX.push_back(10);
		this->mGlobMinX.push_back(18);
		this->mGlobMinX.push_back(24);
		this->mGlobMinX.push_back(28);
		this->mGlobMinX.push_back(30);
		this->mGlobMinX.push_back(30);
		this->mGlobMinX.push_back(28);
		this->mGlobMinX.push_back(24);
		this->mGlobMinX.push_back(18);
		this->mGlobMinX.push_back(10);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class Trid6Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Trid6<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Trid6<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Trid6<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Trid6<IntervalDer<T>>(), box);
	}

	Trid6Benchmark() : Benchmark<T>("Trid 6 function",-50)	{
		this->mGlobMinX.push_back(6);
		this->mGlobMinX.push_back(10);
		this->mGlobMinX.push_back(12);
		this->mGlobMinX.push_back(12);
		this->mGlobMinX.push_back(10);
		this->mGlobMinX.push_back(6);
		this->mBounds.push_back({-36,36});
		this->mBounds.push_back({-36,36});
		this->mBounds.push_back({-36,36});
		this->mBounds.push_back({-36,36});
		this->mBounds.push_back({-36,36});
		this->mBounds.push_back({-36,36});
	}
};
template <class T>
class Trigonometric1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Trigonometric1<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Trigonometric1<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Trigonometric1<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Trigonometric1<IntervalDer<T>>(this->getDim()), box);
	}

	Trigonometric1Benchmark(int size) : Benchmark<T>("Trigonometric 1 function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({0,3.14159});
		}
	}
};
template <class T>
class Trigonometric2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Trigonometric2<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Trigonometric2<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Trigonometric2<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Trigonometric2<IntervalDer<T>>(this->getDim()), box);
	}

	Trigonometric2Benchmark(int size) : Benchmark<T>("Trigonometric 2 function",1)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0.9);
			this->mBounds.push_back({-500,500});
		}
	}
};
template <class T>
class TripodBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Tripod<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Tripod<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Tripod<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Tripod<IntervalDer<T>>(), box);
	}

	TripodBenchmark() : Benchmark<T>("Tripod function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(-50);
		this->mBounds.push_back({-100,100});
		this->mBounds.push_back({-100,100});
	}
};
template <class T>
class Ursem1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Ursem1<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Ursem1<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Ursem1<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Ursem1<IntervalDer<T>>(), box);
	}

	Ursem1Benchmark() : Benchmark<T>("Ursem 1 function",-4.81681)	{
		this->mGlobMinX.push_back(1.7049);
		this->mGlobMinX.push_back(0.0146484);
		this->mBounds.push_back({-2.5,3});
		this->mBounds.push_back({-2,2});
	}
};
template <class T>
class Ursem3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Ursem3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Ursem3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Ursem3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Ursem3<IntervalDer<T>>(), box);
	}

	Ursem3Benchmark() : Benchmark<T>("Ursem 3 function",-2.5)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-2,2});
		this->mBounds.push_back({-1.5,1.5});
	}
};
template <class T>
class Ursem4Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Ursem4<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Ursem4<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Ursem4<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Ursem4<IntervalDer<T>>(), box);
	}

	Ursem4Benchmark() : Benchmark<T>("Ursem 4 function",-1.5)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-2,2});
		this->mBounds.push_back({-2,2});
	}
};
template <class T>
class UrsemWavesBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(UrsemWaves<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(UrsemWaves<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(UrsemWaves<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(UrsemWaves<IntervalDer<T>>(), box);
	}

	UrsemWavesBenchmark() : Benchmark<T>("Ursem Waves function",-7.307)	{
		this->mGlobMinX.push_back(-0.605689);
		this->mGlobMinX.push_back(-1.17756);
		this->mBounds.push_back({-0.9,1.2});
		this->mBounds.push_back({-1.2,1.2});
	}
};
template <class T>
class VenterSobiezcczanskiSobieskiBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(VenterSobiezcczanskiSobieski<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(VenterSobiezcczanskiSobieski<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(VenterSobiezcczanskiSobieski<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(VenterSobiezcczanskiSobieski<IntervalDer<T>>(), box);
	}

	VenterSobiezcczanskiSobieskiBenchmark() : Benchmark<T>("Venter Sobiezcczanski-Sobieski function",-400)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-50,50});
		this->mBounds.push_back({-50,50});
	}
};
template <class T>
class WWavyBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(WWavy<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(WWavy<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(WWavy<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(WWavy<IntervalDer<T>>(this->getDim()), box);
	}

	WWavyBenchmark(int size) : Benchmark<T>("W Wavy function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-3.14159,3.14159});
		}
	}
};
template <class T>
class WatsonBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Watson<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Watson<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Watson<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Watson<IntervalDer<T>>(), box);
	}

	WatsonBenchmark() : Benchmark<T>("Watson function",0.002288)	{
		this->mGlobMinX.push_back(-0.0158);
		this->mGlobMinX.push_back(1.012);
		this->mGlobMinX.push_back(-0.2329);
		this->mGlobMinX.push_back(1.26);
		this->mGlobMinX.push_back(-1.513);
		this->mGlobMinX.push_back(0.9928);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class WayburnSeader1Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(WayburnSeader1<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(WayburnSeader1<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(WayburnSeader1<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(WayburnSeader1<IntervalDer<T>>(), box);
	}

	WayburnSeader1Benchmark() : Benchmark<T>("Wayburn Seader 1 function",0,false)	{
		this->mGlobMinX.push_back(1);
		this->mGlobMinX.push_back(2);
		this->mBounds.push_back({-500,500});
		this->mBounds.push_back({-500,500});
	}
};
template <class T>
class WayburnSeader2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(WayburnSeader2<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(WayburnSeader2<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(WayburnSeader2<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(WayburnSeader2<IntervalDer<T>>(), box);
	}

	WayburnSeader2Benchmark() : Benchmark<T>("Wayburn Seader 2 function",0,false)	{
		this->mGlobMinX.push_back(0.2);
		this->mGlobMinX.push_back(1);
		this->mBounds.push_back({-500,500});
		this->mBounds.push_back({-500,500});
	}
};
template <class T>
class WayburnSeader3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(WayburnSeader3<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(WayburnSeader3<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(WayburnSeader3<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(WayburnSeader3<IntervalDer<T>>(), box);
	}

	WayburnSeader3Benchmark() : Benchmark<T>("Wayburn Seader 3 function",19.1059,false)	{
		this->mGlobMinX.push_back(5.1469);
		this->mGlobMinX.push_back(6.83959);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};
template <class T>
class WeierstrassBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Weierstrass<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Weierstrass<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Weierstrass<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Weierstrass<IntervalDer<T>>(this->getDim()), box);
	}

	WeierstrassBenchmark(int size) : Benchmark<T>("Weierstrass function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-0.5,0.5});
		}
	}
};
template <class T>
class WhitleyBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Whitley<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Whitley<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Whitley<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Whitley<IntervalDer<T>>(this->getDim()), box);
	}

	WhitleyBenchmark(int size) : Benchmark<T>("Whitley function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(1);
			this->mBounds.push_back({-10,10});
		}
	}
};
template <class T>
class WolfeBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Wolfe<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Wolfe<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Wolfe<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Wolfe<IntervalDer<T>>(), box);
	}

	WolfeBenchmark() : Benchmark<T>("Wolfe function",0)	{
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({0,2});
		this->mBounds.push_back({0,2});
		this->mBounds.push_back({0,2});
	}
};
template <class T>
class XinSheYang2Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(XinSheYang2<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(XinSheYang2<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(XinSheYang2<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(XinSheYang2<IntervalDer<T>>(this->getDim()), box);
	}

	XinSheYang2Benchmark(int size) : Benchmark<T>("Xin-She Yang 2 function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-6.28319,6.28319});
		}
	}
};
template <class T>
class XinSheYang3Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(XinSheYang3<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(XinSheYang3<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(XinSheYang3<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(XinSheYang3<IntervalDer<T>>(this->getDim()), box);
	}

	XinSheYang3Benchmark(int size) : Benchmark<T>("Xin-She Yang 3 function",-1)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-20,20});
		}
	}
};
template <class T>
class XinSheYang4Benchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(XinSheYang4<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(XinSheYang4<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(XinSheYang4<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(XinSheYang4<IntervalDer<T>>(this->getDim()), box);
	}

	XinSheYang4Benchmark(int size) : Benchmark<T>("Xin-She Yang 4function",-1)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-10,10});
		}
	}
};
template <class T>
class ZakharovBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Zakharov<T>(this->getDim()), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Zakharov<Interval<T>>(this->getDim()), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Zakharov<ValDer<T>>(this->getDim()), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Zakharov<IntervalDer<T>>(this->getDim()), box);
	}

	ZakharovBenchmark(int size) : Benchmark<T>("Zakharov function",0)	{
		for(int i=0; i < size; ++i) {
			this->mGlobMinX.push_back(0);
			this->mBounds.push_back({-5,10});
		}
	}
};
template <class T>
class ZettlBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Zettl<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Zettl<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Zettl<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Zettl<IntervalDer<T>>(), box);
	}

	ZettlBenchmark() : Benchmark<T>("Zettl function",-0.003791)	{
		this->mGlobMinX.push_back(-0.0299);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-5,10});
		this->mBounds.push_back({-5,10});
	}
};
template <class T>
class ZirilliBenchmark : public Benchmark<T>
{
public:
	virtual T calcFunc(const std::vector<T>& point) const {
		return ::calcFunc(Zirilli<T>(), point);
	}

	virtual Interval<T> calcInterval(const std::vector<Interval<T>>& box) const {
		return ::calcInterval(Zirilli<Interval<T>>(), box);
	}

	virtual ValDer<T> calcGrad(const std::vector<T>& point) const {
		return ::calcGrad(Zirilli<ValDer<T>>(), point);
	}

	virtual IntervalDer<T> calcIntervalGrad(const std::vector<Interval<T>>& box) const {
		return ::calcIntervalGrad(Zirilli<IntervalDer<T>>(), box);
	}

	ZirilliBenchmark() : Benchmark<T>("Zirilli or Aluffi-Pentini’s function",-0.3523,false)	{
		this->mGlobMinX.push_back(-1.0465);
		this->mGlobMinX.push_back(0);
		this->mBounds.push_back({-10,10});
		this->mBounds.push_back({-10,10});
	}
};



template <class T>
using PtrBench = std::shared_ptr<Benchmark<T>>;

template <class T>
class Benchmarks
{
private:
	std::vector<PtrBench<T>> vbm;
public:
	Benchmarks() { fill(); }
	void fill()
	{
		clear();
		add(std::make_shared<Ackley1Benchmark<double>>(3));
		add(std::make_shared<Ackley2Benchmark<double>>(4));
		add(std::make_shared<Ackley3Benchmark<double>>());
		add(std::make_shared<AdjimanBenchmark<double>>());
		add(std::make_shared<Alpine1Benchmark<double>>(3));
		add(std::make_shared<Alpine2Benchmark<double>>());
		add(std::make_shared<BradBenchmark<double>>());
		add(std::make_shared<BartelsConnBenchmark<double>>());
		add(std::make_shared<BealeBenchmark<double>>());
		add(std::make_shared<BiggsEXP2Benchmark<double>>());
		add(std::make_shared<BiggsEXP3Benchmark<double>>());
		add(std::make_shared<BiggsEXP4Benchmark<double>>());
		add(std::make_shared<BiggsEXP5Benchmark<double>>());
		add(std::make_shared<BiggsEXP6Benchmark<double>>());
		add(std::make_shared<BirdBenchmark<double>>());
		add(std::make_shared<Bohachevsky1Benchmark<double>>());
		add(std::make_shared<Bohachevsky2Benchmark<double>>());
		add(std::make_shared<Bohachevsky3Benchmark<double>>());
		add(std::make_shared<BoothBenchmark<double>>());
		add(std::make_shared<BoxBettsQuadraticSumBenchmark<double>>());
		add(std::make_shared<BraninRCOSBenchmark<double>>());
		add(std::make_shared<BraninRCOS2Benchmark<double>>());
		add(std::make_shared<BrentBenchmark<double>>());
		add(std::make_shared<BrownBenchmark<double>>(3));
		add(std::make_shared<Bukin2Benchmark<double>>());
		add(std::make_shared<Bukin4Benchmark<double>>());
		add(std::make_shared<Bukin6Benchmark<double>>());
		add(std::make_shared<CamelSixHumpBenchmark<double>>());
		add(std::make_shared<CamelThreeHumpBenchmark<double>>());
		add(std::make_shared<ChichinadzeBenchmark<double>>());
		add(std::make_shared<ChungReynoldsBenchmark<double>>(3));
		add(std::make_shared<ColvilleBenchmark<double>>());
		add(std::make_shared<ComplexBenchmark<double>>());
		add(std::make_shared<CosineMixtureBenchmark<double>>());
		add(std::make_shared<CrossInTrayBenchmark<double>>());
		add(std::make_shared<CrossLegBenchmark<double>>());
		add(std::make_shared<CubeBenchmark<double>>());
		add(std::make_shared<DavisBenchmark<double>>());
		add(std::make_shared<Deb1Benchmark<double>>(3));
		add(std::make_shared<DeckkersAartsBenchmark<double>>());
		add(std::make_shared<DixonPriceBenchmark<double>>());
		add(std::make_shared<DolanBenchmark<double>>());
		add(std::make_shared<DropWaveBenchmark<double>>());
		add(std::make_shared<EasomBenchmark<double>>());
		add(std::make_shared<EggCrateBenchmark<double>>());
		add(std::make_shared<EggHolderBenchmark<double>>());
		add(std::make_shared<ElAttarVidyasagarDuttBenchmark<double>>());
		add(std::make_shared<EngvallBenchmark<double>>());
		add(std::make_shared<Exp2Benchmark<double>>());
		add(std::make_shared<ExponentialBenchmark<double>>(3));
		add(std::make_shared<FreudensteinRothBenchmark<double>>());
		add(std::make_shared<GoldsteinPriceBenchmark<double>>());
		add(std::make_shared<GramacyLee2Benchmark<double>>());
		add(std::make_shared<GramacyLee3Benchmark<double>>());
		add(std::make_shared<GriewankBenchmark<double>>(3));
		add(std::make_shared<HansenBenchmark<double>>());
		add(std::make_shared<Hartman3Benchmark<double>>());
		add(std::make_shared<Hartman6Benchmark<double>>());
		add(std::make_shared<HelicalValleyBenchmark<double>>());
		add(std::make_shared<HimmelblauBenchmark<double>>());
		add(std::make_shared<HosakiBenchmark<double>>());
		add(std::make_shared<JennrichSampsonBenchmark<double>>());
		add(std::make_shared<KeaneBenchmark<double>>());
		add(std::make_shared<Langerman5Benchmark<double>>());
		add(std::make_shared<LeonBenchmark<double>>());
		add(std::make_shared<MatyasBenchmark<double>>());
		add(std::make_shared<McCormickBenchmark<double>>());
		add(std::make_shared<MieleCantrellBenchmark<double>>());
		add(std::make_shared<Mishra3Benchmark<double>>());
		add(std::make_shared<Mishra4Benchmark<double>>());
		add(std::make_shared<Mishra5Benchmark<double>>());
		add(std::make_shared<Mishra6Benchmark<double>>());
		add(std::make_shared<Mishra7Benchmark<double>>());
		add(std::make_shared<Mishra8Benchmark<double>>());
		add(std::make_shared<Mishra9Benchmark<double>>());
		add(std::make_shared<ParsopoulosBenchmark<double>>());
		add(std::make_shared<PeriodicBenchmark<double>>());
		add(std::make_shared<PinterBenchmark<double>>(3));
		add(std::make_shared<PowellSingular2Benchmark<double>>(8));
		add(std::make_shared<PowellSingular2Benchmark<double>>(3));
		add(std::make_shared<Price1Benchmark<double>>());
		add(std::make_shared<Price2Benchmark<double>>());
		add(std::make_shared<Price3Benchmark<double>>());
		add(std::make_shared<Price4Benchmark<double>>());
		add(std::make_shared<Problem02Benchmark<double>>());
		add(std::make_shared<Problem04Benchmark<double>>());
		add(std::make_shared<Problem05Benchmark<double>>());
		add(std::make_shared<Problem06Benchmark<double>>());
		add(std::make_shared<QingBenchmark<double>>());
		add(std::make_shared<QuadraticBenchmark<double>>());
		add(std::make_shared<QuinticBenchmark<double>>(3));
		add(std::make_shared<RosenbrockBenchmark<double>>(3));
		add(std::make_shared<RosenbrockModifiedBenchmark<double>>());
		add(std::make_shared<RotatedEllipseBenchmark<double>>());
		add(std::make_shared<RotatedEllipse2Benchmark<double>>());
		add(std::make_shared<Scahffer1Benchmark<double>>());
		add(std::make_shared<Scahffer3Benchmark<double>>());
		add(std::make_shared<Scahffer4Benchmark<double>>());
		add(std::make_shared<Scahffer2_6Benchmark<double>>());
		add(std::make_shared<SchafferF6Benchmark<double>>(3));
		add(std::make_shared<SchmidtVettersBenchmark<double>>());
		add(std::make_shared<SchumerSteiglitzBenchmark<double>>(3));
		add(std::make_shared<SchwefelBenchmark<double>>(3));
		add(std::make_shared<Schwefel1_2Benchmark<double>>(3));
		add(std::make_shared<Schwefel2_20Benchmark<double>>(3));
        add(std::make_shared<Schwefel2_22Benchmark<double>>(3));
		add(std::make_shared<Schwefel2_23Benchmark<double>>(3));
		add(std::make_shared<Schwefel2_26Benchmark<double>>());
		add(std::make_shared<Schwefel2_36Benchmark<double>>());
		add(std::make_shared<Schwefel2_4Benchmark<double>>(3));
		add(std::make_shared<Shekel10Benchmark<double>>());
		add(std::make_shared<Shekel5Benchmark<double>>());
		add(std::make_shared<Shekel7Benchmark<double>>());
		add(std::make_shared<ShubertBenchmark<double>>());
		add(std::make_shared<Shubert2Benchmark<double>>());
		add(std::make_shared<Shubert3Benchmark<double>>());
		add(std::make_shared<SolomonBenchmark<double>>());
		add(std::make_shared<SphereBenchmark<double>>(3));
		add(std::make_shared<StrechedVSineWaveBenchmark<double>>(3));
		add(std::make_shared<StyblinskiTangBenchmark<double>>());
		add(std::make_shared<SumSquaresBenchmark<double>>(3));
		add(std::make_shared<Table1HolderTable1Benchmark<double>>());
		add(std::make_shared<Table2HolderTable2Benchmark<double>>());
		add(std::make_shared<Table3CarromBenchmark<double>>());
		add(std::make_shared<TesttubeHolderBenchmark<double>>());
		add(std::make_shared<TrecanniBenchmark<double>>());
		add(std::make_shared<TrefethenBenchmark<double>>());
		add(std::make_shared<Trid10Benchmark<double>>());
		add(std::make_shared<Trid6Benchmark<double>>());
		add(std::make_shared<Trigonometric1Benchmark<double>>(3));
		add(std::make_shared<Trigonometric2Benchmark<double>>(3));
		add(std::make_shared<TripodBenchmark<double>>());
		add(std::make_shared<Ursem1Benchmark<double>>());
		add(std::make_shared<Ursem3Benchmark<double>>());
		add(std::make_shared<Ursem4Benchmark<double>>());
		add(std::make_shared<UrsemWavesBenchmark<double>>());
		add(std::make_shared<VenterSobiezcczanskiSobieskiBenchmark<double>>());
		add(std::make_shared<WWavyBenchmark<double>>(3));
		add(std::make_shared<WayburnSeader1Benchmark<double>>());
		add(std::make_shared<WayburnSeader2Benchmark<double>>());
		add(std::make_shared<WayburnSeader3Benchmark<double>>());
		add(std::make_shared<WeierstrassBenchmark<double>>(3));
		add(std::make_shared<WhitleyBenchmark<double>>(3));
		add(std::make_shared<WolfeBenchmark<double>>());
		add(std::make_shared<XinSheYang2Benchmark<double>>(3));
        add(std::make_shared<XinSheYang3Benchmark<double>>(3));
		add(std::make_shared<XinSheYang4Benchmark<double>>(3));
		add(std::make_shared<ZakharovBenchmark<double>>(3));
		add(std::make_shared<ZettlBenchmark<double>>());
		add(std::make_shared<ZirilliBenchmark<double>>());
	}
	void clear()
	{
		vbm.clear();
	}
	void add(const PtrBench<T> &ptrBenchmark)
	{
		vbm.push_back(ptrBenchmark);
	}
	
	void print() const
	{
		for(PtrBench<T> ptrBench : *this)
			std::cout << *ptrBench;
	}	

	typename std::vector<PtrBench<T>>::iterator begin()
	{
        	return vbm.begin();
    	}
    	typename std::vector<PtrBench<T>>::iterator end()
	{
        	return vbm.end();
    	}
    	typename std::vector<PtrBench<T>>::const_iterator begin() const 
	{
        	return vbm.begin();
    	}
    	typename std::vector<PtrBench<T>>::const_iterator end() const 
	{
        	return vbm.end();
    	}
};




#endif


