#ifndef ALGORITHM__HPP
#define ALGORITHM__HPP

#include <iostream> 
#include <memory>
#include <cmath>
#include <vector>

#include "interval/interval_air.hpp"
#include "interval/enums.h"


using namespace snowgoose::interval;

namespace snowgoose {
namespace expression {

	template<class T> class Algorithm;
	template<typename T> using ptrAlg = std::shared_ptr<Algorithm<T>>;
	template<typename T> using vPtrAlg = std::vector<ptrAlg<T>>;
	/**
	* Base class for algorithms
	*/
	template<class T>
	class Algorithm
	{
	public:
		virtual T Plus(const T& left, const T& right) const = 0;
		virtual T Minus(const T& left, const T& right) const = 0;
		virtual T Mul(const T& left, const T& right) const = 0;
		virtual T Div(const T& left, const T& right) const = 0;
		virtual T Sin(const T& t) const = 0;
		virtual T Cos(const T& t) const = 0;
		virtual T Tan(const T& t) const = 0;
		virtual T Ctg(const T& t) const = 0;
		virtual T ArcCos(const T& t) const = 0;
		virtual T ArcSin(const T& t) const = 0;
		virtual T ArcTan(const T& t) const = 0;
		virtual T ArcCtg(const T& t) const = 0;
		virtual T Exp(const T& t) const = 0;
		virtual T Sqrt(const T& t) const = 0;
		virtual T Sqr(const T& t) const = 0;
		virtual T Pow(const T& base, int exp) const = 0;
		virtual T Pow(const T& base, bool isBaseVar, const T& exp, bool isExpVar) const = 0;
        	virtual T PowDouble(const T& base, double exp) const = 0;
		virtual T Abs(const T& t) const = 0;		
		virtual T Ln(const T& t) const = 0;
		virtual T Log(const T& t, double base) const = 0;
		virtual T Min(const T& left, const T& right) const = 0;
		virtual T Max(const T& left, const T& right) const = 0;
		virtual T IfTrue(IntervalBool ib, const T& left, const T& right) const = 0;
		virtual IntervalBool Condition(Conditions condition, const T& left, const T& right) const = 0;
		virtual T CreateVar(int index) const = 0;
		virtual T CreateConst(double cnst) const = 0;
		virtual vPtrAlg<T> GetNewAlgorithm(Conditions condition, int index, double cnst) const { return vPtrAlg<T>(); }//empty vector
	};
	/**
	* Algorithm to calculate value of a function
	*/
	template<class T=double>
	class FuncAlg : public Algorithm<T>
	{
    private:
        std::vector<T> m_v;
	public:
        FuncAlg(const std::vector<T> &v) : m_v(v) {}
		T Plus(const T& left, const T& right) const { return left + right; }
		T Minus(const T& left, const T& right) const { return left - right; }
		T Mul(const T& left, const T& right) const { return left * right; }
		T Div(const T& left, const T& right) const 
		{ 
			if (right == 0.0)
				throw std::invalid_argument("Invalid operation. Divide by 0.0.");
			return left / right; 
		}
		T Sin(const T& t) const { return std::sin(t); };
		T Cos(const T& t) const { return std::cos(t); };
		T Tan(const T& t) const { return std::tan(t); };
		T Ctg(const T& t) const { return 1.0 / std::tan(t);  };

		T ArcCos(const T& t) const 
		{
			if(t < -1.0 || t > 1.0)
				throw std::invalid_argument("Invalid argument in arccos function. The argument is out of this interval [-1,1].");
			return std::acos(t); 
		};
		T ArcSin(const T& t) const
		{
			if (t < -1.0 || t > 1.0)
				throw std::invalid_argument("Invalid argument in arcsin function. The argument is out of this interval [-1,1].");
			return std::asin(t);
		};
		T ArcTan(const T& t) const { return std::atan(t); }
		T ArcCtg(const T& t) const { return M_PI_2 - std::atan(t); }
		T Exp(const T& t) const { return std::exp(t); };
		T Sqrt(const T& t) const 
		{
			if(t < 0.0)
				throw std::invalid_argument("The function sqrt is not defined for negative numbers");
			return std::sqrt(t); 
		};
		T Sqr(const T& t) const { return t*t; };
		T Pow(const T& base, int exp) const { return std::pow(base, exp); };
		T Pow(const T& base, bool isBaseVar, const T& exp, bool isExpVar) const
		{
			if(base < 0.0)
				throw std::invalid_argument("The function FuncAlg::pow is not define for negative base");
			return std::pow(base, exp);
		};
		T PowDouble(const T& base, double exp) const 
		{
		    return Pow(base, true, T(exp), false); 
		}
			T Abs(const T& t) const 
		{ 
		    return std::abs(t); 
		};
		T Ln(const T& t) const
		{
			if (t <= 0)
				throw std::invalid_argument("The function Ln is not define for negative numbers and 0.0");
			return std::log(t);
		};
		T Log(const T& t, double base) const
		{
			if (t <= 0)
				throw std::invalid_argument("The function Log is not define for negative numbers and 0.0");
			return std::log(t) / std::log(base);
		};
		T Min(const T& left, const T& right) const { return left < right ? left : right; }
		T Max(const T& left, const T& right) const { return left > right ? left : right; }
		T IfTrue(IntervalBool condition, const T& left, const T& right) const { return condition == IntervalBool::True ? left : right; }
		IntervalBool Condition(Conditions condition, const T& left, const T& right) const
		{
			switch (condition)
			{
				case Conditions::More :
					return left > right ? IntervalBool::True : IntervalBool::False;
					break;
				case Conditions::Less :
					return left < right ? IntervalBool::True : IntervalBool::False;
					break;
				case Conditions::LessEqual :
					return left <= right ? IntervalBool::True : IntervalBool::False;
					break;
				case Conditions::MoreEqual :
					return left >= right ? IntervalBool::True : IntervalBool::False;
					break;
			}
			throw std::invalid_argument("Invalid condition.");
		}
		T CreateVar(int index) const { return m_v[index]; }      
		T CreateConst(double cnst) const { return cnst; };
	};


	/**
	* Algorithm to calculate interval estimation of the function
	*/
	template<class T=double>
	class InterEvalAlg : public Algorithm<Interval<T>>
	{
    private:
        std::vector<Interval<T>> m_v;
	public:
        InterEvalAlg(const std::vector<Interval<T>> &v) : m_v(v) {}
		Interval<T> Plus(const Interval<T>& left, const Interval<T>& right) const { return left + right; }
		Interval<T> Minus(const Interval<T>& left, const Interval<T>& right) const { return left - right; }
		Interval<T> Mul(const Interval<T>& left, const Interval<T>& right) const { return left * right; }
		Interval<T> Div(const Interval<T>& left, const Interval<T>& right) const { return left / right; }
		Interval<T> Sin(const Interval<T>& t) const { return sin(t); };
		Interval<T> Cos(const Interval<T>& t) const { return cos(t); };
		Interval<T> Tan(const Interval<T>& t) const { return tg(t);  };
		Interval<T> Ctg(const Interval<T>& t) const { return ctg(t);  };
		Interval<T> ArcCos(const Interval<T>& t) const { return acos(t); };
		Interval<T> ArcSin(const Interval<T>& t) const { return asin(t); };
		Interval<T> ArcTan(const Interval<T>& t) const { return atg(t); };
		Interval<T> ArcCtg(const Interval<T>& t) const { return actg(t); };
		Interval<T> Exp(const Interval<T>& t) const { return exp(t); };
		Interval<T> Sqrt(const Interval<T>& t) const { return sqrt(t); };
		Interval<T> Sqr(const Interval<T>& t) const { return sqr(t); };
		Interval<T> Pow(const Interval<T>& base, int exp) const { return base ^ exp; };
		Interval<T> Pow(const Interval<T>& base, bool isBaseVar, const Interval<T>& exp, bool isExpVar) const { return base ^ exp; };
        	Interval<T> PowDouble(const Interval<T>& base, double exp) const { return base ^ Interval<T>(exp);}
		Interval<T> Abs(const Interval<T>& t) const { return abs(t); };
		Interval<T> Ln(const Interval<T>& t) const { return ln(t); };
		Interval<T> Log(const Interval<T>& t, double base) const { return log(t, base); };
		Interval<T> Min(const Interval<T>& left, const Interval<T>& right) const { return min(IL<T>({ left,right }));}
		Interval<T> Max(const Interval<T>& left, const Interval<T>& right) const { return max(IL<T>({ left,right }));};
		Interval<T> IfTrue(IntervalBool condition, const Interval<T>& left, const Interval<T>& right) const { return ifThen(condition, left, right); }
		IntervalBool Condition(Conditions condition, const Interval<T>& left, const Interval<T>& right) const
		{
			switch (condition)
			{
			case Conditions::More:
				return left > right;
				break;
			case Conditions::Less:
				return left < right;
				break;
			case Conditions::LessEqual:
				return left <= right;
				break;
			case Conditions::MoreEqual:
				return left >= right;
				break;
			}
			throw std::invalid_argument("Invalid condition.");
		}
		Interval<T> CreateVar(int index) const { return m_v[index]; }      
		Interval<T> CreateConst(double cnst) const { return cnst; };
		vPtrAlg<Interval<T>> GetNewAlgorithm(Conditions cond, int index, double cnst) const
		{ 
			std::vector<Interval<T>> leftVec(m_v);
			std::vector<Interval<T>> rightVec(m_v);
			if(cond == Conditions::More || cond == Conditions::MoreEqual) {	
				leftVec[index] = Interval<T>(cnst, m_v[index].rb());				
				rightVec[index] = Interval<T>(m_v[index].lb(), cnst);
			}
			else { // Less or LessEqual
				leftVec[index] = Interval<T>(m_v[index].lb(), cnst);				
				rightVec[index] = Interval<T>(cnst, m_v[index].rb());
			}
			return vPtrAlg<Interval<T>>({ ptrAlg<Interval<T>>(new InterEvalAlg<T>(leftVec)), ptrAlg<Interval<T>>(new InterEvalAlg<T>(rightVec))});
		}
            
	};
}
}

#endif /* ALGORITHM__HPP */
