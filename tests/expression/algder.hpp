#ifndef ALGDER__HPP
#define ALGDER__HPP

#include <iostream> 
#include <memory>
#include <cmath>
#include <vector>

#include "algorithm.hpp"
#include "derivatives/valder.hpp"
#include "derivatives/intervalder.hpp"

using namespace snowgoose::interval;
using namespace snowgoose::derivative;

namespace snowgoose {
namespace expression {

	template<class T>
	class ValDerAlg : public Algorithm<ValDer<T>>
	{
    	private:
        std::vector<T> m_v;      
	public:
        ValDerAlg(const std::vector<T> &v) : m_v(v) {}
		ValDer<T> Plus(const ValDer<T>& left, const ValDer<T>& right) const { return left + right; }
		ValDer<T> Minus(const ValDer<T>& left, const ValDer<T>& right) const { return left - right; }
		ValDer<T> Mul(const ValDer<T>& left, const ValDer<T>& right) const { return left * right; }
		ValDer<T> Div(const ValDer<T>& left, const ValDer<T>& right) const { return left / right; }
		ValDer<T> Sin(const ValDer<T>& t) const { return sin(t); }
		ValDer<T> Cos(const ValDer<T>& t) const { return cos(t); }
		ValDer<T> Tan(const ValDer<T>& t) const { return tg(t);  }
		ValDer<T> Ctg(const ValDer<T>& t) const { return ctg(t);  }
		ValDer<T> ArcCos(const ValDer<T>& t) const { return acos(t); }
		ValDer<T> ArcSin(const ValDer<T>& t) const { return asin(t); }
		ValDer<T> ArcTan(const ValDer<T>& t) const { return atg(t); }
		ValDer<T> ArcCtg(const ValDer<T>& t) const { return actg(t); }
		ValDer<T> Exp(const ValDer<T>& t) const { return exp(t); }
		ValDer<T> Sqrt(const ValDer<T>& t) const { return sqrt(t); }
		ValDer<T> Sqr(const ValDer<T>& t) const { return sqr(t); }
		ValDer<T> Pow(const ValDer<T>& base, int exp) const { return base ^ exp; }
		ValDer<T> Pow(const ValDer<T>& base, bool isBaseVar, const ValDer<T>& exp, bool isExpVar) const { 
		    if(isBaseVar && isExpVar) 
		        throw std::invalid_argument("Invalid argument in DerAlg::Pow. The arguments base and exp depend on var.");
		    if(!isBaseVar)
		        return pow(base.value(), exp);
		    else
		        return base ^ exp.value();
        	};
        	ValDer<T> PowDouble(const ValDer<T>& base, double exp) const { return base ^ (T)exp; }
		ValDer<T> Abs(const ValDer<T>& t) const { return abs(t); }
		ValDer<T> Ln(const ValDer<T>& t) const { return ln(t); }
		ValDer<T> Log(const ValDer<T>& t, double base) const { return log(t, base); }
		ValDer<T> Min(const ValDer<T>& left, const ValDer<T>& right) const { return min(left, right); }
		ValDer<T> Max(const ValDer<T>& left, const ValDer<T>& right) const { return max(left, right); }
		ValDer<T> IfTrue(IntervalBool condition, const ValDer<T>& left, const ValDer<T>& right) const { return ifThen(condition, left, right);  }
		ValDer<T> IfTrue(const ptrCNode<ValDer<T>>& conditionNode, const ptrNode<ValDer<T>>& left, const ptrNode<ValDer<T>>& right, MapIterator &map_iterator) const 
		{
			IntervalBool condition = conditionNode->calcCondition(*this, map_iterator);
			return ifThen(condition, left->calc(*this, map_iterator), right->calc(*this, map_iterator));
		}
		IntervalBool Condition(Conditions condition, const ValDer<T>& left, const ValDer<T>& right) const { return cond(condition, left, right); }
		ValDer<T> CreateVar(int index) const 
		{ 
		    std::vector<T> grad(m_v.size(), 0.0);
		    grad[index] = 1.0;
		    return ValDer<T>(m_v[index], Grad<T>(grad));
		}      
		ValDer<T> CreateConst(double cnst) const { return ValDer<T>(cnst, Grad<T>(m_v.size(), 0.0)); }
	};

	template<class T>
	class IntervalDerAlg : public Algorithm<IntervalDer<T>>
	{
    private:
        std::vector<Interval<T>> m_v;      
	public:
        IntervalDerAlg(const std::vector<Interval<T>> &v) : m_v(v) {}
		IntervalDer<T> Plus(const IntervalDer<T>& left, const IntervalDer<T>& right) const { return left + right; }
		IntervalDer<T> Minus(const IntervalDer<T>& left, const IntervalDer<T>& right) const { return left - right; }
		IntervalDer<T> Mul(const IntervalDer<T>& left, const IntervalDer<T>& right) const { return left * right; }
		IntervalDer<T> Div(const IntervalDer<T>& left, const IntervalDer<T>& right) const { return left / right; }
		IntervalDer<T> Sin(const IntervalDer<T>& t) const { return sin(t); }
		IntervalDer<T> Cos(const IntervalDer<T>& t) const { return cos(t); }
		IntervalDer<T> Tan(const IntervalDer<T>& t) const { return tg(t);  }
		IntervalDer<T> Ctg(const IntervalDer<T>& t) const { return ctg(t);  }
		IntervalDer<T> ArcCos(const IntervalDer<T>& t) const { return acos(t); }
		IntervalDer<T> ArcSin(const IntervalDer<T>& t) const { return asin(t); }
		IntervalDer<T> ArcTan(const IntervalDer<T>& t) const { return atg(t); }
		IntervalDer<T> ArcCtg(const IntervalDer<T>& t) const { return actg(t); }
		IntervalDer<T> Exp(const IntervalDer<T>& t) const { return exp(t); }
		IntervalDer<T> Sqrt(const IntervalDer<T>& t) const { return sqrt(t); }
		IntervalDer<T> Sqr(const IntervalDer<T>& t) const { return sqr(t); }
		IntervalDer<T> Pow(const IntervalDer<T>& base, int exp) const { return base ^ exp; }
		IntervalDer<T> Pow(const IntervalDer<T>& base, bool isBaseVar, const IntervalDer<T>& exp, bool isExpVar) const { 
		    if(isBaseVar && isExpVar) 
		        throw std::invalid_argument("Invalid argument in DerAlg::Pow. The arguments base and exp depend on var.");
		    if(!isBaseVar)
		        return pow(base.value(), exp);
		    else
		        return base ^ exp.value();
        	}
        	IntervalDer<T> PowDouble(const IntervalDer<T>& base, double exp) const { return base ^ Interval<T>(exp); }
		IntervalDer<T> Abs(const IntervalDer<T>& t) const { return abs(t); }
		IntervalDer<T> Ln(const IntervalDer<T>& t) const { return ln(t); }
		IntervalDer<T> Log(const IntervalDer<T>& t, double base) const { return log(t, base); }
		IntervalDer<T> Min(const IntervalDer<T>& left, const IntervalDer<T>& right) const { return min(left, right); }
		IntervalDer<T> Max(const IntervalDer<T>& left, const IntervalDer<T>& right) const { return max(left, right); }
		IntervalDer<T> IfTrue(IntervalBool condition, const IntervalDer<T>& left, const IntervalDer<T>& right) const { return ifThen(condition, left, right);  }
		IntervalDer<T> IfTrue(const ptrCNode<IntervalDer<T>>& conditionNode, const ptrNode<IntervalDer<T>>& left, const ptrNode<IntervalDer<T>>& right, MapIterator &map_iterator) const 
		{
			IntervalBool condition = conditionNode->calcCondition(*this, map_iterator);
			auto var = std::dynamic_pointer_cast<Var<IntervalDer<T>>>(left);
			auto cnst = std::dynamic_pointer_cast<Const<IntervalDer<T>>>(right);
			if(condition == IntervalBool::Intermadiate && var && cnst)
			{
				Conditions cond = conditionNode->getCondition();
				int index = var->getIndex();
				std::vector<Interval<T>> v(m_v);
				if(cond == Conditions::More || cond == Conditions::MoreEqual)					
					v[index] = Interval<T>(cnst->getConst(), v[index].rb());
				else
					v[index] = Interval<T>(v[index].lb(), cnst->getConst());
				IntervalDerAlg<T> alg(v);
				return ifThen(condition, left->calc(alg, map_iterator), right->calc(alg, map_iterator));
			}
			else
				return ifThen(condition, left->calc(*this, map_iterator), right->calc(*this, map_iterator));
		}
		IntervalBool Condition(Conditions condition, const IntervalDer<T>& left, const IntervalDer<T>& right) const { return cond(condition, left, right); }
		IntervalDer<T> CreateVar(int index) const 
		{ 
		    std::vector<Interval<T>> grad(m_v.size(), 0.0);
		    grad[index] = 1.0;
		    return IntervalDer<T>(m_v[index], Grad<Interval<T>>(grad));
		}      
		IntervalDer<T> CreateConst(double cnst) const { return IntervalDer<T>(cnst, Grad<Interval<T>>(m_v.size(), 0.0)); }
		vPtrAlg<IntervalDer<T>> GetNewAlgorithm(Conditions cond, int index, double cnst) const
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
			return vPtrAlg<IntervalDer<T>>({ ptrAlg<IntervalDer<T>>(new IntervalDerAlg<T>(leftVec)), ptrAlg<IntervalDer<T>>(new IntervalDerAlg<T>(rightVec))});
		}
	};    
}
}

#endif /* ALGDER__HPP */
