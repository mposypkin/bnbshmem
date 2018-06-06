#ifndef EXPR__HPP
#define EXPR__HPP

#include <iostream> 
#include <memory>
#include <functional>
#include <math.h>
#include "node.hpp"
#include "algorithm.hpp"
#include "algder.hpp"
#include "derivatives/valder.hpp"
#include "derivatives/intervalder.hpp"

using namespace snowgoose::derivative;

namespace snowgoose {
namespace expression {

	/**
	* This class implements work with mathematical expressions.
	* @param T is real numbers or intervals
	*/
	template <class T>
	class Expr
	{
	public:
		/**
		* Default constructor
		*/
		Expr();
		/**
		* Constructor
		* @param value is constant expression
		*/
		Expr(double value);
		/**
		* Adds two expressions
		* @param value is left expression
		* @return expression
		*/
		Expr operator+(const Expr& value);
		/**
		* Adds real number with expression
		* @param lv is left number
		* @param rv is right expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> operator+(double lv, const Expr<T2>& rv);
		/**
		* Adds expression to given one
		* @param rv is right expression
		* @return expression
		*/
		Expr& operator+=(const Expr& value);
		/**
		* Subtracts the expression from given one
		* @param value is left expression
		* @return expression
		*/
		Expr operator-(const Expr& value);
		/**
		* Subtracts the expression from the real number
		* @param lv is left number
		* @param rv is right expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> operator-(double lv, const Expr<T2>& rv);
		/**
		* Allows to write minus before expression
		* @return expression
		*/
		Expr operator-();
		/**
		* Multiplies expressions
		* @param value is right expression
		* @return expression
		*/
		Expr operator*(const Expr& value);
		/**
		* Multiplies the real number with the expression
		* @param lv is left expression
		* @param rv is right expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> operator*(double lv, const Expr<T2>& rv);
		/**
		* Divides expressions
		* @param value is right expression
		* @return expression
		*/
		Expr operator/(const Expr& value);
		/**
		* Divides a given number by the expression
		* @param lv is left number
		* @param rv is rigth expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> operator/(double lv, const Expr<T2>& rv);
		/**
		* Sinus of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> sin(const Expr<T2>& value);
		/**
		* Cosinus of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> cos(const Expr<T2>& value);
		/**
		* Tangent of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> tg(const Expr<T2>& value);
		/**
		* Cotangent of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> ctg(const Expr<T2>& value);
		/**
		* Arc cosinus of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> acos(const Expr<T2>& value);
		/**
		* Arc sinus of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> asin(const Expr<T2>& value);
		/**
		* Arc tangent of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> atg(const Expr<T2>& value);
		/**
		* Arc cotangent of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> actg(const Expr<T2>& value);
		/**
		* Raises exp 2.718 to the power value 
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> exp(const Expr<T2>& value);
		/**
		* Square root of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> sqrt(const Expr<T2>& value);
		/**
		* Raising to a square power
		* @param value is base
		* @return expression
		*/
		template <class T2> friend Expr<T2> sqr(const Expr<T2>& value);
		/**
		* Raising to a power
		* @param exp is integer exponent
		* @return expression
		*/
		Expr<T> operator^(int exp);
		/**
		* Raising to a power
		* @param exp is real exponent
		* @return expression
		*/
		Expr<T> operator^(double exp);
		/**
		* Raising to a power
		* @param exp is calculated by expression exponent
		* @return expression
		*/
		Expr<T> operator^(const Expr<T>& exp);
		/**
		* Raising to a power
		* @param exp is integer exponent
		* @param base
		* @return expression
		*/
		template <class T2> friend Expr<T2> pow(const Expr<T2>& base, int exp);
		/**
		* Raising to a power
		* @param exp is double exponent
		* @param base
		* @return expression
		*/
		template <class T2> friend Expr<T2> pow(const Expr<T2>& base, double exp);
		/**
		* Raising to a power
		* @param exp is calculated by expression exponent
		* @param base
		* @return expression
		*/
		template <class T2> friend Expr<T2> pow(const Expr<T2>& base, const Expr<T2>& exp);
		/**
		* The absolute value of of the expression.
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> abs(const Expr<T2>& value);
		/**
		* The natural logarithm of the expression
		* @param value is expression
		* @return expression
		*/
		template <class T2> friend Expr<T2> ln(const Expr<T2>& value);
		/**
		* Logarithm to the base of the expression
		* @param value is expression
		* @param base is log's base
		* @return expression
		*/
		template <class T2> friend Expr<T2> log(const Expr<T2>& value, double base);
		/**
		* Returns min expression from two ones.
		* @param lv is left expression
		* @param rv is right expression
		* @return interval
		*/
		template <class T2> friend Expr<T2> min(const Expr<T2>& lv, const Expr<T2>& rv);
		/**
		* Returns max expression from two ones.
		* @param lv is left expression
		* @param rv is right expression
		* @return interval
		*/
		template <class T2> friend Expr<T2> max(const Expr<T2>& lv, const Expr<T2>& rv);
		/**
		* The conditional expression
		* @param condition is expression which returns >,<,<=,>= operations. Don't use other type of expressions here.
		* @param lv is left expression
		* @param rv is right expression
		* @return left expression if true, right if false or union of expressions if Intermadiate
		*/
		template <class T2> friend Expr<T2> ifThen(const Expr<T2>& condition, const Expr<T2>& lv, const Expr<T2>& rv);
		/**
		* Is left expression more than right one
		* @param value is right expression
		* @return expression
		*/
		Expr<T> operator>(const Expr<T>& value);
		/**
		* Is left expression less than right one
		* @param value is right expression
		* @return expression
		*/
		Expr<T> operator<(const Expr<T>& value);
		/**
		* Is left expression more than or equal to right one
		* @param value is right expression
		* @return expression
		*/
		Expr<T> operator>=(const Expr<T>& value);
		/**
		* Is left expression less than or equal to right one
		* @param value is right expression
		* @return expression
		*/
		Expr<T> operator<=(const Expr<T>& value);

		/**
		* Index operator for vector variables
		* @param i is integer index
		* @return expression
		*/
		Expr<T> operator[](int i);
		/**
		* Index operator for vector variables
		* @param iterator
		* @return expression
		*/
		Expr<T> operator[](const Iterator &iterator);
		/**
		* Index operator for vector variables
		* @param expr is calculated expression
		* @return expression
		*/
		Expr<T> operator[](const Expr<T> &expr);
		/**
		* Index operator for vector variables
		* @param f is functor which returns integer index
		* @return expression
		*/
		Expr<T> operator[](const std::function<int()> &f);
		/**
		* Loop of adding
		* @param value is calculated expression
		* @param i is iterator
		* @return expression
		*/
		template <class T2> friend Expr<T2> loopSum(const Expr<T2>& value, const Iterator &i);
		/**
		* Loop of multiplication
		* @param value is calculated expression
		* @param i is iterator
		* @return expression
		*/
		template <class T2> friend Expr<T2> loopMul(const Expr<T2>& value, const Iterator &i);
		/**
		* calculates expressions according to algorithm
		* @param v is vector of variables. Notes, it can be vector of real numbers or intervals.
		* @param alg is algorithm. It allows to calculate either value of function or interval estimation of function.
		* @return real number or interval
		*/
		T calc(const Algorithm<T> & alg) const;
		/**
		* Output the expression
		* @param out is output stream
		* @param v is expression
		* @return output stream
		*/
		template <class T2> friend std::ostream& operator<<(std::ostream & out, const Expr<T2>& v);
		friend class Iterator;
		template <class T3> friend class ExprIndex;
	private:
		/**
		* It is root of tree of the expression
		*/
		ptrNode<T> node;
		/**
		* Private constructor
		* @param n is node 
		*/
		Expr(const ptrNode<T>& n);
	};

	template <class T> Expr<T>::Expr() {}
	template <class T> Expr<T>::Expr(double value) : node(ptrNode<T>(new Const<T>(value)))
	{
	}
	template <class T> Expr<T>::Expr(const ptrNode<T>& n) : node(n) 
	{
	}
	template <class T> Expr<T> Expr<T>::operator+(const Expr& value)
	{
		return Expr(ptrNode<T>(new Plus<T>(node, value.node)));
	}
	template <class T> Expr<T> operator+(double lv, const Expr<T>& rv)
	{
		ptrNode<T> pNode(new Const<T>(lv));
		return Expr<T>(ptrNode<T>(new Plus<T>(pNode, rv.node)));
	}
	template <class T> Expr<T>& Expr<T>::operator+=(const Expr<T>& value)
	{
		*this = *this + value;
		return *this;
	}
	template <class T> Expr<T> Expr<T>::operator-(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Minus<T>(node, value.node)));
	}
	template <class T> Expr<T> operator-(double lv, const Expr<T>& rv)
	{
		ptrNode<T> pNode(new Const<T>(lv));
		return Expr<T>(ptrNode<T>(new Minus<T>(pNode, rv.node)));
	}
	template <class T> Expr<T> Expr<T>::operator-()
	{
		return -1*(*this);
	}
	template <class T> Expr<T> Expr<T>::operator*(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Mul<T>(node, value.node)));
	}
	template <class T> Expr<T> operator*(double lv, const Expr<T>& rv)
	{
		ptrNode<T> pNode(new Const<T>(lv));
		return Expr<T>(ptrNode<T>(new Mul<T>(pNode, rv.node)));
	}
	template <class T> Expr<T> Expr<T>::operator/(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Div<T>(node, value.node)));
	}

	template <class T> Expr<T> operator/(double lv, const Expr<T>& rv)
	{
		ptrNode<T> pNode(new Const<T>(lv));
		return Expr<T>(ptrNode<T>(new Div<T>(pNode, rv.node)));
	}
	template <class T> Expr<T> sin(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Sin<T>(value.node)));
	}

	template <class T> Expr<T> cos(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Cos<T>(value.node)));
	}

	template <class T> Expr<T> tg(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Tg<T>(value.node)));
	}

	template <class T> Expr<T> ctg(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Ctg<T>(value.node)));
	}

	template <class T> Expr<T> acos(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new ArcCos<T>(value.node)));
	}

	template <class T> Expr<T> asin(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new ArcSin<T>(value.node)));
	}

	template <class T> Expr<T> atg(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new ArcTg<T>(value.node)));
	}

	template <class T> Expr<T> actg(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new ArcCtg<T>(value.node)));
	}

	template <class T> Expr<T> exp(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Exp<T>(value.node)));
	}

	template <class T> Expr<T> sqrt(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Sqrt<T>(value.node)));
	}

	template <class T> Expr<T> sqr(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Sqr<T>(value.node)));
	}

	template <class T> Expr<T> Expr<T>::operator^(int exp)
	{
		return pow(*this, exp);
	}

	template <class T> Expr<T> Expr<T>::operator^(double exp)
	{
		return pow(*this, exp);
	}

	template <class T> Expr<T> Expr<T>::operator^(const Expr<T>& exp)
	{
		return pow(*this, exp);
	}

	template <class T> Expr<T> pow(const Expr<T>& base, int exp)
	{
		return Expr<T>(ptrNode<T>(new PowInt<T>(base.node, exp)));
	}

	template <class T> Expr<T> pow(const Expr<T>& base, double exp)
	{
		return Expr<T>(ptrNode<T>(new Pow<T>(base.node, exp)));
	}

	template <class T> Expr<T> pow(const Expr<T>& base, const Expr<T>& exp)
	{
		return Expr<T>(ptrNode<T>(new PowExpr<T>(base.node, exp.node)));
	}

	template <class T> Expr<T> abs(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Abs<T>(value.node)));
	}

	template <class T> Expr<T> ln(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new Ln<T>(value.node)));
	}

	template <class T> Expr<T> log(const Expr<T>& value, double base)
	{
		return Expr<T>(ptrNode<T>(new Log<T>(value.node, base)));
	}

	template <class T> Expr<T> min(const Expr<T>& lv, const Expr<T>& rv)
	{
		return Expr<T>(ptrNode<T>(new Min<T>(lv.node, rv.node)));
	}

	template <class T> Expr<T> max(const Expr<T>& lv, const Expr<T>& rv)
	{
		return Expr<T>(ptrNode<T>(new Max<T>(lv.node, rv.node)));
	}

	template <class T> Expr<T> ifThen(const Expr<T>& condition, const Expr<T>& lv, const Expr<T>& rv)
	{
		ptrCNode<T> cond = std::dynamic_pointer_cast<ConditionNode<T>>(condition.node);
		if (!cond) throw std::invalid_argument("Invalid argument in ifThen function.");
		return Expr<T>(ptrNode<T>(new IfTrue<T>(cond, lv.node, rv.node)));
	}

	template <class T> Expr<T> Expr<T>::operator>(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new ConditionNode<T>(Conditions::More, node, value.node)));
	}
	template <class T> Expr<T> Expr<T>::operator<(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new ConditionNode<T>(Conditions::Less, node, value.node)));
	}
	template <class T> Expr<T> Expr<T>::operator>=(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new ConditionNode<T>(Conditions::MoreEqual, node, value.node)));
	}
	template <class T> Expr<T> Expr<T>::operator<=(const Expr<T>& value)
	{
		return Expr<T>(ptrNode<T>(new ConditionNode<T>(Conditions::LessEqual, node, value.node)));
	}

	template <class T> Expr<T> Expr<T>::operator[](int i)
	{
		return Expr<T>(ptrNode<T>(new Var<T>(i)));
	}

	template <class T> Expr<T> Expr<T>::operator[](const Iterator &iterator)
	{
		return Expr<T>(ptrNode<T>(new Index<T>(iterator)));
	}
	template <class T> Expr<T> Expr<T>::operator[](const Expr<T> &expr)
	{
		return Expr<T>(ptrNode<T>(new ExprIndex<T>(expr)));
	}
	template <class T> Expr<T> Expr<T>::operator[](const std::function<int()> &f)
	{
		return Expr<T>(ptrNode<T>(new CalcIndex<T>(f)));
	}

	template <class T> Expr<T> loopSum(const Expr<T>& value, const Iterator &i)
	{
		return Expr<T>(ptrNode<T>(new CycleSum<T>(value.node, i)));
	}
	template <class T> Expr<T> loopMul(const Expr<T>& value, const Iterator &i)
	{
		return Expr<T>(ptrNode<T>(new CycleMul<T>(value.node, i)));
	}

	template <class T> T Expr<T>::calc(const Algorithm<T> & alg) const
	{
		MapIterator map_iterator;
		return node->calc(alg, map_iterator);
	}

	template <class T> T calcFunc(const Expr<T>& exp, const std::vector<T>& point)
	{
		return exp.calc(FuncAlg<T>(point));
	}

        template <class T> Interval<T> calcInterval(const Expr<Interval<T>>& exp, const std::vector<Interval<T>>& box)
        {
		return exp.calc(InterEvalAlg<T>(box));
        }

	template <class T> ValDer<T> calcGrad(const Expr<ValDer<T>>& exp, const std::vector<T>& point)
	{
		return exp.calc(ValDerAlg<T>(point));
	}

	template <class T> IntervalDer<T> calcIntervalGrad(const Expr<IntervalDer<T>>& exp, const std::vector<Interval<T>>& box)
	{
		return exp.calc(IntervalDerAlg<T>(box));
	}

	template <class T> std::ostream& operator<<(std::ostream & out, const Expr<T>& v)
	{
		return out << *v.node;
	}
}
}

#endif /* EXPR__HPP */
