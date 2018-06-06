/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   interval_.hpp
 * Author: Alexander Usov
 *
 * Created on January 19, 2017, 10:51 AM
 */

#ifndef INTERVAL__HPP
#define INTERVAL__HPP

#include <iostream> 
#include <limits>
#define _USE_MATH_DEFINES
#include <cmath>
#include <exception>
#include <memory>
#include <algorithm> 
#include <limits>
#include "common/utilmacro.hpp"
#include "enums.h"

namespace snowgoose {
  namespace interval {

	template<class T> class Interval;
	template<class T> using IL = std::initializer_list<Interval<T>>;

    /**
    * This class implements interval arithmetic
    */
    template<class T> class Interval
    {
    public:
			/**
			* Constructor
			* @param lb low bound
			* @param rb upper bound
			*/
			Interval(const T lb, const T rb);
			/**
			* Constructor
			* @param pair with low and upper bounds
			*/
			Interval(const std::pair<T, T> &pair);
			/**
			* Constructor
			* @param t with equal low and upper bounds
			*/
			Interval(const T t);
			/**
			* Default constructor. Sets low and upper bound to 0.0.
			* 
			*/
            Interval();
			/**
			* Adds two intervals
			*/
			Interval operator+(const Interval &y) const;
			/**
			* Subtracts two intervals
			* @return interval
			*/
			Interval operator-(const Interval &y) const;
			/**
			* Raising to a power
			* @param exp is integer 
			* @return interval
			*/
			Interval operator^(int exp) const;
			/**
			* Raising to a power
			* @param exp is real exponent
			* @return interval
			*/
			Interval operator^(T exp) const;
			/**
			* Raising to a power
			* @param exp is interval exponent
			* @return interval
			*/
			Interval operator^(const Interval &exp) const;
			/**
			* Multiplies intervals
			* @param y is right interval
			* @return interval
			*/
			Interval operator*(const Interval &y) const;
			/**
			* Divides intervals
			* @param y is right interval
			* @return interval
			*/
			Interval operator/(const Interval &y) const;
			/**
			* Is left interval less than right one
			* @param y is right interval
			* @return IntervalBool 
			*/
			IntervalBool operator < (const Interval &y) const;
			/**
			* Is left interval less than right number
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator < (T y) const; 
			/**
			* Is left interval more than right one
			* @param y is right interval
			* @return IntervalBool
			*/
			IntervalBool operator > (const Interval &y) const;
			/**
			* Is left interval more then right the number
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator > (T y) const;
			/**
			* Is left interval less then or equal to the right interval
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator <= (const Interval &y) const;
			/**
			* Is left interval less then or equals to the right number
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator <= (T y) const;
			/**
			* Is left interval more then or equals to the right interval
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator >= (const Interval &y) const;
			/**
			* Is left interval more then or equal to the right number
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator >= (T y) const;
            /**
			* Is left interval equals to the right interval
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator == (const Interval &y) const;
			/**
			* Is left interval equals to the right number
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator == (T y) const;
            /**
			* Is left interval doesn't equal to the right interval
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator != (const Interval &y) const;
			/**
			* Is left interval doesn't equal to the right number
			* @param y is right number
			* @return IntervalBool
			*/
			IntervalBool operator != (T y) const;
			/**
			* Is left number less than right interval
			* @param x is left number
			* @param y is right number
			* @return IntervalBool
			*/
			template<class T2> friend IntervalBool operator < (T2 x, const Interval<T2> &y);
			/**
			* Is left number more than right interval
			* @param x is left number
			* @param y is right number
			* @return IntervalBool
			*/
			template<class T2> friend IntervalBool operator > (T2 x, const Interval<T2> &y);
			/**
			* Is left number less than or equal to right interval
			* @param x is left number
			* @param y is right number
			* @return IntervalBool
			*/
			template<class T2> friend IntervalBool operator <= (T2 x, const Interval<T2> &y);
			/**
			* Is left number more than or equal to right interval
			* @param x is left number
			* @param y is right number
			* @return IntervalBool
			*/
			template<class T2> friend IntervalBool operator >= (T2 x, const Interval<T2> &y);
			/**
			* Returns min interval from array.
			* @param list array of intervals
			* @return interval
			*/
			template<class T2> friend Interval<T2> min(const IL<T2>& list);
			/**
			* Returns max interval from array.
			* @param list array of intervals
			* @return interval
			*/
			template<class T2> friend Interval<T2> max(const IL<T2>& list);
			/**
			* Multiplies the number with the interval
			* @param t is left number
			* @param y is rigth interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> operator*(T2 t, const Interval<T2> &y);
			/**
			* Adds the number with the interval
			* @param t is left number
			* @param y is rigth interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> operator+(T2 t, const Interval<T2> &y);
			/**
			* Subtracts the interval from the number 
			* @param t is left number
			* @param y is rigth interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> operator-(T2 t, const Interval<T2> &y);
			/**
			* Divides a given number by the interval
			* @param t is left number
			* @param y is rigth interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> operator/(T2 t, const Interval<T2> &y);
			/**
			* Raising to a square power   
			* @param x is base
			* @return interval
			*/
			template<class T2> friend Interval<T2> sqr(const Interval<T2> &x);
			/**
			* Sinus of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> sin(const Interval<T2> &x);
			/**
			* Logarithm to the base of the interval
			* @param x is interval
			* @param base is log's base
			* @return interval
			*/
			template<class T2> friend Interval<T2> log(const Interval<T2> &x, double base);
			/**
			* Arc sinus of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> asin(const Interval<T2> &x);
			/**
			* cosine of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> cos(const Interval<T2> &x);
			/**
			* Arc cosine of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> acos(const Interval<T2> &x);
			/**
			* Tangent of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> tg(const Interval<T2> &x);
			/**
			* �otangent of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> ctg(const Interval<T2> &x);
			/**
			* Arc tangent of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> atg(const Interval<T2> &x);
			/**
			* Arc cotangent of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> actg(const Interval<T2> &x);
			/**
			* Raises the exp 2.718 to the power x
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> exp(const Interval<T2> &x);
			/**
			* Square root of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> sqrt(const Interval<T2> &x);
			/**
			* The absolute value of of the interval.
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> abs(const Interval<T2> &x);
			/**
			* The natural logarithm of the interval
			* @param x is interval
			* @return interval
			*/
			template<class T2> friend Interval<T2> ln(const Interval<T2> &x);
			/**
			* The conditional expression
			* @param ib is IntervalBool
			* @param x is left interval
			* @param y is right interval
			* @return left interval if true, right if false, union of intervals if Intermadiate
			*/
			template<class T2> friend Interval<T2> ifThen(IntervalBool ib, const Interval<T2> &x, const Interval<T2> &y);
			/**
			* Output the interval
			* @param out is output stream
			* @param x is interval
			* @return output stream
			*/
			template<class T2> friend std::ostream& operator<<(std::ostream & out, const Interval<T2> x);
			/**
			* Low bound
			* @return low bound.
			*/
			T lb() const;
			/**
			* Upper bound
			* @return upper bound.
			*/
			T rb() const;
			/**
			*Explicit cast operator. Returns low bound.
			* @return low bound.
			*/
			explicit operator T();
    private:
			/**
			* Raising to a power
			* @param n is unsigned integer exponent
			* @return interval
			*/
			Interval pow(unsigned int n) const;
			/**
			* low bound of the interval
			*/
			T m_lb;
			/**
			* upper bound of the interval
			*/
			T m_rb;
    };   

		template<class T> T Interval<T>::lb() const { return m_lb; }
		template<class T> T Interval<T>::rb() const { return m_rb; }

		template<class T> Interval<T>::Interval(const T lb, const T rb)
		{
			this->m_lb = lb < rb ? lb : rb;
			this->m_rb = rb > lb ? rb : lb;
		}
		template<class T> Interval<T>::Interval(const std::pair<T, T> &pair) : Interval(pair.first, pair.second) {}
		template<class T> Interval<T>::Interval(const T t) : Interval(t, t) {}
		template<class T> Interval<T>::Interval() { m_lb = 0.0; m_rb=0.0; }
		template<class T> Interval<T>::operator T(){ return this->m_lb; }
		template<class T> Interval<T> Interval<T>::operator+(const Interval &y) const
		{
			Interval z(m_lb + y.m_lb, m_rb + y.m_rb);
			return z;
		}
		template<class T> Interval<T> Interval<T>::operator-(const Interval &y) const
		{
			Interval z(m_lb - y.m_rb, m_rb - y.m_lb);
			return z;
		}
		template<class T> Interval<T> Interval<T>::operator^(int exp) const
		{
			if (exp >= 0)
				return pow((unsigned int)exp);
			else
				return 1.0 / pow(::abs(exp));
		}
		template<class T> Interval<T> Interval<T>::operator^(T exp) const
		{
			T lb = m_lb;
			T rb = m_rb;
			if (m_rb < 0)
				throw std::invalid_argument("The function Interval<T>::pow is not define for negative numbers");
			if (m_lb < 0)
				lb = 0.0;
			if (exp < 0) //decreases monotonically
			{
				if (lb == 0.0)
					throw std::invalid_argument("The function Interval<T>::pow is not define for base=0.0 and exponent<0");
				return Interval(std::pow(rb, exp), std::pow(lb, exp));
			}
			else //increases monotonically
				return Interval(std::pow(lb, exp), std::pow(rb, exp));
		}
		template<class T> Interval<T> Interval<T>::operator^(const Interval &y) const
		{
			if (y.m_lb == y.m_rb)
				return (*this) ^ y.m_lb;
			Interval a = (*this) ^ y.m_lb;
			Interval b = (*this) ^ y.m_rb;
			T lb = SGMIN(a.m_lb, b.m_lb);
			T rb = SGMAX(a.m_rb, b.m_rb);
			return Interval(lb, rb);
		}
		template<class T> Interval<T> Interval<T>::operator*(const Interval &y) const
		{
			T t1, t2, t3, t4, lb, rb;
			t1 = m_lb * y.m_lb;
			lb = t1;
			rb = t1;

			t2 = m_lb * y.m_rb;
			lb = SGMIN(lb, t2);
			rb = SGMAX(rb, t2);

			t3 = m_rb * y.m_lb;
			lb = SGMIN(lb, t3);
			rb = SGMAX(rb, t3);

			t4 = m_rb * y.m_rb;
			lb = SGMIN(lb, t4);
			rb = SGMAX(rb, t4);

			return Interval(lb, rb);
		}
		template<class T> Interval<T> Interval<T>::operator/(const Interval &y) const
		{
			if (y.m_lb <= 0.0 && y.m_rb >= 0.0)
				throw std::invalid_argument("Invalid operation. Interval is divided by interval that includes 0.0.");

			T t1, t2, t3, t4, lb, rb;
			t1 = m_lb / y.m_lb;
			lb = t1;
			rb = t1;

			t2 = m_lb / y.m_rb;
			lb = SGMIN(lb, t2);
			rb = SGMAX(rb, t2);

			t3 = m_rb / y.m_lb;
			lb = SGMIN(lb, t3);
			rb = SGMAX(rb, t3);

			t4 = m_rb / y.m_rb;
			lb = SGMIN(lb, t4);
			rb = SGMAX(rb, t4);

			return Interval(lb, rb);
		}
		template<class T> IntervalBool Interval<T>::operator < (const Interval &y) const { return m_rb < y.m_lb ? IntervalBool::True : m_lb > y.m_rb ? IntervalBool::False : IntervalBool::Intermadiate; }
		template<class T> IntervalBool Interval<T>::operator < (T y) const { return m_rb < y ? IntervalBool::True : m_lb > y ? IntervalBool::False : IntervalBool::Intermadiate; }		
		template<class T> IntervalBool operator < (T x, const Interval<T> &y) { return x < y.m_lb ? IntervalBool::True : x > y.m_rb ? IntervalBool::False : IntervalBool::Intermadiate; }
		template<class T> IntervalBool Interval<T>::operator > (const Interval &y) const { return y < (*this); }
		template<class T> IntervalBool Interval<T>::operator > (T y) const { return y < (*this); }
		template<class T> IntervalBool operator > (T x, const Interval<T> &y) { return y < x; }
		template<class T> IntervalBool Interval<T>::operator <= (const Interval &y) const { return m_rb <= y.m_lb ? IntervalBool::True : m_lb >= y.m_rb ? IntervalBool::False : IntervalBool::Intermadiate; }
		template<class T> IntervalBool Interval<T>::operator <= (T y) const { return m_rb <= y ? IntervalBool::True : m_lb >= y ? IntervalBool::False : IntervalBool::Intermadiate; }
		template<class T> IntervalBool operator <= (T x, const Interval<T> &y) { return x <= y.m_lb ? IntervalBool::True : x >= y.m_rb ? IntervalBool::False : IntervalBool::Intermadiate; }
		template<class T> IntervalBool Interval<T>::operator >= (const Interval &y) const { return y <= (*this); }
		template<class T> IntervalBool Interval<T>::operator >= (T y) const { return y <= (*this); }
        template<class T> IntervalBool Interval<T>::operator == (const Interval &y) const { return m_lb == y.m_lb && m_rb == y.m_rb ? IntervalBool::True : IntervalBool::False; }
		template<class T> IntervalBool Interval<T>::operator == (T y) const { return m_lb == y && m_rb == y ? IntervalBool::True : IntervalBool::False; }
        template<class T> IntervalBool Interval<T>::operator != (const Interval &y) const { return m_lb != y.m_lb || m_rb != y.m_rb ? IntervalBool::True : IntervalBool::False; }
		template<class T> IntervalBool Interval<T>::operator != (T y) const { return m_lb != y || m_rb != y ? IntervalBool::True : IntervalBool::False; }
		template<class T> IntervalBool operator >= (T x, const Interval<T> &y) { return y <= x; }

		template<class T> Interval<T> min(const IL<T>& list)
		{
			if (list.size() == 0)
				throw std::invalid_argument("Invalid argument. Empty initializer_list in min function.");
			auto lmin = std::min_element(list.begin(), list.end(), [](const Interval<T> &a, const Interval<T> &b) { return a.m_lb < b.m_lb; });
			auto rmin = std::min_element(list.begin(), list.end(), [](const Interval<T> &a, const Interval<T> &b) { return a.m_rb < b.m_rb; });
			return Interval<T>(lmin->m_lb, rmin->m_rb);
		}
		template<class T> Interval<T> max(const IL<T>& list)
		{
			if (list.size() == 0)
				throw std::invalid_argument("Invalid argument. Empty initializer_list in max function.");

			auto lmin = std::max_element(list.begin(), list.end(), [](const Interval<T> &a, const Interval<T> &b) { return a.m_lb < b.m_lb; });
			auto rmin = std::max_element(list.begin(), list.end(), [](const Interval<T> &a, const Interval<T> &b) { return a.m_rb < b.m_rb; });
			return Interval<T>(lmin->m_lb, rmin->m_rb);
		}
		template<class T> Interval<T> operator*(T t, const Interval<T> &y)
		{
			return Interval<T>(t) * y;
		}
		template<class T> Interval<T> operator+(T t, const Interval<T> &y)
		{
			return Interval<T>(t) + y;
		}
		template<class T> Interval<T> operator-(T t, const Interval<T> &y)
		{
			return Interval<T>(t) - y;
		}
		template<class T> Interval<T> operator/(T t, const Interval<T> &y)
		{
			return Interval<T>(t) / y;
		}
		template<class T> Interval<T> sqr(const Interval<T> &x)
		{
			return x ^ 2;
		}
		template<class T> Interval<T> ln(const Interval<T> &x)
		{
			if (x.m_rb < 0.0)
				throw std::invalid_argument("The function ln is not define for negative numbers");
			if (x.m_lb < 0.0)
				return Interval<T>(-std::numeric_limits<T>::infinity(), std::log(x.m_rb));
			T a = std::log(x.m_lb);
			T b = std::log(x.m_rb);
			T lb = SGMIN(a, b);
			T rb = SGMAX(a, b);
			return Interval<T>(lb, rb);
		}
		template<class T> Interval<T> log(const Interval<T> &x, double base)
		{
			if (x.m_rb <= 0.0)
				throw std::invalid_argument("The function log is not define for negative numbers and 0.0");
			if (base <= 0.0 || base == 1.0)
				throw std::invalid_argument("The function log is not define for negative base and if base equals 1.0");
			if (x.m_lb < 0.0)
				return Interval<T>(-std::numeric_limits<T>::infinity(), std::log(x.m_rb) / std::log(base));
			T a = std::log(x.m_lb) / std::log(base);
			T b = std::log(x.m_rb) / std::log(base);
			T lb = SGMIN(a, b);
			T rb = SGMAX(a, b);
			return Interval<T>(lb, rb);
		}
		template<class T> Interval<T> ifThen(IntervalBool ib, const Interval<T> &x, const Interval<T> &y)
		{
			if (ib == IntervalBool::True)
				return x;
			else if (ib == IntervalBool::False)
				return y;
			else
				return Interval<T>(SGMIN(x.m_lb, y.m_lb), SGMAX(x.m_rb, y.m_rb));
		}
		template<class T> std::ostream& operator<<(std::ostream & out, const Interval<T> x)
		{
			return out << "lb: " << x.m_lb << " rb: " << x.m_rb << '\n';
		}
		template<class T> Interval<T> sin(const Interval<T> &x)
		{
			int k, l;
			bool p2 = false;
			bool p32 = false;
			T lb, rb;

			k = (int)ceil(x.m_lb * M_2_PI);
			l = (int)floor(x.m_rb * M_2_PI);

			for (int i = k; i <= l; i++) {
				int test = i % 4;
				if ((i % 4) == 1 || (i % 4) == -3)
					p2 = true;
				else if ((i % 4) == 3 || (i % 4) == -1)
					p32 = true;
				if (p2 && p32)
					break;
			}
			if (p2)
				rb = 1.;
			else
				rb = SGMAX(std::sin(x.m_lb), std::sin(x.m_rb));
			if (p32)
				lb = -1.;
			else
				lb = SGMIN(std::sin(x.m_lb), std::sin(x.m_rb));

			return Interval<T>(lb, rb);
		}
		template<class T> Interval<T> asin(const Interval<T> &x)
		{
			if (x.m_rb < -1.0 || x.m_lb > 1.0)
				throw std::invalid_argument("Invalid argument in asin function.The argument is out of this interval[-1, 1].");
			T lb = x.m_lb < -1.0 ? -1.0 : x.m_lb;
			T rb = x.m_rb > 1.0 ? 1.0 : x.m_rb;
			return Interval<T>(std::asin(lb), std::asin(rb));
		}
		template<class T> Interval<T> cos(const Interval<T> &x)
		{
			T t1, t2;
			t1 = x.m_lb + M_PI_2;
			t2 = x.m_rb + M_PI_2;
			Interval<T> t(t1, t2);
			return sin(t);
		}
		template<class T> Interval<T> acos(const Interval<T> &x)
		{
			if (x.m_rb < -1.0 || x.m_lb > 1.0)
				throw std::invalid_argument("Invalid argument in acos function.The argument is out of this interval[-1, 1].");
			T lb = x.m_lb < -1.0 ? -1.0 : x.m_lb;
			T rb = x.m_rb > 1.0 ? 1.0 : x.m_rb;
			return Interval<T>(std::acos(rb), std::acos(lb));
		}
		template<class T> Interval<T> tg(const Interval<T> &x)
		{
			int k, l;
			bool p2 = false;

			k = (int)ceil(x.m_lb * M_2_PI);
			l = (int)floor(x.m_rb * M_2_PI);

			for (int i = k; i <= l; i++) {
				if ((i % 2) == 1 || (i % 2) == -1)
					p2 = true;
				if (p2)
					break;
			}
			if (p2)
				throw std::invalid_argument("Invalid argument in tg function. Tg is not defined for PI*N/2.");
			return Interval<T>(std::tan(x.m_lb), std::tan(x.m_rb));
		}
		template<class T> Interval<T> ctg(const Interval<T> &x)
		{
			int k, l;
			bool pi = false;

			k = (int)ceil(x.m_lb * M_2_PI);
			l = (int)floor(x.m_rb * M_2_PI);

			for (int i = k; i <= l; i++) {
				if ((i % 2) == 0)
					pi = true;
				if (pi)
					break;
			}

			if (pi)
				throw std::invalid_argument("Invalid argument in ctg function. Ctg is not defined for PI*N.");

			T lb = 1.0 / std::tan(x.m_rb);
			T rb = 1.0 / std::tan(x.m_lb);
			return Interval<T>(lb, rb);
		}
		template<class T> Interval<T> atg(const Interval<T> &x)
		{
			return Interval<T>(::atan(x.m_lb), std::atan(x.m_rb));
		}
		template<class T> Interval<T> actg(const Interval<T> &x)
		{
			T lb = M_PI_2 - std::atan(x.m_rb);
			T rb = M_PI_2 - std::atan(x.m_lb);
			return Interval<T>(lb, rb);
		}
		template<class T> Interval<T> exp(const Interval<T> &x)
		{
			T t1 = std::exp(x.m_lb);
			T t2 = std::exp(x.m_rb);
			T lb = SGMIN(t1, t2);
			T rb = SGMAX(t1, t2);
			return Interval<T>(lb, rb);
		}
		template<class T> Interval<T> sqrt(const Interval<T> &x)
		{
			if (x.m_rb < 0.0)
				throw std::invalid_argument("The function sqrt is not define for negative numbers");
			if (x.m_lb < 0.0)
				return Interval<T>(0.0, std::sqrt(x.m_rb));
			T a = std::sqrt(x.m_lb);
			T b = std::sqrt(x.m_rb);
			T lb = SGMIN(a, b);
			T rb = SGMAX(a, b);
			return Interval<T>(lb, rb);
		}
		template<class T> Interval<T> abs(const Interval<T> &x)
		{
			T t1 = std::abs(x.m_lb);
			T t2 = std::abs(x.m_rb);
			T lb = (SGMIN(x.m_lb, x.m_rb) <= 0.0 && SGMAX(x.m_lb, x.m_rb) >= 0) ? 0.0 : SGMIN(t1, t2);
			T rb = SGMAX(t1, t2);
			return Interval<T>(lb, rb);
		}
		template<class T> Interval<T> Interval<T>::pow(unsigned int n) const
		{
			T t1 = 1.0, t2 = 1.0, lb, rb;
			for (unsigned int i = 1; i <= n; i++) {
				t1 *= m_lb;
				t2 *= m_rb;
			}
			if (n % 2) {
				lb = t1;
				rb = t2;
			}
			else {
				if ((m_lb <= 0.0) && (m_rb >= 0.0))
					lb = 0.0;
				else
					lb = SGMIN(t1, t2);
				rb = SGMAX(t1, t2);
			}
			return Interval<T>(lb, rb);
		}       
    }
}

#endif /* INTERVAL__HPP */
