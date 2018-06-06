/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   intervalder.hpp
 * Author: alusov
 *
 * Created on June 22, 2017, 11:36 AM
 */

#ifndef INTERVALDER_HPP
#define INTERVALDER_HPP

#include <iostream> 
#include <cmath>
#include "interval/interval_air.hpp"
#include "interval/enums.h"
#include "grad.hpp"


using namespace snowgoose::interval;

namespace snowgoose {
  namespace derivative {
    template<class T> class IntervalDer
    {
        public:
            IntervalDer(const Interval<T> &val, const Grad<Interval<T>> &der) : m_val(val), m_der(der)
            {
            }
            IntervalDer(T val, const Grad<Interval<T>> &der) : m_val(val), m_der(der)
            {
            }
            
            using value_type = Interval<T>;

            IntervalDer operator+(const IntervalDer &y) const;
            IntervalDer operator-(const IntervalDer &y) const;
            IntervalDer operator*(const IntervalDer &y) const;
            IntervalDer operator/(const IntervalDer &y) const;           
            IntervalDer operator^(int exp) const;
            IntervalDer operator^(const Interval<T> &exp) const;
            
            template<class T2, class T3> friend IntervalDer<T2> operator+(T3 t, const IntervalDer<T2> &y);
            template<class T2, class T3> friend IntervalDer<T2> operator-(T3 t, const IntervalDer<T2> &y);
            template<class T2, class T3> friend IntervalDer<T2> operator*(T3 t, const IntervalDer<T2> &y);
            template<class T2, class T3> friend IntervalDer<T2> operator/(T3 t, const IntervalDer<T2> &y);            
            template<class T2, class T3> friend IntervalDer<T2> operator+(const IntervalDer<T2> &y, T3 t);
            template<class T2, class T3> friend IntervalDer<T2> operator-(const IntervalDer<T2> &y, T3 t);
            template<class T2, class T3> friend IntervalDer<T2> operator*(const IntervalDer<T2> &y, T3 t);
            template<class T2, class T3> friend IntervalDer<T2> operator/(const IntervalDer<T2> &y, T3 t);
                       
            template<class T2> friend IntervalDer<T2> sqr(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2> sqrt(const IntervalDer<T2> &x);           
            template<class T2> friend IntervalDer<T2> sin(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2> cos(const IntervalDer<T2>&x);            
            template<class T2> friend IntervalDer<T2> tg(const IntervalDer<T2> &x);   
            template<class T2> friend IntervalDer<T2> ctg(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2> asin(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2> acos(const IntervalDer<T2> &x);     
            template<class T2> friend IntervalDer<T2> atg(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2> actg(const IntervalDer<T2> &x); 
            template<class T2> friend IntervalDer<T2> ln(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2> log(const IntervalDer<T2> &x, double base);
            template<class T2> friend IntervalDer<T2> exp(const IntervalDer<T2> &x);
            template<class T2, class T3> friend IntervalDer<T2> pow(const T3 &base, const IntervalDer<T2> &exp);
            template<class T2> friend IntervalDer<T2> abs(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2> min(const IntervalDer<T2> &x, const IntervalDer<T2> &y);
            template<class T2> friend IntervalDer<T2> max(const IntervalDer<T2> &x, const IntervalDer<T2> &y);
            template<class T2> friend IntervalDer<T2> ifThen(IntervalBool condition, const IntervalDer<T2> &x, const IntervalDer<T2> &y);
            template<class T2> friend IntervalBool cond(Conditions condition, const IntervalDer<T2> &x, const IntervalDer<T2> &y);
            
            template<class T2> friend std::ostream& operator<<(std::ostream & out, const IntervalDer<T2> &x);
            Interval<T>  value() const {return m_val;}
            Grad<Interval<T>> grad() const {return m_der;}
            explicit operator T(){ return m_val.lb(); }
        private:
            Interval<T> m_val;
            Grad<Interval<T>> m_der;
    };
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator+(const IntervalDer &y) const
    {
        return IntervalDer(m_val + y.m_val, m_der + y.m_der);
    }
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator-(const IntervalDer &y) const
    {
        return IntervalDer(m_val - y.m_val, m_der - y.m_der);
    }
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator*(const IntervalDer &y) const
    {
        return IntervalDer(m_val * y.m_val, m_val * y.m_der + y.m_val * m_der);
    }
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator/(const IntervalDer &y) const
    {
        return IntervalDer(m_val / y.m_val, (m_der * y.m_val - m_val * y.m_der)/(y.m_val * y.m_val) );
    }
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator^(int exp) const
    {
        return IntervalDer(m_val^ exp, (double)exp * (m_val ^ (exp - 1.0)) * m_der);
    }
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator^(const Interval<T>  &exp) const
    {
        return IntervalDer(m_val^ exp, exp * (m_val ^ (exp - 1.0)) * m_der);
    }
    
    template<class T2, class T3> IntervalDer<T2>operator+(T3 t, const IntervalDer<T2> &y)
    {
        return IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0)) + y;
    }
    
    template<class T2, class T3> IntervalDer<T2>operator-(T3 t, const IntervalDer<T2> &y)
    {
        return IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0)) - y;
    }
    
    template<class T2, class T3> IntervalDer<T2> operator*(T3 t, const IntervalDer<T2> &y)
    {
        return IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0)) * y;
    }
    
    template<class T2, class T3> IntervalDer<T2> operator/(T3 t, const IntervalDer<T2> &y)
    {
        return IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0)) / y;
    }
    
    template<class T2, class T3> IntervalDer<T2> operator+(const IntervalDer<T2> &y, T3 t)
    {
        return t+y;
    }
    
    template<class T2, class T3> IntervalDer<T2> operator-(const IntervalDer<T2> &y, T3 t)
    {
        return y-IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0));
    }
    
    template<class T2, class T3> IntervalDer<T2> operator*(const IntervalDer<T2> &y, T3 t)
    {
        return t*y;
    }
    
    template<class T2, class T3> IntervalDer<T2> operator/(const IntervalDer<T2> &y, T3 t)
    {
        return y/IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0));
    }
    
    template<class T2> IntervalDer<T2> sqr(const IntervalDer<T2> &x)
    {
        return IntervalDer<T2>(x.m_val * x.m_val, 2.0 * x.m_val * x.m_der);
    }
    
    template<class T2> IntervalDer<T2> sqrt(const IntervalDer<T2> &x)
    {
        auto sq = interval::sqrt(x.m_val);
        return IntervalDer<T2>(sq, x.m_der/(2.0 * sq));
    }
    
    template<class T2> IntervalDer<T2> sin(const IntervalDer<T2> &x)
    {
        return IntervalDer<T2>(interval::sin(x.m_val), interval::cos(x.m_val) * x.m_der);
    }
    
    template<class T2> IntervalDer<T2> cos(const IntervalDer<T2> &x)
    {
        return IntervalDer<T2>(interval::cos(x.m_val), -1.0 * interval::sin(x.m_val) * x.m_der);
    }
      
    template<class T2> IntervalDer<T2> tg(const IntervalDer<T2> &x)
    {
        auto cs = interval::cos(x.m_val);
        return IntervalDer<T2>(interval::tg(x.m_val), x.m_der/(cs * cs));
    }
    
    template<class T2> IntervalDer<T2> ctg(const IntervalDer<T2> &x)
    {
        auto sn = interval::sin(x.m_val);
        return IntervalDer<T2>(interval::ctg(x.m_val), (-1.0 * x.m_der)/(sn * sn));
    }
    
    template<class T2> IntervalDer<T2> asin(const IntervalDer<T2> &x)
    {
        return IntervalDer<T2>(interval::asin(x.m_val), x.m_der/interval::sqrt(1.0 - x.m_val*x.m_val ));
    }
    
    template<class T2> IntervalDer<T2> acos(const IntervalDer<T2> &x)
    {
        return IntervalDer<T2>(interval::acos(x.m_val), (-1.0 * x.m_der)/interval::sqrt(1.0 - x.m_val*x.m_val ));
    }
    
    template<class T2> IntervalDer<T2> atg(const IntervalDer<T2> &x)
    {
        return IntervalDer<T2>(interval::atg(x.m_val), x.m_der/(1.0 + x.m_val*x.m_val ));
    }
    
    template<class T2> IntervalDer<T2> actg(const IntervalDer<T2> &x)
    {
        return IntervalDer<T2>(interval::actg(x.m_val), (-1.0 * x.m_der)/(1.0 + x.m_val*x.m_val ));
    }
    
    template<class T2> IntervalDer<T2> ln(const IntervalDer<T2> &x)
    {
        return IntervalDer<T2>(interval::ln(x.m_val), x.m_der/x.m_val);
    }
    
    template<class T2> IntervalDer<T2> log(const IntervalDer<T2> &x, double base)
    {
        return IntervalDer<T2>(interval::log(x.m_val, base), x.m_der/(x.m_val * std::log(base)));
    }
       
    template<class T2> IntervalDer<T2> exp(const IntervalDer<T2> &x)
    {
        auto ex = interval::exp(x.m_val);
        return IntervalDer<T2>(ex, x.m_der * ex);
    }
    
    template<class T2, class T3> IntervalDer<T2> pow(const T3 &base, const IntervalDer<T2> &exp)
    {
        auto pw = Interval<T2>(base) ^ exp.m_val;
        return IntervalDer<T2>(pw, exp.m_der * pw * interval::ln(base));
    }
    
    template<class T2> IntervalDer<T2> abs(const IntervalDer<T2> &x)
    {
        if(x.m_val==0.0)
            throw std::invalid_argument("Invalid argument in IntervalDer::get_abs. There isn't derivation at zero.");
        
        IntervalBool rez = x.m_val < 0.0;
        if(rez==IntervalBool::True)
            return IntervalDer<T2>(interval::abs(x.m_val), -1.0 * x.m_der);
        else if(rez==IntervalBool::False)
            return IntervalDer<T2>(interval::abs(x.m_val), x.m_der);
        else //IntervalBool::Intermediate case
            throw std::invalid_argument("Invalid argument in IntervalDer::get_abs. There isn't derivation at zero.");
    }
    
    template<class T2> IntervalDer<T2> min(const IntervalDer<T2> &x, const IntervalDer<T2> &y)
    {
        auto condition = (x.m_val <= y.m_val);
        if(condition==IntervalBool::True) 
            return x;
        else if(condition==IntervalBool::False) 
            return y;
        else
        throw std::invalid_argument("IntervalDer::min operation is not defined.");
    }
    
    template<class T2> IntervalDer<T2> max(const IntervalDer<T2> &x, const IntervalDer<T2> &y)
    {
        auto condition = (x.m_val >= y.m_val);
        if(condition==IntervalBool::True) 
            return x;
        else if(condition==IntervalBool::False) 
            return y;
        else
        throw std::invalid_argument("IntervalDer::max operation is not defined.");
    }
    
    template<class T2> IntervalDer<T2> ifThen(IntervalBool condition, const IntervalDer<T2> &x, const IntervalDer<T2> &y)
    {
        if(condition==IntervalBool::True) 
            return x;
        else if(condition==IntervalBool::False) 
            return y;
        else
            throw std::invalid_argument("IntervalDer::ifThen operation is not defined.");
    }
    
    template<class T2> IntervalBool cond(Conditions condition, const IntervalDer<T2> &x, const IntervalDer<T2> &y)
    {
        auto left = x.m_val;
        auto right = y.m_val;
        
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
        throw std::invalid_argument("Invalid condition in IntervalValDer::Condition.");
    }
    
    template<class T2> std::ostream& operator<<(std::ostream & out, const IntervalDer<T2> &x)
    {
        return std::cout << "val: " << x.m_val << " der: " << x.m_der << '\n';
    }
   
  }
}

#endif /* INTERVALDER_HPP */

