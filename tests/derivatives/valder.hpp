/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   valder.hpp
 * Author: alusov
 *
 * Created on June 22, 2017, 11:36 AM
 */

#ifndef VALDER_HPP
#define VALDER_HPP

#include <iostream> 
#include <cmath>
#include "grad.hpp"
#include "interval/enums.h"


namespace snowgoose {
  namespace derivative {
    template<class T> class ValDer
    {
        public:
            ValDer(T val, const Grad<T> &der) : m_val(val), m_der(der)
            {
            }
            
            
            using value_type = T;
            ValDer operator+(const ValDer &y) const;
            ValDer operator-(const ValDer &y) const;
            ValDer operator*(const ValDer &y) const;
            ValDer operator/(const ValDer &y) const;            
            ValDer operator^(int exp) const;
            ValDer operator^(T exp) const;
            
            template<class T2> friend ValDer<T2> operator+(T2 t, const ValDer<T2> &y);
            template<class T2> friend ValDer<T2> operator-(T2 t, const ValDer<T2> &y);
            template<class T2> friend ValDer<T2> operator*(T2 t, const ValDer<T2> &y);
            template<class T2> friend ValDer<T2> operator/(T2 t, const ValDer<T2> &y);           
            template<class T2> friend ValDer<T2> operator+(const ValDer<T2> &y, T2 t);
            template<class T2> friend ValDer<T2> operator-(const ValDer<T2> &y, T2 t);
            template<class T2> friend ValDer<T2> operator*(const ValDer<T2> &y, T2 t);
            template<class T2> friend ValDer<T2> operator/(const ValDer<T2> &y, T2 t);
                       
            template<class T2> friend ValDer<T2> sqr(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2> sqrt(const ValDer<T2> &x);            
            template<class T2> friend ValDer<T2> sin(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2> cos(const ValDer<T2>&x);            
            template<class T2> friend ValDer<T2> tg(const ValDer<T2> &x);  
            template<class T2> friend ValDer<T2> ctg(const ValDer<T2> &x); 
            template<class T2> friend ValDer<T2> asin(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2> acos(const ValDer<T2> &x); 
            template<class T2> friend ValDer<T2> atg(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2> actg(const ValDer<T2> &x); 
            template<class T2> friend ValDer<T2> ln(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2> log(const ValDer<T2> &x, double base);
            template<class T2> friend ValDer<T2> exp(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2> pow(const T2 &base, const ValDer<T2> &exp);
            template<class T2> friend ValDer<T2> abs(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2> min(const ValDer<T2> &x, const ValDer<T2> &y);
            template<class T2> friend ValDer<T2> max(const ValDer<T2> &x, const ValDer<T2> &y);
            template<class T2> friend ValDer<T2> ifThen(IntervalBool condition, const ValDer<T2> &x, const ValDer<T2> &y);
            template<class T2> friend IntervalBool cond(Conditions condition, const ValDer<T2> &x, const ValDer<T2> &y);
            
            template<class T2> friend std::ostream& operator<<(std::ostream & out, const ValDer<T2> &x);
            T value() const {return m_val;}
            Grad<T> grad() const { return m_der; }
            explicit operator T(){ return m_val; }
        private:
            T m_val;
            Grad<T> m_der;
    };
    
    template<class T> ValDer<T> ValDer<T>::operator+(const ValDer &y) const
    {
        return ValDer(m_val + y.m_val, m_der + y.m_der);
    }
    
    template<class T> ValDer<T> ValDer<T>::operator-(const ValDer &y) const
    {
        return ValDer(m_val - y.m_val, m_der - y.m_der);
    }
    
    template<class T> ValDer<T> ValDer<T>::operator*(const ValDer &y) const
    {
        return ValDer(m_val * y.m_val, m_val * y.m_der + y.m_val * m_der);
    }
    
    template<class T> ValDer<T> ValDer<T>::operator/(const ValDer &y) const
    {
        if (y.m_val == 0.0)
            throw std::invalid_argument("Invalid operation. Divide by 0.0.");
        return ValDer(m_val / y.m_val, (m_der * y.m_val - m_val * y.m_der)/(y.m_val * y.m_val) );
    }
    
    template<class T> ValDer<T> ValDer<T>::operator^(int exp) const
    {
        return ValDer(std::pow(m_val, exp), (T)exp * std:: pow(m_val, exp - 1.0) * m_der);
    }
    
    template<class T> ValDer<T> ValDer<T>::operator^(T exp) const
    {
        return ValDer(std::pow(m_val, exp), (T)exp * std:: pow(m_val, exp - 1.0) * m_der);
    }
    
    template<class T2> ValDer<T2>operator+(T2 t, const ValDer<T2>&y)
    {
        return ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0)) + y;
    }
    
    template<class T2> ValDer<T2>operator-(T2 t, const ValDer<T2>&y)
    {
        return ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0)) - y;
    }
    
    template<class T2> ValDer<T2>operator*(T2 t, const ValDer<T2>&y)
    {
        return ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0)) * y;
    }
    
    template<class T2> ValDer<T2>operator/(T2 t, const ValDer<T2>&y)
    {
        return ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0)) / y;;
    }
    
    template<class T2> ValDer<T2>operator+(const ValDer<T2>&y, T2 t)
    {
        return t+y;
    }
    
    template<class T2> ValDer<T2>operator-(const ValDer<T2>&y, T2 t)
    {
        return y-ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0));
    }
    
    template<class T2> ValDer<T2>operator*(const ValDer<T2>&y, T2 t)
    {
        return t*y;
    }
    
    template<class T2> ValDer<T2>operator/(const ValDer<T2>&y, T2 t)
    {
        return y/ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0));
    }
    
    template<class T2> ValDer<T2> sqr(const ValDer<T2>&x)
    {
        return ValDer<T2>(x.m_val * x.m_val, 2.0 * x.m_val * x.m_der);
    }
    
    template<class T2> ValDer<T2> sqrt(const ValDer<T2>&x)
    {
        if(x.m_val < 0.0)
            throw std::invalid_argument("The function ValDer<T>::sqrt is not defined for negative numbers");
        auto sq = std::sqrt(x.m_val);
        return ValDer<T2>(sq, x.m_der/(2.0 * sq));
    }
    
    template<class T2> ValDer<T2> sin(const ValDer<T2>&x)
    {
        return ValDer<T2>(std::sin(x.m_val), std::cos(x.m_val) * x.m_der);
    }
    
    template<class T2> ValDer<T2> cos(const ValDer<T2>&x)
    {
        return ValDer<T2>(std::cos(x.m_val), -1.0 * std::sin(x.m_val) * x.m_der);
    }
      
    template<class T2> ValDer<T2> tg(const ValDer<T2>&x)
    {
        auto cs = std::cos(x.m_val);
        return ValDer<T2>(std::tan(x.m_val), x.m_der/(cs * cs));
    }
    
    template<class T2> ValDer<T2> ctg(const ValDer<T2>&x)
    {
        auto sn = std::sin(x.m_val);
        return ValDer<T2>(1.0/std::tan(x.m_val), (-1.0 * x.m_der)/(sn * sn));
    }
    
    template<class T2> ValDer<T2> asin(const ValDer<T2>&x)
    {
        if(x.m_val < -1.0 || x.m_val > 1.0)
            throw std::invalid_argument("Invalid argument in ValDer::asin. The argument is out of this interval [-1,1]."); 
        return ValDer<T2>(std::asin(x.m_val), x.m_der/std::sqrt(1.0 - x.m_val*x.m_val ));
    }
    
    template<class T2> ValDer<T2> acos(const ValDer<T2>&x)
    {
        if(x.m_val < -1.0 || x.m_val > 1.0)
            throw std::invalid_argument("Invalid argument in ValDer::acos. The argument is out of this interval [-1,1].");        
        return ValDer<T2>(std::acos(x.m_val), (-1.0 * x.m_der)/std::sqrt(1.0 - x.m_val*x.m_val ));
    }
    
    template<class T2> ValDer<T2> atg(const ValDer<T2>&x)
    {
        return ValDer<T2>(std::atan(x.m_val), x.m_der/(1.0 + x.m_val*x.m_val ));            
    }
    
    template<class T2> ValDer<T2> actg(const ValDer<T2>&x)
    {
        return ValDer<T2>(M_PI_2 - std::atan(x.m_val), (-1.0 * x.m_der)/(1.0 + x.m_val*x.m_val ));
    }
    
    template<class T2> ValDer<T2> ln(const ValDer<T2>&x)
    {
        if (x.m_val <= 0)
            throw std::invalid_argument("The function ValDer<T>::ln is not define for negative numbers and 0.0");
        return ValDer<T2>(std::log(x.m_val), x.m_der/x.m_val);
    }
    
    template<class T2> ValDer<T2> log(const ValDer<T2> &x, double base)
    {
        if (x.m_val <= 0)
            throw std::invalid_argument("The function ValDer<T>::ln is not define for negative numbers and 0.0");

        auto lgb = std::log(base);
        return ValDer<T2>(std::log(x.m_val)/lgb, x.m_der/(x.m_val * lgb));
    }
    
    template<class T2> ValDer<T2> exp(const ValDer<T2>&x)
    {
        auto ex = std::exp(x.m_val);
        return ValDer<T2>(ex, x.m_der * ex);
    }
    
    template<class T2> ValDer<T2> pow(const T2 &base, const ValDer<T2> &exp)
    {
        auto pw = std::pow(base, exp.m_val);
        return ValDer<T2>(pw, exp.m_der * pw * std::log(base));
    }
    
    template<class T2> ValDer<T2> abs(const ValDer<T2>&x)
    {
        if(x.m_val==0.0)
            throw std::invalid_argument("Invalid argument in ValDer::abs. There isn't derivation at zero.");
        else if(x.m_val < 0.0)
            return ValDer<T2>(std::abs(x.m_val), -1.0 * x.m_der);
        else
            return ValDer<T2>(std::abs(x.m_val), x.m_der);
    }
    
    template<class T2> ValDer<T2> min(const ValDer<T2> &x, const ValDer<T2> &y)
    {
        return (x.m_val < y.m_val) ? x : y;
    }
    
    template<class T2> ValDer<T2> max(const ValDer<T2> &x, const ValDer<T2> &y)
    {
        return (x.m_val > y.m_val) ? x : y;
    }
            
    template<class T2> ValDer<T2> ifThen(IntervalBool condition, const ValDer<T2> &x, const ValDer<T2> &y)
    {
        return (condition == IntervalBool::True)? x : y;
    }
    
    template<class T2> IntervalBool cond(Conditions condition, const ValDer<T2> &x, const ValDer<T2> &y)
    {
        auto left = x.m_val;
        auto right = y.m_val;
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
        throw std::invalid_argument("Invalid condition in ValDer::Condition.");
    }
        
    template<class T2> std::ostream& operator<<(std::ostream & out, const ValDer<T2> &x)
    {
        return std::cout << "val: " << x.m_val << " der: " << x.m_der << '\n';
    }
   
  }
}

#endif /* VALDER_HPP */

