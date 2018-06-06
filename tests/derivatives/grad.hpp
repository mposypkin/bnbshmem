/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   grad.hpp
 * Author: alusov
 *
 * Created on June 24, 2017, 8:46 PM
 */

#ifndef GRAD_HPP
#define GRAD_HPP

#include <iostream>
#include <initializer_list>
#include <vector>
#include <algorithm>
#include <memory>

namespace snowgoose {
  namespace derivative {

        template <class  T> class Grad
        {
            public:
                Grad(const std::initializer_list<T> &lst);
                Grad(const std::vector<T> &grad) : m_grad(grad) {}
                Grad(std::size_t size, const T &t) : m_grad(size, t) {}
                T operator[](std::size_t i) const;
                Grad operator+(const Grad &y) const;
                Grad operator-(const Grad &y) const;
                Grad operator*(const T &y) const;
                template<class T2, class T3> friend Grad<T2> operator*(const T3 &x, const Grad<T2> &y);
                Grad operator/(const T &y) const;
                template<class T2> friend std::ostream& operator<<(std::ostream & out, const Grad<T2> &y);
                std::size_t size() const {return m_grad.size();}
                bool IsZero() const{
                    for(const T &d : m_grad)
                        if(d!=0.0) 
                            return false;
                    return true;
                }
		const T* getGrad() const { return m_grad.data(); };
            private:
                std::vector<T> m_grad;   
        };

        template<class T> Grad<T>::Grad(const std::initializer_list<T> &lst) : m_grad(lst)
        {
        }
       
        template<class T> T Grad<T>::operator[](std::size_t i) const
        {
            return m_grad[i];
        }
        
        template<class T> Grad<T> Grad<T>::operator+(const Grad &y) const //error control
        {
            std::size_t sz = m_grad.size();
            std::vector<T> grad(sz);
            for(int i=0; i < sz; i++)
                grad[i] = m_grad[i]+y.m_grad[i];
            return Grad<T>(grad);
        }

        template<class T> Grad<T> Grad<T>::operator-(const Grad &y) const //error control
        {
            std::size_t sz = m_grad.size();
            std::vector<T> grad(sz);
            for(int i=0; i < sz; i++)
                grad[i] = m_grad[i]-y.m_grad[i];
            return Grad<T>(grad);
        }

        template<class T> Grad<T> Grad<T>::operator*(const T &y) const
        {
            std::size_t sz = m_grad.size();
            std::vector<T> grad(sz);
            for(int i=0; i < sz; i++)
                grad[i]=m_grad[i]*y;
            return grad;
        }

        template<class T2, class T3> Grad<T2> operator*(const T3 &x, const Grad<T2> &y)
        {
            std::size_t sz = y.m_grad.size();
            std::vector<T2> grad(sz);
            for(int i=0; i < sz; i++)
                grad[i]=y.m_grad[i]*x;
            return grad;
        }

        template<class T> Grad<T> Grad<T>::operator/(const T &y) const //error control
        {
            std::size_t sz = m_grad.size();
            std::vector<T> grad(sz);
            for(int i=0; i < sz; i++)
                grad[i]=m_grad[i]/y;
            return grad;
        }

        template<class T2> std::ostream& operator<<(std::ostream & out, const Grad<T2> &y)
        {
            std::size_t sz = y.m_grad.size();
            for(int i=0; i< sz; i++)
                std::cout << y[i] << ' ';
            return out;
        }
    }
}


#endif /* GRAD_HPP */

