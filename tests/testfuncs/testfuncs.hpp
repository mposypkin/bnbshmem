#ifndef TEST_FUNCS__HPP
#define TEST_FUNCS__HPP

#include <math.h>
#include <iostream>
#include <vector>
#include "expression/expr.hpp"
#include "common/utilmacro.hpp"

namespace OPTITEST {

using namespace snowgoose::expression;

	template <class T>
	Expr<T> Ackley1(int N)
	{
		Expr<T> x;
		Iterator i(0, N - 1);
		Expr<T> y = -20 * exp(-0.2 * sqrt((1.0 / N) * loopSum(sqr(x[i]), i))) - exp((1.0 / N) * loopSum(cos(2 * M_PI * x[i]), i)) + 20.0 + M_E;
		return y;
	}

	template <class T>
	Expr<T> Ackley2(int N)
	{
		Expr<T> x;
		Iterator i(0, N - 1);
		Expr<T> y = -200 * exp(-0.02 * sqrt(loopSum(sqr(x[i]), i)));
		return y;
	}

	template <class T>
	Expr<T> Ackley3()
	{
		Expr<T> x;
		return 200 * exp(-0.02 * sqrt(sqr(x[0]) + sqr(x[1]))) + 5 * exp(cos(3*x[0]) + sin(3*x[1]));
	}

	template <class T>
	Expr<T> Ackley4()
	{
		int n = 2;
		Expr<T> x;
		Iterator i(0, n - 2);
		Expr<T> t = (Expr<T>)i + 1;
		return loopSum(exp(1.0/pow(sqr(x[i]) + sqr(x[t]) ,5.0)) + 3*(cos(2*x[i]) + sin(2*x[t])), i);
	}

	template <class T>
	Expr<T> Adjiman()
	{
		Expr<T> x;
		return cos(x[0]) * sin(x[1]) - x[0] / (sqr(x[1]) + 1);
	}

	template <class T>
	Expr<T> Alpine1(int N)
	{
		Expr<T> x;
		Iterator i(0, N - 1);
		return loopMul(abs(x[i] * sin(x[i] + 0.1 * x[i])), i);
	}

	template <class T>
	Expr<T> Alpine2()
	{
        int N = 3;
		Expr<T> x;
		Iterator i(0, N - 1);
		return loopMul(sqrt(x[i]) * sin(x[i]), i);
	}

	template <class T>
	Expr<T> Brad()
	{
		std::vector<double> vConst = { 0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.10, 4.39 };
		Expr<T> result = 0.0;
		Expr<T> y;
		for (int i = 0; i < 15; i++)
		{
			int I = i + 1;
			int J = 16 - I;
			Expr<T> a = vConst[i] - y[0] - I;
			Expr<T> b = J * y[1] + SGMIN(I, J) * y[2];
			result = result + sqr(a / b);
		}
		return result;
	}

	template <class T>
	Expr<T> BartelsConn()
	{
		Expr<T> x;
		return abs(sqr(x[0]) + sqr(x[1]) + x[0] * x[1]) + abs(sin(x[0])) + abs(cos(x[1]));
	}

	template <class T>
	Expr<T> Beale()
	{
		Expr<T> x;
		return sqr(1.5 - x[0] + x[0] * x[1]) + sqr(2.25 - x[0] + x[0] * sqr(x[1])) + sqr(2.625 - x[0] + x[0] * (x[1] ^ 3));
	}
	

	template <class T>
	Expr<T> BiggsExpr2()
	{
		Expr<T> x;
		Iterator i(1, 10);
		Expr<T> t = 0.1 * (Expr<T>)i;
		Expr<T> y = exp(-t) - 5 * exp(-10 * t);
		return loopSum(sqr(exp(-t * x[0]) - 5 * exp(-t * x[1]) - y), i);
	}

	template <class T>
	Expr<T> BiggsExpr3()
	{
		Expr<T> x;
		Iterator i(1, 10);
		Expr<T> t = 0.1 * (Expr<T>)i;
		Expr<T> y = exp(-t) - 5 * exp(-10 * t);
		return loopSum(sqr(exp(-t * x[0]) - x[2] * exp(-t * x[1]) - y), i);
	}

	template <class T>
	Expr<T> BiggsExpr4()
	{
		Expr<T> x;
		Iterator i(1, 10);
		Expr<T> t = 0.1 * (Expr<T>)i;
		Expr<T> y = exp(-t) - 5 * exp(-10 * t);
		return loopSum(sqr(x[2] * exp(-t * x[0]) - x[3] * exp(-t * x[1]) - y), i);
	}

	template <class T>
	Expr<T> BiggsExpr5()
	{
		Expr<T> x;
		Iterator i(1, 11);
		Expr<T> t = 0.1 * (Expr<T>)i;
		Expr<T> y = exp(-t) - 5 * exp(-10 * t) + 3 * exp(-4 * t);
		return loopSum(sqr(x[2] * exp(-t * x[0]) - x[3] * exp(-t * x[1]) + 3 * exp(-t * x[4]) - y), i);
	}

	template <class T>
	Expr<T> BiggsExpr6()
	{
		Expr<T> x;
		Iterator i(1, 13);
		Expr<T> t = 0.1 * (Expr<T>)i;
		Expr<T> y = exp(-t) - 5 * exp(-10 * t) + 3 * exp(-4 * t);
		return loopSum(sqr(x[2] * exp(-t * x[0]) - x[3] * exp(-t * x[1]) + x[5] * exp(-t * x[4]) - y), i);
	}

	template <class T>
	Expr<T> Bird()
	{
		Expr<T> x;
		return sin(x[0]) * exp(sqr(1 - cos(x[1]))) + cos(x[1]) * exp(sqr(1 - sin(x[0]))) + sqr(x[0] - x[1]);
	}

	template <class T>
	Expr<T> Bohachevsky1()
	{
		Expr<T> x;
		return sqr(x[0]) + 2 * sqr(x[1]) - 0.3 * cos(3 * M_PI * x[0]) - 0.4 * cos(4 * M_PI * x[1]) + 0.7;
	}

	template <class T>
	Expr<T> Bohachevsky2()
	{
		Expr<T> x;
		return sqr(x[0]) + 2 * sqr(x[1]) - 0.3 * cos(3 * M_PI * x[0]) * cos(4 * M_PI * x[1]) + 0.3;
	}

	template <class T>
	Expr<T> Bohachevsky3()
	{
		Expr<T> x;
		return sqr(x[0]) + 2 * sqr(x[1]) - 0.3 * cos(3 * M_PI * x[0] + 4 * M_PI * x[1]) + 0.3;
	}

	template <class T>
	Expr<T>Booth()
	{
		Expr<T> x;
		return sqr(x[0] + 2 * x[1] - 7) + sqr(2 * x[0] + x[1] - 5);
	}

	template <class T>
	Expr<T> BoxBettsQuadraticSum()
	{
		Expr<T> x;
		Iterator i(1, 10);
		Expr<T> t = -0.1 * Expr<T>(i);
                Expr<T> a = exp(t*x[0]) - exp(t*x[1]);
                Expr<T> b = x[2] * (exp(t) - exp(10 * t));
		return loopSum(sqr(a-b), i);
	}

	template <class T>
	Expr<T>BraninRCOS()
	{
		Expr<T> x;
		Expr<T> a = sqr(x[1] - 5.1*sqr(x[0]) / (4 * ::pow(M_PI,2)) + 5 * x[0] / M_PI - 6); 
		Expr<T> b = 10 * (1 - 1.0 / (8 * M_PI)) * cos(x[0]) + 10;
		return a + b;
	}

	template <class T>
	Expr<T>BraninRCOS2()
	{
		Expr<T> x;
		Expr<T> f1 = sqr(x[1] - 5.1 * sqr(x[0]) / (4 * M_PI * M_PI) + 5 * x[0] / M_PI - 6);
		Expr<T> f2 = 10 * (1 - 1.0 / (8 * M_PI)) * cos(x[0]) * cos(x[1]);
                Expr<T> f3 = ln(sqr(x[0]) + sqr(x[1]) + 1);
		return -1.0/(f1+f2+f3+10);
	}


	template <class T>
	Expr<T>Brent()
	{
		Expr<T> x;
		return sqr(x[0] + 10) + sqr(x[1] + 10) + exp(-sqr(x[0]) - sqr(x[1]));
	}


	template <class T>
	Expr<T> Brown(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 2);
		Expr<T> t = (Expr<T>)i + 1;
		return loopSum(pow(sqr(x[i]), sqr(x[t]) + 1) + pow(sqr(x[t]), sqr(x[i]) + 1), i);
	}

	template <class T>
	Expr<T> Bukin2()
	{
		Expr<T> x;
		return 100 * sqr(x[1] - 0.01 * sqr(x[0]) + 1) + 0.01 * sqr(x[0] + 10);
	}

	template <class T>
	Expr<T> Bukin4()
	{
		Expr<T> x;
		return 100 * sqr(x[1])  + 0.01 * abs(x[0] + 10);
	}

	template <class T>
	Expr<T> Bukin6()
	{
		Expr<T> x;
		return 100 * sqrt(abs(x[1] - 0.01 * sqr(x[0]))) + 0.01 * abs(x[0] + 10);
	}

	template <class T>
	Expr<T> CamelThreeHump()
	{
		Expr<T> x;
		return 2 * sqr(x[0]) - 1.05*(x[0] ^ 4) + (x[0] ^ 6) / 6.0 + x[0] * x[1] + sqr(x[1]);
	}

	template <class T>
	Expr<T> CamelSixHump()
	{
		Expr<T> x;
		return (4 - 2.1*sqr(x[0]) + (x[0] ^ 4) / 3)*sqr(x[0]) + x[0] * x[1] + (4 * sqr(x[1]) - 4) * sqr(x[1]);
	}

	template <class T>
	Expr<T> Chichinadze()
	{
		Expr<T> x;
		Expr<T> a = sqr(x[0]) - 12.0 * x[0] + 11.0;
		Expr<T> b = 10 * cos(M_PI_2 * x[0]) + 8*sin(5*x[0] * M_PI_2);
		Expr<T> c = 0.2 * ::sqrt(5.0) / exp(0.5*sqr(x[1] - 0.5));
		return a + b - c;
	}

	template <class T>
	Expr<T> ChungReynolds(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return sqr(loopSum(sqr(x[0]), i));
	}

	template <class T>
	Expr<T> Colville()
	{
		Expr<T> x;
		Expr<T> a = 100 * sqr(x[0] - sqr(x[1])) + sqr(1 - x[0]);
		Expr<T> b = 90 * sqr(x[3] - sqr(x[2])) + sqr(1 - x[2]);
		Expr<T> c = 10.1 * (sqr(x[1] - 1) + sqr(x[3] - 1));
		Expr<T> d = 19.8 * (x[1] - 1) * (x[3] - 1);
		return a + b + c + d;
	}

	template <class T>
	Expr<T> Complex()
	{
		Expr<T> x;
		Expr<T> a = pow(x[0],3);
		Expr<T> b = pow(x[1],3);
		Expr<T> c = 3.0*x[0]*sqr(x[1]);
		return sqr(a-c-1) + sqr(c-b);
	}

	template <class T>
	Expr<T> CosineMixture()
	{
		int n=2;
		Expr<T> x;
		Iterator i(0, n - 1);
		return -0.1* loopSum(cos(5 * M_PI*x[i]), i) + loopSum(sqr(x[i]), i);
	}

        template <class T>
	Expr<T> CrossInTray()
	{
		Expr<T> x;
		Expr<T> a = sqrt(sqr(x[0])+sqr(x[1]));
		Expr<T> b = abs(100.0 - a/M_PI); 
		Expr<T> c = abs(sin(x[0])*sin(x[1])*exp(b));
		Expr<T> d = -0.0001 * pow(c + 1.0, 0.1);
        return d;
	}
	template <class T>
	Expr<T> CrossLeg()
	{
		Expr<T> x;
		Expr<T> a = sqrt(sqr(x[0])+sqr(x[1]));
		Expr<T> b = abs(100.0 - a/M_PI); 
		Expr<T> c = abs(sin(x[0])*sin(x[1])*exp(b));
		Expr<T> d = -1.0/pow(c + 1.0, 0.1);
		return d;
	}

	template <class T>
	Expr<T> Cube()
	{
		Expr<T> x;
		return 100.0 * sqr(x[1] - pow(x[0] , 3)) + sqr(1 - x[0]);
	}

	template <class T>
	Expr<T> Deb1(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return (-1.0 / n) * loopSum(pow(sin(5 * M_PI*x[i]) , 6), i);
	}

	template <class T>
	Expr<T> Davis()
	{
		Expr<T> x;
		Expr<T> a = pow(sqr(x[0])+sqr(x[1]),0.25);
		Expr<T> b = 50*pow(3 * sqr(x[0])+sqr(x[1]) ,0.1);
		Expr<T> c = a*sqr(sin(b)+1);
		return c;
	}

	template <class T>
	Expr<T> DeckkersAarts()
	{
		Expr<T> x;
		return 100000 * sqr(x[0]) + sqr(x[1]) - sqr(sqr(x[0]) + sqr(x[1])) + 0.00001 * pow(sqr(x[0]) + sqr(x[1]), 4);
	}

	template <class T>
	Expr<T> DixonPrice()
	{
		Expr<T> x;
		return sqr(x[0] - 1) + 2.0 * sqr(2 * sqr(x[1]) - x[0]);
	}

	template <class T>
	Expr<T> Dolan()
	{
		Expr<T> x;
		return (x[0] + 1.7 * x[1]) * sin(x[0]) - 1.5 * x[2] - 0.1 * x[3] * cos(x[3] + x[4] - x[0]) + 0.2 * sqr(x[4]) - x[1] - 1.0;
	}

	template <class T>
	Expr<T> DropWave()
	{
		Expr<T> x;
		Expr<T> a = sqr(x[0])+sqr(x[1]);
		Expr<T> b = 1+cos(12*sqrt(a));
		Expr<T> c = 0.5*a + 2;
		return -b/c;
	}

	template <class T>
	Expr<T> Easom()
	{
		Expr<T> x;
		return -cos(x[0]) * cos(x[1]) * exp(-sqr(x[0] - M_PI) - sqr(x[1] - M_PI));
	}

	template <class T>
	Expr<T> ElAttarVidyasagarDutt()
	{
		Expr<T> x;
		Expr<T> a = sqr(sqr(x[0]) + x[1] - 10);
		Expr<T> b = sqr(x[0] + sqr(x[1]) - 7);
		Expr<T> c = sqr(sqr(x[0]) + (x[1] ^ 3) - 1);
		return a + b + c;
	}

	template <class T>
	Expr<T> Engvall()
	{
		Expr<T> x;
		return pow(x[0],4) + pow(x[1],4) + 2.0 * sqr(x[0]) * sqr(x[1]) - 4*x[0] + 3.0;
	}

	template <class T>
	Expr<T> EggCrate()
	{
		Expr<T> x;
		return sqr(x[0]) + sqr(x[1]) + 25 * (sqr(sin(x[0])) + sqr(sin(x[1])));
	}

	template <class T>
	Expr<T> EggHolder()
	{
		int n = 2;
		Expr<T> x;
		Iterator i(0, n - 2);
		Expr<T> t = i;
		return loopSum(-(x[t+1] + 47) * sin(sqrt(abs(x[t+1] + x[i] / 2.0 + 47))) - x[i] * sin(sqrt(abs(x[i] - x[t+1] - 47))), i);
	}

	template <class T>
	Expr<T> Exponential(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return -exp(-0.5 * loopSum( sqr(x[i]) , i));
	}

	template <class T>
	Expr<T> Exp2()
	{
		Expr<T> x;
		Iterator i(0, 9);
		Expr<T> t = -0.1 * (Expr<T>)i;
		return loopSum(sqr(exp(t * x[0]) - 5 * exp(t * x[1]) - exp(t) + 5 * exp(-(Expr<T>)i)), i);
	}

	template <class T>
	Expr<T> F26()
	{
		Expr<T> x;
		return 0.25 * pow(x[0],4) - 0.5 * sqr(x[0]) + 0.1 * x[0]  + 0.5 * sqr(x[1]);
	}

	template <class T>
	Expr<T> FreudensteinRoth()
	{
		Expr<T> x;
		Expr<T> a = sqr(x[0] - 13.0 + ((5 - x[1]) * x[1] - 2.0) * x[1]);
		Expr<T> b = sqr(x[0] - 29.0 + ((x[1] + 1) * x[1] - 14.0) * x[1]);
		return a + b;
	}

	template <class T>
	Expr<T> Giunta()
	{
		Expr<T> x;
		Iterator i(0, 1);
		const double cnst = 16.0 / 15.0;
		Expr<T> a = sin(cnst * x[i] - 1.0);
		Expr<T> b = sqr(sin(cnst * x[i] - 1.0));
		Expr<T> c = (1.0/50.0) * sin(4.0 * (cnst * x[i] - 1.0));
		return 0.6 + loopSum(a + b + c, i);
	}

	template <class T>
	Expr<T> GoldsteinPrice()
	{
		Expr<T> x;
		Expr<T> a = sqr(x[0] + x[1] + 1);
		Expr<T> b = 19 - 14 * x[0] + 3 * sqr(x[0]) - 14 * x[1] + 6 * x[0] * x[1] + 3 * sqr(x[1]);
		Expr<T> c = sqr(2*x[0] - 3*x[1]);
		Expr<T> d = 18 - 32 * x[0] + 12 * sqr(x[0]) + 48 * x[1] - 36 * x[0] * x[1] + 27 * sqr(x[1]);
		return (1 + a * b) * (30 + c * d);
	}

	template <class T>
	Expr<T> GramacyLee2()
	{
		Expr<T> x;
		Expr<T> a = -sqr(x[0])-sqr(x[1]);
		return x[0]*exp(a);
	}

	template <class T>
	Expr<T> GramacyLee3()
	{
		Expr<T> x;
		Iterator i(0, 1);

		Expr<T> a = -sqr(x[i]-1.0);
                Expr<T> b = -0.8 * sqr(x[i]+1.0);
		Expr<T> c = 8.0 * x[i]+0.8;
		
		return -loopMul(exp(a)+exp(b)-0.05*sin(c), i);
	}

	template <class T>
	Expr<T> Griewank(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return loopSum(sqr(x[i]) / 4000.0, i) - loopMul(cos(x[i] / sqrt((Expr<T>)i + 1.0)), i) + 1;
	}

	template <class T>
	Expr<T> Hansen()
	{
		Expr<T> x;
		Iterator i(0, 4);
		Expr<T> t = (Expr<T>)i;
		return loopSum((t + 1) * cos(t * x[0] + t + 1), i) * loopSum((t + 1) * cos((t + 2) * x[1] + t + 1), i);
	}

	template <class T>
	Expr<T> Hartman3()
	{
		Expr<T> x;
		std::vector<std::vector<double>> a = { { 3, 10, 30 } , { 0.1, 10, 35 }, { 3, 10, 30 }, { 0.1, 10, 35 } };
		std::vector<std::vector<double>> p = { { 0.3689, 0.1170, 0.2673 } ,{ 0.4699, 0.4387, 0.7470 },{ 0.1091, 0.8732, 0.5547 }, { 0.03815, 0.5743, 0.8828 } };
		std::vector<double> c = {1.0, 1.2, 3.0, 3.2};
		
		Expr<T> y = 0.0;
		for (int i = 0; i < 4; i++)
		{
			Expr<T> e = 0.0;
			for (int j = 0; j < 3; j++)
				e += a[i][j] * sqr(x[j] - p[i][j]);
			y += c[i] * exp(-e);
		}
		return -y;
	}

	template <class T>
	Expr<T> Hartman6()
	{
		Expr<T> x;
		std::vector<std::vector<double>> a = { { 10, 3, 17, 3.5, 1.7, 8 } ,
		                                       { 0.05, 10, 17, 0.1, 8, 14 },
		                                       { 3, 3.5, 1.7, 10, 17, 8 }, 
						       { 17, 8, 0.05, 10, 0.1, 14 } };
		std::vector<std::vector<double>> p = { { 0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886 },
		                                       { 0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991 },
						       { 0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650 },
						       { 0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381 } };
		std::vector<double> c = { 1.0, 1.2, 3.0, 3.2 };

		Expr<T> y = 0.0;
		for (int i = 0; i < 4; i++)
		{
			Expr<T> e = 0.0;
			for (int j = 0; j < 6; j++)
				e += a[i][j] * sqr(x[j] - p[i][j]);
			y += c[i] * exp(-e);
		}
		return -y;
	}


	template <class T>
	Expr<T> HelicalValley()
	{
		Expr<T> x;
		Expr<T> t = ifThen( x[0]>=0.0, (0.5/M_PI)*atg(x[1]/x[0]), (0.5 / M_PI)*atg(x[1]/x[0] + 0.5));
		Expr<T> a = sqr(x[2] - 10 * t);
		Expr<T> b = sqr(sqrt(sqr(x[0] + x[1])) - 1);
		return 100 * (a + b) + sqr(x[2]);
	}

	template <class T>
	Expr<T> Himmelblau()
	{
		Expr<T> x;
		return sqr(sqr(x[0]) + x[1] - 11) + sqr(x[0] + sqr(x[1]) - 7);
	}

	template <class T>
	Expr<T> Hosaki()
	{
		Expr<T> x;
		return (1 - 8 * x[0] + 7 * sqr(x[0]) - (7.0 / 3.0)*(x[0] ^ 3) + 0.25*(x[0] ^ 4)) * sqr(x[1]) * exp(-x[1]);
	}

	template <class T>
	Expr<T> JennrichSampson()
	{
		Expr<T> x;
		Iterator i(1, 10);
		Expr<T> t = (Expr<T>)i;
		return loopSum(sqr(2 + 2 * t - (exp(t*x[0]) + exp(t*x[1]))) , i);
	}

	template <class T>
	Expr<T> Langerman5()
	{
		int n = 5;
		Expr<T> x;
		std::vector<std::vector<double>> a = {
			{9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
			{9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
			{8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
			{2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
			{8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567}
		};

		std::vector<double> c = { 0.806, 0.517, 0.100, 0.908, 0.965 };

		Expr<T> y = 0.0;
		for (int i = 0; i < 5; i++)
		{ 
			Expr<T> s = 0.0;
			for (int j = 0; j < n; j++)
				s += sqr(x[j] - a[i][j]);
			y += c[i] * exp((-1.0 / M_PI) * s) * cos(M_PI*s);
		}
		return -y;
	}

	template <class T>
	Expr<T> Keane()
	{
		Expr<T> x;
		return -sqr(sin(x[0] - x[1]))*sqr(sin(x[0] + x[1])) / sqrt(sqr(x[0]) + sqr(x[1]));
	}

	template <class T>
	Expr<T> Leon()
	{
		Expr<T> x;
		return 100 * sqr(x[1] - sqr(x[0])) + sqr(1 - x[0]);
	}

	template <class T>
	Expr<T> Matyas()
	{
		Expr<T> x;
		return 0.26 * (sqr(x[0]) + sqr(x[1])) - 0.48*x[0] * x[1];
	}	

	template <class T>
	Expr<T> McCormick()
	{
		Expr<T> x;
		return sin(x[0] + x[1]) + sqr(x[0] - x[1]) - 1.5*x[0] + 2.5*x[1] + 1;
	}

	template <class T>
	Expr<T> MieleCantrell()
	{
		Expr<T> x;
		return (exp(-x[0]) - x[1]) ^ 4;
	}

	template <class T>
	Expr<T> Mishra3()
	{
		Expr<T> x;
		Expr<T> a = abs(sqr(x[0]) + x[1]);
		Expr<T> b = sqrt(a);
		Expr<T> c = sqrt(abs(cos(b)));
		return c + 0.01*(x[0] + x[1]);
	}

	template <class T>
	Expr<T> Mishra4()
	{
		Expr<T> x;
		Expr<T> a = abs(sqr(x[0]) + x[1]);
		Expr<T> b = sqrt(a);
		Expr<T> c = sqrt(abs(sin(b)));

		return c + 0.01*(x[0] + x[1]);
	}

	template <class T>
	Expr<T> Mishra5()
	{
		Expr<T> x;
		Expr<T> a = sqr(sin(sqr( cos(x[0]) + cos(x[1]) )));
		Expr<T> b = sqr(cos(sqr( sin(x[0]) + sin(x[1]) )));
		return sqr(a + b + x[0]) + 0.01*x[0] + 0.1*x[1];
	}

	template <class T>
	Expr<T> Mishra6()
	{
		Expr<T> x;
		Expr<T> a = sqr(sin(sqr( cos(x[0]) + cos(x[1]) )));
		Expr<T> b = sqr(cos(sqr(sin(x[0]) + sin(x[1]))));
	 	Expr<T> c = 0.1*(sqr(x[0] - 1.0) + sqr(x[1] - 1.0));
		return -ln(sqr(a - b + x[0])) + c;
	}

	template <class T>
	Expr<T> Mishra7()
	{
		int n = 2;
		Expr<T> x;
		Iterator i(0, n - 1);
		int N = 1;
		for (int i = 1; i <= n; i++)
			N *= i;
		return sqr(loopMul(x[i], i) - N);
	}

	template <class T>
	Expr<T> Mishra8()
	{
		Expr<T> x;
		Expr<T> a = pow(x[0], 10) - 20*pow(x[0], 9) + 180 * pow(x[0], 8) - 960*pow(x[0], 7) + 3360*pow(x[0], 6);
		Expr<T> b = -8064 * pow(x[0], 5) + 1334 * pow(x[0], 4) - 15360 * pow(x[0], 3) + 11520 * pow(x[0], 2) - 5120*x[0] + 2624;
		Expr<T> c = pow(x[1], 4) + 12 * pow(x[1], 3) + 54 * pow(x[1], 2) + 108 * x[1] + 81;
		return 0.001*sqr(abs(a + b)*abs(c));
	}

	template <class T>
	Expr<T> Mishra9()
	{
		Expr<T> x;
		Expr<T> a = 2 * pow(x[0], 3) + 5 * x[0] * x[1] + 4 * x[2] - 2 * sqr(x[0])*x[2] - 18;
		Expr<T> b = x[0] + pow(x[1], 3) + x[0] * sqr(x[1]) + x[0] * sqr(x[2]) - 22;
		Expr<T> c = 8 * sqr(x[0]) + 2 * x[1] * x[2] + 2 * sqr(x[1]) + 3 * pow(x[1], 3) - 52;
		return sqr(a*sqr(b)*c + a*b*sqr(c) + sqr(b) + sqr(x[0] + x[1] - x[2]));

	}

	template <class T>
	Expr<T> Parsopoulos()
	{
		Expr<T> x;
		return sqr(cos(x[0])) + sqr(sin(x[1]));
	}

	template <class T>
	Expr<T> Pathological(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 2);
		Expr<T> t = i;
		Expr<T> a = sqr(sin(sqrt(100 * sqr(x[i]) + sqr(x[t+1])))) - 0.5;
		Expr<T> b = 1 + 0.001*sqr(sqr(x[i]) - 2 * x[i] * x[t+1] + sqr(x[t+1]));
		return loopSum(0.5 + a / b, i);
	}
	
	template <class T>
	Expr<T> Pinter(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
                Expr<T> t1 = (Expr<T>)i - 1.0;	
                Expr<T> t2 = (Expr<T>)i + 1.0;	
                Expr<T> pred = ifThen(t1 < 0, Expr<T>(0.0), t1);
		Expr<T> sus = ifThen(t2 >= n, Expr<T>(n-1), t2);
		auto A = x[pred] * sin(x[i]) + sin(x[sus]);
		auto B = sqr(x[pred]) - 2 * x[i] + 3 * x[sus] - cos(x[i]) + 1;
		auto t = (Expr<T>)i;
		return loopSum(t*sqr(x[i]), i) + loopSum(20 * t*sqr(sin(A)), i) + loopSum(t*log(1 + t*sqr(B), 10), i);	
	}

	template <class T>
	Expr<T> Periodic()
	{
		Expr<T> x;
		return 1 + sqr(sin(x[0])) + sqr(sin(x[1])) - 0.1*exp(-(sqr(x[0]) + sqr(x[1])));
	}
	
	template <class T>
	Expr<T> PowellSingular2(int n)
	{
		Expr<T> x;
		Iterator i(1, n / 4);
		Expr<T> t = i;
		return loopSum(sqr(x[4*t - 4] + 10 * x[4*t - 3]) + 5 * sqr(x[4*t - 2] - x[4*t - 1]) + pow(x[4*t-3] - x[4*t - 2], 4) + 10 * pow(x[4*t - 4] - x[4*t - 1], 4), i);
	}

	template <class T>
	Expr<T> PowellSum(int n)
	{
		Expr<T> x;
		Iterator i(0, n-1 );
		Expr<T> t = i;
		return loopSum(pow(abs(x[i]), t+2), i);
	}

	template <class T>
	Expr<T> Price1()
	{
		Expr<T> x;
		return sqr(abs(x[0]) - 5) + sqr(abs(x[1]) - 5);
	}

	template <class T>
	Expr<T> Price2()
	{
		Expr<T> x;
		return 1 + sqr(sin(x[0])) + sqr(sin(x[1])) - 0.1*exp(-sqr(x[0]) - sqr(x[1]));
	}

	template <class T>
	Expr<T> Price3()
	{
		Expr<T> x;
		return 100 * sqr(x[1] - sqr(x[0])) + sqr(6.4*sqr(x[1]-0.5) - x[0] - 0.6);
	}

	template <class T>
	Expr<T> Price4()
	{
		Expr<T> x;
		return sqr(2 * pow(x[0], 3) * x[1] - pow(x[1], 3)) + sqr(6*x[0] -sqr(x[1]) + x[1]);
	}

	template <class T>
	Expr<T> Problem02()
	{
		Expr<T> x;
		return sin(x[0])+sin(x[0]*10.0/3.0);
	}

	template <class T>
	Expr<T> Problem04()
	{
		Expr<T> x;
		return -(16.0*sqr(x[0])-24*x[0]+5)*exp(-x[0]);
	}

	template <class T>
	Expr<T> Problem05()
	{
		Expr<T> x;
		return -(1.4-3.0*x[0])*sin(18.0*x[0]);
	}

	template <class T>
	Expr<T> Problem06()
	{
		Expr<T> x;
		return -(x[0]+sin(x[0]))*exp(-sqr(x[0]));
	}

	template <class T>
	Expr<T> Qing()
	{
		Expr<T> x;
		Iterator i(0, 1);
		Expr<T> t = i;
		return loopSum(sqr(sqr(x[i]) - t - 1.0), i);
	}

	template <class T>
	Expr<T> Quadratic()
	{
		Expr<T> x;
		return -3803.84 - 138.08*x[0] - 232.92*x[1] + 128.08*sqr(x[0]) + 203.64*sqr(x[1]) + 182.25 * x[0] * x[1];
	}

	template <class T>
	Expr<T> Quintic(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return loopSum(abs(pow(x[i], 5) - 3*pow(x[i], 4) + 4*pow(x[i], 3) + 2*sqr(x[i]) - 10*x[i] - 4), i);
	}

	template <class T>
	Expr<T> Rosenbrock(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 2);
		Expr<T> t = i;
		return loopSum(100 * sqr(x[t+1] -
		sqr(x[i])) + sqr((x[i]-1)), i);
	}

	template <class T>
	Expr<T> RosenbrockModified()
	{
		Expr<T> x;
		return 74 + 100 * sqr(x[1] - sqr(x[0])) + sqr((1 - x[0])) - 400*exp(-(sqr(x[0]+1)+sqr(x[1]+1))/0.1);
	}

	template <class T>
	Expr<T> RotatedEllipse()
	{
		Expr<T> x;
		return 7 * sqr(x[0]) - 6 * ::sqrt(3)*x[0] * x[1] + 13 * sqr(x[1]);
	}

	template <class T>
	Expr<T> RotatedEllipse2()
	{
		Expr<T> x;
		return sqr(x[0]) - x[0] * x[1] + sqr(x[1]);
	}

	template <class T>
	Expr<T> Solomon()
	{
		int n = 2;
		Expr<T> x;
		Iterator i(0, n - 1);
		return 1.0 - cos(2 * M_PI * sqrt(loopSum(sqr(x[i]), i))) + 0.1*sqrt(loopSum(sqr(x[i]), i));
	}

	template <class T>
	Expr<T> Scahffer1()
	{
		Expr<T> x;
		Expr<T> t = sqr(sqr(x[0]) + sqr(x[1]));
		return 0.5 + (sqr(sin(t)) - 0.5 ) / (1 + 0.001*t);
	}

	template <class T>
	Expr<T> Scahffer3()
	{
		Expr<T> x;
		Expr<T> t = sqr(x[0]) + sqr(x[1]);
		return 0.5 + (sqr(sin(cos(abs(t)))) - 0.5) / (1 + 0.001*sqr(t));
	}

	template <class T>
	Expr<T> Scahffer4()
	{
		Expr<T> x;
		Expr<T> t = sqr(x[0]) + sqr(x[1]);
		return 0.5 + (sqr(cos(sin(abs(t)))) - 0.5) / (1 + 0.001*sqr(t));
	}

	template <class T>
	Expr<T> Scahffer2_6()
	{
		Expr<T> x;
		return max(abs(x[0] + 2*x[1]-7),abs(2*x[0] + x[1]-5));
	}

	template <class T>
	Expr<T> SchmidtVetters()
	{
		Expr<T> x;
		Expr<T> a = 1.0 / (1 + sqr(x[0] - x[1]));
		Expr<T> b = sin((M_PI*x[1] + x[2]) / 2.0);
		Expr<T> c = exp(sqr(((x[0] + x[1]) / x[1]) - 2.0));
		return a + b + c;
	}

	template <class T>
	Expr<T> SchumerSteiglitz(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return loopSum(pow(x[i],4),i);
	}

	template <class T>             
	Expr<T> Schwefel(int n, double e)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return pow(loopSum(sqr(x[i]), i), e);
	}

	template <class T>
	Expr<T> Schwefel1_2(int n)
	{
		Expr<T> x;
		Expr<T> y = 0.0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j <= i; j++)
				y += sqr(x[j]);
		return y;
	}

	template <class T>
	Expr<T> Schwefel2_4(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return loopSum(sqr(x[i] - 1.0), i) + sqr(x[0] - sqr(x[1]));
	}

	template <class T>
	Expr<T> Schwefel2_6()
	{
		Expr<T> x;
		return max(abs(x[0]+2*x[1]-7), abs(2*x[0] + x[1] - 5));
	}

	template <class T>
	Expr<T> Schwefel2_20(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return loopSum(abs(x[i]), i);
	}

	template <class T>
	Expr<T> Schwefel2_22(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return loopSum(abs(x[i]), i) + loopMul(abs(x[i]), i);
	}

	template <class T>
	Expr<T> Schwefel2_23(int n)
	{
		Expr<T> x;
		Iterator i(0, n-1);
		return loopSum(pow(x[i], 10), i);
	}

	template <class T>
	Expr<T> Schwefel2_26()
	{
		Expr<T> x;
		Iterator i(0, 1);
		return -loopSum(x[i]*sin(sqrt(abs(x[i]))), i);
	}

	template <class T>
	Expr<T> Schwefel2_36()
	{
		Expr<T> x;
		return -x[0] * x[1] * (72 - 2 * x[0] - 2 * x[1]);
	}

	template <class T>
	Expr<T> Shekel5()
	{
		Expr<T> x;       
		std::vector<std::vector<double>> a = { { 4,4,4,4 },{ 1,1,1,1 },{ 8,8,8,8 },{ 6,6,6,6 },{3,7,3,7} };
		std::vector<double> c = { 0.1, 0.2, 0.2, 0.4, 0.4 };

		Expr<T> y = 0.0;
		for (int i = 0; i < 5; i++)
		{
			Expr<T> t = 0.0;
			for (int j = 0; j < 4; j++)
				t += sqr(x[j] - a[i][j]);
			y += 1.0 / (t + c[i]);
		}
		return -y;
	}

	template <class T>
	Expr<T> Shekel7()
	{
		Expr<T> x;
		std::vector<std::vector<double>> a = { { 4.0,4.0,4.0,4.0 },{ 1.0,1.0,1.0,1.0 },{ 8.0,8.0,8.0,8.0 },{ 6.0,6.0,6.0,6.0 },{ 3.0,7.0,3.0,7.0 }, {2.0,9.0,2.0,9.0}, {5.0,5.0,3.0,3.0} };
		std::vector<double> c = { 0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3 };

		Expr<T> y = 0.0;
		for (int i = 0; i < 7; i++)
		{
			Expr<T> t = 0.0;
			for (int j = 0; j < 4; j++)
				t += sqr(x[j] - a[i][j]);
			y += 1.0 / (t + c[i]);
		}
		return -y;
	}

	template <class T>
	Expr<T> Shekel10()
	{
		Expr<T> x;
		std::vector<std::vector<double>> a = { { 4,4,4,4 },{ 1,1,1,1 },{ 8,8,8,8 },{ 6,6,6,6 },{ 3,7,3,7 },{ 2,9,2,9 },{ 5,5,3,3 },{8,1,8,1},{6,2,6,2},{7,3.6,7,3.6} };
		std::vector<double> c = { 0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 07,0.5,0.5 }; 

		Expr<T> y = 0.0;
		for (int i = 0; i < 10; i++)
		{
			Expr<T> t = 0.0;
			for (int j = 0; j < 4; j++)
				t += sqr(x[j] - a[i][j]);
			y += 1.0 / (t + c[i]);
		}
		return -y;
	}

	template <class T>
	Expr<T> Shubert()
	{
		Expr<T> x;
		return -sin(2.0*x[0]+1.0)
		-2.0*sin(3.0*x[0]+2.0)
		-3.0*sin(4.0*x[0]+3.0)
		-4.0*sin(5.0*x[0]+4.0)
		-5.0*sin(6.0*x[0]+5.0)
		-sin(2.0*x[1]+1.0)
		-2.0*sin(3.0*x[1]+2.0)
		-3.0*sin(4.0*x[1]+3.0)
		-4.0*sin(5.0*x[1]+4.0)
		-5.0*sin(6.0*x[1]+5.0);
	}

	template <class T>
	Expr<T> Shubert2()
	{
		Expr<T> x;
		Iterator i(1, 5);
		Expr<T> t = (Expr<T>)i;
		Expr<T> s1 = loopSum(t*cos((t + 1.0)*x[0] + t), i);
                Expr<T> s2 = loopSum(t*cos((t + 1.0)*x[1] + t), i);
                return s1 * s2;
	}

	template <class T>
	Expr<T> Shubert3()
	{
		Expr<T> x;
		Iterator i(0, 1), j(0, 4);
		Expr<T> t = (Expr<T>)j;
		return -loopSum(loopSum((t+1)*sin((t + 2.0)*x[i] + t + 1.0), j), i);
	}



	template <class T>
	Expr<T> SchafferF6(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 2);
		Expr<T> t = (Expr<T>)i;
		Expr<T> y = sqr(x[i]) + sqr(x[t + 1]);
		return loopSum(0.5 + (sqr(sin(y)) - 0.5) / sqr(1 + 0.001*y), i);
	}

	template <class T>
	Expr<T> Sphere(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return loopSum(sqr(x[i]), i);
	}

	template <class T>
	Expr<T> StrechedVSineWave(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 2);
		Expr<T> t = i;
		Expr<T> y = sqr(x[t + 1] + sqr(x[i]));
		return loopSum(pow(y, (Expr<T>)0.25)*(sqr(sin(50*pow(y, (Expr<T>)0.1))) + 0.1), i);
	}

	template <class T>
	Expr<T> SumSquares(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		Expr<T> t = (Expr<T>)i + 1;
		return loopSum(t*sqr(x[i]), i);
	}

	template <class T>
	Expr<T> StyblinskiTang()
	{
		int n = 2;
		Expr<T> x;
		Iterator i(0, n - 1);
		return 0.5*loopSum(pow(x[i], 4) - 16 * sqr(x[i]) + 5 * x[i], i);
	}

	template <class T>
	Expr<T> Table1HolderTable1()
	{
		Expr<T> x;
		Expr<T> t = abs(1-sqrt(sqr(x[0])+sqr(x[1]))/M_PI);
		return -abs(cos(x[0])*cos(x[1])*exp(t));
	}

	template <class T>
	Expr<T> Table2HolderTable2()
	{
		Expr<T> x;
		Expr<T> t = abs(1.0 - sqrt(sqr(x[0]) + sqr(x[1])) / M_PI);
		return -abs(sin(x[0])*cos(x[1])*exp(t));
	}

	template <class T>
	Expr<T> Table3Carrom()
	{
		Expr<T> x;
		Expr<T> t = abs(1 - sqrt(sqr(x[0]) + sqr(x[1])) / M_PI);
		return -sqr(cos(x[0])*cos(x[1])*exp(t))/30.0;
	}

	template <class T>
	Expr<T> TesttubeHolder()
	{
		Expr<T> x;
		Expr<T> t = abs(cos((sqr(x[0]) + sqr(x[1])) / 200.0));
		return -4.0*abs(sin(x[0])*cos(x[1])*exp(t));
	}

	template <class T>
	Expr<T> Trecanni()
	{
		Expr<T> x;
		return pow(x[0], 4) + 4 * pow(x[0], 3) + 4 * sqr(x[0]) + sqr(x[1]);
	}

	template <class T>
	Expr<T> Trid6()
	{
		int n = 6;
		Expr<T> x;
		Iterator i(0, n - 1), i2(1, n - 1);
		Expr<T> t = (Expr<T>)i2 - 1;
		return loopSum(sqr(x[i] - 1), i) - loopSum(x[i2] * x[t], i2);
	}

	template <class T>
	Expr<T> Trid10()
	{
		int n = 10;
		Expr<T> x;
		Iterator i(0, n - 1), i2(1, n - 1);
		Expr<T> t = (Expr<T>)i2 - 1;
		return loopSum(sqr(x[i] - 1), i) - loopSum(x[i2] * x[t], i2);
	}

	template <class T>
	Expr<T> Tripod()
	{
		Expr<T> x;
		Expr<T> p0 = ifThen(x[0] >=0, Expr<T>(1.0), Expr<T>(0.0));
		Expr<T> p1 = ifThen(x[1] >=0, Expr<T>(1.0), Expr<T>(0.0));
		return p1*(1.0 + p0) + abs(x[0] + 50.0*p1*(1.0-2.0*p0)) + abs(x[1] + 50.0*(1.0-2.0*p1));
	}

	template <class T>
	Expr<T> Trefethen()
	{
		Expr<T> x;
		Expr<T> a = exp(sin(50 * x[0])) + sin(60 * exp(x[1]));
		Expr<T> b = sin(70 * sin(x[0])) + sin(sin(80 * x[1]));
		Expr<T> c = -sin(10 * (x[0] + x[1])) + 0.25*(sqr(x[0]) + sqr(x[1]));
		return a + b + c;
	}

	template <class T>
	Expr<T> Trigonometric1(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1), j(0, n - 1);
		Expr<T> t = (Expr<T>)i + 1;
		return loopSum(sqr(n - loopSum(cos(x[j]) + t*(1 - cos(x[i]) - sin(x[i])), j)), i);
	}

	template <class T>
	Expr<T> Trigonometric2(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		Expr<T> a = 8 * sqr(sin(7 * sqr(x[i] - 0.9)));
		Expr<T> b = 6 * sqr(sin(14 * sqr(x[i] - 0.9)));
		Expr<T> c = sqr(x[i] - 0.9);
		return 1.0 + loopSum(a + b + c, i);
	}

	template <class T>
	Expr<T> Ursem1()
	{
		Expr<T> x;
		return -sin(2 * x[0] - M_PI_2) - 3 * cos(x[1]) - 0.5*x[0];
	}

	template <class T>
	Expr<T> Ursem3()
	{
		Expr<T> x;
		Expr<T> a = -sin(2.2*M_PI*x[0] + M_PI_2);
		Expr<T> b = -sin(M_PI_2*sqr(x[1]) + M_PI_2);
		Expr<T> c = (2 - abs(x[1])) / 2;
		Expr<T> d = (3 - abs(x[0])) / 2;
		Expr<T> e = (2 - abs(x[0])) / 2;

		return a*c*d + b*c*e;
	}

	template <class T>
	Expr<T> Ursem4()
	{
		Expr<T> x;
		Expr<T> a = -3*sin(2.2*M_PI*x[0] + M_PI_2);
		Expr<T> b = (2 - sqrt(sqr(x[0])+ sqr(x[1]))) / 4;

		return a*b;
	}

	template <class T>
	Expr<T> UrsemWaves()
	{
		Expr<T> x;
		Expr<T> a = -pow(0.3*x[0], 3) + (sqr(x[1]) - 4.5*sqr(x[1]))*x[0] * x[1];
		Expr<T> b = 4.7*cos(3 * x[0] - sqr(x[1])*(2 + x[0]))*sin(2.5*M_PI*x[0]);
		return a + b;
	}

	template <class T>
	Expr<T> VenterSobiezcczanskiSobieski()
	{
		Expr<T> x;
		Expr<T> a = sqr(x[0]) - 100 * sqr(cos(x[0]));
		Expr<T> b = -100 * cos(sqr(x[0]) / 30) + sqr(x[1]);
		Expr<T> c = -100 * sqr(cos(x[1])) - 100 * cos(sqr(x[1]) / 30);
		return a + b + c;
	}

	template <class T>
	Expr<T> Watson()
	{
		Expr<T> x;
		Iterator i(0, 29), j(0, 4), k(0, 5);

		Expr<T> a = (Expr<T>)i / 29.0;
		Expr<T> y = ((Expr<T>)j - 1)*pow(a, (Expr<T>)j)*x[j];
		Expr<T> z = pow(a, (Expr<T>)k)*x[k];

		return loopSum(sqr( loopSum(y, j) - sqr(loopSum(z, k)) -1 ), i) + sqr(x[0]);

	}

	template <class T>
	Expr<T> WayburnSeader1()
	{
		Expr<T> x;
		return sqr(pow(x[0], 6) + pow(x[1], 4) - 17) + sqr(2 * x[0] + x[1] - 4);
	}

	template <class T>
	Expr<T> WayburnSeader2()
	{
		Expr<T> x;
		return sqr(1.613 - 4 * sqr(x[0] - 0.3125) - 4 * sqr(x[1] - 1.625)) + sqr(x[1] - 1);
	}

	template <class T>
	Expr<T> WayburnSeader3()
	{
		Expr<T> x;
		return (2.0 / 3.0)*pow(x[0], 3) - 8.0 * sqr(x[0]) + 33.0 * x[0] - x[0] * x[1] + 5.0 + sqr(sqr(x[0] - 4.0) + sqr(x[1] - 5.0) - 4.0);
	}

	template <class T>
	Expr<T> WWavy(int n)
	{
		Expr<T> x;
		int k = 10;
		Iterator i(0, n - 1);
		return 1 - (1.0 / n) * loopSum(cos(k*x[i]) * exp(-sqr(x[i] / 2)), i);
	}

	template <class T>
	Expr<T> Weierstrass(int n)
	{
		Expr<T> x;
		Expr<T> a = 0.5, b = 3.0;
		int kmax = 20;
		Iterator i(0, n - 1), j(0, kmax);
		Expr<T> k = (Expr<T>)j;

		return loopSum(loopSum(pow(a, k)*cos(2 * M_PI*pow(b, k)*(x[i] + 0.5)), j), i) - n * loopSum(pow(a, k)*cos(M_PI*pow(b, k)), j);
	}

	template <class T>
	Expr<T> Whitley(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1), j(0, n-1);
		Expr<T> a = 100*sqr(sqr(x[i]-x[j])) + sqr(1-x[j]);
		return loopSum(loopSum(sqr(a) / 4000.0 - cos(a) + 1.0, j), i);
	}

	template <class T>
	Expr<T> Wolfe()
	{
		Expr<T> x;
		return (4.0 / 3.0)*(sqr(x[0]) + sqr(x[1]) - x[0] * x[1])*0.75 + x[2];
	}

	template <class T>
	Expr<T> XinSheYang2(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		return loopSum(abs(x[i]), i)*exp(1 - loopSum(sin(sqr(x[i])), i));
	}

	template <class T>
	Expr<T> XinSheYang3(int n)
	{
		Expr<T> x, b=15;
		Iterator i(0, n - 1);
		return exp(-loopSum(pow(x[i] / b, 2 * 5), i)) - 2 * exp(-loopSum(sqr(x[i]), i)) * loopMul(sqr(cos(x[i])), i);
	}

	template <class T>
	Expr<T> XinSheYang4(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1), j(0, n - 1), k(0, n-1);
		Expr<T> a = loopSum(sqr(sin(x[i])), i);
		Expr<T> b = exp(-loopSum(sqr(x[j]), j));
		Expr<T> c = exp(-loopSum(sqr(sin(sqrt(abs(x[k])))), k));
		return (a-b)*c;
	}
	
	template <class T>
	Expr<T> Zakharov(int n)
	{
		Expr<T> x;
		Iterator i(0, n - 1);
		Expr<T> k = (Expr<T>)i;
		return loopSum(sqr(x[i]),i) + sqr(0.5*loopSum(k*x[i], i)) + pow(0.5*loopSum(k*x[i], i), 4);
	}


	template <class T>
	Expr<T> Zettl()
	{
		Expr<T> x;
		return sqr(sqr(x[0]) + sqr(x[1]) - 2 * x[0]) + 0.25*x[0];
	}

	template <class T>
	Expr<T> Zirilli()
	{
		Expr<T> x;
		return 0.25*pow(x[0], 4) - 0.5*sqr(x[0]) + 0.1*x[0] + 0.5*sqr(x[1]);
	}

	template <class T>
	Expr<T> Factorial(int n)
	{
		if (n == 0.0)
			return 1.0;
		Iterator i(1, n);
		Expr<T> result = loopMul((Expr<T>)i, i);
		return result;
	}

}
#endif
