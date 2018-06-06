#ifndef NODE__HPP
#define NODE__HPP

#include <iostream> 
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>


#include "algorithm.hpp"
#include "interval/interval_air.hpp"
#include "interval/enums.h"
#include "expr.hpp"
#include "utils.h"
//#include "mapiterator.hpp"




namespace snowgoose {
namespace expression {

	template<class T> class Node;
	template<class T> class ConditionNode;
	template<typename T> using ptrNode = std::shared_ptr<Node<T>>;
	template<typename T> using ptrCNode = std::shared_ptr<ConditionNode<T>>;
	template<typename T> using vPtrNode = std::vector<ptrNode<T>>;
	class MapIterator;

	template<class T>
	class Node
	{
	public:
		Node(const vPtrNode<T> &childs) : m_childs(childs) {}
		Node() {}
		virtual T calc(const Algorithm<T> &, MapIterator &) const = 0;
		friend std::ostream& operator<<(std::ostream & out, const Node<T>& v) { return v.prn(out);}
		virtual std::ostream& prn(std::ostream & out) const = 0;
        virtual bool IsVar() const {return false;}
        bool IsVarInTree() const
        {
            if(IsVar()) 
                return true;
            else
                for(auto &node : m_childs)
                    if(node->IsVarInTree())
                        return true;
            return false;
        }
	protected:
		vPtrNode<T> m_childs;
	};	

	template<class T>
	class Var : public Node<T>
	{
	private:
		const int index;
	public:
		Var(int i) : index(i) {}
		T calc(const Algorithm<T> &alg, MapIterator &map_iterator) const { return alg.CreateVar(index); };
		std::ostream& prn(std::ostream & out) const { return out << "x[" << index << "]"; };
        	bool IsVar() const { return true; }
		int getIndex() { return index; }
	};


	template <class T>
	class Const : public Node<T>
	{
	private:
		const double m_const;
	public:
		Const(double value) : m_const(value) {}
		T calc(const Algorithm<T> &alg, MapIterator &map_iterator) const { return alg.CreateConst(m_const); }
		std::ostream& prn(std::ostream & out) const { return out << m_const; }
		double getConst() { return m_const; }
	};

	template <class T>
	class Plus : public Node<T>
	{
	public:
		Plus(const ptrNode<T> &left, const ptrNode<T> &right) : Node<T>({ left, right }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Plus(this->m_childs[0]->calc(alg, map_iterator), this->m_childs[1]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "(" << *this->m_childs[0] << " + " << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Minus : public Node<T>
	{
	public:
		Minus(const ptrNode<T> &left, const ptrNode<T> &right) : Node<T>({ left, right }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Minus(this->m_childs[0]->calc(alg, map_iterator), this->m_childs[1]->calc(alg, map_iterator)); }
		std::ostream& prn(std::ostream & out) const { return out << "(" << *this->m_childs[0] << " - " << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Mul : public Node<T>
	{
	public:
		Mul(const ptrNode<T> &left,const ptrNode<T> &right) : Node<T>({ left, right }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Mul(this->m_childs[0]->calc(alg, map_iterator), this->m_childs[1]->calc(alg, map_iterator)); }
		std::ostream& prn(std::ostream & out) const { return out << *this->m_childs[0] << " * " << *this->m_childs[1]; }
	};

	template <class T>
	class Div : public Node<T>
	{
	public:
		Div(const ptrNode<T> &left, const ptrNode<T> &right) : Node<T>({ left, right }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Div(this->m_childs[0]->calc(alg, map_iterator), this->m_childs[1]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "(" << *this->m_childs[0] << "/" << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Sin : public Node<T>
	{
	public:
		Sin(const ptrNode<T> &node) : Node<T>({node}) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Sin(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "sin(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Cos : public Node<T>
	{
	public:
		Cos(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Cos(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "cos(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Tg : public Node<T>
	{
	public:
		Tg(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Tan(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "tg(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Ctg : public Node<T>
	{
	public:
		Ctg(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Ctg(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "ctg(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class ArcCos : public Node<T>
	{
	public:
		ArcCos(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.ArcCos(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "acos(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class ArcSin : public Node<T>
	{
	public:
		ArcSin(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.ArcSin(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "asin(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class ArcTg : public Node<T>
	{
	public:
		ArcTg(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.ArcTan(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "atg(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class ArcCtg : public Node<T>
	{
	public:
		ArcCtg(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.ArcCtg(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "actg(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Exp : public Node<T>
	{
	public:
		Exp(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Exp(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "exp(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Sqrt : public Node<T>
	{
	public:
		Sqrt(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Sqrt(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "sqrt(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Sqr : public Node<T>
	{
	public:
		Sqr(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Sqr(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "sqr(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class PowInt : public Node<T>
	{
	private:
		const int exponent;
	public:
		PowInt(const ptrNode<T> &node, int exp) : Node<T>({ node }), exponent(exp) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { 
			return alg.Pow(this->m_childs[0]->calc(alg, map_iterator), exponent); 
		};
		std::ostream& prn(std::ostream & out) const { return out << "pow(" << *this->m_childs[0] << "," << exponent << ")"; }
	};

	template <class T>
	class Pow : public Node<T>
	{
	private:
		const double exponent;
	public:
		Pow(const ptrNode<T> &node, double exp) : Node<T>({ node }), exponent(exp) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const {
			return alg.PowDouble(this->m_childs[0]->calc(alg, map_iterator), exponent);
		};
		std::ostream& prn(std::ostream & out) const { return out << "pow(" << *this->m_childs[0] << "," << exponent << ")"; }
	};

	template <class T>
	class PowExpr : public Node<T>
	{
	public:
		PowExpr(const ptrNode<T> &base, const ptrNode<T> &exp) : Node<T>({ base, exp }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { 
            return alg.Pow(this->m_childs[0]->calc(alg, map_iterator), 
                           this->m_childs[0]->IsVarInTree(),
                           this->m_childs[1]->calc(alg, map_iterator),
                           this->m_childs[1]->IsVarInTree()); };
		std::ostream& prn(std::ostream & out) const { return out << "pow(" << *this->m_childs[0] << "," << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Abs : public Node<T>
	{
	public:
		Abs(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Abs(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "abs(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Ln : public Node<T>
	{
	public:
		Ln(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Ln(this->m_childs[0]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "ln(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Log : public Node<T>
	{
	private:
		const double base;
	public:
		Log(const ptrNode<T> &node, double b) : Node<T>({ node }), base(b) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Log(this->m_childs[0]->calc(alg, map_iterator), base); };
		std::ostream& prn(std::ostream & out) const { return out << "log(" << *this->m_childs[0] << "," << base << ")"; }
	};

	template <class T>
	class Min : public Node<T>
	{
	public:
		Min(const ptrNode<T> &lv, const ptrNode<T> &rv) : Node<T>({ lv,  rv }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Min(this->m_childs[0]->calc(alg, map_iterator), this->m_childs[1]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "min(" << *this->m_childs[0] << "," << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Max : public Node<T>
	{
	public:
		Max(const ptrNode<T> &lv, const ptrNode<T> &rv) : Node<T>({ lv,  rv }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const { return alg.Max(this->m_childs[0]->calc(alg, map_iterator), this->m_childs[1]->calc(alg, map_iterator)); };
		std::ostream& prn(std::ostream & out) const { return out << "max(" << *this->m_childs[0] << "," << *this->m_childs[1] << ")"; }
	};

	template <class T> class IfTrue;

	template <class T>
	class ConditionNode : public Node<T>
	{
	private:
		const Conditions condition;
	public:
		ConditionNode(Conditions cond, const ptrNode<T> &lv, const ptrNode<T> &rv) : condition(cond), Node<T>({ lv,  rv }) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const 
		{ 
			throw "Invalid operation in ConditionNode.";
		}
		IntervalBool calcCondition(const Algorithm<T> & alg, MapIterator &map_iterator) const
		{
			return alg.Condition(condition, this->m_childs[0]->calc(alg, map_iterator), this->m_childs[1]->calc(alg, map_iterator));
		}
		std::ostream& prn(std::ostream & out) const { return out << " " << *this->m_childs[0] << condition << *this->m_childs[1] << " "; }
		Conditions getCondition() { return condition; }
		friend class IfTrue<T>;
	};

	template <class T>
	class IfTrue : public Node<T>
	{
	private:
		const ptrCNode<T> conditionNode;
	public:
		IfTrue(const ptrCNode<T> &cond, const ptrNode<T> &lv, const ptrNode<T> &rv) : conditionNode(cond), Node<T>({ lv,  rv }) {}

		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const
		{
			IntervalBool condition = conditionNode->calcCondition(alg, map_iterator);
			auto pVar = std::dynamic_pointer_cast<Var<T>>(conditionNode->m_childs[0]);
			auto pCnst = std::dynamic_pointer_cast<Const<T>>(conditionNode->m_childs[1]);
			if(condition == IntervalBool::Intermadiate && pVar && pCnst)
			{
				Conditions cond = conditionNode->getCondition();
				int index = pVar->getIndex();
				double cnst = pCnst->getConst();
				auto vPtrAlg = alg.GetNewAlgorithm(cond, index, cnst);
				return alg.IfTrue(condition, this->m_childs[0]->calc(*vPtrAlg[0], map_iterator), this->m_childs[1]->calc(*vPtrAlg[1], map_iterator));
			}
			else
				return alg.IfTrue(condition, this->m_childs[0]->calc(alg, map_iterator), this->m_childs[1]->calc(alg, map_iterator));
			
		}
		std::ostream& prn(std::ostream & out) const { return out << "ifThen(" << *conditionNode << " , " << *this->m_childs[0] << " , " << *this->m_childs[1] << ")"; }
	};

	template <class T> class Expr;
	template <class T> class IteratorNode;

	/**
	* Iterator is used  for loopSum, loopMul expressions.
	*/
	template <class T> class IteratorNode;
	class Iterator
	{
	public:
		Iterator(int start, int end) : startIterator(start), endIterator(end)
		{
			uid = Utils::getUid();
		}
		int Start() const { return startIterator; }
		int End() const { return endIterator; }
		int Uid() const { return uid; }		
		/**
		* Allows to conver iterator to expression
		* @return expression
		*/
		template <class T>
		operator Expr<T>() const
		{
			ptrNode<T> pNode(new IteratorNode<T>(*this));
			return Expr<T>(pNode);
		}
	private:
		int startIterator;
		int endIterator;
		int uid;
	};

	class MapIterator {
	private :
	    	std::unordered_map<int,int> map_iterator;
	public:
		int Current(const Iterator &i) {
			if(map_iterator.find(i.Uid())==map_iterator.end())
				Reset(i);	 
			return map_iterator[i.Uid()]; 
		}
		bool CanIterate(const Iterator &i) { return map_iterator[i.Uid()] <= i.End(); }
		void Next(const Iterator &i) { map_iterator[i.Uid()]++; }
		void Reset(const Iterator &i) { map_iterator[i.Uid()] = i.Start(); }
	};

	template <class T>
	class IteratorNode : public Node<T>
	{
		const Iterator i;
	public:
		IteratorNode(const Iterator &it) : i(it) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const
		{			
			int index = map_iterator.Current(i);
            		return alg.CreateConst(index);
			
		}
		std::ostream& prn(std::ostream & out) const { return out << "i"; }
	};

	template <class T>
	class Index : public Node<T>
	{
		const Iterator i;
	public:
		Index(const Iterator &it) : i(it) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const
		{
			int index = map_iterator.Current(i);
			return alg.CreateVar(index);
		}
        	bool IsVar() const { return true; }
		std::ostream& prn(std::ostream & out) const { return out << "x[i]"; };
	};

	template <class T>
	class CalcIndex : public Node<T>
	{
		const std::function<int()> func;
	public:
		CalcIndex(const std::function<int()> &f) : func(f) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const
		{
			return alg.CreateVar(func());
		}
        	bool IsVar() const { return true; }
		std::ostream& prn(std::ostream & out) const { return out << "x[calc i]"; };
	};

	template <class T>
	class ExprIndex : public Node<T>
	{
		const Expr<T> expr;
	public:
		ExprIndex(const Expr<T>& e) : expr(e) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const
		{
			size_t index = (size_t)((double)expr.node->calc(alg, map_iterator));
			return alg.CreateVar(index);
		}
        	bool IsVar() const { return true; }
		std::ostream& prn(std::ostream & out) const { return out << "x[" << expr << "]"; };
	};

	template <class T>
	class CycleSum : public Node<T>
	{
	private:
		const Iterator i;
	public:
		CycleSum(const ptrNode<T> &node, const Iterator &it) : Node<T>({ node }), i(it) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const {
            		T result = alg.CreateConst(0.0);
			map_iterator.Reset(i);
			while (map_iterator.CanIterate(i))
			{
				result = alg.Plus(result, this->m_childs[0]->calc(alg, map_iterator));
				map_iterator.Next(i);
			}
			return result;
		};
		std::ostream& prn(std::ostream & out) const { return out << "loopSum(" << *this->m_childs[0] << ",i)"; }
	};

	template <class T>
	class CycleMul : public Node<T>
	{
	private:
		const Iterator i;
	public:
		CycleMul(const ptrNode<T> &node, const Iterator &it) : Node<T>({ node }), i(it) {}
		T calc(const Algorithm<T> & alg, MapIterator &map_iterator) const {		
			T result = alg.CreateConst(1.0);
			map_iterator.Reset(i);
			if (!map_iterator.CanIterate(i))
				return alg.CreateConst(0.0);
			while (map_iterator.CanIterate(i))
			{
				result = alg.Mul(result, this->m_childs[0]->calc(alg, map_iterator));
				map_iterator.Next(i);
			}
			return result;
		};
		std::ostream& prn(std::ostream & out) const { return out << "loopMul(" << *this->m_childs[0] << ",i)"; }
	};

}}

#endif /* NODE__HPP */

