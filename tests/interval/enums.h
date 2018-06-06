#ifndef ENUMS__HPP
#define ENUMS__HPP

enum IntervalBool
{
	True,
	False,
	Intermadiate
};

enum Conditions
{
	More,
	Less,
	MoreEqual,
	LessEqual
};

std::ostream& operator<<(std::ostream & out, Conditions condition) 
{
	switch (condition)
	{
	case Conditions::More:
		return out << ">";
	case Conditions::Less:
		return out << "<";
	case Conditions::LessEqual:
		return out << "<=";
	case Conditions::MoreEqual:
		return out << ">=";
	}
	return out;
}

#endif
