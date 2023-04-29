//-------------------------------------------------------------------------------------------
// SimpleFraction.h
// Description: Класс, реализующий простые дроби
// Date: 31.01.2023
// Authors: Хьюго М.А.
// Ⓒ Sibsutis university
//-------------------------------------------------------------------------------------------

#include <iostream>
#include <string>

#pragma once

namespace Fraction {
class SimpleFraction {
public:
	SimpleFraction();
	SimpleFraction(long long divisible, long long divider);

	void Reduction();
	
	void SetDivisible(long long divisible) { m_divisible = divisible; }
	long long  GetDivisible() const { return m_divisible; }

	void SetDivider(long long divider) { m_divider = divider; }
	long long  GetDivider() const { return m_divider; }

	bool IsOne() const { return m_divisible == m_divider; }
	bool IsZero() const { return m_divisible == 0; }
	bool IsMax() const { return m_divisible == LLONG_MAX && m_divider == 1;}
	bool IsMin() const { return m_divisible == LLONG_MIN && m_divider == 1; }


	static void ToDivider(SimpleFraction& first, SimpleFraction& second);
	SimpleFraction ToAbs();
	std::string ToString() { return std::to_string(m_divisible) + "/" + std::to_string(m_divider); }

	SimpleFraction operator+(SimpleFraction second);
	SimpleFraction& operator+=(SimpleFraction second);
	SimpleFraction operator-(SimpleFraction second);
	SimpleFraction& operator-=(SimpleFraction second);
	SimpleFraction operator*(SimpleFraction second);
	SimpleFraction& operator*=(SimpleFraction second);
	SimpleFraction operator/(SimpleFraction second);
	SimpleFraction& operator/=(SimpleFraction second);

	friend inline std::ostream& operator<< (std::ostream& out, const Fraction::SimpleFraction& fraction);

	friend bool operator== (SimpleFraction first, SimpleFraction second);
	friend bool operator!= (SimpleFraction first, SimpleFraction second);
	friend bool operator> (SimpleFraction first, SimpleFraction second);
	friend bool operator< (SimpleFraction first, SimpleFraction second);

protected:
	long long m_divisible;	// Делимое
	long long m_divider;	// Делитель

	static long long NOD(long long a, long long b);
	static long long NOK(long long a, long long b);
};

inline std::ostream& operator<< (std::ostream& out, const SimpleFraction& fraction) {
	if (fraction.m_divider == 1) {
		out << fraction.m_divisible;
	}
	else {
		out << fraction.m_divisible << "/" << fraction.m_divider;
	}
	return out;
}

inline bool operator== (SimpleFraction first, SimpleFraction second)
{
	return first.m_divisible == second.m_divisible && first.m_divider == second.m_divider;
}

inline bool operator!= (SimpleFraction first, SimpleFraction second)
{
	return first.m_divisible != second.m_divisible || first.m_divider != second.m_divider;
}

inline bool operator> (SimpleFraction first, SimpleFraction second)
{
	SimpleFraction tempFirst = first;
	SimpleFraction tempSecond = second;

	SimpleFraction::ToDivider(tempFirst, tempSecond);
	
	return tempFirst.m_divisible > tempSecond.m_divisible;
}

inline bool operator< (SimpleFraction first, SimpleFraction second)
{
	SimpleFraction tempFirst = first;
	SimpleFraction tempSecond = second;

	SimpleFraction::ToDivider(tempFirst, tempSecond);

	return tempFirst.m_divisible < tempSecond.m_divisible;
}
}