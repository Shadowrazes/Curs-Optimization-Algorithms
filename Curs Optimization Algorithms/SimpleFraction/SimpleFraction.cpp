#include "SimpleFraction.h"
#include <algorithm>

namespace Fraction {
//-----------------------------------------------------------------------------------
SimpleFraction::SimpleFraction() {
	m_divisible = 1;
	m_divider = 1;
}
//-----------------------------------------------------------------------------------
SimpleFraction::SimpleFraction(long long divisible, long long divider) {
	m_divisible = divisible;
	m_divider = divider;
}
//-----------------------------------------------------------------------------------
long long SimpleFraction::NOD(long long a, long long b) {
    if (a < b) {
        std::swap(a, b);
    }
    while (a % b != 0) {
        a = a % b;
        std::swap(a, b);
    }
    return b;
}
//-----------------------------------------------------------------------------------
long long SimpleFraction::NOK(long long a, long long b) {
    return a * b / NOD(a, b);
}
//-----------------------------------------------------------------------------------
void SimpleFraction::Reduction() {
    if (m_divider < 0) {
        m_divisible *= -1;
        m_divider *= -1;
    }

    if (m_divider == 0) {
        return;
    }

    if (m_divisible == 0) {
        m_divider = 1;
        return;
    }

    long long nod = NOD(m_divisible, m_divider);

    if (m_divisible % nod != 0 || m_divider % nod != 0) {
        return;
    }
    else {
        m_divisible /= nod;
        m_divider /= nod;

        if (m_divider < 0) {
            m_divisible *= -1;
            m_divider *= -1;
        }
    }
}
//-----------------------------------------------------------------------------------
void SimpleFraction::ToDivider(SimpleFraction& first, SimpleFraction& second) {
    long long destination = NOK(first.m_divider, second.m_divider);

    long long firstCoef = destination / first.m_divider;
    long long secondCoef = destination / second.m_divider;

    first.m_divisible *= firstCoef;
    first.m_divider *= firstCoef;

    second.m_divisible *= secondCoef;
    second.m_divider *= secondCoef;
}
//-----------------------------------------------------------------------------------
SimpleFraction SimpleFraction::ToAbs() {
    if (m_divisible < 0)
        return SimpleFraction(m_divisible * -1, m_divider);
    return *this;
}
//-----------------------------------------------------------------------------------
SimpleFraction SimpleFraction::operator+(SimpleFraction second) {
    SimpleFraction tempThis = *this;
    ToDivider(tempThis, second);
    tempThis.m_divisible += second.m_divisible;
    tempThis.Reduction();
    return tempThis;
}
//-----------------------------------------------------------------------------------
SimpleFraction& SimpleFraction::operator+=(SimpleFraction second) {
    ToDivider(*this, second);
    m_divisible += second.m_divisible;
    Reduction();
    return *this;
}
//-----------------------------------------------------------------------------------
SimpleFraction SimpleFraction::operator-(SimpleFraction second) {
    SimpleFraction tempThis = *this;
    ToDivider(tempThis, second);
    tempThis.m_divisible -= second.m_divisible;
    tempThis.Reduction();
    return tempThis;
}
//-----------------------------------------------------------------------------------
SimpleFraction& SimpleFraction::operator-=(SimpleFraction second) {
    ToDivider(*this, second);
    m_divisible -= second.m_divisible;
    Reduction();
    return *this;
}
//-----------------------------------------------------------------------------------
SimpleFraction SimpleFraction::operator*(SimpleFraction second) {
    SimpleFraction answer(m_divisible * second.m_divisible, m_divider * second.m_divider);
    answer.Reduction();
    return answer;
}
//-----------------------------------------------------------------------------------
SimpleFraction& SimpleFraction::operator*=(SimpleFraction second) {
    m_divisible *= second.m_divisible;
    m_divider *= second.m_divider;
    Reduction();
    return *this;
}
//-----------------------------------------------------------------------------------
SimpleFraction SimpleFraction::operator/(SimpleFraction second) {
    SimpleFraction answer(m_divisible * second.m_divider, m_divider * second.m_divisible);
    answer.Reduction();
    return answer;
}
//-----------------------------------------------------------------------------------
SimpleFraction& SimpleFraction::operator/=(SimpleFraction second) {
    m_divisible *= second.m_divider;
    m_divider *= second.m_divisible;
    Reduction();
    return *this;
}
//-----------------------------------------------------------------------------------
}