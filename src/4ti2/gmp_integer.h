/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2026 4ti2 team.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef _4ti2__gmp_integer_h_
#define _4ti2__gmp_integer_h_

#include <gmp.h>
#include <cstring>
#include <istream>
#include <ostream>
#include <string>

namespace _4ti2_gmp_
{

class Integer
{
public:
    Integer()
    {
        mpz_init(value_);
    }

    Integer(long long value)
    {
        mpz_init_set_si(value_, static_cast<long>(value));
    }

    Integer(const Integer& other)
    {
        mpz_init_set(value_, other.value_);
    }

    Integer& operator=(const Integer& other)
    {
        if (this != &other) {
            mpz_set(value_, other.value_);
        }
        return *this;
    }

    ~Integer()
    {
        mpz_clear(value_);
    }

    void set_mpz(mpz_srcptr value)
    {
        mpz_set(value_, value);
    }

    mpz_ptr get_mpz_t()
    {
        return value_;
    }

    mpz_srcptr get_mpz_t() const
    {
        return value_;
    }

    double get_d() const
    {
        return mpz_get_d(value_);
    }

    long get_si() const
    {
        return mpz_get_si(value_);
    }

    explicit operator double() const
    {
        return mpz_get_d(value_);
    }

    Integer operator-() const
    {
        Integer result;
        mpz_neg(result.value_, value_);
        return result;
    }

    Integer& operator+=(const Integer& rhs)
    {
        mpz_add(value_, value_, rhs.value_);
        return *this;
    }

    Integer& operator-=(const Integer& rhs)
    {
        mpz_sub(value_, value_, rhs.value_);
        return *this;
    }

    Integer& operator*=(const Integer& rhs)
    {
        mpz_mul(value_, value_, rhs.value_);
        return *this;
    }

    Integer& operator/=(const Integer& rhs)
    {
        mpz_tdiv_q(value_, value_, rhs.value_);
        return *this;
    }

    Integer& operator%=(const Integer& rhs)
    {
        mpz_tdiv_r(value_, value_, rhs.value_);
        return *this;
    }

    friend Integer operator+(Integer lhs, const Integer& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    friend Integer operator+(Integer lhs, long long rhs)
    {
        lhs += Integer(rhs);
        return lhs;
    }

    friend Integer operator+(long long lhs, Integer rhs)
    {
        rhs += Integer(lhs);
        return rhs;
    }

    friend Integer operator-(Integer lhs, const Integer& rhs)
    {
        lhs -= rhs;
        return lhs;
    }

    friend Integer operator-(Integer lhs, long long rhs)
    {
        lhs -= Integer(rhs);
        return lhs;
    }

    friend Integer operator-(long long lhs, const Integer& rhs)
    {
        Integer result(lhs);
        result -= rhs;
        return result;
    }

    friend Integer operator*(Integer lhs, const Integer& rhs)
    {
        lhs *= rhs;
        return lhs;
    }

    friend Integer operator*(Integer lhs, long long rhs)
    {
        lhs *= Integer(rhs);
        return lhs;
    }

    friend Integer operator*(long long lhs, Integer rhs)
    {
        rhs *= Integer(lhs);
        return rhs;
    }

    friend Integer operator/(Integer lhs, const Integer& rhs)
    {
        lhs /= rhs;
        return lhs;
    }

    friend double operator/(const Integer& lhs, double rhs)
    {
        return lhs.get_d() / rhs;
    }

    friend double operator/(double lhs, const Integer& rhs)
    {
        return lhs / rhs.get_d();
    }

    friend Integer operator%(Integer lhs, const Integer& rhs)
    {
        lhs %= rhs;
        return lhs;
    }

    friend bool operator==(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) == 0;
    }

    friend bool operator==(const Integer& lhs, long long rhs)
    {
        return mpz_cmp_si(lhs.value_, static_cast<long>(rhs)) == 0;
    }

    friend bool operator==(long long lhs, const Integer& rhs)
    {
        return mpz_cmp_si(rhs.value_, static_cast<long>(lhs)) == 0;
    }

    friend bool operator!=(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) != 0;
    }

    friend bool operator!=(const Integer& lhs, long long rhs)
    {
        return mpz_cmp_si(lhs.value_, static_cast<long>(rhs)) != 0;
    }

    friend bool operator!=(long long lhs, const Integer& rhs)
    {
        return mpz_cmp_si(rhs.value_, static_cast<long>(lhs)) != 0;
    }

    friend bool operator<(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) < 0;
    }

    friend bool operator<(const Integer& lhs, long long rhs)
    {
        return mpz_cmp_si(lhs.value_, static_cast<long>(rhs)) < 0;
    }

    friend bool operator<(long long lhs, const Integer& rhs)
    {
        return mpz_cmp_si(rhs.value_, static_cast<long>(lhs)) > 0;
    }

    friend bool operator<=(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) <= 0;
    }

    friend bool operator<=(const Integer& lhs, long long rhs)
    {
        return mpz_cmp_si(lhs.value_, static_cast<long>(rhs)) <= 0;
    }

    friend bool operator<=(long long lhs, const Integer& rhs)
    {
        return mpz_cmp_si(rhs.value_, static_cast<long>(lhs)) >= 0;
    }

    friend bool operator>(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) > 0;
    }

    friend bool operator>(const Integer& lhs, long long rhs)
    {
        return mpz_cmp_si(lhs.value_, static_cast<long>(rhs)) > 0;
    }

    friend bool operator>(long long lhs, const Integer& rhs)
    {
        return mpz_cmp_si(rhs.value_, static_cast<long>(lhs)) < 0;
    }

    friend bool operator>=(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) >= 0;
    }

    friend bool operator>=(const Integer& lhs, long long rhs)
    {
        return mpz_cmp_si(lhs.value_, static_cast<long>(rhs)) >= 0;
    }

    friend bool operator>=(long long lhs, const Integer& rhs)
    {
        return mpz_cmp_si(rhs.value_, static_cast<long>(lhs)) <= 0;
    }

    friend std::ostream& operator<<(std::ostream& out, const Integer& value)
    {
        char* text = mpz_get_str(0, 10, value.value_);
        out << text;
        void (*freefunc)(void*, size_t);
        mp_get_memory_functions(0, 0, &freefunc);
        freefunc(text, std::strlen(text) + 1);
        return out;
    }

    friend std::istream& operator>>(std::istream& in, Integer& value)
    {
        std::string text;
        in >> text;
        if (!in) {
            return in;
        }
        if (mpz_set_str(value.value_, text.c_str(), 10) != 0) {
            in.setstate(std::ios::failbit);
        }
        return in;
    }

private:
    mpz_t value_;
};

inline void gcd(Integer& result, const Integer& a, const Integer& b)
{
    mpz_gcd(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
}

inline int calc_precision(const Integer& value)
{
    return mpz_sgn(value.get_mpz_t()) == 0 ? 0 : mpz_sizeinbase(value.get_mpz_t(), 2);
}

} // namespace _4ti2_gmp_

#endif