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
#include <cstdint>
#include <cstddef>

// Portable helper: set mpz_t from int64_t using mpz_import (avoids sizeof(long) assumption).
static inline void mpz_set_int64(mpz_ptr z, int64_t v)
{
    // Two's-complement trick: avoids signed overflow when v == INT64_MIN.
    // -(v+1)+1 == -v for all v, but -(v+1) never overflows since v+1 >= INT64_MIN+1.
    uint64_t uv = v < 0 ? (uint64_t)(-(v + 1)) + 1 : (uint64_t)v;
    mpz_import(z, 1, 1, sizeof(uv), 0, 0, &uv);
    if (v < 0) mpz_neg(z, z);
}

// Portable helper: check if mpz_t value fits in int64_t (replaces mpz_fits_slong_p when sizeof(long) < 8).
static inline int mpz_fits_int64_p(mpz_srcptr v)
{
    if (mpz_sgn(v) == 0) return 1;
    // mpz_sizeinbase may overestimate by at most 1.  Threshold is 65 (not 64) so that a
    // value with exactly 64 real bits is not rejected before the exact mpz_export check below.
    // Values reported as > 65 bits have at least 65 real bits, which exceeds int64_t range.
    if (mpz_sizeinbase(v, 2) > 65) return 0;
    // At most 65 reported bits means at most 65 real bits, requiring at most 2 uint64_t words.
    // Use a 2-element buffer so mpz_export cannot write beyond it for these bounded values.
    uint64_t buf[2];
    size_t count;
    mpz_export(buf, &count, 1, sizeof(uint64_t), 0, 0, v);
    if (count == 0) return 1;
    if (count > 1) return 0;
    // Positive: must be <= INT64_MAX (< 2^63); negative: magnitude must be <= 2^63.
    return mpz_sgn(v) > 0 ? buf[0] <= (uint64_t)INT64_MAX
                           : buf[0] <= (uint64_t)INT64_MAX + 1;
}

// Portable helper: extract int64_t from mpz_t (caller must ensure mpz_fits_int64_p is true).
static inline int64_t mpz_get_int64(mpz_srcptr v)
{
    uint64_t limb;
    size_t count;
    mpz_export(&limb, &count, 1, sizeof(limb), 0, 0, v);
    if (count == 0) return 0;
    if (mpz_sgn(v) > 0) return (int64_t)limb;
    // Two's-complement negation without overflow: -(limb-1)-1 handles limb==2^63 (INT64_MIN).
    return -(int64_t)(limb - 1) - 1;
}

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
        mpz_init(value_);
        mpz_set_int64(value_, static_cast<int64_t>(value));
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
        Integer tmp(rhs);
        return mpz_cmp(lhs.value_, tmp.value_) == 0;
    }

    friend bool operator==(long long lhs, const Integer& rhs)
    {
        Integer tmp(lhs);
        return mpz_cmp(tmp.value_, rhs.value_) == 0;
    }

    friend bool operator!=(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) != 0;
    }

    friend bool operator!=(const Integer& lhs, long long rhs)
    {
        Integer tmp(rhs);
        return mpz_cmp(lhs.value_, tmp.value_) != 0;
    }

    friend bool operator!=(long long lhs, const Integer& rhs)
    {
        Integer tmp(lhs);
        return mpz_cmp(tmp.value_, rhs.value_) != 0;
    }

    friend bool operator<(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) < 0;
    }

    friend bool operator<(const Integer& lhs, long long rhs)
    {
        Integer tmp(rhs);
        return mpz_cmp(lhs.value_, tmp.value_) < 0;
    }

    friend bool operator<(long long lhs, const Integer& rhs)
    {
        Integer tmp(lhs);
        return mpz_cmp(tmp.value_, rhs.value_) < 0;
    }

    friend bool operator<=(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) <= 0;
    }

    friend bool operator<=(const Integer& lhs, long long rhs)
    {
        Integer tmp(rhs);
        return mpz_cmp(lhs.value_, tmp.value_) <= 0;
    }

    friend bool operator<=(long long lhs, const Integer& rhs)
    {
        Integer tmp(lhs);
        return mpz_cmp(tmp.value_, rhs.value_) <= 0;
    }

    friend bool operator>(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) > 0;
    }

    friend bool operator>(const Integer& lhs, long long rhs)
    {
        Integer tmp(rhs);
        return mpz_cmp(lhs.value_, tmp.value_) > 0;
    }

    friend bool operator>(long long lhs, const Integer& rhs)
    {
        Integer tmp(lhs);
        return mpz_cmp(tmp.value_, rhs.value_) > 0;
    }

    friend bool operator>=(const Integer& lhs, const Integer& rhs)
    {
        return mpz_cmp(lhs.value_, rhs.value_) >= 0;
    }

    friend bool operator>=(const Integer& lhs, long long rhs)
    {
        Integer tmp(rhs);
        return mpz_cmp(lhs.value_, tmp.value_) >= 0;
    }

    friend bool operator>=(long long lhs, const Integer& rhs)
    {
        Integer tmp(lhs);
        return mpz_cmp(tmp.value_, rhs.value_) >= 0;
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