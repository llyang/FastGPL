#include <algorithm>
#include <complex>
#include <vector>

#include "FastGPL.h"

#include "FastGPL_internal.h"

#include <iostream>

using std::complex;
using std::vector;

using FastGPL_internal::factorial;
using FastGPL_internal::positive_int_power;

void clean_floating_point_vector(
    vector<complex<double>>::iterator begin, vector<complex<double>>::iterator end)
{
    constexpr double threshold { 1e-15 };

    auto last { std::prev(end) };

    if (begin == last)
        return;

    clean_floating_point_vector(begin, last);
    for (auto it = begin; it != last; ++it) {
        if (abs(*last - *it) < abs(*it) * threshold) {
            *last = *it;
            return;
        }
    }
}

complex<double> FastGPL::G(const vector<complex<double>>& a, double x, bool is_clean)
{
    vector<int> s(a.size(), 1);
    std::transform(
        a.cbegin(), a.cend(), s.begin(), [](auto z) { return z.imag() < 0.0 ? -1 : 1; });
    return G(a, s, x, is_clean);
}

complex<double> FastGPL::G(
    const vector<complex<double>>& b, vector<int> s, double x, bool is_clean)
{
    auto a = b;

    auto len { a.size() };

    if (len != s.size())
        throw FastGPL_error { "GPL: a.size() != s.size()" };

    if (len == 0)
        return 1;

    if (x < 0)
        throw FastGPL_error { "GPL: x < 0" };

    if ((a.back() == 0.0 && x == 0.0) || (a[0] != 0.0 && a[0] == x))
        throw FastGPL_error { "GPL: divergent" };

    if (!is_clean) {
        clean_floating_point_vector(a.begin(), a.end());
    }

    // s respect a
    std::transform(a.cbegin(), a.cend(), s.cbegin(), s.begin(),
        [](auto z, auto sz) { return z.imag() == 0.0 ? sz : (z.imag() < 0.0 ? -1 : 1); });

    if (a.back() != 0.0 && x == 0.0)
        return 0.0;

    if (len == 1)
        return FastGPL_internal::G1(a[0], s[0], x);

    // Check if it's a power of log
    bool is_log { std::all_of(
        ++a.cbegin(), a.cend(), [&a](auto i) { return i == a[0]; }) };
    if (is_log) {
        if (std::all_of(++s.cbegin(), s.cend(), [&s](auto i) { return i == s[0]; }))
            return positive_int_power(FastGPL_internal::G1(a[0], s[0], x), len)
                / factorial(len);
        else
            throw FastGPL_error { "GPL: is_log but s not the same" };
    }

    // Check if it's a PolyLog
    bool is_polylog { len > 1 && a.back() != 0.0
        && std::all_of(a.cbegin(), --a.cend(), [](auto i) { return i == 0.0; }) };
    if (is_polylog && abs(x / a.back()) > 0.2)
        return -PolyLog(len, x / a.back(), -s.back());

    // Check if it's well convergent
    if (FastGPL_internal::is_well_convergent(a, x))
        return FastGPL_internal::G_as_mLi(a, x);

    if (a[len - 1] == 0.0)
        return FastGPL_internal::G_remove_trailing_zeros(a, s, x);

    // Now there are no trailing zeros and no divergences

    // Explicit expression for weight 2 GPLs
    if (len == 2)
        return FastGPL_internal::G2_explicit(a, s, x);

    if (len == 3)
        return FastGPL_internal::G3_dispatch(a, s, x);

    if (len == 4)
        return FastGPL_internal::G4_dispatch(a, s, x);

    if (FastGPL_internal::is_convergent(a, x))
        return FastGPL_internal::G_Hoelder(a, s, x);

    if (len == 5)
        return FastGPL_internal::G5_dispatch(a, s, x);

    if (len == 6)
        return FastGPL_internal::G6_dispatch(a, s, x);

    throw FastGPL_error { "GPL: not implemented" };
}
