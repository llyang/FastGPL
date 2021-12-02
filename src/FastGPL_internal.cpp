#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include <iostream>

#include "FastGPL.h"

#include "FastGPL_internal.h"

using std::complex;
using std::vector;

using FastGPL::FastGPL_error;
using FastGPL::Log;
using FastGPL::PolyLog;

// Check whether the series representation is convergent
bool FastGPL_internal::is_convergent(const vector<complex<double>>& a, double x)
{
    return a[a.size() - 1] != 0.0 && std::all_of(a.cbegin(), a.cend(), [x](auto c) {
        return c == 0.0 || abs(c) >= std::abs(x);
    });
}

// Check whether the series representation is well convergent
bool FastGPL_internal::is_well_convergent(const vector<complex<double>>& a, double x)
{
    return a[a.size() - 1] != 0.0 && std::all_of(a.cbegin(), a.cend(), [x](auto c) {
        return c == 0.0 || abs(c) >= std::abs(x) * FastGPL::Hoelder_threshold;
    });
}

// Convert convergent G (with x>0 and without trailing zeros) to mLi
complex<double> FastGPL_internal::G_as_mLi(const vector<complex<double>>& a, double x)
{
    vector<int> m;
    vector<complex<double>> z;

    for (auto it = a.cbegin(); it != a.cend();) {
        auto next { std::find_if(it, a.cend(), [](auto c) { return c != 0.0; }) };
        if (next != a.cend()) {
            m.push_back(next - it + 1);
            z.push_back(*next / x);
            ++next;
        }
        it = next;
    }

    double sign { m.size() % 2 == 0 ? 1.0 : -1.0 };
    vector<complex<double>> y;
    bool good { true };
    for (auto [p, it] = std::tuple { complex<double> { 1.0 }, z.cbegin() };
         it != z.cend(); p = *it++) {
        y.push_back(p / *it);
        if (abs(y.back()) > 2.) // magic number
            good = false;
    }

    if (good) {
        return sign * mLi_series<double>(m, y);
    } else {
        return sign * static_cast<complex<double>>(mLi_series<long double>(m, y));
    }
}

// Series representation of multiple polylogarithms
// It should be guaranteed by the caller that it is convergent
template <class T>
complex<T> FastGPL_internal::mLi_series(
    const vector<int>& m, const vector<complex<double>>& x)
{
    const int n { static_cast<int>(m.size()) };

    vector<complex<T>> t(n + 1, 0.0);
    t[n] = 1.0;

    complex<T> res;

    int q { 0 };

    vector<complex<T>> x_power(n);
    for (int k = n - 1; k >= 0; --k)
        x_power[k] = positive_int_power(x[k], n - 1 - k);

    constexpr int q_slow { 1000 }; // magic number

    do {
        if (q > q_slow) // magic number
            throw FastGPL_error { "GPL: mLi_series converges too slowly" };

        res = t[0];

        for (int i = 0; i < 2; ++i) {
            ++q;
            for (int k = n - 1; k >= 0; --k) {
                x_power[k] *= x[k]; // positive_int_power(x[k], q + n - 1 - k)
                t[k] += t[k + 1] * x_power[k]
                    / static_cast<T>(
                        positive_int_power(static_cast<double>(q + n - 1 - k), m[k]));
                if (!(std::isfinite(t[k].real()) && std::isfinite(t[k].imag())))
                    throw FastGPL_error { "GPL: Inf or NaN encountered in mLi_series" };
            }
        }
    } while (abs(t[0] - res) > abs(res) * 1e-16);

    return t[0];
}
template complex<double> FastGPL_internal::mLi_series<double>(
    const vector<int>& m, const vector<complex<double>>& x);
template complex<long double> FastGPL_internal::mLi_series<long double>(
    const vector<int>& m, const vector<complex<double>>& x);

// Compute G({a}, {s}, x) without divergences
complex<double> FastGPL_internal::G1(complex<double> a, int s, double x)
{
    if (a == 0.0)
        return log(x);

    if (x == 0.0)
        return 0.0;

    complex<double> z { x / a };

    constexpr double series_threshold { 0.2 };
    if (abs(z) < series_threshold) {
        complex<double> res { -z };
        complex<double> power { z * z };
        complex<double> term { 0 };
        int n { 2 };
        do {
            if (n > 100)
                throw FastGPL_error { "GPL: G1 converges too slowly" };

            term = power / static_cast<double>(n);
            res -= term;

            power *= z;
            ++n;
        } while (abs(term) > abs(res) * 1e-16);

        return res;
    }

    constexpr double unity_threshold { 0.1 };
    if (abs(1.0 - z) < unity_threshold)
        return Log((a - x) / a, s);

    return Log(1.0 - z, s);
}

// Hoelder convolution
complex<double> FastGPL_internal::G_Hoelder(
    const vector<complex<double>>& a, const vector<int>& s, double x)
{
    vector<int> s_left;
    vector<complex<double>> a_left;
    vector<int> s_right { s };
    vector<complex<double>> a_right(a.size(), 0.0);
    std::transform(
        a.cbegin(), a.cend(), a_right.begin(), [x](auto i) { return i * 2.0 / x; });

    complex<double> res { G(a_right, s_right, 1) };
    double sign = 1;
    while (!a_right.empty()) {
        sign = -sign;
        complex<double> ta { a_right[0] };
        int ts { s_right[0] };
        a_right.erase(a_right.begin());
        a_left.insert(a_left.begin(), 2.0 - ta);
        s_right.erase(s_right.begin());
        s_left.insert(s_left.begin(), -ts);
        res += sign * G(a_left, s_left, 1) * G(a_right, s_right, 1);
    }

    return res;
}

// Remove trailing zeros
// The caller must ensure that a[a.size()-1]==0 && x>0
complex<double> FastGPL_internal::G_remove_trailing_zeros(
    vector<complex<double>> a, vector<int> s, double x)
{
    int n { static_cast<int>(a.size()) };

    int n_zeros { static_cast<int>(std::find_if(a.crbegin(), a.crend(), [](auto i) {
        return i != 0.0;
    }) - a.crbegin()) };

    // Normally this shouldn't be true
    if (n_zeros == n)
        return positive_int_power(log(x), n_zeros) / factorial(n_zeros);

    a.pop_back();
    s.pop_back();
    complex<double> res { log(x) * G(a, s, x) };

    a.insert(a.begin(), 0.0);
    s.insert(s.begin(), 1);
    for (int i = 0; i < n - n_zeros;) {
        complex<double> temp { G(a, s, x) };
        res -= temp;

        // If a[++i]==0, then the next term is the same
        for (++i; i < n - n_zeros && a[i] == 0.0; ++i)
            res -= temp;

        if (i < n - n_zeros) { // Now we hit a non-zero index
            a[i - 1] = a[i];
            a[i] = 0.0;
            s[i - 1] = s[i];
            s[i] = 1;
        }
    }

    res /= n_zeros;

    return res;
}

// Explicit expressions for weight 2 GPLs
// (note that simpler cases have already be dealt with earlier)
complex<double> FastGPL_internal::G2_explicit(
    const vector<complex<double>>& a, const vector<int>& s, double x)
{
    double a0r { a[0].real() };
    double a0i { a[0].imag() };
    double a1r { a[1].real() };
    double a1i { a[1].imag() };

    constexpr double unity_threshold { 0.1 };
    complex<double> a01 { 1.0 - a[0] / a[1] };
    if (abs(a01) < unity_threshold)
        a01 = (a[1] - a[0]) / a[1];

    constexpr double Hoelder_low_threshold { 1.01 };

    if (a[1] == x) {
        return G({ 0, a01 }, { 1, -s[0] }, 1);
    } else if (is_convergent(a, x) && abs(a[0] / x) > Hoelder_low_threshold
        && abs(a[1] / x) > Hoelder_low_threshold) {
        return G_Hoelder(a, s, x);
    } else if (abs(a[0]) < abs(a[1]) && abs(x / a[0]) > Hoelder_low_threshold
        && abs(a[1] / a[0]) > Hoelder_low_threshold) {
        complex<double> Gba { G({ a[1] / a[0] }, 1) };
        return G({ 0, a[1] }, { 1, s[1] }, x) + G({ a[1] / a[0], x / a[0] }, 1)
            + G({ 0, a[1] / a[0] }, 1) + Gba * (Log(-x, s[0]) - Log(a[0], s[0]))
            + G1(a[1], s[1], x) * (G({ x / a[0] }, 1) - Gba);
    } else if (abs(a[1]) < abs(a[0]) && abs(x / a[1]) > Hoelder_low_threshold
        && abs(a[0] / a[1]) > Hoelder_low_threshold) {
        return G1(a[0], s[0], x) * G1(a[1], s[1], x)
            - G({ a[1], a[0] }, { s[1], s[0] }, x);
    } else if (fabs(a0i) > fabs(a1i) || (a0i == 0.0 && a1i == 0.0 && a0r < a1r)) {
        return PolyLog(2, (a[1] - x) / (a[1] - a[0]), s[0])
            + G({ 0, a01 }, { 1, -s[0] }, 1)
            + G1(a[1], s[1], x) * Log((x - a[0]) / (a[1] - a[0]), -s[0]);
    } else {
        return G1(a[0], s[0], x) * Log(a01, s[1])
            - G({ 0, (a[0] - a[1]) / a[0] }, { 1, -s[1] }, 1)
            - PolyLog(2, (a[0] - x) / (a[0] - a[1]), s[1]);
    }
}
