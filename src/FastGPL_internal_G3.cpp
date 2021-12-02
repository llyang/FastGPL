#include <algorithm>
#include <complex>
#include <vector>

#include <iostream>

#include "FastGPL.h"

#include "FastGPL_internal.h"

using std::complex;
using std::vector;

using FastGPL::FastGPL_error;
using FastGPL::Log;
using FastGPL::PolyLog;

using FastGPL_internal::G;
using FastGPL_internal::G1;
using FastGPL_internal::Zeta;

// G({0,a,a},x)
complex<double> G3_explicit_0aa(complex<double> a, int s, double x)
{
    if (x == a)
        return Zeta(3);

    const complex<double> xa { x / a };
    const complex<double> xa1 { 1. - xa };

    const complex<double> l1 { Log(xa1, s) };

    if (xa1.real() > 1.) {
        return PolyLog(3, xa, -s) + PolyLog(3, x / (x - a), s)
            - l1 * (PolyLog(2, xa, -s) + l1 * l1 / 6.);
    } else {
        return l1 * (PolyLog(2, xa1, s) + Log(xa, -s) * l1 / 2.) - PolyLog(3, xa1, s)
            + Zeta(3);
    }
}

// G({a,0,a},x)
complex<double> G3_explicit_a0a(complex<double> a, int s, double x)
{
    const complex<double> xa { x / a };
    const complex<double> xa1 { 1. - xa };

    if (xa1.real() > 1.) {
        const complex<double> l1 { Log(xa1, s) };
        return -2. * (PolyLog(3, xa, -s) + PolyLog(3, x / (x - a), s))
            + l1 * (PolyLog(2, xa, -s) + l1 * l1 / 3.);
    } else {
        return 2. * (PolyLog(3, xa1, s) - Zeta(3))
            - Log(xa1, s) * (PolyLog(2, xa1, s) + Zeta(2));
    }
}

// G({a,0,b},x)
// This should only be called if abs(a) < x and/or abs(b) < x
complex<double> G3_a0b(complex<double> a, complex<double> b, int sa, int sb, double x)
{
    if (abs(b) < abs(a)) {
        const complex<double> xb { x / b };
        const complex<double> ab { a / b };
        const complex<double> Lxb { Log(-xb, sb) };
        const complex<double> G0ab { G({ 0, ab }, 1) };
        complex<double> res { G({ 0, 0, a }, { 1, 1, sa }, x) + 2. * G({ 0, 0, ab }, 1)
            + G({ 0, ab, xb }, 1)
            + G1(a, sa, x) * (G({ 0, b }, { 1, sb }, x) + G({ 0, xb }, 1) - G0ab)
            + Lxb * (G0ab - G({ 0, a }, { 1, sa }, x)) };
        return res;
    } else { // abs(a) < abs(b)
        if (b == x) {
            const complex<double> xa { x / a };
            const complex<double> Lxa { Log(-xa, sa) };
            complex<double> res { 2. * G({ 0, 0, xa }, 1) + G({ 0, xa, xa }, 1) - Zeta(3)
                - G({ xa }, 1) * Zeta(2) + G({ 0, xa }, 1) * Lxa };
            return res;
        } else {
            const complex<double> xa { x / a };
            const complex<double> ba { b / a };
            const complex<double> Lxa { Log(-xa, sa) };
            const complex<double> G0ba { G({ 0, ba }, 1) };
            complex<double> res { G({ 0, 0, b }, { 1, 1, sb }, x)
                + 2. * G({ 0, 0, ba }, 1) + G({ 0, ba, xa }, 1)
                + G({ xa }, 1) * G({ 0, b }, { 1, sb }, x) + G0ba * Lxa
                + G1(b, sb, x) * (G({ 0, xa }, 1) - G0ba) };
            return res;
        }
    }
}

// G({0,a,b},x)
// This should only be called if abs(a) < x and/or abs(b) < x
complex<double> G3_0ab(complex<double> a, complex<double> b, int sa, int sb, double x)
{
    if (abs(b) < abs(a)) {
        if (a == x) {
            const complex<double> xb { x / b };
            const complex<double> Lxb { Log(-xb, sb) };
            const complex<double> res { 2. * Zeta(3) - G({ 0, 0, xb }, 1)
                + G({ xb, 0, xb }, 1) + G({ xb }, 1) * (-Lxb * Lxb / 2.0 - 2. * Zeta(2))
                - (Zeta(2) + G({ 0, xb }, 1)) * Lxb };
            return res;
        } else {
            const complex<double> xb { x / b };
            const complex<double> ab { a / b };
            const complex<double> Lxb { Log(-xb, sb) };
            const complex<double> G0ax { G({ 0, a }, { 1, sa }, x) };
            complex<double> res { G({ ab, 0, xb }, 1) - G({ 0, 0, ab }, 1)
                - 2. * G({ 0, 0, a }, { 1, 1, sa }, x)
                + G({ ab }, 1) * (G0ax - Lxb * Lxb / 2.0 - Zeta(2))
                + (G0ax - G({ 0, ab }, 1)) * Lxb };
            return res;
        }
    } else { // abs(a) < abs(b)
        const complex<double> xa { x / a };
        const complex<double> ba { b / a };
        const complex<double> Lxa { Log(-xa, sa) };
        const complex<double> G0xa { G({ 0, xa }, 1) };
        complex<double> res { G({ 0, 0, b }, { 1, 1, sb }, x) - G({ 0, 0, ba }, 1)
            - G({ ba, 0, xa }, 1) - G({ 0, ba, xa }, 1)
            + G({ ba }, 1) * (Lxa * Lxa / 2.0 + Zeta(2) - G({ 0, b }, { 1, sb }, x)) };
        if (b != x) {
            res += G1(b, sb, x) * (G({ 0, ba }, 1) - G0xa);
        }
        return res;
    }
}

// G({a,b,c},x) with a smallest
complex<double> G3_abc_a(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    if (a == b) {
        // take care of G({a,a,c},x)
        if (sa != sb)
            throw FastGPL::FastGPL_error { "G3_abc: a==b but sa!=sb" };
        else {
            const complex<double> xa { x / a };
            const complex<double> ca { c / a };
            const complex<double> Gcaa { G({ c / a, 1 }, 1) };
            complex<double> res { G({ 0, a, c }, { 1, sa, sc }, x) - G({ 0, ca, 1 }, 1)
                - G({ ca, 1, xa }, 1) - G({ ca, xa, 1 }, 1)
                + G({ xa }, 1) * G({ a, c }, { sa, sc }, x) - Gcaa * Log(-xa, sa) };
            if (c != x) {
                res += G1(c, sc, x) * (Gcaa - G({ xa, 1 }, 1));
            }
            return res;
        }
    } else if (a == c) {
        // take care of G({a,b,a},x)
        if (sa != sc)
            throw FastGPL::FastGPL_error { "G3_abc: a==c but sa!=sc" };
        else {
            const complex<double> xa { x / a };
            const complex<double> ba { b / a };
            const complex<double> Gba { G({ ba }, 1) };
            const complex<double> Gbca { G({ ba, 1 }, 1) };
            complex<double> res { G({ 0, b, a }, { 1, sb, sa }, x) + G({ ba, 1, xa }, 1)
                - G({ 0, 1, ba }, 1) + Gbca * Log(-xa, sa)
                + Gba * (G({ 0, a }, { 1, sa }, x) - Zeta(2))
                + G1(a, sa, x) * (G({ ba, xa }, 1) - Gbca) };
            if (b != x) {
                res += G({ b, c }, { sb, sc }, x) * (G({ xa }, 1) - Gba);
            }
            return res;
        }
    } else {
        const complex<double> xa { x / a };
        const complex<double> ba { b / a };
        const complex<double> ca { c / a };
        const complex<double> Gba { G({ ba }, 1) };
        const complex<double> Gbca { G({ ba, ca }, 1) };
        complex<double> res { G({ 0, b, c }, { 1, sb, sc }, x) + G({ ba, ca, xa }, 1)
            - G({ 0, ca, ba }, 1) + Gbca * Log(-xa, sa)
            + Gba * (G({ 0, c }, { 1, sc }, x) + G({ 0, ca }, 1)) };
        if (b != x) {
            res += G({ b, c }, { sb, sc }, x) * (G({ xa }, 1) - Gba);
        }
        if (c != x) {
            res += G1(c, sc, x) * (G({ ba, xa }, 1) - Gbca);
        }
        return res;
    }
}

// G({a,b,c},x) with b smallest
complex<double> G3_abc_b(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{

    if (b == c) {
        // take care of G({a,b,b},x)
        if (sb != sc)
            throw FastGPL::FastGPL_error { "G3_abc: b==c but sb!=sc" };
        else {
            const complex<double> xb { x / b };
            const complex<double> ab { a / b };
            const complex<double> Gax { G1(a, sa, x) };
            const complex<double> Gbx { G1(b, sb, x) };
            return G({ a, 0, b }, { sa, 1, sb }, x) + G({ ab, xb, 1 }, 1)
                - G({ ab, 0, 1 }, 1) + G({ ab, 1 }, 1) * (Gbx - Gax)
                + G({ ab }, 1) * (G({ a, b }, { sa, sb }, x) - G({ 0, b }, { 1, sb }, x))
                - Gbx * G({ ab, xb }, 1) - Gax * Zeta(2);
        }
    } else {
        const complex<double> xb { x / b };
        const complex<double> ab { a / b };
        const complex<double> cb { c / b };
        const complex<double> Gax { G1(a, sa, x) };
        const complex<double> Gacx { G({ a, c }, { sa, sc }, x) };
        const complex<double> Gaxb { G({ ab, xb }, 1) };
        const complex<double> Gab { G({ ab }, 1) };
        const complex<double> Gcb { G({ cb }, 1) };

        complex<double> res { G({ a, 0, c }, { sa, 1, sc }, x) + G({ ab, xb, cb }, 1)
            - G({ ab, 0, cb }, 1) + (Gax - Gab) * Gcb * Log(-xb, sb)
            + Gax * (G({ 0, cb }, 1) + G({ cb, ab }, 1))
            + Gab * (Gacx - G({ 0, c }, { 1, sc }, x))
            - Gcb * (Gacx + G({ 0, a }, { 1, sa }, x) + G({ 0, ab }, 1) + Gaxb) };
        if (c != x) {
            res += G1(c, sc, x) * (G({ ab, cb }, 1) - Gaxb);
        }
        return res;
    }
}

// G({a,b,c},x) with c smallest, c!=a and c!=b
complex<double> G3_abc_c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    const complex<double> xc { x / c };
    const complex<double> ac { a / c };
    const complex<double> bc { b / c };
    const complex<double> Lxc { Log(-xc, sc) };
    const complex<double> Gax { G1(a, sa, x) };
    const complex<double> Gbc { G({ bc }, 1) };
    const complex<double> Gbac { G({ bc, ac }, 1) };
    const complex<double> Gabx { G({ a, b }, { sa, sb }, x) };

    return G({ bc, ac, xc }, 1) + G({ bc, 0, ac }, 1) + G({ 0, bc, ac }, 1)
        - G({ 0, a, b }, { 1, sa, sb }, x) - G({ a, 0, b }, { sa, 0, sb }, x)
        - Gax * (G({ 0, bc }, 1) + Gbac) + Gbc * (Gabx + G({ 0, a }, { 1, sa }, x))
        + (Gabx + Gbac - Gbc * Gax) * Lxc;
}

// G(a,b,c,x)
// This should only be called if abs(a) < x and/or abs(b) < x and/or abs(c) < x
complex<double> G3_abc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G3_abc_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G3_abc_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G3_abc_b(a, b, c, sa, sb, sc, x);
}

complex<double> FastGPL_internal::G3_dispatch(
    const vector<complex<double>>& a, const vector<int>& s, double x)
{
    // take care of G({0,a,a},x)
    if (a[0] == 0.0 && a[1] == a[2]) {
        if (s[1] == s[2])
            return G3_explicit_0aa(a[2], s[2], x);
        else
            throw FastGPL_error { "GPL: a[1]==a[2] but s[1]!=s[2]" };
    }

    // take care of G({a,0,a},x)
    if (a[1] == 0.0 && a[0] == a[2]) {
        if (s[0] == s[2])
            return G3_explicit_a0a(a[2], s[2], x);
        else
            throw FastGPL_error { "GPL: a[0]==a[2] but s[0]!=s[2]" };
    }

    if (is_convergent(a, x))
        return G_Hoelder(a, s, x);

    if (a[0] == 0.0) {
        return G3_0ab(a[1], a[2], s[1], s[2], x);
    } else if (a[1] == 0.0) {
        return G3_a0b(a[0], a[2], s[0], s[2], x);
    } else {
        return G3_abc(a[0], a[1], a[2], s[0], s[1], s[2], x);
    }
}
