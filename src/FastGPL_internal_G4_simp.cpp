#include <algorithm>
#include <complex>
#include <vector>

#include "FastGPL.h"

#include "FastGPL_internal.h"

using std::complex;
using std::vector;

using FastGPL::FastGPL_error;
using FastGPL::Log;
using FastGPL::PolyLog;

using FastGPL_internal::G;
using FastGPL_internal::G1;
using FastGPL_internal::G_Hoelder;
using FastGPL_internal::positive_int_power;
using FastGPL_internal::Zeta;

// G({a,0,0,a},x)
complex<double> G4_explicit_a00a(complex<double> a, int sa, double x)
{
    const complex<double> xa { x / a };
    const complex<double> xa1 { 1. - xa };

    const complex<double> Li2x { PolyLog(2, xa, -sa) };
    return -Li2x * Li2x / 2. - PolyLog(3, xa, -sa) * Log(xa1, sa);
}

// G({0,a,0,a},x)
complex<double> G4_explicit_0a0a(complex<double> a, int sa, double x)
{
    if (x == a) {
        return 3. * Zeta(4) / 4.;
    }

    const complex<double> Li2x { PolyLog(2, x / a, -sa) };
    return Li2x * Li2x / 2. - 2. * FastGPL::S22(x / a, -sa);
}

// G({0,0,a,a},x)
complex<double> G4_explicit_00aa(complex<double> a, int sa, double x)
{
    if (x == a) {
        return Zeta(4) / 4.;
    }

    return FastGPL::S22(x / a, -sa);
}

// G({0,a,a,a},x)
complex<double> G4_explicit_0aaa(complex<double> a, int sa, double x)
{
    if (x == a) {
        return -Zeta(4);
    }

    const complex<double> xa { x / a };
    const complex<double> xa1 { 1. - xa };

    const complex<double> l1 { Log(xa1, sa) };

    return PolyLog(4, xa1, sa) - Zeta(4)
        + l1
        * (Log(xa, -sa) * l1 * l1 / 6. + l1 * PolyLog(2, xa1, sa) / 2.
            - PolyLog(3, xa1, sa));
}

// G({a,0,a,a},x)
complex<double> G4_explicit_a0aa(complex<double> a, int sa, double x)
{
    const complex<double> xa { x / a };
    const complex<double> xa1 { 1. - xa };

    const complex<double> l1 { Log(xa1, sa) };

    return 3. * (Zeta(4) - PolyLog(4, xa1, sa))
        + l1 * (2. * PolyLog(3, xa1, sa) + Zeta(3) - l1 * PolyLog(2, xa1, sa) / 2.);
}

// G({a,a,0,a},x)
complex<double> G4_explicit_aa0a(complex<double> a, int sa, double x)
{
    const complex<double> xa { x / a };
    const complex<double> xa1 { 1. - xa };

    const complex<double> l1 { Log(xa1, sa) };

    return 3. * (PolyLog(4, xa1, sa) - Zeta(4))
        - l1 * (PolyLog(3, xa1, sa) + l1 * Zeta(2) / 2. + 2. * Zeta(3));
}

complex<double> G4_00ab(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G4_0a0b(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G4_a00b(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G4_0abc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_0abc_a(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_0abc_b(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_0abc_c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_a0bc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_a0bc_a(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_a0bc_b(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_a0bc_c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_ab0c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_ab0c_a(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_ab0c_b(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_ab0c_c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G4_abcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G4_abcd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G4_abcd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G4_abcd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G4_abcd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> FastGPL_internal::G4_dispatch(
    const vector<complex<double>>& a, const vector<int>& s, double x)
{
    // take care of G({0,0,a,a},x)
    if (a[0] == 0.0 && a[1] == 0.0 && a[2] == a[3]) {
        if (s[2] == s[3])
            return G4_explicit_00aa(a[3], s[3], x);
        else
            throw FastGPL_error { "GPL: a[2]==a[3] but s[2]!=s[3]" };
    }

    // take care of G({0,a,0,a},x)
    if (a[0] == 0.0 && a[2] == 0.0 && a[1] == a[3]) {
        if (s[1] == s[3])
            return G4_explicit_0a0a(a[3], s[3], x);
        else
            throw FastGPL_error { "GPL: a[1]==a[3] but s[1]!=s[3]" };
    }

    // take care of G({a,0,0,a},x)
    if (a[1] == 0.0 && a[2] == 0.0 && a[0] == a[3]) {
        if (s[0] == s[3])
            return G4_explicit_a00a(a[3], s[3], x);
        else
            throw FastGPL_error { "GPL: a[0]==a[3] but s[0]!=s[3]" };
    }

    // take care of G({0,a,a,a},x)
    if (a[0] == 0.0 && a[1] == a[2] && a[1] == a[3]) {
        if (s[1] == s[2] && s[1] == s[3])
            return G4_explicit_0aaa(a[3], s[3], x);
        else
            throw FastGPL_error { "GPL: a[1]==a[2]==a[3] but s[1]!=s[2]!=s[3]" };
    }

    // take care of G({a,0,a,a},x)
    if (a[1] == 0.0 && a[0] == a[2] && a[0] == a[3]) {
        if (s[0] == s[2] && s[0] == s[3])
            return G4_explicit_a0aa(a[3], s[3], x);
        else
            throw FastGPL_error { "GPL: a[0]==a[2]==a[3] but s[0]!=s[2]!=s[3]" };
    }

    // take care of G({a,a,0,a},x)
    if (a[2] == 0.0 && a[0] == a[1] && a[0] == a[3]) {
        if (s[0] == s[1] && s[0] == s[3])
            return G4_explicit_aa0a(a[3], s[3], x);
        else
            throw FastGPL_error { "GPL: a[0]==a[2]==a[3] but s[0]!=s[1]!=s[3]" };
    }

    if (is_convergent(a, x)) {
        return G_Hoelder(a, s, x);
    }

    if (a[0] == 0.0) {
        if (a[1] == 0.0) // take care of G({0,0,a,b},x)
            return G4_00ab(a[2], a[3], s[2], s[3], x);
        else if (a[2] == 0.0) // take care of G({0,a,0,b},x)
            return G4_0a0b(a[1], a[3], s[1], s[3], x);
        else // take care of G({0,a,b,c},x)
            return G4_0abc(a[1], a[2], a[3], s[1], s[2], s[3], x);
    }

    if (a[1] == 0.0) {
        if (a[2] == 0.0) // take care of G({a,0,0,b},x)
            return G4_a00b(a[0], a[3], s[0], s[3], x);
        else // take care of G({a,0,b,c},x)
            return G4_a0bc(a[0], a[2], a[3], s[0], s[2], s[3], x);
    }

    if (a[2] == 0.0) // take care of G({a,b,0,c},x)
        return G4_ab0c(a[0], a[1], a[3], s[0], s[1], s[3], x);
    else // take care of G({a,b,c,d},x)
        return G4_abcd(a[0], a[1], a[2], a[3], s[0], s[1], s[2], s[3], x);
}
