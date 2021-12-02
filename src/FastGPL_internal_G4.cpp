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

// G({0,0,a,b},x)
// This should only be called if abs(a) < x and/or abs(b) < x
complex<double> G4_00ab(complex<double> a, complex<double> b, int sa, int sb, double x)
{
    if (abs(b) < abs(a)) {
        const complex<double> ab { a / b };
        const complex<double> xb { x / b };
        const complex<double> Lxb { Log(-xb, sb) };

        complex<double> res = Lxb * (G({ 0, 0, a }, { 1, 1, sa }, x) - G({ 0, 0, ab }, 1))
            - G({ ab, 0, 0, xb }, 1) - G({ 0, 0, 0, ab }, 1)
            - 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
            - G({ 0, ab }, 1) * (Lxb * Lxb / 2. + Zeta(2))
            + G({ ab }, 1)
                * (G({ 0, 0, a }, { 1, 1, sa }, x) - Lxb * (Lxb * Lxb / 6. + Zeta(2)));
        return res;
    }

    // abs(a) < abs(b)
    const complex<double> ba { b / a };
    const complex<double> xa { x / a };
    const complex<double> Lxa { Log(-xa, sa) };
    complex<double> res { G({ 0, 0, ba, xa }, 1) + G({ 0, ba, 0, xa }, 1)
        + G({ ba, 0, 0, xa }, 1) + G({ 0, ba }, 1) * G({ 0, b }, { 1, sb }, x)
        + G({ 0, 0, 0, ba }, 1) + G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
        + G({ ba }, 1)
            * (Lxa * (Lxa * Lxa / 6. + Zeta(2)) - G({ 0, 0, b }, { 1, 1, sb }, x)) };
    if (b != x) {
        res += G1(b, sb, x) * (G({ 0, 0, xa }, 1) - G({ 0, 0, ba }, 1));
    }
    return res;
}

// G({0,a,0,b},x)
// This should only be called if abs(a) < x and/or abs(b) < x
complex<double> G4_0a0b(complex<double> a, complex<double> b, int sa, int sb, double x)
{
    if (abs(b) < abs(a)) {
        const complex<double> ab { a / b };
        const complex<double> xb { x / b };
        const complex<double> Lxb { Log(-xb, sb) };
        const complex<double> G0ax { G({ 0, a }, { 1, sa }, x) };
        complex<double> res { 3. * G({ 0, 0, 0, ab }, 1)
            + 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) - G({ 0, ab, 0, xb }, 1)
            + G0ax * (G({ 0, xb }, 1) + G({ 0, b }, { 1, sb }, x))
            + 2. * Lxb * (G({ 0, 0, ab }, 1) - G({ 0, 0, a }, { 1, 1, sa }, x))
            + G({ 0, ab }, 1) * (Lxb * Lxb / 2. + Zeta(2) - G0ax) };
        return res;
    }

    // abs(a) < abs(b)
    const complex<double> ba { b / a };
    const complex<double> xa { x / a };
    const complex<double> Lxa { Log(-xa, sa) };
    const complex<double> G0bx { G({ 0, b }, { 1, sb }, x) };
    complex<double> res { G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
        - 2. * G({ 0, 0, ba, xa }, 1) - G({ 0, ba, 0, xa }, 1) - G({ 0, xa }, 1) * G0bx
        - 3. * G({ 0, 0, 0, ba }, 1)
        + G({ 0, ba }, 1) * (Lxa * Lxa / 2. + Zeta(2) - G0bx) };

    if (b != x) {
        res += 2. * G1(b, sb, x) * (G({ 0, 0, ba }, 1) - G({ 0, 0, xa }, 1));
    }
    return res;
}

// G({a,0,0,b},x)
// This should only be called if abs(a) < x and/or abs(b) < x
complex<double> G4_a00b(complex<double> a, complex<double> b, int sa, int sb, double x)
{
    if (abs(b) < abs(a)) {
        const complex<double> ab { a / b };
        const complex<double> xb { x / b };
        const complex<double> Lxb { Log(-xb, sb) };
        const complex<double> Gax { G1(a, sa, x) };
        complex<double> res { (Gax - Lxb) * G({ 0, 0, ab }, 1)
            + Lxb * G({ 0, 0, a }, { 1, 1, sa }, x)
            + Gax * (G({ 0, 0, b }, { 1, 1, sb }, x) - G({ 0, 0, xb }, 1))
            - G({ 0, 0, ab, xb }, 1)
            - (G({ 0, b }, { 1, sb }, x) + G({ 0, xb }, 1)) * G({ 0, a }, { 1, sa }, x)
            - 3. * G({ 0, 0, 0, ab }, 1) - G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) };
        return res;
    }

    // abs(a) < abs(b)
    const complex<double> ba { b / a };
    const complex<double> xa { x / a };
    const complex<double> Lxa { Log(-xa, sa) };
    const complex<double> G00ba { G({ 0, 0, ba }, 1) };
    complex<double> res { G({ 0, 0, ba, xa }, 1)
        + G({ 0, xa }, 1) * G({ 0, b }, { 1, sb }, x) + Lxa * G00ba
        + G({ xa }, 1) * G({ 0, 0, b }, { 1, 1, sb }, x) + 3. * G({ 0, 0, 0, ba }, 1)
        + G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x) };
    if (b != x) {
        res += G1(b, sb, x) * (G({ 0, 0, xa }, 1) - G00ba);
    }
    return res;
}

// G({0,a,b,c},x) with c smallest, c!=a and c!=b
complex<double> G4_0abc_c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    const complex<double> ac { a / c };
    const complex<double> bc { b / c };
    const complex<double> xc { x / c };
    const complex<double> Lxc { Log(-xc, sc) };
    const complex<double> Gbc { G({ bc }, 1) };
    const complex<double> G0ax { G({ 0, a }, { 1, sa }, x) };

    complex<double> res { G({ 0, 0, bc, ac }, 1) + G({ 0, bc, 0, ac }, 1)
        + G({ bc, 0, 0, ac }, 1) - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
        - 2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x) - G({ bc, ac, 0, xc }, 1)
        + Lxc * (G({ 0, bc, ac }, 1) + G({ bc, 0, ac }, 1))
        + G({ 0, a, b }, { 1, sa, sb }, x) * (Gbc + Lxc)
        - (Lxc * Gbc + G({ 0, bc }, 1)) * G0ax
        + 2. * Gbc * G({ 0, 0, a }, { 1, 1, sa }, x)
        + G({ bc, ac }, 1) * (Lxc * Lxc / 2. + Zeta(2) - G0ax) };
    return res;
}

// G({0,a,b,c},x) with b smallest, a != b
complex<double> G4_0abc_b(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    if (b == c) // take care of G({0,a,b,b},x)
    {
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G4_0abc: b==c but sb!=sc" };
        }

        const complex<double> ab { a / b };
        const complex<double> xb { x / b };
        const complex<double> Lxb { Log(-xb, sb) };
        const complex<double> Gab { G({ ab }, 1) };
        complex<double> res { -G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
            - G({ 0, 0, ab, 1 }, 1) / 2. + G({ ab, 0, 0, 1 }, 1) / 2.
            + G({ ab, 0, 1, xb }, 1) / 2. + G({ ab, 1, 0, xb }, 1) / 2.
            - G({ 0, ab, 1 }, 1) * Lxb / 2.
            + G({ 0, a, b }, { 1, sa, sb }, x) * (Gab + Lxb) / 2.
            + (G({ ab, 0, xb }, 1) - G({ ab, 0, 1 }, 1)) * G1(b, sb, x) / 2.
            - Gab * G({ 0, 0, b }, { 1, 1, sb }, x) / 2.
            + G({ ab, 1 }, 1) * (G({ 0, b }, { 1, sb }, x) - Lxb * Lxb / 2. - Zeta(2))
                / 2. };
        return res;

    } else {
        const complex<double> ab { a / b };
        const complex<double> cb { c / b };
        const complex<double> xb { x / b };
        const complex<double> Lxb { Log(-xb, sb) };
        const complex<double> Gab { G({ ab }, 1) };
        const complex<double> Gcb { G({ cb }, 1) };
        const complex<double> G0ax { G({ 0, a }, { 1, sa }, x) };

        complex<double> res { G({ 0, a, 0, c }, { 1, sa, 1, sc }, x)
            + G({ ab, 0, 0, cb }, 1) - G({ ab, 0, xb, cb }, 1)
            + G({ 0, a, c }, { 1, sa, sc }, x) * (Gab - Gcb) + G({ 0, cb }, 1) * G0ax
            + Gcb
                * (G({ ab, 0, xb }, 1) - G({ 0, 0, ab }, 1)
                    - 2. * G({ 0, 0, a }, { 1, 1, sa }, x)
                    + (G0ax - G({ 0, ab }, 1)) * Lxb)
            - Gab * G({ 0, 0, c }, { 1, 1, sc }, x)
            - G({ cb, ab }, 1) * (Lxb * Lxb / 2. + Zeta(2) - G0ax)
            - G({ ab, cb }, 1) * (Lxb * Lxb / 2. + Zeta(2) - G({ 0, c }, { 1, sc }, x)) };
        if (c != x) {
            res += (G({ ab, 0, xb }, 1) - G({ ab, 0, cb }, 1)) * G1(c, sc, x);
        }
        return res;
    }
}

// G({0,a,b,c},x) with a smallest
complex<double> G4_0abc_a(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    if (a == b) // take care of G({0,a,a,c},x)
    {
        if (sa != sb) {
            throw FastGPL::FastGPL_error { "G4_0abc: a==b but sa!=sb" };
        }

        const complex<double> ca { c / a };
        const complex<double> xa { x / a };
        const complex<double> Lxa { Log(-xa, sa) };
        const complex<double> Gca { G({ ca }, 1) };
        const complex<double> Gax { G1(a, sa, x) };
        const complex<double> Gxa { G({ xa }, 1) };
        const complex<double> G0ax { G({ 0, a }, { 1, sa }, x) };
        const complex<double> G0cx { G({ 0, c }, { 1, sc }, x) };
        const complex<double> G0aca { G({ 0, 1, ca }, 1) };
        const complex<double> G00ca { G({ 0, 0, ca }, 1) };
        complex<double> res {
            (-3. * G({ 0, 0, 1, ca }, 1) - G({ 0, 0, ca, xa }, 1) - G({ 0, 1, 0, ca }, 1)
                - G({ 0, 1, ca, xa }, 1) - G({ 0, ca, 0, xa }, 1) + G({ 0, ca, 1, xa }, 1)
                + G({ ca, 0, 1, xa }, 1) + G({ ca, 1, 0, xa }, 1) - G0aca * Lxa
                + G({ 0, c, a }, { 1, sc, sa }, x) * Gca
                + (-G({ 0, ca, 1 }, 1) + G({ 0, ca, xa }, 1) - G({ ca, 0, 1 }, 1)
                      + G({ ca, 0, xa }, 1))
                    * Gax
                + (G({ 0, a, c }, { 1, sa, sc }, x) - G({ 0, 0, c }, { 1, 1, sc }, x))
                    * (Gax - Gxa)
                - (G({ a, c }, { sa, sc }, x) * (G({ 0, xa }, 1) + Zeta(2)))
                - G({ ca, 1 }, 1) * (Lxa * Lxa / 2. + Zeta(2) - G0ax)
                + Lxa * G({ 0, 0, ca }, 1) + 2. * G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                + G0cx * Zeta(2)
                + G({ 0, ca }, 1)
                    * (Lxa * Lxa / 2. + Zeta(2) + G({ 0, a }, { 1, sa }, x) - G0cx)
                - Gca * (G({ 0, 0, a }, { 1, 1, sa }, x) + Zeta(3)))
            / 2.
        };
        if (c != x) {
            res += G({ c, a }, { sc, sa }, x) * (G({ 0, xa }, 1) - G({ 0, ca }, 1)) / 2.
                + (G1(c, sc, x)
                      * (G0aca - G({ 0, 1, xa }, 1) + G00ca - G({ 0, 0, xa }, 1)))
                    / 2.;
        }
        return res;
    }

    if (a == c) // take care of G({0,a,b,a},x)
    {
        if (sa != sc) {
            throw FastGPL::FastGPL_error { "G4_0abc: a==c but sa!=sc" };
        }

        const complex<double> ba { b / a };
        const complex<double> xa { x / a };
        const complex<double> Lxa { Log(-xa, sa) };
        const complex<double> Gax { G1(a, sa, x) };
        const complex<double> G0ba { G({ 0, ba }, 1) };
        const complex<double> G0ax { G({ 0, a }, { 1, sa }, x) };
        complex<double> res { G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
            - G({ 0, 0, ba, 1 }, 1) - G({ 0, ba, 0, 1 }, 1) - G({ 0, ba, 1, xa }, 1)
            - G({ ba, 0, 0, 1 }, 1) - G({ ba, 0, 1, xa }, 1) - G({ ba, 1, 0, xa }, 1)
            + (G({ 0, ba, 1 }, 1) - G({ 0, ba, xa }, 1) + G({ ba, 0, 1 }, 1)
                  - G({ ba, 0, xa }, 1))
                * Gax
            - G0ba * G0ax
            + G({ ba }, 1)
                * (G({ 0, 0, a }, { 1, 1, sa }, x) - G({ 0, b, a }, { 1, sb, sa }, x))
            + G({ ba, 1 }, 1) * (Lxa * Lxa / 2. + Zeta(2) - G0ax) };
        if (b != x) {
            res += G({ b, a }, { sb, sa }, x) * (G0ba - G({ 0, xa }, 1));
        }
        return res;
    }

    // G({0,a,b,c},x)
    const complex<double> ba { b / a };
    const complex<double> ca { c / a };
    const complex<double> xa { x / a };
    const complex<double> Lxa { Log(-xa, sa) };
    const complex<double> G0cx { G({ 0, c }, { 1, sc }, x) };
    complex<double> res { G({ 0, 0, b, c }, { 1, 1, sb, sc }, x) - G({ 0, 0, ba, ca }, 1)
        - G({ 0, ba, 0, ca }, 1) - G({ 0, ba, ca, xa }, 1) - G({ ba, 0, 0, ca }, 1)
        - G({ ba, 0, ca, xa }, 1) - G({ ba, ca, 0, xa }, 1) - G({ 0, ba }, 1) * G0cx
        + G({ ba }, 1)
            * (G({ 0, 0, c }, { 1, 1, sc }, x) - G({ 0, b, c }, { 1, sb, sc }, x))
        + G({ ba, ca }, 1) * (Lxa * Lxa / 2. + Zeta(2) - G0cx) };
    if (b != x) {
        res += G({ b, c }, { sb, sc }, x) * (G({ 0, ba }, 1) - G({ 0, xa }, 1));
    }
    if (c != x) {
        res += (G({ 0, ba, ca }, 1) - G({ 0, ba, xa }, 1) + G({ ba, 0, ca }, 1)
                   - G({ ba, 0, xa }, 1))
            * G1(c, sc, x);
    }
    return res;
}

// G({a,0,b,c},x) with c smallest, c!=a and c!=b
complex<double> G4_a0bc_c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    const complex<double> ac { a / c };
    const complex<double> bc { b / c };
    const complex<double> xc { x / c };
    const complex<double> Lxc { Log(-xc, sc) };
    const complex<double> G0ax { G({ 0, a }, { 1, sa }, x) };
    const complex<double> Gcx { G1(c, sc, x) };
    const complex<double> Gbc { G({ bc }, 1) };
    const complex<double> Gb0a { G({ bc, 0, ac }, 1) };

    complex<double> res { -G({ 0, b, 0, a }, { 1, sb, 1, sa }, x) - G({ 0, bc, 0, ac }, 1)
        - 2. * G({ bc, 0, 0, ac }, 1) - G({ bc, 0, ac, xc }, 1) - Gb0a * Lxc
        + (-G({ 0, a, b }, { 1, sa, sb }, x) - G({ 0, b, a }, { 1, sb, sa }, x)) * Gcx
        - G({ bc, xc }, 1) * G0ax + G({ c, b }, { sc, sb }, x) * G0ax
        + G1(a, sa, x)
            * (-G({ 0, c, b }, { 1, sc, sb }, x) + Gb0a - G({ bc, 0, xc }, 1)
                - G({ c, 0, b }, { sc, 1, sb }, x) + Gcx * G({ 0, b }, { 1, sb }, x))
        - Gbc * G({ 0, 0, a }, { 1, 1, sa }, x) };
    if (b != x) {
        res += G({ b, 0, a }, { sb, 1, sa }, x) * (Gbc - G({ xc }, 1));
    }
    return res;
}

// G({a,0,b,c},x) with b smallest
complex<double> G4_a0bc_b(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    if (b == c) // take care of G({a,0,b,b},x)
    {
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G4_a0bc: b==c but sb!=sc" };
        }

        const complex<double> ab { a / b };
        const complex<double> xb { x / b };
        const complex<double> Gax { G1(a, sa, x) };
        const complex<double> Lxb { Log(-xb, sb) };
        const complex<double> Gbx { G1(b, sb, x) };
        const complex<double> G0ab { G({ 0, ab }, 1) };
        complex<double> res {
            (G({ 0, 1, 0, ab }, 1) + G({ 0, 1, ab, xb }, 1)
                - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x) - G({ 0, ab, 0, 1 }, 1)
                + G({ 0, ab, 1, xb }, 1) - G({ a, 0, 0, b }, { sa, 1, 1, sb }, x)
                - 2. * G({ ab, 0, 0, 1 }, 1) + G({ 0, 1, ab }, 1) * (Lxb - Gax)
                + G({ 0, ab, 1 }, 1) * (Lxb - Gbx) + G({ 0, ab, xb }, 1) * Gbx
                + G({ a, 0, b }, { sa, 1, sb }, x) * (Gbx - G({ xb }, 1))
                + G({ a, b }, { sa, sb }, x) * (-G0ab - Zeta(2))
                - 2. * G({ ab }, 1) * Zeta(3) + G0ab * G({ 0, b }, { 1, sb }, x)
                - G({ 0, a }, { 1, sa }, x) * Zeta(2)
                + Gax * (Lxb * Zeta(2) + 2. * Zeta(3)))
            / 2.
        };
        return res;
    }

    const complex<double> ab { a / b };
    const complex<double> xb { x / b };
    const complex<double> cb { c / b };
    const complex<double> Gax { G1(a, sa, x) };
    const complex<double> Gcb { G({ cb }, 1) };
    const complex<double> G0ab { G({ 0, ab }, 1) };
    const complex<double> G0cb { G({ 0, cb }, 1) };
    const complex<double> G0c { G({ 0, c }, { 1, sc }, x) };
    const complex<double> G0ax { G({ 0, ab, xb }, 1) };

    complex<double> res { G({ 0, ab, 0, cb }, 1) - G({ 0, ab, xb, cb }, 1)
        + G({ a, 0, 0, c }, { sa, 1, 1, sc }, x)
        + G({ a, c }, { sa, sc }, x) * (G0cb - G0ab) + G0ab * G0c
        + Gcb
            * (G0ax - G({ a, 0, c }, { sa, 1, sc }, x)
                + Log(-xb, sb) * (G0ab - G({ 0, a }, { 1, sa }, x))
                + 2. * G({ 0, 0, ab }, 1) + G({ 0, 0, a }, { 1, 1, sa }, x))
        + Gax
            * (-G({ 0, c, b }, { 1, sc, sb }, x) - G({ 0, cb, ab }, 1)
                + G({ 0, cb, xb }, 1) - G({ b, 0, c }, { sb, 1, sc }, x)
                - G({ cb, 0, ab }, 1) + G({ cb, 0, xb }, 1) + G0c * (Gcb + G1(b, sb, x))
                - G({ 0, 0, c }, { 1, 1, sc }, x)) };

    if (c != x) {
        res += G1(c, sc, x)
            * (-G({ 0, ab, cb }, 1) + G0ax + Gax * (G({ 0, xb }, 1) - G0cb));
    }
    return res;
}

// G({a,0,b,c},x) with a smallest
complex<double> G4_a0bc_a(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    const complex<double> ca { c / a };
    const complex<double> ba { b / a };
    const complex<double> xa { x / a };
    const complex<double> G0bc { G({ 0, ba, ca }, 1) };
    const complex<double> G0ba { G({ 0, ba }, 1) };
    complex<double> res { G({ 0, 0, b, c }, { 1, 1, sb, sc }, x)
        + 2. * G({ 0, 0, ba, ca }, 1) + G({ 0, ba, 0, ca }, 1) + G({ 0, ba, ca, xa }, 1)
        + G0bc * Log(-xa, sa) + G({ 0, b, c }, { 1, sb, sc }, x) * G({ xa }, 1)
        + G0ba * G({ 0, c }, { 1, sc }, x) };
    if (b != x) {
        res += G({ b, c }, { sb, sc }, x) * (-G0ba + G({ 0, xa }, 1));
    }
    if (c != x) {
        res += (-G0bc + G({ 0, ba, xa }, 1)) * G1(c, sc, x);
    }
    return res;
}

//-----------------------------------------------------------------------

// G({a,b,0,c},x) with c smallest, c!=a and c!=b
complex<double> G4_ab0c_c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    const complex<double> ac { a / c };
    const complex<double> bc { b / c };
    const complex<double> xc { x / c };
    const complex<double> Lxc { Log(-xc, sc) };
    const complex<double> Gax { G1(a, sa, x) };
    const complex<double> Gcx { G1(c, sc, x) };
    const complex<double> Gxc { G({ xc }, 1) };
    const complex<double> G0bc { G({ 0, bc }, 1) };
    const complex<double> Gab { G({ a, b }, { sa, sb }, x) };
    complex<double> res { -G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
        - 2. * G({ 0, 0, bc, ac }, 1) - G({ 0, bc, 0, ac }, 1) - G({ 0, bc, ac, xc }, 1)
        + G({ 0, bc, ac }, 1) * (-Lxc + Gax)
        + G({ 0, b, a }, { 1, sb, sa }, x) * (-Gxc + Gcx) + Lxc * Gax * G0bc
        + (-G0bc + G({ 0, xc }, 1)) * Gab - G0bc * G({ 0, a }, { 1, sa }, x)
        + Gab * G({ 0, c }, { 1, sc }, x)
        - Gax
            * (-G({ 0, b }, { 1, sb }, x) * (Gxc - Gcx) - 2. * G({ 0, 0, bc }, 1)
                - G({ 0, 0, b }, { 1, 1, sb }, x)) };

    return res;
}

// G({a,b,0,c},x) with b smallest
complex<double> G4_ab0c_b(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    const complex<double> ab { a / b };
    const complex<double> xb { x / b };
    const complex<double> cb { c / b };
    const complex<double> Gax { G1(a, sa, x) };
    const complex<double> Gab { G({ ab }, 1) };
    const complex<double> G0ab { G({ 0, ab }, 1) };
    const complex<double> G0cb { G({ 0, cb }, 1) };
    complex<double> res { -G({ 0, cb, 0, ab }, 1) + G({ 0, cb, xb, ab }, 1)
        + G({ a, 0, 0, c }, { sa, 1, 1, sc }, x) + G({ 0, cb, ab }, 1) * Gax
        + G({ a, c }, { sa, sc }, x) * (G0ab - G0cb)
        - G0cb * (Log(-xb, sb) * (Gab - Gax) + G({ 0, a }, { 1, sa }, x))
        - (G({ ab, xb }, 1) + G0ab) * G({ 0, c }, { 1, sc }, x)
        - 2. * (Gab - Gax) * G({ 0, 0, cb }, 1)
        + Gab
            * (-G({ 0, cb, xb }, 1) + G({ a, 0, c }, { sa, 1, sc }, x)
                - G({ 0, 0, c }, { 1, 1, sc }, x)) };
    if (c != x) {
        res += (G({ 0, ab, cb }, 1) - G({ 0, ab, xb }, 1) + G({ ab, 0, cb }, 1)
                   - G({ ab, 0, xb }, 1))
            * G1(c, sc, x);
    }
    return res;
}

// G({a,b,0,c},x) with a smallest
complex<double> G4_ab0c_a(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    if (a == b) // take care of G({a,a,0,c},x)
    {
        if (sa != sb) {
            throw FastGPL::FastGPL_error { "G4_ab0c: a==b but sa!=sb" };
        }

        const complex<double> ca { c / a };
        const complex<double> xa { x / a };
        const complex<double> G01c { G({ 0, 1, ca }, 1) };
        const complex<double> Lxa { Log(-xa, sa) };
        const complex<double> Gax { G1(a, sa, x) };
        const complex<double> G0xa { G({ 0, xa }, 1) };
        const complex<double> G0ca { G({ 0, ca }, 1) };
        complex<double> res { -G({ ca }, 1) * G({ 0, 0, 1 }, 1)
            - G({ 0, 0, a, c }, { 1, 1, sa, sc }, x) / 2.
            - G({ 0, 0, c, a }, { 1, 1, sc, sa }, x) / 2. - G({ 0, 1, 0, ca }, 1) / 2.
            - G({ 0, 1, ca, xa }, 1) / 2. + G({ 0, ca, 0, 1 }, 1) / 2.
            - G({ 0, ca, 1, xa }, 1) / 2. + G({ ca, 0, 0, 1 }, 1) + (-G01c * Lxa) / 2.
            + ((-G({ 0, a, c }, { 1, sa, sc }, x) - G({ 0, c, a }, { 1, sc, sa }, x))
                  * G({ xa }, 1))
                / 2.
            + ((-G({ 0, ca, xa }, 1) + G({ a, 0, c }, { sa, 1, sc }, x)) * Gax) / 2.
            + (G({ 0, ca, 1 }, 1) * (-Lxa + Gax)) / 2.
            + (G({ a, c }, { sa, sc }, x) * (-G0xa - Zeta(2))) / 2.
            + (-G0ca * G({ 0, a }, { 1, sa }, x) + G({ 0, c }, { 1, sc }, x) * Zeta(2))
                / 2. };
        if (c != x) {
            res += (G01c - G({ 0, 1, xa }, 1)) * G1(c, sc, x) / 2.
                + G({ c, a }, { sc, sa }, x) * (G0ca - G0xa) / 2.;
        }
        return res;
    }

    const complex<double> ca { c / a };
    const complex<double> ba { b / a };
    const complex<double> xa { x / a };
    const complex<double> Gb0c { G({ ba, 0, ca }, 1) };
    const complex<double> Gba { G({ ba }, 1) };
    complex<double> res { G({ 0, b, 0, c }, { 1, sb, 1, sc }, x) + G({ 0, ba, 0, ca }, 1)
        + 2. * G({ ba, 0, 0, ca }, 1) + G({ ba, 0, ca, xa }, 1) + Gb0c * Log(-xa, sa)
        + G({ ba, xa }, 1) * G({ 0, c }, { 1, sc }, x)
        + Gba * G({ 0, 0, c }, { 1, 1, sc }, x) };
    if (b != x) {
        res += G({ b, 0, c }, { sb, 1, sc }, x) * (-Gba + G({ xa }, 1));
    }
    if (c != x) {
        res += (-Gb0c + G({ ba, 0, xa }, 1)) * G1(c, sc, x);
    }
    return res;
}

complex<double> G4_abcd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{

    if (a == b) {
        if (a == c) {
            const complex<double> da { d / a };
            const complex<double> xa { x / a };
            const complex<double> Gda { G({ da, 1, 1 }, 1) };
            complex<double> res { -(G({ a, d }, { sa, sd }, x) * G({ xa, 1 }, 1))
                + G({ 0, a, a, d }, { 1, sa, sa, sd }, x) + G({ 0, da, 1, 1 }, 1)
                + G({ da, 1, 1, xa }, 1) + G({ da, 1, xa, 1 }, 1) + G({ da, xa, 1, 1 }, 1)
                + Gda * Log(-xa, sa) + G({ a, a, d }, { sa, sa, sd }, x) * G({ xa }, 1) };
            if (d != x) {
                res += (-Gda + G({ xa, 1, 1 }, 1)) * G1(d, sd, x);
            }
            return res;
        }

        const complex<double> ca { c / a };
        const complex<double> da { d / a };
        const complex<double> xa { x / a };
        const complex<double> Gca1 { G({ ca, 1 }, 1) };
        const complex<double> Gcx { G({ ca, xa }, 1) };
        const complex<double> Gc1d { G({ ca, 1, da }, 1) };
        const complex<double> Lxa { Log(-xa, sa) };
        const complex<double> Gca { G({ ca }, 1) };
        const complex<double> Gax { G1(a, sa, x) };
        const complex<double> Gcd { G({ ca, da }, 1) };

        complex<double> res { G({ a, d }, { sa, sd }, x) * (Gca1 / 2. - Gcx / 2.)
            + G({ 0, 1, ca, da }, 1) / 2. - G({ 0, c, a, d }, { 1, sc, sa, sd }, x) / 2.
            - G({ 0, c, d, a }, { 1, sc, sd, sa }, x) / 2. - G({ ca, 1, 0, da }, 1) / 2.
            - G({ ca, 1, da, xa }, 1) / 2. - G({ ca, da, 1, xa }, 1) / 2.
            - Gc1d * Lxa / 2.
            + (-G({ 0, a, d }, { 1, sa, sd }, x) - G({ 0, d, a }, { 1, sd, sa }, x)) * Gca
                / 2.
            + (G({ a, c, d }, { sa, sc, sd }, x) - G({ ca, da, xa }, 1)) * Gax / 2.
            + G({ ca, da, 1 }, 1) * (-Lxa + Gax) / 2.
            - Gca1 * G({ 0, d }, { 1, sd }, x) / 2.
            + Gcd * (-G({ 0, a }, { 1, sa }, x) + Zeta(2)) / 2. };
        if (c != x) {
            res += (G({ c, a, d }, { sc, sa, sd }, x) + G({ c, d, a }, { sc, sd, sa }, x))
                * (Gca - G({ xa }, 1)) / 2.;
        }
        if (d != x) {
            res += (Gcd - Gcx) * G({ d, a }, { sd, sa }, x) / 2.
                + (Gc1d - G({ ca, 1, xa }, 1)) * G1(d, sd, x) / 2.;
        }
        return res;
    } // end of a==b

    const complex<double> ba { b / a };
    const complex<double> ca { c / a };
    const complex<double> da { d / a };
    const complex<double> xa { x / a };
    const complex<double> Gbcd { G({ ba, ca, da }, 1) };
    const complex<double> Gba { G({ ba }, 1) };
    const complex<double> Gbc { G({ ba, ca }, 1) };
    complex<double> res { G({ 0, b, c, d }, { 1, sb, sc, sd }, x)
        + G({ 0, ba, ca, da }, 1) + G({ ba, 0, ca, da }, 1) + G({ ba, ca, 0, da }, 1)
        + G({ ba, ca, da, xa }, 1) + Gbcd * Log(-xa, sa)
        + G({ 0, c, d }, { 1, sc, sd }, x) * Gba + Gbc * G({ 0, d }, { 1, sd }, x) };
    if (b != x) {
        res += G({ b, c, d }, { sb, sc, sd }, x) * (-Gba + G({ xa }, 1));
    }
    if (c != x) {
        res += (-Gbc + G({ ba, xa }, 1)) * G({ c, d }, { sc, sd }, x);
    }
    if (d != x) {
        res += (-Gbcd + G({ ba, ca, xa }, 1)) * G1(d, sd, x);
    }
    return res;
}

complex<double> G4_abcd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{

    if (b == c) {
        if (b == d) {
            const complex<double> ab { a / b };
            const complex<double> xb { x / b };
            const complex<double> Gxb { G({ xb }, 1) };
            const complex<double> Lxb { Log(-xb, sb) };
            const complex<double> Gax { G1(a, sa, x) };
            const complex<double> Gbx { G1(b, sb, x) };
            const complex<double> Gxb1 { G({ xb, 1 }, 1) };

            complex<double> res { -G({ 0, ab, 1, 1 }, 1)
                - G({ 0, b, b, a }, { 1, sb, sb, sa }, x) - G({ ab, 1, 1, xb }, 1)
                - G({ ab, 1, xb, 1 }, 1) - G({ ab, xb, 1, 1 }, 1)
                - G({ b, b, a }, { sb, sb, sa }, x) * Gxb
                + G({ ab, 1, 1 }, 1) * (-Lxb + Gax)
                + (-G({ 0, ab, 1 }, 1) + G({ 0, b, a }, { 1, sb, sa }, x)
                      - G({ ab, 1, xb }, 1) - G({ ab, xb, 1 }, 1)
                      + G({ ab, 1 }, 1) * (-Lxb + Gax))
                    * Gbx
                - (G({ ab, xb }, 1) * Gbx * Gbx) / 2.
                + Gax * (-G({ xb, 1, 1 }, 1) - Gxb1 * Gbx)
                + G({ b, a }, { sb, sa }, x) * (Gxb1 + Gxb * Gbx)
                + Gbx * Gbx
                    * (3. * G({ ab }, 1) * (-Lxb + Gax) + Gax * (-3. * Gxb + Gbx)
                        - 3. * G({ 0, ab }, 1) - 3. * G({ 0, a }, { 1, sa }, x))
                    / 6. };
            return res;
        }

        const complex<double> ab { a / b };
        const complex<double> db { d / b };
        const complex<double> xb { x / b };
        const complex<double> Gxb { G({ xb }, 1) };
        const complex<double> Lxb { Log(-xb, sb) };
        const complex<double> Gax { G1(a, sa, x) };
        const complex<double> Gab { G({ ab, 1 }, 1) };
        const complex<double> Gbd { G({ b, d }, { sb, sd }, x) };
        const complex<double> Gdb { G({ db, 1 }, 1) };
        const complex<double> Gad { G({ ab, 1, db }, 1) };
        const complex<double> Gab1 { G({ ab }, 1) };
        const complex<double> Gdb1 { G({ db }, 1) };
        complex<double> res { (Gab / 2. - G({ ab, xb }, 1) / 2.) * Gbd
            + G({ a, b }, { sa, sb }, x) * (Gbd / 2. + Gdb / 2. - G({ db, xb }, 1) / 2.)
            - G({ 0, a, b, d }, { 1, sa, sb, sd }, x) / 2.
            - G({ 0, a, d, b }, { 1, sa, sd, sb }, x) / 2. - G({ 0, ab, 1, db }, 1) / 2.
            - G({ 0, d, a, b }, { 1, sd, sa, sb }, x) / 2. + G({ 0, db, 1, ab }, 1) / 2.
            - G({ a, 0, d, b }, { sa, 1, sd, sb }, x) / 2. - G({ ab, 0, 1, db }, 1) / 2.
            - G({ ab, 1, 0, db }, 1) / 2. - G({ ab, 1, db, xb }, 1) / 2.
            + G({ db, 0, 1, ab }, 1) / 2. + G({ db, 1, 0, ab }, 1) / 2.
            + G({ db, 1, ab, xb }, 1) / 2. - Gad * Lxb / 2.
            - G({ 0, b, d }, { 1, sb, sd }, x) * Gab1 / 2.
            - G({ 0, a, b }, { 1, sa, sb }, x) * Gdb1 / 2.
            - G({ a, 0, b }, { sa, 1, sb }, x) * Gdb1 / 2.
            + G({ a, b, d }, { sa, sb, sd }, x) * (Gab1 - Gxb) / 2.
            + G({ a, d, b }, { sa, sd, sb }, x) * (Gdb1 - Gxb) / 2.
            + G({ db, 1, ab }, 1) * (Lxb - Gax) / 2. - G({ 0, db, 1 }, 1) * Gax / 2.
            - G({ db, 0, 1 }, 1) * Gax / 2.
            + Gdb * (-Lxb * Gax + G({ 0, a }, { 1, sa }, x)) / 2.
            - Gab * G({ 0, d }, { 1, sd }, x) / 2. };
        if (d != x) {
            res += (G({ d, a, b }, { sd, sa, sb }, x) * (Gdb1 - Gxb)) / 2.
                + (Gad - G({ ab, 1, xb }, 1)) * G1(d, sd, x) / 2.;
        }
        return res;
    } // end of b==c

    if (b == d) {
        const complex<double> ab { a / b };
        const complex<double> cb { c / b };
        const complex<double> xb { x / b };
        const complex<double> Gax { G1(a, sa, x) };
        const complex<double> Gab1 { G({ ab }, 1) };
        const complex<double> Gcb1 { G({ cb, 1 }, 1) };
        const complex<double> Gca { G({ cb, ab }, 1) };
        const complex<double> Gcb { G({ cb }, 1) };
        const complex<double> Gbx { G1(b, sb, x) };
        const complex<double> Gac { G({ ab, cb }, 1) };
        const complex<double> G0bx { G({ 0, b }, { 1, sb }, x) };
        complex<double> res { G({ a, b }, { sa, sb }, x) * (-Gcb1 + Gca)
            + G({ 0, 1, cb, ab }, 1) - G({ 0, ab, cb, 1 }, 1)
            + G({ a, 0, c, b }, { sa, 1, sc, sb }, x) - G({ ab, 0, cb, 1 }, 1)
            - G({ ab, cb, 0, 1 }, 1) - G({ cb, 1, 0, ab }, 1) + G({ cb, 1, xb, ab }, 1)
            - G({ 0, c, b }, { 1, sc, sb }, x) * Gab1 - G({ cb, 1, xb }, 1) * Gab1
            + G({ a, c, b }, { sa, sc, sb }, x) * (Gab1 - Gcb)
            + G({ a, 0, b }, { sa, 1, sb }, x) * Gcb + G({ 0, cb, 1 }, 1) * Gax
            + G({ cb, 0, 1 }, 1) * Gax + G({ cb, 1, ab }, 1) * Gax
            + G({ ab, cb, 1 }, 1) * Gbx - G({ ab, cb, xb }, 1) * Gbx
            + G({ cb, ab, 1 }, 1) * Gbx - G({ cb, ab, xb }, 1) * Gbx
            + Gcb1 * (Log(-xb, sb) * (-Gab1 + Gax) - G({ 0, a }, { 1, sa }, x))
            - Gac * G0bx + Gca * (-G0bx + Zeta(2)) };
        if (c != x) {
            res += (Gac - G({ ab, xb }, 1)) * G({ c, b }, { sc, sb }, x);
        }
        return res;
    } // end of b==d

    const complex<double> ab { a / b };
    const complex<double> cb { c / b };
    const complex<double> db { d / b };
    const complex<double> xb { x / b };
    const complex<double> Gax { G1(a, sa, x) };
    const complex<double> Gab1 { G({ ab }, 1) };
    const complex<double> Gca { G({ cb, ab }, 1) };
    const complex<double> Gcb { G({ cb }, 1) };
    const complex<double> Gac { G({ ab, cb }, 1) };
    const complex<double> Gcd { G({ cb, db }, 1) };
    const complex<double> G0dx { G({ 0, d }, { 1, sd }, x) };
    complex<double> res { G({ a, d }, { sa, sd }, x) * (Gca - Gcd)
        - G({ 0, ab, cb, db }, 1) + G({ 0, db, cb, ab }, 1)
        + G({ a, 0, c, d }, { sa, 1, sc, sd }, x) - G({ ab, 0, cb, db }, 1)
        - G({ ab, cb, 0, db }, 1) - G({ cb, db, 0, ab }, 1) + G({ cb, db, xb, ab }, 1)
        + (-G({ 0, c, d }, { 1, sc, sd }, x) - G({ cb, db, xb }, 1)) * Gab1
        + G({ a, c, d }, { sa, sc, sd }, x) * (Gab1 - Gcb)
        + G({ a, 0, d }, { sa, 1, sd }, x) * Gcb
        + (G({ 0, cb, db }, 1) + G({ cb, 0, db }, 1) + G({ cb, db, ab }, 1)) * Gax
        + Gcd * (Log(-xb, sb) * (-Gab1 + Gax) - G({ 0, a }, { 1, sa }, x)) - Gac * G0dx
        + Gca * (-G({ 0, db }, 1) - G0dx) };
    if (c != x) {
        res += (Gac - G({ ab, xb }, 1)) * G({ c, d }, { sc, sd }, x);
    }
    if (d != x) {
        res += (G({ ab, cb, db }, 1) - G({ ab, cb, xb }, 1) + G({ cb, ab, db }, 1)
                   - G({ cb, ab, xb }, 1))
            * G1(d, sd, x);
    }
    return res;
}

complex<double> G4_abcd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    if (c == d) {
        const complex<double> ac { a / c };
        const complex<double> bc { b / c };
        const complex<double> xc { x / c };
        const complex<double> Gbc1 { G({ bc, 1 }, 1) };
        const complex<double> Gba { G({ bc, ac }, 1) };
        const complex<double> Gbc { G({ bc }, 1) };
        const complex<double> Lxc { Log(-xc, sc) };
        const complex<double> Gax { G1(a, sa, x) };
        const complex<double> Gcx { G1(c, sc, x) };

        const complex<double> Gbax { G({ bc, ac, xc }, 1) };
        const complex<double> G0ax { G({ 0, a }, { 1, sa }, x) };
        complex<double> res { G({ a, c }, { sa, sc }, x) * (Gbc1 / 2. - Gba / 2.)
            + G({ 0, 1 }, 1) * Gba / 2. - G({ 0, 1, bc, ac }, 1) / 2.
            - G({ 0, a, b, c }, { 1, sa, sb, sc }, x) / 2.
            - G({ a, 0, b, c }, { sa, 1, sb, sc }, x) / 2. + G({ bc, 1, 0, ac }, 1) / 2.
            + G({ bc, 1, ac, xc }, 1) / 2. + G({ bc, ac, 1, xc }, 1) / 2.
            - G({ a, 0, c }, { sa, 1, sc }, x) * Gbc / 2.
            + G({ bc, 1, ac }, 1) * (Lxc - Gax) / 2.
            + (-G({ 0, bc, 1 }, 1) - G({ bc, 0, 1 }, 1)) * Gax / 2.
            + G({ bc, ac, 1 }, 1) * (Lxc - Gcx) / 2. + (Gbax * Gcx) / 2.
            + G({ a, b, c }, { sa, sb, sc }, x) * (Gbc - G({ xc }, 1) + Gcx) / 2.
            + (Gbc1 * (-Lxc * Gax + G0ax)) / 2.
            + (Gba * G({ 0, c }, { 1, sc }, x)) / 2. };
        return res;
    }

    const complex<double> ac { a / c };
    const complex<double> bc { b / c };
    const complex<double> dc { d / c };
    const complex<double> xc { x / c };
    const complex<double> Gba { G({ bc, ac }, 1) };
    const complex<double> Gbc { G({ bc }, 1) };
    const complex<double> Lxc { Log(-xc, sc) };
    const complex<double> Gax { G1(a, sa, x) };
    const complex<double> Gbax { G({ bc, ac, xc }, 1) };
    const complex<double> G0ax { G({ 0, a }, { 1, sa }, x) };
    const complex<double> Gdc { G({ dc }, 1) };
    const complex<double> G0dc { G({ 0, dc }, 1) };
    complex<double> res { G({ a, d }, { sa, sd }, x) * (-Gba + G({ bc, dc }, 1))
        + G({ a, b, 0, d }, { sa, sb, 1, sd }, x) + G({ bc, ac, 0, dc }, 1)
        - G({ bc, ac, xc, dc }, 1) + G({ a, b, d }, { sa, sb, sd }, x) * (Gbc - Gdc)
        + (-G({ 0, a, b }, { 1, sa, sb }, x) + G({ 0, bc, ac }, 1)
              - G({ a, 0, b }, { sa, 1, sb }, x) + G({ bc, 0, ac }, 1) + Gbax)
            * Gdc
        + (-G({ bc, dc, ac }, 1) - G({ dc, 0, bc }, 1) - G({ dc, bc, ac }, 1)) * Gax
        + G({ a, b }, { sa, sb }, x) * (G({ dc, bc }, 1) + Lxc * Gdc + G0dc)
        + Gbc
            * (-G({ a, 0, d }, { sa, 1, sd }, x) + Gax * (-G0dc)
                + Gdc * (-Lxc * Gax + G0ax))
        + Gba * (Lxc * Gdc + G({ 0, d }, { 1, sd }, x)) };
    if (d != x) {
        res += (-G({ bc, ac, dc }, 1) + Gbax) * G1(d, sd, x);
    }
    return res;
}

complex<double> G4_abcd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    const complex<double> ad { a / d };
    const complex<double> bd { b / d };
    const complex<double> cd { c / d };
    const complex<double> xd { x / d };
    const complex<double> Gax { G1(a, sa, x) };
    const complex<double> Lxd { Log(-xd, sd) };
    const complex<double> Gcb { G({ cd, bd }, 1) };
    const complex<double> Gcd { G({ cd }, 1) };
    const complex<double> Gabx { G({ a, b }, { sa, sb }, x) };
    const complex<double> Gcx { G({ cd, xd }, 1) };
    complex<double> res { -G({ 0, c, b, a }, { 1, sc, sb, sa }, x)
        - G({ 0, cd, bd, ad }, 1) - G({ cd, 0, bd, ad }, 1) - G({ cd, bd, 0, ad }, 1)
        - G({ cd, bd, ad, xd }, 1)
        + (G({ 0, c, b }, { 1, sc, sb }, x) + G({ 0, cd, bd }, 1) + G({ cd, 0, bd }, 1))
            * Gax
        + G({ cd, bd, ad }, 1) * (-Lxd + Gax)
        + G({ a, b, c }, { sa, sb, sc }, x) * G1(d, sd, x)
        + Gcb * (Lxd * Gax - G({ 0, a }, { 1, sa }, x))
        + Gcd * (-G({ 0, b, a }, { 1, sb, sa }, x) + Gax * G({ 0, b }, { 1, sb }, x))
        + Gabx * (-Gcx - Lxd * Gcd - G({ 0, cd }, 1) - G({ 0, c }, { 1, sc }, x)) };
    if (b != x) {
        res += G({ b, a }, { sb, sa }, x) * (Gcb - Gcx)
            + (-Gcb + Gcx) * Gax * G1(b, sb, x);
    }
    if (c != x) {
        res += G({ c, b, a }, { sc, sb, sa }, x) * (Gcd - G({ xd }, 1))
            + G({ c, b }, { sc, sb }, x) * (-Gcd + G({ xd }, 1)) * Gax
            + Gabx * (Gcd - G({ xd }, 1)) * G1(c, sc, x);
    }
    return res;
}

// G(0,a,b,c,x)
// This should only be called if abs(a) < x and/or abs(b) < x and/or abs(c) < x
complex<double> G4_0abc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G4_0abc_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G4_0abc_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G4_0abc_b(a, b, c, sa, sb, sc, x);
}

// G(a,0,b,c,x)
// This should only be called if abs(a) < x and/or abs(b) < x and/or abs(c) < x
complex<double> G4_a0bc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G4_a0bc_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G4_a0bc_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G4_a0bc_b(a, b, c, sa, sb, sc, x);
}

// G(a,b,0,c,x)
// This should only be called if abs(a) < x and/or abs(b) < x and/or abs(c) < x
complex<double> G4_ab0c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G4_ab0c_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G4_ab0c_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G4_ab0c_b(a, b, c, sa, sb, sc, x);
}

// G(a,b,c,d,x)
// This should only be called if abs(a) < x and/or abs(b) < x and/or abs(c) < x and/or
// abs(d) < x
complex<double> G4_abcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G4_abcd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G4_abcd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G4_abcd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G4_abcd_b(a, b, c, d, sa, sb, sc, sd, x);
}
