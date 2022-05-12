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
//#include<iostream>

// G({0,0,0,a,a},x)
complex<double> G5_explicit_000aa(complex<double> a, int sa, double x)
{
    if (x == a) {
        return -(Zeta(2) * Zeta(3)) + 2. * Zeta(5);
    }

    const complex<double> sy1 = G({ a }, { sa }, x);
    return -(sy1 * G({ 0, 0, 0, x / a }, 1)) + G({ 0, 0, 0, x / a, 1 }, 1)
        + G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
        - G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(2) - sy1 * Zeta(4)
        + G({ 0, a }, { 1, sa }, x) * Zeta(3) + Zeta(5);
}
// G({a,0,0,0,a},x)
complex<double> G5_explicit_a000a(complex<double> a, int sa, double x)
{

    const complex<double> sy1 = G({ a }, { sa }, x);
    return sy1 * G({ 0, 0, 0, x / a }, 1) + G({ 0, 0, 0, 1, x / a }, 1)
        + G({ 0, 0, x / a }, 1) * G({ 0, a }, { 1, sa }, x)
        + G({ 0, x / a }, 1) * G({ 0, 0, a }, { 1, 1, sa }, x)
        + G({ x / a }, 1) * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
        + G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x) + sy1 * Zeta(4)
        + Log(a, sa) * Zeta(4) - Log(-x, sa) * Zeta(4) - 4. * Zeta(5);
}

// G({0,0,a,0,a},x)
complex<double> G5_explicit_00a0a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return 3 * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;

    const complex<double> sy1 = G({ 0, a }, { 1, sa }, x);
    const complex<double> sy2 = G({ a }, { sa }, x);
    const complex<double> sy3 = Log(a, sa);
    const complex<double> sy4 = G({ 0, 0 }, { 1, 1 }, x);

    return sy1 * G({ 0, 0, x / a }, 1) + 3. * sy2 * G({ 0, 0, 0, x / a }, 1)
        + 3. * G({ 0, 0, 0, 1, x / a }, 1) + 2. * G({ 0, 0, 1, 0, x / a }, 1)
        + G({ 0, 1, 0, 0, x / a }, 1) + G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
        - G({ 0, 0, 0 }, { 1, 1, 1 }, x) * Zeta(2)
        + G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(2)
        - (Log(-x, sa) * pow(sy3, 2.) * Zeta(2)) / 2. + (pow(sy3, 3.) * Zeta(2)) / 6.
        - sy3 * Zeta(2) * (sy4 + 2. * Zeta(2))
        + G({ 0 }, { 1 }, x) * Zeta(2) * (sy4 + 2. * Zeta(2)) + 3. * sy2 * Zeta(4)
        - 2. * sy1 * Zeta(3) - 4. * Zeta(5);
}

// G({0,a,0,0,a},x)
complex<double> G5_explicit_0a00a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return -2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;

    const complex<double> sy1 = G({ 0, a }, { 1, sa }, x);
    const complex<double> sy2 = G({ a }, { sa }, x);
    const complex<double> sy3 = Log(a, sa);

    return -2. * sy1 * G({ 0, 0, x / a }, 1) - 3. * sy2 * G({ 0, 0, 0, x / a }, 1)
        - 3. * G({ 0, 0, 0, 1, x / a }, 1) - G({ 0, 0, 1, 0, x / a }, 1)
        - G({ 0, x / a }, 1) * G({ 0, 0, a }, { 1, 1, sa }, x)
        + G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x) - 3. * sy2 * Zeta(4) + sy1 * Zeta(3)
        + G({ 0, 0 }, { 1, 1 }, x) * Zeta(3) + sy3 * Log(-x, sa) * Zeta(3)
        - (pow(sy3, 2.) * Zeta(3)) / 2. + 2. * Zeta(2) * Zeta(3) + 6. * Zeta(5);
}

// G({0,0,a,a,a},x)
complex<double> G5_explicit_00aaa(complex<double> a, int sa, double x)
{
    if (x == a) {
        return Zeta(2) * Zeta(3) - 2. * Zeta(5);
    }

    const complex<double> sy1 = G({ a }, { sa }, x);
    return -(sy1 * G({ 0, 0, x / a, 1 }, 1)) + G({ 0, 0, x / a, 1, 1 }, 1)
        + G({ 0, 0, 0, a, a }, { 1, 1, 1, sa, sa }, x)
        + G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(2)
        - G({ 0, a, a }, { 1, sa, sa }, x) * Zeta(2) + (sy1 * Zeta(4)) / 4.
        - 2. * G({ 0, a }, { 1, sa }, x) * Zeta(3) + Zeta(2) * Zeta(3)
        + G({ a, a }, { sa, sa }, x) * (G({ 0, 0, x / a }, 1) + Zeta(3)) - 2. * Zeta(5);
}

// G({a,0,0,a,a},x)
complex<double> G5_explicit_a00aa(complex<double> a, int sa, double x)
{
    const complex<double> sy1 = G({ a }, { sa }, x);
    return sy1 * G({ 0, 0, 1, x / a }, 1) + G({ 0, 0, 1, 1, x / a }, 1)
        + G({ 0, x / a }, 1) * G({ 0, a, a }, { 1, sa, sa }, x)
        + G({ x / a }, 1) * G({ 0, 0, a, a }, { 1, 1, sa, sa }, x)
        + G({ 0, 0, 0, a, a }, { 1, 1, 1, sa, sa }, x) - (sy1 * Zeta(4)) / 4.
        - (Log(a, sa) * Zeta(4)) / 4. + (Log(-x, sa) * Zeta(4)) / 4.
        - G({ 0, a }, { 1, sa }, x) * Zeta(3) + 3. * Zeta(2) * Zeta(3)
        + G({ a, a }, { sa, sa }, x) * (G({ 0, 0, x / a }, 1) + Zeta(3))
        - (11. * Zeta(5)) / 2. + 3. * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5));
}

// G({a,a,0,0,a},x)
complex<double> G5_explicit_aa00a(complex<double> a, int sa, double x)
{
    const complex<double> sy1 = G({ a }, { sa }, x);
    return -(sy1 * G({ 0, 0, 1, x / a }, 1)) - sy1 * G({ 0, 0, x / a, 1 }, 1)
        - sy1 * G({ 0, 1, 0, x / a }, 1) - 2. * G({ 0, 0, 1, 1, x / a }, 1)
        - G({ 0, 0, 1, x / a, 1 }, 1) - G({ 0, 1, 0, 1, x / a }, 1)
        - (G({ 0, 1, x / a }, 1) + G({ 0, x / a, 1 }, 1)) * G({ 0, a }, { 1, sa }, x)
        - G({ x / a, 1 }, 1) * G({ 0, 0, a }, { 1, 1, sa }, x)
        + G({ x / a }, 1) * G({ a, 0, 0, a }, { sa, 1, 1, sa }, x)
        + G({ 0, a, 0, 0, a }, { 1, sa, 1, 1, sa }, x) + (5. * sy1 * Zeta(4)) / 4.
        + (5. * Log(a, sa) * Zeta(4)) / 4. - (5. * Log(-x, sa) * Zeta(4)) / 4.
        - 3. * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        - 6. * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5))
        - 2. * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.);
}

// G({0,a,0,a,a},x)
complex<double> G5_explicit_0a0aa(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return -3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.;
    }
    complex<double> a = a1 / x1;
    double x = 1.0;
    const complex<double> sy1 = G({ 0, a, a }, { 1, sa, sa }, x);
    const complex<double> sy2 = G({ a }, { sa }, x);
    const complex<double> sy3 = Log(a, sa);
    return -(sy1 * G({ 0, x / a }, 1)) - 2. * sy2 * G({ 0, 0, 1, x / a }, 1)
        - sy2 * G({ 0, 1, 0, x / a }, 1) - 2. * G({ 0, 0, 1, 1, x / a }, 1)
        - G({ 0, 1, 0, 1, x / a }, 1) - G({ 0, 1, 1, 0, x / a }, 1)
        + G({ 0, 0, 0, a, a }, { 1, 1, 1, sa, sa }, x) + sy1 * Zeta(2)
        - G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(2) + (5. * sy2 * Zeta(4)) / 4.
        + 2. * G({ a, a }, { sa, sa }, x) * (-G({ 0, 0, x / a }, 1) - Zeta(3))
        - G({ 0, 0 }, { 1, 1 }, x) * Zeta(3) + G({ 0, a }, { 1, sa }, x) * Zeta(3)
        - sy3 * Log(-x, sa) * Zeta(3) + (pow(sy3, 2.) * Zeta(3)) / 2.
        - 2. * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.) - (9. * Zeta(5)) / 2.
        - 3. * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5));
}

// G({0,a,a,0,a},x)
complex<double> G5_explicit_0aa0a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return 2 * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.;
    }
    complex<double> a = a1 / x1;
    double x = 1.0;
    const complex<double> sy1 = G({ 0, a }, { 1, sa }, x);
    const complex<double> sy2 = G({ a, 0, a }, { sa, 1, sa }, x);
    const complex<double> sy3 = G({ a }, { sa }, x);
    const complex<double> sy4 = Log(a, sa);
    return -(sy2 * G({ 0, x / a }, 1)) + sy1 * G({ 0, x / a, 1 }, 1)
        + 2. * sy3 * G({ 0, 0, 1, x / a }, 1) + 2. * sy3 * G({ 0, 0, x / a, 1 }, 1)
        + sy3 * G({ 0, 1, 0, x / a }, 1) + 4. * G({ 0, 0, 1, 1, x / a }, 1)
        + 2. * G({ 0, 0, 1, x / a, 1 }, 1) + 2. * G({ 0, 1, 0, 1, x / a }, 1)
        + G({ 0, 1, 0, x / a, 1 }, 1) + 2. * G({ 0, 1, 1, 0, x / a }, 1)
        + G({ 0, 0, a, 0, a }, { 1, 1, sa, 1, sa }, x) - sy2 * Zeta(2)
        + G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(2) - (7. * sy3 * Zeta(4)) / 4.
        + 2. * sy1 * Zeta(3) + 2. * sy4 * Log(-x, sa) * Zeta(3) - pow(sy4, 2.) * Zeta(3)
        - 2. * Zeta(2) * Zeta(3)
        + 2. * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2)) * Zeta(3)
        + 2. * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.) + (9. * Zeta(5)) / 2.
        + 6. * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5));
}

// G({a,0,a,0,a},x)
complex<double> G5_explicit_a0a0a(complex<double> a, int sa, double x)
{
    const complex<double> sy1 = G({ a, 0, a }, { sa, 1, sa }, x);
    const complex<double> sy2 = G({ a }, { sa }, x);
    return sy1 * G({ 0, x / a }, 1) + sy2 * G({ 0, 1, 0, x / a }, 1)
        + G({ 0, 1, 0, 1, x / a }, 1) + G({ 0, 1, x / a }, 1) * G({ 0, a }, { 1, sa }, x)
        + G({ x / a }, 1) * G({ 0, a, 0, a }, { 1, sa, 1, sa }, x)
        + G({ 0, 0, a, 0, a }, { 1, 1, sa, 1, sa }, x)
        - (-sy1 + G({ 0, 0, a }, { 1, 1, sa }, x)) * Zeta(2) - (3. * sy2 * Zeta(4)) / 4.
        - (3. * Log(a, sa) * Zeta(4)) / 4. + (3. * Log(-x, sa) * Zeta(4)) / 4.
        + 2. * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        + 2. * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.);
}
// G({0,a,a,a,a},x)
complex<double> G5_explicit_0aaaa(complex<double> a, int sa, double x)
{
    if (x == a) {
        return Zeta(5);
    }

    const complex<double> sy1 = G({ a, a, a }, { sa, sa, sa }, x);
    const complex<double> sy2 = G({ a }, { sa }, x);
    return -(sy1 * G({ 0, x / a }, 1)) - sy2 * G({ 0, x / a, 1, 1 }, 1)
        + G({ 0, x / a, 1, 1, 1 }, 1) + G({ 0, 0, a, a, a }, { 1, 1, sa, sa, sa }, x)
        - sy1 * Zeta(2) + G({ 0, a, a }, { 1, sa, sa }, x) * Zeta(2) - sy2 * Zeta(4)
        + G({ a, a }, { sa, sa }, x) * (G({ 0, x / a, 1 }, 1) - Zeta(3))
        + G({ 0, a }, { 1, sa }, x) * Zeta(3) - Zeta(2) * Zeta(3) + 2. * Zeta(5);
}
// G({a,0,a,a,a},x)
complex<double> G5_explicit_a0aaa(complex<double> a, int sa, double x)
{
    const complex<double> sy1 = G({ a, a, a }, { sa, sa, sa }, x);
    const complex<double> sy2 = G({ a }, { sa }, x);
    return sy1 * G({ 0, x / a }, 1) + sy2 * G({ 0, 1, 1, x / a }, 1)
        + G({ 0, 1, 1, 1, x / a }, 1)
        + G({ x / a }, 1) * G({ 0, a, a, a }, { 1, sa, sa, sa }, x)
        + G({ 0, 0, a, a, a }, { 1, 1, sa, sa, sa }, x) + sy1 * Zeta(2)
        - G({ 0, a, a }, { 1, sa, sa }, x) * Zeta(2) + sy2 * Zeta(4)
        + Log(a, sa) * Zeta(4) - Log(-x, sa) * Zeta(4)
        + G({ a, a }, { sa, sa }, x) * (G({ 0, 1, x / a }, 1) - Zeta(3))
        + G({ 0, a }, { 1, sa }, x) * Zeta(3) - Zeta(2) * Zeta(3)
        + 2. * (Zeta(2) * Zeta(3) - 2. * Zeta(5)) + Zeta(5);
}

// G({a,a,0,a,a},x)
complex<double> G5_explicit_aa0aa(complex<double> a, int sa, double x)
{
    const complex<double> sy1 = G({ a }, { sa }, x);
    return -2. * sy1 * G({ 0, 1, 1, x / a }, 1) - sy1 * G({ 0, 1, x / a, 1 }, 1)
        - 3. * G({ 0, 1, 1, 1, x / a }, 1) - G({ 0, 1, 1, x / a, 1 }, 1)
        - G({ x / a, 1 }, 1) * G({ 0, a, a }, { 1, sa, sa }, x)
        + G({ x / a }, 1) * G({ a, 0, a, a }, { sa, 1, sa, sa }, x)
        + G({ 0, a, 0, a, a }, { 1, sa, 1, sa, sa }, x) - 3. * sy1 * Zeta(4)
        - 3. * Log(a, sa) * Zeta(4) + 3. * Log(-x, sa) * Zeta(4)
        - 2. * G({ 0, a }, { 1, sa }, x) * Zeta(3)
        + G({ a, a }, { sa, sa }, x)
        * (-G({ 0, 1, x / a }, 1) - G({ 0, x / a, 1 }, 1) + 2. * Zeta(3))
        - 2. * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        - 6. * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        - 3. * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.);
}

// G({a,a,a,0,a},x)
complex<double> G5_explicit_aaa0a(complex<double> a, int sa, double x)
{
    const complex<double> sy1 = G({ a }, { sa }, x);
    return sy1 * G({ 0, 1, 1, x / a }, 1) + sy1 * G({ 0, 1, x / a, 1 }, 1)
        + sy1 * G({ 0, x / a, 1, 1 }, 1) + 3. * G({ 0, 1, 1, 1, x / a }, 1)
        + 2. * G({ 0, 1, 1, x / a, 1 }, 1) + G({ 0, 1, x / a, 1, 1 }, 1)
        + G({ x / a, 1, 1 }, 1) * G({ 0, a }, { 1, sa }, x)
        - G({ x / a, 1 }, 1) * G({ a, 0, a }, { sa, 1, sa }, x)
        + G({ x / a }, 1) * G({ a, a, 0, a }, { sa, sa, 1, sa }, x)
        + G({ 0, a, a, 0, a }, { 1, sa, sa, 1, sa }, x) + 3. * sy1 * Zeta(4)
        + 3. * Log(a, sa) * Zeta(4) - 3. * Log(-x, sa) * Zeta(4) + 2. * Zeta(2) * Zeta(3)
        + 6. * (Zeta(2) * Zeta(3) - 2. * Zeta(5)) - (9. * Zeta(5)) / 2.
        + 2. * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.);
}

complex<double> G5_000ab(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G5_00a0b(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G5_0a00b(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G5_a000b(complex<double> a, complex<double> b, int sa, int sb, double x);

complex<double> G5_00abc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G5_00abc_c(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_00abc_b(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_00abc_a(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);

complex<double> G5_0a0bc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G5_0a0bc_c(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_0a0bc_b(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_0a0bc_a(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);

complex<double> G5_0ab0c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G5_0ab0c_c(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_0ab0c_b(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_0ab0c_a(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);

complex<double> G5_a00bc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G5_a00bc_c(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_a00bc_b(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_a00bc_a(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);

complex<double> G5_a0b0c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G5_a0b0c_c(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_a0b0c_b(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_a0b0c_a(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);

complex<double> G5_ab00c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G5_ab00c_c(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_ab00c_b(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);
complex<double> G5_ab00c_a(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x);

complex<double> G5_0abcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_0abcd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_0abcd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_0abcd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_0abcd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G5_a0bcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_a0bcd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_a0bcd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_a0bcd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_a0bcd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G5_ab0cd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_ab0cd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_ab0cd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_ab0cd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_ab0cd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G5_abc0d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_abc0d_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_abc0d_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_abc0d_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G5_abc0d_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G5_abcde(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G5_abcde_e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G5_abcde_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G5_abcde_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G5_abcde_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G5_abcde_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);

complex<double> FastGPL_internal::G5_dispatch(
    const vector<complex<double>>& a, const vector<int>& s, double x)
{
    if (a[0] == 0.0) {
        if (a[1] == 0.0) {
            if (a[2] == 0.0) { // 000ab
                if (a[3] == a[4]) {
                    if (s[3] == s[4])
                        return G5_explicit_000aa(a[4], s[4], x);
                    else
                        throw FastGPL_error { "GPL: a[3]==a[4] but s[3]!=s[4]" };
                } else {
                    // if (is_convergent(a, x)){
                    //     return G_Hoelder(a, s, x);
                    // }
                    return G5_000ab(a[3], a[4], s[3], s[4], x);
                }
            } else if (a[3] == 0.0) { // 00a0b
                if (a[2] == a[4]) {
                    if (s[2] == s[4])
                        return G5_explicit_00a0a(a[4], s[4], x);
                    else
                        throw FastGPL_error { "GPL: a[2]==a[4] but s[2]!=s[4]" };
                } else {
                    // if (is_convergent(a, x)){
                    //     return G_Hoelder(a, s, x);
                    // }
                    return G5_00a0b(a[2], a[4], s[2], s[4], x);
                }
            } else { // 00abc
                if (a[2] == a[3] && a[2] == a[4]) {
                    if (s[2] == s[4] && s[2] == s[4])
                        return G5_explicit_00aaa(a[4], s[4], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[2]==a[3]==a[4] but s[2]!=s[3]!=s[4]"
                        };
                } else {
                    // if (is_convergent(a, x)){
                    //     return G_Hoelder(a, s, x);
                    // }
                    return G5_00abc(a[2], a[3], a[4], s[2], s[3], s[4], x);
                }
            }
        } else if (a[2] == 0.0) {
            if (a[3] == 0.0) { // 0a00b
                if (a[1] == a[4]) {
                    if (s[1] == s[4])
                        return G5_explicit_0a00a(a[4], s[4], x);
                    else
                        throw FastGPL_error { "GPL: a[1]==a[4] but s[1]!=s[4]" };
                } else {
                    // if (is_convergent(a, x)){
                    //     return G_Hoelder(a, s, x);
                    // }
                    return G5_0a00b(a[1], a[4], s[1], s[4], x);
                }
            } else { // 0a0bc
                if (a[1] == a[3] && a[1] == a[4]) {
                    if (s[1] == s[3] && s[1] == s[4])
                        return G5_explicit_0a0aa(a[4], s[4], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[1]==a[3]==a[4] but s[1]!=s[3]!=s[4]"
                        };
                } else {

                    // if (is_convergent(a, x)){
                    //     return G_Hoelder(a, s, x);
                    // }
                    return G5_0a0bc(a[1], a[3], a[4], s[1], s[3], s[4], x);
                }
            }
        } else if (a[3] == 0.0) { // 0ab0c
            if (a[1] == a[2] && a[1] == a[4]) {
                if (s[1] == s[2] && s[1] == s[4])
                    return G5_explicit_0aa0a(a[4], s[4], x);
                else
                    throw FastGPL_error { "GPL: a[1]==a[2]==a[4] but s[1]!=s[2]!=s[4]" };
            } else {

                // if (is_convergent(a, x)){
                //     return G_Hoelder(a, s, x);
                // }
                return G5_0ab0c(a[1], a[2], a[4], s[1], s[2], s[4], x);
            }
        } else { // 0abcd
            if (a[1] == a[2] && a[1] == a[3] && a[1] == a[4]) {
                if (s[1] == s[2] && s[1] == s[3] && s[1] == s[4])
                    return G5_explicit_0aaaa(a[4], s[4], x);
                else
                    throw FastGPL_error {
                        "GPL: a[1]=a[2]==a[3]==a[4] but s[1]!=s[2]!=s[3]!=s[4]"
                    };
            } else {

                // if (is_convergent(a, x)){
                //     return G_Hoelder(a, s, x);
                // }
                return G5_0abcd(a[1], a[2], a[3], a[4], s[1], s[2], s[3], s[4], x);
            }
        }
    }

    if (a[1] == 0.0) {
        if (a[2] == 0.0) {
            if (a[3] == 0.0) { // a000b
                if (a[0] == a[4]) {
                    if (s[0] == s[4])
                        return G5_explicit_a000a(a[4], s[4], x);
                    else
                        throw FastGPL_error { "GPL: a[0]==a[4] but s[0]!=s[4]" };
                } else {
                    // if (is_convergent(a, x)){
                    //     return G_Hoelder(a, s, x);
                    // }
                    return G5_a000b(a[0], a[4], s[0], s[4], x);
                }
            } else { // a00bc
                if (a[0] == a[3] && a[0] == a[4]) {
                    if (s[0] == s[3] && s[0] == s[4])
                        return G5_explicit_a00aa(a[4], s[4], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[0]==a[3]==a[4] but s[0]!=s[3]!=s[4]"
                        };
                } else {

                    // if (is_convergent(a, x)){
                    //     return G_Hoelder(a, s, x);
                    // }
                    return G5_a00bc(a[0], a[3], a[4], s[0], s[3], s[4], x);
                }
            }
        } else if (a[3] == 0.0) { // a0b0c
            if (a[0] == a[2] && a[0] == a[4]) {
                if (s[0] == s[2] && s[0] == s[4])
                    return G5_explicit_a0a0a(a[4], s[4], x);
                else
                    throw FastGPL_error { "GPL: a[0]==a[2]==a[4] but s[0]!=s[2]!=s[4]" };
            } else {
                // if (is_convergent(a, x)){
                //     return G_Hoelder(a, s, x);
                // }
                return G5_a0b0c(a[0], a[2], a[4], s[0], s[2], s[4], x);
            }
        } else { // a0bcd
            if (a[0] == a[2] && a[0] == a[3] && a[0] == a[4]) {
                if (s[0] == s[2] && s[0] == s[3] && s[0] == s[4])
                    return G5_explicit_a0aaa(a[4], s[4], x);
                else
                    throw FastGPL_error {
                        "GPL: a[0]==a[2]==a[3]==a[4] but s[0]!=s[2]!=s[3]!=s[4]"
                    };
            } else {
                // if (is_convergent(a, x)){
                //     return G_Hoelder(a, s, x);
                // }
                return G5_a0bcd(a[0], a[2], a[3], a[4], s[0], s[2], s[3], s[4], x);
            }
        }
    }

    if (a[2] == 0.0) {
        if (a[3] == 0.0) { // ab00c
            if (a[0] == a[1] && a[0] == a[4]) {
                if (s[0] == s[1] && s[0] == s[4])
                    return G5_explicit_aa00a(a[4], s[4], x);
                else
                    throw FastGPL_error { "GPL: a[0]==a[1]==a[4] but s[0]!=s[1]!=s[4]" };
            } else {
                // if (is_convergent(a, x)){
                //     return G_Hoelder(a, s, x);
                // }
                return G5_ab00c(a[0], a[1], a[4], s[0], s[1], s[4], x);
            }
        } else { // ab0cd
            if (a[0] == a[1] && a[0] == a[3] && a[0] == a[4]) {
                if (s[0] == s[1] && s[0] == s[3] && s[0] == s[4])
                    return G5_explicit_aa0aa(a[4], s[4], x);
                else
                    throw FastGPL_error {
                        "GPL: a[0]==a[1]==a[3]==a[4] but s[0]!=s[1]!=s[3]!=s[4]"
                    };
            } else {
                // if (is_convergent(a, x)){
                //     return G_Hoelder(a, s, x);
                // }
                return G5_ab0cd(a[0], a[1], a[3], a[4], s[0], s[1], s[3], s[4], x);
            }
        }
    }

    if (a[3] == 0.0) { // abc0d
        if (a[0] == a[1] && a[0] == a[2] && a[0] == a[4]) {
            if (s[0] == s[1] && s[0] == s[2] && s[0] == s[4])
                return G5_explicit_aaa0a(a[4], s[4], x);
            else
                throw FastGPL_error {
                    "GPL: a[0]==a[1]==a[2]==a[4] but s[0]!=s[1]!=s[2]!=s[4]"
                };
        } else {
            // if (is_convergent(a, x)){
            //     return G_Hoelder(a, s, x);
            // }
            return G5_abc0d(a[0], a[1], a[2], a[4], s[0], s[1], s[2], s[4], x);
        }
    } else { // abcde
        return G5_abcde(a[0], a[1], a[2], a[3], a[4], s[0], s[1], s[2], s[3], s[4], x);
    }
}
