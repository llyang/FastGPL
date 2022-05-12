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

// take care of G({0,0,0,a,b},x)
// take care of G({0,0,a,0,b},x)
// take care of G({0,a,0,0,b},x)
// take care of G({a,0,0,0,b},x)

complex<double> G5_000ab(
    complex<double> a1, complex<double> b1, int sa, int sb, double x1)
{

    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    double x = 1.;
    if (abs(b) < abs(a)) { // b smallest
        const complex<double> sy1 = Log(b, sb);
        const complex<double> sy2 = G({ a / b }, 1);
        const complex<double> sy3 = G({ 0, a / b }, 1);
        const complex<double> sy4 = G({ 0, 0, 0 }, { 1, 1, 1 }, x);
        const complex<double> sy5 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy6 = G({ 0 }, { 1 }, x);
        const complex<double> sy7 = G({ 0, 0, a / b }, 1);
        const complex<double> sy8 = G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x);
        const complex<double> sy9 = G({ 0, 0, 0, a / b }, 1);

        complex<double> res = -(sy3 * sy4) + sy5 * (sy3 * sy6 + sy7)
            - G({ 0, 0, 0, 0, a / b }, 1) + G({ a / b, 0, 0, 0, x / b }, 1)
            + sy8 * (-G({ x / b }, 1) + G({ b }, { sb }, x))
            - 4. * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x) + (sy3 * pow(sy1, 3.)) / 6.
            + Log(-x, sb)
                * (sy1 * sy7 - sy9 - (sy3 * pow(sy1, 2.)) / 2.
                    + (sy2 * pow(sy1, 3.)) / 6.)
            - (sy2 * pow(sy1, 4.)) / 24. + 2. * sy3 * sy6 * Zeta(2) + 2. * sy7 * Zeta(2)
            + pow(sy1, 2.) * (-sy7 / 2. + sy2 * (sy5 / 2. + Zeta(2)))
            + sy1
                * (-(sy3 * sy5) + sy9 - 2. * sy3 * Zeta(2)
                    + sy2 * (sy4 - sy5 * sy6 - 2. * sy6 * Zeta(2)))
            + sy2
                * (-(sy4 * sy6) + sy8 + G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x)
                    - pow(sy5, 2.) + sy5 * (pow(sy6, 2.) - 2. * Zeta(2))
                    + 2. * pow(sy6, 2.) * Zeta(2) + 2. * Zeta(4));
        return res;
    }
    // abs(a) <= abs(b)
    const complex<double> sy1 = Log(a, sa);
    const complex<double> sy2 = G({ b / a }, 1);
    const complex<double> sy3 = G({ 0, 0 }, { 1, 1 }, x);
    const complex<double> sy4 = G({ 0 }, { 1 }, x);
    const complex<double> sy5 = G({ 0, 0, 0 }, { 1, 1, 1 }, x);
    complex<double> res { -G({ 0, 0, 0, 0, b / a }, 1) - G({ 0, 0, 0, b / a, x / a }, 1)
        - G({ 0, 0, b / a, 0, x / a }, 1) - G({ 0, b / a, 0, 0, x / a }, 1)
        - G({ b / a, 0, 0, 0, x / a }, 1)
        - G({ 0, 0, b / a }, 1) * G({ 0, b }, { 1, sb }, x)
        + G({ 0, b / a }, 1) * G({ 0, 0, b }, { 1, 1, sb }, x)
        + G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x)
        - (sy2 * Log(-x, sa) * pow(sy1, 3.)) / 6. + (sy2 * pow(sy1, 4.)) / 24.
        + sy2 * pow(sy1, 2.) * (-sy3 / 2. - Zeta(2))
        + sy1 * sy2 * (sy3 * sy4 - sy5 + 2. * sy4 * Zeta(2))
        + sy2
            * (sy4 * sy5 - G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x)
                - G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x) + pow(sy3, 2.)
                - 2. * pow(sy4, 2.) * Zeta(2) + sy3 * (-pow(sy4, 2.) + 2. * Zeta(2))
                - 2. * Zeta(4)) };
    if (b != x) {
        res += (G({ 0, 0, 0, b / a }, 1) - G({ 0, 0, 0, x / a }, 1))
            * G({ b }, { sb }, x);
    }
    return res;
}
complex<double> G5_00a0b(
    complex<double> a1, complex<double> b1, int sa, int sb, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    double x = 1.;
    if (abs(b) < abs(a)) { // b smallest

        const complex<double> sy1 = Log(b, sb);
        const complex<double> sy2 = G({ 0, a / b }, 1);
        const complex<double> sy3 = G({ 0, 0, a }, { 1, 1, sa }, x);
        const complex<double> sy4 = G({ 0, 0, a / b }, 1);
        const complex<double> sy5 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy6 = G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x);
        const complex<double> sy7 = G({ 0, 0, 0, a / b }, 1);
        const complex<double> sy8 = G({ 0 }, { 1 }, x);
        complex<double> res = -2. * sy4 * sy5 + 3. * sy6 * G({ x / b }, 1)
            + 4. * G({ 0, 0, 0, 0, a / b }, 1) + G({ 0, a / b, 0, 0, x / b }, 1)
            - 3. * sy6 * G({ b }, { sb }, x)
            + sy3 * (G({ 0, x / b }, 1) + G({ 0, b }, { 1, sb }, x))
            + 6. * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x) + sy4 * pow(sy1, 2.)
            + Log(-x, sb) * (-2. * sy1 * sy4 + 3. * sy7 + (sy2 * pow(sy1, 2.)) / 2.)
            - (sy2 * pow(sy1, 3.)) / 6. - 4. * sy4 * Zeta(2)
            + sy2
                * (-sy3 - sy5 * sy8 + G({ 0, 0, 0 }, { 1, 1, 1 }, x) - 2. * sy8 * Zeta(2))
            + sy1 * (-3. * sy7 + sy2 * (sy5 + 2. * Zeta(2)));
        return res;
    }
    // abs(a) <= abs(b)
    const complex<double> sy1 = Log(a, sa);
    const complex<double> sy2 = G({ 0, b / a }, 1);
    const complex<double> sy3 = G({ 0, 0 }, { 1, 1 }, x);
    const complex<double> sy4 = G({ 0 }, { 1 }, x);
    complex<double> res { 4. * G({ 0, 0, 0, 0, b / a }, 1)
        + 3. * G({ 0, 0, 0, b / a, x / a }, 1) + 2. * G({ 0, 0, b / a, 0, x / a }, 1)
        + G({ 0, b / a, 0, 0, x / a }, 1)
        + (2. * G({ 0, 0, b / a }, 1) + G({ 0, 0, x / a }, 1)) * G({ 0, b }, { 1, sb }, x)
        + G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x)
        + (sy2 * Log(-x, sa) * pow(sy1, 2.)) / 2. - (sy2 * pow(sy1, 3.)) / 6.
        + sy1 * sy2 * (sy3 + 2. * Zeta(2))
        + sy2
            * (-(sy3 * sy4) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                - G({ 0, 0, b }, { 1, 1, sb }, x) - 2. * sy4 * Zeta(2)) };
    if (b != x) {
        res += (-3. * G({ 0, 0, 0, b / a }, 1) + 3. * G({ 0, 0, 0, x / a }, 1))
            * G({ b }, { sb }, x);
    }
    return res;
}
complex<double> G5_0a00b(
    complex<double> a1, complex<double> b1, int sa, int sb, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    double x = 1.;
    if (abs(b) < abs(a)) { // b smallest
        const complex<double> sy1 = G({ 0, 0, a }, { 1, 1, sa }, x);
        const complex<double> sy2 = Log(b, sb);
        const complex<double> sy3 = G({ 0, 0, a / b }, 1);
        const complex<double> sy4 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy5 = G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x);
        const complex<double> sy6 = G({ 0, 0, 0, a / b }, 1);

        complex<double> res = 3. * sy2 * sy6 - 3. * sy5 * G({ x / b }, 1)
            - 2. * sy1 * G({ 0, x / b }, 1) - 6. * G({ 0, 0, 0, 0, a / b }, 1)
            + G({ 0, 0, a / b, 0, x / b }, 1) + 3. * sy5 * G({ b }, { sb }, x)
            - 2. * sy1 * G({ 0, b }, { 1, sb }, x)
            + sy4 * (-G({ 0, 0, x / b }, 1) + G({ 0, 0, b }, { 1, 1, sb }, x))
            - 4. * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
            + (sy2 * sy3 - 3. * sy6) * Log(-x, sb) - (sy3 * pow(sy2, 2.)) / 2.
            + sy3 * (sy4 + G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2));
        return res;
    }
    // abs(a) <= abs(b)
    const complex<double> sy1 = Log(a, sa);
    const complex<double> sy2 = G({ 0, 0, b / a }, 1);
    const complex<double> sy3 = G({ 0, b }, { 1, sb }, x);
    complex<double> res { -2. * sy3 * G({ 0, 0, x / a }, 1)
        - 6. * G({ 0, 0, 0, 0, b / a }, 1) - 3. * G({ 0, 0, 0, b / a, x / a }, 1)
        - G({ 0, 0, b / a, 0, x / a }, 1)
        - G({ 0, x / a }, 1) * G({ 0, 0, b }, { 1, 1, sb }, x)
        + G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x) - sy1 * sy2 * Log(-x, sa)
        + (sy2 * pow(sy1, 2.)) / 2.
        + sy2 * (-sy3 - G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
    if (b != x) {
        res += (3. * G({ 0, 0, 0, b / a }, 1) - 3. * G({ 0, 0, 0, x / a }, 1))
            * G({ b }, { sb }, x);
    }
    return res;
}
complex<double> G5_a000b(complex<double> a, complex<double> b, int sa, int sb, double x)
{

    if (abs(b) < abs(a)) { // b smallest
        const complex<double> sy1 = G({ 0, 0, a }, { 1, 1, sa }, x);
        const complex<double> sy2 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy3 = G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x);
        const complex<double> sy4 = G({ 0, 0, 0, a / b }, 1);
        const complex<double> sy5 = G({ a }, { sa }, x);

        complex<double> res = -(sy4 * sy5) + sy3 * G({ x / b }, 1)
            + sy1 * G({ 0, x / b }, 1) + sy2 * G({ 0, 0, x / b }, 1)
            + sy5 * G({ 0, 0, 0, x / b }, 1) + 4. * G({ 0, 0, 0, 0, a / b }, 1)
            + G({ 0, 0, 0, a / b, x / b }, 1) - sy3 * G({ b }, { sb }, x)
            + sy1 * G({ 0, b }, { 1, sb }, x) - sy2 * G({ 0, 0, b }, { 1, 1, sb }, x)
            + sy5 * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
            + G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x) - sy4 * Log(b, sb)
            + sy4 * Log(-x, sb);
        return res;
    }
    // abs(a) <= abs(b)
    const complex<double> sy1 = G({ 0, 0, 0, b / a }, 1);
    complex<double> res { 4. * G({ 0, 0, 0, 0, b / a }, 1)
        + G({ 0, 0, 0, b / a, x / a }, 1)
        + G({ 0, 0, x / a }, 1) * G({ 0, b }, { 1, sb }, x)
        + G({ 0, x / a }, 1) * G({ 0, 0, b }, { 1, 1, sb }, x)
        + G({ x / a }, 1) * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
        + G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x) - sy1 * Log(a, sa)
        + sy1 * Log(-x, sa) };
    if (b != x) {
        res += (-sy1 + G({ 0, 0, 0, x / a }, 1)) * G({ b }, { sb }, x);
    }
    return res;
}

complex<double> G5_00abc_c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;

    const complex<double> sy1 = Log(c, sc);
    const complex<double> sy2 = G({ b / c, a / c }, 1);
    const complex<double> sy3 = G({ 0, 0, a }, { 1, 1, sa }, x);
    const complex<double> sy4 = G({ 0, 0 }, { 1, 1 }, x);
    const complex<double> sy5 = G({ 0, b / c, a / c }, 1);
    const complex<double> sy6 = G({ b / c, 0, a / c }, 1);
    const complex<double> sy7 = G({ 0, 0, a, b }, { 1, 1, sa, sb }, x);
    const complex<double> sy8 = G({ b / c }, 1);
    const complex<double> sy9 = G({ 0, 0, b / c, a / c }, 1);
    const complex<double> sy10 = G({ 0, b / c, 0, a / c }, 1);
    const complex<double> sy11 = G({ b / c, 0, 0, a / c }, 1);
    const complex<double> sy12 = G({ 0 }, { 1 }, x);

    complex<double> res { sy4 * (-sy5 - sy6) - sy3 * G({ 0, b / c }, 1)
        + G({ 0, 0, 0, b / c, a / c }, 1) + G({ 0, 0, b / c, 0, a / c }, 1)
        + G({ 0, b / c, 0, 0, a / c }, 1) + G({ b / c, 0, 0, 0, a / c }, 1)
        + G({ b / c, a / c, 0, 0, x / c }, 1)
        + sy7 * (-G({ x / c }, 1) + G({ c }, { sc }, x))
        + sy8 * (sy7 + 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
        - 3. * G({ 0, 0, 0, a, b }, { 1, 1, 1, sa, sb }, x)
        - G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x)
        + (sy5 / 2. + sy6 / 2.) * pow(sy1, 2.)
        + Log(-x, sc)
            * (sy10 + sy11 + sy1 * (-sy5 - sy6) - sy3 * sy8 + sy9
                + (sy2 * pow(sy1, 2.)) / 2.)
        - (sy2 * pow(sy1, 3.)) / 6. - 2. * sy5 * Zeta(2) - 2. * sy6 * Zeta(2)
        + sy2 * (-sy3 - sy12 * sy4 + G({ 0, 0, 0 }, { 1, 1, 1 }, x) - 2. * sy12 * Zeta(2))
        + sy1 * (-sy10 - sy11 + sy3 * sy8 - sy9 + sy2 * (sy4 + 2. * Zeta(2))) };
    return res;
}
complex<double> G5_00abc_b(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (b == c) {
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G5_00abc: b==c but sb!=sc" };
        }
        const complex<double> sy1 = G({ 0, 0, a }, { 1, 1, sa }, x);
        complex<double> res { -G({ a / b, 0, 0, 0, 1 }, 1)
            + G({ a / b, 0, 0, x / b, 1 }, 1)
            + (G({ a / b, 0, 0, 1 }, 1) - G({ a / b, 0, 0, x / b }, 1))
                * G({ b }, { sb }, x)
            - G({ a / b, 0, 1 }, 1) * G({ 0, b }, { 1, sb }, x)
            + G({ a / b, 1 }, 1) * (-sy1 + G({ 0, 0, b }, { 1, 1, sb }, x))
            + G({ a / b }, 1)
                * (-G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
                    + G({ 0, 0, a, b }, { 1, 1, sa, sb }, x))
            + G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x) - sy1 * Zeta(2) };
        return res;

    } else {
        const complex<double> sy1 = Log(b, sb);
        const complex<double> sy2 = G({ a / b, c / b }, 1);
        const complex<double> sy3 = G({ c / b, a / b }, 1);
        const complex<double> sy4 = G({ 0, 0, a }, { 1, 1, sa }, x);
        const complex<double> sy5 = G({ 0, 0, 0 }, { 1, 1, 1 }, x);
        const complex<double> sy6 = G({ 0, a / b, c / b }, 1);
        const complex<double> sy7 = G({ 0, c / b, a / b }, 1);
        const complex<double> sy8 = G({ c / b, 0, a / b }, 1);
        const complex<double> sy9 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy10 = G({ 0 }, { 1 }, x);
        const complex<double> sy11 = G({ c / b }, 1);
        const complex<double> sy12 = G({ 0, 0, a, c }, { 1, 1, sa, sc }, x);
        const complex<double> sy13 = G({ 0, 0, a / b, c / b }, 1);
        const complex<double> sy14 = G({ 0, 0, c / b, a / b }, 1);
        const complex<double> sy15 = G({ 0, c / b, 0, a / b }, 1);
        const complex<double> sy16 = G({ c / b, 0, 0, a / b }, 1);

        complex<double> res { sy3 * (sy4 - sy5) + (sy10 * sy3 + sy6 + sy7 + sy8) * sy9
            + sy4 * G({ 0, c / b }, 1) - G({ 0, 0, 0, a / b, c / b }, 1)
            - G({ 0, 0, 0, c / b, a / b }, 1) - G({ 0, 0, c / b, 0, a / b }, 1)
            - G({ 0, c / b, 0, 0, a / b }, 1) - G({ a / b, 0, 0, 0, c / b }, 1)
            - G({ a / b, 0, 0, c / b, x / b }, 1) - G({ a / b, 0, c / b, 0, x / b }, 1)
            - G({ a / b, c / b, 0, 0, x / b }, 1) - G({ c / b, 0, 0, 0, a / b }, 1)
            - G({ c / b, a / b, 0, 0, x / b }, 1)
            - G({ a / b, 0, c / b }, 1) * G({ 0, c }, { 1, sc }, x)
            + sy11 * (-sy12 - 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
            + G({ a / b }, 1) * (sy12 - G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x))
            + G({ 0, 0, a, 0, c }, { 1, 1, sa, 1, sc }, x)
            + (-sy6 / 2. - sy7 / 2. - sy8 / 2.) * pow(sy1, 2.)
            + Log(-x, sb)
                * (-sy13 - sy14 - sy15 - sy16 + sy11 * sy4 + sy1 * (sy6 + sy7 + sy8)
                    + (-sy2 / 2. - sy3 / 2.) * pow(sy1, 2.))
            + (sy2 / 6. + sy3 / 6.) * pow(sy1, 3.) + 2. * sy10 * sy3 * Zeta(2)
            + 2. * sy6 * Zeta(2) + 2. * sy7 * Zeta(2) + 2. * sy8 * Zeta(2)
            + sy2
                * (-sy5 + sy10 * sy9 + G({ 0, 0, c }, { 1, 1, sc }, x)
                    + 2. * sy10 * Zeta(2))
            + sy1
                * (sy13 + sy14 + sy15 + sy16 - sy11 * sy4 - sy3 * sy9
                    + sy2 * (-sy9 - 2. * Zeta(2)) - 2. * sy3 * Zeta(2)) };
        if (c != x) {
            res += (G({ a / b, 0, 0, c / b }, 1) - G({ a / b, 0, 0, x / b }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G5_00abc_a(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (a == b) {
        if (sa != sb) {
            throw FastGPL::FastGPL_error { "G5_00abc: a==b but sa!=sb" };
        }
        const complex<double> sy1 = Log(a, sa);
        const complex<double> sy2 = G({ c / a, 1 }, 1);
        const complex<double> sy3 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy4 = G({ 0, 0, c }, { 1, 1, sc }, x);
        const complex<double> sy5 = G({ 0 }, { 1 }, x);
        complex<double> res { -G({ 0, 0, 0, c / a, 1 }, 1)
            - G({ 0, 0, c / a, 1, x / a }, 1) - G({ 0, 0, c / a, x / a, 1 }, 1)
            - G({ 0, c / a, 0, 1, x / a }, 1) - G({ 0, c / a, 0, x / a, 1 }, 1)
            - G({ 0, c / a, 1, 0, x / a }, 1) - G({ c / a, 0, 0, 1, x / a }, 1)
            - G({ c / a, 0, 0, x / a, 1 }, 1) - G({ c / a, 0, 1, 0, x / a }, 1)
            - G({ c / a, 1, 0, 0, x / a }, 1)
            + G({ 0, 0, 0, a, c }, { 1, 1, 1, sa, sc }, x)
            - (sy2 * Log(-x, sa) * pow(sy1, 2.)) / 2. + (sy2 * pow(sy1, 3.)) / 6.
            + sy1 * sy2 * (-sy3 - 2. * Zeta(2)) + sy4 * Zeta(2)
            - G({ 0, a, c }, { 1, sa, sc }, x) * Zeta(2)
            + sy2
                * (sy4 + sy3 * sy5 - G({ 0, 0, 0 }, { 1, 1, 1 }, x) + 2. * sy5 * Zeta(2))
            + G({ 0, c }, { 1, sc }, x) * (-G({ 0, c / a, 1 }, 1) - Zeta(3))
            + G({ a, c }, { sa, sc }, x) * (G({ 0, 0, x / a }, 1) + Zeta(3)) };
        if (c != x) {
            res += (G({ 0, 0, c / a, 1 }, 1) - G({ 0, 0, x / a, 1 }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    if (a == c) {
        if (sa != sc) {
            throw FastGPL::FastGPL_error { "G5_00abc: a==c but sa!=sc" };
        }
        const complex<double> sy1 = Log(a, sa);
        const complex<double> sy2 = G({ b / a, 1 }, 1);
        const complex<double> sy3 = G({ 0, 0, a }, { 1, 1, sa }, x);
        const complex<double> sy4 = G({ 0, 0, b / a }, 1);
        const complex<double> sy5 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy6 = G({ 0 }, { 1 }, x);
        complex<double> res { G({ 0, 0, 0, b / a, 1 }, 1) + G({ 0, 0, b / a, 0, 1 }, 1)
            + G({ 0, 0, b / a, 1, x / a }, 1) + G({ 0, b / a, 0, 0, 1 }, 1)
            + G({ 0, b / a, 0, 1, x / a }, 1) + G({ 0, b / a, 1, 0, x / a }, 1)
            + G({ b / a, 0, 0, 0, 1 }, 1) + G({ b / a, 0, 0, 1, x / a }, 1)
            + G({ b / a, 0, 1, 0, x / a }, 1) + G({ b / a, 1, 0, 0, x / a }, 1)
            + (-G({ 0, 0, b / a, 1 }, 1) + G({ 0, 0, b / a, x / a }, 1)
                  - G({ 0, b / a, 0, 1 }, 1) + G({ 0, b / a, 0, x / a }, 1)
                  - G({ b / a, 0, 0, 1 }, 1) + G({ b / a, 0, 0, x / a }, 1))
                * G({ a }, { sa }, x)
            + (sy4 + G({ 0, b / a, 1 }, 1) + G({ b / a, 0, 1 }, 1))
                * G({ 0, a }, { 1, sa }, x)
            + G({ 0, b / a }, 1) * (-sy3 + G({ 0, b, a }, { 1, sb, sa }, x))
            + G({ b / a }, 1)
                * (G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                    - G({ 0, 0, b, a }, { 1, 1, sb, sa }, x))
            + G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x)
            + (sy2 * Log(-x, sa) * pow(sy1, 2.)) / 2. - (sy2 * pow(sy1, 3.)) / 6.
            + sy1 * sy2 * (sy5 + 2. * Zeta(2))
            + sy2
                * (-sy3 - sy5 * sy6 + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                    - 2. * sy6 * Zeta(2)) };
        if (b != x) {
            res += (-sy4 + G({ 0, 0, x / a }, 1)) * G({ b, a }, { sb, sa }, x);
        }
        return res;
    }
    // G({0,a,b,c},x)
    const complex<double> sy1 = Log(a, sa);
    const complex<double> sy2 = G({ b / a, c / a }, 1);
    const complex<double> sy3 = G({ 0, 0, c }, { 1, 1, sc }, x);
    const complex<double> sy4 = G({ 0, 0, b / a }, 1);
    const complex<double> sy5 = G({ 0, 0 }, { 1, 1 }, x);
    const complex<double> sy6 = G({ 0 }, { 1 }, x);
    complex<double> res {
        G({ 0, 0, 0, b / a, c / a }, 1) + G({ 0, 0, b / a, 0, c / a }, 1)
        + G({ 0, 0, b / a, c / a, x / a }, 1) + G({ 0, b / a, 0, 0, c / a }, 1)
        + G({ 0, b / a, 0, c / a, x / a }, 1) + G({ 0, b / a, c / a, 0, x / a }, 1)
        + G({ b / a, 0, 0, 0, c / a }, 1) + G({ b / a, 0, 0, c / a, x / a }, 1)
        + G({ b / a, 0, c / a, 0, x / a }, 1) + G({ b / a, c / a, 0, 0, x / a }, 1)
        + (sy4 + G({ 0, b / a, c / a }, 1) + G({ b / a, 0, c / a }, 1))
            * G({ 0, c }, { 1, sc }, x)
        + G({ 0, b / a }, 1) * (-sy3 + G({ 0, b, c }, { 1, sb, sc }, x))
        + G({ b / a }, 1)
            * (G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                - G({ 0, 0, b, c }, { 1, 1, sb, sc }, x))
        + G({ 0, 0, 0, b, c }, { 1, 1, 1, sb, sc }, x)
        + (sy2 * Log(-x, sa) * pow(sy1, 2.)) / 2. - (sy2 * pow(sy1, 3.)) / 6.
        + sy1 * sy2 * (sy5 + 2. * Zeta(2))
        + sy2 * (-sy3 - sy5 * sy6 + G({ 0, 0, 0 }, { 1, 1, 1 }, x) - 2. * sy6 * Zeta(2))
    };
    if (b != x) {
        res += (-sy4 + G({ 0, 0, x / a }, 1)) * G({ b, c }, { sb, sc }, x);
    }
    if (c != x) {
        res += (-G({ 0, 0, b / a, c / a }, 1) + G({ 0, 0, b / a, x / a }, 1)
                   - G({ 0, b / a, 0, c / a }, 1) + G({ 0, b / a, 0, x / a }, 1)
                   - G({ b / a, 0, 0, c / a }, 1) + G({ b / a, 0, 0, x / a }, 1))
            * G({ c }, { sc }, x);
    }
    return res;
}

complex<double> G5_0a0bc_c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    const complex<double> sy1 = G({ 0, b / c }, 1);
    const complex<double> sy2 = G({ 0, 0, a }, { 1, 1, sa }, x);
    const complex<double> sy3 = Log(c, sc);
    const complex<double> sy4 = G({ b / c }, 1);
    const complex<double> sy5 = G({ 0, a }, { 1, sa }, x);
    const complex<double> sy6 = G({ b / c, 0, a / c }, 1);
    const complex<double> sy7 = G({ 0, 0 }, { 1, 1 }, x);
    const complex<double> sy8 = G({ 0, a, 0, b }, { 1, sa, 1, sb }, x);
    const complex<double> sy9 = G({ 0, b / c, 0, a / c }, 1);
    const complex<double> sy10 = G({ b / c, 0, 0, a / c }, 1);
    complex<double> res { 2. * sy1 * sy2 + sy6 * sy7
        + sy3 * (2. * sy10 - 2. * sy2 * sy4 + sy1 * sy5 + sy9)
        + sy5 * (sy6 - G({ 0, 0, b / c }, 1)) - G({ 0, 0, b / c, 0, a / c }, 1)
        - 2. * G({ 0, b / c, 0, 0, a / c }, 1) - 3. * G({ b / c, 0, 0, 0, a / c }, 1)
        + G({ b / c, 0, a / c, 0, x / c }, 1)
        + sy8 * (-G({ x / c }, 1) + G({ c }, { sc }, x))
        - 2. * G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x)
        - 2. * G({ 0, a, 0, 0, b }, { 1, sa, 1, 1, sb }, x)
        + (-2. * sy10 + 2. * sy2 * sy4 - sy1 * sy5 + sy3 * (sy4 * sy5 + sy6) - sy9)
            * Log(-x, sc)
        + (-(sy4 * sy5) / 2. - sy6 / 2.) * pow(sy3, 2.) + 2. * sy6 * Zeta(2)
        + sy4
            * (sy8 - 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                + sy5 * (sy7 + 2. * Zeta(2))) };
    return res;
}
complex<double> G5_0a0bc_b(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (b == c) {
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G5_0a0bc_b: b==c but sb!=sc" };
        }

        const complex<double> sy1 = G({ 0, a, b }, { 1, sa, sb }, x);
        const complex<double> sy2 = G({ 0, a / b, 1 }, 1);
        complex<double> res { -G({ 0, a / b, 0, 0, 1 }, 1)
            + G({ 0, a / b, 0, x / b, 1 }, 1)
            + (G({ 0, a / b, 0, 1 }, 1) - G({ 0, a / b, 0, x / b }, 1))
                * G({ b }, { sb }, x)
            - sy2 * G({ 0, b }, { 1, sb }, x)
            + G({ 0, a / b }, 1) * (-sy1 + G({ 0, 0, b }, { 1, 1, sb }, x))
            + G({ 0, a, 0, 0, b }, { 1, sa, 1, 1, sb }, x) - sy1 * Zeta(2)
            + G({ 0, a }, { 1, sa }, x) * (sy2 + Zeta(3)) };
        return res;

    } else {
        const complex<double> sy1 = G({ 0, a, c }, { 1, sa, sc }, x);
        const complex<double> sy2 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy3 = G({ 0, c / b, a / b }, 1);
        const complex<double> sy4 = G({ c / b, 0, a / b }, 1);
        const complex<double> sy5 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy6 = G({ 0, a / b, c / b }, 1);
        const complex<double> sy7 = Log(b, sb);
        const complex<double> sy8 = G({ c / b }, 1);
        const complex<double> sy9 = G({ 0, 0, a }, { 1, 1, sa }, x);
        const complex<double> sy10 = G({ 0, 0, a / b, c / b }, 1);
        const complex<double> sy11 = G({ 0, 0, c / b, a / b }, 1);
        const complex<double> sy12 = G({ 0, c / b, 0, a / b }, 1);
        const complex<double> sy13 = G({ c / b, 0, 0, a / b }, 1);

        complex<double> res { sy5 * (-sy3 - sy4 - sy6)
            + sy7 * (-2. * sy10 - 2. * sy11 - 2. * sy12 - 2. * sy13 + 2. * sy8 * sy9)
            + sy1 * G({ 0, c / b }, 1) + sy2 * (-sy3 - sy4 - G({ 0, 0, c / b }, 1))
            + 3. * G({ 0, 0, 0, a / b, c / b }, 1) + 3. * G({ 0, 0, 0, c / b, a / b }, 1)
            + 3. * G({ 0, 0, c / b, 0, a / b }, 1) - G({ 0, a / b, 0, 0, c / b }, 1)
            - G({ 0, a / b, 0, c / b, x / b }, 1) - G({ 0, a / b, c / b, 0, x / b }, 1)
            + 3. * G({ 0, c / b, 0, 0, a / b }, 1) - G({ 0, c / b, a / b, 0, x / b }, 1)
            + 3. * G({ c / b, 0, 0, 0, a / b }, 1) - G({ c / b, 0, a / b, 0, x / b }, 1)
            + G({ 0, a / b }, 1) * (-sy1 + G({ 0, 0, c }, { 1, 1, sc }, x))
            + G({ 0, a, 0, 0, c }, { 1, sa, 1, 1, sc }, x)
            + (2. * sy10 + 2. * sy11 + 2. * sy12 + 2. * sy13
                  + sy7 * (-sy3 - sy4 - sy6 - sy2 * sy8) - 2. * sy8 * sy9)
                * Log(-x, sb)
            + (sy3 / 2. + sy4 / 2. + sy6 / 2. + (sy2 * sy8) / 2.) * pow(sy7, 2.)
            + sy8
                * (3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                    - G({ 0, a, 0, c }, { 1, sa, 1, sc }, x)
                    + sy2 * (-sy5 - 2. * Zeta(2)))
            + sy6 * (-G({ 0, c }, { 1, sc }, x) - 2. * Zeta(2)) - 2. * sy3 * Zeta(2)
            - 2. * sy4 * Zeta(2) };
        if (c != x) {
            res += (G({ 0, a / b, 0, c / b }, 1) - G({ 0, a / b, 0, x / b }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G5_0a0bc_a(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (a == b) {
        if (sa != sb) {
            throw FastGPL::FastGPL_error { "G5_0a0bc: a==b but sa!=sb" };
        }
        const complex<double> sy1 = Log(a, sa);
        const complex<double> sy2 = G({ 0, 1, c / a }, 1);
        const complex<double> sy3 = G({ 0, a, c }, { 1, sa, sc }, x);
        const complex<double> sy4 = G({ 0, c }, { 1, sc }, x);
        complex<double> res { -(sy3 * G({ 0, x / a }, 1))
            - 3. * G({ 0, 0, 0, 1, c / a }, 1) - 2. * G({ 0, 0, 1, 0, c / a }, 1)
            - 2. * G({ 0, 0, 1, c / a, x / a }, 1) - G({ 0, 1, 0, 0, c / a }, 1)
            - G({ 0, 1, 0, c / a, x / a }, 1) - G({ 0, 1, c / a, 0, x / a }, 1)
            + G({ 0, 0, 0, a, c }, { 1, 1, 1, sa, sc }, x) - sy1 * sy2 * Log(-x, sa)
            + (sy2 * pow(sy1, 2.)) / 2.
            + sy2 * (-sy4 - G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) + sy3 * Zeta(2)
            - G({ 0, 0, c }, { 1, 1, sc }, x) * Zeta(2)
            + G({ a, c }, { sa, sc }, x) * (-2. * G({ 0, 0, x / a }, 1) - 2. * Zeta(3))
            + 2. * sy4 * Zeta(3) };
        if (c != x) {
            res += (2. * G({ 0, 0, 1, c / a }, 1) - 2. * G({ 0, 0, 1, x / a }, 1)
                       + G({ 0, 1, 0, c / a }, 1) - G({ 0, 1, 0, x / a }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    if (a == c) {
        if (sa != sc) {
            throw FastGPL::FastGPL_error { "G5_0a0bc: a==c but sa!=sc" };
        }
        const complex<double> sy1 = G({ 0, b, a }, { 1, sb, sa }, x);
        const complex<double> sy2 = G({ 0, 0, b / a }, 1);
        const complex<double> sy3 = G({ 0, b / a, 1 }, 1);
        const complex<double> sy4 = Log(a, sa);
        complex<double> res { -(sy1 * G({ 0, x / a }, 1))
            - 3. * G({ 0, 0, 0, b / a, 1 }, 1) - 2. * G({ 0, 0, b / a, 0, 1 }, 1)
            - 2. * G({ 0, 0, b / a, 1, x / a }, 1) - G({ 0, b / a, 0, 0, 1 }, 1)
            - G({ 0, b / a, 0, 1, x / a }, 1) - G({ 0, b / a, 1, 0, x / a }, 1)
            + (2. * G({ 0, 0, b / a, 1 }, 1) - 2. * G({ 0, 0, b / a, x / a }, 1)
                  + G({ 0, b / a, 0, 1 }, 1) - G({ 0, b / a, 0, x / a }, 1))
                * G({ a }, { sa }, x)
            + (-2. * sy2 - sy3) * G({ 0, a }, { 1, sa }, x)
            + G({ 0, b / a }, 1) * (-sy1 + G({ 0, 0, a }, { 1, 1, sa }, x))
            + G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x) - sy3 * sy4 * Log(-x, sa)
            + (sy3 * pow(sy4, 2.)) / 2.
            + sy3 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
        if (b != x) {
            res += (2. * sy2 - 2. * G({ 0, 0, x / a }, 1)) * G({ b, a }, { sb, sa }, x);
        }
        return res;
    }
    // G({0,a,b,c},x)
    const complex<double> sy1 = G({ 0, b, c }, { 1, sb, sc }, x);
    const complex<double> sy2 = G({ 0, 0, b / a }, 1);
    const complex<double> sy3 = G({ 0, b / a, c / a }, 1);
    const complex<double> sy4 = Log(a, sa);
    complex<double> res { -(sy1 * G({ 0, x / a }, 1))
        - 3. * G({ 0, 0, 0, b / a, c / a }, 1) - 2. * G({ 0, 0, b / a, 0, c / a }, 1)
        - 2. * G({ 0, 0, b / a, c / a, x / a }, 1) - G({ 0, b / a, 0, 0, c / a }, 1)
        - G({ 0, b / a, 0, c / a, x / a }, 1) - G({ 0, b / a, c / a, 0, x / a }, 1)
        + (-2. * sy2 - sy3) * G({ 0, c }, { 1, sc }, x)
        + G({ 0, b / a }, 1) * (-sy1 + G({ 0, 0, c }, { 1, 1, sc }, x))
        + G({ 0, 0, 0, b, c }, { 1, 1, 1, sb, sc }, x) - sy3 * sy4 * Log(-x, sa)
        + (sy3 * pow(sy4, 2.)) / 2. + sy3 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
    if (b != x) {
        res += 2. * (sy2 - G({ 0, 0, x / a }, 1)) * G({ b, c }, { sb, sc }, x);
    }
    if (c != x) {
        res += (2. * G({ 0, 0, b / a, c / a }, 1) - 2. * G({ 0, 0, b / a, x / a }, 1)
                   + G({ 0, b / a, 0, c / a }, 1) - G({ 0, b / a, 0, x / a }, 1))
            * G({ c }, { sc }, x);
    }
    return res;
}

complex<double> G5_0ab0c_c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;

    const complex<double> sy1 = G({ 0, x / c }, 1);
    const complex<double> sy2 = G({ 0, b, a }, { 1, sb, sa }, x);
    const complex<double> sy3 = G({ 0, b / c }, 1);
    const complex<double> sy4 = Log(c, sc);
    const complex<double> sy5 = G({ 0, b / c, a / c }, 1);
    const complex<double> sy6 = G({ 0, a }, { 1, sa }, x);
    const complex<double> sy7 = G({ 0, 0, b / c, a / c }, 1);
    const complex<double> sy8 = G({ 0, b / c, 0, a / c }, 1);

    complex<double> res { -(sy1 * sy2) + sy4 * (2. * sy7 + sy8)
        - 3. * G({ 0, 0, 0, b / c, a / c }, 1) - 2. * G({ 0, 0, b / c, 0, a / c }, 1)
        - G({ 0, b / c, 0, 0, a / c }, 1) + G({ 0, b / c, a / c, 0, x / c }, 1)
        + sy3 * (sy2 - 2. * G({ 0, 0, a }, { 1, 1, sa }, x))
        + sy6
            * (-(sy3 * sy4) + sy5 + 2. * G({ 0, 0, b / c }, 1)
                + G({ 0, 0, b }, { 1, 1, sb }, x))
        + G({ 0, c }, { 1, sc }, x) * G({ 0, a, b }, { 1, sa, sb }, x)
        + G({ c }, { sc }, x)
            * (-2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
        + G({ x / c }, 1)
            * (sy6 * G({ 0, b }, { 1, sb }, x)
                - 2. * G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                - G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
        - 3. * G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x)
        - G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x)
        + (sy4 * sy5 + sy3 * sy6 - 2. * sy7 - sy8) * Log(-x, sc)
        - (sy5 * pow(sy4, 2.)) / 2. + sy5 * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2)) };
    if (b != x) {
        res += (sy1 * sy6 - sy3 * sy6) * G({ b }, { sb }, x)
            + (-sy1 + sy3) * G({ b, 0, a }, { sb, 1, sa }, x);
    }
    return res;
}
complex<double> G5_0ab0c_b(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (b == c) {
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G5_0ab0c: b==c but sb!=sc" };
        }

        const complex<double> sy1 = G({ 0, a, b }, { 1, sa, sb }, x);
        const complex<double> sy2 = Log(b, sb);
        const complex<double> sy3 = G({ 0, 1, a / b }, 1);
        const complex<double> sy4 = G({ 0, a / b, 1 }, 1);
        const complex<double> sy5 = G({ a / b, 0, 1 }, 1);
        const complex<double> sy6 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy7 = G({ 0, b }, { 1, sb }, x);
        const complex<double> sy8 = G({ b }, { sb }, x);
        const complex<double> sy9 = G({ 0, a / b, 0, 1 }, 1);
        const complex<double> sy10 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy11 = G({ 0, 0, 1, a / b }, 1);
        const complex<double> sy12 = G({ 0, 0, a / b, 1 }, 1);
        const complex<double> sy13 = G({ 0, 1, 0, a / b }, 1);

        complex<double> res { (sy4 + sy5) * sy6 - sy8 * sy9
            + sy7 * (sy5 + G({ a / b, 0, x / b }, 1))
            + sy8
                * (G({ 0, a / b, 0, x / b }, 1) - 2. * G({ a / b, 0, 0, 1 }, 1)
                    + 2. * G({ a / b, 0, 0, x / b }, 1))
            - 3. * G({ 0, 0, 0, 1, a / b }, 1) - 3. * G({ 0, 0, 0, a / b, 1 }, 1)
            - 2. * G({ 0, 0, 1, 0, a / b }, 1) - G({ 0, 0, a / b, 0, 1 }, 1)
            - G({ 0, 1, 0, 0, a / b }, 1) + G({ 0, 1, a / b, 0, x / b }, 1)
            + G({ 0, a / b, 0, 0, 1 }, 1) + G({ 0, a / b, 0, 1, x / b }, 1)
            + G({ 0, a / b, 1, 0, x / b }, 1) + 3. * G({ a / b, 0, 0, 0, 1 }, 1)
            + 2. * G({ a / b, 0, 0, 1, x / b }, 1) + G({ a / b, 0, 1, 0, x / b }, 1)
            + G({ 0, a / b }, 1) * (sy1 - G({ 0, 0, b }, { 1, 1, sb }, x))
            + G({ a / b }, 1)
                * (-G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
                    + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
            + G({ 0, a, 0, 0, b }, { 1, sa, 1, 1, sb }, x)
            + (-sy3 / 2. - sy4 / 2. - sy5 / 2.) * pow(sy2, 2.) + sy1 * Zeta(2)
            + 2. * sy5 * Zeta(2) + 2. * G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(2)
            + sy3 * (sy10 + sy6 + 2. * Zeta(2)) + sy4 * (sy7 + 2. * Zeta(2))
            + Log(-x, sb)
                * (-2. * sy11 - 2. * sy12 - sy13 + sy2 * (sy3 + sy4 + sy5) - sy9
                    - sy10 * Zeta(2))
            + sy2 * (2. * sy11 + 2. * sy12 + sy13 + sy9 + sy10 * Zeta(2))
            - 2. * sy10 * Zeta(3) };
        return res;

    } else {

        const complex<double> sy1 = G({ 0, c / b }, 1);
        const complex<double> sy2 = G({ 0, a, c }, { 1, sa, sc }, x);
        const complex<double> sy3 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy4 = G({ 0, c / b, a / b }, 1);
        const complex<double> sy5 = Log(b, sb);
        const complex<double> sy6 = G({ 0, a / b, c / b }, 1);
        const complex<double> sy7 = G({ a / b, 0, c / b }, 1);
        const complex<double> sy8 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy9 = G({ 0, c }, { 1, sc }, x);
        const complex<double> sy10 = G({ 0, 0, a / b, c / b }, 1);
        const complex<double> sy11 = G({ 0, 0, c / b, a / b }, 1);
        const complex<double> sy12 = G({ 0, a / b, 0, c / b }, 1);
        const complex<double> sy13 = G({ 0, c / b, 0, a / b }, 1);
        complex<double> res { (2. * sy10 + 2. * sy11 + sy12 + sy13 - sy1 * sy3) * sy5
            + (sy4 + sy7) * sy8 + sy3 * (sy4 + 2. * G({ 0, 0, c / b }, 1))
            + sy9 * (sy7 + G({ a / b, 0, x / b }, 1))
            - 3. * G({ 0, 0, 0, a / b, c / b }, 1) - 3. * G({ 0, 0, 0, c / b, a / b }, 1)
            - G({ 0, 0, a / b, 0, c / b }, 1) - 2. * G({ 0, 0, c / b, 0, a / b }, 1)
            + G({ 0, a / b, 0, 0, c / b }, 1) + G({ 0, a / b, 0, c / b, x / b }, 1)
            + G({ 0, a / b, c / b, 0, x / b }, 1) - G({ 0, c / b, 0, 0, a / b }, 1)
            + G({ 0, c / b, a / b, 0, x / b }, 1) + 3. * G({ a / b, 0, 0, 0, c / b }, 1)
            + 2. * G({ a / b, 0, 0, c / b, x / b }, 1)
            + G({ a / b, 0, c / b, 0, x / b }, 1)
            + sy1 * (-sy2 - 2. * G({ 0, 0, a }, { 1, 1, sa }, x))
            + G({ 0, a / b }, 1) * (sy2 - G({ 0, 0, c }, { 1, 1, sc }, x))
            + G({ a / b }, 1)
                * (-G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                    + G({ 0, a, 0, c }, { 1, sa, 1, sc }, x))
            + G({ 0, a, 0, 0, c }, { 1, sa, 1, 1, sc }, x)
            + (-2. * sy10 - 2. * sy11 - sy12 - sy13 + sy1 * sy3 + sy5 * (sy4 + sy6 + sy7))
                * Log(-x, sb)
            + (-sy4 / 2. - sy6 / 2. - sy7 / 2.) * pow(sy5, 2.) + 2. * sy4 * Zeta(2)
            + 2. * sy7 * Zeta(2) + sy6 * (sy8 + sy9 + 2. * Zeta(2)) };
        if (c != x) {
            res += (-sy12 + G({ 0, a / b, 0, x / b }, 1)
                       - 2. * G({ a / b, 0, 0, c / b }, 1)
                       + 2. * G({ a / b, 0, 0, x / b }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G5_0ab0c_a(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (a == b) {
        if (sa != sb) {
            throw FastGPL::FastGPL_error { "G5_0ab0c: a==b but sa!=sb" };
        }
        const complex<double> sy1 = Log(a, sa);
        const complex<double> sy2 = G({ 0, 1, c / a }, 1);
        const complex<double> sy3 = G({ 0, c / a, 1 }, 1);
        const complex<double> sy4 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy5 = G({ 0, c }, { 1, sc }, x);
        const complex<double> sy6 = G({ a, 0, c }, { sa, 1, sc }, x);
        complex<double> res {
            sy3 * sy4 - sy6 * G({ 0, x / a }, 1) + sy5 * (sy3 + G({ 0, x / a, 1 }, 1))
            + 3. * G({ 0, 0, 0, 1, c / a }, 1) + 3. * G({ 0, 0, 0, c / a, 1 }, 1)
            + 2. * G({ 0, 0, 1, 0, c / a }, 1) + 2. * G({ 0, 0, 1, c / a, x / a }, 1)
            + 2. * G({ 0, 0, c / a, 1, x / a }, 1) + 2. * G({ 0, 0, c / a, x / a, 1 }, 1)
            + G({ 0, 1, 0, 0, c / a }, 1) + G({ 0, 1, 0, c / a, x / a }, 1)
            + G({ 0, 1, c / a, 0, x / a }, 1) + G({ 0, c / a, 0, 1, x / a }, 1)
            + G({ 0, c / a, 0, x / a, 1 }, 1) + G({ 0, c / a, 1, 0, x / a }, 1)
            + G({ 0, 0, a, 0, c }, { 1, 1, sa, 1, sc }, x)
            + sy1 * (sy2 + sy3) * Log(-x, sa) + (-sy2 / 2. - sy3 / 2.) * pow(sy1, 2.)
            + 2. * sy3 * Zeta(2) - sy6 * Zeta(2)
            + G({ 0, 0, c }, { 1, 1, sc }, x) * Zeta(2) + sy2 * (sy4 + sy5 + 2. * Zeta(2))
        };
        if (c != x) {
            res += (-2. * G({ 0, 0, 1, c / a }, 1) + 2. * G({ 0, 0, 1, x / a }, 1)
                       - 2. * G({ 0, 0, c / a, 1 }, 1) + 2. * G({ 0, 0, x / a, 1 }, 1)
                       - G({ 0, 1, 0, c / a }, 1) + G({ 0, 1, 0, x / a }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    if (a == c) {
        if (sa != sc) {
            throw FastGPL::FastGPL_error { "G5_0ab0c: a==c but sa!=sc" };
        }

        const complex<double> sy1 = G({ 0, b / a }, 1);
        const complex<double> sy2 = Log(a, sa);
        const complex<double> sy3 = G({ b / a, 0, 1 }, 1);
        complex<double> res {
            -G({ 0, 0, b / a, 0, 1 }, 1) - 2. * G({ 0, b / a, 0, 0, 1 }, 1)
            - G({ 0, b / a, 0, 1, x / a }, 1) - 3. * G({ b / a, 0, 0, 0, 1 }, 1)
            - 2. * G({ b / a, 0, 0, 1, x / a }, 1) - G({ b / a, 0, 1, 0, x / a }, 1)
            + (G({ 0, b / a, 0, 1 }, 1) - G({ 0, b / a, 0, x / a }, 1)
                  + 2. * G({ b / a, 0, 0, 1 }, 1) - 2. * G({ b / a, 0, 0, x / a }, 1))
                * G({ a }, { sa }, x)
            + (-sy3 - G({ 0, b / a, x / a }, 1) - G({ b / a, 0, x / a }, 1))
                * G({ 0, a }, { 1, sa }, x)
            - sy1 * G({ 0, 0, a }, { 1, 1, sa }, x)
            + G({ b / a }, 1)
                * (G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                    - G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
            + G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x) - sy2 * sy3 * Log(-x, sa)
            + (sy3 * pow(sy2, 2.)) / 2. + sy3 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))
        };
        if (b != x) {
            res += (sy1 - G({ 0, x / a }, 1)) * G({ b, 0, a }, { sb, 1, sa }, x);
        }
        return res;
    }
    // G({0,a,b,c},x)
    const complex<double> sy1 = G({ 0, b / a }, 1);
    const complex<double> sy2 = Log(a, sa);
    const complex<double> sy3 = G({ b / a, 0, c / a }, 1);
    complex<double> res { -G({ 0, 0, b / a, 0, c / a }, 1)
        - 2. * G({ 0, b / a, 0, 0, c / a }, 1) - G({ 0, b / a, 0, c / a, x / a }, 1)
        - 3. * G({ b / a, 0, 0, 0, c / a }, 1) - 2. * G({ b / a, 0, 0, c / a, x / a }, 1)
        - G({ b / a, 0, c / a, 0, x / a }, 1)
        + (-sy3 - G({ 0, b / a, x / a }, 1) - G({ b / a, 0, x / a }, 1))
            * G({ 0, c }, { 1, sc }, x)
        - sy1 * G({ 0, 0, c }, { 1, 1, sc }, x)
        + G({ b / a }, 1)
            * (G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                - G({ 0, b, 0, c }, { 1, sb, 1, sc }, x))
        + G({ 0, 0, b, 0, c }, { 1, 1, sb, 1, sc }, x) - sy2 * sy3 * Log(-x, sa)
        + (sy3 * pow(sy2, 2.)) / 2. + sy3 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
    if (b != x) {
        res += (sy1 - G({ 0, x / a }, 1)) * G({ b, 0, c }, { sb, 1, sc }, x);
    }
    if (c != x) {
        res += (G({ 0, b / a, 0, c / a }, 1) - G({ 0, b / a, 0, x / a }, 1)
                   + 2. * G({ b / a, 0, 0, c / a }, 1)
                   - 2. * G({ b / a, 0, 0, x / a }, 1))
            * G({ c }, { sc }, x);
    }
    return res;
}

complex<double> G5_a00bc_c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;

    const complex<double> sy1 = Log(c, sc);
    const complex<double> sy2 = G({ a }, { sa }, x);
    const complex<double> sy3 = G({ b / c }, 1);
    const complex<double> sy4 = G({ 0, a }, { 1, sa }, x);
    const complex<double> sy5 = G({ 0, b / c }, 1);
    const complex<double> sy6 = G({ 0, 0, a }, { 1, 1, sa }, x);
    const complex<double> sy7 = G({ 0, 0, b / c }, 1);
    const complex<double> sy8 = G({ a, 0, 0, b }, { sa, 1, 1, sb }, x);
    const complex<double> sy9 = G({ b / c, 0, 0, a / c }, 1);
    const complex<double> sy10 = G({ 0, 0 }, { 1, 1 }, x);
    const complex<double> sy11 = G({ 0 }, { 1 }, x);
    complex<double> res { -(sy5 * sy6) + sy4 * sy7 + G({ 0, b / c, 0, 0, a / c }, 1)
        + 3. * G({ b / c, 0, 0, 0, a / c }, 1) + G({ b / c, 0, 0, a / c, x / c }, 1)
        + sy8 * (-G({ x / c }, 1) + G({ c }, { sc }, x))
        - G({ 0, a, 0, 0, b }, { 1, sa, 1, 1, sb }, x)
        - 3. * G({ a, 0, 0, 0, b }, { sa, 1, 1, 1, sb }, x)
        + ((sy3 * sy4) / 2. - (sy2 * sy5) / 2.) * pow(sy1, 2.)
        + Log(-x, sc)
            * (sy4 * sy5 + sy1 * (-(sy3 * sy4) + sy2 * sy5) - sy3 * sy6 - sy2 * sy7 + sy9
                - (sy2 * sy3 * pow(sy1, 2.)) / 2.)
        + (sy2 * sy3 * pow(sy1, 3.)) / 6.
        + sy1
            * (-(sy4 * sy5) + sy3 * sy6 - sy9
                + sy2 * (sy7 + sy3 * (-sy10 - 2. * Zeta(2))))
        + sy3
            * (-(sy10 * sy4) + sy8 + G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                - 2. * sy4 * Zeta(2))
        + sy2
            * (sy10 * sy5 - sy9 - G({ 0, 0, 0, b / c }, 1) + 2. * sy5 * Zeta(2)
                + sy3
                    * (sy10 * sy11 - G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                        + 2. * sy11 * Zeta(2))) };
    return res;
}
complex<double> G5_a00bc_b(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (b == c) {
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G5_a00bc: b==c but sb!=sc" };
        }

        const complex<double> sy1 = G({ 0, 0, a / b }, 1);
        const complex<double> sy2 = G({ a }, { sa }, x);
        const complex<double> sy3 = G({ 0, 0, a / b, 1 }, 1);
        const complex<double> sy4 = G({ b }, { sb }, x);
        complex<double> res { -(sy2 * sy3) + sy3 * sy4
            - sy4 * G({ 0, 0, a / b, x / b }, 1) - G({ 0, 0, a / b, 0, 1 }, 1)
            + G({ 0, 0, a / b, x / b, 1 }, 1) - sy1 * G({ 0, b }, { 1, sb }, x)
            + G({ a, 0, 0, 0, b }, { sa, 1, 1, 1, sb }, x)
            - G({ a, 0, b }, { sa, 1, sb }, x) * Zeta(2) - sy2 * Zeta(4)
            + G({ a, b }, { sa, sb }, x) * (sy1 + Zeta(3)) };
        return res;

    } else {
        const complex<double> sy1 = Log(b, sb);
        const complex<double> sy2 = G({ a }, { sa }, x);
        const complex<double> sy3 = G({ c / b }, 1);
        const complex<double> sy4 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy5 = G({ 0, 0, a / b }, 1);
        const complex<double> sy6 = G({ a, c }, { sa, sc }, x);
        const complex<double> sy7 = G({ 0, 0, a }, { 1, 1, sa }, x);
        const complex<double> sy8 = G({ 0, 0, a / b, c / b }, 1);
        const complex<double> sy9 = G({ 0, 0, c / b, a / b }, 1);
        const complex<double> sy10 = G({ 0, c / b, 0, a / b }, 1);
        const complex<double> sy11 = G({ c / b, 0, 0, a / b }, 1);
        const complex<double> sy12 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy13 = G({ 0 }, { 1 }, x);

        complex<double> res { sy5 * sy6 - sy6 * G({ 0, 0, c / b }, 1)
            - 3. * G({ 0, 0, 0, a / b, c / b }, 1) - 3. * G({ 0, 0, 0, c / b, a / b }, 1)
            - G({ 0, 0, a / b, 0, c / b }, 1) - G({ 0, 0, a / b, c / b, x / b }, 1)
            - 3. * G({ 0, 0, c / b, 0, a / b }, 1) - G({ 0, 0, c / b, a / b, x / b }, 1)
            - 3. * G({ 0, c / b, 0, 0, a / b }, 1) - G({ 0, c / b, 0, a / b, x / b }, 1)
            - 3. * G({ c / b, 0, 0, 0, a / b }, 1) - G({ c / b, 0, 0, a / b, x / b }, 1)
            - sy5 * G({ 0, c }, { 1, sc }, x)
            + G({ 0, c / b }, 1) * G({ a, 0, c }, { sa, 1, sc }, x)
            + G({ a, 0, 0, 0, c }, { sa, 1, 1, 1, sc }, x)
            - (sy3 * sy4 * pow(sy1, 2.)) / 2.
            + Log(-x, sb)
                * (-sy10 - sy11 + sy1 * sy3 * sy4 + sy3 * sy7 - sy8 - sy9
                    + (sy2 * sy3 * pow(sy1, 2.)) / 2.)
            - (sy2 * sy3 * pow(sy1, 3.)) / 6.
            + sy3
                * (sy12 * sy4 - G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                    - G({ a, 0, 0, c }, { sa, 1, 1, sc }, x) + 2. * sy4 * Zeta(2))
            + sy1
                * (sy10 + sy11 - sy3 * sy7 + sy8 + sy9
                    + sy2 * sy3 * (sy12 + 2. * Zeta(2)))
            + sy2
                * (sy10 + sy11 + sy9 + G({ 0, 0, 0, c / b }, 1)
                    + sy3
                        * (-(sy12 * sy13) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                            - 2. * sy13 * Zeta(2))) };
        if (c != x) {
            res += (sy8 - G({ 0, 0, a / b, x / b }, 1)) * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G5_a00bc_a(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (a == b) {
        if (sa != sb) {
            throw FastGPL::FastGPL_error { "G5_a00bc: a==b but sa!=sb" };
        }
        const complex<double> sy1 = G({ 0, 0, 1, c / a }, 1);
        complex<double> res { 3. * G({ 0, 0, 0, 1, c / a }, 1)
            + G({ 0, 0, 1, 0, c / a }, 1) + G({ 0, 0, 1, c / a, x / a }, 1)
            + G({ 0, x / a }, 1) * G({ 0, a, c }, { 1, sa, sc }, x)
            + G({ x / a }, 1) * G({ 0, 0, a, c }, { 1, 1, sa, sc }, x)
            + G({ 0, 0, 0, a, c }, { 1, 1, 1, sa, sc }, x) - sy1 * Log(a, sa)
            + sy1 * Log(-x, sa) - G({ 0, c }, { 1, sc }, x) * Zeta(3)
            + G({ a, c }, { sa, sc }, x) * (G({ 0, 0, x / a }, 1) + Zeta(3)) };
        if (c != x) {
            res += (-sy1 + G({ 0, 0, 1, x / a }, 1)) * G({ c }, { sc }, x);
        }
        return res;
    }

    if (a == c) {
        if (sa != sc) {
            throw FastGPL::FastGPL_error { "G5_a00bc: a==c but sa!=sc" };
        }

        const complex<double> sy1 = G({ 0, 0, b / a }, 1);
        const complex<double> sy2 = G({ 0, 0, b / a, 1 }, 1);
        const complex<double> sy3 = G({ a }, { sa }, x);
        complex<double> res { -(sy2 * sy3) + sy3 * G({ 0, 0, b / a, x / a }, 1)
            + 3. * G({ 0, 0, 0, b / a, 1 }, 1) + G({ 0, 0, b / a, 0, 1 }, 1)
            + G({ 0, 0, b / a, 1, x / a }, 1) + sy1 * G({ 0, a }, { 1, sa }, x)
            + G({ 0, x / a }, 1) * G({ 0, b, a }, { 1, sb, sa }, x)
            + G({ x / a }, 1) * G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
            + G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x) - sy2 * Log(a, sa)
            + sy2 * Log(-x, sa) };
        if (b != x) {
            res += (-sy1 + G({ 0, 0, x / a }, 1)) * G({ b, a }, { sb, sa }, x);
        }
        return res;
    }
    // G({0,a,b,c},x)

    const complex<double> sy1 = G({ 0, 0, b / a }, 1);
    const complex<double> sy2 = G({ 0, 0, b / a, c / a }, 1);
    complex<double> res { 3. * G({ 0, 0, 0, b / a, c / a }, 1)
        + G({ 0, 0, b / a, 0, c / a }, 1) + G({ 0, 0, b / a, c / a, x / a }, 1)
        + sy1 * G({ 0, c }, { 1, sc }, x)
        + G({ 0, x / a }, 1) * G({ 0, b, c }, { 1, sb, sc }, x)
        + G({ x / a }, 1) * G({ 0, 0, b, c }, { 1, 1, sb, sc }, x)
        + G({ 0, 0, 0, b, c }, { 1, 1, 1, sb, sc }, x) - sy2 * Log(a, sa)
        + sy2 * Log(-x, sa) };
    if (b != x) {
        res += (-sy1 + G({ 0, 0, x / a }, 1)) * G({ b, c }, { sb, sc }, x);
    }
    if (c != x) {
        res += (-sy2 + G({ 0, 0, b / a, x / a }, 1)) * G({ c }, { sc }, x);
    }
    return res;
}

complex<double> G5_a0b0c_c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;

    const complex<double> sy1 = Log(c, sc);
    const complex<double> sy2 = G({ a }, { sa }, x);
    const complex<double> sy3 = G({ 0, b / c }, 1);
    const complex<double> sy4 = G({ 0, b }, { 1, sb }, x);
    const complex<double> sy5 = G({ 0, c }, { 1, sc }, x);
    const complex<double> sy6 = G({ 0, x / c }, 1);
    const complex<double> sy7 = G({ 0, a }, { 1, sa }, x);
    const complex<double> sy8 = G({ 0, 0, b }, { 1, 1, sb }, x);
    const complex<double> sy9 = G({ 0, 0, b / c }, 1);
    const complex<double> sy10 = G({ c }, { sc }, x);
    const complex<double> sy11 = G({ 0, b / c, 0, a / c }, 1);

    complex<double> res { sy2 * sy4 * (-sy3 + sy5 + sy6)
        + sy7 * (sy1 * sy3 - sy8 - 2. * sy9) + sy1 * (-sy11 - 2. * sy2 * sy9)
        + 2. * G({ 0, 0, b / c, 0, a / c }, 1) + 2. * G({ 0, b / c, 0, 0, a / c }, 1)
        + G({ 0, b / c, 0, a / c, x / c }, 1) + sy3 * G({ 0, 0, a }, { 1, 1, sa }, x)
        + sy5 * (-G({ 0, a, b }, { 1, sa, sb }, x) - G({ 0, b, a }, { 1, sb, sa }, x))
        + sy10
            * (2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                + 2. * G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
        + G({ x / c }, 1)
            * (-(sy4 * sy7) + 2. * sy2 * sy8 + G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
        + G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x)
        + (sy11 - sy1 * sy2 * sy3 - sy3 * sy7 + 2. * sy2 * sy9) * Log(-x, sc)
        + (sy2 * sy3 * pow(sy1, 2.)) / 2.
        + sy2
            * (-sy11 - 2. * sy10 * sy8 + 3. * G({ 0, 0, 0, b / c }, 1)
                + 3. * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
                + sy3 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
    if (b != x) {
        res += (sy3 * sy7 - sy6 * sy7) * G({ b }, { sb }, x)
            + (-sy3 + sy6) * G({ b, 0, a }, { sb, 1, sa }, x);
    }
    return res;
}
complex<double> G5_a0b0c_b(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (b == c) {
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G5_a0b0c: b==c but sb!=sc" };
        }

        const complex<double> sy1 = G({ 0, b }, { 1, sb }, x);
        const complex<double> sy2 = G({ 0, 0, a / b }, 1);
        const complex<double> sy3 = G({ a, b }, { sa, sb }, x);
        const complex<double> sy4 = G({ 0, a / b }, 1);
        const complex<double> sy5 = G({ a, 0, b }, { sa, 1, sb }, x);
        const complex<double> sy6 = G({ b }, { sb }, x);
        const complex<double> sy7 = G({ 0, 0, a / b, 1 }, 1);
        const complex<double> sy8 = G({ 0, a / b, 0, 1 }, 1);
        const complex<double> sy9 = G({ 0, 0, 1, a / b }, 1);
        const complex<double> sy10 = G({ 0, 1, 0, a / b }, 1);
        const complex<double> sy11 = Log(b, sb);
        const complex<double> sy12 = G({ a }, { sa }, x);
        const complex<double> sy13 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy14 = Log(-x, sb);

        complex<double> res { sy10 * sy14 + 2. * sy1 * sy2 - 2. * sy2 * sy3 - sy4 * sy5
            + 2. * sy14 * sy7 - 2. * sy6 * sy7 + sy14 * sy8 - sy6 * sy8 + 2. * sy14 * sy9
            + sy1 * G({ 0, a / b, x / b }, 1) + 2. * sy6 * G({ 0, 0, a / b, x / b }, 1)
            + sy6 * G({ 0, a / b, 0, x / b }, 1) + 6. * G({ 0, 0, 0, 1, a / b }, 1)
            + 6. * G({ 0, 0, 0, a / b, 1 }, 1) + 4. * G({ 0, 0, 1, 0, a / b }, 1)
            + 2. * G({ 0, 0, 1, a / b, x / b }, 1) + 4. * G({ 0, 0, a / b, 0, 1 }, 1)
            + 2. * G({ 0, 0, a / b, 1, x / b }, 1) + 2. * G({ 0, 1, 0, 0, a / b }, 1)
            + G({ 0, 1, 0, a / b, x / b }, 1) + 2. * G({ 0, a / b, 0, 0, 1 }, 1)
            + G({ 0, a / b, 0, 1, x / b }, 1) + sy4 * G({ 0, 0, b }, { 1, 1, sb }, x)
            + G({ a, 0, 0, 0, b }, { sa, 1, 1, 1, sb }, x) + sy13 * sy14 * Zeta(2)
            + sy5 * Zeta(2) - G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(2)
            - (sy12 * pow(sy11, 2.) * Zeta(2)) / 2.
            - sy11 * (sy10 + 2. * sy7 + sy8 + 2. * sy9 - (-sy13 + sy12 * sy14) * Zeta(2))
            - sy12
                * (sy10 + 2. * sy9 - 2. * pow(Zeta(2), 2.)
                    - G({ 0, 0 }, { 1, 1 }, x) * Zeta(2) - 3. * Zeta(4))
            - 2. * sy3 * Zeta(3) };
        return res;

    } else {

        const complex<double> sy1 = Log(b, sb);
        const complex<double> sy2 = G({ a }, { sa }, x);
        const complex<double> sy3 = G({ 0, c / b }, 1);
        const complex<double> sy4 = G({ a, c }, { sa, sc }, x);
        const complex<double> sy5 = G({ 0, 0, a / b }, 1);
        const complex<double> sy6 = G({ a, 0, c }, { sa, 1, sc }, x);
        const complex<double> sy7 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy8 = G({ 0, 0, a / b, c / b }, 1);
        const complex<double> sy9 = G({ 0, 0, c / b, a / b }, 1);
        const complex<double> sy10 = G({ 0, a / b, 0, c / b }, 1);
        const complex<double> sy11 = G({ 0, c / b, 0, a / b }, 1);

        complex<double> res { -2. * sy4 * sy5
            + sy1 * (-sy10 - sy11 + sy3 * sy7 - 2. * sy8 - 2. * sy9)
            + 2. * sy4 * G({ 0, 0, c / b }, 1) + 6. * G({ 0, 0, 0, a / b, c / b }, 1)
            + 6. * G({ 0, 0, 0, c / b, a / b }, 1) + 4. * G({ 0, 0, a / b, 0, c / b }, 1)
            + 2. * G({ 0, 0, a / b, c / b, x / b }, 1)
            + 4. * G({ 0, 0, c / b, 0, a / b }, 1)
            + 2. * G({ 0, 0, c / b, a / b, x / b }, 1)
            + 2. * G({ 0, a / b, 0, 0, c / b }, 1) + G({ 0, a / b, 0, c / b, x / b }, 1)
            + 2. * G({ 0, c / b, 0, 0, a / b }, 1) + G({ 0, c / b, 0, a / b, x / b }, 1)
            + (2. * sy5 + G({ 0, a / b, x / b }, 1)) * G({ 0, c }, { 1, sc }, x)
            + sy3 * (-sy6 + G({ 0, 0, a }, { 1, 1, sa }, x))
            + G({ 0, a / b }, 1) * (-sy6 + G({ 0, 0, c }, { 1, 1, sc }, x))
            + G({ a, 0, 0, 0, c }, { sa, 1, 1, 1, sc }, x)
            + (sy10 + sy11 - sy1 * sy2 * sy3 - sy3 * sy7 + 2. * sy8 + 2. * sy9)
                * Log(-x, sb)
            + (sy2 * sy3 * pow(sy1, 2.)) / 2.
            + sy2
                * (-sy11 - 2. * sy9 - 3. * G({ 0, 0, 0, c / b }, 1)
                    + sy3 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
        if (c != x) {
            res += (-sy10 - 2. * sy8 + 2. * G({ 0, 0, a / b, x / b }, 1)
                       + G({ 0, a / b, 0, x / b }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G5_a0b0c_a(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (a == b) {
        if (sa != sb) {
            throw FastGPL::FastGPL_error { "G5_a0b0c: a==b but sa!=sb" };
        }
        const complex<double> sy1 = G({ a, 0, c }, { sa, 1, sc }, x);
        const complex<double> sy2 = G({ 0, 1, 0, c / a }, 1);
        complex<double> res { sy1 * G({ 0, x / a }, 1) + 2. * G({ 0, 0, 1, 0, c / a }, 1)
            + 2. * G({ 0, 1, 0, 0, c / a }, 1) + G({ 0, 1, 0, c / a, x / a }, 1)
            + G({ 0, 1, x / a }, 1) * G({ 0, c }, { 1, sc }, x)
            + G({ x / a }, 1) * G({ 0, a, 0, c }, { 1, sa, 1, sc }, x)
            + G({ 0, 0, a, 0, c }, { 1, 1, sa, 1, sc }, x) - sy2 * Log(a, sa)
            + sy2 * Log(-x, sa) + sy1 * Zeta(2)
            - G({ 0, 0, c }, { 1, 1, sc }, x) * Zeta(2) };
        if (c != x) {
            res += (-sy2 + G({ 0, 1, 0, x / a }, 1)) * G({ c }, { sc }, x);
        }
        return res;
    }

    if (a == c) {
        if (sa != sc) {
            throw FastGPL::FastGPL_error { "G5_a0b0c: a==c but sa!=sc" };
        }

        const complex<double> sy1 = G({ 0, b / a }, 1);
        const complex<double> sy2 = G({ 0, b / a, 0, 1 }, 1);
        const complex<double> sy3 = G({ a }, { sa }, x);
        complex<double> res { -(sy2 * sy3) + sy3 * G({ 0, b / a, 0, x / a }, 1)
            + 2. * G({ 0, 0, b / a, 0, 1 }, 1) + 2. * G({ 0, b / a, 0, 0, 1 }, 1)
            + G({ 0, b / a, 0, 1, x / a }, 1)
            + G({ 0, b / a, x / a }, 1) * G({ 0, a }, { 1, sa }, x)
            + sy1 * G({ 0, 0, a }, { 1, 1, sa }, x)
            + G({ x / a }, 1) * G({ 0, b, 0, a }, { 1, sb, 1, sa }, x)
            + G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x) - sy2 * Log(a, sa)
            + sy2 * Log(-x, sa) };
        if (b != x) {
            res += (-sy1 + G({ 0, x / a }, 1)) * G({ b, 0, a }, { sb, 1, sa }, x);
        }
        return res;
    }
    // G({0,a,b,c},x)

    const complex<double> sy1 = G({ 0, b / a }, 1);
    const complex<double> sy2 = G({ 0, b / a, 0, c / a }, 1);
    complex<double> res { 2. * G({ 0, 0, b / a, 0, c / a }, 1)
        + 2. * G({ 0, b / a, 0, 0, c / a }, 1) + G({ 0, b / a, 0, c / a, x / a }, 1)
        + G({ 0, b / a, x / a }, 1) * G({ 0, c }, { 1, sc }, x)
        + sy1 * G({ 0, 0, c }, { 1, 1, sc }, x)
        + G({ x / a }, 1) * G({ 0, b, 0, c }, { 1, sb, 1, sc }, x)
        + G({ 0, 0, b, 0, c }, { 1, 1, sb, 1, sc }, x) - sy2 * Log(a, sa)
        + sy2 * Log(-x, sa) };
    if (b != x) {
        res += (-sy1 + G({ 0, x / a }, 1)) * G({ b, 0, c }, { sb, 1, sc }, x);
    }
    if (c != x) {
        res += (-sy2 + G({ 0, b / a, 0, x / a }, 1)) * G({ c }, { sc }, x);
    }
    return res;
}

complex<double> G5_ab00c_c(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x)
{

    const complex<double> sy1 = G({ 0, 0, b / c }, 1);
    const complex<double> sy2 = G({ 0, c }, { 1, sc }, x);
    const complex<double> sy3 = G({ 0, b, a }, { 1, sb, sa }, x);
    const complex<double> sy4 = G({ 0, x / c }, 1);
    const complex<double> sy5 = G({ c }, { sc }, x);
    const complex<double> sy6 = G({ 0, 0, b, a }, { 1, 1, sb, sa }, x);
    const complex<double> sy7 = G({ x / c }, 1);
    const complex<double> sy8 = G({ a }, { sa }, x);
    const complex<double> sy9 = G({ 0, 0, b }, { 1, 1, sb }, x);
    const complex<double> sy10 = Log(c, sc);
    const complex<double> sy11 = G({ 0, 0, b / c, a / c }, 1);
    const complex<double> sy12 = G({ 0, 0, x / c }, 1);
    complex<double> res { -(sy10 * sy11) + sy2 * sy3 + sy3 * sy4 - sy5 * sy6 + sy6 * sy7
        + 3. * G({ 0, 0, 0, b / c, a / c }, 1) + G({ 0, 0, b / c, 0, a / c }, 1)
        + G({ 0, 0, b / c, a / c, x / c }, 1) + sy1 * G({ 0, a }, { 1, sa }, x)
        + G({ a, b }, { sa, sb }, x) * G({ 0, 0, c }, { 1, 1, sc }, x)
        + sy8
            * (sy1 * sy10 - sy11 + sy5 * sy9 - sy7 * sy9 - 3. * G({ 0, 0, 0, b / c }, 1)
                + (-sy2 - sy4) * G({ 0, b }, { 1, sb }, x)
                - G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x))
        + G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x)
        + (sy11 - sy1 * sy8) * Log(-x, sc) };
    if (b != x) {
        res += (sy1 * sy8 - sy12 * sy8) * G({ b }, { sb }, x)
            + (-sy1 + sy12) * G({ b, a }, { sb, sa }, x);
    }
    return res;
}
complex<double> G5_ab00c_b(complex<double> a, complex<double> b, complex<double> c,
    int sa, int sb, int sc, double x)
{
    if (b == c) {
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G5_ab00c: b==c but sb!=sc" };
        }

        const complex<double> sy1 = G({ 0, b }, { 1, sb }, x);
        const complex<double> sy2 = G({ 0, 0, a / b }, 1);
        const complex<double> sy3 = G({ a, b }, { sa, sb }, x);
        const complex<double> sy4 = G({ 0, a / b }, 1);
        const complex<double> sy5 = G({ 0, 0, b }, { 1, 1, sb }, x);
        const complex<double> sy6 = G({ a / b }, 1);
        const complex<double> sy7 = G({ a }, { sa }, x);
        const complex<double> sy8 = G({ 0, 0, 1, a / b }, 1);
        const complex<double> sy9 = G({ b }, { sb }, x);
        const complex<double> sy10 = G({ 0, 0, a / b, 1 }, 1);
        const complex<double> sy11 = G({ 0, a / b, 0, 1 }, 1);
        const complex<double> sy12 = G({ a / b, 0, 0, 1 }, 1);
        const complex<double> sy13 = Log(-x, sb);
        complex<double> res { -(sy10 * sy13) - sy11 * sy13 - sy12 * sy13 - sy1 * sy2
            + sy2 * sy3 - sy4 * sy5 - sy13 * sy8 + sy7 * sy8 + sy10 * sy9 + sy11 * sy9
            + sy12 * sy9 - sy5 * G({ a / b, x / b }, 1) - sy1 * G({ 0, a / b, x / b }, 1)
            - sy1 * G({ a / b, 0, x / b }, 1) - sy9 * G({ 0, 0, a / b, x / b }, 1)
            - sy9 * G({ 0, a / b, 0, x / b }, 1) - sy9 * G({ a / b, 0, 0, x / b }, 1)
            - 3. * G({ 0, 0, 0, 1, a / b }, 1) - 3. * G({ 0, 0, 0, a / b, 1 }, 1)
            - G({ 0, 0, 1, 0, a / b }, 1) - G({ 0, 0, 1, a / b, x / b }, 1)
            - 3. * G({ 0, 0, a / b, 0, 1 }, 1) - G({ 0, 0, a / b, 1, x / b }, 1)
            - 3. * G({ 0, a / b, 0, 0, 1 }, 1) - G({ 0, a / b, 0, 1, x / b }, 1)
            - 3. * G({ a / b, 0, 0, 0, 1 }, 1) - G({ a / b, 0, 0, 1, x / b }, 1)
            + sy4 * G({ a, 0, b }, { sa, 1, sb }, x)
            - sy6 * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
            + sy6 * G({ a, 0, 0, b }, { sa, 1, 1, sb }, x)
            + G({ a, 0, 0, 0, b }, { sa, 1, 1, 1, sb }, x) - 3. * sy7 * Zeta(4)
            + sy3 * Zeta(3) - sy13 * sy7 * Zeta(3) + G({ 0, a }, { 1, sa }, x) * Zeta(3)
            + Log(b, sb) * (sy10 + sy11 + sy12 + sy8 + sy7 * Zeta(3)) };
        return res;

    } else {

        const complex<double> sy1 = G({ a, c }, { sa, sc }, x);
        const complex<double> sy2 = G({ 0, 0, a / b }, 1);
        const complex<double> sy3 = G({ 0, 0, c }, { 1, 1, sc }, x);
        const complex<double> sy4 = G({ 0, 0, c / b }, 1);
        const complex<double> sy5 = G({ a }, { sa }, x);
        const complex<double> sy6 = G({ 0, 0, c / b, a / b }, 1);
        const complex<double> sy7 = G({ 0, 0, a / b, c / b }, 1);
        const complex<double> sy8 = G({ 0, a / b, 0, c / b }, 1);
        const complex<double> sy9 = G({ a / b, 0, 0, c / b }, 1);

        complex<double> res { sy1 * sy2 - sy1 * sy4 - sy3 * G({ a / b, x / b }, 1)
            + sy5 * (sy6 + 3. * G({ 0, 0, 0, c / b }, 1))
            - 3. * G({ 0, 0, 0, a / b, c / b }, 1) - 3. * G({ 0, 0, 0, c / b, a / b }, 1)
            - 3. * G({ 0, 0, a / b, 0, c / b }, 1) - G({ 0, 0, a / b, c / b, x / b }, 1)
            - G({ 0, 0, c / b, 0, a / b }, 1) - G({ 0, 0, c / b, a / b, x / b }, 1)
            - 3. * G({ 0, a / b, 0, 0, c / b }, 1) - G({ 0, a / b, 0, c / b, x / b }, 1)
            - 3. * G({ a / b, 0, 0, 0, c / b }, 1) - G({ a / b, 0, 0, c / b, x / b }, 1)
            - sy4 * G({ 0, a }, { 1, sa }, x)
            + (-sy2 - G({ 0, a / b, x / b }, 1) - G({ a / b, 0, x / b }, 1))
                * G({ 0, c }, { 1, sc }, x)
            + G({ 0, a / b }, 1) * (-sy3 + G({ a, 0, c }, { sa, 1, sc }, x))
            + G({ a / b }, 1)
                * (-G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                    + G({ a, 0, 0, c }, { sa, 1, 1, sc }, x))
            + G({ a, 0, 0, 0, c }, { sa, 1, 1, 1, sc }, x)
            + (-(sy4 * sy5) + sy6 + sy7 + sy8 + sy9) * Log(b, sb)
            + (sy4 * sy5 - sy6 - sy7 - sy8 - sy9) * Log(-x, sb) };
        if (c != x) {
            res += (sy7 + sy8 + sy9 - G({ 0, 0, a / b, x / b }, 1)
                       - G({ 0, a / b, 0, x / b }, 1) - G({ a / b, 0, 0, x / b }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G5_ab00c_a(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (a == b) {
        if (sa != sb) {
            throw FastGPL::FastGPL_error { "G5_ab00c: a==b but sa!=sb" };
        }
        const complex<double> sy1 = G({ 0, 0, 1, c / a }, 1);
        const complex<double> sy2 = G({ 0, 0, c / a, 1 }, 1);
        const complex<double> sy3 = G({ 0, 1, 0, c / a }, 1);
        complex<double> res { -3. * G({ 0, 0, 0, 1, c / a }, 1)
            - 3. * G({ 0, 0, 0, c / a, 1 }, 1) - 3. * G({ 0, 0, 1, 0, c / a }, 1)
            - G({ 0, 0, 1, c / a, x / a }, 1) - G({ 0, 0, c / a, 1, x / a }, 1)
            - G({ 0, 0, c / a, x / a, 1 }, 1) - 2. * G({ 0, 1, 0, 0, c / a }, 1)
            - G({ 0, 1, 0, c / a, x / a }, 1)
            + (-G({ 0, 1, x / a }, 1) - G({ 0, x / a, 1 }, 1)) * G({ 0, c }, { 1, sc }, x)
            - G({ x / a, 1 }, 1) * G({ 0, 0, c }, { 1, 1, sc }, x)
            + G({ x / a }, 1) * G({ a, 0, 0, c }, { sa, 1, 1, sc }, x)
            + G({ 0, a, 0, 0, c }, { 1, sa, 1, 1, sc }, x)
            + (sy1 + sy2 + sy3) * Log(a, sa) + (-sy1 - sy2 - sy3) * Log(-x, sa) };
        if (c != x) {
            res += (sy1 + sy2 + sy3 - G({ 0, 0, 1, x / a }, 1) - G({ 0, 0, x / a, 1 }, 1)
                       - G({ 0, 1, 0, x / a }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    if (a == c) {
        if (sa != sc) {
            throw FastGPL::FastGPL_error { "G5_ab00c: a==c but sa!=sc" };
        }

        const complex<double> sy1 = G({ b / a }, 1);
        const complex<double> sy2 = G({ b / a, 0, 0, 1 }, 1);
        const complex<double> sy3 = G({ a }, { sa }, x);
        complex<double> res { -(sy2 * sy3) + sy3 * G({ b / a, 0, 0, x / a }, 1)
            + G({ 0, b / a, 0, 0, 1 }, 1) + 3. * G({ b / a, 0, 0, 0, 1 }, 1)
            + G({ b / a, 0, 0, 1, x / a }, 1)
            + G({ b / a, 0, x / a }, 1) * G({ 0, a }, { 1, sa }, x)
            + G({ b / a, x / a }, 1) * G({ 0, 0, a }, { 1, 1, sa }, x)
            + sy1 * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
            + G({ 0, b, 0, 0, a }, { 1, sb, 1, 1, sa }, x) - sy2 * Log(a, sa)
            + sy2 * Log(-x, sa) };
        if (b != x) {
            res += (-sy1 + G({ x / a }, 1)) * G({ b, 0, 0, a }, { sb, 1, 1, sa }, x);
        }
        return res;
    }
    // G({0,a,b,c},x)

    const complex<double> sy1 = G({ b / a }, 1);
    const complex<double> sy2 = G({ b / a, 0, 0, c / a }, 1);
    complex<double> res { G({ 0, b / a, 0, 0, c / a }, 1)
        + 3. * G({ b / a, 0, 0, 0, c / a }, 1) + G({ b / a, 0, 0, c / a, x / a }, 1)
        + G({ b / a, 0, x / a }, 1) * G({ 0, c }, { 1, sc }, x)
        + G({ b / a, x / a }, 1) * G({ 0, 0, c }, { 1, 1, sc }, x)
        + sy1 * G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
        + G({ 0, b, 0, 0, c }, { 1, sb, 1, 1, sc }, x) - sy2 * Log(a, sa)
        + sy2 * Log(-x, sa) };
    if (b != x) {
        res += (-sy1 + G({ x / a }, 1)) * G({ b, 0, 0, c }, { sb, 1, 1, sc }, x);
    }
    if (c != x) {
        res += (-sy2 + G({ b / a, 0, 0, x / a }, 1)) * G({ c }, { sc }, x);
    }
    return res;
}

complex<double> G5_0abcd_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;

    const complex<double> sy1 = G({ c / d, b / d }, 1);
    const complex<double> sy2 = G({ 0, a, b }, { 1, sa, sb }, x);
    const complex<double> sy3 = Log(d, sd);
    const complex<double> sy4 = G({ c / d, b / d, a / d }, 1);
    const complex<double> sy5 = G({ 0, a }, { 1, sa }, x);
    const complex<double> sy6 = G({ 0, a, b, c }, { 1, sa, sb, sc }, x);
    const complex<double> sy7 = G({ c / d }, 1);
    const complex<double> sy8 = G({ 0, c / d, b / d, a / d }, 1);
    const complex<double> sy9 = G({ c / d, 0, b / d, a / d }, 1);
    const complex<double> sy10 = G({ c / d, b / d, 0, a / d }, 1);
    complex<double> res { sy3 * (sy10 - sy1 * sy5 + sy2 * sy7 + sy8 + sy9)
        - sy2 * G({ 0, c / d }, 1)
        + sy5 * (sy4 + G({ 0, c / d, b / d }, 1) + G({ c / d, 0, b / d }, 1))
        - G({ 0, 0, c / d, b / d, a / d }, 1) - G({ 0, c / d, 0, b / d, a / d }, 1)
        - G({ 0, c / d, b / d, 0, a / d }, 1) - G({ c / d, 0, 0, b / d, a / d }, 1)
        - G({ c / d, 0, b / d, 0, a / d }, 1) - G({ c / d, b / d, 0, 0, a / d }, 1)
        + G({ c / d, b / d, a / d, 0, x / d }, 1)
        + sy6 * (-G({ x / d }, 1) + G({ d }, { sd }, x))
        + sy1 * (-sy2 - 2. * G({ 0, 0, a }, { 1, 1, sa }, x))
        + sy7
            * (sy6 + 2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
        - 2. * G({ 0, 0, a, b, c }, { 1, 1, sa, sb, sc }, x)
        - G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x)
        - G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x)
        + (-sy10 + sy3 * sy4 + sy1 * sy5 - sy2 * sy7 - sy8 - sy9) * Log(-x, sd)
        - (sy4 * pow(sy3, 2.)) / 2. + sy4 * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2)) };

    return res;
}
complex<double> G5_0abcd_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;

    if (c == d) {
        const complex<double> sy1 = G({ 0, a, c }, { 1, sa, sc }, x);
        const complex<double> sy2 = G({ 0, a, b }, { 1, sa, sb }, x);
        const complex<double> sy3 = G({ b / c, a / c, 1 }, 1);
        complex<double> tmp { (sy1 - sy2) * G({ b / c, 1 }, 1)
            - G({ b / c, a / c, 0, 0, 1 }, 1) + G({ b / c, a / c, 0, x / c, 1 }, 1)
            + (G({ b / c, a / c, 0, 1 }, 1) - G({ b / c, a / c, 0, x / c }, 1))
                * G({ c }, { sc }, x)
            + (sy3 - G({ b / c, 0, 1 }, 1)) * G({ 0, a }, { 1, sa }, x)
            - sy3 * G({ 0, c }, { 1, sc }, x)
            + G({ b / c, a / c }, 1) * (-sy1 + G({ 0, 0, c }, { 1, 1, sc }, x))
            + G({ b / c }, 1)
                * (-G({ 0, a, 0, c }, { 1, sa, 1, sc }, x)
                    + G({ 0, a, b, c }, { 1, sa, sb, sc }, x))
            + G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x) - sy2 * Zeta(2) };
        return tmp;
    } else {

        const complex<double> sy1 = G({ 0, a, b }, { 1, sa, sb }, x);
        const complex<double> sy2 = G({ d / c, b / c }, 1);
        const complex<double> sy3 = G({ 0, 0, a }, { 1, 1, sa }, x);
        const complex<double> sy4 = G({ 0, a, d }, { 1, sa, sd }, x);
        const complex<double> sy5 = G({ b / c, d / c }, 1);
        const complex<double> sy6 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy7 = G({ b / c, d / c, a / c }, 1);
        const complex<double> sy8 = G({ d / c, b / c, a / c }, 1);
        const complex<double> sy9 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy10 = Log(c, sc);
        const complex<double> sy11 = G({ b / c, a / c, d / c }, 1);
        const complex<double> sy12 = G({ d / c }, 1);
        const complex<double> sy13 = G({ 0, a, b, d }, { 1, sa, sb, sd }, x);
        const complex<double> sy14 = G({ 0, b / c, a / c, d / c }, 1);
        const complex<double> sy15 = G({ 0, b / c, d / c, a / c }, 1);
        const complex<double> sy16 = G({ 0, d / c, b / c, a / c }, 1);
        const complex<double> sy17 = G({ b / c, 0, a / c, d / c }, 1);
        const complex<double> sy18 = G({ b / c, 0, d / c, a / c }, 1);
        const complex<double> sy19 = G({ b / c, d / c, 0, a / c }, 1);
        const complex<double> sy20 = G({ d / c, 0, b / c, a / c }, 1);
        const complex<double> sy21 = G({ d / c, b / c, 0, a / c }, 1);
        complex<double> res { sy2 * (sy1 + 2. * sy3) + (2. * sy3 + sy4) * sy5
            + sy6 * (-sy7 - sy8)
            + sy10
                * (-(sy1 * sy12) - sy14 - sy15 - sy16 - sy17 - sy18 - sy19 - sy20 - sy21
                    + (sy2 + sy5) * sy9)
            + sy1 * G({ 0, d / c }, 1)
            + sy9
                * (-sy7 - sy8 - G({ 0, b / c, d / c }, 1) - G({ 0, d / c, b / c }, 1)
                    - G({ b / c, 0, d / c }, 1) - G({ d / c, 0, b / c }, 1))
            + G({ 0, 0, b / c, a / c, d / c }, 1) + G({ 0, 0, b / c, d / c, a / c }, 1)
            + G({ 0, 0, d / c, b / c, a / c }, 1) + G({ 0, b / c, 0, a / c, d / c }, 1)
            + G({ 0, b / c, 0, d / c, a / c }, 1) + G({ 0, b / c, d / c, 0, a / c }, 1)
            + G({ 0, d / c, 0, b / c, a / c }, 1) + G({ 0, d / c, b / c, 0, a / c }, 1)
            + G({ b / c, 0, 0, a / c, d / c }, 1) + G({ b / c, 0, 0, d / c, a / c }, 1)
            + G({ b / c, 0, d / c, 0, a / c }, 1) - G({ b / c, a / c, 0, 0, d / c }, 1)
            - G({ b / c, a / c, 0, d / c, x / c }, 1)
            - G({ b / c, a / c, d / c, 0, x / c }, 1)
            + G({ b / c, d / c, 0, 0, a / c }, 1)
            - G({ b / c, d / c, a / c, 0, x / c }, 1)
            + G({ d / c, 0, 0, b / c, a / c }, 1) + G({ d / c, 0, b / c, 0, a / c }, 1)
            + G({ d / c, b / c, 0, 0, a / c }, 1)
            - G({ d / c, b / c, a / c, 0, x / c }, 1)
            + G({ b / c, a / c }, 1) * (-sy4 + G({ 0, 0, d }, { 1, 1, sd }, x))
            + sy12
                * (-sy13 - 2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                    - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
            + G({ b / c }, 1) * (sy13 - G({ 0, a, 0, d }, { 1, sa, 1, sd }, x))
            + G({ 0, a, b, 0, d }, { 1, sa, sb, 1, sd }, x)
            + (sy1 * sy12 + sy14 + sy15 + sy16 + sy17 + sy18 + sy19 + sy20 + sy21
                  + sy10 * (-sy11 - sy7 - sy8) + (-sy2 - sy5) * sy9)
                * Log(-x, sc)
            + (sy11 / 2. + sy7 / 2. + sy8 / 2.) * pow(sy10, 2.)
            + sy11 * (-sy6 - G({ 0, d }, { 1, sd }, x) - 2. * Zeta(2))
            - 2. * sy7 * Zeta(2) - 2. * sy8 * Zeta(2) };
        if (d != x) {
            res += (G({ b / c, a / c, 0, d / c }, 1) - G({ b / c, a / c, 0, x / c }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G5_0abcd_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;

    if (b == c) {
        if (b == d) { // G0abbb

            const complex<double> sy1 = G({ b, b }, { sb, sb }, x);
            const complex<double> sy2 = G({ a / b, 0, 1 }, 1);
            const complex<double> sy3 = G({ a / b, 1, 1 }, 1);

            complex<double> tmp { -(sy1 * sy2) + sy1 * G({ a / b, 0, x / b }, 1)
                - G({ a / b, 0, 0, 1, 1 }, 1) + G({ a / b, 0, x / b, 1, 1 }, 1)
                + (G({ a / b, 0, 1, 1 }, 1) - G({ a / b, 0, x / b, 1 }, 1))
                    * G({ b }, { sb }, x)
                + (sy2 - sy3) * G({ 0, b }, { 1, sb }, x)
                + G({ a / b, 1 }, 1)
                    * (-G({ 0, a, b }, { 1, sa, sb }, x)
                        + G({ 0, b, b }, { 1, sb, sb }, x))
                + G({ a / b }, 1)
                    * (-G({ 0, 0, b, b }, { 1, 1, sb, sb }, x)
                        + G({ 0, a, b, b }, { 1, sa, sb, sb }, x))
                + G({ 0, a, 0, b, b }, { 1, sa, 1, sb, sb }, x)
                + G({ 0, a }, { 1, sa }, x) * (sy3 - Zeta(3)) };
            return tmp;
        } else { // G0abbd
            const complex<double> sy1 = G({ d / b, 1 }, 1);
            const complex<double> sy2 = G({ 0, a, d }, { 1, sa, sd }, x);
            const complex<double> sy3 = G({ b, d }, { sb, sd }, x);
            const complex<double> sy4 = G({ a / b, 0, 1 }, 1);
            const complex<double> sy5 = G({ a / b, d / b, 1 }, 1);
            const complex<double> sy6 = G({ 0, 0 }, { 1, 1 }, x);
            const complex<double> sy7 = G({ d / b, 1, a / b }, 1);
            const complex<double> sy8 = G({ d / b, a / b, 1 }, 1);
            const complex<double> sy9 = G({ 0, a }, { 1, sa }, x);
            const complex<double> sy10 = Log(b, sb);
            const complex<double> sy11 = G({ 0, a / b, d / b, 1 }, 1);
            const complex<double> sy12 = G({ 0, d / b, 1, a / b }, 1);
            const complex<double> sy13 = G({ 0, d / b, a / b, 1 }, 1);
            const complex<double> sy14 = G({ d / b, 0, 1, a / b }, 1);
            const complex<double> sy15 = G({ d / b, 0, a / b, 1 }, 1);
            const complex<double> sy16 = G({ d / b, 1, 0, a / b }, 1);

            complex<double> res { -(sy3 * sy4) + sy6 * (-sy7 - sy8)
                + sy10 * (-sy11 - sy12 - sy13 - sy14 - sy15 - sy16 + sy1 * sy9)
                + sy9 * (-sy7 - sy8 - G({ 0, d / b, 1 }, 1))
                + sy3 * G({ a / b, 0, x / b }, 1) + G({ 0, 0, a / b, d / b, 1 }, 1)
                + G({ 0, 0, d / b, 1, a / b }, 1) + G({ 0, 0, d / b, a / b, 1 }, 1)
                + G({ 0, d / b, 0, 1, a / b }, 1) + G({ 0, d / b, 0, a / b, 1 }, 1)
                + G({ 0, d / b, 1, 0, a / b }, 1) - G({ a / b, 0, 0, d / b, 1 }, 1)
                - G({ a / b, 0, d / b, 1, x / b }, 1)
                - G({ a / b, 0, d / b, x / b, 1 }, 1)
                - G({ a / b, d / b, 0, 1, x / b }, 1)
                - G({ a / b, d / b, 0, x / b, 1 }, 1)
                - G({ a / b, d / b, 1, 0, x / b }, 1) + G({ d / b, 0, 0, 1, a / b }, 1)
                + G({ d / b, 0, 0, a / b, 1 }, 1) + G({ d / b, 0, 1, 0, a / b }, 1)
                + G({ d / b, 1, 0, 0, a / b }, 1) - G({ d / b, 1, a / b, 0, x / b }, 1)
                - G({ d / b, a / b, 0, 1, x / b }, 1)
                - G({ d / b, a / b, 0, x / b, 1 }, 1)
                - G({ d / b, a / b, 1, 0, x / b }, 1)
                + (sy4 - sy5) * G({ 0, d }, { 1, sd }, x)
                + sy1 * (sy2 + 2. * G({ 0, 0, a }, { 1, 1, sa }, x))
                + G({ a / b, 1 }, 1) * (-sy2 + G({ 0, b, d }, { 1, sb, sd }, x))
                + G({ a / b }, 1)
                    * (-G({ 0, 0, b, d }, { 1, 1, sb, sd }, x)
                        + G({ 0, a, b, d }, { 1, sa, sb, sd }, x))
                + G({ 0, a, 0, b, d }, { 1, sa, 1, sb, sd }, x)
                + (sy11 + sy12 + sy13 + sy14 + sy15 + sy16 + sy10 * (-sy5 - sy7 - sy8)
                      - sy1 * sy9)
                    * Log(-x, sb)
                + (sy5 / 2. + sy7 / 2. + sy8 / 2.) * pow(sy10, 2.)
                + sy5 * (-sy6 - 2. * Zeta(2)) - 2. * sy7 * Zeta(2) - 2. * sy8 * Zeta(2) };
            if (d != x) {
                res += (G({ a / b, 0, d / b, 1 }, 1) - G({ a / b, 0, x / b, 1 }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        }
    } else { // G0abcd can contain G0abcb.
        const complex<double> sy1 = G({ c / b, a / b }, 1);
        const complex<double> sy2 = G({ 0, 0, d }, { 1, 1, sd }, x);
        const complex<double> sy3 = G({ c / b, d / b }, 1);
        const complex<double> sy4 = G({ 0, a, d }, { 1, sa, sd }, x);
        const complex<double> sy5 = G({ a / b, 0, c / b }, 1);
        const complex<double> sy6 = G({ a / b, c / b, d / b }, 1);
        const complex<double> sy7 = G({ c / b, a / b, d / b }, 1);
        const complex<double> sy8 = Log(b, sb);
        const complex<double> sy9 = G({ c / b, d / b, a / b }, 1);
        const complex<double> sy10 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy11 = G({ 0, 0 }, { 1, 1 }, x);
        const complex<double> sy12 = G({ 0, a, c, d }, { 1, sa, sc, sd }, x);
        const complex<double> sy13 = G({ 0, a / b, c / b, d / b }, 1);
        const complex<double> sy14 = G({ 0, c / b, a / b, d / b }, 1);
        const complex<double> sy15 = G({ 0, c / b, d / b, a / b }, 1);
        const complex<double> sy16 = G({ c / b, 0, a / b, d / b }, 1);
        const complex<double> sy17 = G({ c / b, 0, d / b, a / b }, 1);
        const complex<double> sy18 = G({ c / b, d / b, 0, a / b }, 1);
        complex<double> res { -(sy1 * sy2) + sy1 * sy4
            + (sy13 + sy14 + sy15 + sy16 + sy17 + sy18 - sy10 * sy3) * sy8
            + sy11 * (sy7 + sy9)
            + sy10 * (sy9 + G({ 0, c / b, d / b }, 1) + G({ c / b, 0, d / b }, 1))
            - G({ 0, 0, a / b, c / b, d / b }, 1) - G({ 0, 0, c / b, a / b, d / b }, 1)
            - G({ 0, 0, c / b, d / b, a / b }, 1) - G({ 0, c / b, 0, a / b, d / b }, 1)
            - G({ 0, c / b, 0, d / b, a / b }, 1) - G({ 0, c / b, d / b, 0, a / b }, 1)
            + G({ a / b, 0, 0, c / b, d / b }, 1) + G({ a / b, 0, c / b, 0, d / b }, 1)
            + G({ a / b, 0, c / b, d / b, x / b }, 1)
            + G({ a / b, c / b, 0, 0, d / b }, 1)
            + G({ a / b, c / b, 0, d / b, x / b }, 1)
            + G({ a / b, c / b, d / b, 0, x / b }, 1)
            - G({ c / b, 0, 0, a / b, d / b }, 1) - G({ c / b, 0, 0, d / b, a / b }, 1)
            - G({ c / b, 0, d / b, 0, a / b }, 1) + G({ c / b, a / b, 0, 0, d / b }, 1)
            + G({ c / b, a / b, 0, d / b, x / b }, 1)
            + G({ c / b, a / b, d / b, 0, x / b }, 1)
            - G({ c / b, d / b, 0, 0, a / b }, 1)
            + G({ c / b, d / b, a / b, 0, x / b }, 1)
            + (sy5 + sy6 + sy7) * G({ 0, d }, { 1, sd }, x)
            + sy3 * (-sy4 - 2. * G({ 0, 0, a }, { 1, 1, sa }, x))
            + G({ a / b, c / b }, 1) * (-sy2 + G({ 0, c, d }, { 1, sc, sd }, x))
            + G({ a / b }, 1) * (sy12 - G({ 0, 0, c, d }, { 1, 1, sc, sd }, x))
            + G({ c / b }, 1) * (-sy12 + G({ 0, a, 0, d }, { 1, sa, 1, sd }, x))
            + G({ 0, a, 0, c, d }, { 1, sa, 1, sc, sd }, x)
            + (-sy13 - sy14 - sy15 - sy16 - sy17 - sy18 + sy10 * sy3
                  + sy8 * (sy6 + sy7 + sy9))
                * Log(-x, sb)
            + (-sy6 / 2. - sy7 / 2. - sy9 / 2.) * pow(sy8, 2.) + 2. * sy7 * Zeta(2)
            + 2. * sy9 * Zeta(2) + sy6 * (sy11 + 2. * Zeta(2)) };
        if (c != x) {
            res += (-sy5 + G({ a / b, 0, x / b }, 1)) * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (-G({ a / b, 0, c / b, d / b }, 1) + G({ a / b, 0, c / b, x / b }, 1)
                       - G({ a / b, c / b, 0, d / b }, 1)
                       + G({ a / b, c / b, 0, x / b }, 1)
                       - G({ c / b, a / b, 0, d / b }, 1)
                       + G({ c / b, a / b, 0, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G5_0abcd_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;

    if (a == b) {
        if (a == c) {

            const complex<double> sy1 = G({ a, a, d }, { sa, sa, sd }, x);
            const complex<double> sy2 = Log(a, sa);
            const complex<double> sy3 = G({ d / a, 1, 1 }, 1);
            const complex<double> sy4 = G({ 0, d }, { 1, sd }, x);
            complex<double> res { -(sy1 * G({ 0, x / a }, 1))
                - G({ 0, 0, d / a, 1, 1 }, 1) - G({ 0, d / a, 1, 1, x / a }, 1)
                - G({ 0, d / a, 1, x / a, 1 }, 1) - G({ 0, d / a, x / a, 1, 1 }, 1)
                - G({ d / a, 0, 1, 1, x / a }, 1) - G({ d / a, 0, 1, x / a, 1 }, 1)
                - G({ d / a, 0, x / a, 1, 1 }, 1) - G({ d / a, 1, 0, 1, x / a }, 1)
                - G({ d / a, 1, 0, x / a, 1 }, 1) - G({ d / a, 1, 1, 0, x / a }, 1)
                + G({ 0, 0, a, a, d }, { 1, 1, sa, sa, sd }, x) - sy2 * sy3 * Log(-x, sa)
                + (sy3 * pow(sy2, 2.)) / 2.
                + sy3 * (-sy4 - G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) - sy1 * Zeta(2)
                + G({ 0, a, d }, { 1, sa, sd }, x) * Zeta(2)
                + G({ a, d }, { sa, sd }, x) * (G({ 0, x / a, 1 }, 1) - Zeta(3))
                + sy4 * Zeta(3) };
            if (d != x) {
                res += (G({ 0, d / a, 1, 1 }, 1) - G({ 0, x / a, 1, 1 }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        } else {

            const complex<double> sy1 = G({ 0, c, d }, { 1, sc, sd }, x);
            const complex<double> sy2 = G({ a, c, d }, { sa, sc, sd }, x);
            const complex<double> sy3 = Log(a, sa);
            const complex<double> sy4 = G({ c / a, 1, d / a }, 1);
            const complex<double> sy5 = G({ c / a, d / a, 1 }, 1);
            const complex<double> sy6 = G({ 0, 0 }, { 1, 1 }, x);
            const complex<double> sy7 = G({ 0, c / a, 1 }, 1);
            complex<double> res { sy5 * sy6 - sy2 * G({ 0, x / a }, 1)
                + G({ 0, 0, c / a, 1, d / a }, 1) + G({ 0, 0, c / a, d / a, 1 }, 1)
                + G({ 0, c / a, 0, 1, d / a }, 1) + G({ 0, c / a, 0, d / a, 1 }, 1)
                + G({ 0, c / a, 1, 0, d / a }, 1) + G({ 0, c / a, 1, d / a, x / a }, 1)
                + G({ 0, c / a, d / a, 1, x / a }, 1)
                + G({ 0, c / a, d / a, x / a, 1 }, 1) + G({ c / a, 0, 0, 1, d / a }, 1)
                + G({ c / a, 0, 0, d / a, 1 }, 1) + G({ c / a, 0, 1, 0, d / a }, 1)
                + G({ c / a, 0, 1, d / a, x / a }, 1)
                + G({ c / a, 0, d / a, 1, x / a }, 1)
                + G({ c / a, 0, d / a, x / a, 1 }, 1) + G({ c / a, 1, 0, 0, d / a }, 1)
                + G({ c / a, 1, 0, d / a, x / a }, 1)
                + G({ c / a, 1, d / a, 0, x / a }, 1)
                + G({ c / a, d / a, 0, 1, x / a }, 1)
                + G({ c / a, d / a, 0, x / a, 1 }, 1)
                + G({ c / a, d / a, 1, 0, x / a }, 1)
                + (sy4 + sy5 + sy7) * G({ 0, d }, { 1, sd }, x)
                + G({ c / a, 1 }, 1) * (sy1 - G({ 0, 0, d }, { 1, 1, sd }, x))
                + G({ 0, 0, a, c, d }, { 1, 1, sa, sc, sd }, x)
                + sy3 * (sy4 + sy5) * Log(-x, sa) + (-sy4 / 2. - sy5 / 2.) * pow(sy3, 2.)
                + sy1 * Zeta(2) - sy2 * Zeta(2) + 2. * sy5 * Zeta(2)
                + sy4 * (sy6 + 2. * Zeta(2)) };
            if (c != x) {
                res += (-sy7 + G({ 0, x / a, 1 }, 1)) * G({ c, d }, { sc, sd }, x);
            }
            if (d != x) {
                res += (-G({ 0, c / a, 1, d / a }, 1) + G({ 0, c / a, 1, x / a }, 1)
                           - G({ 0, c / a, d / a, 1 }, 1) + G({ 0, c / a, x / a, 1 }, 1)
                           - G({ c / a, 0, 1, d / a }, 1) + G({ c / a, 0, 1, x / a }, 1)
                           - G({ c / a, 0, d / a, 1 }, 1) + G({ c / a, 0, x / a, 1 }, 1)
                           - G({ c / a, 1, 0, d / a }, 1) + G({ c / a, 1, 0, x / a }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        }
    } else {
        const complex<double> sy1 = G({ 0, c, d }, { 1, sc, sd }, x);
        const complex<double> sy2 = G({ 0, b / a }, 1);
        const complex<double> sy3 = G({ 0, b / a, c / a }, 1);
        const complex<double> sy4 = G({ b / a, 0, c / a }, 1);
        const complex<double> sy5 = G({ b / a, c / a, d / a }, 1);
        const complex<double> sy6 = Log(a, sa);
        complex<double> res { -(sy1 * sy2) - G({ 0, 0, b / a, c / a, d / a }, 1)
            - G({ 0, b / a, 0, c / a, d / a }, 1) - G({ 0, b / a, c / a, 0, d / a }, 1)
            - G({ 0, b / a, c / a, d / a, x / a }, 1)
            - G({ b / a, 0, 0, c / a, d / a }, 1) - G({ b / a, 0, c / a, 0, d / a }, 1)
            - G({ b / a, 0, c / a, d / a, x / a }, 1)
            - G({ b / a, c / a, 0, 0, d / a }, 1)
            - G({ b / a, c / a, 0, d / a, x / a }, 1)
            - G({ b / a, c / a, d / a, 0, x / a }, 1)
            + (-sy3 - sy4 - sy5) * G({ 0, d }, { 1, sd }, x)
            + G({ b / a, c / a }, 1) * (-sy1 + G({ 0, 0, d }, { 1, 1, sd }, x))
            + G({ b / a }, 1)
                * (G({ 0, 0, c, d }, { 1, 1, sc, sd }, x)
                    - G({ 0, b, c, d }, { 1, sb, sc, sd }, x))
            + G({ 0, 0, b, c, d }, { 1, 1, sb, sc, sd }, x) - sy5 * sy6 * Log(-x, sa)
            + (sy5 * pow(sy6, 2.)) / 2.
            + sy5 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
        if (b != x) {
            res += (sy2 - G({ 0, x / a }, 1)) * G({ b, c, d }, { sb, sc, sd }, x);
        }
        if (c != x) {
            res += (sy3 + sy4 - G({ 0, b / a, x / a }, 1) - G({ b / a, 0, x / a }, 1))
                * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (G({ 0, b / a, c / a, d / a }, 1) - G({ 0, b / a, c / a, x / a }, 1)
                       + G({ b / a, 0, c / a, d / a }, 1)
                       - G({ b / a, 0, c / a, x / a }, 1)
                       + G({ b / a, c / a, 0, d / a }, 1)
                       - G({ b / a, c / a, 0, x / a }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G5_a0bcd_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;

    const complex<double> sy1 = Log(d, sd);
    const complex<double> sy2 = G({ a }, { sa }, x);
    const complex<double> sy3 = G({ c / d, b / d }, 1);
    const complex<double> sy4 = G({ a, 0, b }, { sa, 1, sb }, x);
    const complex<double> sy5 = G({ 0, a }, { 1, sa }, x);
    const complex<double> sy6 = G({ 0, c / d, b / d }, 1);
    const complex<double> sy7 = G({ c / d, 0, b / d }, 1);
    const complex<double> sy8 = G({ a, 0, b, c }, { sa, 1, sb, sc }, x);
    const complex<double> sy9 = G({ c / d }, 1);
    const complex<double> sy10 = G({ c / d, b / d, 0, a / d }, 1);

    complex<double> res { sy5 * (-sy6 - sy7)
        + sy1 * (-sy10 + sy3 * sy5 + sy2 * (-sy6 - sy7) + sy4 * sy9)
        - sy4 * G({ 0, c / d }, 1) + G({ 0, c / d, b / d, 0, a / d }, 1)
        + G({ c / d, 0, b / d, 0, a / d }, 1) + 2. * G({ c / d, b / d, 0, 0, a / d }, 1)
        + G({ c / d, b / d, 0, a / d, x / d }, 1)
        + sy8 * (-G({ x / d }, 1) + G({ d }, { sd }, x))
        + sy3 * (-sy4 + G({ 0, 0, a }, { 1, 1, sa }, x))
        + sy9
            * (sy8 + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                + 2. * G({ a, 0, 0, b }, { sa, 1, 1, sb }, x))
        - G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x)
        - 2. * G({ a, 0, 0, b, c }, { sa, 1, 1, sb, sc }, x)
        - G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x)
        + (sy10 - sy1 * sy2 * sy3 - sy3 * sy5 + sy2 * (sy6 + sy7) - sy4 * sy9)
            * Log(-x, sd)
        + (sy2 * sy3 * pow(sy1, 2.)) / 2.
        + sy2
            * (-sy10 + G({ 0, 0, c / d, b / d }, 1) + G({ 0, c / d, 0, b / d }, 1)
                + G({ c / d, 0, 0, b / d }, 1)
                + sy3 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };

    return res;
}
complex<double> G5_a0bcd_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (c == d) {

        const complex<double> sy1 = G({ a, 0, b }, { sa, 1, sb }, x);
        const complex<double> sy2 = G({ b / c, 0, a / c }, 1);
        const complex<double> sy3 = G({ b / c, 0, a / c, 1 }, 1);
        const complex<double> sy4 = G({ c }, { sc }, x);
        complex<double> tmp { sy3 * sy4 - sy4 * G({ b / c, 0, a / c, x / c }, 1)
            - G({ b / c, 0, a / c, 0, 1 }, 1) + G({ b / c, 0, a / c, x / c, 1 }, 1)
            + (-sy3 + G({ b / c, 0, 0, 1 }, 1)) * G({ a }, { sa }, x)
            - sy2 * G({ 0, c }, { 1, sc }, x)
            + (sy2 - G({ b / c, 0, 1 }, 1)) * G({ a, c }, { sa, sc }, x)
            + G({ b / c, 1 }, 1) * (-sy1 + G({ a, 0, c }, { sa, 1, sc }, x))
            + G({ b / c }, 1)
                * (-G({ a, 0, 0, c }, { sa, 1, 1, sc }, x)
                    + G({ a, 0, b, c }, { sa, 1, sb, sc }, x))
            + G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x) - sy1 * Zeta(2) };
        return tmp;
    } else {

        const complex<double> sy1 = Log(c, sc);
        const complex<double> sy2 = G({ a }, { sa }, x);
        const complex<double> sy3 = G({ b / c, d / c }, 1);
        const complex<double> sy4 = G({ d / c, b / c }, 1);
        const complex<double> sy5 = G({ a, 0, b }, { sa, 1, sb }, x);
        const complex<double> sy6 = G({ 0, 0, a }, { 1, 1, sa }, x);
        const complex<double> sy7 = G({ b / c, 0, a / c }, 1);
        const complex<double> sy8 = G({ a, d }, { sa, sd }, x);
        const complex<double> sy9 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy10 = G({ 0, b / c, d / c }, 1);
        const complex<double> sy11 = G({ 0, d / c, b / c }, 1);
        const complex<double> sy12 = G({ d / c, 0, b / c }, 1);
        const complex<double> sy13 = G({ d / c }, 1);
        const complex<double> sy14 = G({ a, 0, b, d }, { sa, 1, sb, sd }, x);
        const complex<double> sy15 = G({ b / c, 0, a / c, d / c }, 1);
        const complex<double> sy16 = G({ b / c, 0, d / c, a / c }, 1);
        const complex<double> sy17 = G({ b / c, d / c, 0, a / c }, 1);
        const complex<double> sy18 = G({ d / c, b / c, 0, a / c }, 1);
        const complex<double> sy19 = G({ 0, 0 }, { 1, 1 }, x);
        complex<double> res { sy4 * (sy5 - sy6) + sy7 * sy8 + (sy10 + sy11 + sy12) * sy9
            + sy1
                * (sy15 + sy16 + sy17 + sy18 + (sy10 + sy11 + sy12) * sy2 - sy13 * sy5
                    - sy3 * sy9 - sy4 * sy9)
            + sy5 * G({ 0, d / c }, 1) - sy8 * G({ b / c, 0, d / c }, 1)
            - G({ 0, b / c, 0, a / c, d / c }, 1) - G({ 0, b / c, 0, d / c, a / c }, 1)
            - G({ 0, b / c, d / c, 0, a / c }, 1) - G({ 0, d / c, b / c, 0, a / c }, 1)
            - 2. * G({ b / c, 0, 0, a / c, d / c }, 1)
            - 2. * G({ b / c, 0, 0, d / c, a / c }, 1)
            - G({ b / c, 0, a / c, 0, d / c }, 1)
            - G({ b / c, 0, a / c, d / c, x / c }, 1)
            - 2. * G({ b / c, 0, d / c, 0, a / c }, 1)
            - G({ b / c, 0, d / c, a / c, x / c }, 1)
            - 2. * G({ b / c, d / c, 0, 0, a / c }, 1)
            - G({ b / c, d / c, 0, a / c, x / c }, 1)
            - G({ d / c, 0, b / c, 0, a / c }, 1)
            - 2. * G({ d / c, b / c, 0, 0, a / c }, 1)
            - G({ d / c, b / c, 0, a / c, x / c }, 1) - sy7 * G({ 0, d }, { 1, sd }, x)
            + sy3 * (-sy6 + G({ a, 0, d }, { sa, 1, sd }, x))
            + sy13
                * (-sy14 - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                    - 2. * G({ a, 0, 0, b }, { sa, 1, 1, sb }, x))
            + G({ b / c }, 1) * (sy14 - G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
            + G({ a, 0, b, 0, d }, { sa, 1, sb, 1, sd }, x)
            + (-sy15 - sy16 - sy17 - sy18 + (-sy10 - sy11 - sy12) * sy2
                  + sy1 * sy2 * (sy3 + sy4) + sy13 * sy5 + sy3 * sy9 + sy4 * sy9)
                * Log(-x, sc)
            + sy2 * (-sy3 / 2. - sy4 / 2.) * pow(sy1, 2.)
            + sy2
                * (sy16 + sy17 + sy18 + sy19 * sy4 - G({ 0, 0, b / c, d / c }, 1)
                    - G({ 0, 0, d / c, b / c }, 1) - G({ 0, d / c, 0, b / c }, 1)
                    + G({ b / c, 0, 0, d / c }, 1) - G({ d / c, 0, 0, b / c }, 1)
                    + 2. * sy4 * Zeta(2) + sy3 * (sy19 + 2. * Zeta(2))) };
        if (d != x) {
            res += (sy15 - G({ b / c, 0, a / c, x / c }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G5_a0bcd_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;

    if (b == c) {
        if (b == d) {

            const complex<double> sy1 = G({ 0, a / b, 1 }, 1);
            const complex<double> sy2 = G({ 0, a / b }, 1);
            const complex<double> sy3 = G({ a, b, b }, { sa, sb, sb }, x);
            const complex<double> sy4 = G({ a }, { sa }, x);
            const complex<double> sy5 = G({ 0, a / b, 1, 1 }, 1);
            const complex<double> sy6 = G({ b }, { sb }, x);

            complex<double> tmp { -(sy2 * sy3) - sy4 * sy5 + sy5 * sy6
                - sy6 * G({ 0, a / b, x / b, 1 }, 1) - G({ 0, a / b, 0, 1, 1 }, 1)
                + G({ 0, a / b, x / b, 1, 1 }, 1)
                + (-sy1 + G({ 0, a / b, x / b }, 1)) * G({ b, b }, { sb, sb }, x)
                + sy2 * G({ 0, b, b }, { 1, sb, sb }, x)
                + G({ a, 0, 0, b, b }, { sa, 1, 1, sb, sb }, x) - sy3 * Zeta(2)
                + G({ a, 0, b }, { sa, 1, sb }, x) * Zeta(2) + (sy4 * Zeta(4)) / 4.
                + G({ a, b }, { sa, sb }, x) * (sy1 - Zeta(3)) };
            return tmp;
        } else {

            const complex<double> sy1 = Log(b, sb);
            const complex<double> sy2 = G({ a }, { sa }, x);
            const complex<double> sy3 = G({ d / b, 1 }, 1);
            const complex<double> sy4 = G({ b, d }, { sb, sd }, x);
            const complex<double> sy5 = G({ 0, a / b, 1 }, 1);
            const complex<double> sy6 = G({ a, 0, d }, { sa, 1, sd }, x);
            const complex<double> sy7 = G({ a, b, d }, { sa, sb, sd }, x);
            const complex<double> sy8 = G({ 0, a }, { 1, sa }, x);
            const complex<double> sy9 = G({ 0, a / b, d / b, 1 }, 1);
            const complex<double> sy10 = G({ 0, d / b, 1, a / b }, 1);
            const complex<double> sy11 = G({ 0, d / b, a / b, 1 }, 1);
            const complex<double> sy12 = G({ d / b, 0, 1, a / b }, 1);
            const complex<double> sy13 = G({ d / b, 0, a / b, 1 }, 1);
            const complex<double> sy14 = G({ d / b, 1, 0, a / b }, 1);
            complex<double> res { -(sy4 * sy5)
                + sy1 * (sy10 + sy11 + sy12 + sy13 + sy14 - sy3 * sy8 + sy9)
                + sy4 * G({ 0, a / b, x / b }, 1) - 2. * G({ 0, 0, a / b, d / b, 1 }, 1)
                - 2. * G({ 0, 0, d / b, 1, a / b }, 1)
                - 2. * G({ 0, 0, d / b, a / b, 1 }, 1) - G({ 0, a / b, 0, d / b, 1 }, 1)
                - G({ 0, a / b, d / b, 1, x / b }, 1)
                - G({ 0, a / b, d / b, x / b, 1 }, 1)
                - 2. * G({ 0, d / b, 0, 1, a / b }, 1)
                - 2. * G({ 0, d / b, 0, a / b, 1 }, 1)
                - 2. * G({ 0, d / b, 1, 0, a / b }, 1)
                - G({ 0, d / b, 1, a / b, x / b }, 1)
                - G({ 0, d / b, a / b, 1, x / b }, 1)
                - G({ 0, d / b, a / b, x / b, 1 }, 1)
                - 2. * G({ d / b, 0, 0, 1, a / b }, 1)
                - 2. * G({ d / b, 0, 0, a / b, 1 }, 1)
                - 2. * G({ d / b, 0, 1, 0, a / b }, 1)
                - G({ d / b, 0, 1, a / b, x / b }, 1)
                - G({ d / b, 0, a / b, 1, x / b }, 1)
                - G({ d / b, 0, a / b, x / b, 1 }, 1)
                - 2. * G({ d / b, 1, 0, 0, a / b }, 1)
                - G({ d / b, 1, 0, a / b, x / b }, 1)
                + (sy5 - G({ 0, d / b, 1 }, 1)) * G({ a, d }, { sa, sd }, x)
                + sy3 * (sy6 - G({ 0, 0, a }, { 1, 1, sa }, x))
                + G({ 0, a / b }, 1) * (-sy7 + G({ 0, b, d }, { 1, sb, sd }, x))
                + G({ a, 0, 0, b, d }, { sa, 1, 1, sb, sd }, x)
                + (-sy10 - sy11 - sy12 - sy13 - sy14 + sy1 * sy2 * sy3 + sy3 * sy8 - sy9)
                    * Log(-x, sb)
                - (sy2 * sy3 * pow(sy1, 2.)) / 2. + sy6 * Zeta(2) - sy7 * Zeta(2)
                + sy2
                    * (sy10 + sy11 + sy12 + sy13 + sy14 + G({ 0, 0, d / b, 1 }, 1)
                        + sy3 * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2))) };
            if (d != x) {
                res += (sy9 - G({ 0, a / b, x / b, 1 }, 1)) * G({ d }, { sd }, x);
            }
            return res;
        }
    } else {

        const complex<double> sy1 = Log(b, sb);
        const complex<double> sy2 = G({ a }, { sa }, x);
        const complex<double> sy3 = G({ c / b, d / b }, 1);
        const complex<double> sy4 = G({ a, d }, { sa, sd }, x);
        const complex<double> sy5 = G({ 0, c / b, a / b }, 1);
        const complex<double> sy6 = G({ a, 0, d }, { sa, 1, sd }, x);
        const complex<double> sy7 = G({ a, c, d }, { sa, sc, sd }, x);
        const complex<double> sy8 = G({ 0, a / b, c / b }, 1);
        const complex<double> sy9 = G({ c / b, 0, a / b }, 1);
        const complex<double> sy10 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy11 = G({ 0, a / b, c / b, d / b }, 1);
        const complex<double> sy12 = G({ 0, c / b, a / b, d / b }, 1);
        const complex<double> sy13 = G({ 0, c / b, d / b, a / b }, 1);
        const complex<double> sy14 = G({ c / b, 0, a / b, d / b }, 1);
        const complex<double> sy15 = G({ c / b, 0, d / b, a / b }, 1);
        const complex<double> sy16 = G({ c / b, d / b, 0, a / b }, 1);
        complex<double> res { sy1
                * (-sy11 - sy12 - sy13 - sy14 - sy15 - sy16 + sy10 * sy3)
            - sy4 * sy5 + (-sy6 + sy7) * G({ 0, c / b }, 1)
            + sy4 * (-sy9 + G({ 0, c / b, d / b }, 1) + G({ c / b, 0, d / b }, 1))
            + 2. * G({ 0, 0, a / b, c / b, d / b }, 1)
            + 2. * G({ 0, 0, c / b, a / b, d / b }, 1)
            + 2. * G({ 0, 0, c / b, d / b, a / b }, 1)
            + G({ 0, a / b, 0, c / b, d / b }, 1) + G({ 0, a / b, c / b, 0, d / b }, 1)
            + G({ 0, a / b, c / b, d / b, x / b }, 1)
            + 2. * G({ 0, c / b, 0, a / b, d / b }, 1)
            + 2. * G({ 0, c / b, 0, d / b, a / b }, 1)
            + G({ 0, c / b, a / b, 0, d / b }, 1)
            + G({ 0, c / b, a / b, d / b, x / b }, 1)
            + 2. * G({ 0, c / b, d / b, 0, a / b }, 1)
            + G({ 0, c / b, d / b, a / b, x / b }, 1)
            + 2. * G({ c / b, 0, 0, a / b, d / b }, 1)
            + 2. * G({ c / b, 0, 0, d / b, a / b }, 1)
            + G({ c / b, 0, a / b, 0, d / b }, 1)
            + G({ c / b, 0, a / b, d / b, x / b }, 1)
            + 2. * G({ c / b, 0, d / b, 0, a / b }, 1)
            + G({ c / b, 0, d / b, a / b, x / b }, 1)
            + 2. * G({ c / b, d / b, 0, 0, a / b }, 1)
            + G({ c / b, d / b, 0, a / b, x / b }, 1)
            + (sy5 + sy8 + sy9) * G({ 0, d }, { 1, sd }, x)
            + sy3 * (-sy6 + G({ 0, 0, a }, { 1, 1, sa }, x))
            + G({ 0, a / b }, 1) * (-sy7 + G({ 0, c, d }, { 1, sc, sd }, x))
            + G({ c / b }, 1)
                * (G({ a, 0, 0, d }, { sa, 1, 1, sd }, x)
                    - G({ a, 0, c, d }, { sa, 1, sc, sd }, x))
            + G({ a, 0, 0, c, d }, { sa, 1, 1, sc, sd }, x)
            + (sy11 + sy12 + sy13 + sy14 + sy15 + sy16 - sy10 * sy3 - sy1 * sy2 * sy3)
                * Log(-x, sb)
            + (sy2 * sy3 * pow(sy1, 2.)) / 2.
            + sy2
                * (-sy13 - sy15 - sy16 - G({ 0, 0, c / b, d / b }, 1)
                    - G({ 0, c / b, 0, d / b }, 1) - G({ c / b, 0, 0, d / b }, 1)
                    + sy3 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
        if (c != x) {
            res += (-sy8 + G({ 0, a / b, x / b }, 1)) * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (-sy11 - sy12 - sy14 + G({ 0, a / b, c / b, x / b }, 1)
                       + G({ 0, c / b, a / b, x / b }, 1)
                       + G({ c / b, 0, a / b, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G5_a0bcd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{

    if (a == b) {
        if (a == c) {

            const complex<double> sy1 = G({ a, a, d }, { sa, sa, sd }, x);
            const complex<double> sy2 = G({ 0, 1, 1, d / a }, 1);
            complex<double> res { sy1 * G({ 0, x / a }, 1)
                + 2. * G({ 0, 0, 1, 1, d / a }, 1) + G({ 0, 1, 0, 1, d / a }, 1)
                + G({ 0, 1, 1, 0, d / a }, 1) + G({ 0, 1, 1, d / a, x / a }, 1)
                + G({ x / a }, 1) * G({ 0, a, a, d }, { 1, sa, sa, sd }, x)
                + G({ 0, 0, a, a, d }, { 1, 1, sa, sa, sd }, x) - sy2 * Log(a, sa)
                + sy2 * Log(-x, sa) + sy1 * Zeta(2)
                - G({ 0, a, d }, { 1, sa, sd }, x) * Zeta(2)
                + G({ a, d }, { sa, sd }, x) * (G({ 0, 1, x / a }, 1) - Zeta(3))
                + G({ 0, d }, { 1, sd }, x) * Zeta(3) };
            if (d != x) {
                res += (-sy2 + G({ 0, 1, 1, x / a }, 1)) * G({ d }, { sd }, x);
            }
            return res;
        } else {

            const complex<double> sy1 = G({ 0, 1, c / a }, 1);
            const complex<double> sy2 = G({ a, c, d }, { sa, sc, sd }, x);
            const complex<double> sy3 = G({ 0, 1, c / a, d / a }, 1);
            complex<double> res { sy2 * G({ 0, x / a }, 1)
                + 2. * G({ 0, 0, 1, c / a, d / a }, 1) + G({ 0, 1, 0, c / a, d / a }, 1)
                + G({ 0, 1, c / a, 0, d / a }, 1) + G({ 0, 1, c / a, d / a, x / a }, 1)
                + sy1 * G({ 0, d }, { 1, sd }, x)
                + G({ x / a }, 1) * G({ 0, a, c, d }, { 1, sa, sc, sd }, x)
                + G({ 0, 0, a, c, d }, { 1, 1, sa, sc, sd }, x) - sy3 * Log(a, sa)
                + sy3 * Log(-x, sa) + sy2 * Zeta(2)
                - G({ 0, c, d }, { 1, sc, sd }, x) * Zeta(2) };
            if (c != x) {
                res += (-sy1 + G({ 0, 1, x / a }, 1)) * G({ c, d }, { sc, sd }, x);
            }
            if (d != x) {
                res += (-sy3 + G({ 0, 1, c / a, x / a }, 1)) * G({ d }, { sd }, x);
            }
            return res;
        }
    } else {

        const complex<double> sy1 = G({ 0, b / a, c / a }, 1);
        const complex<double> sy2 = G({ 0, b / a }, 1);
        const complex<double> sy3 = G({ 0, b / a, c / a, d / a }, 1);
        complex<double> res { 2. * G({ 0, 0, b / a, c / a, d / a }, 1)
            + G({ 0, b / a, 0, c / a, d / a }, 1) + G({ 0, b / a, c / a, 0, d / a }, 1)
            + G({ 0, b / a, c / a, d / a, x / a }, 1) + sy1 * G({ 0, d }, { 1, sd }, x)
            + sy2 * G({ 0, c, d }, { 1, sc, sd }, x)
            + G({ x / a }, 1) * G({ 0, b, c, d }, { 1, sb, sc, sd }, x)
            + G({ 0, 0, b, c, d }, { 1, 1, sb, sc, sd }, x) - sy3 * Log(a, sa)
            + sy3 * Log(-x, sa) };
        if (b != x) {
            res += (-sy2 + G({ 0, x / a }, 1)) * G({ b, c, d }, { sb, sc, sd }, x);
        }
        if (c != x) {
            res += (-sy1 + G({ 0, b / a, x / a }, 1)) * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (-sy3 + G({ 0, b / a, c / a, x / a }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G5_ab0cd_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;

    const complex<double> sy1 = Log(d, sd);
    const complex<double> sy2 = G({ c / d }, 1);
    const complex<double> sy3 = G({ a, b }, { sa, sb }, x);
    const complex<double> sy4 = G({ 0, c / d }, 1);
    const complex<double> sy5 = G({ 0, a, b }, { 1, sa, sb }, x);
    const complex<double> sy6 = G({ a, 0, b }, { sa, 1, sb }, x);
    const complex<double> sy7 = G({ c / d, x / d }, 1);
    const complex<double> sy8 = G({ c / d, 0, b / d }, 1);
    const complex<double> sy9 = G({ c / d, 0, x / d }, 1);
    const complex<double> sy10 = G({ a, b, 0, c }, { sa, sb, 1, sc }, x);
    const complex<double> sy11 = G({ a }, { sa }, x);
    const complex<double> sy12 = G({ c / d, 0, b / d, a / d }, 1);

    complex<double> res { sy4 * (sy5 + sy6)
        + sy1 * (-sy12 + sy3 * sy4 + sy2 * (-sy5 - sy6) + sy11 * sy8)
        + sy3 * (sy9 - G({ 0, 0, c / d }, 1)) + G({ 0, c / d, 0, b / d, a / d }, 1)
        + 2. * G({ c / d, 0, 0, b / d, a / d }, 1) + G({ c / d, 0, b / d, 0, a / d }, 1)
        + G({ c / d, 0, b / d, a / d, x / d }, 1)
        + sy10 * (-G({ x / d }, 1) + G({ d }, { sd }, x))
        + sy8 * G({ 0, a }, { 1, sa }, x)
        + sy11
            * (-sy12 - G({ 0, c / d, 0, b / d }, 1) - 2. * G({ c / d, 0, 0, b / d }, 1)
                - sy7 * G({ 0, b }, { 1, sb }, x))
        + sy7 * (sy5 + sy6 + G({ 0, b, a }, { 1, sb, sa }, x))
        - G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x)
        - G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x)
        - 2. * G({ a, b, 0, 0, c }, { sa, sb, 1, 1, sc }, x)
        + (sy12 + sy1 * sy2 * sy3 - sy3 * sy4 + sy2 * (sy5 + sy6) - sy11 * sy8)
            * Log(-x, sd)
        - (sy2 * sy3 * pow(sy1, 2.)) / 2.
        + sy2
            * (sy10 - sy11 * G({ 0, 0, b }, { 1, 1, sb }, x)
                + G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                + sy3 * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2))) };
    if (b != x) {
        res += (sy11 * sy8 - sy11 * sy9) * G({ b }, { sb }, x)
            + (-sy8 + sy9) * G({ b, a }, { sb, sa }, x);
    }

    return res;
}
complex<double> G5_ab0cd_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (c == d) {

        const complex<double> sy1 = G({ 0, b / c, a / c }, 1);
        const complex<double> sy2 = G({ 0, b / c, 1 }, 1);
        const complex<double> sy3 = G({ x / c, 1 }, 1);
        const complex<double> sy4 = G({ a, b, c }, { sa, sb, sc }, x);
        const complex<double> sy5 = G({ a }, { sa }, x);
        const complex<double> sy6 = G({ 0, b / c, a / c, 1 }, 1);
        const complex<double> sy7 = G({ c }, { sc }, x);
        const complex<double> sy8 = G({ 0, 1, b / c }, 1);
        const complex<double> sy9 = G({ 0, 1, x / c }, 1);
        const complex<double> sy10 = G({ 0, x / c, 1 }, 1);
        complex<double> tmp { sy6 * sy7 - sy7 * G({ 0, b / c, a / c, x / c }, 1)
            - G({ 0, b / c, a / c, 0, 1 }, 1) + G({ 0, b / c, a / c, x / c, 1 }, 1)
            + sy5 * (-sy6 + G({ 0, b / c, 0, 1 }, 1) - sy3 * G({ 0, b }, { 1, sb }, x))
            - sy1 * G({ 0, c }, { 1, sc }, x) + (sy1 - sy2) * G({ a, c }, { sa, sc }, x)
            + sy3
                * (G({ 0, a, b }, { 1, sa, sb }, x) + G({ 0, b, a }, { 1, sb, sa }, x)
                    + G({ a, 0, b }, { sa, 1, sb }, x))
            + G({ 0, b / c }, 1) * (-sy4 + G({ a, 0, c }, { sa, 1, sc }, x))
            + G({ a, b, 0, 0, c }, { sa, sb, 1, 1, sc }, x) - sy4 * Zeta(2)
            + G({ a, b }, { sa, sb }, x) * (sy10 - sy8 + sy9 + Zeta(3)) };
        if (b != x) {
            tmp += (-(sy10 * sy5) + sy2 * sy5 + sy5 * sy8 - sy5 * sy9)
                    * G({ b }, { sb }, x)
                + (sy10 - sy2 - sy8 + sy9) * G({ b, a }, { sb, sa }, x);
        }
        return tmp;
    } else {

        const complex<double> sy1 = Log(c, sc);
        const complex<double> sy2 = G({ d / c }, 1);
        const complex<double> sy3 = G({ a, b }, { sa, sb }, x);
        const complex<double> sy4 = G({ 0, b / c, a / c }, 1);
        const complex<double> sy5 = G({ a, d }, { sa, sd }, x);
        const complex<double> sy6 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy7 = G({ 0, b / c, d / c }, 1);
        const complex<double> sy8 = G({ d / c, x / c }, 1);
        const complex<double> sy9 = G({ 0, a, b }, { 1, sa, sb }, x);
        const complex<double> sy10 = G({ a, 0, b }, { sa, 1, sb }, x);
        const complex<double> sy11 = G({ a, b, d }, { sa, sb, sd }, x);
        const complex<double> sy12 = G({ 0, d / c, b / c }, 1);
        const complex<double> sy13 = G({ d / c, 0, b / c }, 1);
        const complex<double> sy14 = G({ d / c, 0, x / c }, 1);
        const complex<double> sy15 = G({ a }, { sa }, x);
        const complex<double> sy16 = G({ 0, b / c, a / c, d / c }, 1);
        const complex<double> sy17 = G({ 0, b / c, d / c, a / c }, 1);
        const complex<double> sy18 = G({ 0, d / c, b / c, a / c }, 1);
        const complex<double> sy19 = G({ d / c, 0, b / c, a / c }, 1);
        complex<double> res { sy4 * sy5 + (-sy12 - sy13) * sy6 - sy5 * sy7 - sy6 * sy7
            + sy1
                * (sy16 + sy17 + sy18 + sy19 + sy15 * (-sy12 - sy13 - sy7)
                    + sy2 * (sy10 + sy9))
            + sy11 * G({ 0, d / c }, 1) + sy3 * (-sy12 - sy14 - G({ 0, 0, d / c }, 1))
            - 2. * G({ 0, 0, b / c, a / c, d / c }, 1)
            - 2. * G({ 0, 0, b / c, d / c, a / c }, 1)
            - 2. * G({ 0, 0, d / c, b / c, a / c }, 1)
            - G({ 0, b / c, 0, a / c, d / c }, 1) - G({ 0, b / c, 0, d / c, a / c }, 1)
            - G({ 0, b / c, a / c, 0, d / c }, 1)
            - G({ 0, b / c, a / c, d / c, x / c }, 1)
            - G({ 0, b / c, d / c, 0, a / c }, 1)
            - G({ 0, b / c, d / c, a / c, x / c }, 1)
            - 2. * G({ 0, d / c, 0, b / c, a / c }, 1)
            - G({ 0, d / c, b / c, 0, a / c }, 1)
            - G({ 0, d / c, b / c, a / c, x / c }, 1)
            - 2. * G({ d / c, 0, 0, b / c, a / c }, 1)
            - G({ d / c, 0, b / c, 0, a / c }, 1)
            - G({ d / c, 0, b / c, a / c, x / c }, 1)
            + sy15
                * (sy17 + sy18 + sy19 + 2. * G({ 0, 0, b / c, d / c }, 1)
                    + 2. * G({ 0, 0, d / c, b / c }, 1) + G({ 0, b / c, 0, d / c }, 1)
                    + 2. * G({ 0, d / c, 0, b / c }, 1)
                    + 2. * G({ d / c, 0, 0, b / c }, 1) + sy8 * G({ 0, b }, { 1, sb }, x))
            - sy4 * G({ 0, d }, { 1, sd }, x)
            + sy8 * (-sy10 - sy9 - G({ 0, b, a }, { 1, sb, sa }, x))
            + G({ 0, b / c }, 1) * (-sy11 + G({ a, 0, d }, { sa, 1, sd }, x))
            + G({ a, b, 0, 0, d }, { sa, sb, 1, 1, sd }, x)
            + (-sy16 - sy17 - sy18 - sy19 - sy1 * sy2 * sy3 + sy15 * (sy12 + sy13 + sy7)
                  + sy2 * (-sy10 - sy9))
                * Log(-x, sc)
            + (sy2 * sy3 * pow(sy1, 2.)) / 2.
            + sy2
                * (sy15 * G({ 0, 0, b }, { 1, 1, sb }, x)
                    - G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                    - G({ a, b, 0, d }, { sa, sb, 1, sd }, x)
                    + sy3 * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
        if (b != x) {
            res += (-(sy13 * sy15) + sy14 * sy15) * G({ b }, { sb }, x)
                + (sy13 - sy14) * G({ b, a }, { sb, sa }, x);
        }
        if (d != x) {
            res += (sy16 - G({ 0, b / c, a / c, x / c }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G5_ab0cd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    const complex<double> sy1 = G({ 0, c, d }, { 1, sc, sd }, x);
    const complex<double> sy2 = G({ a, d }, { sa, sd }, x);
    const complex<double> sy3 = G({ 0, c / b, a / b }, 1);
    const complex<double> sy4 = G({ 0, c / b, d / b }, 1);
    const complex<double> sy5 = G({ a, c, d }, { sa, sc, sd }, x);
    const complex<double> sy6 = G({ 0, a / b, c / b }, 1);
    const complex<double> sy7 = G({ a / b, 0, c / b }, 1);
    const complex<double> sy8 = G({ a }, { sa }, x);
    const complex<double> sy9 = G({ 0, c / b, d / b, a / b }, 1);
    const complex<double> sy10 = G({ 0, a / b, c / b, d / b }, 1);
    const complex<double> sy11 = G({ 0, c / b, a / b, d / b }, 1);
    const complex<double> sy12 = G({ a / b, 0, c / b, d / b }, 1);
    complex<double> res { sy2 * sy3 - sy2 * sy4 + (-sy1 + sy5) * G({ 0, a / b }, 1)
        - sy1 * G({ a / b, x / b }, 1)
        + sy8 * (sy9 + 2. * G({ 0, 0, c / b, d / b }, 1) + G({ 0, c / b, 0, d / b }, 1))
        - 2. * G({ 0, 0, a / b, c / b, d / b }, 1)
        - 2. * G({ 0, 0, c / b, a / b, d / b }, 1)
        - 2. * G({ 0, 0, c / b, d / b, a / b }, 1)
        - 2. * G({ 0, a / b, 0, c / b, d / b }, 1) - G({ 0, a / b, c / b, 0, d / b }, 1)
        - G({ 0, a / b, c / b, d / b, x / b }, 1) - G({ 0, c / b, 0, a / b, d / b }, 1)
        - G({ 0, c / b, 0, d / b, a / b }, 1) - G({ 0, c / b, a / b, 0, d / b }, 1)
        - G({ 0, c / b, a / b, d / b, x / b }, 1) - G({ 0, c / b, d / b, 0, a / b }, 1)
        - G({ 0, c / b, d / b, a / b, x / b }, 1)
        - 2. * G({ a / b, 0, 0, c / b, d / b }, 1) - G({ a / b, 0, c / b, 0, d / b }, 1)
        - G({ a / b, 0, c / b, d / b, x / b }, 1) - sy4 * G({ 0, a }, { 1, sa }, x)
        + (-sy3 - sy6 - sy7) * G({ 0, d }, { 1, sd }, x)
        + G({ 0, c / b }, 1) * (-sy5 + G({ a, 0, d }, { sa, 1, sd }, x))
        + G({ a / b }, 1)
            * (-G({ 0, 0, c, d }, { 1, 1, sc, sd }, x)
                + G({ a, 0, c, d }, { sa, 1, sc, sd }, x))
        + G({ a, 0, 0, c, d }, { sa, 1, 1, sc, sd }, x)
        + (sy10 + sy11 + sy12 - sy4 * sy8 + sy9) * Log(b, sb)
        + (-sy10 - sy11 - sy12 + sy4 * sy8 - sy9) * Log(-x, sb) };
    if (c != x) {
        res += (sy6 + sy7 - G({ 0, a / b, x / b }, 1) - G({ a / b, 0, x / b }, 1))
            * G({ c, d }, { sc, sd }, x);
    }
    if (d != x) {
        res += (sy10 + sy11 + sy12 - G({ 0, a / b, c / b, x / b }, 1)
                   - G({ 0, c / b, a / b, x / b }, 1) - G({ a / b, 0, c / b, x / b }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}
complex<double> G5_ab0cd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{

    if (a == b) {

        const complex<double> sy1 = G({ 0, 1, c / a }, 1);
        const complex<double> sy2 = G({ 0, c / a, 1 }, 1);
        const complex<double> sy3 = G({ 0, 1, c / a, d / a }, 1);
        const complex<double> sy4 = G({ 0, c / a, 1, d / a }, 1);
        const complex<double> sy5 = G({ 0, c / a, d / a, 1 }, 1);
        complex<double> res { -2. * G({ 0, 0, 1, c / a, d / a }, 1)
            - 2. * G({ 0, 0, c / a, 1, d / a }, 1) - 2. * G({ 0, 0, c / a, d / a, 1 }, 1)
            - G({ 0, 1, 0, c / a, d / a }, 1) - G({ 0, 1, c / a, 0, d / a }, 1)
            - G({ 0, 1, c / a, d / a, x / a }, 1) - G({ 0, c / a, 0, 1, d / a }, 1)
            - G({ 0, c / a, 0, d / a, 1 }, 1) - G({ 0, c / a, 1, 0, d / a }, 1)
            - G({ 0, c / a, 1, d / a, x / a }, 1) - G({ 0, c / a, d / a, 1, x / a }, 1)
            - G({ 0, c / a, d / a, x / a, 1 }, 1)
            + (-sy1 - sy2) * G({ 0, d }, { 1, sd }, x)
            - G({ x / a, 1 }, 1) * G({ 0, c, d }, { 1, sc, sd }, x)
            + G({ x / a }, 1) * G({ a, 0, c, d }, { sa, 1, sc, sd }, x)
            + G({ 0, a, 0, c, d }, { 1, sa, 1, sc, sd }, x)
            + (sy3 + sy4 + sy5) * Log(a, sa) + (-sy3 - sy4 - sy5) * Log(-x, sa) };
        if (c != x) {
            res += (sy1 + sy2 - G({ 0, 1, x / a }, 1) - G({ 0, x / a, 1 }, 1))
                * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (sy3 + sy4 + sy5 - G({ 0, 1, c / a, x / a }, 1)
                       - G({ 0, c / a, 1, x / a }, 1) - G({ 0, c / a, x / a, 1 }, 1))
                * G({ d }, { sd }, x);
        }
        return res;

    } else {

        const complex<double> sy1 = G({ b / a, 0, c / a }, 1);
        const complex<double> sy2 = G({ b / a }, 1);
        const complex<double> sy3 = G({ b / a, 0, c / a, d / a }, 1);
        complex<double> res { G({ 0, b / a, 0, c / a, d / a }, 1)
            + 2. * G({ b / a, 0, 0, c / a, d / a }, 1)
            + G({ b / a, 0, c / a, 0, d / a }, 1)
            + G({ b / a, 0, c / a, d / a, x / a }, 1) + sy1 * G({ 0, d }, { 1, sd }, x)
            + G({ b / a, x / a }, 1) * G({ 0, c, d }, { 1, sc, sd }, x)
            + sy2 * G({ 0, 0, c, d }, { 1, 1, sc, sd }, x)
            + G({ 0, b, 0, c, d }, { 1, sb, 1, sc, sd }, x) - sy3 * Log(a, sa)
            + sy3 * Log(-x, sa) };
        if (b != x) {
            res += (-sy2 + G({ x / a }, 1)) * G({ b, 0, c, d }, { sb, 1, sc, sd }, x);
        }
        if (c != x) {
            res += (-sy1 + G({ b / a, 0, x / a }, 1)) * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (-sy3 + G({ b / a, 0, c / a, x / a }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G5_abc0d_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{

    const complex<double> sy1 = G({ x / d }, 1);
    const complex<double> sy2 = G({ 0, c }, { 1, sc }, x);
    const complex<double> sy3 = G({ a, b }, { sa, sb }, x);
    const complex<double> sy4 = G({ 0, c / d }, 1);
    const complex<double> sy5 = Log(d, sd);
    const complex<double> sy6 = G({ 0, c / d, b / d }, 1);
    const complex<double> sy7 = G({ 0, c / d, x / d }, 1);
    const complex<double> sy8 = G({ d }, { sd }, x);
    const complex<double> sy9 = G({ 0, c, b, a }, { 1, sc, sb, sa }, x);
    const complex<double> sy10 = G({ a }, { sa }, x);
    const complex<double> sy11 = G({ 0, c, b }, { 1, sc, sb }, x);
    const complex<double> sy12 = G({ 0, c / d, b / d, a / d }, 1);
    const complex<double> sy13 = G({ 0, x / d }, 1);

    complex<double> res { sy1 * sy2 * sy3 - sy12 * sy5 + sy8 * (-(sy2 * sy3) - sy9)
        + sy1 * sy9 + 2. * G({ 0, 0, c / d, b / d, a / d }, 1)
        + G({ 0, c / d, 0, b / d, a / d }, 1) + G({ 0, c / d, b / d, 0, a / d }, 1)
        + G({ 0, c / d, b / d, a / d, x / d }, 1) + sy6 * G({ 0, a }, { 1, sa }, x)
        + sy3 * (sy7 + 2. * G({ 0, 0, c / d }, 1) + G({ 0, 0, c }, { 1, 1, sc }, x))
        + sy4 * (-(sy3 * sy5) + G({ 0, b, a }, { 1, sb, sa }, x))
        + G({ 0, d }, { 1, sd }, x) * G({ a, b, c }, { sa, sb, sc }, x)
        + sy10
            * (-(sy1 * sy11) - sy12 + sy5 * sy6 + sy11 * sy8
                - 2. * G({ 0, 0, c / d, b / d }, 1) - G({ 0, c / d, 0, b / d }, 1)
                - sy4 * G({ 0, b }, { 1, sb }, x)
                - G({ 0, 0, c, b }, { 1, 1, sc, sb }, x))
        + G({ 0, 0, c, b, a }, { 1, 1, sc, sb, sa }, x)
        + (sy12 + sy3 * sy4 - sy10 * sy6) * Log(-x, sd) };
    if (b != x) {
        res += (sy10 * sy6 - sy10 * sy7) * G({ b }, { sb }, x)
            + (-sy6 + sy7) * G({ b, a }, { sb, sa }, x);
    }
    if (c != x) {
        res += (sy13 * sy3 - sy3 * sy4) * G({ c }, { sc }, x)
            + (-(sy10 * sy13) + sy10 * sy4) * G({ c, b }, { sc, sb }, x)
            + (sy13 - sy4) * G({ c, b, a }, { sc, sb, sa }, x);
    }
    return res;
}
complex<double> G5_abc0d_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    if (c == d) {

        const complex<double> sy1 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy2 = G({ a, b }, { sa, sb }, x);
        const complex<double> sy3 = G({ 0, 1, b / c }, 1);
        const complex<double> sy4 = G({ a, c }, { sa, sc }, x);
        const complex<double> sy5 = G({ 0, b / c, 1 }, 1);
        const complex<double> sy6 = G({ a, 0, c }, { sa, 1, sc }, x);
        const complex<double> sy7 = G({ a, b, c }, { sa, sb, sc }, x);
        const complex<double> sy8 = G({ b / c, 0, 1 }, 1);
        const complex<double> sy9 = G({ 0, b / c, a / c }, 1);
        const complex<double> sy10 = G({ b / c, 0, a / c }, 1);
        const complex<double> sy11 = G({ c }, { sc }, x);
        const complex<double> sy12 = G({ 0, b / c, a / c, 1 }, 1);
        const complex<double> sy13 = G({ a }, { sa }, x);
        const complex<double> sy14 = G({ 0, 1, b / c, a / c }, 1);
        const complex<double> sy15 = G({ 0, b / c, 1, a / c }, 1);
        const complex<double> sy16 = G({ b / c, 0, 1, a / c }, 1);
        const complex<double> sy17 = G({ b / c, 0, a / c, 1 }, 1);
        const complex<double> sy18 = G({ b / c, a / c, 0, 1 }, 1);

        complex<double> tmp { -(sy11 * sy12) + (sy1 + sy2) * sy3 + sy4 * sy5
            + sy1 * (sy5 + sy8) + sy4 * (-sy10 + sy8 - sy9)
            + (-sy6 + sy7) * G({ 0, b / c }, 1)
            + sy13
                * (-sy14 - sy15 - sy16 - 2. * G({ 0, 0, 1, b / c }, 1)
                    - 2. * G({ 0, 0, b / c, 1 }, 1) - G({ 0, 1, 0, b / c }, 1)
                    - 2. * G({ 0, b / c, 0, 1 }, 1) - 2. * G({ b / c, 0, 0, 1 }, 1))
            + sy11
                * (-sy17 - sy18 + G({ 0, b / c, a / c, x / c }, 1)
                    + G({ b / c, 0, a / c, x / c }, 1) + G({ b / c, a / c, 0, x / c }, 1))
            + 2. * G({ 0, 0, 1, b / c, a / c }, 1) + 2. * G({ 0, 0, b / c, 1, a / c }, 1)
            + 2. * G({ 0, 0, b / c, a / c, 1 }, 1) + G({ 0, 1, 0, b / c, a / c }, 1)
            + G({ 0, 1, b / c, 0, a / c }, 1) + G({ 0, 1, b / c, a / c, x / c }, 1)
            + 2. * G({ 0, b / c, 0, 1, a / c }, 1) + 2. * G({ 0, b / c, 0, a / c, 1 }, 1)
            + G({ 0, b / c, 1, 0, a / c }, 1) + G({ 0, b / c, 1, a / c, x / c }, 1)
            + 2. * G({ 0, b / c, a / c, 0, 1 }, 1) + G({ 0, b / c, a / c, 1, x / c }, 1)
            + 2. * G({ b / c, 0, 0, 1, a / c }, 1) + 2. * G({ b / c, 0, 0, a / c, 1 }, 1)
            + G({ b / c, 0, 1, 0, a / c }, 1) + G({ b / c, 0, 1, a / c, x / c }, 1)
            + 2. * G({ b / c, 0, a / c, 0, 1 }, 1) + G({ b / c, 0, a / c, 1, x / c }, 1)
            + 2. * G({ b / c, a / c, 0, 0, 1 }, 1) + G({ b / c, a / c, 0, 1, x / c }, 1)
            + (sy10 + sy9 + G({ b / c, a / c, x / c }, 1)) * G({ 0, c }, { 1, sc }, x)
            + G({ b / c, a / c }, 1) * (-sy6 + G({ 0, 0, c }, { 1, 1, sc }, x))
            + G({ b / c }, 1)
                * (-G({ a, 0, 0, c }, { sa, 1, 1, sc }, x)
                    + G({ a, b, 0, c }, { sa, sb, 1, sc }, x))
            + G({ a, b, 0, 0, c }, { sa, sb, 1, 1, sc }, x) + sy7 * Zeta(2)
            + G({ 0, a, b }, { 1, sa, sb }, x) * Zeta(2)
            + G({ a, 0, b }, { sa, 1, sb }, x) * Zeta(2)
            + Log(-x, sc)
                * (sy12 + sy14 + sy15 + sy16 + sy17 + sy18 + sy13 * (-sy3 - sy5 - sy8)
                    - sy2 * Zeta(2))
            + Log(c, sc)
                * (-sy12 - sy14 - sy15 - sy16 - sy17 - sy18 + sy13 * (sy3 + sy5 + sy8)
                    + sy2 * Zeta(2))
            - 2. * sy2 * Zeta(3) };
        return tmp;
    } else {

        const complex<double> sy1 = G({ a, d }, { sa, sd }, x);
        const complex<double> sy2 = G({ 0, b / c, a / c }, 1);
        const complex<double> sy3 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy4 = G({ 0, b / c, d / c }, 1);
        const complex<double> sy5 = G({ a, b }, { sa, sb }, x);
        const complex<double> sy6 = G({ 0, d / c, b / c }, 1);
        const complex<double> sy7 = G({ a, 0, d }, { sa, 1, sd }, x);
        const complex<double> sy8 = G({ 0, d / c }, 1);
        const complex<double> sy9 = G({ a, b, d }, { sa, sb, sd }, x);
        const complex<double> sy10 = G({ b / c, 0, d / c }, 1);
        const complex<double> sy11 = G({ b / c, 0, a / c }, 1);
        const complex<double> sy12 = G({ a }, { sa }, x);
        const complex<double> sy13 = G({ 0, b / c, d / c, a / c }, 1);
        const complex<double> sy14 = G({ 0, d / c, b / c, a / c }, 1);
        const complex<double> sy15 = G({ b / c, 0, d / c, a / c }, 1);
        const complex<double> sy16 = G({ 0, b / c, a / c, d / c }, 1);
        const complex<double> sy17 = G({ b / c, 0, a / c, d / c }, 1);
        const complex<double> sy18 = G({ b / c, a / c, 0, d / c }, 1);
        complex<double> res { -(sy1 * sy2) + sy3 * sy4 + sy1 * (sy10 - sy11 + sy4)
            + sy3 * (sy10 + sy6) + (-sy7 + sy9) * G({ 0, b / c }, 1)
            + sy5 * (sy6 + 2. * G({ 0, 0, d / c }, 1))
            + sy12
                * (-sy13 - sy14 - sy15 - 2. * G({ 0, 0, b / c, d / c }, 1)
                    - 2. * G({ 0, 0, d / c, b / c }, 1)
                    - 2. * G({ 0, b / c, 0, d / c }, 1) - G({ 0, d / c, 0, b / c }, 1)
                    - 2. * G({ b / c, 0, 0, d / c }, 1))
            + 2. * G({ 0, 0, b / c, a / c, d / c }, 1)
            + 2. * G({ 0, 0, b / c, d / c, a / c }, 1)
            + 2. * G({ 0, 0, d / c, b / c, a / c }, 1)
            + 2. * G({ 0, b / c, 0, a / c, d / c }, 1)
            + 2. * G({ 0, b / c, 0, d / c, a / c }, 1)
            + 2. * G({ 0, b / c, a / c, 0, d / c }, 1)
            + G({ 0, b / c, a / c, d / c, x / c }, 1)
            + G({ 0, b / c, d / c, 0, a / c }, 1)
            + G({ 0, b / c, d / c, a / c, x / c }, 1)
            + G({ 0, d / c, 0, b / c, a / c }, 1) + G({ 0, d / c, b / c, 0, a / c }, 1)
            + G({ 0, d / c, b / c, a / c, x / c }, 1)
            + 2. * G({ b / c, 0, 0, a / c, d / c }, 1)
            + 2. * G({ b / c, 0, 0, d / c, a / c }, 1)
            + 2. * G({ b / c, 0, a / c, 0, d / c }, 1)
            + G({ b / c, 0, a / c, d / c, x / c }, 1)
            + G({ b / c, 0, d / c, 0, a / c }, 1)
            + G({ b / c, 0, d / c, a / c, x / c }, 1)
            + 2. * G({ b / c, a / c, 0, 0, d / c }, 1)
            + G({ b / c, a / c, 0, d / c, x / c }, 1)
            + (sy11 + sy2 + G({ b / c, a / c, x / c }, 1)) * G({ 0, d }, { 1, sd }, x)
            + G({ b / c, a / c }, 1) * (-sy7 + G({ 0, 0, d }, { 1, 1, sd }, x))
            + sy8
                * (-sy9 - G({ 0, a, b }, { 1, sa, sb }, x)
                    - G({ a, 0, b }, { sa, 1, sb }, x))
            + G({ b / c }, 1)
                * (-G({ a, 0, 0, d }, { sa, 1, 1, sd }, x)
                    + G({ a, b, 0, d }, { sa, sb, 1, sd }, x))
            + G({ a, b, 0, 0, d }, { sa, sb, 1, 1, sd }, x)
            + (-sy13 - sy14 - sy15 - sy16 - sy17 - sy18 + sy12 * (sy10 + sy4 + sy6)
                  - sy5 * sy8)
                * Log(c, sc)
            + (sy13 + sy14 + sy15 + sy16 + sy17 + sy18 + sy12 * (-sy10 - sy4 - sy6)
                  + sy5 * sy8)
                * Log(-x, sc) };
        if (d != x) {
            res += (-sy16 - sy17 - sy18 + G({ 0, b / c, a / c, x / c }, 1)
                       + G({ b / c, 0, a / c, x / c }, 1)
                       + G({ b / c, a / c, 0, x / c }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G5_abc0d_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{

    if (b == c) {

        const complex<double> sy1 = G({ a, d }, { sa, sd }, x);
        const complex<double> sy2 = G({ 0, 1, a / b }, 1);
        const complex<double> sy3 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy4 = G({ 0, 1, d / b }, 1);
        const complex<double> sy5 = G({ 0, d / b, 1 }, 1);
        const complex<double> sy6 = G({ 0, a / b, 1 }, 1);
        const complex<double> sy7 = G({ b, 0, d }, { sb, 1, sd }, x);
        const complex<double> sy8 = G({ a }, { sa }, x);
        const complex<double> sy9 = G({ 0, 1, d / b, a / b }, 1);
        const complex<double> sy10 = G({ 0, d / b, 1, a / b }, 1);
        const complex<double> sy11 = G({ 0, d / b, a / b, 1 }, 1);
        const complex<double> sy12 = G({ 0, 1, a / b, d / b }, 1);
        const complex<double> sy13 = G({ 0, a / b, 1, d / b }, 1);
        const complex<double> sy14 = G({ 0, a / b, d / b, 1 }, 1);
        const complex<double> sy15 = G({ a / b, 0, 1, d / b }, 1);
        const complex<double> sy16 = G({ a / b, 0, d / b, 1 }, 1);
        complex<double> res { -(sy1 * sy2) + sy3 * sy4 + sy3 * sy5
            + sy1 * (sy4 + sy5 - sy6) - sy7 * G({ a / b, x / b }, 1)
            + sy8
                * (-sy10 - sy11 - sy9 - 2. * G({ 0, 0, 1, d / b }, 1)
                    - 2. * G({ 0, 0, d / b, 1 }, 1) - G({ 0, 1, 0, d / b }, 1))
            + 2. * G({ 0, 0, 1, a / b, d / b }, 1) + 2. * G({ 0, 0, 1, d / b, a / b }, 1)
            + 2. * G({ 0, 0, a / b, 1, d / b }, 1) + 2. * G({ 0, 0, a / b, d / b, 1 }, 1)
            + 2. * G({ 0, 0, d / b, 1, a / b }, 1) + 2. * G({ 0, 0, d / b, a / b, 1 }, 1)
            + G({ 0, 1, 0, a / b, d / b }, 1) + G({ 0, 1, 0, d / b, a / b }, 1)
            + G({ 0, 1, a / b, 0, d / b }, 1) + G({ 0, 1, a / b, d / b, x / b }, 1)
            + G({ 0, 1, d / b, 0, a / b }, 1) + G({ 0, 1, d / b, a / b, x / b }, 1)
            + 2. * G({ 0, a / b, 0, 1, d / b }, 1) + 2. * G({ 0, a / b, 0, d / b, 1 }, 1)
            + G({ 0, a / b, 1, 0, d / b }, 1) + G({ 0, a / b, 1, d / b, x / b }, 1)
            + G({ 0, a / b, d / b, 1, x / b }, 1) + G({ 0, a / b, d / b, x / b, 1 }, 1)
            + G({ 0, d / b, 0, 1, a / b }, 1) + G({ 0, d / b, 0, a / b, 1 }, 1)
            + G({ 0, d / b, 1, 0, a / b }, 1) + G({ 0, d / b, 1, a / b, x / b }, 1)
            + G({ 0, d / b, a / b, 1, x / b }, 1) + G({ 0, d / b, a / b, x / b, 1 }, 1)
            + 2. * G({ a / b, 0, 0, 1, d / b }, 1) + 2. * G({ a / b, 0, 0, d / b, 1 }, 1)
            + G({ a / b, 0, 1, 0, d / b }, 1) + G({ a / b, 0, 1, d / b, x / b }, 1)
            + G({ a / b, 0, d / b, 1, x / b }, 1) + G({ a / b, 0, d / b, x / b, 1 }, 1)
            + (sy2 + sy6 + G({ a / b, x / b, 1 }, 1)) * G({ 0, d }, { 1, sd }, x)
            + G({ a / b, 1 }, 1) * (sy7 - G({ a, 0, d }, { sa, 1, sd }, x))
            + G({ a / b }, 1)
                * (-G({ 0, b, 0, d }, { 1, sb, 1, sd }, x)
                    + G({ a, b, 0, d }, { sa, sb, 1, sd }, x))
            + G({ a, 0, b, 0, d }, { sa, 1, sb, 1, sd }, x)
            + (-sy10 - sy11 - sy12 - sy13 - sy14 - sy15 - sy16 + (sy4 + sy5) * sy8 - sy9)
                * Log(b, sb)
            + (sy10 + sy11 + sy12 + sy13 + sy14 + sy15 + sy16 + (-sy4 - sy5) * sy8 + sy9)
                * Log(-x, sb) };
        if (d != x) {
            res += (-sy12 - sy13 - sy14 - sy15 - sy16 + G({ 0, 1, a / b, x / b }, 1)
                       + G({ 0, a / b, 1, x / b }, 1) + G({ 0, a / b, x / b, 1 }, 1)
                       + G({ a / b, 0, 1, x / b }, 1) + G({ a / b, 0, x / b, 1 }, 1))
                * G({ d }, { sd }, x);
        }
        return res;

    } else {

        const complex<double> sy1 = G({ a / b, c / b }, 1);
        const complex<double> sy2 = G({ 0, 0, d }, { 1, 1, sd }, x);
        const complex<double> sy3 = G({ c / b, a / b }, 1);
        const complex<double> sy4 = G({ a, d }, { sa, sd }, x);
        const complex<double> sy5 = G({ c / b, 0, a / b }, 1);
        const complex<double> sy6 = G({ c / b, 0, d / b }, 1);
        const complex<double> sy7 = G({ a, c, 0, d }, { sa, sc, 1, sd }, x);
        const complex<double> sy8 = G({ a }, { sa }, x);
        const complex<double> sy9 = G({ c / b, 0, d / b, a / b }, 1);
        const complex<double> sy10 = G({ a / b, c / b, 0, d / b }, 1);
        const complex<double> sy11 = G({ c / b, 0, a / b, d / b }, 1);
        const complex<double> sy12 = G({ c / b, a / b, 0, d / b }, 1);
        complex<double> res { -(sy1 * sy2) - sy2 * sy3 + sy4 * sy5 - sy4 * sy6
            + sy8
                * (sy9 + G({ 0, c / b, 0, d / b }, 1) + 2. * G({ c / b, 0, 0, d / b }, 1))
            - G({ 0, a / b, c / b, 0, d / b }, 1) - G({ 0, c / b, 0, a / b, d / b }, 1)
            - G({ 0, c / b, 0, d / b, a / b }, 1) - G({ 0, c / b, a / b, 0, d / b }, 1)
            - G({ a / b, 0, c / b, 0, d / b }, 1)
            - 2. * G({ a / b, c / b, 0, 0, d / b }, 1)
            - G({ a / b, c / b, 0, d / b, x / b }, 1)
            - 2. * G({ c / b, 0, 0, a / b, d / b }, 1)
            - 2. * G({ c / b, 0, 0, d / b, a / b }, 1)
            - 2. * G({ c / b, 0, a / b, 0, d / b }, 1)
            - G({ c / b, 0, a / b, d / b, x / b }, 1)
            - G({ c / b, 0, d / b, 0, a / b }, 1)
            - G({ c / b, 0, d / b, a / b, x / b }, 1)
            - 2. * G({ c / b, a / b, 0, 0, d / b }, 1)
            - G({ c / b, a / b, 0, d / b, x / b }, 1) - sy6 * G({ 0, a }, { 1, sa }, x)
            + (-sy5 - G({ a / b, c / b, x / b }, 1) - G({ c / b, a / b, x / b }, 1))
                * G({ 0, d }, { 1, sd }, x)
            + sy3 * G({ a, 0, d }, { sa, 1, sd }, x)
            + G({ a / b }, 1) * (sy7 - G({ 0, c, 0, d }, { 1, sc, 1, sd }, x))
            + G({ c / b }, 1) * (-sy7 + G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
            + G({ a, 0, c, 0, d }, { sa, 1, sc, 1, sd }, x)
            + (sy10 + sy11 + sy12 - sy6 * sy8 + sy9) * Log(b, sb)
            + (-sy10 - sy11 - sy12 + sy6 * sy8 - sy9) * Log(-x, sb) };
        if (c != x) {
            res += (sy1 - G({ a / b, x / b }, 1)) * G({ c, 0, d }, { sc, 1, sd }, x);
        }
        if (d != x) {
            res += (sy10 + sy11 + sy12 - G({ a / b, c / b, 0, x / b }, 1)
                       - G({ c / b, 0, a / b, x / b }, 1)
                       - G({ c / b, a / b, 0, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G5_abc0d_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{

    if (a == b) {
        if (a == c) {
            const complex<double> sy1 = G({ 0, 1, 1, d / a }, 1);
            const complex<double> sy2 = G({ 0, 1, d / a, 1 }, 1);
            const complex<double> sy3 = G({ 0, d / a, 1, 1 }, 1);
            complex<double> res { 2. * G({ 0, 0, 1, 1, d / a }, 1)
                + 2. * G({ 0, 0, 1, d / a, 1 }, 1) + 2. * G({ 0, 0, d / a, 1, 1 }, 1)
                + G({ 0, 1, 0, 1, d / a }, 1) + G({ 0, 1, 0, d / a, 1 }, 1)
                + G({ 0, 1, 1, 0, d / a }, 1) + G({ 0, 1, 1, d / a, x / a }, 1)
                + G({ 0, 1, d / a, 1, x / a }, 1) + G({ 0, 1, d / a, x / a, 1 }, 1)
                + G({ 0, d / a, 1, 1, x / a }, 1) + G({ 0, d / a, 1, x / a, 1 }, 1)
                + G({ 0, d / a, x / a, 1, 1 }, 1)
                + G({ x / a, 1, 1 }, 1) * G({ 0, d }, { 1, sd }, x)
                - G({ x / a, 1 }, 1) * G({ a, 0, d }, { sa, 1, sd }, x)
                + G({ x / a }, 1) * G({ a, a, 0, d }, { sa, sa, 1, sd }, x)
                + G({ 0, a, a, 0, d }, { 1, sa, sa, 1, sd }, x)
                + (-sy1 - sy2 - sy3) * Log(a, sa) + (sy1 + sy2 + sy3) * Log(-x, sa) };
            if (d != x) {
                res += (-sy1 - sy2 - sy3 + G({ 0, 1, 1, x / a }, 1)
                           + G({ 0, 1, x / a, 1 }, 1) + G({ 0, x / a, 1, 1 }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        } else {

            const complex<double> sy1 = G({ c / a, 1 }, 1);
            const complex<double> sy2 = G({ c / a, 0, 1, d / a }, 1);
            const complex<double> sy3 = G({ c / a, 0, d / a, 1 }, 1);
            const complex<double> sy4 = G({ c / a, 1, 0, d / a }, 1);
            complex<double> res { -G({ 0, c / a, 0, 1, d / a }, 1)
                - G({ 0, c / a, 0, d / a, 1 }, 1) - G({ 0, c / a, 1, 0, d / a }, 1)
                - 2. * G({ c / a, 0, 0, 1, d / a }, 1)
                - 2. * G({ c / a, 0, 0, d / a, 1 }, 1)
                - 2. * G({ c / a, 0, 1, 0, d / a }, 1)
                - G({ c / a, 0, 1, d / a, x / a }, 1)
                - G({ c / a, 0, d / a, 1, x / a }, 1)
                - G({ c / a, 0, d / a, x / a, 1 }, 1)
                - 2. * G({ c / a, 1, 0, 0, d / a }, 1)
                - G({ c / a, 1, 0, d / a, x / a }, 1)
                + (-G({ c / a, 1, x / a }, 1) - G({ c / a, x / a, 1 }, 1))
                    * G({ 0, d }, { 1, sd }, x)
                - sy1 * G({ 0, 0, d }, { 1, 1, sd }, x)
                + G({ x / a }, 1) * G({ a, c, 0, d }, { sa, sc, 1, sd }, x)
                + G({ 0, a, c, 0, d }, { 1, sa, sc, 1, sd }, x)
                + (sy2 + sy3 + sy4) * Log(a, sa) + (-sy2 - sy3 - sy4) * Log(-x, sa) };
            if (c != x) {
                res += (sy1 - G({ x / a, 1 }, 1)) * G({ c, 0, d }, { sc, 1, sd }, x);
            }
            if (d != x) {
                res += (sy2 + sy3 + sy4 - G({ c / a, 0, 1, x / a }, 1)
                           - G({ c / a, 0, x / a, 1 }, 1) - G({ c / a, 1, 0, x / a }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        }
    } else {

        const complex<double> sy1 = G({ b / a, c / a }, 1);
        const complex<double> sy2 = G({ b / a }, 1);
        const complex<double> sy3 = G({ b / a, c / a, 0, d / a }, 1);
        complex<double> res { G({ 0, b / a, c / a, 0, d / a }, 1)
            + G({ b / a, 0, c / a, 0, d / a }, 1)
            + 2. * G({ b / a, c / a, 0, 0, d / a }, 1)
            + G({ b / a, c / a, 0, d / a, x / a }, 1)
            + G({ b / a, c / a, x / a }, 1) * G({ 0, d }, { 1, sd }, x)
            + sy1 * G({ 0, 0, d }, { 1, 1, sd }, x)
            + sy2 * G({ 0, c, 0, d }, { 1, sc, 1, sd }, x)
            + G({ 0, b, c, 0, d }, { 1, sb, sc, 1, sd }, x) - sy3 * Log(a, sa)
            + sy3 * Log(-x, sa) };
        if (b != x) {
            res += (-sy2 + G({ x / a }, 1)) * G({ b, c, 0, d }, { sb, sc, 1, sd }, x);
        }
        if (c != x) {
            res += (-sy1 + G({ b / a, x / a }, 1)) * G({ c, 0, d }, { sc, 1, sd }, x);
        }
        if (d != x) {
            res += (-sy3 + G({ b / a, c / a, 0, x / a }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G5_abcde_e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{
    const complex<double> sy1 = G({ d / e, c / e }, 1);
    const complex<double> sy2 = G({ a, b, c }, { sa, sb, sc }, x);
    const complex<double> sy3 = G({ d / e, c / e, b / e }, 1);
    const complex<double> sy4 = G({ a, b }, { sa, sb }, x);
    const complex<double> sy5 = G({ a, b, c, d }, { sa, sb, sc, sd }, x);
    const complex<double> sy6 = G({ d / e }, 1);
    const complex<double> sy7 = G({ a }, { sa }, x);
    const complex<double> sy8 = G({ d / e, c / e, b / e, a / e }, 1);
    complex<double> res { -(sy2 * G({ 0, d / e }, 1))
        + sy4 * (sy3 + G({ 0, d / e, c / e }, 1) + G({ d / e, 0, c / e }, 1))
        + sy7
            * (-sy8 - G({ 0, d / e, c / e, b / e }, 1) - G({ d / e, 0, c / e, b / e }, 1)
                - G({ d / e, c / e, 0, b / e }, 1))
        + G({ 0, d / e, c / e, b / e, a / e }, 1)
        + G({ d / e, 0, c / e, b / e, a / e }, 1)
        + G({ d / e, c / e, 0, b / e, a / e }, 1)
        + G({ d / e, c / e, b / e, 0, a / e }, 1)
        + G({ d / e, c / e, b / e, a / e, x / e }, 1)
        + sy5 * (-G({ x / e }, 1) + G({ e }, { se }, x)) + sy3 * G({ 0, a }, { 1, sa }, x)
        + sy1
            * (-sy2 - G({ 0, a, b }, { 1, sa, sb }, x) - G({ a, 0, b }, { sa, 1, sb }, x))
        + sy6
            * (sy5 + G({ 0, a, b, c }, { 1, sa, sb, sc }, x)
                + G({ a, 0, b, c }, { sa, 1, sb, sc }, x)
                + G({ a, b, 0, c }, { sa, sb, 1, sc }, x))
        - G({ 0, a, b, c, d }, { 1, sa, sb, sc, sd }, x)
        - G({ a, 0, b, c, d }, { sa, 1, sb, sc, sd }, x)
        - G({ a, b, 0, c, d }, { sa, sb, 1, sc, sd }, x)
        - G({ a, b, c, 0, d }, { sa, sb, sc, 1, sd }, x)
        + (-(sy1 * sy4) + sy2 * sy6 + sy3 * sy7 - sy8) * Log(e, se)
        + (sy1 * sy4 - sy2 * sy6 - sy3 * sy7 + sy8) * Log(-x, se) };

    return res;
}
complex<double> G5_abcde_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{
    if (d == e) {

        const complex<double> sy1 = G({ a, b, d }, { sa, sb, sd }, x);
        const complex<double> sy2 = G({ a, b, c }, { sa, sb, sc }, x);
        const complex<double> sy3 = G({ a, d }, { sa, sd }, x);
        const complex<double> sy4 = G({ c / d, b / d, 1 }, 1);
        const complex<double> sy5 = G({ c / d, b / d, a / d }, 1);
        const complex<double> sy6 = G({ c / d, b / d, a / d, 1 }, 1);
        const complex<double> sy7 = G({ d }, { sd }, x);
        complex<double> tmp { -(sy3 * sy4) + sy3 * sy5 + sy6 * sy7
            + (sy1 - sy2) * G({ c / d, 1 }, 1)
            - sy7 * G({ c / d, b / d, a / d, x / d }, 1)
            - G({ c / d, b / d, a / d, 0, 1 }, 1)
            + G({ c / d, b / d, a / d, x / d, 1 }, 1)
            + (-sy6 + G({ c / d, b / d, 0, 1 }, 1)) * G({ a }, { sa }, x)
            - sy5 * G({ 0, d }, { 1, sd }, x)
            + (sy4 - G({ c / d, 0, 1 }, 1)) * G({ a, b }, { sa, sb }, x)
            + G({ c / d, b / d }, 1) * (-sy1 + G({ a, 0, d }, { sa, 1, sd }, x))
            + G({ c / d }, 1)
                * (-G({ a, b, 0, d }, { sa, sb, 1, sd }, x)
                    + G({ a, b, c, d }, { sa, sb, sc, sd }, x))
            + G({ a, b, c, 0, d }, { sa, sb, sc, 1, sd }, x) - sy2 * Zeta(2) };
        return tmp;
    } else {

        const complex<double> sy1 = G({ a, b, c }, { sa, sb, sc }, x);
        const complex<double> sy2 = G({ e / d, c / d }, 1);
        const complex<double> sy3 = G({ 0, a, b }, { 1, sa, sb }, x);
        const complex<double> sy4 = G({ a, 0, b }, { sa, 1, sb }, x);
        const complex<double> sy5 = G({ a, b, e }, { sa, sb, se }, x);
        const complex<double> sy6 = G({ c / d, e / d }, 1);
        const complex<double> sy7 = G({ c / d, b / d, a / d }, 1);
        const complex<double> sy8 = G({ a, e }, { sa, se }, x);
        const complex<double> sy9 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy10 = G({ c / d, b / d, e / d }, 1);
        const complex<double> sy11 = G({ c / d, e / d, b / d }, 1);
        const complex<double> sy12 = G({ e / d, c / d, b / d }, 1);
        const complex<double> sy13 = G({ a, b }, { sa, sb }, x);
        const complex<double> sy14 = G({ e / d }, 1);
        const complex<double> sy15 = G({ a, b, c, e }, { sa, sb, sc, se }, x);
        const complex<double> sy16 = G({ a }, { sa }, x);
        const complex<double> sy17 = G({ c / d, b / d, a / d, e / d }, 1);
        const complex<double> sy18 = G({ c / d, b / d, e / d, a / d }, 1);
        const complex<double> sy19 = G({ c / d, e / d, b / d, a / d }, 1);
        const complex<double> sy20 = G({ e / d, c / d, b / d, a / d }, 1);

        complex<double> res { sy2 * (sy1 + sy3 + sy4) + (sy3 + sy4 + sy5) * sy6
            - sy10 * sy8 + sy7 * sy8 - sy10 * sy9 + (-sy11 - sy12) * sy9
            + sy1 * G({ 0, e / d }, 1)
            + sy13
                * (-sy11 - sy12 - G({ 0, c / d, e / d }, 1) - G({ 0, e / d, c / d }, 1)
                    - G({ c / d, 0, e / d }, 1) - G({ e / d, 0, c / d }, 1))
            + sy16
                * (sy18 + sy19 + sy20 + G({ 0, c / d, b / d, e / d }, 1)
                    + G({ 0, c / d, e / d, b / d }, 1) + G({ 0, e / d, c / d, b / d }, 1)
                    + G({ c / d, 0, b / d, e / d }, 1) + G({ c / d, 0, e / d, b / d }, 1)
                    + G({ c / d, b / d, 0, e / d }, 1) + G({ c / d, e / d, 0, b / d }, 1)
                    + G({ e / d, 0, c / d, b / d }, 1) + G({ e / d, c / d, 0, b / d }, 1))
            - G({ 0, c / d, b / d, a / d, e / d }, 1)
            - G({ 0, c / d, b / d, e / d, a / d }, 1)
            - G({ 0, c / d, e / d, b / d, a / d }, 1)
            - G({ 0, e / d, c / d, b / d, a / d }, 1)
            - G({ c / d, 0, b / d, a / d, e / d }, 1)
            - G({ c / d, 0, b / d, e / d, a / d }, 1)
            - G({ c / d, 0, e / d, b / d, a / d }, 1)
            - G({ c / d, b / d, 0, a / d, e / d }, 1)
            - G({ c / d, b / d, 0, e / d, a / d }, 1)
            - G({ c / d, b / d, a / d, 0, e / d }, 1)
            - G({ c / d, b / d, a / d, e / d, x / d }, 1)
            - G({ c / d, b / d, e / d, 0, a / d }, 1)
            - G({ c / d, b / d, e / d, a / d, x / d }, 1)
            - G({ c / d, e / d, 0, b / d, a / d }, 1)
            - G({ c / d, e / d, b / d, 0, a / d }, 1)
            - G({ c / d, e / d, b / d, a / d, x / d }, 1)
            - G({ e / d, 0, c / d, b / d, a / d }, 1)
            - G({ e / d, c / d, 0, b / d, a / d }, 1)
            - G({ e / d, c / d, b / d, 0, a / d }, 1)
            - G({ e / d, c / d, b / d, a / d, x / d }, 1)
            - sy7 * G({ 0, e }, { 1, se }, x)
            + G({ c / d, b / d }, 1) * (-sy5 + G({ a, 0, e }, { sa, 1, se }, x))
            + sy14
                * (-sy15 - G({ 0, a, b, c }, { 1, sa, sb, sc }, x)
                    - G({ a, 0, b, c }, { sa, 1, sb, sc }, x)
                    - G({ a, b, 0, c }, { sa, sb, 1, sc }, x))
            + G({ c / d }, 1) * (sy15 - G({ a, b, 0, e }, { sa, sb, 1, se }, x))
            + G({ a, b, c, 0, e }, { sa, sb, sc, 1, se }, x)
            + (-(sy1 * sy14) + (-sy10 - sy11 - sy12) * sy16 + sy17 + sy18 + sy19 + sy20
                  + sy13 * (sy2 + sy6))
                * Log(d, sd)
            + (sy1 * sy14 + (sy10 + sy11 + sy12) * sy16 - sy17 - sy18 - sy19 - sy20
                  + sy13 * (-sy2 - sy6))
                * Log(-x, sd) };
        if (e != x) {
            res += (sy17 - G({ c / d, b / d, a / d, x / d }, 1)) * G({ e }, { se }, x);
        }
        return res;
    }
}

// complex<double> G5_abcde_c(complex<double> a, complex<double> b, complex<double> c,
//     complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
//     double x)
// {

//     return 1.0;
// }

complex<double> G5_abcde_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{

    if (c == d) {
        if (c == e) {
            const complex<double> sy1 = G({ a, c, c }, { sa, sc, sc }, x);
            const complex<double> sy2 = G({ a, c }, { sa, sc }, x);
            const complex<double> sy3 = G({ b / c, 1, 1 }, 1);
            const complex<double> sy4 = G({ b / c, a / c, 1 }, 1);
            const complex<double> sy5 = G({ c, c }, { sc, sc }, x);
            const complex<double> sy6 = G({ b / c, a / c, 1, 1 }, 1);
            const complex<double> sy7 = G({ c }, { sc }, x);
            complex<double> tmp { -(sy2 * sy3) + sy2 * sy4 - sy4 * sy5 + sy6 * sy7
                + sy5 * G({ b / c, a / c, x / c }, 1)
                - sy7 * G({ b / c, a / c, x / c, 1 }, 1) - G({ b / c, a / c, 0, 1, 1 }, 1)
                + G({ b / c, a / c, x / c, 1, 1 }, 1)
                + (-sy6 + G({ b / c, 0, 1, 1 }, 1)) * G({ a }, { sa }, x)
                + G({ b / c, a / c }, 1) * (-sy1 + G({ 0, c, c }, { 1, sc, sc }, x))
                + G({ b / c, 1 }, 1) * (sy1 - G({ a, b, c }, { sa, sb, sc }, x))
                + G({ b / c }, 1)
                    * (-G({ a, 0, c, c }, { sa, 1, sc, sc }, x)
                        + G({ a, b, c, c }, { sa, sb, sc, sc }, x))
                + G({ a, b, 0, c, c }, { sa, sb, 1, sc, sc }, x)
                + G({ a, b }, { sa, sb }, x) * (sy3 - Zeta(3)) };
            return tmp;
        } else {

            const complex<double> sy1 = G({ e / c, 1 }, 1);
            const complex<double> sy2 = G({ a, b, e }, { sa, sb, se }, x);
            const complex<double> sy3 = G({ a, c, e }, { sa, sc, se }, x);
            const complex<double> sy4 = G({ c, e }, { sc, se }, x);
            const complex<double> sy5 = G({ b / c, a / c, 1 }, 1);
            const complex<double> sy6 = G({ b / c, e / c, 1 }, 1);
            const complex<double> sy7 = G({ 0, a }, { 1, sa }, x);
            const complex<double> sy8 = G({ e / c, 1, b / c }, 1);
            const complex<double> sy9 = G({ e / c, b / c, 1 }, 1);
            const complex<double> sy10 = G({ a, b }, { sa, sb }, x);
            const complex<double> sy11 = G({ a }, { sa }, x);
            const complex<double> sy12 = G({ b / c, a / c, e / c, 1 }, 1);
            const complex<double> sy13 = G({ b / c, e / c, 1, a / c }, 1);
            const complex<double> sy14 = G({ b / c, e / c, a / c, 1 }, 1);
            const complex<double> sy15 = G({ e / c, 1, b / c, a / c }, 1);
            const complex<double> sy16 = G({ e / c, b / c, 1, a / c }, 1);
            const complex<double> sy17 = G({ e / c, b / c, a / c, 1 }, 1);
            complex<double> res { -(sy4 * sy5) - sy6 * sy7 + sy7 * (-sy8 - sy9)
                + (-sy2 + sy3) * G({ b / c, 1 }, 1)
                + sy10 * (-sy8 - sy9 - G({ 0, e / c, 1 }, 1))
                + sy4 * G({ b / c, a / c, x / c }, 1)
                + sy11
                    * (sy13 + sy14 + sy15 + sy16 + sy17 + G({ 0, b / c, e / c, 1 }, 1)
                        + G({ 0, e / c, 1, b / c }, 1) + G({ 0, e / c, b / c, 1 }, 1)
                        + G({ b / c, 0, e / c, 1 }, 1) + G({ e / c, 0, 1, b / c }, 1)
                        + G({ e / c, 0, b / c, 1 }, 1) + G({ e / c, 1, 0, b / c }, 1))
                - G({ 0, b / c, a / c, e / c, 1 }, 1)
                - G({ 0, b / c, e / c, 1, a / c }, 1)
                - G({ 0, b / c, e / c, a / c, 1 }, 1)
                - G({ 0, e / c, 1, b / c, a / c }, 1)
                - G({ 0, e / c, b / c, 1, a / c }, 1)
                - G({ 0, e / c, b / c, a / c, 1 }, 1)
                - G({ b / c, 0, a / c, e / c, 1 }, 1)
                - G({ b / c, 0, e / c, 1, a / c }, 1)
                - G({ b / c, 0, e / c, a / c, 1 }, 1)
                - G({ b / c, a / c, 0, e / c, 1 }, 1)
                - G({ b / c, a / c, e / c, 1, x / c }, 1)
                - G({ b / c, a / c, e / c, x / c, 1 }, 1)
                - G({ b / c, e / c, 0, 1, a / c }, 1)
                - G({ b / c, e / c, 0, a / c, 1 }, 1)
                - G({ b / c, e / c, 1, 0, a / c }, 1)
                - G({ b / c, e / c, 1, a / c, x / c }, 1)
                - G({ b / c, e / c, a / c, 1, x / c }, 1)
                - G({ b / c, e / c, a / c, x / c, 1 }, 1)
                - G({ e / c, 0, 1, b / c, a / c }, 1)
                - G({ e / c, 0, b / c, 1, a / c }, 1)
                - G({ e / c, 0, b / c, a / c, 1 }, 1)
                - G({ e / c, 1, 0, b / c, a / c }, 1)
                - G({ e / c, 1, b / c, 0, a / c }, 1)
                - G({ e / c, 1, b / c, a / c, x / c }, 1)
                - G({ e / c, b / c, 0, 1, a / c }, 1)
                - G({ e / c, b / c, 0, a / c, 1 }, 1)
                - G({ e / c, b / c, 1, 0, a / c }, 1)
                - G({ e / c, b / c, 1, a / c, x / c }, 1)
                - G({ e / c, b / c, a / c, 1, x / c }, 1)
                - G({ e / c, b / c, a / c, x / c, 1 }, 1)
                + (sy5 - sy6) * G({ a, e }, { sa, se }, x)
                + G({ b / c, a / c }, 1) * (-sy3 + G({ 0, c, e }, { 1, sc, se }, x))
                + sy1
                    * (sy2 + G({ 0, a, b }, { 1, sa, sb }, x)
                        + G({ a, 0, b }, { sa, 1, sb }, x))
                + G({ b / c }, 1)
                    * (-G({ a, 0, c, e }, { sa, 1, sc, se }, x)
                        + G({ a, b, c, e }, { sa, sb, sc, se }, x))
                + G({ a, b, 0, c, e }, { sa, sb, 1, sc, se }, x)
                + (sy1 * sy10 + sy12 + sy13 + sy14 + sy15 + sy16 + sy17
                      + sy11 * (-sy6 - sy8 - sy9))
                    * Log(c, sc)
                + (-(sy1 * sy10) - sy12 - sy13 - sy14 - sy15 - sy16 - sy17
                      + sy11 * (sy6 + sy8 + sy9))
                    * Log(-x, sc) };
            if (e != x) {
                res += (sy12 - G({ b / c, a / c, x / c, 1 }, 1)) * G({ e }, { se }, x);
            }
            return res;
        }
    } else {

        const complex<double> sy1 = G({ d / c, b / c }, 1);
        const complex<double> sy2 = G({ a, 0, e }, { sa, 1, se }, x);
        const complex<double> sy3 = G({ d / c, e / c }, 1);
        const complex<double> sy4 = G({ a, b, e }, { sa, sb, se }, x);
        const complex<double> sy5 = G({ a, d, e }, { sa, sd, se }, x);
        const complex<double> sy6 = G({ a, e }, { sa, se }, x);
        const complex<double> sy7 = G({ b / c, d / c, a / c }, 1);
        const complex<double> sy8 = G({ 0, a }, { 1, sa }, x);
        const complex<double> sy9 = G({ b / c, d / c, e / c }, 1);
        const complex<double> sy10 = G({ b / c, a / c, d / c }, 1);
        const complex<double> sy11 = G({ d / c, b / c, a / c }, 1);
        const complex<double> sy12 = G({ d / c, b / c, e / c }, 1);
        const complex<double> sy13 = G({ a, b }, { sa, sb }, x);
        const complex<double> sy14 = G({ d / c, e / c, b / c }, 1);
        const complex<double> sy15 = G({ a, b, d, e }, { sa, sb, sd, se }, x);
        const complex<double> sy16 = G({ a }, { sa }, x);
        const complex<double> sy17 = G({ b / c, a / c, d / c, e / c }, 1);
        const complex<double> sy18 = G({ b / c, d / c, a / c, e / c }, 1);
        const complex<double> sy19 = G({ b / c, d / c, e / c, a / c }, 1);
        const complex<double> sy20 = G({ d / c, b / c, a / c, e / c }, 1);
        const complex<double> sy21 = G({ d / c, b / c, e / c, a / c }, 1);
        const complex<double> sy22 = G({ d / c, e / c, b / c, a / c }, 1);
        complex<double> res { -(sy1 * sy2) + sy1 * sy4 - sy6 * sy7 + (sy12 + sy14) * sy8
            + sy8 * sy9 + sy6 * (-sy11 + sy12 + sy9)
            + (-sy2 + sy5) * G({ b / c, d / c }, 1)
            + sy13 * (sy14 + G({ 0, d / c, e / c }, 1) + G({ d / c, 0, e / c }, 1))
            + sy16
                * (-sy19 - sy21 - sy22 - G({ 0, b / c, d / c, e / c }, 1)
                    - G({ 0, d / c, b / c, e / c }, 1) - G({ 0, d / c, e / c, b / c }, 1)
                    - G({ b / c, 0, d / c, e / c }, 1) - G({ b / c, d / c, 0, e / c }, 1)
                    - G({ d / c, 0, b / c, e / c }, 1) - G({ d / c, 0, e / c, b / c }, 1)
                    - G({ d / c, b / c, 0, e / c }, 1) - G({ d / c, e / c, 0, b / c }, 1))
            + G({ 0, b / c, a / c, d / c, e / c }, 1)
            + G({ 0, b / c, d / c, a / c, e / c }, 1)
            + G({ 0, b / c, d / c, e / c, a / c }, 1)
            + G({ 0, d / c, b / c, a / c, e / c }, 1)
            + G({ 0, d / c, b / c, e / c, a / c }, 1)
            + G({ 0, d / c, e / c, b / c, a / c }, 1)
            + G({ b / c, 0, a / c, d / c, e / c }, 1)
            + G({ b / c, 0, d / c, a / c, e / c }, 1)
            + G({ b / c, 0, d / c, e / c, a / c }, 1)
            + G({ b / c, a / c, 0, d / c, e / c }, 1)
            + G({ b / c, a / c, d / c, 0, e / c }, 1)
            + G({ b / c, a / c, d / c, e / c, x / c }, 1)
            + G({ b / c, d / c, 0, a / c, e / c }, 1)
            + G({ b / c, d / c, 0, e / c, a / c }, 1)
            + G({ b / c, d / c, a / c, 0, e / c }, 1)
            + G({ b / c, d / c, a / c, e / c, x / c }, 1)
            + G({ b / c, d / c, e / c, 0, a / c }, 1)
            + G({ b / c, d / c, e / c, a / c, x / c }, 1)
            + G({ d / c, 0, b / c, a / c, e / c }, 1)
            + G({ d / c, 0, b / c, e / c, a / c }, 1)
            + G({ d / c, 0, e / c, b / c, a / c }, 1)
            + G({ d / c, b / c, 0, a / c, e / c }, 1)
            + G({ d / c, b / c, 0, e / c, a / c }, 1)
            + G({ d / c, b / c, a / c, 0, e / c }, 1)
            + G({ d / c, b / c, a / c, e / c, x / c }, 1)
            + G({ d / c, b / c, e / c, 0, a / c }, 1)
            + G({ d / c, b / c, e / c, a / c, x / c }, 1)
            + G({ d / c, e / c, 0, b / c, a / c }, 1)
            + G({ d / c, e / c, b / c, 0, a / c }, 1)
            + G({ d / c, e / c, b / c, a / c, x / c }, 1)
            + (sy10 + sy11 + sy7) * G({ 0, e }, { 1, se }, x)
            + G({ b / c, a / c }, 1) * (-sy5 + G({ 0, d, e }, { 1, sd, se }, x))
            + sy3
                * (-sy4 - G({ 0, a, b }, { 1, sa, sb }, x)
                    - G({ a, 0, b }, { sa, 1, sb }, x))
            + G({ b / c }, 1) * (sy15 - G({ a, 0, d, e }, { sa, 1, sd, se }, x))
            + G({ d / c }, 1) * (-sy15 + G({ a, b, 0, e }, { sa, sb, 1, se }, x))
            + G({ a, b, 0, d, e }, { sa, sb, 1, sd, se }, x)
            + (-sy17 - sy18 - sy19 - sy20 - sy21 - sy22 - sy13 * sy3
                  + sy16 * (sy12 + sy14 + sy9))
                * Log(c, sc)
            + (sy17 + sy18 + sy19 + sy20 + sy21 + sy22 + sy13 * sy3
                  + sy16 * (-sy12 - sy14 - sy9))
                * Log(-x, sc) };
        if (d != x) {
            res += (-sy10 + G({ b / c, a / c, x / c }, 1)) * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (-sy17 - sy18 - sy20 + G({ b / c, a / c, d / c, x / c }, 1)
                       + G({ b / c, d / c, a / c, x / c }, 1)
                       + G({ d / c, b / c, a / c, x / c }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G5_abcde_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{
    if (b == c) {
        if (b == d) {
            if (b == e) { // abbbb

                const complex<double> sy1 = G({ a / b, 1, 1 }, 1);
                const complex<double> sy2 = G({ b, b }, { sb, sb }, x);
                const complex<double> sy3 = G({ b, b, b }, { sb, sb, sb }, x);
                const complex<double> sy4 = G({ a / b }, 1);
                const complex<double> sy5 = G({ a }, { sa }, x);
                const complex<double> sy6 = G({ a / b, 1, 1, 1 }, 1);
                const complex<double> sy7 = G({ b }, { sb }, x);
                complex<double> res { -(sy1 * sy2) - sy5 * sy6 + sy6 * sy7
                    - sy3 * G({ a / b, x / b }, 1) + sy2 * G({ a / b, x / b, 1 }, 1)
                    - sy7 * G({ a / b, x / b, 1, 1 }, 1) - G({ a / b, 0, 1, 1, 1 }, 1)
                    + G({ a / b, x / b, 1, 1, 1 }, 1) + sy1 * G({ a, b }, { sa, sb }, x)
                    + G({ a / b, 1 }, 1) * (sy3 - G({ a, b, b }, { sa, sb, sb }, x))
                    - sy4 * G({ 0, b, b, b }, { 1, sb, sb, sb }, x)
                    + sy4 * G({ a, b, b, b }, { sa, sb, sb, sb }, x)
                    + G({ a, 0, b, b, b }, { sa, 1, sb, sb, sb }, x) - sy5 * Zeta(4) };
                return res;
            } else { // abbbe

                const complex<double> sy1 = G({ b, e }, { sb, se }, x);
                const complex<double> sy2 = G({ a / b, 1, 1 }, 1);
                const complex<double> sy3 = G({ b, b, e }, { sb, sb, se }, x);
                const complex<double> sy4 = G({ e / b, 1, 1 }, 1);
                const complex<double> sy5 = G({ a }, { sa }, x);
                const complex<double> sy6 = G({ e / b, 1, 1, a / b }, 1);
                const complex<double> sy7 = G({ e / b, 1, a / b, 1 }, 1);
                const complex<double> sy8 = G({ e / b, a / b, 1, 1 }, 1);
                const complex<double> sy9 = G({ a / b, e / b, 1, 1 }, 1);

                complex<double> res { -(sy1 * sy2) - sy3 * G({ a / b, x / b }, 1)
                    + sy1 * G({ a / b, x / b, 1 }, 1)
                    + sy5 * (sy6 + sy7 + sy8 + G({ 0, e / b, 1, 1 }, 1))
                    - G({ 0, a / b, e / b, 1, 1 }, 1) - G({ 0, e / b, 1, 1, a / b }, 1)
                    - G({ 0, e / b, 1, a / b, 1 }, 1) - G({ 0, e / b, a / b, 1, 1 }, 1)
                    - G({ a / b, 0, e / b, 1, 1 }, 1)
                    - G({ a / b, e / b, 1, 1, x / b }, 1)
                    - G({ a / b, e / b, 1, x / b, 1 }, 1)
                    - G({ a / b, e / b, x / b, 1, 1 }, 1)
                    - G({ e / b, 0, 1, 1, a / b }, 1) - G({ e / b, 0, 1, a / b, 1 }, 1)
                    - G({ e / b, 0, a / b, 1, 1 }, 1) - G({ e / b, 1, 0, 1, a / b }, 1)
                    - G({ e / b, 1, 0, a / b, 1 }, 1) - G({ e / b, 1, 1, 0, a / b }, 1)
                    - G({ e / b, 1, 1, a / b, x / b }, 1)
                    - G({ e / b, 1, a / b, 1, x / b }, 1)
                    - G({ e / b, 1, a / b, x / b, 1 }, 1)
                    - G({ e / b, a / b, 1, 1, x / b }, 1)
                    - G({ e / b, a / b, 1, x / b, 1 }, 1)
                    - G({ e / b, a / b, x / b, 1, 1 }, 1)
                    - sy4 * G({ 0, a }, { 1, sa }, x)
                    + (sy2 - sy4) * G({ a, e }, { sa, se }, x)
                    + G({ a / b, 1 }, 1) * (sy3 - G({ a, b, e }, { sa, sb, se }, x))
                    + G({ a / b }, 1)
                        * (-G({ 0, b, b, e }, { 1, sb, sb, se }, x)
                            + G({ a, b, b, e }, { sa, sb, sb, se }, x))
                    + G({ a, 0, b, b, e }, { sa, 1, sb, sb, se }, x)
                    + (-(sy4 * sy5) + sy6 + sy7 + sy8 + sy9) * Log(b, sb)
                    + (sy4 * sy5 - sy6 - sy7 - sy8 - sy9) * Log(-x, sb) };

                if (e != x) {
                    res += (sy9 - G({ a / b, x / b, 1, 1 }, 1)) * G({ e }, { se }, x);
                }
                return res;
            }
        } else { // abbde

            const complex<double> sy1 = G({ a, d, e }, { sa, sd, se }, x);
            const complex<double> sy2 = G({ b, d, e }, { sb, sd, se }, x);
            const complex<double> sy3 = G({ a, e }, { sa, se }, x);
            const complex<double> sy4 = G({ d / b, 1, a / b }, 1);
            const complex<double> sy5 = G({ 0, a }, { 1, sa }, x);
            const complex<double> sy6 = G({ d / b, 1, e / b }, 1);
            const complex<double> sy7 = G({ a / b, d / b, 1 }, 1);
            const complex<double> sy8 = G({ d / b, a / b, 1 }, 1);
            const complex<double> sy9 = G({ d / b, e / b, 1 }, 1);
            const complex<double> sy10 = G({ a }, { sa }, x);
            const complex<double> sy11 = G({ d / b, 1, e / b, a / b }, 1);
            const complex<double> sy12 = G({ d / b, e / b, 1, a / b }, 1);
            const complex<double> sy13 = G({ d / b, e / b, a / b, 1 }, 1);
            const complex<double> sy14 = G({ a / b, d / b, 1, e / b }, 1);
            const complex<double> sy15 = G({ a / b, d / b, e / b, 1 }, 1);
            const complex<double> sy16 = G({ d / b, 1, a / b, e / b }, 1);
            const complex<double> sy17 = G({ d / b, a / b, 1, e / b }, 1);
            const complex<double> sy18 = G({ d / b, a / b, e / b, 1 }, 1);
            complex<double> res { -(sy3 * sy4) + sy5 * sy6 + sy5 * sy9
                + sy3 * (sy6 - sy8 + sy9) + (-sy1 + sy2) * G({ a / b, 1 }, 1)
                - sy2 * G({ a / b, x / b }, 1)
                + sy10
                    * (-sy11 - sy12 - sy13 - G({ 0, d / b, 1, e / b }, 1)
                        - G({ 0, d / b, e / b, 1 }, 1) - G({ d / b, 0, 1, e / b }, 1)
                        - G({ d / b, 0, e / b, 1 }, 1) - G({ d / b, 1, 0, e / b }, 1))
                + G({ 0, a / b, d / b, 1, e / b }, 1)
                + G({ 0, a / b, d / b, e / b, 1 }, 1)
                + G({ 0, d / b, 1, a / b, e / b }, 1)
                + G({ 0, d / b, 1, e / b, a / b }, 1)
                + G({ 0, d / b, a / b, 1, e / b }, 1)
                + G({ 0, d / b, a / b, e / b, 1 }, 1)
                + G({ 0, d / b, e / b, 1, a / b }, 1)
                + G({ 0, d / b, e / b, a / b, 1 }, 1)
                + G({ a / b, 0, d / b, 1, e / b }, 1)
                + G({ a / b, 0, d / b, e / b, 1 }, 1)
                + G({ a / b, d / b, 0, 1, e / b }, 1)
                + G({ a / b, d / b, 0, e / b, 1 }, 1)
                + G({ a / b, d / b, 1, 0, e / b }, 1)
                + G({ a / b, d / b, 1, e / b, x / b }, 1)
                + G({ a / b, d / b, e / b, 1, x / b }, 1)
                + G({ a / b, d / b, e / b, x / b, 1 }, 1)
                + G({ d / b, 0, 1, a / b, e / b }, 1)
                + G({ d / b, 0, 1, e / b, a / b }, 1)
                + G({ d / b, 0, a / b, 1, e / b }, 1)
                + G({ d / b, 0, a / b, e / b, 1 }, 1)
                + G({ d / b, 0, e / b, 1, a / b }, 1)
                + G({ d / b, 0, e / b, a / b, 1 }, 1)
                + G({ d / b, 1, 0, a / b, e / b }, 1)
                + G({ d / b, 1, 0, e / b, a / b }, 1)
                + G({ d / b, 1, a / b, 0, e / b }, 1)
                + G({ d / b, 1, a / b, e / b, x / b }, 1)
                + G({ d / b, 1, e / b, 0, a / b }, 1)
                + G({ d / b, 1, e / b, a / b, x / b }, 1)
                + G({ d / b, a / b, 0, 1, e / b }, 1)
                + G({ d / b, a / b, 0, e / b, 1 }, 1)
                + G({ d / b, a / b, 1, 0, e / b }, 1)
                + G({ d / b, a / b, 1, e / b, x / b }, 1)
                + G({ d / b, a / b, e / b, 1, x / b }, 1)
                + G({ d / b, a / b, e / b, x / b, 1 }, 1)
                + G({ d / b, e / b, 0, 1, a / b }, 1)
                + G({ d / b, e / b, 0, a / b, 1 }, 1)
                + G({ d / b, e / b, 1, 0, a / b }, 1)
                + G({ d / b, e / b, 1, a / b, x / b }, 1)
                + G({ d / b, e / b, a / b, 1, x / b }, 1)
                + G({ d / b, e / b, a / b, x / b, 1 }, 1)
                + (sy4 + sy7 + sy8) * G({ 0, e }, { 1, se }, x)
                + G({ d / b, 1 }, 1) * (sy1 - G({ a, 0, e }, { sa, 1, se }, x))
                + G({ a / b }, 1)
                    * (-G({ 0, b, d, e }, { 1, sb, sd, se }, x)
                        + G({ a, b, d, e }, { sa, sb, sd, se }, x))
                + G({ a, 0, b, d, e }, { sa, 1, sb, sd, se }, x)
                + (-sy11 - sy12 - sy13 - sy14 - sy15 - sy16 - sy17 - sy18
                      + sy10 * (sy6 + sy9))
                    * Log(b, sb)
                + (sy11 + sy12 + sy13 + sy14 + sy15 + sy16 + sy17 + sy18
                      + sy10 * (-sy6 - sy9))
                    * Log(-x, sb) };
            if (d != x) {
                res += (-sy7 + G({ a / b, x / b, 1 }, 1)) * G({ d, e }, { sd, se }, x);
            }
            if (e != x) {
                res += (-sy14 - sy15 - sy16 - sy17 - sy18
                           + G({ a / b, d / b, 1, x / b }, 1)
                           + G({ a / b, d / b, x / b, 1 }, 1)
                           + G({ d / b, 1, a / b, x / b }, 1)
                           + G({ d / b, a / b, 1, x / b }, 1)
                           + G({ d / b, a / b, x / b, 1 }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else { // abcde

        const complex<double> sy1 = G({ a / b, c / b }, 1);
        const complex<double> sy2 = G({ 0, d, e }, { 1, sd, se }, x);
        const complex<double> sy3 = G({ c / b, a / b }, 1);
        const complex<double> sy4 = G({ a, d, e }, { sa, sd, se }, x);
        const complex<double> sy5 = G({ a / b, c / b, d / b }, 1);
        const complex<double> sy6 = G({ c / b, a / b, d / b }, 1);
        const complex<double> sy7 = G({ c / b, d / b, a / b }, 1);
        const complex<double> sy8 = G({ a, e }, { sa, se }, x);
        const complex<double> sy9 = G({ c / b, d / b, e / b }, 1);
        const complex<double> sy10 = G({ a, c, d, e }, { sa, sc, sd, se }, x);
        const complex<double> sy11 = G({ a }, { sa }, x);
        const complex<double> sy12 = G({ c / b, d / b, e / b, a / b }, 1);
        const complex<double> sy13 = G({ a / b, c / b, d / b, e / b }, 1);
        const complex<double> sy14 = G({ c / b, a / b, d / b, e / b }, 1);
        const complex<double> sy15 = G({ c / b, d / b, a / b, e / b }, 1);
        complex<double> res { -(sy1 * sy2) - sy2 * sy3 + sy3 * sy4 + sy7 * sy8 - sy8 * sy9
            + sy11
                * (sy12 + G({ 0, c / b, d / b, e / b }, 1)
                    + G({ c / b, 0, d / b, e / b }, 1) + G({ c / b, d / b, 0, e / b }, 1))
            - G({ 0, a / b, c / b, d / b, e / b }, 1)
            - G({ 0, c / b, a / b, d / b, e / b }, 1)
            - G({ 0, c / b, d / b, a / b, e / b }, 1)
            - G({ 0, c / b, d / b, e / b, a / b }, 1)
            - G({ a / b, 0, c / b, d / b, e / b }, 1)
            - G({ a / b, c / b, 0, d / b, e / b }, 1)
            - G({ a / b, c / b, d / b, 0, e / b }, 1)
            - G({ a / b, c / b, d / b, e / b, x / b }, 1)
            - G({ c / b, 0, a / b, d / b, e / b }, 1)
            - G({ c / b, 0, d / b, a / b, e / b }, 1)
            - G({ c / b, 0, d / b, e / b, a / b }, 1)
            - G({ c / b, a / b, 0, d / b, e / b }, 1)
            - G({ c / b, a / b, d / b, 0, e / b }, 1)
            - G({ c / b, a / b, d / b, e / b, x / b }, 1)
            - G({ c / b, d / b, 0, a / b, e / b }, 1)
            - G({ c / b, d / b, 0, e / b, a / b }, 1)
            - G({ c / b, d / b, a / b, 0, e / b }, 1)
            - G({ c / b, d / b, a / b, e / b, x / b }, 1)
            - G({ c / b, d / b, e / b, 0, a / b }, 1)
            - G({ c / b, d / b, e / b, a / b, x / b }, 1)
            - sy9 * G({ 0, a }, { 1, sa }, x)
            + (-sy5 - sy6 - sy7) * G({ 0, e }, { 1, se }, x)
            + G({ c / b, d / b }, 1) * (-sy4 + G({ a, 0, e }, { sa, 1, se }, x))
            + G({ a / b }, 1) * (sy10 - G({ 0, c, d, e }, { 1, sc, sd, se }, x))
            + G({ c / b }, 1) * (-sy10 + G({ a, 0, d, e }, { sa, 1, sd, se }, x))
            + G({ a, 0, c, d, e }, { sa, 1, sc, sd, se }, x)
            + (sy12 + sy13 + sy14 + sy15 - sy11 * sy9) * Log(b, sb)
            + (-sy12 - sy13 - sy14 - sy15 + sy11 * sy9) * Log(-x, sb) };
        if (c != x) {
            res += (sy1 - G({ a / b, x / b }, 1)) * G({ c, d, e }, { sc, sd, se }, x);
        }
        if (d != x) {
            res += (sy5 + sy6 - G({ a / b, c / b, x / b }, 1)
                       - G({ c / b, a / b, x / b }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (sy13 + sy14 + sy15 - G({ a / b, c / b, d / b, x / b }, 1)
                       - G({ c / b, a / b, d / b, x / b }, 1)
                       - G({ c / b, d / b, a / b, x / b }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G5_abcde_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{

    if (a == b) {
        if (a == c) {
            if (a == d) { // aaaae
                const complex<double> sy1 = G({ e / a, 1, 1, 1 }, 1);
                complex<double> res { -G({ 0, e / a, 1, 1, 1 }, 1)
                    - G({ e / a, 1, 1, 1, x / a }, 1) - G({ e / a, 1, 1, x / a, 1 }, 1)
                    - G({ e / a, 1, x / a, 1, 1 }, 1) - G({ e / a, x / a, 1, 1, 1 }, 1)
                    + G({ x / a, 1, 1 }, 1) * G({ a, e }, { sa, se }, x)
                    - G({ x / a, 1 }, 1) * G({ a, a, e }, { sa, sa, se }, x)
                    + G({ x / a }, 1) * G({ a, a, a, e }, { sa, sa, sa, se }, x)
                    + G({ 0, a, a, a, e }, { 1, sa, sa, sa, se }, x) + sy1 * Log(a, sa)
                    - sy1 * Log(-x, sa) };
                if (e != x) {
                    res += (sy1 - G({ x / a, 1, 1, 1 }, 1)) * G({ e }, { se }, x);
                }
                return res;
            } else { // aaade
                const complex<double> sy1 = G({ d / a, 1, 1 }, 1);
                const complex<double> sy2 = G({ d / a, 1, 1, e / a }, 1);
                const complex<double> sy3 = G({ d / a, 1, e / a, 1 }, 1);
                const complex<double> sy4 = G({ d / a, e / a, 1, 1 }, 1);
                complex<double> res { G({ 0, d / a, 1, 1, e / a }, 1)
                    + G({ 0, d / a, 1, e / a, 1 }, 1) + G({ 0, d / a, e / a, 1, 1 }, 1)
                    + G({ d / a, 0, 1, 1, e / a }, 1) + G({ d / a, 0, 1, e / a, 1 }, 1)
                    + G({ d / a, 0, e / a, 1, 1 }, 1) + G({ d / a, 1, 0, 1, e / a }, 1)
                    + G({ d / a, 1, 0, e / a, 1 }, 1) + G({ d / a, 1, 1, 0, e / a }, 1)
                    + G({ d / a, 1, 1, e / a, x / a }, 1)
                    + G({ d / a, 1, e / a, 1, x / a }, 1)
                    + G({ d / a, 1, e / a, x / a, 1 }, 1)
                    + G({ d / a, e / a, 1, 1, x / a }, 1)
                    + G({ d / a, e / a, 1, x / a, 1 }, 1)
                    + G({ d / a, e / a, x / a, 1, 1 }, 1)
                    + sy1 * G({ 0, e }, { 1, se }, x)
                    - G({ x / a, 1 }, 1) * G({ a, d, e }, { sa, sd, se }, x)
                    + G({ x / a }, 1) * G({ a, a, d, e }, { sa, sa, sd, se }, x)
                    + G({ 0, a, a, d, e }, { 1, sa, sa, sd, se }, x)
                    + (-sy2 - sy3 - sy4) * Log(a, sa) + (sy2 + sy3 + sy4) * Log(-x, sa) };
                if (d != x) {
                    res += (-sy1 + G({ x / a, 1, 1 }, 1)) * G({ d, e }, { sd, se }, x);
                }
                if (e != x) {
                    res += (-sy2 - sy3 - sy4 + G({ d / a, 1, 1, x / a }, 1)
                               + G({ d / a, 1, x / a, 1 }, 1)
                               + G({ d / a, x / a, 1, 1 }, 1))
                        * G({ e }, { se }, x);
                }
                return res;
            }
        } else { // aacde
            const complex<double> sy1 = G({ c / a, 1 }, 1);
            const complex<double> sy2 = G({ c / a, 1, d / a }, 1);
            const complex<double> sy3 = G({ c / a, d / a, 1 }, 1);
            const complex<double> sy4 = G({ c / a, 1, d / a, e / a }, 1);
            const complex<double> sy5 = G({ c / a, d / a, 1, e / a }, 1);
            const complex<double> sy6 = G({ c / a, d / a, e / a, 1 }, 1);
            complex<double> res { -G({ 0, c / a, 1, d / a, e / a }, 1)
                - G({ 0, c / a, d / a, 1, e / a }, 1)
                - G({ 0, c / a, d / a, e / a, 1 }, 1)
                - G({ c / a, 0, 1, d / a, e / a }, 1)
                - G({ c / a, 0, d / a, 1, e / a }, 1)
                - G({ c / a, 0, d / a, e / a, 1 }, 1)
                - G({ c / a, 1, 0, d / a, e / a }, 1)
                - G({ c / a, 1, d / a, 0, e / a }, 1)
                - G({ c / a, 1, d / a, e / a, x / a }, 1)
                - G({ c / a, d / a, 0, 1, e / a }, 1)
                - G({ c / a, d / a, 0, e / a, 1 }, 1)
                - G({ c / a, d / a, 1, 0, e / a }, 1)
                - G({ c / a, d / a, 1, e / a, x / a }, 1)
                - G({ c / a, d / a, e / a, 1, x / a }, 1)
                - G({ c / a, d / a, e / a, x / a, 1 }, 1)
                + (-sy2 - sy3) * G({ 0, e }, { 1, se }, x)
                - sy1 * G({ 0, d, e }, { 1, sd, se }, x)
                + G({ x / a }, 1) * G({ a, c, d, e }, { sa, sc, sd, se }, x)
                + G({ 0, a, c, d, e }, { 1, sa, sc, sd, se }, x)
                + (sy4 + sy5 + sy6) * Log(a, sa) + (-sy4 - sy5 - sy6) * Log(-x, sa) };
            if (c != x) {
                res += (sy1 - G({ x / a, 1 }, 1)) * G({ c, d, e }, { sc, sd, se }, x);
            }
            if (d != x) {
                res += (sy2 + sy3 - G({ c / a, 1, x / a }, 1) - G({ c / a, x / a, 1 }, 1))
                    * G({ d, e }, { sd, se }, x);
            }
            if (e != x) {
                res += (sy4 + sy5 + sy6 - G({ c / a, 1, d / a, x / a }, 1)
                           - G({ c / a, d / a, 1, x / a }, 1)
                           - G({ c / a, d / a, x / a, 1 }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else { // abcde

        const complex<double> sy1 = G({ b / a, c / a }, 1);
        const complex<double> sy2 = G({ b / a, c / a, d / a }, 1);
        const complex<double> sy3 = G({ b / a }, 1);
        const complex<double> sy4 = G({ b / a, c / a, d / a, e / a }, 1);
        complex<double> res { G({ 0, b / a, c / a, d / a, e / a }, 1)
            + G({ b / a, 0, c / a, d / a, e / a }, 1)
            + G({ b / a, c / a, 0, d / a, e / a }, 1)
            + G({ b / a, c / a, d / a, 0, e / a }, 1)
            + G({ b / a, c / a, d / a, e / a, x / a }, 1)
            + sy2 * G({ 0, e }, { 1, se }, x) + sy1 * G({ 0, d, e }, { 1, sd, se }, x)
            + sy3 * G({ 0, c, d, e }, { 1, sc, sd, se }, x)
            + G({ 0, b, c, d, e }, { 1, sb, sc, sd, se }, x) - sy4 * Log(a, sa)
            + sy4 * Log(-x, sa) };
        if (b != x) {
            res += (-sy3 + G({ x / a }, 1)) * G({ b, c, d, e }, { sb, sc, sd, se }, x);
        }
        if (c != x) {
            res += (-sy1 + G({ b / a, x / a }, 1)) * G({ c, d, e }, { sc, sd, se }, x);
        }
        if (d != x) {
            res += (-sy2 + G({ b / a, c / a, x / a }, 1)) * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (-sy4 + G({ b / a, c / a, d / a, x / a }, 1)) * G({ e }, { se }, x);
        }
        return res;
    }
}

// G(0,0,a,b,c,x)
// This should only be called if abs(a) < x and/or abs(b) < x and/or abs(c) < x
complex<double> G5_00abc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_00abc_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_00abc_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G5_00abc_b(a, b, c, sa, sb, sc, x);
}
// take care of G({0,a,0,b,c},x)
complex<double> G5_0a0bc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_0a0bc_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_0a0bc_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G5_0a0bc_b(a, b, c, sa, sb, sc, x);
}

// take care of G({0,a,b,0,c},x)
complex<double> G5_0ab0c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_0ab0c_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_0ab0c_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G5_0ab0c_b(a, b, c, sa, sb, sc, x);
}
// take care of G({a,0,0,b,c},x)
complex<double> G5_a00bc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_a00bc_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_a00bc_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G5_a00bc_b(a, b, c, sa, sb, sc, x);
}
// take care of G({a,0,b,0,c},x)
complex<double> G5_a0b0c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_a0b0c_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_a0b0c_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G5_a0b0c_b(a, b, c, sa, sb, sc, x);
}
// take care of G({a,b,0,0,c},x)
complex<double> G5_ab00c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x)
{
    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_ab00c_c(a, b, c, sa, sb, sc, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_ab00c_a(a, b, c, sa, sb, sc, x);

    // b is smallest
    return G5_ab00c_b(a, b, c, sa, sb, sc, x);
}

// take care of G({0,a,b,c,d},x)
complex<double> G5_0abcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G5_0abcd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_0abcd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_0abcd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G5_0abcd_b(a, b, c, d, sa, sb, sc, sd, x);
}
// take care of G({a,0,b,c,d},x)
complex<double> G5_a0bcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G5_a0bcd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_a0bcd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_a0bcd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G5_a0bcd_b(a, b, c, d, sa, sb, sc, sd, x);
}
// take care of G({a,b,0,c,d},x)
complex<double> G5_ab0cd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G5_ab0cd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_ab0cd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_ab0cd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G5_ab0cd_b(a, b, c, d, sa, sb, sc, sd, x);
}
// take care of G({a,b,c,0,d},x)
complex<double> G5_abc0d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G5_abc0d_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_abc0d_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_abc0d_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G5_abc0d_b(a, b, c, d, sa, sb, sc, sd, x);
}
// take care of G({a,b,c,d,e},x)
complex<double> G5_abcde(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{
    // e is smallest
    if (abs(e) < abs(d) && abs(e) < abs(c) && abs(e) < abs(b) && abs(e) < abs(a))
        return G5_abcde_e(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G5_abcde_d(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G5_abcde_c(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G5_abcde_a(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // b is smallest
    return G5_abcde_b(a, b, c, d, e, sa, sb, sc, sd, se, x);
}
