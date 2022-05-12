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

complex<double> G6_explicit_0000aa(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return -pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a }, { sa }, x) };

    return sy[0] * G({ 0, 0, 0, 0, x / a }, 1) - G({ 0, 0, 0, 0, x / a, 1 }, 1)
        + G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x)
        - G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) * Zeta(2)
        - G({ 0, a }, { 1, sa }, x) * Zeta(4) - Zeta(6)
        + G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(3) + sy[0] * Zeta(5);
}
complex<double> G6_explicit_a0000a(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a }, { sa }, x) };
    return sy[0] * G({ 0, 0, 0, 0, x / a }, 1) + G({ 0, 0, 0, 0, 1, x / a }, 1)
        + G({ 0, 0, 0, x / a }, 1) * G({ 0, a }, { 1, sa }, x)
        + G({ 0, 0, x / a }, 1) * G({ 0, 0, a }, { 1, 1, sa }, x)
        + G({ 0, x / a }, 1) * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
        + G({ x / a }, 1) * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
        + G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x) - 5. * Zeta(6)
        + sy[0] * Zeta(5) + Log(a, sa) * Zeta(5) - Log(-x, sa) * Zeta(5);
}
complex<double> G6_explicit_000a0a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ 0, a }, { 1, sa }, x), G({ a }, { sa }, x),
        Log(a, sa), G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x),
        G({ 0, 0, 0 }, { 1, 1, 1 }, x) };
    return -(sy[0] * G({ 0, 0, 0, x / a }, 1)) - 4. * sy[1] * G({ 0, 0, 0, 0, x / a }, 1)
        - 4. * G({ 0, 0, 0, 0, 1, x / a }, 1) - 3. * G({ 0, 0, 0, 1, 0, x / a }, 1)
        - 2. * G({ 0, 0, 1, 0, 0, x / a }, 1) - G({ 0, 1, 0, 0, 0, x / a }, 1)
        + G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x)
        - 2. * sy[3] * pow(Zeta(2), 2.) - sy[4] * sy[5] * Zeta(2)
        + G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x) * Zeta(2)
        + G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) * Zeta(2)
        + (Log(-x, sa) * pow(sy[2], 3.) * Zeta(2)) / 6. - (pow(sy[2], 4.) * Zeta(2)) / 24.
        - pow(sy[3], 2.) * Zeta(2)
        + (pow(sy[2], 2.) * Zeta(2) * (sy[3] + 2. * Zeta(2))) / 2.
        + pow(sy[4], 2.) * Zeta(2) * (sy[3] + 2. * Zeta(2))
        - sy[2] * Zeta(2) * (-sy[5] + sy[4] * (sy[3] + 2. * Zeta(2)))
        + 3. * sy[0] * Zeta(4) + 2. * Zeta(2) * Zeta(4) + 5. * Zeta(6)
        - 2. * G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(3) - 4. * sy[1] * Zeta(5);
}
complex<double> G6_explicit_0a000a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return -pow(Zeta(3), 2.) + (25. * Zeta(6)) / 12.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ 0, a }, { 1, sa }, x), G({ a }, { sa }, x), Log(a, sa) };
    return -3. * sy[0] * G({ 0, 0, 0, x / a }, 1)
        - 4. * sy[1] * G({ 0, 0, 0, 0, x / a }, 1) - 4. * G({ 0, 0, 0, 0, 1, x / a }, 1)
        - G({ 0, 0, 0, 1, 0, x / a }, 1)
        - 2. * G({ 0, 0, x / a }, 1) * G({ 0, 0, a }, { 1, 1, sa }, x)
        - G({ 0, x / a }, 1) * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
        + G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x) + sy[0] * Zeta(4)
        + G({ 0, 0 }, { 1, 1 }, x) * Zeta(4) + sy[2] * Log(-x, sa) * Zeta(4)
        - (pow(sy[2], 2.) * Zeta(4)) / 2. + 2. * Zeta(2) * Zeta(4) + 10. * Zeta(6)
        - 4. * sy[1] * Zeta(5);
}
complex<double> G6_explicit_00a00a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, a }, { 1, sa }, x),
              G({ a }, { sa }, x), Log(a, sa), G({ 0, 0 }, { 1, 1 }, x) };
    return sy[0] * G({ 0, 0, x / a }, 1) + 3. * sy[1] * G({ 0, 0, 0, x / a }, 1)
        + 6. * sy[2] * G({ 0, 0, 0, 0, x / a }, 1) + 6. * G({ 0, 0, 0, 0, 1, x / a }, 1)
        + 3. * G({ 0, 0, 0, 1, 0, x / a }, 1) + G({ 0, 0, 1, 0, 0, x / a }, 1)
        + G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x) - 3. * sy[1] * Zeta(4)
        - 10. * Zeta(6) + sy[0] * Zeta(3) - G({ 0, 0, 0 }, { 1, 1, 1 }, x) * Zeta(3)
        - (Log(-x, sa) * pow(sy[3], 2.) * Zeta(3)) / 2. + (pow(sy[3], 3.) * Zeta(3)) / 6.
        - sy[3] * (sy[4] + 2. * Zeta(2)) * Zeta(3)
        + G({ 0 }, { 1 }, x) * (sy[4] + 2. * Zeta(2)) * Zeta(3) + 6. * sy[2] * Zeta(5);
}

complex<double> G6_explicit_000aaa(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return pow(Zeta(3), 2.) - (23. * Zeta(6)) / 16.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a, a }, { sa, sa }, x), G({ a }, { sa }, x),
        G({ 0, 0, a }, { 1, 1, sa }, x) };
    return -(sy[0] * G({ 0, 0, 0, x / a }, 1)) + sy[1] * G({ 0, 0, 0, x / a, 1 }, 1)
        - G({ 0, 0, 0, x / a, 1, 1 }, 1)
        + G({ 0, 0, 0, 0, a, a }, { 1, 1, 1, 1, sa, sa }, x) - pow(Zeta(3), 2.) / 2.
        + G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) * Zeta(2)
        - G({ 0, 0, a, a }, { 1, 1, sa, sa }, x) * Zeta(2) - sy[0] * Zeta(4)
        + (5. * G({ 0, a }, { 1, sa }, x) * Zeta(4)) / 4. + (3. * Zeta(6)) / 4.
        - sy[2] * Zeta(3) - (sy[2] - G({ 0, a, a }, { 1, sa, sa }, x)) * Zeta(3)
        - sy[1] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5));
}
complex<double> G6_explicit_a000aa(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ a, a }, { sa, sa }, x), G({ a }, { sa }, x) };
    return sy[0] * G({ 0, 0, 0, x / a }, 1) + sy[1] * G({ 0, 0, 0, 1, x / a }, 1)
        + G({ 0, 0, 0, 1, 1, x / a }, 1)
        + G({ 0, 0, x / a }, 1) * G({ 0, a, a }, { 1, sa, sa }, x)
        + G({ 0, x / a }, 1) * G({ 0, 0, a, a }, { 1, 1, sa, sa }, x)
        + G({ x / a }, 1) * G({ 0, 0, 0, a, a }, { 1, 1, 1, sa, sa }, x)
        + G({ 0, 0, 0, 0, a, a }, { 1, 1, 1, 1, sa, sa }, x) + pow(Zeta(3), 2.)
        + sy[0] * Zeta(4) - G({ 0, a }, { 1, sa }, x) * Zeta(4)
        + 4. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.) - (4. * Zeta(6)) / 3.
        - sy[1] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5))
        - Log(a, sa) * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5))
        + Log(-x, sa) * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5));
}
complex<double> G6_explicit_aa000a(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ 0, a }, { 1, sa }, x), G({ a }, { sa }, x), Log(a, sa), Log(-x, sa) };
    return -(sy[0] * G({ 0, 0, 1, x / a }, 1)) - sy[0] * G({ 0, 0, x / a, 1 }, 1)
        - sy[0] * G({ 0, 1, 0, x / a }, 1) - sy[1] * G({ 0, 0, 0, 1, x / a }, 1)
        - sy[1] * G({ 0, 0, 0, x / a, 1 }, 1) - sy[1] * G({ 0, 0, 1, 0, x / a }, 1)
        - sy[1] * G({ 0, 1, 0, 0, x / a }, 1) - 2. * G({ 0, 0, 0, 1, 1, x / a }, 1)
        - G({ 0, 0, 0, 1, x / a, 1 }, 1) - G({ 0, 0, 1, 0, 1, x / a }, 1)
        - G({ 0, 1, 0, 0, 1, x / a }, 1)
        - (G({ 0, 1, x / a }, 1) + G({ 0, x / a, 1 }, 1))
        * G({ 0, 0, a }, { 1, 1, sa }, x)
        - G({ x / a, 1 }, 1) * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
        + G({ x / a }, 1) * G({ a, 0, 0, 0, a }, { sa, 1, 1, 1, sa }, x)
        + G({ 0, a, 0, 0, 0, a }, { 1, sa, 1, 1, 1, sa }, x)
        - 4. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        - 4. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        - 8. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.)
        - 3. * (-pow(Zeta(3), 2.) + (25. * Zeta(6)) / 12.)
        + sy[1] * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        + sy[2] * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        - sy[3] * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        + 2. * sy[1] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5))
        + 2. * sy[2] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5))
        - 2. * sy[3] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5))
        + sy[1] * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.)
        + sy[2] * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.)
        - sy[3] * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_0a00aa(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return (3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ 0, a, a }, { 1, sa, sa }, x),
        G({ a, a }, { sa, sa }, x), G({ a }, { sa }, x), Log(a, sa) };
    return -2. * sy[0] * G({ 0, 0, x / a }, 1) - 3. * sy[1] * G({ 0, 0, 0, x / a }, 1)
        - 3. * sy[2] * G({ 0, 0, 0, 1, x / a }, 1) - sy[2] * G({ 0, 0, 1, 0, x / a }, 1)
        - 3. * G({ 0, 0, 0, 1, 1, x / a }, 1) - G({ 0, 0, 1, 0, 1, x / a }, 1)
        - G({ 0, 0, 1, 1, 0, x / a }, 1)
        - G({ 0, x / a }, 1) * G({ 0, 0, a, a }, { 1, 1, sa, sa }, x)
        + G({ 0, 0, 0, 0, a, a }, { 1, 1, 1, 1, sa, sa }, x) - pow(Zeta(3), 2.) / 2.
        - 3. * sy[1] * Zeta(4) - (G({ 0, 0 }, { 1, 1 }, x) * Zeta(4)) / 4.
        + (11. * G({ 0, a }, { 1, sa }, x) * Zeta(4)) / 4.
        - (sy[3] * Log(-x, sa) * Zeta(4)) / 4. + (pow(sy[3], 2.) * Zeta(4)) / 8.
        - (Zeta(2) * Zeta(4)) / 2. - 3. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        - 6. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.) + Zeta(6) / 2.
        - (-sy[0] + G({ 0, 0, a }, { 1, 1, sa }, x)) * Zeta(3)
        + sy[2] * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        + 3. * sy[2] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5));
}
complex<double> G6_explicit_00a0aa(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return -3. * pow(Zeta(3), 2.) + (203. * Zeta(6)) / 48.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ 0, a, a }, { 1, sa, sa }, x),
        G({ a, a }, { sa, sa }, x), G({ a }, { sa }, x), Log(a, sa), G({ 0 }, { 1 }, x),
        G({ 0, 0 }, { 1, 1 }, x), G({ 0, 0, a }, { 1, 1, sa }, x) };

    return sy[0] * G({ 0, 0, x / a }, 1) + 3. * sy[1] * G({ 0, 0, 0, x / a }, 1)
        + 3. * sy[2] * G({ 0, 0, 0, 1, x / a }, 1)
        + 2. * sy[2] * G({ 0, 0, 1, 0, x / a }, 1) + sy[2] * G({ 0, 1, 0, 0, x / a }, 1)
        + 3. * G({ 0, 0, 0, 1, 1, x / a }, 1) + 2. * G({ 0, 0, 1, 0, 1, x / a }, 1)
        + 2. * G({ 0, 0, 1, 1, 0, x / a }, 1) + G({ 0, 1, 0, 0, 1, x / a }, 1)
        + G({ 0, 1, 0, 1, 0, x / a }, 1) + G({ 0, 1, 1, 0, 0, x / a }, 1)
        + G({ 0, 0, 0, 0, a, a }, { 1, 1, 1, 1, sa, sa }, x) - pow(Zeta(3), 2.)
        - G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) * Zeta(2)
        + G({ 0, 0, a, a }, { 1, 1, sa, sa }, x) * Zeta(2) + 3. * sy[1] * Zeta(4)
        - (7. * G({ 0, a }, { 1, sa }, x) * Zeta(4)) / 4.
        + 3. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        + 2. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        + 4. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.) + (25. * Zeta(6)) / 12.
        - sy[4] * sy[5] * Zeta(3) - sy[6] * Zeta(3) + 2. * (-sy[0] + sy[6]) * Zeta(3)
        + G({ 0, 0, 0 }, { 1, 1, 1 }, x) * Zeta(3)
        + (Log(-x, sa) * pow(sy[3], 2.) * Zeta(3)) / 2. - (pow(sy[3], 3.) * Zeta(3)) / 6.
        - 2. * sy[4] * Zeta(2) * Zeta(3) + sy[3] * (sy[5] + 2. * Zeta(2)) * Zeta(3)
        - 2. * sy[2] * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        - 3. * sy[2] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5))
        - sy[2] * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_0aa00a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return -pow(Zeta(3), 2.) + (13. * Zeta(6)) / 16.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ 0, a }, { 1, sa }, x),
        G({ a, 0, 0, a }, { sa, 1, 1, sa }, x), G({ a }, { sa }, x), Log(a, sa) };
    return -(sy[1] * G({ 0, x / a }, 1)) + 2. * sy[0] * G({ 0, 0, 1, x / a }, 1)
        + 2. * sy[0] * G({ 0, 0, x / a, 1 }, 1) + sy[0] * G({ 0, 1, 0, x / a }, 1)
        + 3. * sy[2] * G({ 0, 0, 0, 1, x / a }, 1)
        + 3. * sy[2] * G({ 0, 0, 0, x / a, 1 }, 1)
        + 3. * sy[2] * G({ 0, 0, 1, 0, x / a }, 1)
        + 2. * sy[2] * G({ 0, 1, 0, 0, x / a }, 1) + 6. * G({ 0, 0, 0, 1, 1, x / a }, 1)
        + 3. * G({ 0, 0, 0, 1, x / a, 1 }, 1) + 4. * G({ 0, 0, 1, 0, 1, x / a }, 1)
        + G({ 0, 0, 1, 0, x / a, 1 }, 1) + 2. * G({ 0, 0, 1, 1, 0, x / a }, 1)
        + 2. * G({ 0, 1, 0, 0, 1, x / a }, 1) + G({ 0, 1, 0, 1, 0, x / a }, 1)
        + G({ 0, x / a, 1 }, 1) * G({ 0, 0, a }, { 1, 1, sa }, x)
        + G({ 0, 0, a, 0, 0, a }, { 1, 1, sa, 1, 1, sa }, x) - sy[1] * Zeta(2)
        + G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) * Zeta(2) + (5. * sy[0] * Zeta(4)) / 4.
        + (5. * sy[3] * Log(-x, sa) * Zeta(4)) / 4. - (5. * pow(sy[3], 2.) * Zeta(4)) / 8.
        + (5. * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2)) * Zeta(4)) / 4.
        + 6. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        + 5. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        + 12. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.)
        + 3. * (-pow(Zeta(3), 2.) + (25. * Zeta(6)) / 12.)
        - 3. * sy[2] * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        - 6. * sy[2] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5))
        - 2. * sy[2] * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_00aa0a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return (3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a, 0, a }, { sa, 1, sa }, x),
        G({ 0, a }, { 1, sa }, x), G({ a }, { sa }, x), Log(a, sa), G({ 0 }, { 1 }, x),
        G({ 0, 0 }, { 1, 1 }, x), G({ 0, 0, a }, { 1, 1, sa }, x) };
    return sy[0] * G({ 0, 0, x / a }, 1) - sy[1] * G({ 0, 0, x / a, 1 }, 1)
        - 3. * sy[2] * G({ 0, 0, 0, 1, x / a }, 1)
        - 3. * sy[2] * G({ 0, 0, 0, x / a, 1 }, 1)
        - 2. * sy[2] * G({ 0, 0, 1, 0, x / a }, 1) - sy[2] * G({ 0, 1, 0, 0, x / a }, 1)
        - 6. * G({ 0, 0, 0, 1, 1, x / a }, 1) - 3. * G({ 0, 0, 0, 1, x / a, 1 }, 1)
        - 4. * G({ 0, 0, 1, 0, 1, x / a }, 1) - 2. * G({ 0, 0, 1, 0, x / a, 1 }, 1)
        - 4. * G({ 0, 0, 1, 1, 0, x / a }, 1) - 2. * G({ 0, 1, 0, 0, 1, x / a }, 1)
        - G({ 0, 1, 0, 0, x / a, 1 }, 1) - 2. * G({ 0, 1, 0, 1, 0, x / a }, 1)
        - 2. * G({ 0, 1, 1, 0, 0, x / a }, 1)
        + G({ 0, 0, 0, a, 0, a }, { 1, 1, 1, sa, 1, sa }, x) + pow(Zeta(3), 2.)
        + G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) * Zeta(2)
        - G({ 0, a, 0, a }, { 1, sa, 1, sa }, x) * Zeta(2) - (7. * sy[1] * Zeta(4)) / 4.
        - 3. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        - 2. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        - 8. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.) - (25. * Zeta(6)) / 12.
        + 2. * sy[4] * sy[5] * Zeta(3) + 2. * sy[6] * Zeta(3) - (-sy[0] + sy[6]) * Zeta(3)
        - 2. * G({ 0, 0, 0 }, { 1, 1, 1 }, x) * Zeta(3)
        - Log(-x, sa) * pow(sy[3], 2.) * Zeta(3) + (pow(sy[3], 3.) * Zeta(3)) / 3.
        + 4. * sy[4] * Zeta(2) * Zeta(3) - 2. * sy[3] * (sy[5] + 2. * Zeta(2)) * Zeta(3)
        + 2. * sy[2] * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        + 6. * sy[2] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5))
        + sy[2] * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_a0a00a(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ a, 0, 0, a }, { sa, 1, 1, sa }, x), G({ a }, { sa }, x) };
    return sy[0] * G({ 0, x / a }, 1) + sy[1] * G({ 0, 1, 0, 0, x / a }, 1)
        + G({ 0, 1, 0, 0, 1, x / a }, 1)
        + G({ 0, 1, 0, x / a }, 1) * G({ 0, a }, { 1, sa }, x)
        + G({ 0, 1, x / a }, 1) * G({ 0, 0, a }, { 1, 1, sa }, x)
        + G({ x / a }, 1) * G({ 0, a, 0, 0, a }, { 1, sa, 1, 1, sa }, x)
        + G({ 0, 0, a, 0, 0, a }, { 1, 1, sa, 1, 1, sa }, x)
        - (-sy[0] + G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)) * Zeta(2)
        + 2. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        + 3. * (-pow(Zeta(3), 2.) + (25. * Zeta(6)) / 12.)
        - sy[1] * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.)
        - Log(a, sa) * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.)
        + Log(-x, sa) * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_a00a0a(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ a, 0, a }, { sa, 1, sa }, x), G({ a }, { sa }, x) };
    return sy[0] * G({ 0, 0, x / a }, 1) + sy[1] * G({ 0, 0, 1, 0, x / a }, 1)
        + G({ 0, 0, 1, 0, 1, x / a }, 1)
        + G({ 0, 0, 1, x / a }, 1) * G({ 0, a }, { 1, sa }, x)
        + G({ 0, x / a }, 1) * G({ 0, a, 0, a }, { 1, sa, 1, sa }, x)
        + G({ x / a }, 1) * G({ 0, 0, a, 0, a }, { 1, 1, sa, 1, sa }, x)
        + G({ 0, 0, 0, a, 0, a }, { 1, 1, 1, sa, 1, sa }, x)
        + 3. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        + 2. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        - (-sy[0] + G({ 0, 0, a }, { 1, 1, sa }, x)) * Zeta(3)
        - sy[1] * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        - Log(a, sa) * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        + Log(-x, sa) * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_0a0a0a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return (-3. * Zeta(6)) / 16.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ a, 0, a }, { sa, 1, sa }, x), G({ 0, a }, { 1, sa }, x),
              G({ 0, a, 0, a }, { 1, sa, 1, sa }, x), G({ a }, { sa }, x), Log(a, sa) };
    return -(sy[2] * G({ 0, x / a }, 1)) - 2. * sy[0] * G({ 0, 0, x / a }, 1)
        - 2. * sy[1] * G({ 0, 0, 1, x / a }, 1) - sy[1] * G({ 0, 1, 0, x / a }, 1)
        - 2. * sy[3] * G({ 0, 0, 1, 0, x / a }, 1)
        - 2. * sy[3] * G({ 0, 1, 0, 0, x / a }, 1) - 2. * G({ 0, 0, 1, 0, 1, x / a }, 1)
        - 2. * G({ 0, 1, 0, 0, 1, x / a }, 1) - G({ 0, 1, 0, 1, 0, x / a }, 1)
        + G({ 0, 0, 0, a, 0, a }, { 1, 1, 1, sa, 1, sa }, x) + sy[2] * Zeta(2)
        - G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) * Zeta(2) - (3. * sy[1] * Zeta(4)) / 4.
        - (3. * G({ 0, 0 }, { 1, 1 }, x) * Zeta(4)) / 4.
        - (3. * sy[4] * Log(-x, sa) * Zeta(4)) / 4. + (3. * pow(sy[4], 2.) * Zeta(4)) / 8.
        - (3. * Zeta(2) * Zeta(4)) / 2. - 3. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        - 4. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        - 3. * (-pow(Zeta(3), 2.) + (25. * Zeta(6)) / 12.)
        + 2. * (-sy[0] + G({ 0, 0, a }, { 1, 1, sa }, x)) * Zeta(3)
        + 2. * sy[3] * (3. * Zeta(2) * Zeta(3) - (11. * Zeta(5)) / 2.)
        + 2. * sy[3] * (-2. * Zeta(2) * Zeta(3) + (9. * Zeta(5)) / 2.);
}

complex<double> G6_explicit_00aaaa(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return -pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ a, a, a }, { sa, sa, sa }, x), G({ a, a }, { sa, sa }, x),
              G({ a }, { sa }, x), G({ 0, a, a }, { 1, sa, sa }, x) };
    return sy[0] * G({ 0, 0, x / a }, 1) - sy[1] * G({ 0, 0, x / a, 1 }, 1)
        + sy[2] * G({ 0, 0, x / a, 1, 1 }, 1) - G({ 0, 0, x / a, 1, 1, 1 }, 1)
        + G({ 0, 0, 0, a, a, a }, { 1, 1, 1, sa, sa, sa }, x) + pow(Zeta(3), 2.)
        + G({ 0, 0, a, a }, { 1, 1, sa, sa }, x) * Zeta(2)
        - G({ 0, a, a, a }, { 1, sa, sa, sa }, x) * Zeta(2) + (sy[1] * Zeta(4)) / 4.
        - (5. * G({ 0, a }, { 1, sa }, x) * Zeta(4)) / 4. - (23. * Zeta(6)) / 16.
        - sy[3] * Zeta(3) - (-sy[0] + sy[3]) * Zeta(3)
        + G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(3)
        - sy[2] * (Zeta(2) * Zeta(3) - 2. * Zeta(5));
}
complex<double> G6_explicit_a00aaa(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a, a, a }, { sa, sa, sa }, x),
        G({ a, a }, { sa, sa }, x), G({ a }, { sa }, x) };
    return sy[0] * G({ 0, 0, x / a }, 1) + sy[1] * G({ 0, 0, 1, x / a }, 1)
        + sy[2] * G({ 0, 0, 1, 1, x / a }, 1) + G({ 0, 0, 1, 1, 1, x / a }, 1)
        + G({ 0, x / a }, 1) * G({ 0, a, a, a }, { 1, sa, sa, sa }, x)
        + G({ x / a }, 1) * G({ 0, 0, a, a, a }, { 1, 1, sa, sa, sa }, x)
        + G({ 0, 0, 0, a, a, a }, { 1, 1, 1, sa, sa, sa }, x)
        - (3. * pow(Zeta(3), 2.)) / 2. - (sy[1] * Zeta(4)) / 4.
        + (G({ 0, a }, { 1, sa }, x) * Zeta(4)) / 4.
        + 3. * (pow(Zeta(3), 2.) - (23. * Zeta(6)) / 16.) + (97. * Zeta(6)) / 48.
        - (-sy[0] + G({ 0, a, a }, { 1, sa, sa }, x)) * Zeta(3)
        - sy[2] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        - Log(a, sa) * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        + Log(-x, sa) * (Zeta(2) * Zeta(3) - 2. * Zeta(5));
}
complex<double> G6_explicit_aa00aa(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ 0, a, a }, { 1, sa, sa }, x),
        G({ a, a }, { sa, sa }, x), G({ a }, { sa }, x), Log(a, sa), Log(-x, sa) };
    return -(sy[0] * G({ 0, 1, x / a }, 1)) - sy[0] * G({ 0, x / a, 1 }, 1)
        - sy[1] * G({ 0, 0, 1, x / a }, 1) - sy[1] * G({ 0, 0, x / a, 1 }, 1)
        - sy[1] * G({ 0, 1, 0, x / a }, 1) - 2. * sy[2] * G({ 0, 0, 1, 1, x / a }, 1)
        - sy[2] * G({ 0, 0, 1, x / a, 1 }, 1) - sy[2] * G({ 0, 1, 0, 1, x / a }, 1)
        - 3. * G({ 0, 0, 1, 1, 1, x / a }, 1) - G({ 0, 0, 1, 1, x / a, 1 }, 1)
        - G({ 0, 1, 0, 1, 1, x / a }, 1)
        - G({ x / a, 1 }, 1) * G({ 0, 0, a, a }, { 1, 1, sa, sa }, x)
        + G({ x / a }, 1) * G({ a, 0, 0, a, a }, { sa, 1, 1, sa, sa }, x)
        + G({ 0, a, 0, 0, a, a }, { 1, sa, 1, 1, sa, sa }, x)
        + (5. * sy[1] * Zeta(4)) / 4. - (5. * G({ 0, a }, { 1, sa }, x) * Zeta(4)) / 4.
        - 4. * ((3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.)
        - 9. * (pow(Zeta(3), 2.) - (23. * Zeta(6)) / 16.) + (3. * Zeta(6)) / 16.
        - 5. * (-3. * pow(Zeta(3), 2.) + (203. * Zeta(6)) / 48.)
        + 3. * sy[2] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        + 3. * sy[3] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        - 3. * sy[4] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        + sy[2] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.)
        + sy[3] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.)
        - sy[4] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_aaa00a(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a }, { sa }, x), Log(a, sa), Log(-x, sa) };
    return sy[0] * G({ 0, 0, 1, 1, x / a }, 1) + sy[0] * G({ 0, 0, 1, x / a, 1 }, 1)
        + sy[0] * G({ 0, 0, x / a, 1, 1 }, 1) + sy[0] * G({ 0, 1, 0, 1, x / a }, 1)
        + sy[0] * G({ 0, 1, 0, x / a, 1 }, 1) + sy[0] * G({ 0, 1, 1, 0, x / a }, 1)
        + 3. * G({ 0, 0, 1, 1, 1, x / a }, 1) + 2. * G({ 0, 0, 1, 1, x / a, 1 }, 1)
        + G({ 0, 0, 1, x / a, 1, 1 }, 1) + 2. * G({ 0, 1, 0, 1, 1, x / a }, 1)
        + G({ 0, 1, 0, 1, x / a, 1 }, 1) + G({ 0, 1, 1, 0, 1, x / a }, 1)
        + (G({ 0, 1, 1, x / a }, 1) + G({ 0, 1, x / a, 1 }, 1) + G({ 0, x / a, 1, 1 }, 1))
        * G({ 0, a }, { 1, sa }, x)
        + G({ x / a, 1, 1 }, 1) * G({ 0, 0, a }, { 1, 1, sa }, x)
        - G({ x / a, 1 }, 1) * G({ a, 0, 0, a }, { sa, 1, 1, sa }, x)
        + G({ x / a }, 1) * G({ a, a, 0, 0, a }, { sa, sa, 1, 1, sa }, x)
        + G({ 0, a, a, 0, 0, a }, { 1, sa, sa, 1, 1, sa }, x)
        + 7. * ((3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.)
        + 9. * (pow(Zeta(3), 2.) - (23. * Zeta(6)) / 16.)
        + 2. * (-pow(Zeta(3), 2.) + (13. * Zeta(6)) / 16.) - (3. * Zeta(6)) / 8.
        + 6. * (-3. * pow(Zeta(3), 2.) + (203. * Zeta(6)) / 48.)
        - sy[0] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        - sy[1] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        + sy[2] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        - 3. * sy[0] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        - 3. * sy[1] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        + 3. * sy[2] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        - 2. * sy[0] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.)
        - 2. * sy[1] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.)
        + 2. * sy[2] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_0a0aaa(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a, a, a }, { sa, sa, sa }, x),
        G({ a, a }, { sa, sa }, x), G({ 0, a, a, a }, { 1, sa, sa, sa }, x),
        G({ a }, { sa }, x), Log(a, sa), G({ 0, a, a }, { 1, sa, sa }, x) };
    return -(sy[2] * G({ 0, x / a }, 1)) - 2. * sy[0] * G({ 0, 0, x / a }, 1)
        - 2. * sy[1] * G({ 0, 0, 1, x / a }, 1) - sy[1] * G({ 0, 1, 0, x / a }, 1)
        - 2. * sy[3] * G({ 0, 0, 1, 1, x / a }, 1) - sy[3] * G({ 0, 1, 0, 1, x / a }, 1)
        - sy[3] * G({ 0, 1, 1, 0, x / a }, 1) - 2. * G({ 0, 0, 1, 1, 1, x / a }, 1)
        - G({ 0, 1, 0, 1, 1, x / a }, 1) - G({ 0, 1, 1, 0, 1, x / a }, 1)
        - G({ 0, 1, 1, 1, 0, x / a }, 1)
        + G({ 0, 0, 0, a, a, a }, { 1, 1, 1, sa, sa, sa }, x) - pow(Zeta(3), 2.) / 2.
        + sy[2] * Zeta(2) - G({ 0, 0, a, a }, { 1, 1, sa, sa }, x) * Zeta(2)
        + (5. * sy[1] * Zeta(4)) / 4. + G({ 0, 0 }, { 1, 1 }, x) * Zeta(4)
        - (G({ 0, a }, { 1, sa }, x) * Zeta(4)) / 4. + sy[4] * Log(-x, sa) * Zeta(4)
        - (pow(sy[4], 2.) * Zeta(4)) / 2. + 2. * Zeta(2) * Zeta(4)
        - 2. * ((3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.)
        - 3. * (pow(Zeta(3), 2.) - (23. * Zeta(6)) / 16.) + (19. * Zeta(6)) / 12.
        - 2. * (-3. * pow(Zeta(3), 2.) + (203. * Zeta(6)) / 48.) - sy[5] * Zeta(3)
        + 2. * (-sy[0] + sy[5]) * Zeta(3) + G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(3)
        + sy[3] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        + 2. * sy[3] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        + sy[3] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_0aa0aa(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ 0, a, a }, { 1, sa, sa }, x), G({ a, a }, { sa, sa }, x),
              G({ a, 0, a, a }, { sa, 1, sa, sa }, x), G({ a }, { sa }, x), Log(a, sa) };
    return -(sy[2] * G({ 0, x / a }, 1)) + sy[0] * G({ 0, x / a, 1 }, 1)
        + 2. * sy[1] * G({ 0, 0, 1, x / a }, 1) + 2. * sy[1] * G({ 0, 0, x / a, 1 }, 1)
        + sy[1] * G({ 0, 1, 0, x / a }, 1) + 4. * sy[3] * G({ 0, 0, 1, 1, x / a }, 1)
        + 2. * sy[3] * G({ 0, 0, 1, x / a, 1 }, 1)
        + 2. * sy[3] * G({ 0, 1, 0, 1, x / a }, 1) + sy[3] * G({ 0, 1, 0, x / a, 1 }, 1)
        + 2. * sy[3] * G({ 0, 1, 1, 0, x / a }, 1) + 6. * G({ 0, 0, 1, 1, 1, x / a }, 1)
        + 2. * G({ 0, 0, 1, 1, x / a, 1 }, 1) + 3. * G({ 0, 1, 0, 1, 1, x / a }, 1)
        + G({ 0, 1, 0, 1, x / a, 1 }, 1) + 3. * G({ 0, 1, 1, 0, 1, x / a }, 1)
        + G({ 0, 1, 1, 0, x / a, 1 }, 1) + 3. * G({ 0, 1, 1, 1, 0, x / a }, 1)
        + G({ 0, 0, a, 0, a, a }, { 1, 1, sa, 1, sa, sa }, x) - sy[2] * Zeta(2)
        + G({ 0, 0, a, a }, { 1, 1, sa, sa }, x) * Zeta(2) - (7. * sy[1] * Zeta(4)) / 4.
        - (5. * G({ 0, a }, { 1, sa }, x) * Zeta(4)) / 4.
        - 3. * sy[4] * Log(-x, sa) * Zeta(4) + (3. * pow(sy[4], 2.) * Zeta(4)) / 2.
        - 3. * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2)) * Zeta(4)
        + 7. * ((3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.)
        + 9. * (pow(Zeta(3), 2.) - (23. * Zeta(6)) / 16.)
        + 2. * (-pow(Zeta(3), 2.) + (13. * Zeta(6)) / 16.) - (3. * Zeta(6)) / 8.
        + 6. * (-3. * pow(Zeta(3), 2.) + (203. * Zeta(6)) / 48.) + 2. * sy[0] * Zeta(3)
        - 2. * G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(3)
        - 2. * sy[3] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        - 6. * sy[3] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        - 3. * sy[3] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_0aaa0a(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return -pow(Zeta(3), 2.) + (25. * Zeta(6)) / 12.;
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ a, 0, a }, { sa, 1, sa }, x), G({ 0, a }, { 1, sa }, x),
              G({ a, a, 0, a }, { sa, sa, 1, sa }, x), G({ a }, { sa }, x), Log(a, sa) };
    return -(sy[2] * G({ 0, x / a }, 1)) + sy[0] * G({ 0, x / a, 1 }, 1)
        - sy[1] * G({ 0, x / a, 1, 1 }, 1) - 2. * sy[3] * G({ 0, 0, 1, 1, x / a }, 1)
        - 2. * sy[3] * G({ 0, 0, 1, x / a, 1 }, 1)
        - 2. * sy[3] * G({ 0, 0, x / a, 1, 1 }, 1) - sy[3] * G({ 0, 1, 0, 1, x / a }, 1)
        - sy[3] * G({ 0, 1, 0, x / a, 1 }, 1) - sy[3] * G({ 0, 1, 1, 0, x / a }, 1)
        - 6. * G({ 0, 0, 1, 1, 1, x / a }, 1) - 4. * G({ 0, 0, 1, 1, x / a, 1 }, 1)
        - 2. * G({ 0, 0, 1, x / a, 1, 1 }, 1) - 3. * G({ 0, 1, 0, 1, 1, x / a }, 1)
        - 2. * G({ 0, 1, 0, 1, x / a, 1 }, 1) - G({ 0, 1, 0, x / a, 1, 1 }, 1)
        - 3. * G({ 0, 1, 1, 0, 1, x / a }, 1) - 2. * G({ 0, 1, 1, 0, x / a, 1 }, 1)
        - 3. * G({ 0, 1, 1, 1, 0, x / a }, 1)
        + G({ 0, 0, a, a, 0, a }, { 1, 1, sa, sa, 1, sa }, x) + pow(Zeta(3), 2.)
        - sy[2] * Zeta(2) + G({ 0, a, 0, a }, { 1, sa, 1, sa }, x) * Zeta(2)
        + 3. * sy[1] * Zeta(4) + 3. * sy[4] * Log(-x, sa) * Zeta(4)
        - (3. * pow(sy[4], 2.) * Zeta(4)) / 2.
        + 3. * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2)) * Zeta(4)
        - 4. * ((3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.)
        - 9. * (pow(Zeta(3), 2.) - (23. * Zeta(6)) / 16.) - (5. * Zeta(6)) / 8.
        - 4. * (-3. * pow(Zeta(3), 2.) + (203. * Zeta(6)) / 48.) - sy[0] * Zeta(3)
        + G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(3)
        + sy[3] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        + 6. * sy[3] * (Zeta(2) * Zeta(3) - 2. * Zeta(5))
        + 2. * sy[3] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_a0a0aa(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a, a }, { sa, sa }, x),
        G({ a, 0, a, a }, { sa, 1, sa, sa }, x), G({ a }, { sa }, x) };
    return sy[1] * G({ 0, x / a }, 1) + sy[0] * G({ 0, 1, 0, x / a }, 1)
        + sy[2] * G({ 0, 1, 0, 1, x / a }, 1) + G({ 0, 1, 0, 1, 1, x / a }, 1)
        + G({ 0, 1, x / a }, 1) * G({ 0, a, a }, { 1, sa, sa }, x)
        + G({ x / a }, 1) * G({ 0, a, 0, a, a }, { 1, sa, 1, sa, sa }, x)
        + G({ 0, 0, a, 0, a, a }, { 1, 1, sa, 1, sa, sa }, x)
        - (-sy[1] + G({ 0, 0, a, a }, { 1, 1, sa, sa }, x)) * Zeta(2)
        - (3. * sy[0] * Zeta(4)) / 4. + (3. * G({ 0, a }, { 1, sa }, x) * Zeta(4)) / 4.
        + 2. * ((3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.)
        - (3. * Zeta(6)) / 16. + 2. * (-3. * pow(Zeta(3), 2.) + (203. * Zeta(6)) / 48.)
        - sy[2] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.)
        - Log(a, sa) * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.)
        + Log(-x, sa) * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_a0aa0a(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a, 0, a }, { sa, 1, sa }, x),
        G({ a, a, 0, a }, { sa, sa, 1, sa }, x), G({ a }, { sa }, x) };
    return sy[1] * G({ 0, x / a }, 1) + sy[0] * G({ 0, 1, x / a }, 1)
        + sy[2] * G({ 0, 1, 1, 0, x / a }, 1) + G({ 0, 1, 1, 0, 1, x / a }, 1)
        + G({ 0, 1, 1, x / a }, 1) * G({ 0, a }, { 1, sa }, x)
        + G({ x / a }, 1) * G({ 0, a, a, 0, a }, { 1, sa, sa, 1, sa }, x)
        + G({ 0, 0, a, a, 0, a }, { 1, 1, sa, sa, 1, sa }, x) + sy[1] * Zeta(2)
        - G({ 0, a, 0, a }, { 1, sa, 1, sa }, x) * Zeta(2)
        + 2. * ((3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.)
        + 2. * (-pow(Zeta(3), 2.) + (13. * Zeta(6)) / 16.) - (3. * Zeta(6)) / 16.
        - sy[0] * Zeta(3) + G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(3)
        - sy[2] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        - Log(a, sa) * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        + Log(-x, sa) * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.);
}
complex<double> G6_explicit_aa0a0a(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a, 0, a }, { sa, 1, sa }, x),
        G({ 0, a }, { 1, sa }, x), G({ a }, { sa }, x), Log(a, sa), Log(-x, sa) };

    return -(sy[0] * G({ 0, 1, x / a }, 1)) - sy[0] * G({ 0, x / a, 1 }, 1)
        - 2. * sy[1] * G({ 0, 1, 1, x / a }, 1) - sy[1] * G({ 0, 1, x / a, 1 }, 1)
        - sy[2] * G({ 0, 1, 0, 1, x / a }, 1) - sy[2] * G({ 0, 1, 0, x / a, 1 }, 1)
        - 2. * sy[2] * G({ 0, 1, 1, 0, x / a }, 1) - 2. * G({ 0, 1, 0, 1, 1, x / a }, 1)
        - G({ 0, 1, 0, 1, x / a, 1 }, 1) - 2. * G({ 0, 1, 1, 0, 1, x / a }, 1)
        - G({ x / a, 1 }, 1) * G({ 0, a, 0, a }, { 1, sa, 1, sa }, x)
        + G({ x / a }, 1) * G({ a, 0, a, 0, a }, { sa, 1, sa, 1, sa }, x)
        + G({ 0, a, 0, a, 0, a }, { 1, sa, 1, sa, 1, sa }, x)
        - 8. * ((3. * pow(Zeta(3), 2.)) / 2. - (53. * Zeta(6)) / 24.)
        - 4. * (-pow(Zeta(3), 2.) + (13. * Zeta(6)) / 16.) + (9. * Zeta(6)) / 16.
        - 4. * (-3. * pow(Zeta(3), 2.) + (203. * Zeta(6)) / 48.) + 2. * sy[0] * Zeta(3)
        - 2. * G({ 0, 0, a }, { 1, 1, sa }, x) * Zeta(3)
        + 2. * sy[2] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        + 2. * sy[3] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        - 2. * sy[4] * (2. * Zeta(2) * Zeta(3) - (9. * Zeta(5)) / 2.)
        + 2. * sy[2] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.)
        + 2. * sy[3] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.)
        - 2. * sy[4] * (-3. * Zeta(2) * Zeta(3) + (11. * Zeta(5)) / 2.);
}

complex<double> G6_explicit_0aaaaa(complex<double> a1, int sa, double x1)
{
    if (x1 == a1) {
        return -Zeta(6);
    }

    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ a, a, a }, { sa, sa, sa }, x), G({ a, a }, { sa, sa }, x),
              G({ a, a, a, a }, { sa, sa, sa, sa }, x), G({ a }, { sa }, x) };
    return -(sy[2] * G({ 0, x / a }, 1)) + sy[0] * G({ 0, x / a, 1 }, 1)
        - sy[1] * G({ 0, x / a, 1, 1 }, 1) + sy[3] * G({ 0, x / a, 1, 1, 1 }, 1)
        - G({ 0, x / a, 1, 1, 1, 1 }, 1)
        + G({ 0, 0, a, a, a, a }, { 1, 1, sa, sa, sa, sa }, x) - pow(Zeta(3), 2.) / 2.
        - sy[2] * Zeta(2) + G({ 0, a, a, a }, { 1, sa, sa, sa }, x) * Zeta(2)
        - sy[1] * Zeta(4) + G({ 0, a }, { 1, sa }, x) * Zeta(4) + (3. * Zeta(6)) / 4.
        + (-sy[0] + G({ 0, a, a }, { 1, sa, sa }, x)) * Zeta(3) - sy[3] * Zeta(5);
}
complex<double> G6_explicit_a0aaaa(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ a, a, a }, { sa, sa, sa }, x), G({ a, a }, { sa, sa }, x),
              G({ a, a, a, a }, { sa, sa, sa, sa }, x), G({ a }, { sa }, x) };

    return sy[2] * G({ 0, x / a }, 1) + sy[0] * G({ 0, 1, x / a }, 1)
        + sy[1] * G({ 0, 1, 1, x / a }, 1) + sy[3] * G({ 0, 1, 1, 1, x / a }, 1)
        + G({ 0, 1, 1, 1, 1, x / a }, 1)
        + G({ x / a }, 1) * G({ 0, a, a, a, a }, { 1, sa, sa, sa, sa }, x)
        + G({ 0, 0, a, a, a, a }, { 1, 1, sa, sa, sa, sa }, x) + pow(Zeta(3), 2.) / 2.
        + sy[2] * Zeta(2) - G({ 0, a, a, a }, { 1, sa, sa, sa }, x) * Zeta(2)
        + sy[1] * Zeta(4) - G({ 0, a }, { 1, sa }, x) * Zeta(4)
        + 2. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.) + Zeta(6) / 4.
        + (-sy[0] + G({ 0, a, a }, { 1, sa, sa }, x)) * Zeta(3) - sy[3] * Zeta(5)
        - Log(a, sa) * Zeta(5) + Log(-x, sa) * Zeta(5);
}
complex<double> G6_explicit_aa0aaa(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a, a, a }, { sa, sa, sa }, x),
        G({ a, a }, { sa, sa }, x), G({ a }, { sa }, x) };

    return -(sy[0] * G({ 0, 1, x / a }, 1)) - sy[0] * G({ 0, x / a, 1 }, 1)
        - 2. * sy[1] * G({ 0, 1, 1, x / a }, 1) - sy[1] * G({ 0, 1, x / a, 1 }, 1)
        - 3. * sy[2] * G({ 0, 1, 1, 1, x / a }, 1) - sy[2] * G({ 0, 1, 1, x / a, 1 }, 1)
        - 4. * G({ 0, 1, 1, 1, 1, x / a }, 1) - G({ 0, 1, 1, 1, x / a, 1 }, 1)
        - G({ x / a, 1 }, 1) * G({ 0, a, a, a }, { 1, sa, sa, sa }, x)
        + G({ x / a }, 1) * G({ a, 0, a, a, a }, { sa, 1, sa, sa, sa }, x)
        + G({ 0, a, 0, a, a, a }, { 1, sa, 1, sa, sa, sa }, x) - 3. * sy[1] * Zeta(4)
        + 3. * G({ 0, a }, { 1, sa }, x) * Zeta(4)
        - 4. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        - 4. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        - 8. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.)
        - 3. * (-pow(Zeta(3), 2.) + (25. * Zeta(6)) / 12.)
        - 2. * (-sy[0] + G({ 0, a, a }, { 1, sa, sa }, x)) * Zeta(3)
        + 4. * sy[2] * Zeta(5) + 4. * Log(a, sa) * Zeta(5) - 4. * Log(-x, sa) * Zeta(5);
}
complex<double> G6_explicit_aaa0aa(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy
        = { G({ a, a }, { sa, sa }, x), G({ a }, { sa }, x) };

    return sy[0] * G({ 0, 1, 1, x / a }, 1) + sy[0] * G({ 0, 1, x / a, 1 }, 1)
        + sy[0] * G({ 0, x / a, 1, 1 }, 1) + 3. * sy[1] * G({ 0, 1, 1, 1, x / a }, 1)
        + 2. * sy[1] * G({ 0, 1, 1, x / a, 1 }, 1) + sy[1] * G({ 0, 1, x / a, 1, 1 }, 1)
        + 6. * G({ 0, 1, 1, 1, 1, x / a }, 1) + 3. * G({ 0, 1, 1, 1, x / a, 1 }, 1)
        + G({ 0, 1, 1, x / a, 1, 1 }, 1)
        + G({ x / a, 1, 1 }, 1) * G({ 0, a, a }, { 1, sa, sa }, x)
        - G({ x / a, 1 }, 1) * G({ a, 0, a, a }, { sa, 1, sa, sa }, x)
        + G({ x / a }, 1) * G({ a, a, 0, a, a }, { sa, sa, 1, sa, sa }, x)
        + G({ 0, a, a, 0, a, a }, { 1, sa, sa, 1, sa, sa }, x) + 3. * sy[0] * Zeta(4)
        - 3. * G({ 0, a }, { 1, sa }, x) * Zeta(4)
        + 6. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        + 5. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        + 12. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.)
        + 3. * (-pow(Zeta(3), 2.) + (25. * Zeta(6)) / 12.) - 6. * sy[1] * Zeta(5)
        - 6. * Log(a, sa) * Zeta(5) + 6. * Log(-x, sa) * Zeta(5);
}
complex<double> G6_explicit_aaaa0a(complex<double> a1, int sa, double x1)
{
    complex<double> a = a1 / x1;
    double x = 1.0;
    const vector<complex<double>> sy = { G({ a }, { sa }, x) };
    return -(sy[0] * G({ 0, 1, 1, 1, x / a }, 1)) - sy[0] * G({ 0, 1, 1, x / a, 1 }, 1)
        - sy[0] * G({ 0, 1, x / a, 1, 1 }, 1) - sy[0] * G({ 0, x / a, 1, 1, 1 }, 1)
        - 4. * G({ 0, 1, 1, 1, 1, x / a }, 1) - 3. * G({ 0, 1, 1, 1, x / a, 1 }, 1)
        - 2. * G({ 0, 1, 1, x / a, 1, 1 }, 1) - G({ 0, 1, x / a, 1, 1, 1 }, 1)
        - G({ x / a, 1, 1, 1 }, 1) * G({ 0, a }, { 1, sa }, x)
        + G({ x / a, 1, 1 }, 1) * G({ a, 0, a }, { sa, 1, sa }, x)
        - G({ x / a, 1 }, 1) * G({ a, a, 0, a }, { sa, sa, 1, sa }, x)
        + G({ x / a }, 1) * G({ a, a, a, 0, a }, { sa, sa, sa, 1, sa }, x)
        + G({ 0, a, a, a, 0, a }, { 1, sa, sa, sa, 1, sa }, x) + pow(Zeta(3), 2.)
        - 3. * (pow(Zeta(3), 2.) - (4. * Zeta(6)) / 3.)
        - 2. * (pow(Zeta(3), 2.) / 2. - Zeta(6) / 2.)
        - 8. * (-pow(Zeta(3), 2.) / 2. + (3. * Zeta(6)) / 4.) - (25. * Zeta(6)) / 12.
        + 4. * sy[0] * Zeta(5) + 4. * Log(a, sa) * Zeta(5) - 4. * Log(-x, sa) * Zeta(5);
}

complex<double> G6_0000ab(
    complex<double> a1, complex<double> b1, int sa, int sb, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    double x = 1.;
    if (abs(b) < abs(a)) { // b smallest
        const vector<complex<double>> sy = { Log(b, sb), G({ a / b }, 1),
            G({ 0, a / b }, 1), G({ 0, 0 }, { 1, 1 }, x), G({ 0, 0, 0 }, { 1, 1, 1 }, x),
            G({ 0, 0, a / b }, 1), G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x),
            G({ 0, 0, 0, a / b }, 1), G({ 0, 0, 0, 0, a / b }, 1), G({ 0 }, { 1 }, x),
            G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x) };
        complex<double> res { -(sy[4] * sy[5]) - G({ 0, 0, 0, 0, 0, a / b }, 1)
            - G({ a / b, 0, 0, 0, 0, x / b }, 1)
            + sy[6] * (-G({ x / b }, 1) + G({ b }, { sb }, x))
            - 5. * G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x)
            - (sy[2] * pow(sy[0], 4.)) / 24.
            + Log(-x, sb)
                * (sy[0] * sy[7] - sy[8] - (sy[5] * pow(sy[0], 2.)) / 2.
                    + (sy[2] * pow(sy[0], 3.)) / 6. - (sy[1] * pow(sy[0], 4.)) / 24.)
            + (sy[1] * pow(sy[0], 5.)) / 120. - sy[2] * pow(sy[3], 2.)
            + pow(sy[0], 3.) * (sy[5] / 6. + sy[1] * (-sy[3] / 6. - Zeta(2) / 3.))
            + 2. * sy[7] * Zeta(2) + 2. * sy[2] * pow(sy[9], 2.) * Zeta(2)
            + sy[3]
                * (sy[9] * sy[5] + sy[7] + sy[2] * pow(sy[9], 2.) - 2. * sy[2] * Zeta(2))
            + sy[9] * (-(sy[2] * sy[4]) + 2. * sy[5] * Zeta(2))
            + pow(sy[0], 2.)
                * ((sy[2] * sy[3]) / 2. - sy[7] / 2. + sy[2] * Zeta(2)
                    + sy[1] * ((sy[9] * sy[3]) / 2. - sy[4] / 2. + sy[9] * Zeta(2)))
            + sy[0]
                * (sy[2] * sy[4] + sy[3] * (-(sy[9] * sy[2]) - sy[5]) + sy[8]
                    - 2. * sy[9] * sy[2] * Zeta(2) - 2. * sy[5] * Zeta(2)
                    + sy[1]
                        * (-sy[10] + sy[9] * sy[4] + pow(sy[3], 2.)
                            - 2. * pow(sy[9], 2.) * Zeta(2)
                            + sy[3] * (-pow(sy[9], 2.) + 2. * Zeta(2)) - 2. * Zeta(4)))
            + sy[2] * (sy[10] + 2. * Zeta(4))
            + sy[1]
                * (sy[6] - G({ 0, 0, 0, 0, 0 }, { 1, 1, 1, 1, 1 }, x)
                    - sy[4] * pow(sy[9], 2.) - 2. * sy[9] * pow(sy[3], 2.)
                    + 2. * sy[4] * Zeta(2) + 2. * pow(sy[9], 3.) * Zeta(2)
                    + sy[3] * (2. * sy[4] + pow(sy[9], 3.) - 4. * sy[9] * Zeta(2))
                    + sy[9] * (sy[10] + 2. * Zeta(4))) };
        return res;
    }

    // a smallest
    const vector<complex<double>> sy
        = { Log(a, sa), G({ b / a }, 1), G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x),
              G({ 0, 0, 0 }, { 1, 1, 1 }, x), G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x) };
    complex<double> res { G({ 0, 0, 0, 0, 0, b / a }, 1)
        + G({ 0, 0, 0, 0, b / a, x / a }, 1) + G({ 0, 0, 0, b / a, 0, x / a }, 1)
        + G({ 0, 0, b / a, 0, 0, x / a }, 1) + G({ 0, b / a, 0, 0, 0, x / a }, 1)
        + G({ b / a, 0, 0, 0, 0, x / a }, 1)
        + G({ 0, 0, 0, b / a }, 1) * G({ 0, b }, { 1, sb }, x)
        - G({ 0, 0, b / a }, 1) * G({ 0, 0, b }, { 1, 1, sb }, x)
        + G({ 0, b / a }, 1) * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
        + G({ 0, 0, 0, 0, 0, b }, { 1, 1, 1, 1, 1, sb }, x)
        + (sy[1] * Log(-x, sa) * pow(sy[0], 4.)) / 24. - (sy[1] * pow(sy[0], 5.)) / 120.
        + sy[1] * pow(sy[0], 3.) * (sy[2] / 6. + Zeta(2) / 3.)
        + sy[1] * pow(sy[0], 2.) * (-(sy[2] * sy[3]) / 2. + sy[4] / 2. - sy[3] * Zeta(2))
        + sy[1]
            * (G({ 0, 0, 0, 0, 0 }, { 1, 1, 1, 1, 1 }, x)
                - G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x)
                + 2. * sy[3] * pow(sy[2], 2.) + sy[4] * pow(sy[3], 2.)
                - 2. * sy[4] * Zeta(2) - 2. * pow(sy[3], 3.) * Zeta(2)
                + sy[2] * (-2. * sy[4] - pow(sy[3], 3.) + 4. * sy[3] * Zeta(2))
                + sy[3] * (-sy[5] - 2. * Zeta(4)))
        + sy[0] * sy[1]
            * (-(sy[3] * sy[4]) + sy[5] - pow(sy[2], 2.)
                + sy[2] * (pow(sy[3], 2.) - 2. * Zeta(2)) + 2. * pow(sy[3], 2.) * Zeta(2)
                + 2. * Zeta(4)) };

    if (b != x) {
        res += (-G({ 0, 0, 0, 0, b / a }, 1) + G({ 0, 0, 0, 0, x / a }, 1))
            * G({ b }, { sb }, x);
    }
    return res;
}
complex<double> G6_a0000b(
    complex<double> a1, complex<double> b1, int sa, int sb, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    double x = 1.;
    if (abs(b) < abs(a)) { // b smallest
        const vector<complex<double>> sy = { G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x),
            G({ 0, a }, { 1, sa }, x), G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x),
            G({ 0, 0, 0, 0, a / b }, 1), G({ a }, { sa }, x) };
        complex<double> res { sy[3] * sy[4] - sy[2] * G({ x / b }, 1)
            - sy[0] * G({ 0, x / b }, 1) - sy[1] * G({ 0, 0, 0, x / b }, 1)
            - sy[4] * G({ 0, 0, 0, 0, x / b }, 1) - 5. * G({ 0, 0, 0, 0, 0, a / b }, 1)
            - G({ 0, 0, 0, 0, a / b, x / b }, 1) + sy[2] * G({ b }, { sb }, x)
            - sy[0] * G({ 0, b }, { 1, sb }, x)
            + G({ 0, 0, a }, { 1, 1, sa }, x)
                * (-G({ 0, 0, x / b }, 1) + G({ 0, 0, b }, { 1, 1, sb }, x))
            - sy[1] * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
            + sy[4] * G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x)
            - G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x) + sy[3] * Log(b, sb)
            - sy[3] * Log(-x, sb) };
        return res;
    }

    // a smallest
    const vector<complex<double>> sy = { G({ 0, 0, 0, 0, b / a }, 1) };
    complex<double> res { 5. * G({ 0, 0, 0, 0, 0, b / a }, 1)
        + G({ 0, 0, 0, 0, b / a, x / a }, 1)
        + G({ 0, 0, 0, x / a }, 1) * G({ 0, b }, { 1, sb }, x)
        + G({ 0, 0, x / a }, 1) * G({ 0, 0, b }, { 1, 1, sb }, x)
        + G({ 0, x / a }, 1) * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
        + G({ x / a }, 1) * G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x)
        + G({ 0, 0, 0, 0, 0, b }, { 1, 1, 1, 1, 1, sb }, x) - sy[0] * Log(a, sa)
        + sy[0] * Log(-x, sa) };

    if (b != x) {
        res += (-sy[0] + G({ 0, 0, 0, 0, x / a }, 1)) * G({ b }, { sb }, x);
    }
    return res;
}
complex<double> G6_000a0b(
    complex<double> a1, complex<double> b1, int sa, int sb, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    double x = 1.;
    if (abs(b) < abs(a)) { // b smallest
        const vector<complex<double>> sy = { Log(b, sb), G({ 0, a / b }, 1),
            G({ 0, 0, a / b }, 1), G({ 0, 0, 0 }, { 1, 1, 1 }, x),
            G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x), G({ 0, 0 }, { 1, 1 }, x),
            G({ 0 }, { 1 }, x), G({ 0, 0, 0, a / b }, 1),
            G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x), G({ 0, 0, 0, 0, a / b }, 1) };
        complex<double> res { 2. * sy[2] * sy[3]
            + sy[5] * (-2. * sy[2] * sy[6] - 3. * sy[7]) + 4. * sy[8] * G({ x / b }, 1)
            + 5. * G({ 0, 0, 0, 0, 0, a / b }, 1) - G({ 0, a / b, 0, 0, 0, x / b }, 1)
            - 4. * sy[8] * G({ b }, { sb }, x)
            + sy[4] * (G({ 0, x / b }, 1) + G({ 0, b }, { 1, sb }, x))
            + 10. * G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x)
            - (sy[2] * pow(sy[0], 3.)) / 3.
            + Log(-x, sb)
                * (4. * sy[9] - 3. * sy[0] * sy[7] + sy[2] * pow(sy[0], 2.)
                    - (sy[1] * pow(sy[0], 3.)) / 6.)
            + (sy[1] * pow(sy[0], 4.)) / 24.
            + pow(sy[0], 2.) * ((3. * sy[7]) / 2. + sy[1] * (-sy[5] / 2. - Zeta(2)))
            - 4. * sy[2] * sy[6] * Zeta(2) - 6. * sy[7] * Zeta(2)
            + sy[0]
                * (-4. * sy[9] + 2. * sy[2] * sy[5] + 4. * sy[2] * Zeta(2)
                    + sy[1] * (-sy[3] + sy[5] * sy[6] + 2. * sy[6] * Zeta(2)))
            + sy[1]
                * (-sy[4] + sy[3] * sy[6] - G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x)
                    + pow(sy[5], 2.) - 2. * pow(sy[6], 2.) * Zeta(2)
                    + sy[5] * (-pow(sy[6], 2.) + 2. * Zeta(2)) - 2. * Zeta(4)) };
        return res;
    }

    // a smallest
    const vector<complex<double>> sy = { Log(a, sa), G({ 0, b / a }, 1),
        G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x), G({ 0, 0, 0 }, { 1, 1, 1 }, x) };
    complex<double> res { -5. * G({ 0, 0, 0, 0, 0, b / a }, 1)
        - 4. * G({ 0, 0, 0, 0, b / a, x / a }, 1)
        - 3. * G({ 0, 0, 0, b / a, 0, x / a }, 1)
        - 2. * G({ 0, 0, b / a, 0, 0, x / a }, 1) - G({ 0, b / a, 0, 0, 0, x / a }, 1)
        + (-3. * G({ 0, 0, 0, b / a }, 1) - G({ 0, 0, 0, x / a }, 1))
            * G({ 0, b }, { 1, sb }, x)
        + 2. * G({ 0, 0, b / a }, 1) * G({ 0, 0, b }, { 1, 1, sb }, x)
        + G({ 0, 0, 0, 0, 0, b }, { 1, 1, 1, 1, 1, sb }, x)
        - (sy[1] * Log(-x, sa) * pow(sy[0], 3.)) / 6. + (sy[1] * pow(sy[0], 4.)) / 24.
        + sy[1] * pow(sy[0], 2.) * (-sy[2] / 2. - Zeta(2))
        + sy[0] * sy[1] * (sy[2] * sy[3] - sy[4] + 2. * sy[3] * Zeta(2))
        + sy[1]
            * (sy[3] * sy[4] - G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x)
                - G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x) + pow(sy[2], 2.)
                - 2. * pow(sy[3], 2.) * Zeta(2) + sy[2] * (-pow(sy[3], 2.) + 2. * Zeta(2))
                - 2. * Zeta(4)) };

    if (b != x) {
        res += (4. * G({ 0, 0, 0, 0, b / a }, 1) - 4. * G({ 0, 0, 0, 0, x / a }, 1))
            * G({ b }, { sb }, x);
    }
    return res;
}
complex<double> G6_0a000b(
    complex<double> a1, complex<double> b1, int sa, int sb, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    double x = 1.;
    if (abs(b) < abs(a)) { // b smallest
        const vector<complex<double>> sy = { G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x),
            Log(b, sb), G({ 0, 0, 0, a / b }, 1), G({ 0, a }, { 1, sa }, x),
            G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x), G({ 0, 0, 0, 0, a / b }, 1) };
        complex<double> res { -4. * sy[1] * sy[5] + 4. * sy[4] * G({ x / b }, 1)
            + 3. * sy[0] * G({ 0, x / b }, 1) + 10. * G({ 0, 0, 0, 0, 0, a / b }, 1)
            - G({ 0, 0, 0, a / b, 0, x / b }, 1) - 4. * sy[4] * G({ b }, { sb }, x)
            + 3. * sy[0] * G({ 0, b }, { 1, sb }, x)
            + G({ 0, 0, a }, { 1, 1, sa }, x)
                * (2. * G({ 0, 0, x / b }, 1) - 2. * G({ 0, 0, b }, { 1, 1, sb }, x))
            + sy[3] * (G({ 0, 0, 0, x / b }, 1) + G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x))
            + 5. * G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x)
            + (-(sy[1] * sy[2]) + 4. * sy[5]) * Log(-x, sb)
            + (sy[2] * pow(sy[1], 2.)) / 2.
            + sy[2] * (-sy[3] - G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
        return res;
    }

    // a smallest
    const vector<complex<double>> sy
        = { Log(a, sa), G({ 0, 0, 0, b / a }, 1), G({ 0, b }, { 1, sb }, x) };
    complex<double> res { -3. * sy[2] * G({ 0, 0, 0, x / a }, 1)
        - 10. * G({ 0, 0, 0, 0, 0, b / a }, 1) - 4. * G({ 0, 0, 0, 0, b / a, x / a }, 1)
        - G({ 0, 0, 0, b / a, 0, x / a }, 1)
        - 2. * G({ 0, 0, x / a }, 1) * G({ 0, 0, b }, { 1, 1, sb }, x)
        - G({ 0, x / a }, 1) * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
        + G({ 0, 0, 0, 0, 0, b }, { 1, 1, 1, 1, 1, sb }, x) - sy[0] * sy[1] * Log(-x, sa)
        + (sy[1] * pow(sy[0], 2.)) / 2.
        + sy[1] * (-sy[2] - G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };

    if (b != x) {
        res += (4. * G({ 0, 0, 0, 0, b / a }, 1) - 4. * G({ 0, 0, 0, 0, x / a }, 1))
            * G({ b }, { sb }, x);
    }
    return res;
}
complex<double> G6_00a00b(
    complex<double> a1, complex<double> b1, int sa, int sb, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    double x = 1.;
    if (abs(b) < abs(a)) { // b smallest
        const vector<complex<double>> sy
            = { Log(b, sb), G({ 0, 0, a / b }, 1), G({ 0, 0, a }, { 1, 1, sa }, x),
                  G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x), G({ 0, 0, 0, a / b }, 1),
                  G({ 0, 0 }, { 1, 1 }, x), G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x),
                  G({ 0, 0, 0, 0, a / b }, 1), G({ 0 }, { 1 }, x) };
        complex<double> res { 3. * sy[4] * sy[5] - 6. * sy[6] * G({ x / b }, 1)
            - 3. * sy[3] * G({ 0, x / b }, 1) - 10. * G({ 0, 0, 0, 0, 0, a / b }, 1)
            - G({ 0, 0, a / b, 0, 0, x / b }, 1) + 6. * sy[6] * G({ b }, { sb }, x)
            - 3. * sy[3] * G({ 0, b }, { 1, sb }, x)
            + sy[2] * (-G({ 0, 0, x / b }, 1) + G({ 0, 0, b }, { 1, 1, sb }, x))
            - 10. * G({ 0, 0, 0, 0, 0, a }, { 1, 1, 1, 1, 1, sa }, x)
            - (3. * sy[4] * pow(sy[0], 2.)) / 2.
            + Log(-x, sb)
                * (3. * sy[0] * sy[4] - 6. * sy[7] - (sy[1] * pow(sy[0], 2.)) / 2.)
            + (sy[1] * pow(sy[0], 3.)) / 6.
            + sy[0] * (6. * sy[7] + sy[1] * (-sy[5] - 2. * Zeta(2)))
            + 6. * sy[4] * Zeta(2)
            + sy[1]
                * (sy[2] + sy[5] * sy[8] - G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                    + 2. * sy[8] * Zeta(2)) };
        return res;
    }

    // a smallest
    const vector<complex<double>> sy = { Log(a, sa), G({ 0, 0, b / a }, 1),
        G({ 0, 0, b }, { 1, 1, sb }, x), G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x) };
    complex<double> res { sy[2] * G({ 0, 0, x / a }, 1)
        + 10. * G({ 0, 0, 0, 0, 0, b / a }, 1) + 6. * G({ 0, 0, 0, 0, b / a, x / a }, 1)
        + 3. * G({ 0, 0, 0, b / a, 0, x / a }, 1) + G({ 0, 0, b / a, 0, 0, x / a }, 1)
        + (3. * G({ 0, 0, 0, b / a }, 1) + 3. * G({ 0, 0, 0, x / a }, 1))
            * G({ 0, b }, { 1, sb }, x)
        + G({ 0, 0, 0, 0, 0, b }, { 1, 1, 1, 1, 1, sb }, x)
        + (sy[1] * Log(-x, sa) * pow(sy[0], 2.)) / 2. - (sy[1] * pow(sy[0], 3.)) / 6.
        + sy[0] * sy[1] * (sy[3] + 2. * Zeta(2))
        + sy[1]
            * (-sy[2] - sy[3] * sy[4] + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                - 2. * sy[4] * Zeta(2)) };

    if (b != x) {
        res += (-6. * G({ 0, 0, 0, 0, b / a }, 1) + 6. * G({ 0, 0, 0, 0, x / a }, 1))
            * G({ b }, { sb }, x);
    }
    return res;
}

complex<double> G6_000abc(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy = { Log(c, sc), G({ b / c, a / c }, 1),
            G({ 0, b / c, a / c }, 1), G({ b / c, 0, a / c }, 1),
            G({ 0, 0, 0 }, { 1, 1, 1 }, x), G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x), G({ 0, 0, b / c, a / c }, 1),
            G({ 0, b / c, 0, a / c }, 1), G({ b / c, 0, 0, a / c }, 1),
            G({ 0, 0, 0, a, b }, { 1, 1, 1, sa, sb }, x), G({ b / c }, 1),
            G({ 0, 0, 0, b / c, a / c }, 1), G({ 0, 0, b / c, 0, a / c }, 1),
            G({ 0, b / c, 0, 0, a / c }, 1), G({ b / c, 0, 0, 0, a / c }, 1) };
        complex<double> res { (sy[2] + sy[3]) * sy[4]
            + sy[6] * (-sy[9] - sy[10] + (-sy[2] - sy[3]) * sy[7] - sy[8])
            - sy[5] * G({ 0, b / c }, 1) + G({ 0, 0, 0, 0, b / c, a / c }, 1)
            + G({ 0, 0, 0, b / c, 0, a / c }, 1) + G({ 0, 0, b / c, 0, 0, a / c }, 1)
            + G({ 0, b / c, 0, 0, 0, a / c }, 1) + G({ b / c, 0, 0, 0, 0, a / c }, 1)
            - G({ b / c, a / c, 0, 0, 0, x / c }, 1)
            + sy[11] * (-G({ x / c }, 1) + G({ c }, { sc }, x))
            + sy[12] * (sy[11] + 4. * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x))
            - 4. * G({ 0, 0, 0, 0, a, b }, { 1, 1, 1, 1, sa, sb }, x)
            - G({ 0, 0, 0, a, 0, b }, { 1, 1, 1, sa, 1, sb }, x)
            + (-sy[2] / 6. - sy[3] / 6.) * pow(sy[0], 3.)
            + Log(-x, sc)
                * (sy[13] + sy[14] + sy[15] + sy[16] - sy[12] * sy[5]
                    + sy[0] * (-sy[9] - sy[10] - sy[8])
                    + (sy[2] / 2. + sy[3] / 2.) * pow(sy[0], 2.)
                    - (sy[1] * pow(sy[0], 3.)) / 6.)
            + (sy[1] * pow(sy[0], 4.)) / 24.
            + pow(sy[0], 2.)
                * (sy[9] / 2. + sy[10] / 2. + sy[8] / 2.
                    + sy[1] * (-sy[6] / 2. - Zeta(2)))
            - 2. * sy[9] * Zeta(2) - 2. * sy[10] * Zeta(2) - 2. * sy[8] * Zeta(2)
            + sy[7] * (-2. * sy[2] * Zeta(2) - 2. * sy[3] * Zeta(2))
            + sy[0]
                * (-sy[13] - sy[14] - sy[15] - sy[16] + sy[12] * sy[5]
                    + (sy[2] + sy[3]) * sy[6] + 2. * sy[2] * Zeta(2)
                    + 2. * sy[3] * Zeta(2)
                    + sy[1] * (-sy[4] + sy[6] * sy[7] + 2. * sy[7] * Zeta(2)))
            + sy[1]
                * (-sy[5] + sy[4] * sy[7] - G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x)
                    + pow(sy[6], 2.) - 2. * pow(sy[7], 2.) * Zeta(2)
                    + sy[6] * (-pow(sy[7], 2.) + 2. * Zeta(2)) - 2. * Zeta(4)) };
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        if (a == b) // aaC
        {
            if (sa != sb) {
                throw FastGPL::FastGPL_error { "G6_000abc: a==b but sa!=sb" };
            }
            const vector<complex<double>> sy = { Log(a, sa), G({ c / a, 1 }, 1),
                G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x), G({ 0, 0 }, { 1, 1 }, x),
                G({ 0 }, { 1 }, x), G({ 0, 0, 0 }, { 1, 1, 1 }, x) };
            complex<double> res { G({ 0, 0, 0, 0, c / a, 1 }, 1)
                + G({ 0, 0, 0, c / a, 1, x / a }, 1) + G({ 0, 0, 0, c / a, x / a, 1 }, 1)
                + G({ 0, 0, c / a, 0, 1, x / a }, 1) + G({ 0, 0, c / a, 0, x / a, 1 }, 1)
                + G({ 0, 0, c / a, 1, 0, x / a }, 1) + G({ 0, c / a, 0, 0, 1, x / a }, 1)
                + G({ 0, c / a, 0, 0, x / a, 1 }, 1) + G({ 0, c / a, 0, 1, 0, x / a }, 1)
                + G({ 0, c / a, 1, 0, 0, x / a }, 1) + G({ c / a, 0, 0, 0, 1, x / a }, 1)
                + G({ c / a, 0, 0, 0, x / a, 1 }, 1) + G({ c / a, 0, 0, 1, 0, x / a }, 1)
                + G({ c / a, 0, 1, 0, 0, x / a }, 1) + G({ c / a, 1, 0, 0, 0, x / a }, 1)
                + G({ 0, 0, 0, 0, a, c }, { 1, 1, 1, 1, sa, sc }, x)
                + (sy[1] * Log(-x, sa) * pow(sy[0], 3.)) / 6.
                - (sy[1] * pow(sy[0], 4.)) / 24. + sy[2] * Zeta(2)
                - G({ 0, 0, a, c }, { 1, 1, sa, sc }, x) * Zeta(2)
                + sy[1] * pow(sy[0], 2.) * (sy[3] / 2. + Zeta(2))
                + sy[0] * sy[1] * (-(sy[3] * sy[4]) + sy[5] - 2. * sy[4] * Zeta(2))
                + G({ a, c }, { sa, sc }, x) * (-G({ 0, 0, 0, x / a }, 1) - Zeta(4))
                + G({ 0, c }, { 1, sc }, x) * (G({ 0, 0, c / a, 1 }, 1) + Zeta(4))
                + sy[1]
                    * (sy[2] - sy[4] * sy[5] + G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x)
                        - pow(sy[3], 2.) + sy[3] * (pow(sy[4], 2.) - 2. * Zeta(2))
                        + 2. * pow(sy[4], 2.) * Zeta(2) + 2. * Zeta(4))
                + G({ 0, 0, c }, { 1, 1, sc }, x) * (-G({ 0, c / a, 1 }, 1) - Zeta(3))
                + G({ 0, a, c }, { 1, sa, sc }, x) * Zeta(3) };
            if (c != x) {
                res += (-G({ 0, 0, 0, c / a, 1 }, 1) + G({ 0, 0, 0, x / a, 1 }, 1))
                    * G({ c }, { sc }, x);
            }
            return res;
        }
        const vector<complex<double>> sy = { Log(a, sa), G({ b / a, c / a }, 1),
            G({ 0, 0, c }, { 1, 1, sc }, x), G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x),
            G({ 0, 0, 0, b / a }, 1), G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x),
            G({ 0, 0, 0 }, { 1, 1, 1 }, x) };
        complex<double> res { sy[2]
                * (G({ 0, b / a, c / a }, 1) + G({ b / a, 0, c / a }, 1))
            - G({ 0, 0, 0, 0, b / a, c / a }, 1) - G({ 0, 0, 0, b / a, 0, c / a }, 1)
            - G({ 0, 0, 0, b / a, c / a, x / a }, 1) - G({ 0, 0, b / a, 0, 0, c / a }, 1)
            - G({ 0, 0, b / a, 0, c / a, x / a }, 1)
            - G({ 0, 0, b / a, c / a, 0, x / a }, 1) - G({ 0, b / a, 0, 0, 0, c / a }, 1)
            - G({ 0, b / a, 0, 0, c / a, x / a }, 1)
            - G({ 0, b / a, 0, c / a, 0, x / a }, 1)
            - G({ 0, b / a, c / a, 0, 0, x / a }, 1) - G({ b / a, 0, 0, 0, 0, c / a }, 1)
            - G({ b / a, 0, 0, 0, c / a, x / a }, 1)
            - G({ b / a, 0, 0, c / a, 0, x / a }, 1)
            - G({ b / a, 0, c / a, 0, 0, x / a }, 1)
            - G({ b / a, c / a, 0, 0, 0, x / a }, 1)
            + (-sy[4] - G({ 0, 0, b / a, c / a }, 1) - G({ 0, b / a, 0, c / a }, 1)
                  - G({ b / a, 0, 0, c / a }, 1))
                * G({ 0, c }, { 1, sc }, x)
            + G({ 0, 0, b / a }, 1) * (sy[2] - G({ 0, b, c }, { 1, sb, sc }, x))
            + G({ 0, b / a }, 1) * (-sy[3] + G({ 0, 0, b, c }, { 1, 1, sb, sc }, x))
            + G({ b / a }, 1)
                * (G({ 0, 0, 0, 0, c }, { 1, 1, 1, 1, sc }, x)
                    - G({ 0, 0, 0, b, c }, { 1, 1, 1, sb, sc }, x))
            + G({ 0, 0, 0, 0, b, c }, { 1, 1, 1, 1, sb, sc }, x)
            - (sy[1] * Log(-x, sa) * pow(sy[0], 3.)) / 6. + (sy[1] * pow(sy[0], 4.)) / 24.
            + sy[1] * pow(sy[0], 2.) * (-sy[5] / 2. - Zeta(2))
            + sy[0] * sy[1] * (sy[5] * sy[6] - sy[7] + 2. * sy[6] * Zeta(2))
            + sy[1]
                * (-sy[3] + sy[6] * sy[7] - G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x)
                    + pow(sy[5], 2.) - 2. * pow(sy[6], 2.) * Zeta(2)
                    + sy[5] * (-pow(sy[6], 2.) + 2. * Zeta(2)) - 2. * Zeta(4)) }; // aBC
        if (b != x) {
            res += (sy[4] - G({ 0, 0, 0, x / a }, 1)) * G({ b, c }, { sb, sc }, x);
        }
        if (c != x) {
            res += (G({ 0, 0, 0, b / a, c / a }, 1) - G({ 0, 0, 0, b / a, x / a }, 1)
                       + G({ 0, 0, b / a, 0, c / a }, 1) - G({ 0, 0, b / a, 0, x / a }, 1)
                       + G({ 0, b / a, 0, 0, c / a }, 1) - G({ 0, b / a, 0, 0, x / a }, 1)
                       + G({ b / a, 0, 0, 0, c / a }, 1)
                       - G({ b / a, 0, 0, 0, x / a }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest
    if (b == c) { // Abb
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G6_000abc: b==c but sb!=sc" };
        }
        const vector<complex<double>> sy = { G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x) };
        complex<double> res { G({ a / b, 0, 0, 0, 0, 1 }, 1)
            - G({ a / b, 0, 0, 0, x / b, 1 }, 1)
            + (-G({ a / b, 0, 0, 0, 1 }, 1) + G({ a / b, 0, 0, 0, x / b }, 1))
                * G({ b }, { sb }, x)
            + G({ a / b, 0, 0, 1 }, 1) * G({ 0, b }, { 1, sb }, x)
            - G({ a / b, 0, 1 }, 1) * G({ 0, 0, b }, { 1, 1, sb }, x)
            + G({ a / b, 1 }, 1) * (-sy[0] + G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x))
            + G({ a / b }, 1)
                * (-G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x)
                    + G({ 0, 0, 0, a, b }, { 1, 1, 1, sa, sb }, x))
            + G({ 0, 0, 0, a, 0, b }, { 1, 1, 1, sa, 1, sb }, x) - sy[0] * Zeta(2) };
        return res;
    } else { // AbC
        const vector<complex<double>> sy = { Log(b, sb), G({ a / b, c / b }, 1),
            G({ c / b, a / b }, 1), G({ 0, 0 }, { 1, 1 }, x),
            G({ 0, 0, 0 }, { 1, 1, 1 }, x), G({ 0, a / b, c / b }, 1),
            G({ 0, c / b, a / b }, 1), G({ c / b, 0, a / b }, 1),
            G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x), G({ c / b }, 1),
            G({ 0, 0, 0, a, c }, { 1, 1, 1, sa, sc }, x), G({ 0, 0, a / b, c / b }, 1),
            G({ 0, 0, c / b, a / b }, 1), G({ 0, c / b, 0, a / b }, 1),
            G({ c / b, 0, 0, a / b }, 1), G({ 0, 0, 0, a / b, c / b }, 1),
            G({ 0, 0, 0, c / b, a / b }, 1), G({ 0, 0, c / b, 0, a / b }, 1),
            G({ 0, c / b, 0, 0, a / b }, 1), G({ c / b, 0, 0, 0, a / b }, 1),
            G({ 0 }, { 1 }, x), G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x) };
        complex<double> res { sy[4] * (-sy[5] - sy[6] - sy[7])
            + sy[8] * G({ 0, c / b }, 1) - G({ 0, 0, 0, 0, a / b, c / b }, 1)
            - G({ 0, 0, 0, 0, c / b, a / b }, 1) - G({ 0, 0, 0, c / b, 0, a / b }, 1)
            - G({ 0, 0, c / b, 0, 0, a / b }, 1) - G({ 0, c / b, 0, 0, 0, a / b }, 1)
            + G({ a / b, 0, 0, 0, 0, c / b }, 1) + G({ a / b, 0, 0, 0, c / b, x / b }, 1)
            + G({ a / b, 0, 0, c / b, 0, x / b }, 1)
            + G({ a / b, 0, c / b, 0, 0, x / b }, 1)
            + G({ a / b, c / b, 0, 0, 0, x / b }, 1) - G({ c / b, 0, 0, 0, 0, a / b }, 1)
            + G({ c / b, a / b, 0, 0, 0, x / b }, 1)
            + G({ a / b, 0, 0, c / b }, 1) * G({ 0, c }, { 1, sc }, x)
            - G({ a / b, 0, c / b }, 1) * G({ 0, 0, c }, { 1, 1, sc }, x)
            + sy[9] * (-sy[10] - 4. * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x))
            + G({ a / b }, 1) * (sy[10] - G({ 0, 0, 0, 0, c }, { 1, 1, 1, 1, sc }, x))
            + G({ 0, 0, 0, a, 0, c }, { 1, 1, 1, sa, 1, sc }, x)
            + (sy[5] / 6. + sy[6] / 6. + sy[7] / 6.) * pow(sy[0], 3.)
            + Log(-x, sb)
                * (sy[0] * (sy[11] + sy[12] + sy[13] + sy[14]) - sy[15] - sy[16] - sy[17]
                    - sy[18] - sy[19] + sy[9] * sy[8]
                    + (-sy[5] / 2. - sy[6] / 2. - sy[7] / 2.) * pow(sy[0], 2.)
                    + (sy[1] / 6. + sy[2] / 6.) * pow(sy[0], 3.))
            + (-sy[1] / 24. - sy[2] / 24.) * pow(sy[0], 4.) - sy[2] * pow(sy[3], 2.)
            + 2. * sy[11] * Zeta(2) + 2. * sy[12] * Zeta(2) + 2. * sy[13] * Zeta(2)
            + 2. * sy[14] * Zeta(2) + 2. * sy[2] * pow(sy[20], 2.) * Zeta(2)
            + sy[3]
                * (sy[11] + sy[12] + sy[13] + sy[14] + sy[20] * (sy[5] + sy[6] + sy[7])
                    + sy[2] * pow(sy[20], 2.) - 2. * sy[2] * Zeta(2))
            + sy[20]
                * (-(sy[2] * sy[4]) + 2. * sy[5] * Zeta(2) + 2. * sy[6] * Zeta(2)
                    + 2. * sy[7] * Zeta(2))
            + pow(sy[0], 2.)
                * (-sy[11] / 2. - sy[12] / 2. - sy[13] / 2. - sy[14] / 2.
                    + (sy[2] * sy[3]) / 2. + sy[2] * Zeta(2)
                    + sy[1] * (sy[3] / 2. + Zeta(2)))
            + sy[0]
                * (sy[15] + sy[16] + sy[17] + sy[18] + sy[19] + sy[2] * sy[4]
                    + sy[3] * (-(sy[20] * sy[2]) - sy[5] - sy[6] - sy[7]) - sy[9] * sy[8]
                    - 2. * sy[20] * sy[2] * Zeta(2) - 2. * sy[5] * Zeta(2)
                    - 2. * sy[6] * Zeta(2) - 2. * sy[7] * Zeta(2)
                    + sy[1] * (-(sy[20] * sy[3]) + sy[4] - 2. * sy[20] * Zeta(2)))
            + sy[2] * (sy[21] + sy[8] + 2. * Zeta(4))
            + sy[1]
                * (sy[21] - sy[20] * sy[4] + G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                    - pow(sy[3], 2.) + sy[3] * (pow(sy[20], 2.) - 2. * Zeta(2))
                    + 2. * pow(sy[20], 2.) * Zeta(2) + 2. * Zeta(4)) };
        if (c != x) {
            res += (-G({ a / b, 0, 0, 0, c / b }, 1) + G({ a / b, 0, 0, 0, x / b }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G6_a000bc(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy = { Log(c, sc), G({ a }, { sa }, x),
            G({ b / c }, 1), G({ 0, 0 }, { 1, 1 }, x), G({ 0, a }, { 1, sa }, x),
            G({ 0, b / c }, 1), G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, 0, b / c }, 1),
            G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x),
            G({ a, 0, 0, 0, b }, { sa, 1, 1, 1, sb }, x), G({ 0, 0, 0, b / c }, 1),
            G({ b / c, 0, 0, 0, a / c }, 1), G({ 0, 0, 0 }, { 1, 1, 1 }, x),
            G({ 0 }, { 1 }, x) };
        complex<double> res { -(sy[3] * sy[4] * sy[5]) - sy[6] * sy[7] + sy[5] * sy[8]
            - G({ 0, b / c, 0, 0, 0, a / c }, 1) - 4. * G({ b / c, 0, 0, 0, 0, a / c }, 1)
            - G({ b / c, 0, 0, 0, a / c, x / c }, 1)
            + sy[9] * (-G({ x / c }, 1) + G({ c }, { sc }, x))
            - G({ 0, a, 0, 0, 0, b }, { 1, sa, 1, 1, 1, sb }, x)
            - 4. * G({ a, 0, 0, 0, 0, b }, { sa, 1, 1, 1, 1, sb }, x)
            + (-(sy[2] * sy[4]) / 6. + (sy[1] * sy[5]) / 6.) * pow(sy[0], 3.)
            + Log(-x, sc)
                * (-sy[11] - sy[10] * sy[1] - sy[5] * sy[6] + sy[4] * sy[7]
                    + sy[0] * (-(sy[4] * sy[5]) + sy[2] * sy[6] + sy[1] * sy[7])
                    + sy[2] * sy[8]
                    + ((sy[2] * sy[4]) / 2. - (sy[1] * sy[5]) / 2.) * pow(sy[0], 2.)
                    + (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.)
            - (sy[1] * sy[2] * pow(sy[0], 4.)) / 24.
            + sy[4] * (sy[10] - 2. * sy[5] * Zeta(2))
            + sy[2]
                * (sy[9] + sy[12] * sy[4] + sy[3] * (-(sy[13] * sy[4]) + sy[6])
                    - G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
                    - 2. * sy[13] * sy[4] * Zeta(2) + 2. * sy[6] * Zeta(2))
            + pow(sy[0], 2.)
                * ((sy[4] * sy[5]) / 2. - (sy[2] * sy[6]) / 2.
                    + sy[1] * (-sy[7] / 2. + sy[2] * (sy[3] / 2. + Zeta(2))))
            + sy[0]
                * (sy[11] + sy[5] * sy[6] - sy[4] * sy[7]
                    + sy[2] * (sy[3] * sy[4] - sy[8] + 2. * sy[4] * Zeta(2))
                    + sy[1]
                        * (sy[10] - sy[3] * sy[5] - 2. * sy[5] * Zeta(2)
                            + sy[2] * (sy[12] - sy[13] * sy[3] - 2. * sy[13] * Zeta(2))))
            + sy[1]
                * (sy[11] - sy[12] * sy[5] + sy[3] * (sy[13] * sy[5] + sy[7])
                    - G({ 0, 0, 0, 0, b / c }, 1) + 2. * sy[13] * sy[5] * Zeta(2)
                    + 2. * sy[7] * Zeta(2)
                    + sy[2]
                        * (-(sy[12] * sy[13]) + G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x)
                            - pow(sy[3], 2.) + sy[3] * (pow(sy[13], 2.) - 2. * Zeta(2))
                            + 2. * pow(sy[13], 2.) * Zeta(2) + 2. * Zeta(4))) };
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        const vector<complex<double>> sy
            = { G({ 0, 0, 0, b / a }, 1), G({ 0, 0, 0, b / a, c / a }, 1) };
        complex<double> res { 4. * G({ 0, 0, 0, 0, b / a, c / a }, 1)
            + G({ 0, 0, 0, b / a, 0, c / a }, 1) + G({ 0, 0, 0, b / a, c / a, x / a }, 1)
            + sy[0] * G({ 0, c }, { 1, sc }, x)
            + G({ 0, 0, x / a }, 1) * G({ 0, b, c }, { 1, sb, sc }, x)
            + G({ 0, x / a }, 1) * G({ 0, 0, b, c }, { 1, 1, sb, sc }, x)
            + G({ x / a }, 1) * G({ 0, 0, 0, b, c }, { 1, 1, 1, sb, sc }, x)
            + G({ 0, 0, 0, 0, b, c }, { 1, 1, 1, 1, sb, sc }, x) - sy[1] * Log(a, sa)
            + sy[1] * Log(-x, sa) }; // aBC
        if (b != x) {
            res += (-sy[0] + G({ 0, 0, 0, x / a }, 1)) * G({ b, c }, { sb, sc }, x);
        }
        if (c != x) {
            res += (-sy[1] + G({ 0, 0, 0, b / a, x / a }, 1)) * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest
    if (b == c) { // Abb
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G6_a000bc: b==c but sb!=sc" };
        }
        const vector<complex<double>> sy = { G({ 0, 0, 0, a / b }, 1),
            G({ a }, { sa }, x), G({ 0, 0, 0, a / b, 1 }, 1), G({ b }, { sb }, x) };
        complex<double> res { sy[1] * sy[2] - sy[2] * sy[3]
            + sy[3] * G({ 0, 0, 0, a / b, x / b }, 1) + G({ 0, 0, 0, a / b, 0, 1 }, 1)
            - G({ 0, 0, 0, a / b, x / b, 1 }, 1) + sy[0] * G({ 0, b }, { 1, sb }, x)
            + G({ a, 0, 0, 0, 0, b }, { sa, 1, 1, 1, 1, sb }, x)
            - G({ a, 0, 0, b }, { sa, 1, 1, sb }, x) * Zeta(2)
            + G({ a, b }, { sa, sb }, x) * (-sy[0] - Zeta(4))
            + G({ a, 0, b }, { sa, 1, sb }, x) * Zeta(3) + sy[1] * Zeta(5) };
        return res;
    } else { // AbC
        const vector<complex<double>> sy = { Log(b, sb), G({ a }, { sa }, x),
            G({ c / b }, 1), G({ 0, a }, { 1, sa }, x), G({ 0, 0, 0, a / b }, 1),
            G({ a, c }, { sa, sc }, x), G({ 0, 0, a }, { 1, 1, sa }, x),
            G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x), G({ 0, 0, 0, a / b, c / b }, 1),
            G({ 0, 0, 0, c / b, a / b }, 1), G({ 0, 0, c / b, 0, a / b }, 1),
            G({ 0, c / b, 0, 0, a / b }, 1), G({ c / b, 0, 0, 0, a / b }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0, 0, 0 }, { 1, 1, 1 }, x),
            G({ 0 }, { 1 }, x) };
        complex<double> res { -(sy[4] * sy[5]) + sy[5] * G({ 0, 0, 0, c / b }, 1)
            + 4. * G({ 0, 0, 0, 0, a / b, c / b }, 1)
            + 4. * G({ 0, 0, 0, 0, c / b, a / b }, 1) + G({ 0, 0, 0, a / b, 0, c / b }, 1)
            + G({ 0, 0, 0, a / b, c / b, x / b }, 1)
            + 4. * G({ 0, 0, 0, c / b, 0, a / b }, 1)
            + G({ 0, 0, 0, c / b, a / b, x / b }, 1)
            + 4. * G({ 0, 0, c / b, 0, 0, a / b }, 1)
            + G({ 0, 0, c / b, 0, a / b, x / b }, 1)
            + 4. * G({ 0, c / b, 0, 0, 0, a / b }, 1)
            + G({ 0, c / b, 0, 0, a / b, x / b }, 1)
            + 4. * G({ c / b, 0, 0, 0, 0, a / b }, 1)
            + G({ c / b, 0, 0, 0, a / b, x / b }, 1) + sy[4] * G({ 0, c }, { 1, sc }, x)
            - G({ 0, 0, c / b }, 1) * G({ a, 0, c }, { sa, 1, sc }, x)
            + G({ 0, c / b }, 1) * G({ a, 0, 0, c }, { sa, 1, 1, sc }, x)
            + G({ a, 0, 0, 0, 0, c }, { sa, 1, 1, 1, 1, sc }, x)
            + (sy[2] * sy[3] * pow(sy[0], 3.)) / 6.
            + Log(-x, sb)
                * (sy[9] + sy[10] + sy[11] + sy[12] - sy[0] * sy[2] * sy[6]
                    - sy[2] * sy[7] + sy[8] - (sy[2] * sy[3] * pow(sy[0], 2.)) / 2.
                    - (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.)
            + (sy[1] * sy[2] * pow(sy[0], 4.)) / 24.
            + pow(sy[0], 2.)
                * ((sy[2] * sy[6]) / 2. + sy[1] * sy[2] * (-sy[13] / 2. - Zeta(2)))
            + sy[2]
                * (-(sy[14] * sy[3]) + sy[13] * (sy[15] * sy[3] - sy[6])
                    + G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
                    - G({ a, 0, 0, 0, c }, { sa, 1, 1, 1, sc }, x)
                    + 2. * sy[15] * sy[3] * Zeta(2) - 2. * sy[6] * Zeta(2))
            + sy[0]
                * (-sy[9] - sy[10] - sy[11] - sy[12] - sy[8]
                    + sy[1] * sy[2] * (-sy[14] + sy[13] * sy[15] + 2. * sy[15] * Zeta(2))
                    + sy[2] * (-(sy[13] * sy[3]) + sy[7] - 2. * sy[3] * Zeta(2)))
            + sy[1]
                * (-sy[9] - sy[10] - sy[11] - sy[12] - G({ 0, 0, 0, 0, c / b }, 1)
                    + sy[2]
                        * (sy[14] * sy[15] - G({ 0, 0, 0, 0 }, { 1, 1, 1, 1 }, x)
                            + pow(sy[13], 2.) - 2. * pow(sy[15], 2.) * Zeta(2)
                            + sy[13] * (-pow(sy[15], 2.) + 2. * Zeta(2))
                            - 2. * Zeta(4))) };
        if (c != x) {
            res += (-sy[8] + G({ 0, 0, 0, a / b, x / b }, 1)) * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G6_ab000c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy = { G({ 0, 0, c }, { 1, 1, sc }, x),
            G({ 0, b, a }, { 1, sb, sa }, x), G({ 0, 0, x / c }, 1),
            G({ 0, 0, 0, b / c }, 1), G({ 0, c }, { 1, sc }, x),
            G({ 0, 0, b, a }, { 1, 1, sb, sa }, x), G({ 0, x / c }, 1),
            G({ c }, { sc }, x), G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x),
            G({ x / c }, 1), Log(c, sc), G({ 0, 0, 0, b / c, a / c }, 1),
            G({ a }, { sa }, x), G({ 0, 0, b }, { 1, 1, sb }, x),
            G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x), G({ 0, 0, 0, x / c }, 1) };
        complex<double> res { sy[10] * sy[11] + sy[0] * sy[1] - sy[1] * sy[2]
            - sy[4] * sy[5] - sy[5] * sy[6] - sy[9] * sy[8] + sy[7] * sy[8]
            - 4. * G({ 0, 0, 0, 0, b / c, a / c }, 1) - G({ 0, 0, 0, b / c, 0, a / c }, 1)
            - G({ 0, 0, 0, b / c, a / c, x / c }, 1) - sy[3] * G({ 0, a }, { 1, sa }, x)
            + G({ a, b }, { sa, sb }, x) * G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
            + sy[12]
                * (sy[11] + sy[9] * sy[14] - sy[10] * sy[3] + sy[13] * sy[4]
                    + sy[13] * sy[6] - sy[14] * sy[7] + 4. * G({ 0, 0, 0, 0, b / c }, 1)
                    + (-sy[0] + sy[2]) * G({ 0, b }, { 1, sb }, x)
                    + G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x))
            - G({ 0, 0, 0, 0, b, a }, { 1, 1, 1, 1, sb, sa }, x)
            + (-sy[11] + sy[12] * sy[3]) * Log(-x, sc) };
        if (b != x) {
            res += (sy[12] * sy[15] - sy[12] * sy[3]) * G({ b }, { sb }, x)
                + (-sy[15] + sy[3]) * G({ b, a }, { sb, sa }, x);
        }
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        if (a == b) // aaC
        {
            if (sa != sb) {
                throw FastGPL::FastGPL_error { "G6_ab000c: a==b but sa!=sb" };
            }
            const vector<complex<double>> sy
                = { G({ 0, 0, 0, 1, c / a }, 1), G({ 0, 0, 0, c / a, 1 }, 1),
                      G({ 0, 0, 1, 0, c / a }, 1), G({ 0, 1, 0, 0, c / a }, 1) };
            complex<double> res { -4. * G({ 0, 0, 0, 0, 1, c / a }, 1)
                - 4. * G({ 0, 0, 0, 0, c / a, 1 }, 1)
                - 4. * G({ 0, 0, 0, 1, 0, c / a }, 1) - G({ 0, 0, 0, 1, c / a, x / a }, 1)
                - G({ 0, 0, 0, c / a, 1, x / a }, 1) - G({ 0, 0, 0, c / a, x / a, 1 }, 1)
                - 4. * G({ 0, 0, 1, 0, 0, c / a }, 1) - G({ 0, 0, 1, 0, c / a, x / a }, 1)
                - 3. * G({ 0, 1, 0, 0, 0, c / a }, 1) - G({ 0, 1, 0, 0, c / a, x / a }, 1)
                + (-G({ 0, 0, 1, x / a }, 1) - G({ 0, 0, x / a, 1 }, 1)
                      - G({ 0, 1, 0, x / a }, 1))
                    * G({ 0, c }, { 1, sc }, x)
                + (-G({ 0, 1, x / a }, 1) - G({ 0, x / a, 1 }, 1))
                    * G({ 0, 0, c }, { 1, 1, sc }, x)
                - G({ x / a, 1 }, 1) * G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                + G({ x / a }, 1) * G({ a, 0, 0, 0, c }, { sa, 1, 1, 1, sc }, x)
                + G({ 0, a, 0, 0, 0, c }, { 1, sa, 1, 1, 1, sc }, x)
                + (sy[0] + sy[1] + sy[2] + sy[3]) * Log(a, sa)
                + (-sy[0] - sy[1] - sy[2] - sy[3]) * Log(-x, sa) };
            if (c != x) {
                res += (sy[0] + sy[1] + sy[2] + sy[3] - G({ 0, 0, 0, 1, x / a }, 1)
                           - G({ 0, 0, 0, x / a, 1 }, 1) - G({ 0, 0, 1, 0, x / a }, 1)
                           - G({ 0, 1, 0, 0, x / a }, 1))
                    * G({ c }, { sc }, x);
            }
            return res;
        }

        const vector<complex<double>> sy
            = { G({ b / a }, 1), G({ b / a, 0, 0, 0, c / a }, 1) };
        complex<double> res { G({ 0, b / a, 0, 0, 0, c / a }, 1)
            + 4. * G({ b / a, 0, 0, 0, 0, c / a }, 1)
            + G({ b / a, 0, 0, 0, c / a, x / a }, 1)
            + G({ b / a, 0, 0, x / a }, 1) * G({ 0, c }, { 1, sc }, x)
            + G({ b / a, 0, x / a }, 1) * G({ 0, 0, c }, { 1, 1, sc }, x)
            + G({ b / a, x / a }, 1) * G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
            + sy[0] * G({ 0, 0, 0, 0, c }, { 1, 1, 1, 1, sc }, x)
            + G({ 0, b, 0, 0, 0, c }, { 1, sb, 1, 1, 1, sc }, x) - sy[1] * Log(a, sa)
            + sy[1] * Log(-x, sa) }; // aBC
        if (b != x) {
            res += (-sy[0] + G({ x / a }, 1))
                * G({ b, 0, 0, 0, c }, { sb, 1, 1, 1, sc }, x);
        }
        if (c != x) {
            res += (-sy[1] + G({ b / a, 0, 0, 0, x / a }, 1)) * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest
    const vector<complex<double>> sy = { G({ 0, 0, c }, { 1, 1, sc }, x),
        G({ a, c }, { sa, sc }, x), G({ 0, 0, 0, a / b }, 1),
        G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x), G({ 0, 0, 0, c / b }, 1),
        G({ a }, { sa }, x), G({ 0, 0, 0, c / b, a / b }, 1),
        G({ 0, 0, 0, a / b, c / b }, 1), G({ 0, 0, a / b, 0, c / b }, 1),
        G({ 0, a / b, 0, 0, c / b }, 1), G({ a / b, 0, 0, 0, c / b }, 1) };
    complex<double> res { sy[1] * sy[2] - sy[1] * sy[4] - sy[3] * G({ a / b, x / b }, 1)
        + sy[0] * (-G({ 0, a / b, x / b }, 1) - G({ a / b, 0, x / b }, 1))
        + sy[5] * (sy[6] + 4. * G({ 0, 0, 0, 0, c / b }, 1))
        - 4. * G({ 0, 0, 0, 0, a / b, c / b }, 1)
        - 4. * G({ 0, 0, 0, 0, c / b, a / b }, 1)
        - 4. * G({ 0, 0, 0, a / b, 0, c / b }, 1) - G({ 0, 0, 0, a / b, c / b, x / b }, 1)
        - G({ 0, 0, 0, c / b, 0, a / b }, 1) - G({ 0, 0, 0, c / b, a / b, x / b }, 1)
        - 4. * G({ 0, 0, a / b, 0, 0, c / b }, 1) - G({ 0, 0, a / b, 0, c / b, x / b }, 1)
        - 4. * G({ 0, a / b, 0, 0, 0, c / b }, 1) - G({ 0, a / b, 0, 0, c / b, x / b }, 1)
        - 4. * G({ a / b, 0, 0, 0, 0, c / b }, 1) - G({ a / b, 0, 0, 0, c / b, x / b }, 1)
        - sy[4] * G({ 0, a }, { 1, sa }, x)
        + (-sy[2] - G({ 0, 0, a / b, x / b }, 1) - G({ 0, a / b, 0, x / b }, 1)
              - G({ a / b, 0, 0, x / b }, 1))
            * G({ 0, c }, { 1, sc }, x)
        + G({ 0, 0, a / b }, 1) * (-sy[0] + G({ a, 0, c }, { sa, 1, sc }, x))
        + G({ 0, a / b }, 1) * (-sy[3] + G({ a, 0, 0, c }, { sa, 1, 1, sc }, x))
        + G({ a / b }, 1)
            * (-G({ 0, 0, 0, 0, c }, { 1, 1, 1, 1, sc }, x)
                + G({ a, 0, 0, 0, c }, { sa, 1, 1, 1, sc }, x))
        + G({ a, 0, 0, 0, 0, c }, { sa, 1, 1, 1, 1, sc }, x)
        + (sy[9] + sy[10] - sy[4] * sy[5] + sy[6] + sy[7] + sy[8]) * Log(b, sb)
        + (-sy[9] - sy[10] + sy[4] * sy[5] - sy[6] - sy[7] - sy[8]) * Log(-x, sb) };
    if (c != x) {
        res += (sy[9] + sy[10] + sy[7] + sy[8] - G({ 0, 0, 0, a / b, x / b }, 1)
                   - G({ 0, 0, a / b, 0, x / b }, 1) - G({ 0, a / b, 0, 0, x / b }, 1)
                   - G({ a / b, 0, 0, 0, x / b }, 1))
            * G({ c }, { sc }, x);
    }
    return res;
}
complex<double> G6_00a0bc(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy = { Log(c, sc), G({ b / c, 0, a / c }, 1),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, b / c }, 1),
            G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x), G({ b / c }, 1),
            G({ 0, b / c, 0, a / c }, 1), G({ b / c, 0, 0, a / c }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x),
            G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x), G({ 0, 0, b / c, 0, a / c }, 1),
            G({ 0, b / c, 0, 0, a / c }, 1), G({ b / c, 0, 0, 0, a / c }, 1) };
        complex<double> res { 3. * sy[3] * sy[4]
            + (sy[9] * sy[1] + sy[6] + 2. * sy[7]) * sy[8]
            + sy[2] * (sy[1] - G({ 0, 0, b / c }, 1)) - G({ 0, 0, 0, b / c, 0, a / c }, 1)
            - 2. * G({ 0, 0, b / c, 0, 0, a / c }, 1)
            - 3. * G({ 0, b / c, 0, 0, 0, a / c }, 1)
            - 4. * G({ b / c, 0, 0, 0, 0, a / c }, 1)
            - G({ b / c, 0, a / c, 0, 0, x / c }, 1)
            + sy[10] * (-G({ x / c }, 1) + G({ c }, { sc }, x))
            - 3. * G({ 0, 0, 0, a, 0, b }, { 1, 1, 1, sa, 1, sb }, x)
            - 2. * G({ 0, 0, a, 0, 0, b }, { 1, 1, sa, 1, 1, sb }, x)
            + (-(sy[2] * sy[5]) / 2. - sy[6] / 2. - sy[7]) * pow(sy[0], 2.)
            + Log(-x, sc)
                * (-sy[11] - 2. * sy[12] - 3. * sy[13] - sy[2] * sy[3]
                    + 3. * sy[4] * sy[5] + sy[0] * (sy[2] * sy[5] + sy[6] + 2. * sy[7])
                    - (sy[1] * pow(sy[0], 2.)) / 2.)
            + (sy[1] * pow(sy[0], 3.)) / 6. + 2. * sy[6] * Zeta(2) + 4. * sy[7] * Zeta(2)
            + sy[1] * (-G({ 0, 0, 0 }, { 1, 1, 1 }, x) + 2. * sy[9] * Zeta(2))
            + sy[0]
                * (sy[11] + 2. * sy[12] + 3. * sy[13] + sy[2] * sy[3] - 3. * sy[4] * sy[5]
                    - sy[1] * sy[8] - 2. * sy[1] * Zeta(2))
            + sy[5]
                * (sy[10] - 6. * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
                    + sy[2] * (sy[8] + 2. * Zeta(2))) };
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        const vector<complex<double>> sy = { G({ 0, b, c }, { 1, sb, sc }, x),
            G({ 0, 0, c }, { 1, 1, sc }, x), Log(a, sa), G({ 0, b / a, c / a }, 1),
            G({ 0, 0, 0, b / a }, 1), G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x) };
        complex<double> res { -(sy[1] * sy[3])
            + (2. * sy[0] - 2. * sy[1]) * G({ 0, 0, b / a }, 1)
            + sy[0] * G({ 0, 0, x / a }, 1) + 4. * G({ 0, 0, 0, 0, b / a, c / a }, 1)
            + 3. * G({ 0, 0, 0, b / a, 0, c / a }, 1)
            + 3. * G({ 0, 0, 0, b / a, c / a, x / a }, 1)
            + 2. * G({ 0, 0, b / a, 0, 0, c / a }, 1)
            + 2. * G({ 0, 0, b / a, 0, c / a, x / a }, 1)
            + 2. * G({ 0, 0, b / a, c / a, 0, x / a }, 1)
            + G({ 0, b / a, 0, 0, 0, c / a }, 1) + G({ 0, b / a, 0, 0, c / a, x / a }, 1)
            + G({ 0, b / a, 0, c / a, 0, x / a }, 1)
            + G({ 0, b / a, c / a, 0, 0, x / a }, 1)
            + (3. * sy[4] + 2. * G({ 0, 0, b / a, c / a }, 1)
                  + G({ 0, b / a, 0, c / a }, 1))
                * G({ 0, c }, { 1, sc }, x)
            + G({ 0, b / a }, 1)
                * (G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                    - G({ 0, 0, b, c }, { 1, 1, sb, sc }, x))
            + G({ 0, 0, 0, 0, b, c }, { 1, 1, 1, 1, sb, sc }, x)
            + (sy[3] * Log(-x, sa) * pow(sy[2], 2.)) / 2. - (sy[3] * pow(sy[2], 3.)) / 6.
            + sy[2] * sy[3] * (sy[5] + 2. * Zeta(2))
            + sy[3]
                * (-(sy[5] * sy[6]) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                    - 2. * sy[6] * Zeta(2)) }; // aBC
        if (b != x) {
            res += (-3. * sy[4] + 3. * G({ 0, 0, 0, x / a }, 1))
                * G({ b, c }, { sb, sc }, x);
        }
        if (c != x) {
            res += (-3. * G({ 0, 0, 0, b / a, c / a }, 1)
                       + 3. * G({ 0, 0, 0, b / a, x / a }, 1)
                       - 2. * G({ 0, 0, b / a, 0, c / a }, 1)
                       + 2. * G({ 0, 0, b / a, 0, x / a }, 1)
                       - G({ 0, b / a, 0, 0, c / a }, 1)
                       + G({ 0, b / a, 0, 0, x / a }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest
    if (b == c) { // Abb
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G6_00a0bc: b==c but sb!=sc" };
        }
        const vector<complex<double>> sy
            = { G({ 0, a / b, 1 }, 1), G({ 0, 0, a, b }, { 1, 1, sa, sb }, x) };
        complex<double> res { G({ 0, a / b, 0, 0, 0, 1 }, 1)
            - G({ 0, a / b, 0, 0, x / b, 1 }, 1)
            + (-G({ 0, a / b, 0, 0, 1 }, 1) + G({ 0, a / b, 0, 0, x / b }, 1))
                * G({ b }, { sb }, x)
            + G({ 0, a / b, 0, 1 }, 1) * G({ 0, b }, { 1, sb }, x)
            - sy[0] * G({ 0, 0, b }, { 1, 1, sb }, x)
            + G({ 0, a / b }, 1) * (-sy[1] + G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x))
            + G({ 0, 0, a, 0, 0, b }, { 1, 1, sa, 1, 1, sb }, x) - sy[1] * Zeta(2)
            + G({ 0, 0, a }, { 1, 1, sa }, x) * (sy[0] + Zeta(3)) };
        return res;
    } else { // AbC
        const vector<complex<double>> sy = { G({ 0, 0, a }, { 1, 1, sa }, x),
            G({ 0, c / b, a / b }, 1), G({ c / b, 0, a / b }, 1), Log(b, sb),
            G({ 0, a / b, c / b }, 1), G({ 0, 0, 0 }, { 1, 1, 1 }, x),
            G({ 0, 0, a, c }, { 1, 1, sa, sc }, x), G({ 0, 0 }, { 1, 1 }, x),
            G({ 0 }, { 1 }, x), G({ 0, 0, a / b, c / b }, 1),
            G({ 0, 0, c / b, a / b }, 1), G({ 0, c / b, 0, a / b }, 1),
            G({ c / b, 0, 0, a / b }, 1), G({ c / b }, 1),
            G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x), G({ 0, 0, 0, a / b, c / b }, 1),
            G({ 0, 0, 0, c / b, a / b }, 1), G({ 0, 0, c / b, 0, a / b }, 1),
            G({ 0, c / b, 0, 0, a / b }, 1), G({ c / b, 0, 0, 0, a / b }, 1) };
        complex<double> res { (sy[1] + sy[2]) * sy[5]
            + sy[7]
                * (-2. * sy[9] - 2. * sy[10] - 2. * sy[11] - 2. * sy[12]
                    + (-sy[1] - sy[2]) * sy[8] - sy[4] * sy[8])
            + sy[6] * G({ 0, c / b }, 1)
            + sy[0] * (-sy[1] - sy[2] - G({ 0, 0, c / b }, 1))
            + 4. * G({ 0, 0, 0, 0, a / b, c / b }, 1)
            + 4. * G({ 0, 0, 0, 0, c / b, a / b }, 1)
            + 4. * G({ 0, 0, 0, c / b, 0, a / b }, 1)
            + 4. * G({ 0, 0, c / b, 0, 0, a / b }, 1) + G({ 0, a / b, 0, 0, 0, c / b }, 1)
            + G({ 0, a / b, 0, 0, c / b, x / b }, 1)
            + G({ 0, a / b, 0, c / b, 0, x / b }, 1)
            + G({ 0, a / b, c / b, 0, 0, x / b }, 1)
            + 4. * G({ 0, c / b, 0, 0, 0, a / b }, 1)
            + G({ 0, c / b, a / b, 0, 0, x / b }, 1)
            + 4. * G({ c / b, 0, 0, 0, 0, a / b }, 1)
            + G({ c / b, 0, a / b, 0, 0, x / b }, 1)
            + G({ 0, a / b, 0, c / b }, 1) * G({ 0, c }, { 1, sc }, x)
            + G({ 0, a / b }, 1) * (-sy[6] + G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x))
            + G({ 0, 0, a, 0, 0, c }, { 1, 1, sa, 1, 1, sc }, x)
            + (sy[9] + sy[10] + sy[11] + sy[12] + (sy[0] * sy[13]) / 2.) * pow(sy[3], 2.)
            + Log(-x, sb)
                * (-3. * sy[13] * sy[14] + 3. * sy[15] + 3. * sy[16] + 3. * sy[17]
                    + 3. * sy[18] + 3. * sy[19]
                    + (-2. * sy[9] - 2. * sy[10] - 2. * sy[11] - 2. * sy[12]
                          - sy[0] * sy[13])
                        * sy[3]
                    + (sy[1] / 2. + sy[2] / 2. + sy[4] / 2.) * pow(sy[3], 2.))
            + (-sy[1] / 6. - sy[2] / 6. - sy[4] / 6.) * pow(sy[3], 3.)
            + sy[13]
                * (6. * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
                    - G({ 0, 0, a, 0, c }, { 1, 1, sa, 1, sc }, x)
                    + sy[0] * (-sy[7] - 2. * Zeta(2)))
            - 4. * sy[9] * Zeta(2) - 4. * sy[10] * Zeta(2) - 4. * sy[11] * Zeta(2)
            - 4. * sy[12] * Zeta(2)
            + sy[8] * (-2. * sy[1] * Zeta(2) - 2. * sy[2] * Zeta(2))
            + sy[3]
                * (3. * sy[13] * sy[14] - 3. * sy[15] - 3. * sy[16] - 3. * sy[17]
                    - 3. * sy[18] - 3. * sy[19] + (sy[1] + sy[2] + sy[4]) * sy[7]
                    + 2. * sy[1] * Zeta(2) + 2. * sy[2] * Zeta(2) + 2. * sy[4] * Zeta(2))
            + sy[4] * (sy[5] - G({ 0, 0, c }, { 1, 1, sc }, x) - 2. * sy[8] * Zeta(2)) };
        if (c != x) {
            res += (-G({ 0, a / b, 0, 0, c / b }, 1) + G({ 0, a / b, 0, 0, x / b }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G6_0a00bc(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy = { Log(c, sc), G({ b / c }, 1),
            G({ 0, a }, { 1, sa }, x), G({ 0, 0, a }, { 1, 1, sa }, x),
            G({ 0, 0, b / c }, 1), G({ 0, b / c }, 1),
            G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x), G({ b / c, 0, 0, a / c }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0, a, 0, 0, b }, { 1, sa, 1, 1, sb }, x),
            G({ 0, b / c, 0, 0, a / c }, 1), G({ b / c, 0, 0, 0, a / c }, 1),
            G({ 0 }, { 1 }, x) };
        complex<double> res { 2. * sy[3] * sy[4] - 3. * sy[5] * sy[6] - sy[7] * sy[8]
            + G({ 0, 0, b / c, 0, 0, a / c }, 1) + 3. * G({ 0, b / c, 0, 0, 0, a / c }, 1)
            + 6. * G({ b / c, 0, 0, 0, 0, a / c }, 1)
            - G({ b / c, 0, 0, a / c, 0, x / c }, 1)
            + sy[9] * (-G({ x / c }, 1) + G({ c }, { sc }, x))
            - 2. * G({ 0, 0, a, 0, 0, b }, { 1, 1, sa, 1, 1, sb }, x)
            - 3. * G({ 0, a, 0, 0, 0, b }, { 1, sa, 1, 1, 1, sb }, x)
            + (sy[1] * sy[3] - (sy[2] * sy[5]) / 2. + sy[7] / 2.) * pow(sy[0], 2.)
            + Log(-x, sc)
                * (sy[10] + 3. * sy[11] - sy[2] * sy[4] + 2. * sy[3] * sy[5]
                    - 3. * sy[1] * sy[6]
                    + sy[0] * (-2. * sy[1] * sy[3] + sy[2] * sy[5] - sy[7])
                    - (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.)
            + (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.
            + sy[0]
                * (-sy[10] - 3. * sy[11] + sy[2] * sy[4] - 2. * sy[3] * sy[5]
                    + sy[1] * (3. * sy[6] + sy[2] * (-sy[8] - 2. * Zeta(2))))
            - 2. * sy[7] * Zeta(2)
            + sy[2]
                * (-sy[7] + sy[5] * sy[8] - G({ 0, 0, 0, b / c }, 1)
                    + 2. * sy[5] * Zeta(2))
            + sy[1]
                * (sy[9] - 2. * sy[3] * sy[8]
                    + 4. * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
                    - 4. * sy[3] * Zeta(2)
                    + sy[2]
                        * (sy[12] * sy[8] - G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                            + 2. * sy[12] * Zeta(2))) };
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        const vector<complex<double>> sy = { G({ 0, b, c }, { 1, sb, sc }, x),
            G({ 0, 0, 0, b / a }, 1), G({ 0, 0, b / a, c / a }, 1), Log(a, sa) };

        complex<double> res { -2. * sy[0] * G({ 0, 0, x / a }, 1)
            - 6. * G({ 0, 0, 0, 0, b / a, c / a }, 1)
            - 3. * G({ 0, 0, 0, b / a, 0, c / a }, 1)
            - 3. * G({ 0, 0, 0, b / a, c / a, x / a }, 1)
            - G({ 0, 0, b / a, 0, 0, c / a }, 1) - G({ 0, 0, b / a, 0, c / a, x / a }, 1)
            - G({ 0, 0, b / a, c / a, 0, x / a }, 1)
            + (-3. * sy[1] - sy[2]) * G({ 0, c }, { 1, sc }, x)
            + G({ 0, 0, b / a }, 1) * (-sy[0] + G({ 0, 0, c }, { 1, 1, sc }, x))
            - G({ 0, x / a }, 1) * G({ 0, 0, b, c }, { 1, 1, sb, sc }, x)
            + G({ 0, 0, 0, 0, b, c }, { 1, 1, 1, 1, sb, sc }, x)
            - sy[2] * sy[3] * Log(-x, sa) + (sy[2] * pow(sy[3], 2.)) / 2.
            + sy[2] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) }; // aBC
        if (b != x) {
            res += (3. * sy[1] - 3. * G({ 0, 0, 0, x / a }, 1))
                * G({ b, c }, { sb, sc }, x);
        }
        if (c != x) {
            res += (3. * G({ 0, 0, 0, b / a, c / a }, 1)
                       - 3. * G({ 0, 0, 0, b / a, x / a }, 1)
                       + G({ 0, 0, b / a, 0, c / a }, 1)
                       - G({ 0, 0, b / a, 0, x / a }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest
    if (b == c) { // Abb
        if (sb != sc) {
            throw FastGPL::FastGPL_error { "G6_0a00bc: b==c but sb!=sc" };
        }
        const vector<complex<double>> sy
            = { G({ 0, a, b }, { 1, sa, sb }, x), G({ 0, 0, a / b, 1 }, 1) };
        complex<double> res { G({ 0, 0, a / b, 0, 0, 1 }, 1)
            - G({ 0, 0, a / b, 0, x / b, 1 }, 1)
            + (-G({ 0, 0, a / b, 0, 1 }, 1) + G({ 0, 0, a / b, 0, x / b }, 1))
                * G({ b }, { sb }, x)
            + sy[1] * G({ 0, b }, { 1, sb }, x)
            + G({ 0, 0, a / b }, 1) * (sy[0] - G({ 0, 0, b }, { 1, 1, sb }, x))
            + G({ 0, a, 0, 0, 0, b }, { 1, sa, 1, 1, 1, sb }, x)
            - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x) * Zeta(2)
            + G({ 0, a }, { 1, sa }, x) * (-sy[1] - Zeta(4)) + sy[0] * Zeta(3) };
        return res;
    } else { // AbC
        const vector<complex<double>> sy
            = { Log(b, sb), G({ c / b }, 1), G({ 0, a }, { 1, sa }, x),
                  G({ 0, a, c }, { 1, sa, sc }, x), G({ 0, 0, a }, { 1, 1, sa }, x),
                  G({ 0, 0, a / b, c / b }, 1), G({ 0, 0, c / b, a / b }, 1),
                  G({ 0, c / b, 0, a / b }, 1), G({ c / b, 0, 0, a / b }, 1),
                  G({ 0, 0 }, { 1, 1 }, x), G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x),
                  G({ 0, 0, 0, a / b, c / b }, 1), G({ 0, 0, 0, c / b, a / b }, 1),
                  G({ 0, 0, c / b, 0, a / b }, 1), G({ 0, c / b, 0, 0, a / b }, 1),
                  G({ c / b, 0, 0, 0, a / b }, 1), G({ 0 }, { 1 }, x) };
        complex<double> res { sy[9] * (sy[5] + sy[6] + sy[7] + sy[8])
            - sy[3] * G({ 0, 0, c / b }, 1)
            + sy[2] * (sy[6] + sy[7] + sy[8] + G({ 0, 0, 0, c / b }, 1))
            - 6. * G({ 0, 0, 0, 0, a / b, c / b }, 1)
            - 6. * G({ 0, 0, 0, 0, c / b, a / b }, 1)
            - 6. * G({ 0, 0, 0, c / b, 0, a / b }, 1) + G({ 0, 0, a / b, 0, 0, c / b }, 1)
            + G({ 0, 0, a / b, 0, c / b, x / b }, 1)
            + G({ 0, 0, a / b, c / b, 0, x / b }, 1)
            - 6. * G({ 0, 0, c / b, 0, 0, a / b }, 1)
            + G({ 0, 0, c / b, a / b, 0, x / b }, 1)
            - 6. * G({ 0, c / b, 0, 0, 0, a / b }, 1)
            + G({ 0, c / b, 0, a / b, 0, x / b }, 1)
            - 6. * G({ c / b, 0, 0, 0, 0, a / b }, 1)
            + G({ c / b, 0, 0, a / b, 0, x / b }, 1)
            + G({ 0, 0, a / b }, 1) * (sy[3] - G({ 0, 0, c }, { 1, 1, sc }, x))
            + G({ 0, c / b }, 1) * G({ 0, a, 0, c }, { 1, sa, 1, sc }, x)
            + G({ 0, a, 0, 0, 0, c }, { 1, sa, 1, 1, 1, sc }, x)
            + (-(sy[1] * sy[4]) - sy[5] / 2. - sy[6] / 2. - sy[7] / 2. - sy[8] / 2.)
                * pow(sy[0], 2.)
            + Log(-x, sb)
                * (-3. * sy[11] - 3. * sy[12] - 3. * sy[13] - 3. * sy[14] - 3. * sy[15]
                    + 3. * sy[10] * sy[1]
                    + sy[0] * (2. * sy[1] * sy[4] + sy[5] + sy[6] + sy[7] + sy[8])
                    + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.)
            - (sy[1] * sy[2] * pow(sy[0], 3.)) / 6. + 2. * sy[6] * Zeta(2)
            + 2. * sy[7] * Zeta(2) + 2. * sy[8] * Zeta(2)
            + sy[5] * (G({ 0, c }, { 1, sc }, x) + 2. * Zeta(2))
            + sy[1]
                * (2. * sy[9] * sy[4] - 4. * G({ 0, 0, 0, 0, a }, { 1, 1, 1, 1, sa }, x)
                    - G({ 0, a, 0, 0, c }, { 1, sa, 1, 1, sc }, x) + 4. * sy[4] * Zeta(2)
                    + sy[2]
                        * (-(sy[9] * sy[16]) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                            - 2. * sy[16] * Zeta(2)))
            + sy[0]
                * (3. * sy[11] + 3. * sy[12] + 3. * sy[13] + 3. * sy[14] + 3. * sy[15]
                    + sy[1] * (-3. * sy[10] + sy[2] * (sy[9] + 2. * Zeta(2)))) };
        if (c != x) {
            res += (-G({ 0, 0, a / b, 0, c / b }, 1) + G({ 0, 0, a / b, 0, x / b }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }
}
complex<double> G6_00ab0c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy = { Log(c, sc), G({ 0, b / c, a / c }, 1),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, b / c }, 1), G({ 0, x / c }, 1),
            G({ 0, 0, b, a }, { 1, 1, sb, sa }, x),
            G({ 0, b, 0, a }, { 1, sb, 1, sa }, x), G({ 0, 0, b / c, a / c }, 1),
            G({ 0, b / c, 0, a / c }, 1), G({ 0, 0 }, { 1, 1 }, x),
            G({ 0, 0, 0, b / c, a / c }, 1), G({ 0, 0, b / c, 0, a / c }, 1),
            G({ 0, b / c, 0, 0, a / c }, 1), G({ 0 }, { 1 }, x) };
        complex<double> res { -(sy[4] * sy[5]) - sy[4] * sy[6]
            + sy[9] * (2. * sy[7] + sy[8]) - 4. * G({ 0, 0, 0, 0, b / c, a / c }, 1)
            - 3. * G({ 0, 0, 0, b / c, 0, a / c }, 1)
            - 2. * G({ 0, 0, b / c, 0, 0, a / c }, 1) - G({ 0, b / c, 0, 0, 0, a / c }, 1)
            - G({ 0, b / c, a / c, 0, 0, x / c }, 1)
            + sy[2]
                * (sy[1] - sy[0] * sy[3] + 2. * G({ 0, 0, b / c }, 1)
                    + G({ 0, 0, b }, { 1, 1, sb }, x))
            + sy[3] * (sy[5] + sy[6] - 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
            + G({ 0, c }, { 1, sc }, x) * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
            + G({ c }, { sc }, x)
                * (-3. * G({ 0, 0, 0, a, b }, { 1, 1, 1, sa, sb }, x)
                    - G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x))
            + G({ x / c }, 1)
                * (sy[2] * G({ 0, b }, { 1, sb }, x)
                    - 3. * G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x)
                    - 2. * G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x)
                    - G({ 0, b, 0, 0, a }, { 1, sb, 1, 1, sa }, x))
            - 6. * G({ 0, 0, 0, 0, b, a }, { 1, 1, 1, 1, sb, sa }, x)
            - 3. * G({ 0, 0, 0, b, 0, a }, { 1, 1, 1, sb, 1, sa }, x)
            - G({ 0, 0, b, 0, 0, a }, { 1, 1, sb, 1, 1, sa }, x)
            + (-sy[7] - sy[8] / 2.) * pow(sy[0], 2.)
            + Log(-x, sc)
                * (-3. * sy[10] - 2. * sy[11] - sy[12] + sy[2] * sy[3]
                    + sy[0] * (2. * sy[7] + sy[8]) - (sy[1] * pow(sy[0], 2.)) / 2.)
            + (sy[1] * pow(sy[0], 3.)) / 6.
            + sy[0]
                * (3. * sy[10] + 2. * sy[11] + sy[12] + sy[1] * (-sy[9] - 2. * Zeta(2)))
            + 4. * sy[7] * Zeta(2) + 2. * sy[8] * Zeta(2)
            + sy[1]
                * (sy[9] * sy[13] - G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                    + 2. * sy[13] * Zeta(2)) };
        if (b != x) {
            res += (-(sy[2] * sy[3]) + sy[2] * sy[4]) * G({ b }, { sb }, x)
                + (sy[3] - sy[4]) * G({ b, 0, 0, a }, { sb, 1, 1, sa }, x);
        }
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        if (a == b) // aaC
        {
            if (sa != sb) {
                throw FastGPL::FastGPL_error { "G6_00ab0c: a==b but sa!=sb" };
            }
            const vector<complex<double>> sy = { Log(a, sa), G({ 0, 1, c / a }, 1),
                G({ 0, c / a, 1 }, 1), G({ 0 }, { 1 }, x), G({ 0, 0 }, { 1, 1 }, x),
                G({ 0, 0, 0 }, { 1, 1, 1 }, x), G({ a, 0, c }, { sa, 1, sc }, x),
                G({ 0, 0, c }, { 1, 1, sc }, x) };
            complex<double> res { sy[2] * sy[3] * sy[4] - sy[2] * sy[5]
                + sy[6] * G({ 0, 0, x / a }, 1) - 4. * G({ 0, 0, 0, 0, 1, c / a }, 1)
                - 4. * G({ 0, 0, 0, 0, c / a, 1 }, 1)
                - 3. * G({ 0, 0, 0, 1, 0, c / a }, 1)
                - 3. * G({ 0, 0, 0, 1, c / a, x / a }, 1)
                - 3. * G({ 0, 0, 0, c / a, 1, x / a }, 1)
                - 3. * G({ 0, 0, 0, c / a, x / a, 1 }, 1)
                - 2. * G({ 0, 0, 1, 0, 0, c / a }, 1)
                - 2. * G({ 0, 0, 1, 0, c / a, x / a }, 1)
                - 2. * G({ 0, 0, 1, c / a, 0, x / a }, 1)
                - 2. * G({ 0, 0, c / a, 0, 1, x / a }, 1)
                - 2. * G({ 0, 0, c / a, 0, x / a, 1 }, 1)
                - 2. * G({ 0, 0, c / a, 1, 0, x / a }, 1) - G({ 0, 1, 0, 0, 0, c / a }, 1)
                - G({ 0, 1, 0, 0, c / a, x / a }, 1) - G({ 0, 1, 0, c / a, 0, x / a }, 1)
                - G({ 0, 1, c / a, 0, 0, x / a }, 1) - G({ 0, c / a, 0, 0, 1, x / a }, 1)
                - G({ 0, c / a, 0, 0, x / a, 1 }, 1) - G({ 0, c / a, 0, 1, 0, x / a }, 1)
                - G({ 0, c / a, 1, 0, 0, x / a }, 1)
                + (-2. * G({ 0, 0, 1, c / a }, 1) - 2. * G({ 0, 0, c / a, 1 }, 1)
                      - G({ 0, 0, x / a, 1 }, 1) - G({ 0, 1, 0, c / a }, 1))
                    * G({ 0, c }, { 1, sc }, x)
                + G({ 0, 0, 0, a, 0, c }, { 1, 1, 1, sa, 1, sc }, x)
                + (-sy[1] / 2. - sy[2] / 2.) * Log(-x, sa) * pow(sy[0], 2.)
                + (sy[1] / 6. + sy[2] / 6.) * pow(sy[0], 3.)
                + 2. * sy[2] * sy[3] * Zeta(2)
                + G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x) * Zeta(2)
                - G({ 0, a, 0, c }, { 1, sa, 1, sc }, x) * Zeta(2)
                + sy[0]
                    * (-(sy[2] * sy[4]) + sy[1] * (-sy[4] - 2. * Zeta(2))
                        - 2. * sy[2] * Zeta(2))
                + sy[1] * (sy[3] * sy[4] - sy[5] + sy[7] + 2. * sy[3] * Zeta(2))
                + sy[7] * (sy[2] - Zeta(3)) + sy[6] * Zeta(3) };
            if (c != x) {
                res += (3. * G({ 0, 0, 0, 1, c / a }, 1)
                           - 3. * G({ 0, 0, 0, 1, x / a }, 1)
                           + 3. * G({ 0, 0, 0, c / a, 1 }, 1)
                           - 3. * G({ 0, 0, 0, x / a, 1 }, 1)
                           + 2. * G({ 0, 0, 1, 0, c / a }, 1)
                           - 2. * G({ 0, 0, 1, 0, x / a }, 1)
                           + G({ 0, 1, 0, 0, c / a }, 1) - G({ 0, 1, 0, 0, x / a }, 1))
                    * G({ c }, { sc }, x);
            }
            return res;
        }
        const vector<complex<double>> sy = { G({ 0, 0, b / a }, 1),
            G({ 0, 0, c }, { 1, 1, sc }, x), Log(a, sa), G({ b / a, 0, c / a }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x) };

        complex<double> res { sy[0] * sy[1] - sy[1] * sy[3]
            + G({ 0, 0, 0, b / a, 0, c / a }, 1) + 2. * G({ 0, 0, b / a, 0, 0, c / a }, 1)
            + G({ 0, 0, b / a, 0, c / a, x / a }, 1)
            + 3. * G({ 0, b / a, 0, 0, 0, c / a }, 1)
            + 2. * G({ 0, b / a, 0, 0, c / a, x / a }, 1)
            + G({ 0, b / a, 0, c / a, 0, x / a }, 1)
            + 4. * G({ b / a, 0, 0, 0, 0, c / a }, 1)
            + 3. * G({ b / a, 0, 0, 0, c / a, x / a }, 1)
            + 2. * G({ b / a, 0, 0, c / a, 0, x / a }, 1)
            + G({ b / a, 0, c / a, 0, 0, x / a }, 1)
            + (G({ 0, 0, b / a, x / a }, 1) + G({ 0, b / a, 0, c / a }, 1)
                  + G({ 0, b / a, 0, x / a }, 1) + 2. * G({ b / a, 0, 0, c / a }, 1)
                  + G({ b / a, 0, 0, x / a }, 1))
                * G({ 0, c }, { 1, sc }, x)
            + G({ 0, b / a }, 1)
                * (-G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                    + G({ 0, b, 0, c }, { 1, sb, 1, sc }, x))
            + G({ b / a }, 1)
                * (G({ 0, 0, 0, 0, c }, { 1, 1, 1, 1, sc }, x)
                    - G({ 0, 0, b, 0, c }, { 1, 1, sb, 1, sc }, x))
            + G({ 0, 0, 0, b, 0, c }, { 1, 1, 1, sb, 1, sc }, x)
            + (sy[3] * Log(-x, sa) * pow(sy[2], 2.)) / 2. - (sy[3] * pow(sy[2], 3.)) / 6.
            + sy[2] * sy[3] * (sy[4] + 2. * Zeta(2))
            + sy[3]
                * (-(sy[4] * sy[5]) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                    - 2. * sy[5] * Zeta(2)) }; // aBC
        if (b != x) {
            res += (-sy[0] + G({ 0, 0, x / a }, 1)) * G({ b, 0, c }, { sb, 1, sc }, x);
        }
        if (c != x) {
            res += (-G({ 0, 0, b / a, 0, c / a }, 1) + G({ 0, 0, b / a, 0, x / a }, 1)
                       - 2. * G({ 0, b / a, 0, 0, c / a }, 1)
                       + 2. * G({ 0, b / a, 0, 0, x / a }, 1)
                       - 3. * G({ b / a, 0, 0, 0, c / a }, 1)
                       + 3. * G({ b / a, 0, 0, 0, x / a }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest
    const vector<complex<double>> sy = { G({ 0, 0, a }, { 1, 1, sa }, x),
        G({ 0, c / b, a / b }, 1), G({ 0, 0, 0 }, { 1, 1, 1 }, x),
        G({ a / b, 0, c / b }, 1), Log(b, sb), G({ 0, a / b, c / b }, 1),
        G({ 0, 0, c }, { 1, 1, sc }, x), G({ 0, c / b }, 1),
        G({ 0, 0, a, c }, { 1, 1, sa, sc }, x), G({ 0, 0, a / b, c / b }, 1),
        G({ 0, 0, c / b, a / b }, 1), G({ 0, a / b, 0, c / b }, 1),
        G({ 0, c / b, 0, a / b }, 1), G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x),
        G({ 0, c }, { 1, sc }, x), G({ 0, 0, 0, a / b, c / b }, 1),
        G({ 0, 0, 0, c / b, a / b }, 1), G({ 0, 0, a / b, 0, c / b }, 1),
        G({ 0, 0, c / b, 0, a / b }, 1), G({ 0, c / b, 0, 0, a / b }, 1) };
    complex<double> res { sy[2] * (-sy[1] - sy[3])
        + sy[13] * (2. * sy[9] + 2. * sy[10] + sy[11] + sy[12] + sy[14] * (sy[1] + sy[3]))
        + sy[3] * sy[6] + sy[0] * (sy[1] + 2. * G({ 0, 0, c / b }, 1))
        + sy[15] * (-2. * G({ a / b, 0, 0, c / b }, 1) - G({ a / b, 0, 0, x / b }, 1))
        - 4. * G({ 0, 0, 0, 0, a / b, c / b }, 1)
        - 4. * G({ 0, 0, 0, 0, c / b, a / b }, 1) - G({ 0, 0, 0, a / b, 0, c / b }, 1)
        - 3. * G({ 0, 0, 0, c / b, 0, a / b }, 1)
        - 2. * G({ 0, 0, c / b, 0, 0, a / b }, 1) - G({ 0, a / b, 0, 0, 0, c / b }, 1)
        - G({ 0, a / b, 0, 0, c / b, x / b }, 1) - G({ 0, a / b, 0, c / b, 0, x / b }, 1)
        - G({ 0, a / b, c / b, 0, 0, x / b }, 1) - G({ 0, c / b, 0, 0, 0, a / b }, 1)
        - G({ 0, c / b, a / b, 0, 0, x / b }, 1) - 4. * G({ a / b, 0, 0, 0, 0, c / b }, 1)
        - 3. * G({ a / b, 0, 0, 0, c / b, x / b }, 1)
        - 2. * G({ a / b, 0, 0, c / b, 0, x / b }, 1)
        - G({ a / b, 0, c / b, 0, 0, x / b }, 1)
        + sy[7] * (-sy[8] - 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
        + G({ 0, a / b }, 1) * (sy[8] - G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x))
        + G({ a / b }, 1)
            * (-G({ 0, 0, 0, 0, c }, { 1, 1, 1, 1, sc }, x)
                + G({ 0, 0, a, 0, c }, { 1, 1, sa, 1, sc }, x))
        + G({ 0, 0, a, 0, 0, c }, { 1, 1, sa, 1, 1, sc }, x)
        + (-sy[9] - sy[10] - sy[11] / 2. - sy[12] / 2.) * pow(sy[4], 2.)
        + Log(-x, sb)
            * (-3. * sy[16] - 3. * sy[17] - sy[18] - 2. * sy[19] - sy[20]
                + (2. * sy[9] + 2. * sy[10] + sy[11] + sy[12]) * sy[4] + sy[0] * sy[7]
                + (-sy[1] / 2. - sy[3] / 2. - sy[5] / 2.) * pow(sy[4], 2.))
        + (sy[1] / 6. + sy[3] / 6. + sy[5] / 6.) * pow(sy[4], 3.) + 4. * sy[9] * Zeta(2)
        + 4. * sy[10] * Zeta(2) + 2. * sy[12] * Zeta(2)
        + sy[11] * (-sy[15] + 2. * Zeta(2))
        + sy[5] * (sy[13] * sy[14] - sy[2] + sy[6] + 2. * sy[14] * Zeta(2))
        + sy[4]
            * (3. * sy[16] + 3. * sy[17] + sy[18] + 2. * sy[19] + sy[20]
                + sy[13] * (-sy[1] - sy[3]) - sy[0] * sy[7]
                + sy[5] * (-sy[13] - 2. * Zeta(2)) - 2. * sy[1] * Zeta(2)
                - 2. * sy[3] * Zeta(2))
        + sy[14] * (2. * sy[1] * Zeta(2) + 2. * sy[3] * Zeta(2)) };
    if (c != x) {
        res += (G({ 0, a / b, 0, 0, c / b }, 1) - G({ 0, a / b, 0, 0, x / b }, 1)
                   + 3. * G({ a / b, 0, 0, 0, c / b }, 1)
                   - 3. * G({ a / b, 0, 0, 0, x / b }, 1))
            * G({ c }, { sc }, x);
    }
    return res;
}
complex<double> G6_0ab00c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy = { G({ 0, 0, b / c }, 1),
            G({ 0, b, a }, { 1, sb, sa }, x), G({ 0, 0, x / c }, 1),
            G({ 0, a }, { 1, sa }, x), G({ 0, x / c }, 1), G({ c }, { sc }, x),
            G({ 0, 0, b }, { 1, 1, sb }, x), G({ x / c }, 1), Log(c, sc),
            G({ 0, 0, b / c, a / c }, 1), G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x),
            G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x), G({ 0, 0, 0, b / c, a / c }, 1),
            G({ 0, 0, b / c, 0, a / c }, 1) };
        complex<double> res { sy[1] * sy[2] + (-3. * sy[10] - sy[11]) * sy[5]
            + (3. * sy[10] + sy[11]) * sy[7] + (-3. * sy[12] - sy[13]) * sy[8]
            + 6. * G({ 0, 0, 0, 0, b / c, a / c }, 1)
            + 3. * G({ 0, 0, 0, b / c, 0, a / c }, 1) + G({ 0, 0, b / c, 0, 0, a / c }, 1)
            - G({ 0, 0, b / c, a / c, 0, x / c }, 1)
            + sy[0] * (-sy[1] + 2. * G({ 0, 0, a }, { 1, 1, sa }, x))
            + G({ 0, 0, c }, { 1, 1, sc }, x) * G({ 0, a, b }, { 1, sa, sb }, x)
            + sy[3]
                * (-sy[9] + sy[5] * sy[6] - sy[6] * sy[7] + sy[0] * sy[8]
                    - 3. * G({ 0, 0, 0, b / c }, 1) - sy[4] * G({ 0, b }, { 1, sb }, x)
                    - G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x))
            + G({ 0, c }, { 1, sc }, x)
                * (-2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                    - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
            + sy[4]
                * (2. * G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                    + G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
            + 4. * G({ 0, 0, 0, 0, b, a }, { 1, 1, 1, 1, sb, sa }, x)
            + G({ 0, 0, 0, b, 0, a }, { 1, 1, 1, sb, 1, sa }, x)
            + (3. * sy[12] + sy[13] - sy[0] * sy[3] - sy[9] * sy[8]) * Log(-x, sc)
            + (sy[9] * pow(sy[8], 2.)) / 2.
            + sy[9] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
        if (b != x) {
            res += (sy[0] * sy[3] - sy[2] * sy[3]) * G({ b }, { sb }, x)
                + (-sy[0] + sy[2]) * G({ b, 0, a }, { sb, 1, sa }, x);
        }
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        if (a == b) // aaC
        {
            if (sa != sb) {
                throw FastGPL::FastGPL_error { "G6_0ab00c: a==b but sa!=sb" };
            }
            const vector<complex<double>> sy
                = { Log(a, sa), G({ 0, 0, 1, c / a }, 1), G({ 0, 0, c / a, 1 }, 1),
                      G({ 0, 1, 0, c / a }, 1), G({ 0, 0 }, { 1, 1 }, x),
                      G({ 0, c }, { 1, sc }, x), G({ a, 0, 0, c }, { sa, 1, 1, sc }, x) };
            complex<double> res { (sy[2] + sy[3]) * sy[4] - sy[6] * G({ 0, x / a }, 1)
                + sy[5]
                    * (sy[2] + sy[3] + 2. * G({ 0, 0, 1, x / a }, 1)
                        + 2. * G({ 0, 0, x / a, 1 }, 1) + G({ 0, 1, 0, x / a }, 1))
                + 6. * G({ 0, 0, 0, 0, 1, c / a }, 1)
                + 6. * G({ 0, 0, 0, 0, c / a, 1 }, 1)
                + 6. * G({ 0, 0, 0, 1, 0, c / a }, 1)
                + 3. * G({ 0, 0, 0, 1, c / a, x / a }, 1)
                + 3. * G({ 0, 0, 0, c / a, 1, x / a }, 1)
                + 3. * G({ 0, 0, 0, c / a, x / a, 1 }, 1)
                + 5. * G({ 0, 0, 1, 0, 0, c / a }, 1)
                + 3. * G({ 0, 0, 1, 0, c / a, x / a }, 1)
                + G({ 0, 0, 1, c / a, 0, x / a }, 1) + G({ 0, 0, c / a, 0, 1, x / a }, 1)
                + G({ 0, 0, c / a, 0, x / a, 1 }, 1) + G({ 0, 0, c / a, 1, 0, x / a }, 1)
                + 3. * G({ 0, 1, 0, 0, 0, c / a }, 1)
                + 2. * G({ 0, 1, 0, 0, c / a, x / a }, 1)
                + G({ 0, 1, 0, c / a, 0, x / a }, 1)
                + G({ 0, x / a, 1 }, 1) * G({ 0, 0, c }, { 1, 1, sc }, x)
                + G({ 0, 0, a, 0, 0, c }, { 1, 1, sa, 1, 1, sc }, x)
                + sy[0] * (sy[1] + sy[2] + sy[3]) * Log(-x, sa)
                + (-sy[1] / 2. - sy[2] / 2. - sy[3] / 2.) * pow(sy[0], 2.)
                + 2. * sy[2] * Zeta(2) + 2. * sy[3] * Zeta(2) - sy[6] * Zeta(2)
                + G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x) * Zeta(2)
                + sy[1] * (sy[4] + sy[5] + 2. * Zeta(2)) };
            if (c != x) {
                res += (-3. * G({ 0, 0, 0, 1, c / a }, 1)
                           + 3. * G({ 0, 0, 0, 1, x / a }, 1)
                           - 3. * G({ 0, 0, 0, c / a, 1 }, 1)
                           + 3. * G({ 0, 0, 0, x / a, 1 }, 1)
                           - 3. * G({ 0, 0, 1, 0, c / a }, 1)
                           + 3. * G({ 0, 0, 1, 0, x / a }, 1)
                           - 2. * G({ 0, 1, 0, 0, c / a }, 1)
                           + 2. * G({ 0, 1, 0, 0, x / a }, 1))
                    * G({ c }, { sc }, x);
            }
            return res;
        }

        const vector<complex<double>> sy
            = { G({ 0, b / a }, 1), Log(a, sa), G({ b / a, 0, 0, c / a }, 1) };
        complex<double> res { -G({ 0, 0, b / a, 0, 0, c / a }, 1)
            - 3. * G({ 0, b / a, 0, 0, 0, c / a }, 1)
            - G({ 0, b / a, 0, 0, c / a, x / a }, 1)
            - 6. * G({ b / a, 0, 0, 0, 0, c / a }, 1)
            - 3. * G({ b / a, 0, 0, 0, c / a, x / a }, 1)
            - G({ b / a, 0, 0, c / a, 0, x / a }, 1)
            + (-sy[2] - G({ 0, b / a, 0, x / a }, 1) - 2. * G({ b / a, 0, 0, x / a }, 1))
                * G({ 0, c }, { 1, sc }, x)
            + (-G({ 0, b / a, x / a }, 1) - G({ b / a, 0, x / a }, 1))
                * G({ 0, 0, c }, { 1, 1, sc }, x)
            - sy[0] * G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
            + G({ b / a }, 1)
                * (G({ 0, 0, 0, 0, c }, { 1, 1, 1, 1, sc }, x)
                    - G({ 0, b, 0, 0, c }, { 1, sb, 1, 1, sc }, x))
            + G({ 0, 0, b, 0, 0, c }, { 1, 1, sb, 1, 1, sc }, x)
            - sy[1] * sy[2] * Log(-x, sa) + (sy[2] * pow(sy[1], 2.)) / 2.
            + sy[2] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) }; // aBC
        if (b != x) {
            res += (sy[0] - G({ 0, x / a }, 1)) * G({ b, 0, 0, c }, { sb, 1, 1, sc }, x);
        }
        if (c != x) {
            res += (G({ 0, b / a, 0, 0, c / a }, 1) - G({ 0, b / a, 0, 0, x / a }, 1)
                       + 3. * G({ b / a, 0, 0, 0, c / a }, 1)
                       - 3. * G({ b / a, 0, 0, 0, x / a }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest
    const vector<complex<double>> sy = { G({ 0, 0, c / b }, 1),
        G({ 0, a, c }, { 1, sa, sc }, x), G({ 0, 0, c }, { 1, 1, sc }, x),
        G({ 0, a }, { 1, sa }, x), G({ 0, 0, c / b, a / b }, 1), Log(b, sb),
        G({ 0, 0, a / b, c / b }, 1), G({ 0, a / b, 0, c / b }, 1),
        G({ a / b, 0, 0, c / b }, 1), G({ 0, 0 }, { 1, 1 }, x), G({ 0, c }, { 1, sc }, x),
        G({ 0, 0, 0, a / b, c / b }, 1), G({ 0, 0, 0, c / b, a / b }, 1),
        G({ 0, 0, a / b, 0, c / b }, 1), G({ 0, 0, c / b, 0, a / b }, 1),
        G({ 0, a / b, 0, 0, c / b }, 1) };
    complex<double> res { (3. * sy[11] + 3. * sy[12] + 2. * sy[13] + sy[14] + sy[15]
                              - sy[0] * sy[3])
            * sy[5]
        + sy[9] * (sy[4] + sy[7] + sy[8]) + (sy[1] - sy[2]) * G({ 0, 0, a / b }, 1)
        + sy[2] * G({ a / b, 0, x / b }, 1)
        + sy[3] * (sy[4] + 3. * G({ 0, 0, 0, c / b }, 1))
        + sy[10]
            * (sy[7] + sy[8] + G({ 0, a / b, 0, x / b }, 1)
                + 2. * G({ a / b, 0, 0, x / b }, 1))
        - 6. * G({ 0, 0, 0, 0, a / b, c / b }, 1)
        - 6. * G({ 0, 0, 0, 0, c / b, a / b }, 1)
        - 3. * G({ 0, 0, 0, a / b, 0, c / b }, 1)
        - 3. * G({ 0, 0, 0, c / b, 0, a / b }, 1) + G({ 0, 0, a / b, 0, c / b, x / b }, 1)
        + G({ 0, 0, a / b, c / b, 0, x / b }, 1) - G({ 0, 0, c / b, 0, 0, a / b }, 1)
        + G({ 0, 0, c / b, a / b, 0, x / b }, 1) + 3. * G({ 0, a / b, 0, 0, 0, c / b }, 1)
        + 2. * G({ 0, a / b, 0, 0, c / b, x / b }, 1)
        + G({ 0, a / b, 0, c / b, 0, x / b }, 1) + 6. * G({ a / b, 0, 0, 0, 0, c / b }, 1)
        + 3. * G({ a / b, 0, 0, 0, c / b, x / b }, 1)
        + G({ a / b, 0, 0, c / b, 0, x / b }, 1)
        + sy[0] * (-sy[1] - 2. * G({ 0, 0, a }, { 1, 1, sa }, x))
        + G({ 0, a / b }, 1)
            * (-G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                + G({ 0, a, 0, c }, { 1, sa, 1, sc }, x))
        + G({ a / b }, 1)
            * (-G({ 0, 0, 0, 0, c }, { 1, 1, 1, 1, sc }, x)
                + G({ 0, a, 0, 0, c }, { 1, sa, 1, 1, sc }, x))
        + G({ 0, a, 0, 0, 0, c }, { 1, sa, 1, 1, 1, sc }, x)
        + (-3. * sy[11] - 3. * sy[12] - 2. * sy[13] - sy[14] - sy[15] + sy[0] * sy[3]
              + sy[5] * (sy[4] + sy[6] + sy[7] + sy[8]))
            * Log(-x, sb)
        + (-sy[4] / 2. - sy[6] / 2. - sy[7] / 2. - sy[8] / 2.) * pow(sy[5], 2.)
        + 2. * sy[4] * Zeta(2) + 2. * sy[7] * Zeta(2) + 2. * sy[8] * Zeta(2)
        + sy[6] * (sy[9] + sy[10] + 2. * Zeta(2)) };
    if (c != x) {
        res += (-sy[13] - 2. * sy[15] + G({ 0, 0, a / b, 0, x / b }, 1)
                   + 2. * G({ 0, a / b, 0, 0, x / b }, 1)
                   - 3. * G({ a / b, 0, 0, 0, c / b }, 1)
                   + 3. * G({ a / b, 0, 0, 0, x / b }, 1))
            * G({ c }, { sc }, x);
    }
    return res;
}
complex<double> G6_a00b0c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy
            = { Log(c, sc), G({ a }, { sa }, x), G({ 0, b / c }, 1), G({ x / c }, 1),
                  G({ 0, b }, { 1, sb }, x), G({ 0, 0, a }, { 1, 1, sa }, x),
                  G({ 0, 0, b }, { 1, 1, sb }, x), G({ 0, 0, b / c }, 1),
                  G({ 0, a }, { 1, sa }, x), G({ 0, x / c }, 1), G({ c }, { sc }, x),
                  G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x), G({ 0, 0, 0, b / c }, 1),
                  G({ 0, c }, { 1, sc }, x), G({ 0, b / c, 0, 0, a / c }, 1),
                  G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x) };
        complex<double> res { sy[3] * sy[4] * sy[5] + sy[5] * (sy[6] + 2. * sy[7])
            + (-3. * sy[11] - 3. * sy[12] - sy[9] * sy[4] + 2. * sy[10] * sy[6]
                  - 2. * sy[3] * sy[6])
                * sy[8]
            - 2. * G({ 0, 0, b / c, 0, 0, a / c }, 1)
            - 3. * G({ 0, b / c, 0, 0, 0, a / c }, 1)
            - G({ 0, b / c, 0, 0, a / c, x / c }, 1)
            + sy[13]
                * (-G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                    - G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                    - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
            + sy[10]
                * (-3. * G({ 0, 0, 0, a, b }, { 1, 1, 1, sa, sb }, x)
                    - 3. * G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x)
                    - G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x)
                    - 2. * G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x))
            - sy[3] * G({ 0, b, 0, 0, a }, { 1, sb, 1, 1, sa }, x)
            - G({ 0, 0, b, 0, 0, a }, { 1, 1, sb, 1, 1, sa }, x)
            + (sy[1] * sy[7] - (sy[2] * sy[8]) / 2.) * pow(sy[0], 2.)
            + Log(-x, sc)
                * (-sy[14] + 3. * sy[12] * sy[1] + sy[2] * sy[5] - 2. * sy[7] * sy[8]
                    + sy[0] * (-2. * sy[1] * sy[7] + sy[2] * sy[8])
                    + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.)
            - (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.
            + sy[2]
                * (sy[15] * sy[8] - G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                    + sy[8] * (sy[4] + 2. * Zeta(2)))
            + sy[1]
                * (-3. * sy[10] * sy[11] + sy[14] + 3. * sy[11] * sy[3] + sy[9] * sy[6]
                    + sy[13] * sy[6] - 2. * sy[15] * sy[7]
                    + 4. * G({ 0, 0, 0, 0, b / c }, 1)
                    + 6. * G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x)
                    - 4. * sy[7] * Zeta(2)
                    + sy[2]
                        * (-(sy[15] * sy[16]) - sy[6] + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                            - 2. * sy[16] * Zeta(2)))
            + sy[0]
                * (sy[14] - sy[2] * sy[5] + 2. * sy[7] * sy[8]
                    + sy[1] * (-3. * sy[12] + sy[2] * (sy[15] + 2. * Zeta(2)))) };
        if (b != x) {
            res += (sy[9] * sy[5] - sy[2] * sy[5]) * G({ b }, { sb }, x)
                + (-sy[9] + sy[2]) * G({ b, 0, 0, a }, { sb, 1, 1, sa }, x);
        }
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        const vector<complex<double>> sy
            = { G({ 0, 0, b / a }, 1), G({ 0, 0, b / a, 0, c / a }, 1) };
        complex<double> res { 3. * G({ 0, 0, 0, b / a, 0, c / a }, 1)
            + 2. * G({ 0, 0, b / a, 0, 0, c / a }, 1)
            + G({ 0, 0, b / a, 0, c / a, x / a }, 1)
            + G({ 0, 0, b / a, x / a }, 1) * G({ 0, c }, { 1, sc }, x)
            + sy[0] * G({ 0, 0, c }, { 1, 1, sc }, x)
            + G({ 0, x / a }, 1) * G({ 0, b, 0, c }, { 1, sb, 1, sc }, x)
            + G({ x / a }, 1) * G({ 0, 0, b, 0, c }, { 1, 1, sb, 1, sc }, x)
            + G({ 0, 0, 0, b, 0, c }, { 1, 1, 1, sb, 1, sc }, x) - sy[1] * Log(a, sa)
            + sy[1] * Log(-x, sa) }; // aBC
        if (b != x) {
            res += (-sy[0] + G({ 0, 0, x / a }, 1)) * G({ b, 0, c }, { sb, 1, sc }, x);
        }
        if (c != x) {
            res += (-sy[1] + G({ 0, 0, b / a, 0, x / a }, 1)) * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest

    const vector<complex<double>> sy = { Log(b, sb), G({ a }, { sa }, x),
        G({ 0, c / b }, 1), G({ 0, a }, { 1, sa }, x), G({ a, 0, c }, { sa, 1, sc }, x),
        G({ a, c }, { sa, sc }, x), G({ 0, 0, 0, a / b }, 1),
        G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, 0, 0, a / b, c / b }, 1),
        G({ 0, 0, 0, c / b, a / b }, 1), G({ 0, 0, a / b, 0, c / b }, 1),
        G({ 0, 0, c / b, 0, a / b }, 1), G({ 0, c / b, 0, 0, a / b }, 1),
        G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x) };
    complex<double> res { 3. * sy[5] * sy[6] + 2. * sy[4] * G({ 0, 0, c / b }, 1)
        - 3. * sy[5] * G({ 0, 0, 0, c / b }, 1) - 12. * G({ 0, 0, 0, 0, a / b, c / b }, 1)
        - 12. * G({ 0, 0, 0, 0, c / b, a / b }, 1)
        - 6. * G({ 0, 0, 0, a / b, 0, c / b }, 1)
        - 3. * G({ 0, 0, 0, a / b, c / b, x / b }, 1)
        - 9. * G({ 0, 0, 0, c / b, 0, a / b }, 1)
        - 3. * G({ 0, 0, 0, c / b, a / b, x / b }, 1)
        - 2. * G({ 0, 0, a / b, 0, 0, c / b }, 1) - G({ 0, 0, a / b, 0, c / b, x / b }, 1)
        - 6. * G({ 0, 0, c / b, 0, 0, a / b }, 1)
        - 2. * G({ 0, 0, c / b, 0, a / b, x / b }, 1)
        - 3. * G({ 0, c / b, 0, 0, 0, a / b }, 1) - G({ 0, c / b, 0, 0, a / b, x / b }, 1)
        + (-3. * sy[6] - G({ 0, 0, a / b, x / b }, 1)) * G({ 0, c }, { 1, sc }, x)
        + G({ 0, 0, a / b }, 1) * (sy[4] - G({ 0, 0, c }, { 1, 1, sc }, x))
        + G({ a, 0, 0, 0, 0, c }, { sa, 1, 1, 1, 1, sc }, x)
        - (sy[2] * sy[3] * pow(sy[0], 2.)) / 2.
        + Log(-x, sb)
            * (-3. * sy[9] - sy[10] - 2. * sy[11] - sy[12] + sy[0] * sy[2] * sy[3]
                + sy[2] * sy[7] - 3. * sy[8] + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.)
        - (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.
        + sy[2]
            * (sy[13] * sy[3] - G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                - G({ a, 0, 0, c }, { sa, 1, 1, sc }, x) + 2. * sy[3] * Zeta(2))
        + sy[0]
            * (3. * sy[9] + sy[10] + 2. * sy[11] + sy[12] - sy[2] * sy[7] + 3. * sy[8]
                + sy[1] * sy[2] * (sy[13] + 2. * Zeta(2)))
        + sy[1]
            * (3. * sy[9] + 2. * sy[11] + sy[12] + 4. * G({ 0, 0, 0, 0, c / b }, 1)
                + sy[2]
                    * (-(sy[13] * sy[14]) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                        - 2. * sy[14] * Zeta(2))) };
    if (c != x) {
        res += (sy[10] + 3. * sy[8] - 3. * G({ 0, 0, 0, a / b, x / b }, 1)
                   - G({ 0, 0, a / b, 0, x / b }, 1))
            * G({ c }, { sc }, x);
    }
    return res;
}
complex<double> G6_a0b00c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy = { G({ a }, { sa }, x),
            G({ 0, c }, { 1, sc }, x), G({ 0, 0, b }, { 1, 1, sb }, x),
            G({ 0, 0, b / c }, 1), G({ 0, b }, { 1, sb }, x),
            G({ 0, 0, c }, { 1, 1, sc }, x), G({ 0, 0, x / c }, 1),
            G({ 0, a }, { 1, sa }, x), G({ 0, x / c }, 1), G({ c }, { sc }, x),
            G({ x / c }, 1), Log(c, sc), G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x),
            G({ 0, 0, 0, b / c }, 1), G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x),
            G({ 0, 0, b / c, 0, a / c }, 1) };
        complex<double> res { sy[10] * (-3. * sy[0] * sy[12] - sy[14])
            + sy[9] * (3. * sy[0] * sy[12] + sy[14]) + sy[11] * sy[15]
            - 2. * sy[0] * sy[1] * sy[2] + sy[0] * sy[4] * (sy[3] + sy[5] - sy[6])
            + sy[7]
                * (sy[12] + 3. * sy[13] - sy[9] * sy[2] + sy[10] * sy[2] - sy[11] * sy[3]
                    + sy[4] * sy[8])
            - 3. * G({ 0, 0, 0, b / c, 0, a / c }, 1)
            - 2. * G({ 0, 0, b / c, 0, 0, a / c }, 1)
            - G({ 0, 0, b / c, 0, a / c, x / c }, 1)
            - sy[3] * G({ 0, 0, a }, { 1, 1, sa }, x)
            + sy[5]
                * (-G({ 0, a, b }, { 1, sa, sb }, x) - G({ 0, b, a }, { 1, sb, sa }, x))
            + sy[1]
                * (2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                    + 2. * G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                    + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
            + sy[8] * (-2. * sy[0] * sy[2] - G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
            - G({ 0, 0, 0, b, 0, a }, { 1, 1, 1, sb, 1, sa }, x)
            + (-sy[15] + sy[0] * (-3. * sy[13] + sy[11] * sy[3]) + sy[3] * sy[7])
                * Log(-x, sc)
            + sy[0]
                * (3. * sy[11] * sy[13] + sy[15] - 6. * G({ 0, 0, 0, 0, b / c }, 1)
                    - 4. * G({ 0, 0, 0, 0, b }, { 1, 1, 1, 1, sb }, x)
                    - (sy[3] * pow(sy[11], 2.)) / 2.
                    + sy[3] * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2))) };
        if (b != x) {
            res += (-(sy[3] * sy[7]) + sy[6] * sy[7]) * G({ b }, { sb }, x)
                + (sy[3] - sy[6]) * G({ b, 0, a }, { sb, 1, sa }, x);
        }
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        const vector<complex<double>> sy
            = { G({ 0, b / a }, 1), G({ 0, b / a, 0, 0, c / a }, 1) };
        complex<double> res { 2. * G({ 0, 0, b / a, 0, 0, c / a }, 1)
            + 3. * G({ 0, b / a, 0, 0, 0, c / a }, 1)
            + G({ 0, b / a, 0, 0, c / a, x / a }, 1)
            + G({ 0, b / a, 0, x / a }, 1) * G({ 0, c }, { 1, sc }, x)
            + G({ 0, b / a, x / a }, 1) * G({ 0, 0, c }, { 1, 1, sc }, x)
            + sy[0] * G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
            + G({ x / a }, 1) * G({ 0, b, 0, 0, c }, { 1, sb, 1, 1, sc }, x)
            + G({ 0, 0, b, 0, 0, c }, { 1, 1, sb, 1, 1, sc }, x) - sy[1] * Log(a, sa)
            + sy[1] * Log(-x, sa) }; // aBC
        if (b != x) {
            res += (-sy[0] + G({ 0, x / a }, 1)) * G({ b, 0, 0, c }, { sb, 1, 1, sc }, x);
        }
        if (c != x) {
            res += (-sy[1] + G({ 0, b / a, 0, 0, x / a }, 1)) * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest

    const vector<complex<double>> sy
        = { Log(b, sb), G({ a }, { sa }, x), G({ 0, 0, c / b }, 1),
              G({ 0, 0, c }, { 1, 1, sc }, x), G({ a, 0, c }, { sa, 1, sc }, x),
              G({ a, c }, { sa, sc }, x), G({ 0, 0, 0, a / b }, 1),
              G({ 0, a }, { 1, sa }, x), G({ 0, 0, 0, a / b, c / b }, 1),
              G({ 0, 0, 0, c / b, a / b }, 1), G({ 0, 0, a / b, 0, c / b }, 1),
              G({ 0, 0, c / b, 0, a / b }, 1), G({ 0, a / b, 0, 0, c / b }, 1) };
    complex<double> res { -3. * sy[5] * sy[6]
        + sy[0]
            * (-3. * sy[9] - 2. * sy[10] - sy[11] - sy[12] + sy[2] * sy[7] - 3. * sy[8])
        + (2. * sy[3] - 2. * sy[4]) * G({ 0, 0, a / b }, 1)
        + sy[3] * G({ 0, a / b, x / b }, 1) + 3. * sy[5] * G({ 0, 0, 0, c / b }, 1)
        + 12. * G({ 0, 0, 0, 0, a / b, c / b }, 1)
        + 12. * G({ 0, 0, 0, 0, c / b, a / b }, 1)
        + 9. * G({ 0, 0, 0, a / b, 0, c / b }, 1)
        + 3. * G({ 0, 0, 0, a / b, c / b, x / b }, 1)
        + 6. * G({ 0, 0, 0, c / b, 0, a / b }, 1)
        + 3. * G({ 0, 0, 0, c / b, a / b, x / b }, 1)
        + 6. * G({ 0, 0, a / b, 0, 0, c / b }, 1)
        + 2. * G({ 0, 0, a / b, 0, c / b, x / b }, 1)
        + 2. * G({ 0, 0, c / b, 0, 0, a / b }, 1) + G({ 0, 0, c / b, 0, a / b, x / b }, 1)
        + 3. * G({ 0, a / b, 0, 0, 0, c / b }, 1) + G({ 0, a / b, 0, 0, c / b, x / b }, 1)
        + (3. * sy[6] + 2. * G({ 0, 0, a / b, x / b }, 1) + G({ 0, a / b, 0, x / b }, 1))
            * G({ 0, c }, { 1, sc }, x)
        + sy[2] * (-sy[4] + G({ 0, 0, a }, { 1, 1, sa }, x))
        + G({ 0, a / b }, 1)
            * (G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x)
                - G({ a, 0, 0, c }, { sa, 1, 1, sc }, x))
        + G({ a, 0, 0, 0, 0, c }, { sa, 1, 1, 1, 1, sc }, x)
        + (3. * sy[9] + 2. * sy[10] + sy[11] + sy[12] - sy[0] * sy[1] * sy[2]
              - sy[2] * sy[7] + 3. * sy[8])
            * Log(-x, sb)
        + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.
        + sy[1]
            * (-3. * sy[9] - sy[11] - 6. * G({ 0, 0, 0, 0, c / b }, 1)
                + sy[2] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
    if (c != x) {
        res += (-2. * sy[10] - sy[12] - 3. * sy[8] + 3. * G({ 0, 0, 0, a / b, x / b }, 1)
                   + 2. * G({ 0, 0, a / b, 0, x / b }, 1)
                   + G({ 0, a / b, 0, 0, x / b }, 1))
            * G({ c }, { sc }, x);
    }
    return res;
}
complex<double> G6_0a0b0c(complex<double> a1, complex<double> b1, complex<double> c1,
    int sa, int sb, int sc, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    double x = 1.;
    if (abs(c) < abs(b) && abs(c) < abs(a)) { // c is smallest
        const vector<complex<double>> sy = { G({ x / c }, 1), G({ 0, b }, { 1, sb }, x),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, 0, b }, { 1, 1, sb }, x),
            G({ 0, 0, b / c }, 1), G({ 0, b / c }, 1),
            G({ 0, b, 0, a }, { 1, sb, 1, sa }, x), G({ 0, x / c }, 1), Log(c, sc),
            G({ 0, a }, { 1, sa }, x), G({ 0, b / c, 0, a / c }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ c }, { sc }, x),
            G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x), G({ 0, 0, b / c, 0, a / c }, 1),
            G({ 0, b / c, 0, 0, a / c }, 1) };
        complex<double> res { -(sy[10] * sy[11]) - 2. * sy[0] * sy[1] * sy[2]
            + sy[2] * (-2. * sy[3] - 4. * sy[4]) + sy[6] * sy[7]
            + (-2. * sy[14] - 2. * sy[15] - 2. * sy[9] * sy[4] + 2. * sy[2] * sy[5])
                * sy[8]
            + 3. * G({ 0, 0, 0, b / c, 0, a / c }, 1)
            + 4. * G({ 0, 0, b / c, 0, 0, a / c }, 1)
            + 3. * G({ 0, b / c, 0, 0, 0, a / c }, 1)
            - G({ 0, b / c, 0, a / c, 0, x / c }, 1)
            + sy[5] * (-sy[6] + 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
            + G({ 0, c }, { 1, sc }, x) * G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
            + sy[12]
                * (2. * sy[13] + 6. * G({ 0, 0, 0, a, b }, { 1, 1, 1, sa, sb }, x)
                    + 6. * G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x)
                    + 2. * G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x))
            + sy[0] * (2. * sy[13] + 2. * G({ 0, b, 0, 0, a }, { 1, sb, 1, 1, sa }, x))
            + 3. * G({ 0, 0, 0, b, 0, a }, { 1, 1, 1, sb, 1, sa }, x)
            + 2. * G({ 0, 0, b, 0, 0, a }, { 1, 1, sb, 1, 1, sa }, x)
            + (2. * sy[14] + 2. * sy[15] + 2. * sy[9] * sy[4] - 2. * sy[2] * sy[5]
                  + (-sy[10] - sy[9] * sy[5]) * sy[8])
                * Log(-x, sc)
            + (sy[10] / 2. + (sy[9] * sy[5]) / 2.) * pow(sy[8], 2.)
            + sy[9]
                * (-sy[10] + 2. * sy[0] * sy[3] - 2. * sy[12] * sy[3] + sy[1] * sy[7]
                    + 3. * G({ 0, 0, 0, b / c }, 1)
                    + 3. * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
                    + sy[5] * (-sy[11] - sy[1] - 2. * Zeta(2)))
            - 2. * sy[10] * Zeta(2) };
        if (b != x) {
            res += (2. * sy[2] * sy[5] - 2. * sy[2] * sy[7]) * G({ b }, { sb }, x)
                + (-2. * sy[5] + 2. * sy[7]) * G({ b, 0, 0, a }, { sb, 1, 1, sa }, x);
        }
        return res;
    }

    if (abs(a) <= abs(b)) { // a is smallest
        const vector<complex<double>> sy
            = { G({ 0, 0, b / a }, 1), G({ 0, b, 0, c }, { 1, sb, 1, sc }, x), Log(a, sa),
                  G({ 0, b / a, 0, c / a }, 1) };
        complex<double> res { -(sy[1] * G({ 0, x / a }, 1))
            - 3. * G({ 0, 0, 0, b / a, 0, c / a }, 1)
            - 4. * G({ 0, 0, b / a, 0, 0, c / a }, 1)
            - 2. * G({ 0, 0, b / a, 0, c / a, x / a }, 1)
            - 3. * G({ 0, b / a, 0, 0, 0, c / a }, 1)
            - 2. * G({ 0, b / a, 0, 0, c / a, x / a }, 1)
            - G({ 0, b / a, 0, c / a, 0, x / a }, 1)
            + (-sy[3] - 2. * G({ 0, 0, b / a, x / a }, 1) - G({ 0, b / a, 0, x / a }, 1))
                * G({ 0, c }, { 1, sc }, x)
            - 2. * sy[0] * G({ 0, 0, c }, { 1, 1, sc }, x)
            + G({ 0, b / a }, 1) * (-sy[1] + G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x))
            + G({ 0, 0, 0, b, 0, c }, { 1, 1, 1, sb, 1, sc }, x)
            - sy[2] * sy[3] * Log(-x, sa) + (sy[3] * pow(sy[2], 2.)) / 2.
            + sy[3] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) }; // aBC
        if (b != x) {
            res += (2. * sy[0] - 2. * G({ 0, 0, x / a }, 1))
                * G({ b, 0, c }, { sb, 1, sc }, x);
        }
        if (c != x) {
            res += (2. * G({ 0, 0, b / a, 0, c / a }, 1)
                       - 2. * G({ 0, 0, b / a, 0, x / a }, 1)
                       + 2. * G({ 0, b / a, 0, 0, c / a }, 1)
                       - 2. * G({ 0, b / a, 0, 0, x / a }, 1))
                * G({ c }, { sc }, x);
        }
        return res;
    }

    // b is smallest

    const vector<complex<double>> sy = { G({ 0, a, c }, { 1, sa, sc }, x),
        G({ 0, c / b }, 1), G({ 0, a, 0, c }, { 1, sa, 1, sc }, x),
        G({ 0, c }, { 1, sc }, x), G({ 0, a / b, 0, c / b }, 1), G({ 0, 0 }, { 1, 1 }, x),
        G({ 0, 0, a / b, c / b }, 1), G({ 0, 0, c / b, a / b }, 1),
        G({ 0, c / b, 0, a / b }, 1), Log(b, sb), G({ 0, a }, { 1, sa }, x),
        G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, 0, 0, a / b, c / b }, 1),
        G({ 0, 0, 0, c / b, a / b }, 1), G({ 0, 0, a / b, 0, c / b }, 1),
        G({ 0, 0, c / b, 0, a / b }, 1), G({ 0, c / b, 0, 0, a / b }, 1) };
    complex<double> res { sy[9]
            * (-6. * sy[12] - 6. * sy[13] - 2. * sy[14] - 4. * sy[15] - 2. * sy[16]
                + 2. * sy[11] * sy[1])
        + sy[5] * (-sy[4] - 2. * sy[6] - 2. * sy[7] - sy[8])
        + 2. * sy[0] * G({ 0, 0, c / b }, 1)
        + sy[3] * (-sy[4] - G({ 0, a / b, 0, x / b }, 1))
        + 12. * G({ 0, 0, 0, 0, a / b, c / b }, 1)
        + 12. * G({ 0, 0, 0, 0, c / b, a / b }, 1)
        + 3. * G({ 0, 0, 0, a / b, 0, c / b }, 1)
        + 9. * G({ 0, 0, 0, c / b, 0, a / b }, 1)
        - 2. * G({ 0, 0, a / b, 0, 0, c / b }, 1)
        - 2. * G({ 0, 0, a / b, 0, c / b, x / b }, 1)
        - 2. * G({ 0, 0, a / b, c / b, 0, x / b }, 1)
        + 6. * G({ 0, 0, c / b, 0, 0, a / b }, 1)
        - 2. * G({ 0, 0, c / b, a / b, 0, x / b }, 1)
        - 3. * G({ 0, a / b, 0, 0, 0, c / b }, 1)
        - 2. * G({ 0, a / b, 0, 0, c / b, x / b }, 1)
        - G({ 0, a / b, 0, c / b, 0, x / b }, 1) + 3. * G({ 0, c / b, 0, 0, 0, a / b }, 1)
        - G({ 0, c / b, 0, a / b, 0, x / b }, 1)
        + G({ 0, 0, a / b }, 1) * (-2. * sy[0] + 2. * G({ 0, 0, c }, { 1, 1, sc }, x))
        + sy[1] * (-sy[2] + 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
        + G({ 0, a / b }, 1) * (-sy[2] + G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x))
        + G({ 0, a, 0, 0, 0, c }, { 1, sa, 1, 1, 1, sc }, x)
        + (6. * sy[12] + 6. * sy[13] + 2. * sy[14] + 4. * sy[15] + 2. * sy[16]
              - 2. * sy[11] * sy[1]
              + sy[9] * (-(sy[10] * sy[1]) - sy[4] - 2. * sy[6] - 2. * sy[7] - sy[8]))
            * Log(-x, sb)
        + ((sy[10] * sy[1]) / 2. + sy[4] / 2. + sy[6] + sy[7] + sy[8] / 2.)
            * pow(sy[9], 2.)
        + sy[10]
            * (-2. * sy[7] - sy[8] - 3. * G({ 0, 0, 0, c / b }, 1)
                + sy[1] * (-sy[5] - 2. * Zeta(2)))
        + sy[6] * (-2. * sy[3] - 4. * Zeta(2)) - 2. * sy[4] * Zeta(2)
        - 4. * sy[7] * Zeta(2) - 2. * sy[8] * Zeta(2) };
    if (c != x) {
        res += (2. * sy[14] - 2. * G({ 0, 0, a / b, 0, x / b }, 1)
                   + 2. * G({ 0, a / b, 0, 0, c / b }, 1)
                   - 2. * G({ 0, a / b, 0, 0, x / b }, 1))
            * G({ c }, { sc }, x);
    }
    return res;
}
