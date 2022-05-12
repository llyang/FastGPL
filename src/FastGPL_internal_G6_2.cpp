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

complex<double> G6_00abcd_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { Log(d, sd), G({ c / d, b / d, a / d }, 1),
        G({ 0, 0, a }, { 1, 1, sa }, x), G({ c / d, b / d }, 1),
        G({ 0, 0, a, b }, { 1, 1, sa, sb }, x), G({ 0, c / d, b / d, a / d }, 1),
        G({ c / d, 0, b / d, a / d }, 1), G({ c / d, b / d, 0, a / d }, 1),
        G({ 0, 0 }, { 1, 1 }, x), G({ 0, 0, a, b, c }, { 1, 1, sa, sb, sc }, x),
        G({ c / d }, 1), G({ 0, 0, c / d, b / d, a / d }, 1),
        G({ 0, c / d, 0, b / d, a / d }, 1), G({ 0, c / d, b / d, 0, a / d }, 1),
        G({ c / d, 0, 0, b / d, a / d }, 1), G({ c / d, 0, b / d, 0, a / d }, 1),
        G({ c / d, b / d, 0, 0, a / d }, 1), G({ 0 }, { 1 }, x) };
    complex<double> res { (sy[5] + sy[6] + sy[7]) * sy[8] - sy[4] * G({ 0, c / d }, 1)
        + sy[2] * (sy[1] + G({ 0, c / d, b / d }, 1) + G({ c / d, 0, b / d }, 1))
        - G({ 0, 0, 0, c / d, b / d, a / d }, 1) - G({ 0, 0, c / d, 0, b / d, a / d }, 1)
        - G({ 0, 0, c / d, b / d, 0, a / d }, 1) - G({ 0, c / d, 0, 0, b / d, a / d }, 1)
        - G({ 0, c / d, 0, b / d, 0, a / d }, 1) - G({ 0, c / d, b / d, 0, 0, a / d }, 1)
        - G({ c / d, 0, 0, 0, b / d, a / d }, 1) - G({ c / d, 0, 0, b / d, 0, a / d }, 1)
        - G({ c / d, 0, b / d, 0, 0, a / d }, 1) - G({ c / d, b / d, 0, 0, 0, a / d }, 1)
        - G({ c / d, b / d, a / d, 0, 0, x / d }, 1)
        + sy[9] * (-G({ x / d }, 1) + G({ d }, { sd }, x))
        + sy[3] * (-sy[4] - 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
        + sy[10]
            * (sy[9] + 3. * G({ 0, 0, 0, a, b }, { 1, 1, 1, sa, sb }, x)
                + G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x))
        - 3. * G({ 0, 0, 0, a, b, c }, { 1, 1, 1, sa, sb, sc }, x)
        - G({ 0, 0, a, 0, b, c }, { 1, 1, sa, 1, sb, sc }, x)
        - G({ 0, 0, a, b, 0, c }, { 1, 1, sa, sb, 1, sc }, x)
        + (-sy[5] / 2. - sy[6] / 2. - sy[7] / 2.) * pow(sy[0], 2.)
        + Log(-x, sd)
            * (-sy[11] - sy[12] - sy[13] - sy[14] - sy[15] - sy[16] + sy[2] * sy[3]
                - sy[10] * sy[4] + sy[0] * (sy[5] + sy[6] + sy[7])
                - (sy[1] * pow(sy[0], 2.)) / 2.)
        + (sy[1] * pow(sy[0], 3.)) / 6.
        + sy[0]
            * (sy[11] + sy[12] + sy[13] + sy[14] + sy[15] + sy[16] - sy[2] * sy[3]
                + sy[10] * sy[4] + sy[1] * (-sy[8] - 2. * Zeta(2)))
        + 2. * sy[5] * Zeta(2) + 2. * sy[6] * Zeta(2) + 2. * sy[7] * Zeta(2)
        + sy[1]
            * (sy[17] * sy[8] - G({ 0, 0, 0 }, { 1, 1, 1 }, x) + 2. * sy[17] * Zeta(2)) };
    return res;
}
complex<double> G6_00abcd_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;

    if (c == d) { // abcc
        const vector<complex<double>> sy
            = { G({ b / c, a / c, 1 }, 1), G({ 0, 0, a, c }, { 1, 1, sa, sc }, x),
                  G({ 0, 0, a, b }, { 1, 1, sa, sb }, x) };
        complex<double> res { (sy[1] - sy[2]) * G({ b / c, 1 }, 1)
            + G({ b / c, a / c, 0, 0, 0, 1 }, 1) - G({ b / c, a / c, 0, 0, x / c, 1 }, 1)
            + (-G({ b / c, a / c, 0, 0, 1 }, 1) + G({ b / c, a / c, 0, 0, x / c }, 1))
                * G({ c }, { sc }, x)
            + G({ b / c, a / c, 0, 1 }, 1) * G({ 0, c }, { 1, sc }, x)
            + (sy[0] - G({ b / c, 0, 1 }, 1)) * G({ 0, 0, a }, { 1, 1, sa }, x)
            - sy[0] * G({ 0, 0, c }, { 1, 1, sc }, x)
            + G({ b / c, a / c }, 1) * (-sy[1] + G({ 0, 0, 0, c }, { 1, 1, 1, sc }, x))
            + G({ b / c }, 1)
                * (-G({ 0, 0, a, 0, c }, { 1, 1, sa, 1, sc }, x)
                    + G({ 0, 0, a, b, c }, { 1, 1, sa, sb, sc }, x))
            + G({ 0, 0, a, b, 0, c }, { 1, 1, sa, sb, 1, sc }, x) - sy[2] * Zeta(2) };
        return res;
    } else { // abcd
        const vector<complex<double>> sy = { G({ 0, 0, a }, { 1, 1, sa }, x),
            G({ b / c, d / c, a / c }, 1), G({ d / c, b / c, a / c }, 1), Log(c, sc),
            G({ b / c, a / c, d / c }, 1), G({ 0, 0, 0 }, { 1, 1, 1 }, x),
            G({ 0, 0, a, b }, { 1, 1, sa, sb }, x), G({ d / c, b / c }, 1),
            G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x), G({ 0, 0, a, d }, { 1, 1, sa, sd }, x),
            G({ b / c, d / c }, 1), G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x),
            G({ 0, b / c, a / c, d / c }, 1), G({ 0, b / c, d / c, a / c }, 1),
            G({ 0, d / c, b / c, a / c }, 1), G({ b / c, 0, a / c, d / c }, 1),
            G({ b / c, 0, d / c, a / c }, 1), G({ b / c, d / c, 0, a / c }, 1),
            G({ d / c, 0, b / c, a / c }, 1), G({ d / c, b / c, 0, a / c }, 1),
            G({ d / c }, 1), G({ 0, 0, a, b, d }, { 1, 1, sa, sb, sd }, x),
            G({ 0, 0, b / c, a / c, d / c }, 1), G({ 0, 0, b / c, d / c, a / c }, 1),
            G({ 0, 0, d / c, b / c, a / c }, 1), G({ 0, b / c, 0, a / c, d / c }, 1),
            G({ 0, b / c, 0, d / c, a / c }, 1), G({ 0, b / c, d / c, 0, a / c }, 1),
            G({ 0, d / c, 0, b / c, a / c }, 1), G({ 0, d / c, b / c, 0, a / c }, 1),
            G({ b / c, 0, 0, a / c, d / c }, 1), G({ b / c, 0, 0, d / c, a / c }, 1),
            G({ b / c, 0, d / c, 0, a / c }, 1), G({ b / c, d / c, 0, 0, a / c }, 1),
            G({ d / c, 0, 0, b / c, a / c }, 1), G({ d / c, 0, b / c, 0, a / c }, 1),
            G({ d / c, b / c, 0, 0, a / c }, 1) };
        complex<double> res { sy[11]
                * (-sy[13] - sy[14] - sy[15] - sy[16] - sy[17] - sy[18] - sy[19] - sy[20]
                    + sy[12] * (-sy[1] - sy[2]))
            + (sy[1] + sy[2]) * sy[5] + sy[10] * (sy[9] + 3. * sy[8])
            + sy[7] * (sy[6] + 3. * sy[8]) + sy[6] * G({ 0, d / c }, 1)
            + sy[0]
                * (-sy[1] - sy[2] - G({ 0, b / c, d / c }, 1) - G({ 0, d / c, b / c }, 1)
                    - G({ b / c, 0, d / c }, 1) - G({ d / c, 0, b / c }, 1))
            + G({ 0, 0, 0, b / c, a / c, d / c }, 1)
            + G({ 0, 0, 0, b / c, d / c, a / c }, 1)
            + G({ 0, 0, 0, d / c, b / c, a / c }, 1)
            + G({ 0, 0, b / c, 0, a / c, d / c }, 1)
            + G({ 0, 0, b / c, 0, d / c, a / c }, 1)
            + G({ 0, 0, b / c, d / c, 0, a / c }, 1)
            + G({ 0, 0, d / c, 0, b / c, a / c }, 1)
            + G({ 0, 0, d / c, b / c, 0, a / c }, 1)
            + G({ 0, b / c, 0, 0, a / c, d / c }, 1)
            + G({ 0, b / c, 0, 0, d / c, a / c }, 1)
            + G({ 0, b / c, 0, d / c, 0, a / c }, 1)
            + G({ 0, b / c, d / c, 0, 0, a / c }, 1)
            + G({ 0, d / c, 0, 0, b / c, a / c }, 1)
            + G({ 0, d / c, 0, b / c, 0, a / c }, 1)
            + G({ 0, d / c, b / c, 0, 0, a / c }, 1)
            + G({ b / c, 0, 0, 0, a / c, d / c }, 1)
            + G({ b / c, 0, 0, 0, d / c, a / c }, 1)
            + G({ b / c, 0, 0, d / c, 0, a / c }, 1)
            + G({ b / c, 0, d / c, 0, 0, a / c }, 1)
            + G({ b / c, a / c, 0, 0, 0, d / c }, 1)
            + G({ b / c, a / c, 0, 0, d / c, x / c }, 1)
            + G({ b / c, a / c, 0, d / c, 0, x / c }, 1)
            + G({ b / c, a / c, d / c, 0, 0, x / c }, 1)
            + G({ b / c, d / c, 0, 0, 0, a / c }, 1)
            + G({ b / c, d / c, a / c, 0, 0, x / c }, 1)
            + G({ d / c, 0, 0, 0, b / c, a / c }, 1)
            + G({ d / c, 0, 0, b / c, 0, a / c }, 1)
            + G({ d / c, 0, b / c, 0, 0, a / c }, 1)
            + G({ d / c, b / c, 0, 0, 0, a / c }, 1)
            + G({ d / c, b / c, a / c, 0, 0, x / c }, 1)
            + G({ b / c, a / c, 0, d / c }, 1) * G({ 0, d }, { 1, sd }, x)
            + G({ b / c, a / c }, 1) * (-sy[9] + G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x))
            + sy[21]
                * (-sy[22] - 3. * G({ 0, 0, 0, a, b }, { 1, 1, 1, sa, sb }, x)
                    - G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x))
            + G({ b / c }, 1) * (sy[22] - G({ 0, 0, a, 0, d }, { 1, 1, sa, 1, sd }, x))
            + G({ 0, 0, a, b, 0, d }, { 1, 1, sa, sb, 1, sd }, x)
            + (sy[13] / 2. + sy[14] / 2. + sy[15] / 2. + sy[16] / 2. + sy[17] / 2.
                  + sy[18] / 2. + sy[19] / 2. + sy[20] / 2.)
                * pow(sy[3], 2.)
            + Log(-x, sc)
                * (-(sy[0] * sy[10]) + sy[23] + sy[24] + sy[25] + sy[26] + sy[27] + sy[28]
                    + sy[29] + sy[30] + sy[31] + sy[32] + sy[33] + sy[34] + sy[35]
                    + sy[36] + sy[37]
                    + (-sy[13] - sy[14] - sy[15] - sy[16] - sy[17] - sy[18] - sy[19]
                          - sy[20])
                        * sy[3]
                    + sy[21] * sy[6] - sy[0] * sy[7]
                    + (sy[1] / 2. + sy[2] / 2. + sy[4] / 2.) * pow(sy[3], 2.))
            + (-sy[1] / 6. - sy[2] / 6. - sy[4] / 6.) * pow(sy[3], 3.)
            - 2. * sy[13] * Zeta(2) - 2. * sy[14] * Zeta(2) - 2. * sy[15] * Zeta(2)
            - 2. * sy[16] * Zeta(2) - 2. * sy[17] * Zeta(2) - 2. * sy[18] * Zeta(2)
            - 2. * sy[19] * Zeta(2) - 2. * sy[20] * Zeta(2)
            + sy[4]
                * (-(sy[11] * sy[12]) + sy[5] - G({ 0, 0, d }, { 1, 1, sd }, x)
                    - 2. * sy[12] * Zeta(2))
            + sy[12] * (-2. * sy[1] * Zeta(2) - 2. * sy[2] * Zeta(2))
            + sy[3]
                * (sy[0] * sy[10] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27] - sy[28]
                    + sy[11] * (sy[1] + sy[2]) - sy[29] - sy[30] - sy[31] - sy[32]
                    - sy[33] - sy[34] - sy[35] - sy[36] - sy[37] - sy[21] * sy[6]
                    + sy[0] * sy[7] + 2. * sy[1] * Zeta(2) + 2. * sy[2] * Zeta(2)
                    + sy[4] * (sy[11] + 2. * Zeta(2))) };
        if (d != x) {
            res += (-G({ b / c, a / c, 0, 0, d / c }, 1)
                       + G({ b / c, a / c, 0, 0, x / c }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_00abcd_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (b == c) {
        if (b == d) { // abbb
            const vector<complex<double>> sy
                = { G({ a / b, 0, 1 }, 1), G({ a / b, 1, 1 }, 1),
                      G({ b, b }, { sb, sb }, x), G({ a / b, 0, 0, 1 }, 1) };
            complex<double> res { sy[2] * sy[3] - sy[2] * G({ a / b, 0, 0, x / b }, 1)
                + G({ a / b, 0, 0, 0, 1, 1 }, 1) - G({ a / b, 0, 0, x / b, 1, 1 }, 1)
                + (-G({ a / b, 0, 0, 1, 1 }, 1) + G({ a / b, 0, 0, x / b, 1 }, 1))
                    * G({ b }, { sb }, x)
                + (-sy[3] + G({ a / b, 0, 1, 1 }, 1)) * G({ 0, b }, { 1, sb }, x)
                + (sy[0] - sy[1]) * G({ 0, 0, b }, { 1, 1, sb }, x)
                - sy[0] * G({ 0, b, b }, { 1, sb, sb }, x)
                + G({ a / b, 1 }, 1)
                    * (-G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                        + G({ 0, 0, b, b }, { 1, 1, sb, sb }, x))
                + G({ a / b }, 1)
                    * (-G({ 0, 0, 0, b, b }, { 1, 1, 1, sb, sb }, x)
                        + G({ 0, 0, a, b, b }, { 1, 1, sa, sb, sb }, x))
                + G({ 0, 0, a, 0, b, b }, { 1, 1, sa, 1, sb, sb }, x)
                + G({ 0, 0, a }, { 1, 1, sa }, x) * (sy[1] - Zeta(3)) };
            return res;
        } else { // abbd
            const vector<complex<double>> sy = { G({ a / b, 0, 1 }, 1),
                G({ a / b, d / b, 1 }, 1), G({ 0, 0, a }, { 1, 1, sa }, x),
                G({ d / b, 1, a / b }, 1), G({ d / b, a / b, 1 }, 1), Log(b, sb),
                G({ 0, 0, 0 }, { 1, 1, 1 }, x), G({ d / b, 1 }, 1),
                G({ 0, 0, a, d }, { 1, 1, sa, sd }, x), G({ b, d }, { sb, sd }, x),
                G({ a / b, 0, 0, 1 }, 1), G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x),
                G({ 0, a / b, d / b, 1 }, 1), G({ 0, d / b, 1, a / b }, 1),
                G({ 0, d / b, a / b, 1 }, 1), G({ d / b, 0, 1, a / b }, 1),
                G({ d / b, 0, a / b, 1 }, 1), G({ d / b, 1, 0, a / b }, 1),
                G({ 0, 0, a / b, d / b, 1 }, 1), G({ 0, 0, d / b, 1, a / b }, 1),
                G({ 0, 0, d / b, a / b, 1 }, 1), G({ 0, d / b, 0, 1, a / b }, 1),
                G({ 0, d / b, 0, a / b, 1 }, 1), G({ 0, d / b, 1, 0, a / b }, 1),
                G({ d / b, 0, 0, 1, a / b }, 1), G({ d / b, 0, 0, a / b, 1 }, 1),
                G({ d / b, 0, 1, 0, a / b }, 1), G({ d / b, 1, 0, 0, a / b }, 1) };
            complex<double> res { sy[9] * sy[10]
                + sy[11]
                    * (-sy[13] - sy[14] - sy[15] - sy[16] - sy[17] - sy[18]
                        + sy[12] * (-sy[3] - sy[4]))
                + (sy[3] + sy[4]) * sy[6]
                + sy[2] * (-sy[3] - sy[4] - G({ 0, d / b, 1 }, 1))
                - sy[9] * G({ a / b, 0, 0, x / b }, 1)
                + G({ 0, 0, 0, a / b, d / b, 1 }, 1) + G({ 0, 0, 0, d / b, 1, a / b }, 1)
                + G({ 0, 0, 0, d / b, a / b, 1 }, 1) + G({ 0, 0, d / b, 0, 1, a / b }, 1)
                + G({ 0, 0, d / b, 0, a / b, 1 }, 1) + G({ 0, 0, d / b, 1, 0, a / b }, 1)
                + G({ 0, d / b, 0, 0, 1, a / b }, 1) + G({ 0, d / b, 0, 0, a / b, 1 }, 1)
                + G({ 0, d / b, 0, 1, 0, a / b }, 1) + G({ 0, d / b, 1, 0, 0, a / b }, 1)
                + G({ a / b, 0, 0, 0, d / b, 1 }, 1)
                + G({ a / b, 0, 0, d / b, 1, x / b }, 1)
                + G({ a / b, 0, 0, d / b, x / b, 1 }, 1)
                + G({ a / b, 0, d / b, 0, 1, x / b }, 1)
                + G({ a / b, 0, d / b, 0, x / b, 1 }, 1)
                + G({ a / b, 0, d / b, 1, 0, x / b }, 1)
                + G({ a / b, d / b, 0, 0, 1, x / b }, 1)
                + G({ a / b, d / b, 0, 0, x / b, 1 }, 1)
                + G({ a / b, d / b, 0, 1, 0, x / b }, 1)
                + G({ a / b, d / b, 1, 0, 0, x / b }, 1)
                + G({ d / b, 0, 0, 0, 1, a / b }, 1) + G({ d / b, 0, 0, 0, a / b, 1 }, 1)
                + G({ d / b, 0, 0, 1, 0, a / b }, 1) + G({ d / b, 0, 1, 0, 0, a / b }, 1)
                + G({ d / b, 1, 0, 0, 0, a / b }, 1)
                + G({ d / b, 1, a / b, 0, 0, x / b }, 1)
                + G({ d / b, a / b, 0, 0, 1, x / b }, 1)
                + G({ d / b, a / b, 0, 0, x / b, 1 }, 1)
                + G({ d / b, a / b, 0, 1, 0, x / b }, 1)
                + G({ d / b, a / b, 1, 0, 0, x / b }, 1)
                + (-sy[10] + G({ a / b, 0, d / b, 1 }, 1)) * G({ 0, d }, { 1, sd }, x)
                + (sy[0] - sy[1]) * G({ 0, 0, d }, { 1, 1, sd }, x)
                - sy[0] * G({ 0, b, d }, { 1, sb, sd }, x)
                + sy[7] * (sy[8] + 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
                + G({ a / b, 1 }, 1) * (-sy[8] + G({ 0, 0, b, d }, { 1, 1, sb, sd }, x))
                + G({ a / b }, 1)
                    * (-G({ 0, 0, 0, b, d }, { 1, 1, 1, sb, sd }, x)
                        + G({ 0, 0, a, b, d }, { 1, 1, sa, sb, sd }, x))
                + G({ 0, 0, a, 0, b, d }, { 1, 1, sa, 1, sb, sd }, x)
                + (sy[13] / 2. + sy[14] / 2. + sy[15] / 2. + sy[16] / 2. + sy[17] / 2.
                      + sy[18] / 2.)
                    * pow(sy[5], 2.)
                + Log(-x, sb)
                    * (sy[19] + sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25]
                        + sy[26] + sy[27] + sy[28]
                        + (-sy[13] - sy[14] - sy[15] - sy[16] - sy[17] - sy[18]) * sy[5]
                        - sy[2] * sy[7]
                        + (sy[1] / 2. + sy[3] / 2. + sy[4] / 2.) * pow(sy[5], 2.))
                + (-sy[1] / 6. - sy[3] / 6. - sy[4] / 6.) * pow(sy[5], 3.)
                - 2. * sy[13] * Zeta(2) - 2. * sy[14] * Zeta(2) - 2. * sy[15] * Zeta(2)
                - 2. * sy[16] * Zeta(2) - 2. * sy[17] * Zeta(2) - 2. * sy[18] * Zeta(2)
                + sy[1] * (-(sy[11] * sy[12]) + sy[6] - 2. * sy[12] * Zeta(2))
                + sy[12] * (-2. * sy[3] * Zeta(2) - 2. * sy[4] * Zeta(2))
                + sy[5]
                    * (-sy[19] - sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25]
                        - sy[26] - sy[27] - sy[28] + sy[11] * (sy[3] + sy[4])
                        + sy[2] * sy[7] + 2. * sy[3] * Zeta(2) + 2. * sy[4] * Zeta(2)
                        + sy[1] * (sy[11] + 2. * Zeta(2))) };
            if (d != x) {
                res += (-G({ a / b, 0, 0, d / b, 1 }, 1)
                           + G({ a / b, 0, 0, x / b, 1 }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        }
    } else { // abcd
        const vector<complex<double>> sy = { G({ a / b, 0, c / b }, 1),
            G({ a / b, c / b, d / b }, 1), G({ c / b, a / b, d / b }, 1),
            G({ 0, 0, 0 }, { 1, 1, 1 }, x), G({ c / b, d / b, a / b }, 1), Log(b, sb),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ c / b, a / b }, 1),
            G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x), G({ c / b, d / b }, 1),
            G({ 0, 0, a, d }, { 1, 1, sa, sd }, x), G({ a / b, 0, 0, c / b }, 1),
            G({ 0, a / b, c / b, d / b }, 1), G({ 0, c / b, a / b, d / b }, 1),
            G({ 0, c / b, d / b, a / b }, 1), G({ c / b, 0, a / b, d / b }, 1),
            G({ c / b, 0, d / b, a / b }, 1), G({ c / b, d / b, 0, a / b }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x),
            G({ 0, 0, a, c, d }, { 1, 1, sa, sc, sd }, x),
            G({ 0, 0, a / b, c / b, d / b }, 1), G({ 0, 0, c / b, a / b, d / b }, 1),
            G({ 0, 0, c / b, d / b, a / b }, 1), G({ 0, c / b, 0, a / b, d / b }, 1),
            G({ 0, c / b, 0, d / b, a / b }, 1), G({ 0, c / b, d / b, 0, a / b }, 1),
            G({ c / b, 0, 0, a / b, d / b }, 1), G({ c / b, 0, 0, d / b, a / b }, 1),
            G({ c / b, 0, d / b, 0, a / b }, 1), G({ c / b, d / b, 0, 0, a / b }, 1) };
        complex<double> res { sy[3] * (-sy[2] - sy[4])
            + sy[18]
                * (sy[12] + sy[13] + sy[14] + sy[15] + sy[16] + sy[17]
                    + sy[19] * (sy[2] + sy[4]))
            + sy[10] * sy[7] - sy[7] * sy[8]
            + sy[6] * (sy[4] + G({ 0, c / b, d / b }, 1) + G({ c / b, 0, d / b }, 1))
            - G({ 0, 0, 0, a / b, c / b, d / b }, 1)
            - G({ 0, 0, 0, c / b, a / b, d / b }, 1)
            - G({ 0, 0, 0, c / b, d / b, a / b }, 1)
            - G({ 0, 0, c / b, 0, a / b, d / b }, 1)
            - G({ 0, 0, c / b, 0, d / b, a / b }, 1)
            - G({ 0, 0, c / b, d / b, 0, a / b }, 1)
            - G({ 0, c / b, 0, 0, a / b, d / b }, 1)
            - G({ 0, c / b, 0, 0, d / b, a / b }, 1)
            - G({ 0, c / b, 0, d / b, 0, a / b }, 1)
            - G({ 0, c / b, d / b, 0, 0, a / b }, 1)
            - G({ a / b, 0, 0, 0, c / b, d / b }, 1)
            - G({ a / b, 0, 0, c / b, 0, d / b }, 1)
            - G({ a / b, 0, 0, c / b, d / b, x / b }, 1)
            - G({ a / b, 0, c / b, 0, 0, d / b }, 1)
            - G({ a / b, 0, c / b, 0, d / b, x / b }, 1)
            - G({ a / b, 0, c / b, d / b, 0, x / b }, 1)
            - G({ a / b, c / b, 0, 0, 0, d / b }, 1)
            - G({ a / b, c / b, 0, 0, d / b, x / b }, 1)
            - G({ a / b, c / b, 0, d / b, 0, x / b }, 1)
            - G({ a / b, c / b, d / b, 0, 0, x / b }, 1)
            - G({ c / b, 0, 0, 0, a / b, d / b }, 1)
            - G({ c / b, 0, 0, 0, d / b, a / b }, 1)
            - G({ c / b, 0, 0, d / b, 0, a / b }, 1)
            - G({ c / b, 0, d / b, 0, 0, a / b }, 1)
            - G({ c / b, a / b, 0, 0, 0, d / b }, 1)
            - G({ c / b, a / b, 0, 0, d / b, x / b }, 1)
            - G({ c / b, a / b, 0, d / b, 0, x / b }, 1)
            - G({ c / b, a / b, d / b, 0, 0, x / b }, 1)
            - G({ c / b, d / b, 0, 0, 0, a / b }, 1)
            - G({ c / b, d / b, a / b, 0, 0, x / b }, 1)
            + (-sy[11] - G({ a / b, 0, c / b, d / b }, 1)
                  - G({ a / b, c / b, 0, d / b }, 1) - G({ c / b, a / b, 0, d / b }, 1))
                * G({ 0, d }, { 1, sd }, x)
            + (sy[0] + sy[1] + sy[2]) * G({ 0, 0, d }, { 1, 1, sd }, x)
            - sy[0] * G({ 0, c, d }, { 1, sc, sd }, x)
            + sy[9] * (-sy[10] - 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
            + G({ a / b, c / b }, 1) * (-sy[8] + G({ 0, 0, c, d }, { 1, 1, sc, sd }, x))
            + G({ a / b }, 1) * (sy[20] - G({ 0, 0, 0, c, d }, { 1, 1, 1, sc, sd }, x))
            + G({ c / b }, 1) * (-sy[20] + G({ 0, 0, a, 0, d }, { 1, 1, sa, 1, sd }, x))
            + G({ 0, 0, a, 0, c, d }, { 1, 1, sa, 1, sc, sd }, x)
            + (-sy[12] / 2. - sy[13] / 2. - sy[14] / 2. - sy[15] / 2. - sy[16] / 2.
                  - sy[17] / 2.)
                * pow(sy[5], 2.)
            + Log(-x, sb)
                * (-sy[21] - sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27] - sy[28]
                    - sy[29] - sy[30]
                    + (sy[12] + sy[13] + sy[14] + sy[15] + sy[16] + sy[17]) * sy[5]
                    + sy[9] * sy[6]
                    + (-sy[1] / 2. - sy[2] / 2. - sy[4] / 2.) * pow(sy[5], 2.))
            + (sy[1] / 6. + sy[2] / 6. + sy[4] / 6.) * pow(sy[5], 3.)
            + 2. * sy[12] * Zeta(2) + 2. * sy[13] * Zeta(2) + 2. * sy[14] * Zeta(2)
            + 2. * sy[15] * Zeta(2) + 2. * sy[16] * Zeta(2) + 2. * sy[17] * Zeta(2)
            + sy[1] * (sy[18] * sy[19] - sy[3] + 2. * sy[19] * Zeta(2))
            + sy[5]
                * (sy[21] + sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27] + sy[28]
                    + sy[29] + sy[30] + sy[18] * (-sy[2] - sy[4]) - sy[9] * sy[6]
                    + sy[1] * (-sy[18] - 2. * Zeta(2)) - 2. * sy[2] * Zeta(2)
                    - 2. * sy[4] * Zeta(2))
            + sy[19] * (2. * sy[2] * Zeta(2) + 2. * sy[4] * Zeta(2)) };
        if (c != x) {
            res += (sy[11] - G({ a / b, 0, 0, x / b }, 1)) * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (G({ a / b, 0, 0, c / b, d / b }, 1)
                       - G({ a / b, 0, 0, c / b, x / b }, 1)
                       + G({ a / b, 0, c / b, 0, d / b }, 1)
                       - G({ a / b, 0, c / b, 0, x / b }, 1)
                       + G({ a / b, c / b, 0, 0, d / b }, 1)
                       - G({ a / b, c / b, 0, 0, x / b }, 1)
                       + G({ c / b, a / b, 0, 0, d / b }, 1)
                       - G({ c / b, a / b, 0, 0, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_00abcd_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (a == b) {
        if (a == c) { // aaad
            const vector<complex<double>> sy = { G({ a, a, d }, { sa, sa, sd }, x),
                Log(a, sa), G({ d / a, 1, 1 }, 1), G({ 0, 0 }, { 1, 1 }, x),
                G({ 0 }, { 1 }, x), G({ 0, 0, d }, { 1, 1, sd }, x) };
            complex<double> res { sy[0] * G({ 0, 0, x / a }, 1)
                + G({ 0, 0, 0, d / a, 1, 1 }, 1) + G({ 0, 0, d / a, 1, 1, x / a }, 1)
                + G({ 0, 0, d / a, 1, x / a, 1 }, 1) + G({ 0, 0, d / a, x / a, 1, 1 }, 1)
                + G({ 0, d / a, 0, 1, 1, x / a }, 1) + G({ 0, d / a, 0, 1, x / a, 1 }, 1)
                + G({ 0, d / a, 0, x / a, 1, 1 }, 1) + G({ 0, d / a, 1, 0, 1, x / a }, 1)
                + G({ 0, d / a, 1, 0, x / a, 1 }, 1) + G({ 0, d / a, 1, 1, 0, x / a }, 1)
                + G({ d / a, 0, 0, 1, 1, x / a }, 1) + G({ d / a, 0, 0, 1, x / a, 1 }, 1)
                + G({ d / a, 0, 0, x / a, 1, 1 }, 1) + G({ d / a, 0, 1, 0, 1, x / a }, 1)
                + G({ d / a, 0, 1, 0, x / a, 1 }, 1) + G({ d / a, 0, 1, 1, 0, x / a }, 1)
                + G({ d / a, 1, 0, 0, 1, x / a }, 1) + G({ d / a, 1, 0, 0, x / a, 1 }, 1)
                + G({ d / a, 1, 0, 1, 0, x / a }, 1) + G({ d / a, 1, 1, 0, 0, x / a }, 1)
                + G({ 0, 0, 0, a, a, d }, { 1, 1, 1, sa, sa, sd }, x)
                + (sy[2] * Log(-x, sa) * pow(sy[1], 2.)) / 2.
                - (sy[2] * pow(sy[1], 3.)) / 6.
                + G({ 0, 0, a, d }, { 1, 1, sa, sd }, x) * Zeta(2)
                - G({ 0, a, a, d }, { 1, sa, sa, sd }, x) * Zeta(2)
                + sy[1] * sy[2] * (sy[3] + 2. * Zeta(2))
                + sy[2]
                    * (-(sy[3] * sy[4]) - sy[5] + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                        - 2. * sy[4] * Zeta(2))
                + G({ 0, d }, { 1, sd }, x) * (G({ 0, d / a, 1, 1 }, 1) - Zeta(4) / 4.)
                + G({ a, d }, { sa, sd }, x) * (-G({ 0, 0, x / a, 1 }, 1) + Zeta(4) / 4.)
                + sy[0] * Zeta(3) + sy[5] * Zeta(3)
                - 2. * G({ 0, a, d }, { 1, sa, sd }, x) * Zeta(3) };
            if (d != x) {
                res += (-G({ 0, 0, d / a, 1, 1 }, 1) + G({ 0, 0, x / a, 1, 1 }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        } else { // aacd
            const vector<complex<double>> sy = { G({ 0, c, d }, { 1, sc, sd }, x),
                G({ 0, c / a, 1 }, 1), G({ a, c, d }, { sa, sc, sd }, x), Log(a, sa),
                G({ c / a, 1, d / a }, 1), G({ c / a, d / a, 1 }, 1), G({ 0 }, { 1 }, x),
                G({ 0, 0 }, { 1, 1 }, x), G({ 0, 0, 0 }, { 1, 1, 1 }, x),
                G({ 0, 0, c, d }, { 1, 1, sc, sd }, x), G({ 0, 0, c / a, 1 }, 1) };

            complex<double> res { -(sy[0] * sy[1]) + sy[5] * sy[6] * sy[7] - sy[5] * sy[8]
                + sy[2] * G({ 0, 0, x / a }, 1) - G({ 0, 0, 0, c / a, 1, d / a }, 1)
                - G({ 0, 0, 0, c / a, d / a, 1 }, 1) - G({ 0, 0, c / a, 0, 1, d / a }, 1)
                - G({ 0, 0, c / a, 0, d / a, 1 }, 1) - G({ 0, 0, c / a, 1, 0, d / a }, 1)
                - G({ 0, 0, c / a, 1, d / a, x / a }, 1)
                - G({ 0, 0, c / a, d / a, 1, x / a }, 1)
                - G({ 0, 0, c / a, d / a, x / a, 1 }, 1)
                - G({ 0, c / a, 0, 0, 1, d / a }, 1) - G({ 0, c / a, 0, 0, d / a, 1 }, 1)
                - G({ 0, c / a, 0, 1, 0, d / a }, 1)
                - G({ 0, c / a, 0, 1, d / a, x / a }, 1)
                - G({ 0, c / a, 0, d / a, 1, x / a }, 1)
                - G({ 0, c / a, 0, d / a, x / a, 1 }, 1)
                - G({ 0, c / a, 1, 0, 0, d / a }, 1)
                - G({ 0, c / a, 1, 0, d / a, x / a }, 1)
                - G({ 0, c / a, 1, d / a, 0, x / a }, 1)
                - G({ 0, c / a, d / a, 0, 1, x / a }, 1)
                - G({ 0, c / a, d / a, 0, x / a, 1 }, 1)
                - G({ 0, c / a, d / a, 1, 0, x / a }, 1)
                - G({ c / a, 0, 0, 0, 1, d / a }, 1) - G({ c / a, 0, 0, 0, d / a, 1 }, 1)
                - G({ c / a, 0, 0, 1, 0, d / a }, 1)
                - G({ c / a, 0, 0, 1, d / a, x / a }, 1)
                - G({ c / a, 0, 0, d / a, 1, x / a }, 1)
                - G({ c / a, 0, 0, d / a, x / a, 1 }, 1)
                - G({ c / a, 0, 1, 0, 0, d / a }, 1)
                - G({ c / a, 0, 1, 0, d / a, x / a }, 1)
                - G({ c / a, 0, 1, d / a, 0, x / a }, 1)
                - G({ c / a, 0, d / a, 0, 1, x / a }, 1)
                - G({ c / a, 0, d / a, 0, x / a, 1 }, 1)
                - G({ c / a, 0, d / a, 1, 0, x / a }, 1)
                - G({ c / a, 1, 0, 0, 0, d / a }, 1)
                - G({ c / a, 1, 0, 0, d / a, x / a }, 1)
                - G({ c / a, 1, 0, d / a, 0, x / a }, 1)
                - G({ c / a, 1, d / a, 0, 0, x / a }, 1)
                - G({ c / a, d / a, 0, 0, 1, x / a }, 1)
                - G({ c / a, d / a, 0, 0, x / a, 1 }, 1)
                - G({ c / a, d / a, 0, 1, 0, x / a }, 1)
                - G({ c / a, d / a, 1, 0, 0, x / a }, 1)
                + (-sy[10] - G({ 0, c / a, 1, d / a }, 1) - G({ 0, c / a, d / a, 1 }, 1)
                      - G({ c / a, 0, 1, d / a }, 1) - G({ c / a, 0, d / a, 1 }, 1)
                      - G({ c / a, 1, 0, d / a }, 1))
                    * G({ 0, d }, { 1, sd }, x)
                + (sy[1] + sy[4] + sy[5]) * G({ 0, 0, d }, { 1, 1, sd }, x)
                + G({ c / a, 1 }, 1) * (sy[9] - G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x))
                + G({ 0, 0, 0, a, c, d }, { 1, 1, 1, sa, sc, sd }, x)
                + (-sy[4] / 2. - sy[5] / 2.) * Log(-x, sa) * pow(sy[3], 2.)
                + (sy[4] / 6. + sy[5] / 6.) * pow(sy[3], 3.) + sy[9] * Zeta(2)
                + 2. * sy[5] * sy[6] * Zeta(2)
                - G({ 0, a, c, d }, { 1, sa, sc, sd }, x) * Zeta(2)
                + sy[3]
                    * (-(sy[5] * sy[7]) + sy[4] * (-sy[7] - 2. * Zeta(2))
                        - 2. * sy[5] * Zeta(2))
                + sy[4] * (sy[6] * sy[7] - sy[8] + 2. * sy[6] * Zeta(2)) - sy[0] * Zeta(3)
                + sy[2] * Zeta(3) };
            if (c != x) {
                res += (sy[10] - G({ 0, 0, x / a, 1 }, 1)) * G({ c, d }, { sc, sd }, x);
            }
            if (d != x) {
                res += (G({ 0, 0, c / a, 1, d / a }, 1) - G({ 0, 0, c / a, 1, x / a }, 1)
                           + G({ 0, 0, c / a, d / a, 1 }, 1)
                           - G({ 0, 0, c / a, x / a, 1 }, 1)
                           + G({ 0, c / a, 0, 1, d / a }, 1)
                           - G({ 0, c / a, 0, 1, x / a }, 1)
                           + G({ 0, c / a, 0, d / a, 1 }, 1)
                           - G({ 0, c / a, 0, x / a, 1 }, 1)
                           + G({ 0, c / a, 1, 0, d / a }, 1)
                           - G({ 0, c / a, 1, 0, x / a }, 1)
                           + G({ c / a, 0, 0, 1, d / a }, 1)
                           - G({ c / a, 0, 0, 1, x / a }, 1)
                           + G({ c / a, 0, 0, d / a, 1 }, 1)
                           - G({ c / a, 0, 0, x / a, 1 }, 1)
                           + G({ c / a, 0, 1, 0, d / a }, 1)
                           - G({ c / a, 0, 1, 0, x / a }, 1)
                           + G({ c / a, 1, 0, 0, d / a }, 1)
                           - G({ c / a, 1, 0, 0, x / a }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        }
    } else { // abcd
        const vector<complex<double>> sy = { G({ 0, 0, b / a }, 1),
            G({ 0, c, d }, { 1, sc, sd }, x), G({ 0, b / a, c / a }, 1),
            G({ b / a, 0, c / a }, 1), G({ b / a, c / a, d / a }, 1), Log(a, sa),
            G({ 0, 0, c, d }, { 1, 1, sc, sd }, x), G({ 0, 0, b / a, c / a }, 1),
            G({ 0, b / a, 0, c / a }, 1), G({ b / a, 0, 0, c / a }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x) };
        complex<double> res { sy[0] * sy[1] + sy[1] * sy[2] + sy[1] * sy[3]
            + G({ 0, 0, 0, b / a, c / a, d / a }, 1)
            + G({ 0, 0, b / a, 0, c / a, d / a }, 1)
            + G({ 0, 0, b / a, c / a, 0, d / a }, 1)
            + G({ 0, 0, b / a, c / a, d / a, x / a }, 1)
            + G({ 0, b / a, 0, 0, c / a, d / a }, 1)
            + G({ 0, b / a, 0, c / a, 0, d / a }, 1)
            + G({ 0, b / a, 0, c / a, d / a, x / a }, 1)
            + G({ 0, b / a, c / a, 0, 0, d / a }, 1)
            + G({ 0, b / a, c / a, 0, d / a, x / a }, 1)
            + G({ 0, b / a, c / a, d / a, 0, x / a }, 1)
            + G({ b / a, 0, 0, 0, c / a, d / a }, 1)
            + G({ b / a, 0, 0, c / a, 0, d / a }, 1)
            + G({ b / a, 0, 0, c / a, d / a, x / a }, 1)
            + G({ b / a, 0, c / a, 0, 0, d / a }, 1)
            + G({ b / a, 0, c / a, 0, d / a, x / a }, 1)
            + G({ b / a, 0, c / a, d / a, 0, x / a }, 1)
            + G({ b / a, c / a, 0, 0, 0, d / a }, 1)
            + G({ b / a, c / a, 0, 0, d / a, x / a }, 1)
            + G({ b / a, c / a, 0, d / a, 0, x / a }, 1)
            + G({ b / a, c / a, d / a, 0, 0, x / a }, 1)
            + (sy[9] + sy[7] + sy[8] + G({ 0, b / a, c / a, d / a }, 1)
                  + G({ b / a, 0, c / a, d / a }, 1) + G({ b / a, c / a, 0, d / a }, 1))
                * G({ 0, d }, { 1, sd }, x)
            + (-sy[2] - sy[3] - sy[4]) * G({ 0, 0, d }, { 1, 1, sd }, x)
            + G({ b / a, c / a }, 1) * (-sy[6] + G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x))
            + G({ 0, b / a }, 1) * (-sy[6] + G({ 0, b, c, d }, { 1, sb, sc, sd }, x))
            + G({ b / a }, 1)
                * (G({ 0, 0, 0, c, d }, { 1, 1, 1, sc, sd }, x)
                    - G({ 0, 0, b, c, d }, { 1, 1, sb, sc, sd }, x))
            + G({ 0, 0, 0, b, c, d }, { 1, 1, 1, sb, sc, sd }, x)
            + (sy[4] * Log(-x, sa) * pow(sy[5], 2.)) / 2. - (sy[4] * pow(sy[5], 3.)) / 6.
            + sy[4] * sy[5] * (sy[10] + 2. * Zeta(2))
            + sy[4]
                * (-(sy[10] * sy[11]) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                    - 2. * sy[11] * Zeta(2)) };
        if (b != x) {
            res += (-sy[0] + G({ 0, 0, x / a }, 1)) * G({ b, c, d }, { sb, sc, sd }, x);
        }
        if (c != x) {
            res += (-sy[9] - sy[7] - sy[8] + G({ 0, 0, b / a, x / a }, 1)
                       + G({ 0, b / a, 0, x / a }, 1) + G({ b / a, 0, 0, x / a }, 1))
                * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (-G({ 0, 0, b / a, c / a, d / a }, 1)
                       + G({ 0, 0, b / a, c / a, x / a }, 1)
                       - G({ 0, b / a, 0, c / a, d / a }, 1)
                       + G({ 0, b / a, 0, c / a, x / a }, 1)
                       - G({ 0, b / a, c / a, 0, d / a }, 1)
                       + G({ 0, b / a, c / a, 0, x / a }, 1)
                       - G({ b / a, 0, 0, c / a, d / a }, 1)
                       + G({ b / a, 0, 0, c / a, x / a }, 1)
                       - G({ b / a, 0, c / a, 0, d / a }, 1)
                       + G({ b / a, 0, c / a, 0, x / a }, 1)
                       - G({ b / a, c / a, 0, 0, d / a }, 1)
                       + G({ b / a, c / a, 0, 0, x / a }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G6_a00bcd_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { Log(d, sd), G({ a }, { sa }, x),
        G({ c / d, b / d }, 1), G({ 0, a }, { 1, sa }, x), G({ 0, c / d, b / d }, 1),
        G({ c / d, 0, b / d }, 1), G({ 0, 0, a }, { 1, 1, sa }, x),
        G({ a, 0, 0, b }, { sa, 1, 1, sb }, x), G({ 0, 0, c / d, b / d }, 1),
        G({ 0, c / d, 0, b / d }, 1), G({ c / d, 0, 0, b / d }, 1),
        G({ a, 0, 0, b, c }, { sa, 1, 1, sb, sc }, x), G({ c / d }, 1),
        G({ c / d, b / d, 0, 0, a / d }, 1), G({ 0, 0 }, { 1, 1 }, x),
        G({ 0 }, { 1 }, x) };
    complex<double> res { (sy[4] + sy[5]) * sy[6] + sy[3] * (-sy[9] - sy[10] - sy[8])
        - sy[7] * G({ 0, c / d }, 1) - G({ 0, c / d, b / d, 0, 0, a / d }, 1)
        - G({ c / d, 0, b / d, 0, 0, a / d }, 1)
        - 3. * G({ c / d, b / d, 0, 0, 0, a / d }, 1)
        - G({ c / d, b / d, 0, 0, a / d, x / d }, 1)
        + sy[11] * (-G({ x / d }, 1) + G({ d }, { sd }, x))
        + sy[12]
            * (sy[11] + G({ 0, a, 0, 0, b }, { 1, sa, 1, 1, sb }, x)
                + 3. * G({ a, 0, 0, 0, b }, { sa, 1, 1, 1, sb }, x))
        - G({ 0, a, 0, 0, b, c }, { 1, sa, 1, 1, sb, sc }, x)
        - 3. * G({ a, 0, 0, 0, b, c }, { sa, 1, 1, 1, sb, sc }, x)
        - G({ a, 0, 0, b, 0, c }, { sa, 1, 1, sb, 1, sc }, x)
        + (-(sy[2] * sy[3]) / 2. + sy[1] * (sy[4] / 2. + sy[5] / 2.)) * pow(sy[0], 2.)
        + Log(-x, sd)
            * (-sy[13] + sy[0] * (sy[2] * sy[3] + sy[1] * (-sy[4] - sy[5]))
                + sy[3] * (-sy[4] - sy[5]) + sy[2] * sy[6] - sy[12] * sy[7]
                + sy[1] * (sy[9] + sy[10] + sy[8])
                + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.)
        - (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.
        + sy[2]
            * (sy[14] * sy[3] - sy[7] - G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                + 2. * sy[3] * Zeta(2))
        + sy[1]
            * (sy[13] + sy[14] * (-sy[4] - sy[5]) + G({ 0, 0, 0, c / d, b / d }, 1)
                + G({ 0, 0, c / d, 0, b / d }, 1) + G({ 0, c / d, 0, 0, b / d }, 1)
                + G({ c / d, 0, 0, 0, b / d }, 1) - 2. * sy[4] * Zeta(2)
                - 2. * sy[5] * Zeta(2)
                + sy[2]
                    * (-(sy[14] * sy[15]) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                        - 2. * sy[15] * Zeta(2)))
        + sy[0]
            * (sy[13] + sy[3] * (sy[4] + sy[5]) - sy[2] * sy[6] + sy[12] * sy[7]
                + sy[1] * (-sy[9] - sy[10] - sy[8] + sy[2] * (sy[14] + 2. * Zeta(2)))) };
    return res;
}
complex<double> G6_a00bcd_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (c == d) { // abcc
        const vector<complex<double>> sy
            = { G({ a, 0, 0, b }, { sa, 1, 1, sb }, x), G({ b / c, 0, 0, a / c }, 1),
                  G({ c }, { sc }, x), G({ b / c, 0, 0, a / c, 1 }, 1) };
        complex<double> res { -(sy[2] * sy[3])
            + sy[2] * G({ b / c, 0, 0, a / c, x / c }, 1)
            + G({ b / c, 0, 0, a / c, 0, 1 }, 1) - G({ b / c, 0, 0, a / c, x / c, 1 }, 1)
            + (sy[3] - G({ b / c, 0, 0, 0, 1 }, 1)) * G({ a }, { sa }, x)
            + sy[1] * G({ 0, c }, { 1, sc }, x)
            + (-sy[1] + G({ b / c, 0, 0, 1 }, 1)) * G({ a, c }, { sa, sc }, x)
            - G({ b / c, 0, 1 }, 1) * G({ a, 0, c }, { sa, 1, sc }, x)
            + G({ b / c, 1 }, 1) * (-sy[0] + G({ a, 0, 0, c }, { sa, 1, 1, sc }, x))
            + G({ b / c }, 1)
                * (-G({ a, 0, 0, 0, c }, { sa, 1, 1, 1, sc }, x)
                    + G({ a, 0, 0, b, c }, { sa, 1, 1, sb, sc }, x))
            + G({ a, 0, 0, b, 0, c }, { sa, 1, 1, sb, 1, sc }, x) - sy[0] * Zeta(2) };
        return res;
    } else { // abcd
        const vector<complex<double>> sy = { Log(c, sc), G({ a }, { sa }, x),
            G({ b / c, d / c }, 1), G({ d / c, b / c }, 1), G({ 0, 0 }, { 1, 1 }, x),
            G({ 0, a }, { 1, sa }, x), G({ 0, b / c, d / c }, 1),
            G({ 0, d / c, b / c }, 1), G({ d / c, 0, b / c }, 1),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ a, 0, 0, b }, { sa, 1, 1, sb }, x),
            G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x), G({ b / c, 0, 0, a / c }, 1),
            G({ a, d }, { sa, sd }, x), G({ d / c }, 1),
            G({ a, 0, 0, b, d }, { sa, 1, 1, sb, sd }, x), G({ 0, 0, b / c, d / c }, 1),
            G({ 0, 0, d / c, b / c }, 1), G({ 0, d / c, 0, b / c }, 1),
            G({ d / c, 0, 0, b / c }, 1), G({ b / c, 0, 0, a / c, d / c }, 1),
            G({ b / c, 0, 0, d / c, a / c }, 1), G({ b / c, 0, d / c, 0, a / c }, 1),
            G({ b / c, d / c, 0, 0, a / c }, 1), G({ d / c, b / c, 0, 0, a / c }, 1),
            G({ 0, 0, 0 }, { 1, 1, 1 }, x), G({ 0 }, { 1 }, x) };
        complex<double> res { -(sy[12] * sy[13]) + (sy[10] + sy[11]) * sy[3]
            - sy[3] * sy[4] * sy[5] + sy[9] * (-sy[6] - sy[7] - sy[8])
            + sy[10] * G({ 0, d / c }, 1) + sy[13] * G({ b / c, 0, 0, d / c }, 1)
            + G({ 0, b / c, 0, 0, a / c, d / c }, 1)
            + G({ 0, b / c, 0, 0, d / c, a / c }, 1)
            + G({ 0, b / c, 0, d / c, 0, a / c }, 1)
            + G({ 0, b / c, d / c, 0, 0, a / c }, 1)
            + G({ 0, d / c, b / c, 0, 0, a / c }, 1)
            + 3. * G({ b / c, 0, 0, 0, a / c, d / c }, 1)
            + 3. * G({ b / c, 0, 0, 0, d / c, a / c }, 1)
            + G({ b / c, 0, 0, a / c, 0, d / c }, 1)
            + G({ b / c, 0, 0, a / c, d / c, x / c }, 1)
            + 3. * G({ b / c, 0, 0, d / c, 0, a / c }, 1)
            + G({ b / c, 0, 0, d / c, a / c, x / c }, 1)
            + 3. * G({ b / c, 0, d / c, 0, 0, a / c }, 1)
            + G({ b / c, 0, d / c, 0, a / c, x / c }, 1)
            + 3. * G({ b / c, d / c, 0, 0, 0, a / c }, 1)
            + G({ b / c, d / c, 0, 0, a / c, x / c }, 1)
            + G({ d / c, 0, b / c, 0, 0, a / c }, 1)
            + 3. * G({ d / c, b / c, 0, 0, 0, a / c }, 1)
            + G({ d / c, b / c, 0, 0, a / c, x / c }, 1)
            + sy[12] * G({ 0, d }, { 1, sd }, x)
            - G({ b / c, 0, d / c }, 1) * G({ a, 0, d }, { sa, 1, sd }, x)
            + sy[14]
                * (-sy[15] - G({ 0, a, 0, 0, b }, { 1, sa, 1, 1, sb }, x)
                    - 3. * G({ a, 0, 0, 0, b }, { sa, 1, 1, 1, sb }, x))
            + G({ b / c }, 1) * (sy[15] - G({ a, 0, 0, 0, d }, { sa, 1, 1, 1, sd }, x))
            + G({ a, 0, 0, b, 0, d }, { sa, 1, 1, sb, 1, sd }, x)
            + ((sy[2] * sy[5]) / 2. + (sy[3] * sy[5]) / 2.
                  + sy[1] * (-sy[6] / 2. - sy[7] / 2. - sy[8] / 2.))
                * pow(sy[0], 2.)
            + Log(-x, sc)
                * (sy[10] * sy[14] + sy[1] * (-sy[16] - sy[17] - sy[18] - sy[19]) + sy[20]
                    + sy[21] + sy[22] + sy[23] + sy[24] - sy[9] * sy[2] - sy[9] * sy[3]
                    + sy[5] * (sy[6] + sy[7] + sy[8])
                    + sy[0]
                        * (-(sy[2] * sy[5]) - sy[3] * sy[5]
                            + sy[1] * (sy[6] + sy[7] + sy[8]))
                    + sy[1] * (-sy[2] / 2. - sy[3] / 2.) * pow(sy[0], 2.))
            + sy[1] * (sy[2] / 6. + sy[3] / 6.) * pow(sy[0], 3.)
            + sy[5] * (sy[16] + sy[17] + sy[18] + sy[19] - 2. * sy[3] * Zeta(2))
            + sy[2]
                * (sy[11] - sy[4] * sy[5] + G({ a, 0, 0, d }, { sa, 1, 1, sd }, x)
                    - 2. * sy[5] * Zeta(2))
            + sy[1]
                * (-sy[21] - sy[22] - sy[23] - sy[24] - sy[25] * sy[3]
                    + sy[4] * (sy[26] * sy[3] + sy[6] + sy[7] + sy[8])
                    - G({ 0, 0, 0, b / c, d / c }, 1) - G({ 0, 0, 0, d / c, b / c }, 1)
                    - G({ 0, 0, d / c, 0, b / c }, 1) - G({ 0, d / c, 0, 0, b / c }, 1)
                    - G({ b / c, 0, 0, 0, d / c }, 1) - G({ d / c, 0, 0, 0, b / c }, 1)
                    + 2. * sy[26] * sy[3] * Zeta(2) + 2. * sy[6] * Zeta(2)
                    + 2. * sy[7] * Zeta(2) + 2. * sy[8] * Zeta(2)
                    + sy[2] * (-sy[25] + sy[26] * sy[4] + 2. * sy[26] * Zeta(2)))
            + sy[0]
                * (-(sy[10] * sy[14]) - sy[20] - sy[21] - sy[22] - sy[23] - sy[24]
                    + sy[9] * sy[2] + sy[9] * sy[3] + sy[5] * (-sy[6] - sy[7] - sy[8])
                    + sy[1]
                        * (sy[16] + sy[17] + sy[18] + sy[19] - sy[3] * sy[4]
                            + sy[2] * (-sy[4] - 2. * Zeta(2)) - 2. * sy[3] * Zeta(2))) };
        if (d != x) {
            res += (-sy[20] + G({ b / c, 0, 0, a / c, x / c }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_a00bcd_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (b == c) {
        if (b == d) { // abbb
            const vector<complex<double>> sy
                = { G({ a, b, b }, { sa, sb, sb }, x), G({ a, b }, { sa, sb }, x),
                      G({ 0, 0, a / b, 1 }, 1), G({ b, b }, { sb, sb }, x),
                      G({ a }, { sa }, x), G({ 0, 0, a / b, 1, 1 }, 1),
                      G({ b }, { sb }, x), G({ a, 0, b }, { sa, 1, sb }, x) };
            complex<double> res { -(sy[1] * sy[2]) + sy[2] * sy[3] + sy[4] * sy[5]
                - sy[5] * sy[6] - sy[3] * G({ 0, 0, a / b, x / b }, 1)
                + sy[6] * G({ 0, 0, a / b, x / b, 1 }, 1) + G({ 0, 0, a / b, 0, 1, 1 }, 1)
                - G({ 0, 0, a / b, x / b, 1, 1 }, 1)
                + G({ 0, 0, a / b }, 1) * (sy[0] - G({ 0, b, b }, { 1, sb, sb }, x))
                + G({ a, 0, 0, 0, b, b }, { sa, 1, 1, 1, sb, sb }, x)
                + G({ a, 0, 0, b }, { sa, 1, 1, sb }, x) * Zeta(2)
                - G({ a, 0, b, b }, { sa, 1, sb, sb }, x) * Zeta(2)
                + (sy[1] * Zeta(4)) / 4. - sy[7] * Zeta(3) - (-sy[0] + sy[7]) * Zeta(3)
                - sy[4] * (-(Zeta(2) * Zeta(3)) + 2. * Zeta(5)) };
            return res;
        } else { // abbd
            const vector<complex<double>> sy = { Log(b, sb), G({ a }, { sa }, x),
                G({ d / b, 1 }, 1), G({ 0, a }, { 1, sa }, x),
                G({ a, 0, d }, { sa, 1, sd }, x), G({ a, b, d }, { sa, sb, sd }, x),
                G({ b, d }, { sb, sd }, x), G({ 0, 0, a / b, 1 }, 1),
                G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, 0, a / b, d / b, 1 }, 1),
                G({ 0, 0, d / b, 1, a / b }, 1), G({ 0, 0, d / b, a / b, 1 }, 1),
                G({ 0, d / b, 0, 1, a / b }, 1), G({ 0, d / b, 0, a / b, 1 }, 1),
                G({ 0, d / b, 1, 0, a / b }, 1), G({ d / b, 0, 0, 1, a / b }, 1),
                G({ d / b, 0, 0, a / b, 1 }, 1), G({ d / b, 0, 1, 0, a / b }, 1),
                G({ d / b, 1, 0, 0, a / b }, 1), G({ 0, 0 }, { 1, 1 }, x),
                G({ a, 0, 0, d }, { sa, 1, 1, sd }, x), G({ 0 }, { 1 }, x) };
            complex<double> res { sy[6] * sy[7] - sy[4] * G({ 0, d / b, 1 }, 1)
                - sy[6] * G({ 0, 0, a / b, x / b }, 1)
                + 3. * G({ 0, 0, 0, a / b, d / b, 1 }, 1)
                + 3. * G({ 0, 0, 0, d / b, 1, a / b }, 1)
                + 3. * G({ 0, 0, 0, d / b, a / b, 1 }, 1)
                + G({ 0, 0, a / b, 0, d / b, 1 }, 1)
                + G({ 0, 0, a / b, d / b, 1, x / b }, 1)
                + G({ 0, 0, a / b, d / b, x / b, 1 }, 1)
                + 3. * G({ 0, 0, d / b, 0, 1, a / b }, 1)
                + 3. * G({ 0, 0, d / b, 0, a / b, 1 }, 1)
                + 3. * G({ 0, 0, d / b, 1, 0, a / b }, 1)
                + G({ 0, 0, d / b, 1, a / b, x / b }, 1)
                + G({ 0, 0, d / b, a / b, 1, x / b }, 1)
                + G({ 0, 0, d / b, a / b, x / b, 1 }, 1)
                + 3. * G({ 0, d / b, 0, 0, 1, a / b }, 1)
                + 3. * G({ 0, d / b, 0, 0, a / b, 1 }, 1)
                + 3. * G({ 0, d / b, 0, 1, 0, a / b }, 1)
                + G({ 0, d / b, 0, 1, a / b, x / b }, 1)
                + G({ 0, d / b, 0, a / b, 1, x / b }, 1)
                + G({ 0, d / b, 0, a / b, x / b, 1 }, 1)
                + 3. * G({ 0, d / b, 1, 0, 0, a / b }, 1)
                + G({ 0, d / b, 1, 0, a / b, x / b }, 1)
                + 3. * G({ d / b, 0, 0, 0, 1, a / b }, 1)
                + 3. * G({ d / b, 0, 0, 0, a / b, 1 }, 1)
                + 3. * G({ d / b, 0, 0, 1, 0, a / b }, 1)
                + G({ d / b, 0, 0, 1, a / b, x / b }, 1)
                + G({ d / b, 0, 0, a / b, 1, x / b }, 1)
                + G({ d / b, 0, 0, a / b, x / b, 1 }, 1)
                + 3. * G({ d / b, 0, 1, 0, 0, a / b }, 1)
                + G({ d / b, 0, 1, 0, a / b, x / b }, 1)
                + 3. * G({ d / b, 1, 0, 0, 0, a / b }, 1)
                + G({ d / b, 1, 0, 0, a / b, x / b }, 1)
                + (-sy[7] + G({ 0, 0, d / b, 1 }, 1)) * G({ a, d }, { sa, sd }, x)
                + G({ 0, 0, a / b }, 1) * (sy[5] - G({ 0, b, d }, { 1, sb, sd }, x))
                + G({ a, 0, 0, 0, b, d }, { sa, 1, 1, 1, sb, sd }, x)
                + (sy[2] * sy[3] * pow(sy[0], 2.)) / 2.
                + Log(-x, sb)
                    * (sy[9] + sy[10] + sy[11] + sy[12] + sy[13] + sy[14] + sy[15]
                        + sy[16] + sy[17] + sy[18] - sy[0] * sy[2] * sy[3] - sy[2] * sy[8]
                        - (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.)
                + (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.
                + sy[0]
                    * (-sy[9] - sy[10] - sy[11] - sy[12] - sy[13] - sy[14] - sy[15]
                        - sy[16] - sy[17] - sy[18] + sy[2] * sy[8]
                        + sy[1] * sy[2] * (-sy[19] - 2. * Zeta(2)))
                + sy[20] * Zeta(2) - G({ a, 0, b, d }, { sa, 1, sb, sd }, x) * Zeta(2)
                + sy[2]
                    * (sy[20] - sy[19] * sy[3] + G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                        - 2. * sy[3] * Zeta(2))
                + sy[1]
                    * (-sy[10] - sy[11] - sy[12] - sy[13] - sy[14] - sy[15] - sy[16]
                        - sy[17] - sy[18] - G({ 0, 0, 0, d / b, 1 }, 1)
                        + sy[2]
                            * (sy[19] * sy[21] - G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                                + 2. * sy[21] * Zeta(2)))
                - sy[4] * Zeta(3) + sy[5] * Zeta(3) };
            if (d != x) {
                res += (-sy[9] + G({ 0, 0, a / b, x / b, 1 }, 1)) * G({ d }, { sd }, x);
            }
            return res;
        }
    } else { // abcd
        const vector<complex<double>> sy = { Log(b, sb), G({ a }, { sa }, x),
            G({ c / b, d / b }, 1), G({ 0, a }, { 1, sa }, x),
            G({ a, 0, d }, { sa, 1, sd }, x), G({ a, c, d }, { sa, sc, sd }, x),
            G({ a, d }, { sa, sd }, x), G({ 0, 0, c / b, a / b }, 1),
            G({ a, 0, 0, d }, { sa, 1, 1, sd }, x), G({ 0, 0, a / b, c / b }, 1),
            G({ 0, c / b, 0, a / b }, 1), G({ c / b, 0, 0, a / b }, 1),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, 0, a / b, c / b, d / b }, 1),
            G({ 0, 0, c / b, a / b, d / b }, 1), G({ 0, 0, c / b, d / b, a / b }, 1),
            G({ 0, c / b, 0, a / b, d / b }, 1), G({ 0, c / b, 0, d / b, a / b }, 1),
            G({ 0, c / b, d / b, 0, a / b }, 1), G({ c / b, 0, 0, a / b, d / b }, 1),
            G({ c / b, 0, 0, d / b, a / b }, 1), G({ c / b, 0, d / b, 0, a / b }, 1),
            G({ c / b, d / b, 0, 0, a / b }, 1), G({ 0, 0 }, { 1, 1 }, x),
            G({ 0 }, { 1 }, x) };
        complex<double> res { sy[6] * sy[7] + (sy[4] - sy[5]) * G({ 0, 0, c / b }, 1)
            + sy[4] * (G({ 0, c / b, d / b }, 1) + G({ c / b, 0, d / b }, 1))
            + sy[6]
                * (sy[10] + sy[11] - G({ 0, 0, c / b, d / b }, 1)
                    - G({ 0, c / b, 0, d / b }, 1) - G({ c / b, 0, 0, d / b }, 1))
            - 3. * G({ 0, 0, 0, a / b, c / b, d / b }, 1)
            - 3. * G({ 0, 0, 0, c / b, a / b, d / b }, 1)
            - 3. * G({ 0, 0, 0, c / b, d / b, a / b }, 1)
            - G({ 0, 0, a / b, 0, c / b, d / b }, 1)
            - G({ 0, 0, a / b, c / b, 0, d / b }, 1)
            - G({ 0, 0, a / b, c / b, d / b, x / b }, 1)
            - 3. * G({ 0, 0, c / b, 0, a / b, d / b }, 1)
            - 3. * G({ 0, 0, c / b, 0, d / b, a / b }, 1)
            - G({ 0, 0, c / b, a / b, 0, d / b }, 1)
            - G({ 0, 0, c / b, a / b, d / b, x / b }, 1)
            - 3. * G({ 0, 0, c / b, d / b, 0, a / b }, 1)
            - G({ 0, 0, c / b, d / b, a / b, x / b }, 1)
            - 3. * G({ 0, c / b, 0, 0, a / b, d / b }, 1)
            - 3. * G({ 0, c / b, 0, 0, d / b, a / b }, 1)
            - G({ 0, c / b, 0, a / b, 0, d / b }, 1)
            - G({ 0, c / b, 0, a / b, d / b, x / b }, 1)
            - 3. * G({ 0, c / b, 0, d / b, 0, a / b }, 1)
            - G({ 0, c / b, 0, d / b, a / b, x / b }, 1)
            - 3. * G({ 0, c / b, d / b, 0, 0, a / b }, 1)
            - G({ 0, c / b, d / b, 0, a / b, x / b }, 1)
            - 3. * G({ c / b, 0, 0, 0, a / b, d / b }, 1)
            - 3. * G({ c / b, 0, 0, 0, d / b, a / b }, 1)
            - G({ c / b, 0, 0, a / b, 0, d / b }, 1)
            - G({ c / b, 0, 0, a / b, d / b, x / b }, 1)
            - 3. * G({ c / b, 0, 0, d / b, 0, a / b }, 1)
            - G({ c / b, 0, 0, d / b, a / b, x / b }, 1)
            - 3. * G({ c / b, 0, d / b, 0, 0, a / b }, 1)
            - G({ c / b, 0, d / b, 0, a / b, x / b }, 1)
            - 3. * G({ c / b, d / b, 0, 0, 0, a / b }, 1)
            - G({ c / b, d / b, 0, 0, a / b, x / b }, 1)
            + (-sy[9] - sy[10] - sy[11] - sy[7]) * G({ 0, d }, { 1, sd }, x)
            + G({ 0, 0, a / b }, 1) * (sy[5] - G({ 0, c, d }, { 1, sc, sd }, x))
            + G({ 0, c / b }, 1) * (-sy[8] + G({ a, 0, c, d }, { sa, 1, sc, sd }, x))
            + G({ c / b }, 1)
                * (G({ a, 0, 0, 0, d }, { sa, 1, 1, 1, sd }, x)
                    - G({ a, 0, 0, c, d }, { sa, 1, 1, sc, sd }, x))
            + G({ a, 0, 0, 0, c, d }, { sa, 1, 1, 1, sc, sd }, x)
            - (sy[2] * sy[3] * pow(sy[0], 2.)) / 2.
            + Log(-x, sb)
                * (-sy[13] - sy[14] - sy[15] - sy[16] - sy[17] - sy[18] - sy[19] - sy[20]
                    - sy[21] - sy[22] + sy[12] * sy[2] + sy[0] * sy[2] * sy[3]
                    + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.)
            - (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.
            + sy[2]
                * (sy[23] * sy[3] - sy[8] - G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x)
                    + 2. * sy[3] * Zeta(2))
            + sy[0]
                * (sy[13] + sy[14] + sy[15] + sy[16] + sy[17] + sy[18] + sy[19] + sy[20]
                    + sy[21] + sy[22] - sy[12] * sy[2]
                    + sy[1] * sy[2] * (sy[23] + 2. * Zeta(2)))
            + sy[1]
                * (sy[15] + sy[17] + sy[18] + sy[20] + sy[21] + sy[22]
                    + G({ 0, 0, 0, c / b, d / b }, 1) + G({ 0, 0, c / b, 0, d / b }, 1)
                    + G({ 0, c / b, 0, 0, d / b }, 1) + G({ c / b, 0, 0, 0, d / b }, 1)
                    + sy[2]
                        * (-(sy[23] * sy[24]) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                            - 2. * sy[24] * Zeta(2))) };
        if (c != x) {
            res += (sy[9] - G({ 0, 0, a / b, x / b }, 1)) * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (sy[13] + sy[14] + sy[16] + sy[19]
                       - G({ 0, 0, a / b, c / b, x / b }, 1)
                       - G({ 0, 0, c / b, a / b, x / b }, 1)
                       - G({ 0, c / b, 0, a / b, x / b }, 1)
                       - G({ c / b, 0, 0, a / b, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_a00bcd_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ 0, 0, b / a }, 1),
        G({ 0, 0, b / a, c / a }, 1), G({ 0, 0, b / a, c / a, d / a }, 1) };
    complex<double> res { 3. * G({ 0, 0, 0, b / a, c / a, d / a }, 1)
        + G({ 0, 0, b / a, 0, c / a, d / a }, 1) + G({ 0, 0, b / a, c / a, 0, d / a }, 1)
        + G({ 0, 0, b / a, c / a, d / a, x / a }, 1) + sy[1] * G({ 0, d }, { 1, sd }, x)
        + sy[0] * G({ 0, c, d }, { 1, sc, sd }, x)
        + G({ 0, x / a }, 1) * G({ 0, b, c, d }, { 1, sb, sc, sd }, x)
        + G({ x / a }, 1) * G({ 0, 0, b, c, d }, { 1, 1, sb, sc, sd }, x)
        + G({ 0, 0, 0, b, c, d }, { 1, 1, 1, sb, sc, sd }, x) - sy[2] * Log(a, sa)
        + sy[2] * Log(-x, sa) };
    if (b != x) {
        res += (-sy[0] + G({ 0, 0, x / a }, 1)) * G({ b, c, d }, { sb, sc, sd }, x);
    }
    if (c != x) {
        res += (-sy[1] + G({ 0, 0, b / a, x / a }, 1)) * G({ c, d }, { sc, sd }, x);
    }
    if (d != x) {
        res += (-sy[2] + G({ 0, 0, b / a, c / a, x / a }, 1)) * G({ d }, { sd }, x);
    }
    return res;
}

complex<double> G6_ab00cd_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { Log(d, sd), G({ c / d }, 1),
        G({ a, b }, { sa, sb }, x), G({ 0, c / d }, 1), G({ 0, a, b }, { 1, sa, sb }, x),
        G({ a, 0, b }, { sa, 1, sb }, x), G({ 0, 0, c / d }, 1),
        G({ c / d, 0, x / d }, 1), G({ 0, 0, a, b }, { 1, 1, sa, sb }, x),
        G({ 0, a, 0, b }, { 1, sa, 1, sb }, x), G({ a, 0, 0, b }, { sa, 1, 1, sb }, x),
        G({ c / d, x / d }, 1), G({ c / d, 0, 0, b / d }, 1),
        G({ a, b, 0, 0, c }, { sa, sb, 1, 1, sc }, x), G({ a }, { sa }, x),
        G({ c / d, 0, 0, b / d, a / d }, 1), G({ 0, 0 }, { 1, 1 }, x),
        G({ c / d, 0, 0, x / d }, 1), G({ 0 }, { 1 }, x) };
    complex<double> res { (sy[4] + sy[5]) * sy[6] - sy[4] * sy[7] - sy[5] * sy[7]
        + sy[3] * (-sy[9] - sy[10] - sy[8]) - G({ 0, c / d, 0, 0, b / d, a / d }, 1)
        - 3. * G({ c / d, 0, 0, 0, b / d, a / d }, 1)
        - G({ c / d, 0, 0, b / d, 0, a / d }, 1)
        - G({ c / d, 0, 0, b / d, a / d, x / d }, 1)
        + sy[13] * (-G({ x / d }, 1) + G({ d }, { sd }, x))
        - sy[12] * G({ 0, a }, { 1, sa }, x)
        + sy[14]
            * (sy[15] + G({ 0, c / d, 0, 0, b / d }, 1)
                + 3. * G({ c / d, 0, 0, 0, b / d }, 1) + sy[7] * G({ 0, b }, { 1, sb }, x)
                + sy[11] * G({ 0, 0, b }, { 1, 1, sb }, x))
        - sy[7] * G({ 0, b, a }, { 1, sb, sa }, x)
        + sy[11] * (-sy[9] - sy[10] - sy[8] - G({ 0, 0, b, a }, { 1, 1, sb, sa }, x))
        - G({ 0, a, b, 0, 0, c }, { 1, sa, sb, 1, 1, sc }, x)
        - G({ a, 0, b, 0, 0, c }, { sa, 1, sb, 1, 1, sc }, x)
        - 3. * G({ a, b, 0, 0, 0, c }, { sa, sb, 1, 1, 1, sc }, x)
        + (-(sy[2] * sy[3]) / 2. + sy[1] * (sy[4] / 2. + sy[5] / 2.)) * pow(sy[0], 2.)
        + Log(-x, sd)
            * (sy[12] * sy[14] - sy[15]
                + sy[0] * (sy[2] * sy[3] + sy[1] * (-sy[4] - sy[5]))
                + sy[3] * (sy[4] + sy[5]) - sy[2] * sy[6]
                + sy[1] * (-sy[9] - sy[10] - sy[8])
                - (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.)
        + (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.
        + sy[0]
            * (-(sy[12] * sy[14]) + sy[15] + sy[3] * (-sy[4] - sy[5]) + sy[2] * sy[6]
                + sy[1] * (sy[9] + sy[10] + sy[8] + sy[2] * (-sy[16] - 2. * Zeta(2))))
        + sy[2]
            * (-sy[17] + sy[16] * sy[3] - G({ 0, 0, 0, c / d }, 1) + 2. * sy[3] * Zeta(2))
        + sy[1]
            * (sy[13] + sy[16] * (-sy[4] - sy[5])
                + sy[14] * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
                - G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x) - 2. * sy[4] * Zeta(2)
                - 2. * sy[5] * Zeta(2)
                + sy[2]
                    * (sy[16] * sy[18] - G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                        + 2. * sy[18] * Zeta(2))) };
    if (b != x) {
        res += (-(sy[12] * sy[14]) + sy[14] * sy[17]) * G({ b }, { sb }, x)
            + (sy[12] - sy[17]) * G({ b, a }, { sb, sa }, x);
    }
    return res;
}
complex<double> G6_ab00cd_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (c == d) { // abcc
        const vector<complex<double>> sy
            = { G({ 0, x / c, 1 }, 1), G({ a, b, c }, { sa, sb, sc }, x),
                  G({ 0, 0, b / c, 1 }, 1), G({ 0, 0, b / c, a / c }, 1),
                  G({ x / c, 1 }, 1), G({ 0, 0, a, b }, { 1, 1, sa, sb }, x),
                  G({ 0, 0, b, a }, { 1, 1, sb, sa }, x),
                  G({ 0, a, 0, b }, { 1, sa, 1, sb }, x),
                  G({ a, 0, 0, b }, { sa, 1, 1, sb }, x), G({ c }, { sc }, x),
                  G({ 0, 0, b / c, a / c, 1 }, 1), G({ a }, { sa }, x),
                  G({ 0, 0, b }, { 1, 1, sb }, x), G({ 0, 0, 1, b / c }, 1),
                  G({ 0, 0, 1, x / c }, 1), G({ 0, 0, x / c, 1 }, 1) };
        complex<double> res { -(sy[9] * sy[10]) + sy[4] * (-sy[5] - sy[6] - sy[7] - sy[8])
            + sy[9] * G({ 0, 0, b / c, a / c, x / c }, 1)
            + G({ 0, 0, b / c, a / c, 0, 1 }, 1) - G({ 0, 0, b / c, a / c, x / c, 1 }, 1)
            + sy[3] * G({ 0, c }, { 1, sc }, x)
            + (sy[2] - sy[3]) * G({ a, c }, { sa, sc }, x)
            + sy[0]
                * (-G({ 0, a, b }, { 1, sa, sb }, x) - G({ 0, b, a }, { 1, sb, sa }, x)
                    - G({ a, 0, b }, { sa, 1, sb }, x))
            + G({ 0, 0, b / c }, 1) * (sy[1] - G({ a, 0, c }, { sa, 1, sc }, x))
            + G({ a, b, 0, 0, 0, c }, { sa, sb, 1, 1, 1, sc }, x) - sy[5] * Zeta(2)
            - sy[6] * Zeta(2) - sy[7] * Zeta(2) - sy[8] * Zeta(2)
            - G({ a, b, 0, c }, { sa, sb, 1, sc }, x) * Zeta(2)
            + sy[11]
                * (sy[10] + sy[12] * sy[4] - G({ 0, 0, b / c, 0, 1 }, 1)
                    + sy[0] * G({ 0, b }, { 1, sb }, x) + sy[12] * Zeta(2))
            + G({ a, b }, { sa, sb }, x) * (sy[13] - sy[14] - sy[15] - Zeta(4))
            + sy[1] * Zeta(3) };
        if (b != x) {
            res += (-(sy[11] * sy[13]) + sy[11] * sy[14] + sy[11] * sy[15]
                       - sy[11] * sy[2])
                    * G({ b }, { sb }, x)
                + (sy[13] - sy[14] - sy[15] + sy[2]) * G({ b, a }, { sb, sa }, x);
        }
        return res;
    } else { // abcd
        const vector<complex<double>> sy = { Log(c, sc), G({ d / c }, 1),
            G({ a, b }, { sa, sb }, x), G({ 0, a, b }, { 1, sa, sb }, x),
            G({ a, 0, b }, { sa, 1, sb }, x), G({ 0, d / c, x / c }, 1),
            G({ 0, b, a }, { 1, sb, sa }, x), G({ a, b, d }, { sa, sb, sd }, x),
            G({ d / c, 0, x / c }, 1), G({ 0, 0, b / c, a / c }, 1),
            G({ a, d }, { sa, sd }, x), G({ 0, a }, { 1, sa }, x),
            G({ 0, 0, b / c, d / c }, 1), G({ d / c, x / c }, 1),
            G({ 0, 0, a, b }, { 1, 1, sa, sb }, x),
            G({ 0, 0, b, a }, { 1, 1, sb, sa }, x),
            G({ 0, a, 0, b }, { 1, sa, 1, sb }, x),
            G({ a, 0, 0, b }, { sa, 1, 1, sb }, x), G({ 0, d / c }, 1),
            G({ 0, 0, d / c, b / c }, 1), G({ 0, d / c, 0, b / c }, 1),
            G({ d / c, 0, 0, b / c }, 1), G({ 0, d / c, 0, x / c }, 1),
            G({ d / c, 0, 0, x / c }, 1), G({ a }, { sa }, x),
            G({ 0, 0, b }, { 1, 1, sb }, x), G({ 0, 0, b / c, d / c, a / c }, 1),
            G({ 0, 0, d / c, b / c, a / c }, 1), G({ 0, d / c, 0, b / c, a / c }, 1),
            G({ d / c, 0, 0, b / c, a / c }, 1), G({ 0, 0, b / c, a / c, d / c }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0 }, { 1 }, x) };
        complex<double> res {
            -(sy[9] * sy[10]) + sy[10] * sy[12] + sy[11] * sy[12]
            + sy[13] * (sy[14] + sy[15] + sy[16] + sy[17])
            + sy[11] * (sy[19] + sy[20] + sy[21]) + sy[5] * (sy[4] + sy[6])
            + sy[4] * sy[8] + sy[6] * sy[8] + sy[3] * (sy[5] + sy[8])
            - sy[7] * G({ 0, 0, d / c }, 1)
            + sy[2] * (sy[19] + sy[22] + sy[23] + G({ 0, 0, 0, d / c }, 1))
            + 3. * G({ 0, 0, 0, b / c, a / c, d / c }, 1)
            + 3. * G({ 0, 0, 0, b / c, d / c, a / c }, 1)
            + 3. * G({ 0, 0, 0, d / c, b / c, a / c }, 1)
            + G({ 0, 0, b / c, 0, a / c, d / c }, 1)
            + G({ 0, 0, b / c, 0, d / c, a / c }, 1)
            + G({ 0, 0, b / c, a / c, 0, d / c }, 1)
            + G({ 0, 0, b / c, a / c, d / c, x / c }, 1)
            + G({ 0, 0, b / c, d / c, 0, a / c }, 1)
            + G({ 0, 0, b / c, d / c, a / c, x / c }, 1)
            + 3. * G({ 0, 0, d / c, 0, b / c, a / c }, 1)
            + G({ 0, 0, d / c, b / c, 0, a / c }, 1)
            + G({ 0, 0, d / c, b / c, a / c, x / c }, 1)
            + 3. * G({ 0, d / c, 0, 0, b / c, a / c }, 1)
            + G({ 0, d / c, 0, b / c, 0, a / c }, 1)
            + G({ 0, d / c, 0, b / c, a / c, x / c }, 1)
            + 3. * G({ d / c, 0, 0, 0, b / c, a / c }, 1)
            + G({ d / c, 0, 0, b / c, 0, a / c }, 1)
            + G({ d / c, 0, 0, b / c, a / c, x / c }, 1)
            + sy[24]
                * (-(sy[13] * sy[25]) - sy[18] * sy[25] - sy[26] - sy[27] - sy[28]
                    - sy[29] - 3. * G({ 0, 0, 0, b / c, d / c }, 1)
                    - 3. * G({ 0, 0, 0, d / c, b / c }, 1)
                    - G({ 0, 0, b / c, 0, d / c }, 1)
                    - 3. * G({ 0, 0, d / c, 0, b / c }, 1)
                    - 3. * G({ 0, d / c, 0, 0, b / c }, 1)
                    - 3. * G({ d / c, 0, 0, 0, b / c }, 1)
                    + (-sy[5] - sy[8]) * G({ 0, b }, { 1, sb }, x))
            + sy[9] * G({ 0, d }, { 1, sd }, x)
            + G({ 0, 0, b / c }, 1) * (sy[7] - G({ a, 0, d }, { sa, 1, sd }, x))
            + sy[18]
                * (sy[14] + sy[15] + sy[16] + sy[17]
                    + G({ a, b, 0, d }, { sa, sb, 1, sd }, x))
            + G({ a, b, 0, 0, 0, d }, { sa, sb, 1, 1, 1, sd }, x)
            + sy[1] * (-sy[3] / 2. - sy[4] / 2.) * pow(sy[0], 2.)
            + Log(-x, sc)
                * ((sy[14] + sy[16] + sy[17]) * sy[1]
                    + (-sy[12] - sy[19] - sy[20] - sy[21]) * sy[24] + sy[26] + sy[27]
                    + sy[28] + sy[29] + sy[30] + sy[0] * sy[1] * (sy[3] + sy[4])
                    + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.)
            - (sy[1] * sy[2] * pow(sy[0], 3.)) / 6.
            + sy[1]
                * (sy[31] * (sy[3] + sy[4])
                    - sy[24] * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
                    + G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x)
                    - G({ a, b, 0, 0, d }, { sa, sb, 1, 1, sd }, x) + 2. * sy[3] * Zeta(2)
                    + 2. * sy[4] * Zeta(2)
                    + sy[2]
                        * (-(sy[31] * sy[32]) + G({ 0, 0, 0 }, { 1, 1, 1 }, x)
                            - 2. * sy[32] * Zeta(2)))
            + sy[0]
                * ((sy[12] + sy[19] + sy[20] + sy[21]) * sy[24] - sy[26] - sy[27] - sy[28]
                    - sy[29] - sy[30]
                    + sy[1]
                        * (-sy[14] - sy[16] - sy[17] + sy[2] * (sy[31] + 2. * Zeta(2))))
        };
        if (b != x) {
            res += (sy[20] * sy[24] + sy[21] * sy[24] - sy[22] * sy[24] - sy[23] * sy[24])
                    * G({ b }, { sb }, x)
                + (-sy[20] - sy[21] + sy[22] + sy[23]) * G({ b, a }, { sb, sa }, x);
        }
        if (d != x) {
            res += (-sy[30] + G({ 0, 0, b / c, a / c, x / c }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_ab00cd_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ a, c, d }, { sa, sc, sd }, x),
        G({ 0, c, d }, { 1, sc, sd }, x), G({ 0, 0, c, d }, { 1, 1, sc, sd }, x),
        G({ a, d }, { sa, sd }, x), G({ 0, 0, c / b, a / b }, 1),
        G({ 0, 0, c / b, d / b }, 1), G({ 0, 0, a / b, c / b }, 1),
        G({ 0, a / b, 0, c / b }, 1), G({ a / b, 0, 0, c / b }, 1), G({ a }, { sa }, x),
        G({ 0, 0, c / b, d / b, a / b }, 1), G({ 0, 0, a / b, c / b, d / b }, 1),
        G({ 0, 0, c / b, a / b, d / b }, 1), G({ 0, a / b, 0, c / b, d / b }, 1),
        G({ a / b, 0, 0, c / b, d / b }, 1) };
    complex<double> res { sy[3] * sy[4] - sy[3] * sy[5] - sy[2] * G({ a / b, x / b }, 1)
        + (sy[0] - sy[1]) * G({ 0, 0, a / b }, 1)
        + sy[1] * (-G({ 0, a / b, x / b }, 1) - G({ a / b, 0, x / b }, 1))
        + sy[9]
            * (sy[10] + 3. * G({ 0, 0, 0, c / b, d / b }, 1)
                + G({ 0, 0, c / b, 0, d / b }, 1))
        - 3. * G({ 0, 0, 0, a / b, c / b, d / b }, 1)
        - 3. * G({ 0, 0, 0, c / b, a / b, d / b }, 1)
        - 3. * G({ 0, 0, 0, c / b, d / b, a / b }, 1)
        - 3. * G({ 0, 0, a / b, 0, c / b, d / b }, 1)
        - G({ 0, 0, a / b, c / b, 0, d / b }, 1)
        - G({ 0, 0, a / b, c / b, d / b, x / b }, 1)
        - G({ 0, 0, c / b, 0, a / b, d / b }, 1) - G({ 0, 0, c / b, 0, d / b, a / b }, 1)
        - G({ 0, 0, c / b, a / b, 0, d / b }, 1)
        - G({ 0, 0, c / b, a / b, d / b, x / b }, 1)
        - G({ 0, 0, c / b, d / b, 0, a / b }, 1)
        - G({ 0, 0, c / b, d / b, a / b, x / b }, 1)
        - 3. * G({ 0, a / b, 0, 0, c / b, d / b }, 1)
        - G({ 0, a / b, 0, c / b, 0, d / b }, 1)
        - G({ 0, a / b, 0, c / b, d / b, x / b }, 1)
        - 3. * G({ a / b, 0, 0, 0, c / b, d / b }, 1)
        - G({ a / b, 0, 0, c / b, 0, d / b }, 1)
        - G({ a / b, 0, 0, c / b, d / b, x / b }, 1) - sy[5] * G({ 0, a }, { 1, sa }, x)
        + (-sy[4] - sy[6] - sy[7] - sy[8]) * G({ 0, d }, { 1, sd }, x)
        + G({ 0, 0, c / b }, 1) * (-sy[0] + G({ a, 0, d }, { sa, 1, sd }, x))
        + G({ 0, a / b }, 1) * (-sy[2] + G({ a, 0, c, d }, { sa, 1, sc, sd }, x))
        + G({ a / b }, 1)
            * (-G({ 0, 0, 0, c, d }, { 1, 1, 1, sc, sd }, x)
                + G({ a, 0, 0, c, d }, { sa, 1, 1, sc, sd }, x))
        + G({ a, 0, 0, 0, c, d }, { sa, 1, 1, 1, sc, sd }, x)
        + (sy[10] + sy[11] + sy[12] + sy[13] + sy[14] - sy[9] * sy[5]) * Log(b, sb)
        + (-sy[10] - sy[11] - sy[12] - sy[13] - sy[14] + sy[9] * sy[5]) * Log(-x, sb) };
    if (c != x) {
        res += (sy[6] + sy[7] + sy[8] - G({ 0, 0, a / b, x / b }, 1)
                   - G({ 0, a / b, 0, x / b }, 1) - G({ a / b, 0, 0, x / b }, 1))
            * G({ c, d }, { sc, sd }, x);
    }
    if (d != x) {
        res += (sy[11] + sy[12] + sy[13] + sy[14] - G({ 0, 0, a / b, c / b, x / b }, 1)
                   - G({ 0, 0, c / b, a / b, x / b }, 1)
                   - G({ 0, a / b, 0, c / b, x / b }, 1)
                   - G({ a / b, 0, 0, c / b, x / b }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}
complex<double> G6_ab00cd_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (a == b) {
        const vector<complex<double>> sy = { G({ 0, c, d }, { 1, sc, sd }, x),
            G({ 0, 0, 1, c / a }, 1), G({ 0, 0, c / a, 1 }, 1), G({ 0, 1, 0, c / a }, 1),
            G({ 0, 0, 1, c / a, d / a }, 1), G({ 0, 0, c / a, 1, d / a }, 1),
            G({ 0, 0, c / a, d / a, 1 }, 1), G({ 0, 1, 0, c / a, d / a }, 1) };
        complex<double> res { -(sy[0] * G({ 0, 1, x / a }, 1))
            - sy[0] * G({ 0, x / a, 1 }, 1) - 3. * G({ 0, 0, 0, 1, c / a, d / a }, 1)
            - 3. * G({ 0, 0, 0, c / a, 1, d / a }, 1)
            - 3. * G({ 0, 0, 0, c / a, d / a, 1 }, 1)
            - 3. * G({ 0, 0, 1, 0, c / a, d / a }, 1) - G({ 0, 0, 1, c / a, 0, d / a }, 1)
            - G({ 0, 0, 1, c / a, d / a, x / a }, 1) - G({ 0, 0, c / a, 0, 1, d / a }, 1)
            - G({ 0, 0, c / a, 0, d / a, 1 }, 1) - G({ 0, 0, c / a, 1, 0, d / a }, 1)
            - G({ 0, 0, c / a, 1, d / a, x / a }, 1)
            - G({ 0, 0, c / a, d / a, 1, x / a }, 1)
            - G({ 0, 0, c / a, d / a, x / a, 1 }, 1)
            - 2. * G({ 0, 1, 0, 0, c / a, d / a }, 1) - G({ 0, 1, 0, c / a, 0, d / a }, 1)
            - G({ 0, 1, 0, c / a, d / a, x / a }, 1)
            + (-sy[1] - sy[2] - sy[3]) * G({ 0, d }, { 1, sd }, x)
            - G({ x / a, 1 }, 1) * G({ 0, 0, c, d }, { 1, 1, sc, sd }, x)
            + G({ x / a }, 1) * G({ a, 0, 0, c, d }, { sa, 1, 1, sc, sd }, x)
            + G({ 0, a, 0, 0, c, d }, { 1, sa, 1, 1, sc, sd }, x)
            + (sy[4] + sy[5] + sy[6] + sy[7]) * Log(a, sa)
            + (-sy[4] - sy[5] - sy[6] - sy[7]) * Log(-x, sa) };
        if (c != x) {
            res += (sy[1] + sy[2] + sy[3] - G({ 0, 0, 1, x / a }, 1)
                       - G({ 0, 0, x / a, 1 }, 1) - G({ 0, 1, 0, x / a }, 1))
                * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (sy[4] + sy[5] + sy[6] + sy[7] - G({ 0, 0, 1, c / a, x / a }, 1)
                       - G({ 0, 0, c / a, 1, x / a }, 1) - G({ 0, 0, c / a, x / a, 1 }, 1)
                       - G({ 0, 1, 0, c / a, x / a }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    } else { // abcd
        const vector<complex<double>> sy = { G({ b / a, 0, 0, c / a }, 1),
            G({ b / a }, 1), G({ b / a, 0, 0, c / a, d / a }, 1) };
        complex<double> res { G({ 0, b / a, 0, 0, c / a, d / a }, 1)
            + 3. * G({ b / a, 0, 0, 0, c / a, d / a }, 1)
            + G({ b / a, 0, 0, c / a, 0, d / a }, 1)
            + G({ b / a, 0, 0, c / a, d / a, x / a }, 1)
            + sy[0] * G({ 0, d }, { 1, sd }, x)
            + G({ b / a, 0, x / a }, 1) * G({ 0, c, d }, { 1, sc, sd }, x)
            + G({ b / a, x / a }, 1) * G({ 0, 0, c, d }, { 1, 1, sc, sd }, x)
            + sy[1] * G({ 0, 0, 0, c, d }, { 1, 1, 1, sc, sd }, x)
            + G({ 0, b, 0, 0, c, d }, { 1, sb, 1, 1, sc, sd }, x) - sy[2] * Log(a, sa)
            + sy[2] * Log(-x, sa) };
        if (b != x) {
            res += (-sy[1] + G({ x / a }, 1))
                * G({ b, 0, 0, c, d }, { sb, 1, 1, sc, sd }, x);
        }
        if (c != x) {
            res += (-sy[0] + G({ b / a, 0, 0, x / a }, 1)) * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (-sy[2] + G({ b / a, 0, 0, c / a, x / a }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G6_abc00d_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { G({ a, b, c }, { sa, sb, sc }, x),
        G({ 0, 0, c / d }, 1), G({ a, b }, { sa, sb }, x), G({ 0, 0, c / d, b / d }, 1),
        G({ 0, x / d }, 1), G({ 0, a, b, c }, { 1, sa, sb, sc }, x),
        G({ a, 0, b, c }, { sa, 1, sb, sc }, x), G({ a, b, 0, c }, { sa, sb, 1, sc }, x),
        G({ x / d }, 1), G({ 0, 0, a, b, c }, { 1, 1, sa, sb, sc }, x),
        G({ a }, { sa }, x), G({ 0, 0, c / d, b / d, a / d }, 1),
        G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x),
        G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x),
        G({ a, 0, 0, b, c }, { sa, 1, 1, sb, sc }, x),
        G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x),
        G({ a, b, 0, 0, c }, { sa, sb, 1, 1, sc }, x) };
    complex<double> res { -(sy[4] * sy[5]) + sy[4] * (-sy[6] - sy[7]) - sy[9] * sy[8]
        + (-sy[12] - sy[13] - sy[14] - sy[15] - sy[16]) * sy[8]
        + sy[2] * (-sy[3] - 3. * G({ 0, 0, 0, c / d }, 1))
        + sy[10]
            * (sy[11] + 3. * G({ 0, 0, 0, c / d, b / d }, 1)
                + G({ 0, 0, c / d, 0, b / d }, 1))
        - 3. * G({ 0, 0, 0, c / d, b / d, a / d }, 1)
        - G({ 0, 0, c / d, 0, b / d, a / d }, 1) - G({ 0, 0, c / d, b / d, 0, a / d }, 1)
        - G({ 0, 0, c / d, b / d, a / d, x / d }, 1)
        + (sy[9] + sy[12] + sy[13] + sy[14] + sy[15] + sy[16]) * G({ d }, { sd }, x)
        - sy[3] * G({ 0, a }, { 1, sa }, x)
        + (-sy[5] - sy[6] - sy[7]) * G({ 0, d }, { 1, sd }, x)
        + sy[0] * (-G({ 0, 0, x / d }, 1) + G({ 0, 0, d }, { 1, 1, sd }, x))
        + sy[1]
            * (sy[0] + G({ 0, a, b }, { 1, sa, sb }, x)
                + G({ a, 0, b }, { sa, 1, sb }, x))
        - G({ 0, 0, 0, a, b, c }, { 1, 1, 1, sa, sb, sc }, x)
        - G({ 0, 0, a, 0, b, c }, { 1, 1, sa, 1, sb, sc }, x)
        - G({ 0, 0, a, b, 0, c }, { 1, 1, sa, sb, 1, sc }, x)
        - G({ 0, a, 0, 0, b, c }, { 1, sa, 1, 1, sb, sc }, x)
        - G({ 0, a, 0, b, 0, c }, { 1, sa, 1, sb, 1, sc }, x)
        - G({ 0, a, b, 0, 0, c }, { 1, sa, sb, 1, 1, sc }, x)
        - G({ a, 0, 0, 0, b, c }, { sa, 1, 1, 1, sb, sc }, x)
        - G({ a, 0, 0, b, 0, c }, { sa, 1, 1, sb, 1, sc }, x)
        - G({ a, 0, b, 0, 0, c }, { sa, 1, sb, 1, 1, sc }, x)
        - G({ a, b, 0, 0, 0, c }, { sa, sb, 1, 1, 1, sc }, x)
        + (sy[11] + sy[1] * sy[2] - sy[10] * sy[3]) * Log(d, sd)
        + (-sy[11] - sy[1] * sy[2] + sy[10] * sy[3]) * Log(-x, sd) };

    return res;
}
complex<double> G6_abc00d_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ 0, b / c, a / c }, 1),
        G({ a, 0, d }, { sa, 1, sd }, x), G({ 0, 0, d / c }, 1),
        G({ a, b, d }, { sa, sb, sd }, x), G({ b / c, 0, a / c }, 1),
        G({ a, d }, { sa, sd }, x), G({ 0, 0, b / c, a / c }, 1),
        G({ 0, a }, { 1, sa }, x), G({ 0, 0, b / c, d / c }, 1),
        G({ a, b }, { sa, sb }, x), G({ 0, 0, d / c, b / c }, 1),
        G({ a, 0, 0, d }, { sa, 1, 1, sd }, x), G({ 0, b / c, 0, d / c }, 1),
        G({ b / c, 0, 0, d / c }, 1), G({ 0, b / c, 0, a / c }, 1),
        G({ b / c, 0, 0, a / c }, 1), G({ a }, { sa }, x),
        G({ 0, 0, b / c, d / c, a / c }, 1), G({ 0, 0, d / c, b / c, a / c }, 1),
        G({ 0, b / c, 0, d / c, a / c }, 1), G({ b / c, 0, 0, d / c, a / c }, 1),
        G({ 0, 0, b / c, a / c, d / c }, 1), G({ 0, b / c, 0, a / c, d / c }, 1),
        G({ 0, b / c, a / c, 0, d / c }, 1), G({ b / c, 0, 0, a / c, d / c }, 1),
        G({ b / c, 0, a / c, 0, d / c }, 1), G({ b / c, a / c, 0, 0, d / c }, 1) };
    complex<double> res { -(sy[0] * sy[1]) - sy[1] * sy[4] - sy[5] * sy[6]
        + (sy[10] + sy[12] + sy[13]) * sy[7] + sy[7] * sy[8]
        + sy[5] * (sy[12] + sy[13] - sy[14] - sy[15] + sy[8])
        + (-sy[1] + sy[3]) * G({ 0, 0, b / c }, 1)
        + sy[9] * (sy[10] + 3. * G({ 0, 0, 0, d / c }, 1))
        + sy[16]
            * (-sy[17] - sy[18] - sy[19] - sy[20] - 3. * G({ 0, 0, 0, b / c, d / c }, 1)
                - 3. * G({ 0, 0, 0, d / c, b / c }, 1)
                - 3. * G({ 0, 0, b / c, 0, d / c }, 1) - G({ 0, 0, d / c, 0, b / c }, 1)
                - 3. * G({ 0, b / c, 0, 0, d / c }, 1)
                - 3. * G({ b / c, 0, 0, 0, d / c }, 1))
        + 3. * G({ 0, 0, 0, b / c, a / c, d / c }, 1)
        + 3. * G({ 0, 0, 0, b / c, d / c, a / c }, 1)
        + 3. * G({ 0, 0, 0, d / c, b / c, a / c }, 1)
        + 3. * G({ 0, 0, b / c, 0, a / c, d / c }, 1)
        + 3. * G({ 0, 0, b / c, 0, d / c, a / c }, 1)
        + 3. * G({ 0, 0, b / c, a / c, 0, d / c }, 1)
        + G({ 0, 0, b / c, a / c, d / c, x / c }, 1)
        + G({ 0, 0, b / c, d / c, 0, a / c }, 1)
        + G({ 0, 0, b / c, d / c, a / c, x / c }, 1)
        + G({ 0, 0, d / c, 0, b / c, a / c }, 1) + G({ 0, 0, d / c, b / c, 0, a / c }, 1)
        + G({ 0, 0, d / c, b / c, a / c, x / c }, 1)
        + 3. * G({ 0, b / c, 0, 0, a / c, d / c }, 1)
        + 3. * G({ 0, b / c, 0, 0, d / c, a / c }, 1)
        + 3. * G({ 0, b / c, 0, a / c, 0, d / c }, 1)
        + G({ 0, b / c, 0, a / c, d / c, x / c }, 1)
        + G({ 0, b / c, 0, d / c, 0, a / c }, 1)
        + G({ 0, b / c, 0, d / c, a / c, x / c }, 1)
        + 3. * G({ 0, b / c, a / c, 0, 0, d / c }, 1)
        + G({ 0, b / c, a / c, 0, d / c, x / c }, 1)
        + 3. * G({ b / c, 0, 0, 0, a / c, d / c }, 1)
        + 3. * G({ b / c, 0, 0, 0, d / c, a / c }, 1)
        + 3. * G({ b / c, 0, 0, a / c, 0, d / c }, 1)
        + G({ b / c, 0, 0, a / c, d / c, x / c }, 1)
        + G({ b / c, 0, 0, d / c, 0, a / c }, 1)
        + G({ b / c, 0, 0, d / c, a / c, x / c }, 1)
        + 3. * G({ b / c, 0, a / c, 0, 0, d / c }, 1)
        + G({ b / c, 0, a / c, 0, d / c, x / c }, 1)
        + 3. * G({ b / c, a / c, 0, 0, 0, d / c }, 1)
        + G({ b / c, a / c, 0, 0, d / c, x / c }, 1)
        + (sy[14] + sy[15] + sy[6] + G({ 0, b / c, a / c, x / c }, 1)
              + G({ b / c, 0, a / c, x / c }, 1) + G({ b / c, a / c, 0, x / c }, 1))
            * G({ 0, d }, { 1, sd }, x)
        + (sy[0] + sy[4] + G({ b / c, a / c, x / c }, 1))
            * G({ 0, 0, d }, { 1, 1, sd }, x)
        + sy[2]
            * (-sy[3] - G({ 0, a, b }, { 1, sa, sb }, x)
                - G({ a, 0, b }, { sa, 1, sb }, x))
        + G({ b / c, a / c }, 1) * (-sy[11] + G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x))
        + G({ 0, b / c }, 1) * (-sy[11] + G({ a, b, 0, d }, { sa, sb, 1, sd }, x))
        + G({ b / c }, 1)
            * (-G({ a, 0, 0, 0, d }, { sa, 1, 1, 1, sd }, x)
                + G({ a, b, 0, 0, d }, { sa, sb, 1, 1, sd }, x))
        + G({ a, b, 0, 0, 0, d }, { sa, sb, 1, 1, 1, sd }, x)
        + (-sy[17] - sy[18] - sy[19] - sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25]
              - sy[26] - sy[9] * sy[2] + sy[16] * (sy[10] + sy[12] + sy[13] + sy[8]))
            * Log(c, sc)
        + (sy[17] + sy[18] + sy[19] + sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25]
              + sy[26] + sy[9] * sy[2] + sy[16] * (-sy[10] - sy[12] - sy[13] - sy[8]))
            * Log(-x, sc) };
    if (d != x) {
        res += (-sy[21] - sy[22] - sy[23] - sy[24] - sy[25] - sy[26]
                   + G({ 0, 0, b / c, a / c, x / c }, 1)
                   + G({ 0, b / c, 0, a / c, x / c }, 1)
                   + G({ 0, b / c, a / c, 0, x / c }, 1)
                   + G({ b / c, 0, 0, a / c, x / c }, 1)
                   + G({ b / c, 0, a / c, 0, x / c }, 1)
                   + G({ b / c, a / c, 0, 0, x / c }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}
complex<double> G6_abc00d_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (b == c) {
        // abbd
        const vector<complex<double>> sy = { G({ 0, 1, a / b }, 1),
            G({ a, 0, d }, { sa, 1, sd }, x), G({ 0, a / b, 1 }, 1),
            G({ a, d }, { sa, sd }, x), G({ 0, 0, 1, a / b }, 1),
            G({ 0, a }, { 1, sa }, x), G({ 0, 0, 1, d / b }, 1), G({ 0, 0, d / b, 1 }, 1),
            G({ 0, 1, 0, d / b }, 1), G({ 0, 0, a / b, 1 }, 1), G({ 0, 1, 0, a / b }, 1),
            G({ b, 0, 0, d }, { sb, 1, 1, sd }, x), G({ a }, { sa }, x),
            G({ 0, 0, 1, d / b, a / b }, 1), G({ 0, 0, d / b, 1, a / b }, 1),
            G({ 0, 0, d / b, a / b, 1 }, 1), G({ 0, 1, 0, d / b, a / b }, 1),
            G({ 0, 0, 1, a / b, d / b }, 1), G({ 0, 0, a / b, 1, d / b }, 1),
            G({ 0, 0, a / b, d / b, 1 }, 1), G({ 0, 1, 0, a / b, d / b }, 1),
            G({ 0, 1, a / b, 0, d / b }, 1), G({ 0, a / b, 0, 1, d / b }, 1),
            G({ 0, a / b, 0, d / b, 1 }, 1), G({ 0, a / b, 1, 0, d / b }, 1),
            G({ a / b, 0, 0, 1, d / b }, 1), G({ a / b, 0, 0, d / b, 1 }, 1),
            G({ a / b, 0, 1, 0, d / b }, 1) };
        complex<double> res { -(sy[0] * sy[1]) - sy[1] * sy[2] - sy[3] * sy[4]
            + sy[5] * sy[6] + sy[5] * (sy[7] + sy[8])
            + sy[3] * (-sy[9] - sy[10] + sy[6] + sy[7] + sy[8])
            - sy[11] * G({ a / b, x / b }, 1)
            + sy[12]
                * (-sy[13] - sy[14] - sy[15] - sy[16] - 3. * G({ 0, 0, 0, 1, d / b }, 1)
                    - 3. * G({ 0, 0, 0, d / b, 1 }, 1) - 3. * G({ 0, 0, 1, 0, d / b }, 1)
                    - 2. * G({ 0, 1, 0, 0, d / b }, 1))
            + 3. * G({ 0, 0, 0, 1, a / b, d / b }, 1)
            + 3. * G({ 0, 0, 0, 1, d / b, a / b }, 1)
            + 3. * G({ 0, 0, 0, a / b, 1, d / b }, 1)
            + 3. * G({ 0, 0, 0, a / b, d / b, 1 }, 1)
            + 3. * G({ 0, 0, 0, d / b, 1, a / b }, 1)
            + 3. * G({ 0, 0, 0, d / b, a / b, 1 }, 1)
            + 3. * G({ 0, 0, 1, 0, a / b, d / b }, 1)
            + 3. * G({ 0, 0, 1, 0, d / b, a / b }, 1)
            + 3. * G({ 0, 0, 1, a / b, 0, d / b }, 1)
            + G({ 0, 0, 1, a / b, d / b, x / b }, 1) + G({ 0, 0, 1, d / b, 0, a / b }, 1)
            + G({ 0, 0, 1, d / b, a / b, x / b }, 1)
            + 3. * G({ 0, 0, a / b, 0, 1, d / b }, 1)
            + 3. * G({ 0, 0, a / b, 0, d / b, 1 }, 1)
            + 3. * G({ 0, 0, a / b, 1, 0, d / b }, 1)
            + G({ 0, 0, a / b, 1, d / b, x / b }, 1)
            + G({ 0, 0, a / b, d / b, 1, x / b }, 1)
            + G({ 0, 0, a / b, d / b, x / b, 1 }, 1) + G({ 0, 0, d / b, 0, 1, a / b }, 1)
            + G({ 0, 0, d / b, 0, a / b, 1 }, 1) + G({ 0, 0, d / b, 1, 0, a / b }, 1)
            + G({ 0, 0, d / b, 1, a / b, x / b }, 1)
            + G({ 0, 0, d / b, a / b, 1, x / b }, 1)
            + G({ 0, 0, d / b, a / b, x / b, 1 }, 1)
            + 2. * G({ 0, 1, 0, 0, a / b, d / b }, 1)
            + 2. * G({ 0, 1, 0, 0, d / b, a / b }, 1)
            + 2. * G({ 0, 1, 0, a / b, 0, d / b }, 1)
            + G({ 0, 1, 0, a / b, d / b, x / b }, 1) + G({ 0, 1, 0, d / b, 0, a / b }, 1)
            + G({ 0, 1, 0, d / b, a / b, x / b }, 1)
            + 2. * G({ 0, 1, a / b, 0, 0, d / b }, 1)
            + G({ 0, 1, a / b, 0, d / b, x / b }, 1)
            + 3. * G({ 0, a / b, 0, 0, 1, d / b }, 1)
            + 3. * G({ 0, a / b, 0, 0, d / b, 1 }, 1)
            + 3. * G({ 0, a / b, 0, 1, 0, d / b }, 1)
            + G({ 0, a / b, 0, 1, d / b, x / b }, 1)
            + G({ 0, a / b, 0, d / b, 1, x / b }, 1)
            + G({ 0, a / b, 0, d / b, x / b, 1 }, 1)
            + 2. * G({ 0, a / b, 1, 0, 0, d / b }, 1)
            + G({ 0, a / b, 1, 0, d / b, x / b }, 1)
            + 3. * G({ a / b, 0, 0, 0, 1, d / b }, 1)
            + 3. * G({ a / b, 0, 0, 0, d / b, 1 }, 1)
            + 3. * G({ a / b, 0, 0, 1, 0, d / b }, 1)
            + G({ a / b, 0, 0, 1, d / b, x / b }, 1)
            + G({ a / b, 0, 0, d / b, 1, x / b }, 1)
            + G({ a / b, 0, 0, d / b, x / b, 1 }, 1)
            + 2. * G({ a / b, 0, 1, 0, 0, d / b }, 1)
            + G({ a / b, 0, 1, 0, d / b, x / b }, 1)
            + (sy[9] + sy[10] + sy[4] + G({ 0, 1, a / b, x / b }, 1)
                  + G({ 0, a / b, 1, x / b }, 1) + G({ 0, a / b, x / b, 1 }, 1)
                  + G({ a / b, 0, 1, x / b }, 1) + G({ a / b, 0, x / b, 1 }, 1))
                * G({ 0, d }, { 1, sd }, x)
            + (sy[0] + sy[2] + G({ a / b, x / b, 1 }, 1))
                * G({ 0, 0, d }, { 1, 1, sd }, x)
            + G({ a / b, 1 }, 1) * (sy[11] - G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
            + G({ a / b }, 1)
                * (-G({ 0, b, 0, 0, d }, { 1, sb, 1, 1, sd }, x)
                    + G({ a, b, 0, 0, d }, { sa, sb, 1, 1, sd }, x))
            + G({ a, 0, b, 0, 0, d }, { sa, 1, sb, 1, 1, sd }, x)
            + (-sy[13] - sy[14] - sy[15] - sy[16] - sy[17] - sy[18] - sy[19] - sy[20]
                  - sy[21] - sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27]
                  + sy[12] * (sy[6] + sy[7] + sy[8]))
                * Log(b, sb)
            + (sy[13] + sy[14] + sy[15] + sy[16] + sy[17] + sy[18] + sy[19] + sy[20]
                  + sy[21] + sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27]
                  + sy[12] * (-sy[6] - sy[7] - sy[8]))
                * Log(-x, sb) };
        if (d != x) {
            res += (-sy[17] - sy[18] - sy[19] - sy[20] - sy[21] - sy[22] - sy[23] - sy[24]
                       - sy[25] - sy[26] - sy[27] + G({ 0, 0, 1, a / b, x / b }, 1)
                       + G({ 0, 0, a / b, 1, x / b }, 1) + G({ 0, 0, a / b, x / b, 1 }, 1)
                       + G({ 0, 1, 0, a / b, x / b }, 1) + G({ 0, 1, a / b, 0, x / b }, 1)
                       + G({ 0, a / b, 0, 1, x / b }, 1) + G({ 0, a / b, 0, x / b, 1 }, 1)
                       + G({ 0, a / b, 1, 0, x / b }, 1) + G({ a / b, 0, 0, 1, x / b }, 1)
                       + G({ a / b, 0, 0, x / b, 1 }, 1)
                       + G({ a / b, 0, 1, 0, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;

    } else { // abcd
        const vector<complex<double>> sy = { G({ c / b, 0, a / b }, 1),
            G({ a / b, c / b }, 1), G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x),
            G({ c / b, a / b }, 1), G({ a, d }, { sa, sd }, x),
            G({ c / b, 0, 0, a / b }, 1), G({ c / b, 0, 0, d / b }, 1),
            G({ a, c, 0, 0, d }, { sa, sc, 1, 1, sd }, x), G({ a }, { sa }, x),
            G({ c / b, 0, 0, d / b, a / b }, 1), G({ a / b, c / b, 0, 0, d / b }, 1),
            G({ c / b, 0, 0, a / b, d / b }, 1), G({ c / b, 0, a / b, 0, d / b }, 1),
            G({ c / b, a / b, 0, 0, d / b }, 1) };
        complex<double> res {
            -(sy[1] * sy[2]) - sy[2] * sy[3] + sy[4] * sy[5] - sy[4] * sy[6]
            + sy[8]
                * (sy[9] + G({ 0, c / b, 0, 0, d / b }, 1)
                    + 3. * G({ c / b, 0, 0, 0, d / b }, 1))
            - G({ 0, a / b, c / b, 0, 0, d / b }, 1)
            - G({ 0, c / b, 0, 0, a / b, d / b }, 1)
            - G({ 0, c / b, 0, 0, d / b, a / b }, 1)
            - G({ 0, c / b, 0, a / b, 0, d / b }, 1)
            - G({ 0, c / b, a / b, 0, 0, d / b }, 1)
            - G({ a / b, 0, c / b, 0, 0, d / b }, 1)
            - 3. * G({ a / b, c / b, 0, 0, 0, d / b }, 1)
            - G({ a / b, c / b, 0, 0, d / b, x / b }, 1)
            - 3. * G({ c / b, 0, 0, 0, a / b, d / b }, 1)
            - 3. * G({ c / b, 0, 0, 0, d / b, a / b }, 1)
            - 3. * G({ c / b, 0, 0, a / b, 0, d / b }, 1)
            - G({ c / b, 0, 0, a / b, d / b, x / b }, 1)
            - G({ c / b, 0, 0, d / b, 0, a / b }, 1)
            - G({ c / b, 0, 0, d / b, a / b, x / b }, 1)
            - 3. * G({ c / b, 0, a / b, 0, 0, d / b }, 1)
            - G({ c / b, 0, a / b, 0, d / b, x / b }, 1)
            - 3. * G({ c / b, a / b, 0, 0, 0, d / b }, 1)
            - G({ c / b, a / b, 0, 0, d / b, x / b }, 1)
            - sy[6] * G({ 0, a }, { 1, sa }, x)
            + (-sy[5] - G({ a / b, c / b, 0, x / b }, 1)
                  - G({ c / b, 0, a / b, x / b }, 1) - G({ c / b, a / b, 0, x / b }, 1))
                * G({ 0, d }, { 1, sd }, x)
            + (-sy[0] - G({ a / b, c / b, x / b }, 1) - G({ c / b, a / b, x / b }, 1))
                * G({ 0, 0, d }, { 1, 1, sd }, x)
            + sy[0] * G({ a, 0, d }, { sa, 1, sd }, x)
            + sy[3] * G({ a, 0, 0, d }, { sa, 1, 1, sd }, x)
            + G({ a / b }, 1) * (sy[7] - G({ 0, c, 0, 0, d }, { 1, sc, 1, 1, sd }, x))
            + G({ c / b }, 1) * (-sy[7] + G({ a, 0, 0, 0, d }, { sa, 1, 1, 1, sd }, x))
            + G({ a, 0, c, 0, 0, d }, { sa, 1, sc, 1, 1, sd }, x)
            + (sy[9] + sy[10] + sy[11] + sy[12] + sy[13] - sy[6] * sy[8]) * Log(b, sb)
            + (-sy[9] - sy[10] - sy[11] - sy[12] - sy[13] + sy[6] * sy[8]) * Log(-x, sb)
        };
        if (c != x) {
            res += (sy[1] - G({ a / b, x / b }, 1))
                * G({ c, 0, 0, d }, { sc, 1, 1, sd }, x);
        }
        if (d != x) {
            res += (sy[10] + sy[11] + sy[12] + sy[13]
                       - G({ a / b, c / b, 0, 0, x / b }, 1)
                       - G({ c / b, 0, 0, a / b, x / b }, 1)
                       - G({ c / b, 0, a / b, 0, x / b }, 1)
                       - G({ c / b, a / b, 0, 0, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_abc00d_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (a == b) {
        if (a == c) { // aaad
            const vector<complex<double>> sy
                = { G({ 0, 0, 1, 1, d / a }, 1), G({ 0, 0, 1, d / a, 1 }, 1),
                      G({ 0, 0, d / a, 1, 1 }, 1), G({ 0, 1, 0, 1, d / a }, 1),
                      G({ 0, 1, 0, d / a, 1 }, 1), G({ 0, 1, 1, 0, d / a }, 1) };
            complex<double> res { 3. * G({ 0, 0, 0, 1, 1, d / a }, 1)
                + 3. * G({ 0, 0, 0, 1, d / a, 1 }, 1)
                + 3. * G({ 0, 0, 0, d / a, 1, 1 }, 1)
                + 3. * G({ 0, 0, 1, 0, 1, d / a }, 1)
                + 3. * G({ 0, 0, 1, 0, d / a, 1 }, 1)
                + 3. * G({ 0, 0, 1, 1, 0, d / a }, 1) + G({ 0, 0, 1, 1, d / a, x / a }, 1)
                + G({ 0, 0, 1, d / a, 1, x / a }, 1) + G({ 0, 0, 1, d / a, x / a, 1 }, 1)
                + G({ 0, 0, d / a, 1, 1, x / a }, 1) + G({ 0, 0, d / a, 1, x / a, 1 }, 1)
                + G({ 0, 0, d / a, x / a, 1, 1 }, 1) + 2. * G({ 0, 1, 0, 0, 1, d / a }, 1)
                + 2. * G({ 0, 1, 0, 0, d / a, 1 }, 1)
                + 2. * G({ 0, 1, 0, 1, 0, d / a }, 1) + G({ 0, 1, 0, 1, d / a, x / a }, 1)
                + G({ 0, 1, 0, d / a, 1, x / a }, 1) + G({ 0, 1, 0, d / a, x / a, 1 }, 1)
                + 2. * G({ 0, 1, 1, 0, 0, d / a }, 1) + G({ 0, 1, 1, 0, d / a, x / a }, 1)
                + (G({ 0, 1, 1, x / a }, 1) + G({ 0, 1, x / a, 1 }, 1)
                      + G({ 0, x / a, 1, 1 }, 1))
                    * G({ 0, d }, { 1, sd }, x)
                + G({ x / a, 1, 1 }, 1) * G({ 0, 0, d }, { 1, 1, sd }, x)
                - G({ x / a, 1 }, 1) * G({ a, 0, 0, d }, { sa, 1, 1, sd }, x)
                + G({ x / a }, 1) * G({ a, a, 0, 0, d }, { sa, sa, 1, 1, sd }, x)
                + G({ 0, a, a, 0, 0, d }, { 1, sa, sa, 1, 1, sd }, x)
                + (-sy[0] - sy[1] - sy[2] - sy[3] - sy[4] - sy[5]) * Log(a, sa)
                + (sy[0] + sy[1] + sy[2] + sy[3] + sy[4] + sy[5]) * Log(-x, sa) };
            if (d != x) {
                res += (-sy[0] - sy[1] - sy[2] - sy[3] - sy[4] - sy[5]
                           + G({ 0, 0, 1, 1, x / a }, 1) + G({ 0, 0, 1, x / a, 1 }, 1)
                           + G({ 0, 0, x / a, 1, 1 }, 1) + G({ 0, 1, 0, 1, x / a }, 1)
                           + G({ 0, 1, 0, x / a, 1 }, 1) + G({ 0, 1, 1, 0, x / a }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        } else { // aacd

            const vector<complex<double>> sy = { G({ c / a, 1 }, 1),
                G({ c / a, 0, 0, 1, d / a }, 1), G({ c / a, 0, 0, d / a, 1 }, 1),
                G({ c / a, 0, 1, 0, d / a }, 1), G({ c / a, 1, 0, 0, d / a }, 1) };
            complex<double> res { -G({ 0, c / a, 0, 0, 1, d / a }, 1)
                - G({ 0, c / a, 0, 0, d / a, 1 }, 1) - G({ 0, c / a, 0, 1, 0, d / a }, 1)
                - G({ 0, c / a, 1, 0, 0, d / a }, 1)
                - 3. * G({ c / a, 0, 0, 0, 1, d / a }, 1)
                - 3. * G({ c / a, 0, 0, 0, d / a, 1 }, 1)
                - 3. * G({ c / a, 0, 0, 1, 0, d / a }, 1)
                - G({ c / a, 0, 0, 1, d / a, x / a }, 1)
                - G({ c / a, 0, 0, d / a, 1, x / a }, 1)
                - G({ c / a, 0, 0, d / a, x / a, 1 }, 1)
                - 3. * G({ c / a, 0, 1, 0, 0, d / a }, 1)
                - G({ c / a, 0, 1, 0, d / a, x / a }, 1)
                - 3. * G({ c / a, 1, 0, 0, 0, d / a }, 1)
                - G({ c / a, 1, 0, 0, d / a, x / a }, 1)
                + (-G({ c / a, 0, 1, x / a }, 1) - G({ c / a, 0, x / a, 1 }, 1)
                      - G({ c / a, 1, 0, x / a }, 1))
                    * G({ 0, d }, { 1, sd }, x)
                + (-G({ c / a, 1, x / a }, 1) - G({ c / a, x / a, 1 }, 1))
                    * G({ 0, 0, d }, { 1, 1, sd }, x)
                - sy[0] * G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x)
                + G({ x / a }, 1) * G({ a, c, 0, 0, d }, { sa, sc, 1, 1, sd }, x)
                + G({ 0, a, c, 0, 0, d }, { 1, sa, sc, 1, 1, sd }, x)
                + (sy[1] + sy[2] + sy[3] + sy[4]) * Log(a, sa)
                + (-sy[1] - sy[2] - sy[3] - sy[4]) * Log(-x, sa) };
            if (c != x) {
                res += (sy[0] - G({ x / a, 1 }, 1))
                    * G({ c, 0, 0, d }, { sc, 1, 1, sd }, x);
            }
            if (d != x) {
                res += (sy[1] + sy[2] + sy[3] + sy[4] - G({ c / a, 0, 0, 1, x / a }, 1)
                           - G({ c / a, 0, 0, x / a, 1 }, 1)
                           - G({ c / a, 0, 1, 0, x / a }, 1)
                           - G({ c / a, 1, 0, 0, x / a }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        }
    } else { // abcd
        const vector<complex<double>> sy = { G({ b / a, c / a }, 1), G({ b / a }, 1),
            G({ b / a, c / a, 0, 0, d / a }, 1) };
        complex<double> res { G({ 0, b / a, c / a, 0, 0, d / a }, 1)
            + G({ b / a, 0, c / a, 0, 0, d / a }, 1)
            + 3. * G({ b / a, c / a, 0, 0, 0, d / a }, 1)
            + G({ b / a, c / a, 0, 0, d / a, x / a }, 1)
            + G({ b / a, c / a, 0, x / a }, 1) * G({ 0, d }, { 1, sd }, x)
            + G({ b / a, c / a, x / a }, 1) * G({ 0, 0, d }, { 1, 1, sd }, x)
            + sy[0] * G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x)
            + sy[1] * G({ 0, c, 0, 0, d }, { 1, sc, 1, 1, sd }, x)
            + G({ 0, b, c, 0, 0, d }, { 1, sb, sc, 1, 1, sd }, x) - sy[2] * Log(a, sa)
            + sy[2] * Log(-x, sa) };
        if (b != x) {
            res += (-sy[1] + G({ x / a }, 1))
                * G({ b, c, 0, 0, d }, { sb, sc, 1, 1, sd }, x);
        }
        if (c != x) {
            res += (-sy[0] + G({ b / a, x / a }, 1))
                * G({ c, 0, 0, d }, { sc, 1, 1, sd }, x);
        }
        if (d != x) {
            res += (-sy[2] + G({ b / a, c / a, 0, 0, x / a }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G6_0a0bcd_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { G({ 0, 0, a }, { 1, 1, sa }, x),
        G({ 0, c / d, b / d }, 1), G({ c / d, 0, b / d }, 1), G({ c / d, b / d }, 1),
        G({ 0, a, 0, b }, { 1, sa, 1, sb }, x), Log(d, sd), G({ 0, a }, { 1, sa }, x),
        G({ c / d, b / d, 0, a / d }, 1), G({ 0, 0 }, { 1, 1 }, x),
        G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x), G({ c / d }, 1),
        G({ 0, c / d, b / d, 0, a / d }, 1), G({ c / d, 0, b / d, 0, a / d }, 1),
        G({ c / d, b / d, 0, 0, a / d }, 1) };
    complex<double> res { sy[0] * (-2. * sy[1] - 2. * sy[2])
        + sy[5]
            * (-sy[11] - sy[12] - 2. * sy[13] + 2. * sy[0] * sy[3] + sy[10] * sy[4]
                + (-sy[1] - sy[2]) * sy[6])
        - sy[7] * sy[8] - sy[4] * G({ 0, c / d }, 1)
        + G({ 0, 0, c / d, b / d, 0, a / d }, 1) + G({ 0, c / d, 0, b / d, 0, a / d }, 1)
        + 2. * G({ 0, c / d, b / d, 0, 0, a / d }, 1)
        + G({ c / d, 0, 0, b / d, 0, a / d }, 1)
        + 2. * G({ c / d, 0, b / d, 0, 0, a / d }, 1)
        + 3. * G({ c / d, b / d, 0, 0, 0, a / d }, 1)
        - G({ c / d, b / d, 0, a / d, 0, x / d }, 1)
        + sy[9] * (-G({ x / d }, 1) + G({ d }, { sd }, x))
        + sy[3] * (-sy[4] + 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
        + sy[10]
            * (sy[9] + 2. * G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x)
                + 2. * G({ 0, a, 0, 0, b }, { 1, sa, 1, 1, sb }, x))
        - 2. * G({ 0, 0, a, 0, b, c }, { 1, 1, sa, 1, sb, sc }, x)
        - 2. * G({ 0, a, 0, 0, b, c }, { 1, sa, 1, 1, sb, sc }, x)
        - G({ 0, a, 0, b, 0, c }, { 1, sa, 1, sb, 1, sc }, x)
        + (sy[11] + sy[12] + 2. * sy[13] - 2. * sy[0] * sy[3] - sy[10] * sy[4]
              + (sy[1] + sy[2]) * sy[6] + sy[5] * (-(sy[3] * sy[6]) - sy[7]))
            * Log(-x, sd)
        + ((sy[3] * sy[6]) / 2. + sy[7] / 2.) * pow(sy[5], 2.)
        + sy[6]
            * (-sy[7] + G({ 0, 0, c / d, b / d }, 1) + G({ 0, c / d, 0, b / d }, 1)
                + G({ c / d, 0, 0, b / d }, 1) + sy[3] * (-sy[8] - 2. * Zeta(2)))
        - 2. * sy[7] * Zeta(2) };
    return res;
}
complex<double> G6_0a0bcd_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (c == d) { // abcc
        const vector<complex<double>> sy = { G({ b / c, 0, a / c }, 1),
            G({ 0, a, 0, b }, { 1, sa, 1, sb }, x), G({ b / c, 0, a / c, 1 }, 1) };
        complex<double> res { G({ b / c, 0, a / c, 0, 0, 1 }, 1)
            - G({ b / c, 0, a / c, 0, x / c, 1 }, 1)
            + (-G({ b / c, 0, a / c, 0, 1 }, 1) + G({ b / c, 0, a / c, 0, x / c }, 1))
                * G({ c }, { sc }, x)
            + (-sy[2] + G({ b / c, 0, 0, 1 }, 1)) * G({ 0, a }, { 1, sa }, x)
            + sy[2] * G({ 0, c }, { 1, sc }, x) - sy[0] * G({ 0, 0, c }, { 1, 1, sc }, x)
            + (sy[0] - G({ b / c, 0, 1 }, 1)) * G({ 0, a, c }, { 1, sa, sc }, x)
            + G({ b / c, 1 }, 1) * (-sy[1] + G({ 0, a, 0, c }, { 1, sa, 1, sc }, x))
            + G({ b / c }, 1)
                * (-G({ 0, a, 0, 0, c }, { 1, sa, 1, 1, sc }, x)
                    + G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x))
            + G({ 0, a, 0, b, 0, c }, { 1, sa, 1, sb, 1, sc }, x) - sy[1] * Zeta(2) };
        return res;
    } else { // abcd
        const vector<complex<double>> sy = { G({ b / c, 0, a / c }, 1),
            G({ 0, a, d }, { 1, sa, sd }, x), G({ 0, 0, a }, { 1, 1, sa }, x),
            G({ 0, b / c, d / c }, 1), G({ 0, d / c, b / c }, 1),
            G({ d / c, 0, b / c }, 1), G({ 0, a, 0, b }, { 1, sa, 1, sb }, x),
            G({ d / c, b / c }, 1), G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x),
            G({ b / c, d / c }, 1), Log(c, sc), G({ 0, a }, { 1, sa }, x),
            G({ b / c, 0, a / c, d / c }, 1), G({ b / c, 0, d / c, a / c }, 1),
            G({ b / c, d / c, 0, a / c }, 1), G({ d / c, b / c, 0, a / c }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ d / c }, 1),
            G({ 0, a, 0, b, d }, { 1, sa, 1, sb, sd }, x),
            G({ 0, b / c, 0, a / c, d / c }, 1), G({ 0, b / c, 0, d / c, a / c }, 1),
            G({ 0, b / c, d / c, 0, a / c }, 1), G({ 0, d / c, b / c, 0, a / c }, 1),
            G({ b / c, 0, 0, a / c, d / c }, 1), G({ b / c, 0, 0, d / c, a / c }, 1),
            G({ b / c, 0, d / c, 0, a / c }, 1), G({ b / c, d / c, 0, 0, a / c }, 1),
            G({ d / c, 0, b / c, 0, a / c }, 1), G({ d / c, b / c, 0, 0, a / c }, 1) };
        complex<double> res { (sy[12] + sy[13] + sy[14] + sy[15]) * sy[16] + sy[0] * sy[1]
            + sy[2] * (2. * sy[3] + 2. * sy[4] + 2. * sy[5])
            + sy[10]
                * (sy[19] + sy[20] + sy[21] + sy[22] + 2. * sy[23] + 2. * sy[24]
                    + 2. * sy[25] + 2. * sy[26] + sy[27] + 2. * sy[28]
                    - 2. * sy[9] * sy[2] + sy[11] * (sy[3] + sy[4] + sy[5])
                    - sy[17] * sy[6] - 2. * sy[2] * sy[7])
            + sy[7] * (sy[6] - 3. * sy[8]) + sy[6] * G({ 0, d / c }, 1)
            - sy[1] * G({ b / c, 0, d / c }, 1) - G({ 0, 0, b / c, 0, a / c, d / c }, 1)
            - G({ 0, 0, b / c, 0, d / c, a / c }, 1)
            - G({ 0, 0, b / c, d / c, 0, a / c }, 1)
            - G({ 0, 0, d / c, b / c, 0, a / c }, 1)
            - 2. * G({ 0, b / c, 0, 0, a / c, d / c }, 1)
            - 2. * G({ 0, b / c, 0, 0, d / c, a / c }, 1)
            - 2. * G({ 0, b / c, 0, d / c, 0, a / c }, 1)
            - 2. * G({ 0, b / c, d / c, 0, 0, a / c }, 1)
            - G({ 0, d / c, 0, b / c, 0, a / c }, 1)
            - 2. * G({ 0, d / c, b / c, 0, 0, a / c }, 1)
            - 3. * G({ b / c, 0, 0, 0, a / c, d / c }, 1)
            - 3. * G({ b / c, 0, 0, 0, d / c, a / c }, 1)
            - 3. * G({ b / c, 0, 0, d / c, 0, a / c }, 1)
            + G({ b / c, 0, a / c, 0, 0, d / c }, 1)
            + G({ b / c, 0, a / c, 0, d / c, x / c }, 1)
            + G({ b / c, 0, a / c, d / c, 0, x / c }, 1)
            - 3. * G({ b / c, 0, d / c, 0, 0, a / c }, 1)
            + G({ b / c, 0, d / c, a / c, 0, x / c }, 1)
            - 3. * G({ b / c, d / c, 0, 0, 0, a / c }, 1)
            + G({ b / c, d / c, 0, a / c, 0, x / c }, 1)
            - G({ d / c, 0, 0, b / c, 0, a / c }, 1)
            - 2. * G({ d / c, 0, b / c, 0, 0, a / c }, 1)
            - 3. * G({ d / c, b / c, 0, 0, 0, a / c }, 1)
            + G({ d / c, b / c, 0, a / c, 0, x / c }, 1)
            - sy[0] * G({ 0, 0, d }, { 1, 1, sd }, x)
            + sy[9] * (-3. * sy[8] + G({ 0, a, 0, d }, { 1, sa, 1, sd }, x))
            + sy[17]
                * (-sy[18] - 2. * G({ 0, 0, a, 0, b }, { 1, 1, sa, 1, sb }, x)
                    - 2. * G({ 0, a, 0, 0, b }, { 1, sa, 1, 1, sb }, x))
            + G({ b / c }, 1) * (sy[18] - G({ 0, a, 0, 0, d }, { 1, sa, 1, 1, sd }, x))
            + G({ 0, a, 0, b, 0, d }, { 1, sa, 1, sb, 1, sd }, x)
            + (-sy[19] - sy[20] - sy[21] - sy[22] - 2. * sy[23] - 2. * sy[24]
                  - 2. * sy[25] - 2. * sy[26] - sy[27] - 2. * sy[28] + 2. * sy[9] * sy[2]
                  + sy[11] * (-sy[3] - sy[4] - sy[5]) + sy[17] * sy[6]
                  + 2. * sy[2] * sy[7]
                  + sy[10]
                      * (sy[12] + sy[13] + sy[14] + sy[15] + sy[11] * (sy[9] + sy[7])))
                * Log(-x, sc)
            + (-sy[12] / 2. - sy[13] / 2. - sy[14] / 2. - sy[15] / 2.
                  + sy[11] * (-sy[9] / 2. - sy[7] / 2.))
                * pow(sy[10], 2.)
            + 2. * sy[13] * Zeta(2) + 2. * sy[14] * Zeta(2) + 2. * sy[15] * Zeta(2)
            + sy[12] * (G({ 0, d }, { 1, sd }, x) + 2. * Zeta(2))
            + sy[11]
                * (sy[13] + sy[14] + sy[15] + sy[16] * sy[7]
                    - G({ 0, 0, b / c, d / c }, 1) - G({ 0, 0, d / c, b / c }, 1)
                    - G({ 0, d / c, 0, b / c }, 1) + G({ b / c, 0, 0, d / c }, 1)
                    - G({ d / c, 0, 0, b / c }, 1) + 2. * sy[7] * Zeta(2)
                    + sy[9] * (sy[16] + 2. * Zeta(2))) };

        if (d != x) {
            res += (-G({ b / c, 0, a / c, 0, d / c }, 1)
                       + G({ b / c, 0, a / c, 0, x / c }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_0a0bcd_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (b == c) {
        if (b == d) { // abbb
            const vector<complex<double>> sy = { G({ 0, a / b, 1 }, 1),
                G({ 0, a, b, b }, { 1, sa, sb, sb }, x), G({ b, b }, { sb, sb }, x),
                G({ 0, a / b, 0, 1 }, 1), G({ 0, a / b, 1, 1 }, 1) };
            complex<double> res { sy[2] * sy[3] - sy[2] * G({ 0, a / b, 0, x / b }, 1)
                + G({ 0, a / b, 0, 0, 1, 1 }, 1) - G({ 0, a / b, 0, x / b, 1, 1 }, 1)
                + (-G({ 0, a / b, 0, 1, 1 }, 1) + G({ 0, a / b, 0, x / b, 1 }, 1))
                    * G({ b }, { sb }, x)
                + (-sy[3] + sy[4]) * G({ 0, b }, { 1, sb }, x)
                - sy[0] * G({ 0, b, b }, { 1, sb, sb }, x)
                + G({ 0, a / b }, 1) * (-sy[1] + G({ 0, 0, b, b }, { 1, 1, sb, sb }, x))
                + G({ 0, a, 0, 0, b, b }, { 1, sa, 1, 1, sb, sb }, x) - sy[1] * Zeta(2)
                + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x) * Zeta(2)
                + G({ 0, a }, { 1, sa }, x) * (-sy[4] + Zeta(4) / 4.)
                + G({ 0, a, b }, { 1, sa, sb }, x) * (sy[0] - Zeta(3)) };
            return res;
        } else { // abbd
            const vector<complex<double>> sy = { G({ 0, a / b, 1 }, 1),
                G({ d / b, 1 }, 1), G({ 0, a, 0, d }, { 1, sa, 1, sd }, x),
                G({ 0, a, b, d }, { 1, sa, sb, sd }, x), G({ b, d }, { sb, sd }, x),
                G({ 0, a / b, 0, 1 }, 1), G({ 0, a / b, d / b, 1 }, 1), Log(b, sb),
                G({ 0, a }, { 1, sa }, x), G({ 0, d / b, 1, a / b }, 1),
                G({ 0, d / b, a / b, 1 }, 1), G({ d / b, 0, 1, a / b }, 1),
                G({ d / b, 0, a / b, 1 }, 1), G({ d / b, 1, 0, a / b }, 1),
                G({ 0, 0 }, { 1, 1 }, x), G({ 0, 0, a }, { 1, 1, sa }, x),
                G({ 0, 0, a / b, d / b, 1 }, 1), G({ 0, 0, d / b, 1, a / b }, 1),
                G({ 0, 0, d / b, a / b, 1 }, 1), G({ 0, d / b, 0, 1, a / b }, 1),
                G({ 0, d / b, 0, a / b, 1 }, 1), G({ 0, d / b, 1, 0, a / b }, 1),
                G({ d / b, 0, 0, 1, a / b }, 1), G({ d / b, 0, 0, a / b, 1 }, 1),
                G({ d / b, 0, 1, 0, a / b }, 1), G({ d / b, 1, 0, 0, a / b }, 1) };
            complex<double> res { sy[4] * sy[5]
                + sy[14] * (sy[9] + sy[10] + sy[11] + sy[12] + sy[13] + sy[6])
                + (2. * sy[16] + 2. * sy[17] + 2. * sy[18] - 2. * sy[15] * sy[1]
                      + 2. * sy[19] + 2. * sy[20] + 2. * sy[21] + 2. * sy[22]
                      + 2. * sy[23] + 2. * sy[24] + 2. * sy[25])
                    * sy[7]
                - sy[4] * G({ 0, a / b, 0, x / b }, 1)
                - 3. * G({ 0, 0, 0, a / b, d / b, 1 }, 1)
                - 3. * G({ 0, 0, 0, d / b, 1, a / b }, 1)
                - 3. * G({ 0, 0, 0, d / b, a / b, 1 }, 1)
                - 3. * G({ 0, 0, d / b, 0, 1, a / b }, 1)
                - 3. * G({ 0, 0, d / b, 0, a / b, 1 }, 1)
                - 3. * G({ 0, 0, d / b, 1, 0, a / b }, 1)
                + G({ 0, a / b, 0, 0, d / b, 1 }, 1)
                + G({ 0, a / b, 0, d / b, 1, x / b }, 1)
                + G({ 0, a / b, 0, d / b, x / b, 1 }, 1)
                + G({ 0, a / b, d / b, 0, 1, x / b }, 1)
                + G({ 0, a / b, d / b, 0, x / b, 1 }, 1)
                + G({ 0, a / b, d / b, 1, 0, x / b }, 1)
                - 3. * G({ 0, d / b, 0, 0, 1, a / b }, 1)
                - 3. * G({ 0, d / b, 0, 0, a / b, 1 }, 1)
                - 3. * G({ 0, d / b, 0, 1, 0, a / b }, 1)
                - 3. * G({ 0, d / b, 1, 0, 0, a / b }, 1)
                + G({ 0, d / b, 1, a / b, 0, x / b }, 1)
                + G({ 0, d / b, a / b, 0, 1, x / b }, 1)
                + G({ 0, d / b, a / b, 0, x / b, 1 }, 1)
                + G({ 0, d / b, a / b, 1, 0, x / b }, 1)
                - 3. * G({ d / b, 0, 0, 0, 1, a / b }, 1)
                - 3. * G({ d / b, 0, 0, 0, a / b, 1 }, 1)
                - 3. * G({ d / b, 0, 0, 1, 0, a / b }, 1)
                - 3. * G({ d / b, 0, 1, 0, 0, a / b }, 1)
                + G({ d / b, 0, 1, a / b, 0, x / b }, 1)
                + G({ d / b, 0, a / b, 0, 1, x / b }, 1)
                + G({ d / b, 0, a / b, 0, x / b, 1 }, 1)
                + G({ d / b, 0, a / b, 1, 0, x / b }, 1)
                - 3. * G({ d / b, 1, 0, 0, 0, a / b }, 1)
                + G({ d / b, 1, 0, a / b, 0, x / b }, 1)
                + (-sy[5] + sy[6]) * G({ 0, d }, { 1, sd }, x)
                + (sy[0] - G({ 0, d / b, 1 }, 1)) * G({ 0, a, d }, { 1, sa, sd }, x)
                - sy[0] * G({ 0, b, d }, { 1, sb, sd }, x)
                + sy[1] * (sy[2] - 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
                + G({ 0, a / b }, 1) * (-sy[3] + G({ 0, 0, b, d }, { 1, 1, sb, sd }, x))
                + G({ 0, a, 0, 0, b, d }, { 1, sa, 1, 1, sb, sd }, x)
                + (-2. * sy[16] - 2. * sy[17] - 2. * sy[18] + 2. * sy[15] * sy[1]
                      - 2. * sy[19] - 2. * sy[20] - 2. * sy[21] - 2. * sy[22]
                      - 2. * sy[23] - 2. * sy[24] - 2. * sy[25]
                      + sy[7]
                          * (sy[9] + sy[10] + sy[11] + sy[12] + sy[13] + sy[6]
                              + sy[1] * sy[8]))
                    * Log(-x, sb)
                + (-sy[9] / 2. - sy[10] / 2. - sy[11] / 2. - sy[12] / 2. - sy[13] / 2.
                      - sy[6] / 2. - (sy[1] * sy[8]) / 2.)
                    * pow(sy[7], 2.)
                + 2. * sy[9] * Zeta(2) + 2. * sy[10] * Zeta(2) + 2. * sy[11] * Zeta(2)
                + 2. * sy[12] * Zeta(2) + 2. * sy[13] * Zeta(2) + sy[2] * Zeta(2)
                - sy[3] * Zeta(2) + 2. * sy[6] * Zeta(2)
                + sy[8]
                    * (sy[9] + sy[10] + sy[11] + sy[12] + sy[13]
                        + G({ 0, 0, d / b, 1 }, 1) + sy[1] * (sy[14] + 2. * Zeta(2))) };
            if (d != x) {
                res += (-G({ 0, a / b, 0, d / b, 1 }, 1)
                           + G({ 0, a / b, 0, x / b, 1 }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        }
    } else { // abcd
        const vector<complex<double>> sy = { G({ 0, a / b, c / b }, 1),
            G({ 0, a, d }, { 1, sa, sd }, x), G({ 0, c / b, a / b }, 1),
            G({ c / b, 0, a / b }, 1), G({ c / b, d / b }, 1),
            G({ 0, a, 0, d }, { 1, sa, 1, sd }, x),
            G({ 0, a, c, d }, { 1, sa, sc, sd }, x), G({ 0, a / b, 0, c / b }, 1),
            G({ 0, a / b, c / b, d / b }, 1), G({ 0, c / b, a / b, d / b }, 1),
            G({ c / b, 0, a / b, d / b }, 1), G({ 0, 0 }, { 1, 1 }, x),
            G({ 0, c / b, d / b, a / b }, 1), G({ c / b, 0, d / b, a / b }, 1),
            G({ c / b, d / b, 0, a / b }, 1), Log(b, sb), G({ 0, a }, { 1, sa }, x),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ 0, 0, a / b, c / b, d / b }, 1),
            G({ 0, 0, c / b, a / b, d / b }, 1), G({ 0, 0, c / b, d / b, a / b }, 1),
            G({ 0, c / b, 0, a / b, d / b }, 1), G({ 0, c / b, 0, d / b, a / b }, 1),
            G({ 0, c / b, d / b, 0, a / b }, 1), G({ c / b, 0, 0, a / b, d / b }, 1),
            G({ c / b, 0, 0, d / b, a / b }, 1), G({ c / b, 0, d / b, 0, a / b }, 1),
            G({ c / b, d / b, 0, 0, a / b }, 1) };
        complex<double> res { -(sy[1] * sy[2])
            + sy[15]
                * (-2. * sy[18] - 2. * sy[19] - 2. * sy[20] - 2. * sy[21] - 2. * sy[22]
                    - 2. * sy[23] - 2. * sy[24] - 2. * sy[25] - 2. * sy[26] - 2. * sy[27]
                    + 2. * sy[17] * sy[4])
            + sy[11] * (-sy[9] - sy[10] - sy[12] - sy[13] - sy[14] - sy[8])
            + (-sy[5] + sy[6]) * G({ 0, c / b }, 1)
            + sy[1] * (-sy[3] + G({ 0, c / b, d / b }, 1) + G({ c / b, 0, d / b }, 1))
            + 3. * G({ 0, 0, 0, a / b, c / b, d / b }, 1)
            + 3. * G({ 0, 0, 0, c / b, a / b, d / b }, 1)
            + 3. * G({ 0, 0, 0, c / b, d / b, a / b }, 1)
            + 3. * G({ 0, 0, c / b, 0, a / b, d / b }, 1)
            + 3. * G({ 0, 0, c / b, 0, d / b, a / b }, 1)
            + 3. * G({ 0, 0, c / b, d / b, 0, a / b }, 1)
            - G({ 0, a / b, 0, 0, c / b, d / b }, 1)
            - G({ 0, a / b, 0, c / b, 0, d / b }, 1)
            - G({ 0, a / b, 0, c / b, d / b, x / b }, 1)
            - G({ 0, a / b, c / b, 0, 0, d / b }, 1)
            - G({ 0, a / b, c / b, 0, d / b, x / b }, 1)
            - G({ 0, a / b, c / b, d / b, 0, x / b }, 1)
            + 3. * G({ 0, c / b, 0, 0, a / b, d / b }, 1)
            + 3. * G({ 0, c / b, 0, 0, d / b, a / b }, 1)
            + 3. * G({ 0, c / b, 0, d / b, 0, a / b }, 1)
            - G({ 0, c / b, a / b, 0, 0, d / b }, 1)
            - G({ 0, c / b, a / b, 0, d / b, x / b }, 1)
            - G({ 0, c / b, a / b, d / b, 0, x / b }, 1)
            + 3. * G({ 0, c / b, d / b, 0, 0, a / b }, 1)
            - G({ 0, c / b, d / b, a / b, 0, x / b }, 1)
            + 3. * G({ c / b, 0, 0, 0, a / b, d / b }, 1)
            + 3. * G({ c / b, 0, 0, 0, d / b, a / b }, 1)
            + 3. * G({ c / b, 0, 0, d / b, 0, a / b }, 1)
            - G({ c / b, 0, a / b, 0, 0, d / b }, 1)
            - G({ c / b, 0, a / b, 0, d / b, x / b }, 1)
            - G({ c / b, 0, a / b, d / b, 0, x / b }, 1)
            + 3. * G({ c / b, 0, d / b, 0, 0, a / b }, 1)
            - G({ c / b, 0, d / b, a / b, 0, x / b }, 1)
            + 3. * G({ c / b, d / b, 0, 0, 0, a / b }, 1)
            - G({ c / b, d / b, 0, a / b, 0, x / b }, 1)
            + (-sy[9] - sy[10] - sy[7] - sy[8]) * G({ 0, d }, { 1, sd }, x)
            + (sy[0] + sy[2] + sy[3]) * G({ 0, 0, d }, { 1, 1, sd }, x)
            - sy[0] * G({ 0, c, d }, { 1, sc, sd }, x)
            + sy[4] * (-sy[5] + 3. * G({ 0, 0, 0, a }, { 1, 1, 1, sa }, x))
            + G({ 0, a / b }, 1) * (-sy[6] + G({ 0, 0, c, d }, { 1, 1, sc, sd }, x))
            + G({ c / b }, 1)
                * (G({ 0, a, 0, 0, d }, { 1, sa, 1, 1, sd }, x)
                    - G({ 0, a, 0, c, d }, { 1, sa, 1, sc, sd }, x))
            + G({ 0, a, 0, 0, c, d }, { 1, sa, 1, 1, sc, sd }, x)
            + (2. * sy[18] + 2. * sy[19] + 2. * sy[20] + 2. * sy[21] + 2. * sy[22]
                  + 2. * sy[23] + 2. * sy[24] + 2. * sy[25] + 2. * sy[26] + 2. * sy[27]
                  - 2. * sy[17] * sy[4]
                  + sy[15]
                      * (-sy[9] - sy[10] - sy[12] - sy[13] - sy[14] - sy[16] * sy[4]
                          - sy[8]))
                * Log(-x, sb)
            + (sy[9] / 2. + sy[10] / 2. + sy[12] / 2. + sy[13] / 2. + sy[14] / 2.
                  + (sy[16] * sy[4]) / 2. + sy[8] / 2.)
                * pow(sy[15], 2.)
            + sy[16]
                * (-sy[12] - sy[13] - sy[14] - G({ 0, 0, c / b, d / b }, 1)
                    - G({ 0, c / b, 0, d / b }, 1) - G({ c / b, 0, 0, d / b }, 1)
                    + sy[4] * (-sy[11] - 2. * Zeta(2)))
            - 2. * sy[9] * Zeta(2) - 2. * sy[10] * Zeta(2) - 2. * sy[12] * Zeta(2)
            - 2. * sy[13] * Zeta(2) - 2. * sy[14] * Zeta(2) - 2. * sy[8] * Zeta(2) };
        if (c != x) {
            res += (sy[7] - G({ 0, a / b, 0, x / b }, 1)) * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (G({ 0, a / b, 0, c / b, d / b }, 1)
                       - G({ 0, a / b, 0, c / b, x / b }, 1)
                       + G({ 0, a / b, c / b, 0, d / b }, 1)
                       - G({ 0, a / b, c / b, 0, x / b }, 1)
                       + G({ 0, c / b, a / b, 0, d / b }, 1)
                       - G({ 0, c / b, a / b, 0, x / b }, 1)
                       + G({ c / b, 0, a / b, 0, d / b }, 1)
                       - G({ c / b, 0, a / b, 0, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_0a0bcd_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ 0, b / a, c / a }, 1), G({ 0, 0, b / a }, 1),
        G({ 0, c, d }, { 1, sc, sd }, x), G({ 0, b, c, d }, { 1, sb, sc, sd }, x),
        G({ 0, 0, b / a, c / a }, 1), G({ 0, b / a, 0, c / a }, 1),
        G({ 0, b / a, c / a, d / a }, 1), Log(a, sa) };
    complex<double> res { -(sy[0] * sy[2]) - 2. * sy[1] * sy[2]
        - sy[3] * G({ 0, x / a }, 1) - 3. * G({ 0, 0, 0, b / a, c / a, d / a }, 1)
        - 2. * G({ 0, 0, b / a, 0, c / a, d / a }, 1)
        - 2. * G({ 0, 0, b / a, c / a, 0, d / a }, 1)
        - 2. * G({ 0, 0, b / a, c / a, d / a, x / a }, 1)
        - G({ 0, b / a, 0, 0, c / a, d / a }, 1) - G({ 0, b / a, 0, c / a, 0, d / a }, 1)
        - G({ 0, b / a, 0, c / a, d / a, x / a }, 1)
        - G({ 0, b / a, c / a, 0, 0, d / a }, 1)
        - G({ 0, b / a, c / a, 0, d / a, x / a }, 1)
        - G({ 0, b / a, c / a, d / a, 0, x / a }, 1)
        + (-2. * sy[4] - sy[5] - sy[6]) * G({ 0, d }, { 1, sd }, x)
        + sy[0] * G({ 0, 0, d }, { 1, 1, sd }, x)
        + G({ 0, b / a }, 1) * (-sy[3] + G({ 0, 0, c, d }, { 1, 1, sc, sd }, x))
        + G({ 0, 0, 0, b, c, d }, { 1, 1, 1, sb, sc, sd }, x)
        - sy[6] * sy[7] * Log(-x, sa) + (sy[6] * pow(sy[7], 2.)) / 2.
        + sy[6] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
    if (b != x) {
        res += (2. * sy[1] - 2. * G({ 0, 0, x / a }, 1))
            * G({ b, c, d }, { sb, sc, sd }, x);
    }
    if (c != x) {
        res += (2. * sy[4] + sy[5] - 2. * G({ 0, 0, b / a, x / a }, 1)
                   - G({ 0, b / a, 0, x / a }, 1))
            * G({ c, d }, { sc, sd }, x);
    }
    if (d != x) {
        res += (2. * G({ 0, 0, b / a, c / a, d / a }, 1)
                   - 2. * G({ 0, 0, b / a, c / a, x / a }, 1)
                   + G({ 0, b / a, 0, c / a, d / a }, 1)
                   - G({ 0, b / a, 0, c / a, x / a }, 1)
                   + G({ 0, b / a, c / a, 0, d / a }, 1)
                   - G({ 0, b / a, c / a, 0, x / a }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}

complex<double> G6_0ab0cd_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { G({ 0, b, a }, { 1, sb, sa }, x),
        G({ c / d, 0, b / d }, 1), G({ c / d, 0, x / d }, 1),
        G({ 0, a, b }, { 1, sa, sb }, x), G({ 0, c / d }, 1),
        G({ 0, 0, a, b }, { 1, 1, sa, sb }, x), G({ 0, a, 0, b }, { 1, sa, 1, sb }, x),
        G({ c / d, x / d }, 1), G({ 0, a }, { 1, sa }, x), G({ c / d }, 1), Log(d, sd),
        G({ c / d, 0, b / d, a / d }, 1), G({ 0, 0 }, { 1, 1 }, x),
        G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x),
        G({ 0, c / d, 0, b / d, a / d }, 1), G({ c / d, 0, 0, b / d, a / d }, 1),
        G({ c / d, 0, b / d, 0, a / d }, 1) };
    complex<double> res { -(sy[11] * sy[12]) + sy[0] * sy[2]
        + sy[10] * (-sy[14] - 2. * sy[15] - sy[16] + sy[3] * sy[4])
        + sy[4] * (2. * sy[5] + sy[6]) + sy[3] * (sy[2] - G({ 0, 0, c / d }, 1))
        + G({ 0, 0, c / d, 0, b / d, a / d }, 1)
        + 2. * G({ 0, c / d, 0, 0, b / d, a / d }, 1)
        + G({ 0, c / d, 0, b / d, 0, a / d }, 1)
        + 3. * G({ c / d, 0, 0, 0, b / d, a / d }, 1)
        + 2. * G({ c / d, 0, 0, b / d, 0, a / d }, 1)
        + G({ c / d, 0, b / d, 0, 0, a / d }, 1)
        - G({ c / d, 0, b / d, a / d, 0, x / d }, 1)
        + sy[13] * (-G({ x / d }, 1) + G({ d }, { sd }, x))
        + sy[1] * (-sy[0] + 2. * G({ 0, 0, a }, { 1, 1, sa }, x))
        + sy[8]
            * (-sy[11] + sy[10] * sy[1] - G({ 0, c / d, 0, b / d }, 1)
                - 2. * G({ c / d, 0, 0, b / d }, 1) - sy[7] * G({ 0, b }, { 1, sb }, x)
                - sy[9] * G({ 0, 0, b }, { 1, 1, sb }, x))
        + sy[7]
            * (2. * sy[5] + sy[6] + 2. * G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                + G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
        - 2. * G({ 0, 0, a, b, 0, c }, { 1, 1, sa, sb, 1, sc }, x)
        - G({ 0, a, 0, b, 0, c }, { 1, sa, 1, sb, 1, sc }, x)
        - 2. * G({ 0, a, b, 0, 0, c }, { 1, sa, sb, 1, 1, sc }, x)
        + (-(sy[10] * sy[11]) + sy[14] + 2. * sy[15] + sy[16] - sy[3] * sy[4]
              + sy[9] * (sy[10] * sy[3] + 2. * sy[5] + sy[6]) - sy[1] * sy[8])
            * Log(-x, sd)
        + (sy[11] * pow(sy[10], 2.)) / 2. - 2. * sy[11] * Zeta(2)
        + sy[9]
            * (sy[13] + sy[10] * (-2. * sy[5] - sy[6])
                + 3. * G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x)
                + G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x)
                - (sy[3] * pow(sy[10], 2.)) / 2. + sy[3] * (sy[12] + 2. * Zeta(2))) };
    if (b != x) {
        res += (sy[1] * sy[8] - sy[2] * sy[8]) * G({ b }, { sb }, x)
            + (-sy[1] + sy[2]) * G({ b, 0, a }, { sb, 1, sa }, x);
    }
    return res;
}
complex<double> G6_0ab0cd_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (c == d) { // abcc
        const vector<complex<double>> sy = { G({ 0, 1, b / c }, 1),
            G({ 0, a, b }, { 1, sa, sb }, x), G({ 0, b, a }, { 1, sb, sa }, x),
            G({ 0, 1, x / c }, 1), G({ 0, b / c, a / c }, 1), G({ 0, b / c, 1 }, 1),
            G({ 0, x / c, 1 }, 1), G({ 0, a, b, c }, { 1, sa, sb, sc }, x),
            G({ x / c, 1 }, 1), G({ 0, a }, { 1, sa }, x), G({ 0, b / c, a / c, 1 }, 1) };
        complex<double> res { sy[0] * (-sy[1] - sy[2]) + sy[2] * sy[3]
            + sy[2] * (-sy[5] + sy[6]) + G({ 0, b / c, a / c, 0, 0, 1 }, 1)
            - G({ 0, b / c, a / c, 0, x / c, 1 }, 1)
            + (-G({ 0, b / c, a / c, 0, 1 }, 1) + G({ 0, b / c, a / c, 0, x / c }, 1))
                * G({ c }, { sc }, x)
            + sy[9]
                * (-sy[10] + G({ 0, b / c, 0, 1 }, 1) - sy[8] * G({ 0, b }, { 1, sb }, x))
            + sy[10] * G({ 0, c }, { 1, sc }, x) - sy[4] * G({ 0, 0, c }, { 1, 1, sc }, x)
            + (sy[4] - sy[5]) * G({ 0, a, c }, { 1, sa, sc }, x)
            + G({ 0, b / c }, 1) * (-sy[7] + G({ 0, a, 0, c }, { 1, sa, 1, sc }, x))
            + sy[8]
                * (2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                    + 2. * G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                    + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                    + G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
            + G({ 0, a, b, 0, 0, c }, { 1, sa, sb, 1, 1, sc }, x) - sy[7] * Zeta(2)
            + sy[1] * (sy[3] + sy[6] + Zeta(3)) };

        if (b != x) {
            res += (sy[0] * sy[9] - sy[9] * sy[3] + sy[9] * sy[5] - sy[9] * sy[6])
                    * G({ b }, { sb }, x)
                + (-sy[0] + sy[3] - sy[5] + sy[6]) * G({ b, 0, a }, { sb, 1, sa }, x);
        }
        return res;
    } else { // abcd
        const vector<complex<double>> sy = { G({ 0, b / c, a / c }, 1),
            G({ 0, a, d }, { 1, sa, sd }, x), G({ 0, 0, a }, { 1, 1, sa }, x),
            G({ 0, b / c, d / c }, 1), G({ 0, d / c, b / c }, 1),
            G({ d / c, 0, b / c }, 1), G({ 0, b, a }, { 1, sb, sa }, x),
            G({ 0, a, b }, { 1, sa, sb }, x), G({ d / c, 0, x / c }, 1),
            G({ 0, a, b, d }, { 1, sa, sb, sd }, x), G({ d / c, x / c }, 1),
            G({ 0, 0, a, b }, { 1, 1, sa, sb }, x),
            G({ 0, a, 0, b }, { 1, sa, 1, sb }, x), Log(c, sc),
            G({ 0, b / c, a / c, d / c }, 1), G({ 0, b / c, d / c, a / c }, 1),
            G({ 0, d / c, b / c, a / c }, 1), G({ d / c, 0, b / c, a / c }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0, a }, { 1, sa }, x), G({ d / c }, 1),
            G({ 0, 0, b / c, a / c, d / c }, 1), G({ 0, 0, b / c, d / c, a / c }, 1),
            G({ 0, 0, d / c, b / c, a / c }, 1), G({ 0, b / c, 0, a / c, d / c }, 1),
            G({ 0, b / c, 0, d / c, a / c }, 1), G({ 0, b / c, d / c, 0, a / c }, 1),
            G({ 0, d / c, 0, b / c, a / c }, 1), G({ 0, d / c, b / c, 0, a / c }, 1),
            G({ d / c, 0, 0, b / c, a / c }, 1), G({ d / c, 0, b / c, 0, a / c }, 1) };
        complex<double> res { (sy[14] + sy[15] + sy[16] + sy[17]) * sy[18] + sy[0] * sy[1]
            + sy[13]
                * (2. * sy[21] + 2. * sy[22] + 2. * sy[23] + sy[24] + sy[25] + sy[26]
                    + 2. * sy[27] + sy[28] + 2. * sy[29] + sy[30])
            - sy[1] * sy[3] - 2. * sy[2] * sy[3] + sy[2] * (-2. * sy[4] - 2. * sy[5])
            + sy[5] * sy[6] - sy[6] * sy[8] + sy[9] * G({ 0, d / c }, 1)
            + sy[7] * (-sy[4] - sy[8] - G({ 0, 0, d / c }, 1))
            - 3. * G({ 0, 0, 0, b / c, a / c, d / c }, 1)
            - 3. * G({ 0, 0, 0, b / c, d / c, a / c }, 1)
            - 3. * G({ 0, 0, 0, d / c, b / c, a / c }, 1)
            - 2. * G({ 0, 0, b / c, 0, a / c, d / c }, 1)
            - 2. * G({ 0, 0, b / c, 0, d / c, a / c }, 1)
            - 2. * G({ 0, 0, b / c, d / c, 0, a / c }, 1)
            - 3. * G({ 0, 0, d / c, 0, b / c, a / c }, 1)
            - 2. * G({ 0, 0, d / c, b / c, 0, a / c }, 1)
            - G({ 0, b / c, 0, 0, a / c, d / c }, 1)
            - G({ 0, b / c, 0, 0, d / c, a / c }, 1)
            - G({ 0, b / c, 0, d / c, 0, a / c }, 1)
            + G({ 0, b / c, a / c, 0, 0, d / c }, 1)
            + G({ 0, b / c, a / c, 0, d / c, x / c }, 1)
            + G({ 0, b / c, a / c, d / c, 0, x / c }, 1)
            - G({ 0, b / c, d / c, 0, 0, a / c }, 1)
            + G({ 0, b / c, d / c, a / c, 0, x / c }, 1)
            - 3. * G({ 0, d / c, 0, 0, b / c, a / c }, 1)
            - 2. * G({ 0, d / c, 0, b / c, 0, a / c }, 1)
            - G({ 0, d / c, b / c, 0, 0, a / c }, 1)
            + G({ 0, d / c, b / c, a / c, 0, x / c }, 1)
            - 3. * G({ d / c, 0, 0, 0, b / c, a / c }, 1)
            - 2. * G({ d / c, 0, 0, b / c, 0, a / c }, 1)
            - G({ d / c, 0, b / c, 0, 0, a / c }, 1)
            + G({ d / c, 0, b / c, a / c, 0, x / c }, 1)
            + sy[19]
                * (sy[15] + sy[16] + sy[17] + sy[13] * (-sy[3] - sy[4] - sy[5])
                    + 2. * G({ 0, 0, b / c, d / c }, 1)
                    + 2. * G({ 0, 0, d / c, b / c }, 1) + G({ 0, b / c, 0, d / c }, 1)
                    + 2. * G({ 0, d / c, 0, b / c }, 1)
                    + 2. * G({ d / c, 0, 0, b / c }, 1)
                    + sy[10] * G({ 0, b }, { 1, sb }, x)
                    + sy[20] * G({ 0, 0, b }, { 1, 1, sb }, x))
            - sy[0] * G({ 0, 0, d }, { 1, 1, sd }, x)
            + G({ 0, b / c }, 1) * (-sy[9] + G({ 0, a, 0, d }, { 1, sa, 1, sd }, x))
            + sy[10]
                * (-2. * sy[11] - sy[12] - 2. * G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                    - G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
            + G({ 0, a, b, 0, 0, d }, { 1, sa, sb, 1, 1, sd }, x)
            + (sy[13] * (sy[14] + sy[15] + sy[16] + sy[17]) - 2. * sy[21] - 2. * sy[22]
                  - 2. * sy[23] - sy[24] - sy[25] - sy[26] - 2. * sy[27] - sy[28]
                  - 2. * sy[29] - sy[30] + sy[19] * (sy[3] + sy[4] + sy[5])
                  + sy[20] * (-2. * sy[11] - sy[12] - sy[13] * sy[7]))
                * Log(-x, sc)
            + (-sy[14] / 2. - sy[15] / 2. - sy[16] / 2. - sy[17] / 2.) * pow(sy[13], 2.)
            + sy[20]
                * ((2. * sy[11] + sy[12]) * sy[13]
                    - 3. * G({ 0, 0, 0, b, a }, { 1, 1, 1, sb, sa }, x)
                    - G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x)
                    - G({ 0, a, b, 0, d }, { 1, sa, sb, 1, sd }, x)
                    + (sy[7] * pow(sy[13], 2.)) / 2. + sy[7] * (-sy[18] - 2. * Zeta(2)))
            + 2. * sy[15] * Zeta(2) + 2. * sy[16] * Zeta(2) + 2. * sy[17] * Zeta(2)
            + sy[14] * (G({ 0, d }, { 1, sd }, x) + 2. * Zeta(2)) };
        if (b != x) {
            res += (-(sy[19] * sy[5]) + sy[19] * sy[8]) * G({ b }, { sb }, x)
                + (sy[5] - sy[8]) * G({ b, 0, a }, { sb, 1, sa }, x);
        }
        if (d != x) {
            res += (-G({ 0, b / c, a / c, 0, d / c }, 1)
                       + G({ 0, b / c, a / c, 0, x / c }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_0ab0cd_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ 0, a / b, c / b }, 1),
        G({ 0, c, d }, { 1, sc, sd }, x), G({ 0, a, d }, { 1, sa, sd }, x),
        G({ 0, c / b, a / b }, 1), G({ 0, c / b, d / b }, 1), G({ a / b, 0, c / b }, 1),
        G({ 0, a, c, d }, { 1, sa, sc, sd }, x), G({ 0, a }, { 1, sa }, x),
        G({ 0, c / b, d / b, a / b }, 1), Log(b, sb), G({ 0, a / b, c / b, d / b }, 1),
        G({ 0, c / b, a / b, d / b }, 1), G({ a / b, 0, c / b, d / b }, 1),
        G({ 0, 0 }, { 1, 1 }, x), G({ 0, a / b, 0, c / b }, 1),
        G({ a / b, 0, 0, c / b }, 1), G({ 0, 0, a / b, c / b, d / b }, 1),
        G({ 0, 0, c / b, a / b, d / b }, 1), G({ 0, 0, c / b, d / b, a / b }, 1),
        G({ 0, a / b, 0, c / b, d / b }, 1), G({ 0, c / b, 0, a / b, d / b }, 1),
        G({ 0, c / b, 0, d / b, a / b }, 1), G({ 0, c / b, d / b, 0, a / b }, 1) };
    complex<double> res { sy[0] * sy[1] + sy[2] * sy[3] - sy[2] * sy[4]
        + sy[9]
            * (2. * sy[16] + 2. * sy[17] + 2. * sy[18] + sy[19] + sy[20] + sy[21] + sy[22]
                - sy[4] * sy[7])
        + sy[13] * (sy[11] + sy[12] + sy[8]) + sy[1] * (sy[5] + G({ a / b, 0, x / b }, 1))
        + sy[7]
            * (sy[8] + 2. * G({ 0, 0, c / b, d / b }, 1) + G({ 0, c / b, 0, d / b }, 1))
        - 3. * G({ 0, 0, 0, a / b, c / b, d / b }, 1)
        - 3. * G({ 0, 0, 0, c / b, a / b, d / b }, 1)
        - 3. * G({ 0, 0, 0, c / b, d / b, a / b }, 1)
        - G({ 0, 0, a / b, 0, c / b, d / b }, 1)
        - 2. * G({ 0, 0, c / b, 0, a / b, d / b }, 1)
        - 2. * G({ 0, 0, c / b, 0, d / b, a / b }, 1)
        - 2. * G({ 0, 0, c / b, d / b, 0, a / b }, 1)
        + G({ 0, a / b, 0, 0, c / b, d / b }, 1) + G({ 0, a / b, 0, c / b, 0, d / b }, 1)
        + G({ 0, a / b, 0, c / b, d / b, x / b }, 1)
        + G({ 0, a / b, c / b, 0, 0, d / b }, 1)
        + G({ 0, a / b, c / b, 0, d / b, x / b }, 1)
        + G({ 0, a / b, c / b, d / b, 0, x / b }, 1)
        - G({ 0, c / b, 0, 0, a / b, d / b }, 1) - G({ 0, c / b, 0, 0, d / b, a / b }, 1)
        - G({ 0, c / b, 0, d / b, 0, a / b }, 1) + G({ 0, c / b, a / b, 0, 0, d / b }, 1)
        + G({ 0, c / b, a / b, 0, d / b, x / b }, 1)
        + G({ 0, c / b, a / b, d / b, 0, x / b }, 1)
        - G({ 0, c / b, d / b, 0, 0, a / b }, 1)
        + G({ 0, c / b, d / b, a / b, 0, x / b }, 1)
        + 3. * G({ a / b, 0, 0, 0, c / b, d / b }, 1)
        + 2. * G({ a / b, 0, 0, c / b, 0, d / b }, 1)
        + 2. * G({ a / b, 0, 0, c / b, d / b, x / b }, 1)
        + G({ a / b, 0, c / b, 0, 0, d / b }, 1)
        + G({ a / b, 0, c / b, 0, d / b, x / b }, 1)
        + G({ a / b, 0, c / b, d / b, 0, x / b }, 1)
        + (sy[10] + sy[11] + sy[12] + sy[14] + 2. * sy[15]) * G({ 0, d }, { 1, sd }, x)
        - 2. * sy[4] * G({ 0, 0, a }, { 1, 1, sa }, x)
        + (-sy[0] - sy[3] - sy[5]) * G({ 0, 0, d }, { 1, 1, sd }, x)
        + G({ 0, a / b }, 1) * (sy[6] - G({ 0, 0, c, d }, { 1, 1, sc, sd }, x))
        + G({ 0, c / b }, 1) * (-sy[6] + G({ 0, a, 0, d }, { 1, sa, 1, sd }, x))
        + G({ a / b }, 1)
            * (-G({ 0, 0, 0, c, d }, { 1, 1, 1, sc, sd }, x)
                + G({ 0, a, 0, c, d }, { 1, sa, 1, sc, sd }, x))
        + G({ 0, a, 0, 0, c, d }, { 1, sa, 1, 1, sc, sd }, x)
        + (-2. * sy[16] - 2. * sy[17] - 2. * sy[18] - sy[19] - sy[20] - sy[21] - sy[22]
              + sy[4] * sy[7] + sy[9] * (sy[10] + sy[11] + sy[12] + sy[8]))
            * Log(-x, sb)
        + (-sy[10] / 2. - sy[11] / 2. - sy[12] / 2. - sy[8] / 2.) * pow(sy[9], 2.)
        + 2. * sy[11] * Zeta(2) + 2. * sy[12] * Zeta(2) + 2. * sy[8] * Zeta(2)
        + sy[10] * (sy[13] + 2. * Zeta(2)) };
    if (c != x) {
        res += (-sy[14] - 2. * sy[15] + G({ 0, a / b, 0, x / b }, 1)
                   + 2. * G({ a / b, 0, 0, x / b }, 1))
            * G({ c, d }, { sc, sd }, x);
    }
    if (d != x) {
        res += (-sy[19] + G({ 0, a / b, 0, c / b, x / b }, 1)
                   - G({ 0, a / b, c / b, 0, d / b }, 1)
                   + G({ 0, a / b, c / b, 0, x / b }, 1)
                   - G({ 0, c / b, a / b, 0, d / b }, 1)
                   + G({ 0, c / b, a / b, 0, x / b }, 1)
                   - 2. * G({ a / b, 0, 0, c / b, d / b }, 1)
                   + 2. * G({ a / b, 0, 0, c / b, x / b }, 1)
                   - G({ a / b, 0, c / b, 0, d / b }, 1)
                   + G({ a / b, 0, c / b, 0, x / b }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}
complex<double> G6_0ab0cd_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (a == b) { // aacd

        const vector<complex<double>> sy = { G({ 0, 1, c / a }, 1),
            G({ 0, c, d }, { 1, sc, sd }, x), G({ 0, c / a, 1 }, 1), Log(a, sa),
            G({ 0, 1, c / a, d / a }, 1), G({ 0, c / a, 1, d / a }, 1),
            G({ 0, c / a, d / a, 1 }, 1), G({ 0, 0 }, { 1, 1 }, x),
            G({ 0, 0, 1, c / a }, 1), G({ 0, 0, c / a, 1 }, 1), G({ 0, 1, 0, c / a }, 1),
            G({ a, 0, c, d }, { sa, 1, sc, sd }, x) };
        complex<double> res { sy[0] * sy[1] + (sy[5] + sy[6]) * sy[7]
            - sy[11] * G({ 0, x / a }, 1) + sy[1] * (sy[2] + G({ 0, x / a, 1 }, 1))
            + 3. * G({ 0, 0, 0, 1, c / a, d / a }, 1)
            + 3. * G({ 0, 0, 0, c / a, 1, d / a }, 1)
            + 3. * G({ 0, 0, 0, c / a, d / a, 1 }, 1)
            + 2. * G({ 0, 0, 1, 0, c / a, d / a }, 1)
            + 2. * G({ 0, 0, 1, c / a, 0, d / a }, 1)
            + 2. * G({ 0, 0, 1, c / a, d / a, x / a }, 1)
            + 2. * G({ 0, 0, c / a, 0, 1, d / a }, 1)
            + 2. * G({ 0, 0, c / a, 0, d / a, 1 }, 1)
            + 2. * G({ 0, 0, c / a, 1, 0, d / a }, 1)
            + 2. * G({ 0, 0, c / a, 1, d / a, x / a }, 1)
            + 2. * G({ 0, 0, c / a, d / a, 1, x / a }, 1)
            + 2. * G({ 0, 0, c / a, d / a, x / a, 1 }, 1)
            + G({ 0, 1, 0, 0, c / a, d / a }, 1) + G({ 0, 1, 0, c / a, 0, d / a }, 1)
            + G({ 0, 1, 0, c / a, d / a, x / a }, 1) + G({ 0, 1, c / a, 0, 0, d / a }, 1)
            + G({ 0, 1, c / a, 0, d / a, x / a }, 1)
            + G({ 0, 1, c / a, d / a, 0, x / a }, 1) + G({ 0, c / a, 0, 0, 1, d / a }, 1)
            + G({ 0, c / a, 0, 0, d / a, 1 }, 1) + G({ 0, c / a, 0, 1, 0, d / a }, 1)
            + G({ 0, c / a, 0, 1, d / a, x / a }, 1)
            + G({ 0, c / a, 0, d / a, 1, x / a }, 1)
            + G({ 0, c / a, 0, d / a, x / a, 1 }, 1) + G({ 0, c / a, 1, 0, 0, d / a }, 1)
            + G({ 0, c / a, 1, 0, d / a, x / a }, 1)
            + G({ 0, c / a, 1, d / a, 0, x / a }, 1)
            + G({ 0, c / a, d / a, 0, 1, x / a }, 1)
            + G({ 0, c / a, d / a, 0, x / a, 1 }, 1)
            + G({ 0, c / a, d / a, 1, 0, x / a }, 1)
            + (2. * sy[9] + sy[10] + sy[4] + sy[5] + sy[6] + 2. * sy[8])
                * G({ 0, d }, { 1, sd }, x)
            + (-sy[0] - sy[2]) * G({ 0, 0, d }, { 1, 1, sd }, x)
            + G({ 0, 0, a, 0, c, d }, { 1, 1, sa, 1, sc, sd }, x)
            + sy[3] * (sy[4] + sy[5] + sy[6]) * Log(-x, sa)
            + (-sy[4] / 2. - sy[5] / 2. - sy[6] / 2.) * pow(sy[3], 2.) - sy[11] * Zeta(2)
            + 2. * sy[5] * Zeta(2) + 2. * sy[6] * Zeta(2)
            + G({ 0, 0, c, d }, { 1, 1, sc, sd }, x) * Zeta(2)
            + sy[4] * (sy[7] + 2. * Zeta(2)) };
        if (c != x) {
            res += (-2. * sy[9] - sy[10] - 2. * sy[8] + 2. * G({ 0, 0, 1, x / a }, 1)
                       + 2. * G({ 0, 0, x / a, 1 }, 1) + G({ 0, 1, 0, x / a }, 1))
                * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (-2. * G({ 0, 0, 1, c / a, d / a }, 1)
                       + 2. * G({ 0, 0, 1, c / a, x / a }, 1)
                       - 2. * G({ 0, 0, c / a, 1, d / a }, 1)
                       + 2. * G({ 0, 0, c / a, 1, x / a }, 1)
                       - 2. * G({ 0, 0, c / a, d / a, 1 }, 1)
                       + 2. * G({ 0, 0, c / a, x / a, 1 }, 1)
                       - G({ 0, 1, 0, c / a, d / a }, 1) + G({ 0, 1, 0, c / a, x / a }, 1)
                       - G({ 0, 1, c / a, 0, d / a }, 1) + G({ 0, 1, c / a, 0, x / a }, 1)
                       - G({ 0, c / a, 0, 1, d / a }, 1) + G({ 0, c / a, 0, 1, x / a }, 1)
                       - G({ 0, c / a, 0, d / a, 1 }, 1) + G({ 0, c / a, 0, x / a, 1 }, 1)
                       - G({ 0, c / a, 1, 0, d / a }, 1)
                       + G({ 0, c / a, 1, 0, x / a }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    } else { // abcd
        const vector<complex<double>> sy = { G({ 0, c, d }, { 1, sc, sd }, x),
            G({ b / a, 0, c / a }, 1), G({ 0, b / a }, 1), G({ 0, b / a, 0, c / a }, 1),
            G({ b / a, 0, 0, c / a }, 1), G({ b / a, 0, c / a, d / a }, 1), Log(a, sa) };
        complex<double> res { -(sy[0] * G({ 0, b / a, x / a }, 1))
            + sy[0] * (-sy[1] - G({ b / a, 0, x / a }, 1))
            - G({ 0, 0, b / a, 0, c / a, d / a }, 1)
            - 2. * G({ 0, b / a, 0, 0, c / a, d / a }, 1)
            - G({ 0, b / a, 0, c / a, 0, d / a }, 1)
            - G({ 0, b / a, 0, c / a, d / a, x / a }, 1)
            - 3. * G({ b / a, 0, 0, 0, c / a, d / a }, 1)
            - 2. * G({ b / a, 0, 0, c / a, 0, d / a }, 1)
            - 2. * G({ b / a, 0, 0, c / a, d / a, x / a }, 1)
            - G({ b / a, 0, c / a, 0, 0, d / a }, 1)
            - G({ b / a, 0, c / a, 0, d / a, x / a }, 1)
            - G({ b / a, 0, c / a, d / a, 0, x / a }, 1)
            + (-sy[3] - 2. * sy[4] - sy[5]) * G({ 0, d }, { 1, sd }, x)
            + sy[1] * G({ 0, 0, d }, { 1, 1, sd }, x)
            - sy[2] * G({ 0, 0, c, d }, { 1, 1, sc, sd }, x)
            + G({ b / a }, 1)
                * (G({ 0, 0, 0, c, d }, { 1, 1, 1, sc, sd }, x)
                    - G({ 0, b, 0, c, d }, { 1, sb, 1, sc, sd }, x))
            + G({ 0, 0, b, 0, c, d }, { 1, 1, sb, 1, sc, sd }, x)
            - sy[5] * sy[6] * Log(-x, sa) + (sy[5] * pow(sy[6], 2.)) / 2.
            + sy[5] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
        if (b != x) {
            res += (sy[2] - G({ 0, x / a }, 1)) * G({ b, 0, c, d }, { sb, 1, sc, sd }, x);
        }
        if (c != x) {
            res += (sy[3] + 2. * sy[4] - G({ 0, b / a, 0, x / a }, 1)
                       - 2. * G({ b / a, 0, 0, x / a }, 1))
                * G({ c, d }, { sc, sd }, x);
        }
        if (d != x) {
            res += (G({ 0, b / a, 0, c / a, d / a }, 1)
                       - G({ 0, b / a, 0, c / a, x / a }, 1)
                       + 2. * G({ b / a, 0, 0, c / a, d / a }, 1)
                       - 2. * G({ b / a, 0, 0, c / a, x / a }, 1)
                       + G({ b / a, 0, c / a, 0, d / a }, 1)
                       - G({ b / a, 0, c / a, 0, x / a }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G6_0abc0d_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy
        = { G({ 0, c / d, b / d }, 1), G({ 0, a, b }, { 1, sa, sb }, x),
              G({ 0, c / d }, 1), G({ 0, a, b, c }, { 1, sa, sb, sc }, x),
              G({ 0, a }, { 1, sa }, x), G({ 0, c / d, b / d, a / d }, 1), Log(d, sd),
              G({ x / d }, 1), G({ 0, 0, a, b, c }, { 1, 1, sa, sb, sc }, x),
              G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x),
              G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x),
              G({ 0, 0, c / d, b / d, a / d }, 1), G({ 0, c / d, 0, b / d, a / d }, 1),
              G({ 0, c / d, b / d, 0, a / d }, 1) };
    complex<double> res { (-2. * sy[11] - sy[12] - sy[13] - sy[1] * sy[2] + sy[0] * sy[4])
            * sy[6]
        + (sy[9] + sy[10]) * sy[7] + 2. * sy[7] * sy[8]
        + sy[1] * (sy[0] + 2. * G({ 0, 0, c / d }, 1))
        + sy[4]
            * (-sy[5] - 2. * G({ 0, 0, c / d, b / d }, 1) - G({ 0, c / d, 0, b / d }, 1))
        + 3. * G({ 0, 0, 0, c / d, b / d, a / d }, 1)
        + 2. * G({ 0, 0, c / d, 0, b / d, a / d }, 1)
        + 2. * G({ 0, 0, c / d, b / d, 0, a / d }, 1)
        + G({ 0, c / d, 0, 0, b / d, a / d }, 1) + G({ 0, c / d, 0, b / d, 0, a / d }, 1)
        + G({ 0, c / d, b / d, 0, 0, a / d }, 1)
        - G({ 0, c / d, b / d, a / d, 0, x / d }, 1)
        + (-sy[9] - sy[10] - 2. * sy[8]) * G({ d }, { sd }, x)
        + sy[3] * (G({ 0, x / d }, 1) + G({ 0, d }, { 1, sd }, x))
        + 2. * sy[0] * G({ 0, 0, a }, { 1, 1, sa }, x)
        + sy[2]
            * (-sy[3] - 2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
        + 3. * G({ 0, 0, 0, a, b, c }, { 1, 1, 1, sa, sb, sc }, x)
        + 2. * G({ 0, 0, a, 0, b, c }, { 1, 1, sa, 1, sb, sc }, x)
        + 2. * G({ 0, 0, a, b, 0, c }, { 1, 1, sa, sb, 1, sc }, x)
        + G({ 0, a, 0, 0, b, c }, { 1, sa, 1, 1, sb, sc }, x)
        + G({ 0, a, 0, b, 0, c }, { 1, sa, 1, sb, 1, sc }, x)
        + G({ 0, a, b, 0, 0, c }, { 1, sa, sb, 1, 1, sc }, x)
        + (2. * sy[11] + sy[12] + sy[13] + sy[1] * sy[2] - sy[0] * sy[4] - sy[5] * sy[6])
            * Log(-x, sd)
        + (sy[5] * pow(sy[6], 2.)) / 2.
        + sy[5] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
    return res;
}
complex<double> G6_0abc0d_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ 0, a, d }, { 1, sa, sd }, x),
        G({ 0, b / c, a / c }, 1), G({ 0, 0, a }, { 1, 1, sa }, x),
        G({ 0, b / c, d / c }, 1), G({ 0, a, b }, { 1, sa, sb }, x),
        G({ 0, d / c, b / c }, 1), G({ b / c, 0, a / c }, 1), G({ b / c, 0, d / c }, 1),
        G({ 0, a, 0, d }, { 1, sa, 1, sd }, x), G({ 0, d / c }, 1),
        G({ 0, a, b, d }, { 1, sa, sb, sd }, x), G({ 0, a }, { 1, sa }, x),
        G({ 0, b / c, d / c, a / c }, 1), G({ 0, d / c, b / c, a / c }, 1),
        G({ b / c, 0, d / c, a / c }, 1), G({ 0, 0 }, { 1, 1 }, x),
        G({ b / c, 0, a / c, d / c }, 1), G({ b / c, a / c, 0, d / c }, 1), Log(c, sc),
        G({ 0, b / c, a / c, d / c }, 1), G({ 0, d }, { 1, sd }, x),
        G({ 0, 0, b / c, a / c, d / c }, 1), G({ 0, 0, b / c, d / c, a / c }, 1),
        G({ 0, 0, d / c, b / c, a / c }, 1), G({ 0, b / c, 0, a / c, d / c }, 1),
        G({ 0, b / c, 0, d / c, a / c }, 1), G({ 0, b / c, a / c, 0, d / c }, 1),
        G({ 0, b / c, d / c, 0, a / c }, 1), G({ 0, d / c, 0, b / c, a / c }, 1),
        G({ 0, d / c, b / c, 0, a / c }, 1), G({ b / c, 0, 0, a / c, d / c }, 1),
        G({ b / c, 0, 0, d / c, a / c }, 1), G({ b / c, 0, a / c, 0, d / c }, 1),
        G({ b / c, 0, d / c, 0, a / c }, 1) };
    complex<double> res { sy[15] * (-sy[12] - sy[13] - sy[14] - sy[16] - sy[17])
        - sy[0] * sy[1] + 2. * sy[2] * sy[3] + sy[0] * (sy[3] - sy[6] + sy[7])
        + sy[2] * (2. * sy[5] + 2. * sy[7])
        + sy[18]
            * (-2. * sy[21] - 2. * sy[22] - 2. * sy[23] - 2. * sy[24] - 2. * sy[25]
                - sy[26] - sy[27] - sy[28] - sy[29] - 2. * sy[30] - 2. * sy[31] - sy[32]
                - sy[33] - sy[9] * sy[4] + sy[11] * (sy[3] + sy[5] + sy[7]))
        + (sy[10] - sy[8]) * G({ 0, b / c }, 1)
        + sy[4] * (sy[5] + 2. * G({ 0, 0, d / c }, 1))
        + sy[11]
            * (-sy[12] - sy[13] - sy[14] - 2. * G({ 0, 0, b / c, d / c }, 1)
                - 2. * G({ 0, 0, d / c, b / c }, 1) - 2. * G({ 0, b / c, 0, d / c }, 1)
                - G({ 0, d / c, 0, b / c }, 1) - 2. * G({ b / c, 0, 0, d / c }, 1))
        + sy[20] * (-sy[16] - sy[17] - G({ b / c, a / c, 0, x / c }, 1))
        + 3. * G({ 0, 0, 0, b / c, a / c, d / c }, 1)
        + 3. * G({ 0, 0, 0, b / c, d / c, a / c }, 1)
        + 3. * G({ 0, 0, 0, d / c, b / c, a / c }, 1)
        + 3. * G({ 0, 0, b / c, 0, a / c, d / c }, 1)
        + 3. * G({ 0, 0, b / c, 0, d / c, a / c }, 1)
        + G({ 0, 0, b / c, a / c, 0, d / c }, 1)
        + 2. * G({ 0, 0, b / c, d / c, 0, a / c }, 1)
        + 2. * G({ 0, 0, d / c, 0, b / c, a / c }, 1)
        + 2. * G({ 0, 0, d / c, b / c, 0, a / c }, 1)
        + 3. * G({ 0, b / c, 0, 0, a / c, d / c }, 1)
        + 3. * G({ 0, b / c, 0, 0, d / c, a / c }, 1)
        + G({ 0, b / c, 0, a / c, 0, d / c }, 1)
        + 2. * G({ 0, b / c, 0, d / c, 0, a / c }, 1)
        - G({ 0, b / c, a / c, 0, 0, d / c }, 1)
        - G({ 0, b / c, a / c, 0, d / c, x / c }, 1)
        - G({ 0, b / c, a / c, d / c, 0, x / c }, 1)
        + G({ 0, b / c, d / c, 0, 0, a / c }, 1)
        - G({ 0, b / c, d / c, a / c, 0, x / c }, 1)
        + G({ 0, d / c, 0, 0, b / c, a / c }, 1) + G({ 0, d / c, 0, b / c, 0, a / c }, 1)
        + G({ 0, d / c, b / c, 0, 0, a / c }, 1)
        - G({ 0, d / c, b / c, a / c, 0, x / c }, 1)
        + 3. * G({ b / c, 0, 0, 0, a / c, d / c }, 1)
        + 3. * G({ b / c, 0, 0, 0, d / c, a / c }, 1)
        + G({ b / c, 0, 0, a / c, 0, d / c }, 1)
        + 2. * G({ b / c, 0, 0, d / c, 0, a / c }, 1)
        - G({ b / c, 0, a / c, 0, 0, d / c }, 1)
        - G({ b / c, 0, a / c, 0, d / c, x / c }, 1)
        - G({ b / c, 0, a / c, d / c, 0, x / c }, 1)
        + G({ b / c, 0, d / c, 0, 0, a / c }, 1)
        - G({ b / c, 0, d / c, a / c, 0, x / c }, 1)
        - 3. * G({ b / c, a / c, 0, 0, 0, d / c }, 1)
        - 2. * G({ b / c, a / c, 0, 0, d / c, x / c }, 1)
        - G({ b / c, a / c, 0, d / c, 0, x / c }, 1)
        + (sy[1] + sy[6]) * G({ 0, 0, d }, { 1, 1, sd }, x)
        + G({ b / c, a / c }, 1) * (-sy[8] + G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x))
        + sy[9]
            * (-sy[10] - 2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
        + G({ b / c }, 1)
            * (-G({ 0, a, 0, 0, d }, { 1, sa, 1, 1, sd }, x)
                + G({ 0, a, b, 0, d }, { 1, sa, sb, 1, sd }, x))
        + G({ 0, a, b, 0, 0, d }, { 1, sa, sb, 1, 1, sd }, x)
        + (sy[18] * (-sy[12] - sy[13] - sy[14] - sy[16] - sy[17] - sy[19]) + 2. * sy[21]
              + 2. * sy[22] + 2. * sy[23] + 2. * sy[24] + 2. * sy[25] + sy[26] + sy[27]
              + sy[28] + sy[29] + 2. * sy[30] + 2. * sy[31] + sy[32] + sy[33]
              + sy[9] * sy[4] + sy[11] * (-sy[3] - sy[5] - sy[7]))
            * Log(-x, sc)
        + (sy[12] / 2. + sy[13] / 2. + sy[14] / 2. + sy[16] / 2. + sy[17] / 2.
              + sy[19] / 2.)
            * pow(sy[18], 2.)
        + sy[19] * (-sy[15] - sy[20] - 2. * Zeta(2)) - 2. * sy[12] * Zeta(2)
        - 2. * sy[13] * Zeta(2) - 2. * sy[14] * Zeta(2) - 2. * sy[16] * Zeta(2)
        - 2. * sy[17] * Zeta(2) };
    if (d != x) {
        res += (sy[26] + sy[32] - G({ 0, b / c, a / c, 0, x / c }, 1)
                   - G({ b / c, 0, a / c, 0, x / c }, 1)
                   + 2. * G({ b / c, a / c, 0, 0, d / c }, 1)
                   - 2. * G({ b / c, a / c, 0, 0, x / c }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}
complex<double> G6_0abc0d_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (b == c) {
        // abbd
        const vector<complex<double>> sy
            = { G({ 0, 1, a / b }, 1), G({ 0, a, d }, { 1, sa, sd }, x),
                  G({ 0, 1, d / b }, 1), G({ 0, 0, a }, { 1, 1, sa }, x),
                  G({ 0, d / b, 1 }, 1), G({ 0, a / b, 1 }, 1), G({ a / b, 0, 1 }, 1),
                  G({ b, 0, d }, { sb, 1, sd }, x), G({ 0, a }, { 1, sa }, x),
                  G({ 0, 1, d / b, a / b }, 1), G({ 0, d / b, 1, a / b }, 1),
                  G({ 0, d / b, a / b, 1 }, 1), G({ 0, 0 }, { 1, 1 }, x),
                  G({ 0, a / b, 1, d / b }, 1), G({ 0, a / b, d / b, 1 }, 1),
                  G({ a / b, 0, 1, d / b }, 1), G({ a / b, 0, d / b, 1 }, 1), Log(b, sb),
                  G({ 0, 1, a / b, d / b }, 1), G({ 0, d }, { 1, sd }, x),
                  G({ 0, 0, 1, a / b, d / b }, 1), G({ 0, 0, 1, d / b, a / b }, 1),
                  G({ 0, 0, a / b, 1, d / b }, 1), G({ 0, 0, a / b, d / b, 1 }, 1),
                  G({ 0, 0, d / b, 1, a / b }, 1), G({ 0, 0, d / b, a / b, 1 }, 1),
                  G({ 0, 1, 0, a / b, d / b }, 1), G({ 0, 1, 0, d / b, a / b }, 1),
                  G({ 0, 1, d / b, 0, a / b }, 1), G({ 0, a / b, 0, 1, d / b }, 1),
                  G({ 0, a / b, 0, d / b, 1 }, 1), G({ 0, d / b, 0, 1, a / b }, 1),
                  G({ 0, d / b, 0, a / b, 1 }, 1), G({ 0, d / b, 1, 0, a / b }, 1) };
        complex<double> res {
            sy[12] * (-sy[9] - sy[10] - sy[11] - sy[13] - sy[14] - sy[15] - sy[16])
            - sy[0] * sy[1] + sy[2] * (sy[1] + 2. * sy[3]) + 2. * sy[3] * sy[4]
            + sy[1] * (sy[4] - sy[5]) - sy[6] * sy[7]
            + sy[17]
                * (-2. * sy[20] - 2. * sy[21] - 2. * sy[22] - 2. * sy[23] - 2. * sy[24]
                    - 2. * sy[25] - sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31]
                    - sy[32] - sy[33] + (sy[2] + sy[4]) * sy[8])
            + sy[7] * G({ a / b, 0, x / b }, 1)
            + sy[8]
                * (-sy[9] - sy[10] - sy[11] - 2. * G({ 0, 0, 1, d / b }, 1)
                    - 2. * G({ 0, 0, d / b, 1 }, 1) - G({ 0, 1, 0, d / b }, 1))
            + sy[19] * (-sy[13] - sy[14] - sy[15] - sy[16] - G({ a / b, 0, x / b, 1 }, 1))
            + 3. * G({ 0, 0, 0, 1, a / b, d / b }, 1)
            + 3. * G({ 0, 0, 0, 1, d / b, a / b }, 1)
            + 3. * G({ 0, 0, 0, a / b, 1, d / b }, 1)
            + 3. * G({ 0, 0, 0, a / b, d / b, 1 }, 1)
            + 3. * G({ 0, 0, 0, d / b, 1, a / b }, 1)
            + 3. * G({ 0, 0, 0, d / b, a / b, 1 }, 1)
            + 2. * G({ 0, 0, 1, 0, a / b, d / b }, 1)
            + 2. * G({ 0, 0, 1, 0, d / b, a / b }, 1)
            + 2. * G({ 0, 0, 1, d / b, 0, a / b }, 1) + G({ 0, 0, a / b, 0, 1, d / b }, 1)
            + G({ 0, 0, a / b, 0, d / b, 1 }, 1) + 2. * G({ 0, 0, d / b, 0, 1, a / b }, 1)
            + 2. * G({ 0, 0, d / b, 0, a / b, 1 }, 1)
            + 2. * G({ 0, 0, d / b, 1, 0, a / b }, 1) + G({ 0, 1, 0, 0, a / b, d / b }, 1)
            + G({ 0, 1, 0, 0, d / b, a / b }, 1) + G({ 0, 1, 0, d / b, 0, a / b }, 1)
            - G({ 0, 1, a / b, 0, 0, d / b }, 1) - G({ 0, 1, a / b, 0, d / b, x / b }, 1)
            - G({ 0, 1, a / b, d / b, 0, x / b }, 1) + G({ 0, 1, d / b, 0, 0, a / b }, 1)
            - G({ 0, 1, d / b, a / b, 0, x / b }, 1) - G({ 0, a / b, 0, 0, 1, d / b }, 1)
            - G({ 0, a / b, 0, 0, d / b, 1 }, 1) - G({ 0, a / b, 0, 1, 0, d / b }, 1)
            - G({ 0, a / b, 0, 1, d / b, x / b }, 1)
            - G({ 0, a / b, 0, d / b, 1, x / b }, 1)
            - G({ 0, a / b, 0, d / b, x / b, 1 }, 1) - G({ 0, a / b, 1, 0, 0, d / b }, 1)
            - G({ 0, a / b, 1, 0, d / b, x / b }, 1)
            - G({ 0, a / b, 1, d / b, 0, x / b }, 1)
            - G({ 0, a / b, d / b, 0, 1, x / b }, 1)
            - G({ 0, a / b, d / b, 0, x / b, 1 }, 1)
            - G({ 0, a / b, d / b, 1, 0, x / b }, 1) + G({ 0, d / b, 0, 0, 1, a / b }, 1)
            + G({ 0, d / b, 0, 0, a / b, 1 }, 1) + G({ 0, d / b, 0, 1, 0, a / b }, 1)
            + G({ 0, d / b, 1, 0, 0, a / b }, 1) - G({ 0, d / b, 1, a / b, 0, x / b }, 1)
            - G({ 0, d / b, a / b, 0, 1, x / b }, 1)
            - G({ 0, d / b, a / b, 0, x / b, 1 }, 1)
            - G({ 0, d / b, a / b, 1, 0, x / b }, 1)
            - 3. * G({ a / b, 0, 0, 0, 1, d / b }, 1)
            - 3. * G({ a / b, 0, 0, 0, d / b, 1 }, 1)
            - 2. * G({ a / b, 0, 0, 1, 0, d / b }, 1)
            - 2. * G({ a / b, 0, 0, 1, d / b, x / b }, 1)
            - 2. * G({ a / b, 0, 0, d / b, 1, x / b }, 1)
            - 2. * G({ a / b, 0, 0, d / b, x / b, 1 }, 1)
            - G({ a / b, 0, 1, 0, 0, d / b }, 1) - G({ a / b, 0, 1, 0, d / b, x / b }, 1)
            - G({ a / b, 0, 1, d / b, 0, x / b }, 1)
            - G({ a / b, 0, d / b, 0, 1, x / b }, 1)
            - G({ a / b, 0, d / b, 0, x / b, 1 }, 1)
            - G({ a / b, 0, d / b, 1, 0, x / b }, 1)
            + (sy[0] + sy[5] + sy[6]) * G({ 0, 0, d }, { 1, 1, sd }, x)
            + G({ a / b, 1 }, 1)
                * (-G({ 0, a, 0, d }, { 1, sa, 1, sd }, x)
                    + G({ 0, b, 0, d }, { 1, sb, 1, sd }, x))
            + G({ a / b }, 1)
                * (-G({ 0, 0, b, 0, d }, { 1, 1, sb, 1, sd }, x)
                    + G({ 0, a, b, 0, d }, { 1, sa, sb, 1, sd }, x))
            + G({ 0, a, 0, b, 0, d }, { 1, sa, 1, sb, 1, sd }, x)
            + (sy[17]
                      * (-sy[9] - sy[10] - sy[11] - sy[13] - sy[14] - sy[15] - sy[16]
                          - sy[18])
                  + 2. * sy[20] + 2. * sy[21] + 2. * sy[22] + 2. * sy[23] + 2. * sy[24]
                  + 2. * sy[25] + sy[26] + sy[27] + sy[28] + sy[29] + sy[30] + sy[31]
                  + sy[32] + sy[33] + (-sy[2] - sy[4]) * sy[8])
                * Log(-x, sb)
            + (sy[9] / 2. + sy[10] / 2. + sy[11] / 2. + sy[13] / 2. + sy[14] / 2.
                  + sy[15] / 2. + sy[16] / 2. + sy[18] / 2.)
                * pow(sy[17], 2.)
            + sy[18] * (-sy[12] - sy[19] - 2. * Zeta(2)) - 2. * sy[9] * Zeta(2)
            - 2. * sy[10] * Zeta(2) - 2. * sy[11] * Zeta(2) - 2. * sy[13] * Zeta(2)
            - 2. * sy[14] * Zeta(2) - 2. * sy[15] * Zeta(2) - 2. * sy[16] * Zeta(2)
        };
        if (d != x) {
            res += (sy[29] + sy[30] + G({ 0, 1, a / b, 0, d / b }, 1)
                       - G({ 0, 1, a / b, 0, x / b }, 1) - G({ 0, a / b, 0, 1, x / b }, 1)
                       - G({ 0, a / b, 0, x / b, 1 }, 1) + G({ 0, a / b, 1, 0, d / b }, 1)
                       - G({ 0, a / b, 1, 0, x / b }, 1)
                       + 2. * G({ a / b, 0, 0, 1, d / b }, 1)
                       - 2. * G({ a / b, 0, 0, 1, x / b }, 1)
                       + 2. * G({ a / b, 0, 0, d / b, 1 }, 1)
                       - 2. * G({ a / b, 0, 0, x / b, 1 }, 1)
                       + G({ a / b, 0, 1, 0, d / b }, 1)
                       - G({ a / b, 0, 1, 0, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;

    } else { // abcd
        const vector<complex<double>> sy = { G({ a / b, 0, c / b }, 1),
            G({ c / b, 0, a / b }, 1), G({ 0, a, d }, { 1, sa, sd }, x),
            G({ c / b, 0, d / b }, 1), G({ c / b, a / b }, 1),
            G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x), G({ 0, a }, { 1, sa }, x),
            G({ c / b, 0, d / b, a / b }, 1), Log(b, sb),
            G({ a / b, c / b, 0, d / b }, 1), G({ c / b, 0, a / b, d / b }, 1),
            G({ c / b, a / b, 0, d / b }, 1), G({ 0, 0 }, { 1, 1 }, x),
            G({ 0, a, c, 0, d }, { 1, sa, sc, 1, sd }, x),
            G({ 0, a / b, c / b, 0, d / b }, 1), G({ 0, c / b, 0, a / b, d / b }, 1),
            G({ 0, c / b, 0, d / b, a / b }, 1), G({ 0, c / b, a / b, 0, d / b }, 1),
            G({ c / b, 0, 0, a / b, d / b }, 1), G({ c / b, 0, 0, d / b, a / b }, 1),
            G({ c / b, 0, a / b, 0, d / b }, 1), G({ c / b, 0, d / b, 0, a / b }, 1) };
        complex<double> res { sy[1] * sy[2] - sy[2] * sy[3] - sy[4] * sy[5]
            + sy[12] * (sy[10] + sy[11] + sy[7])
            + (sy[14] + sy[15] + sy[16] + sy[17] + 2. * sy[18] + 2. * sy[19] + sy[20]
                  + sy[21] - sy[3] * sy[6])
                * sy[8]
            + sy[6]
                * (sy[7] + G({ 0, c / b, 0, d / b }, 1)
                    + 2. * G({ c / b, 0, 0, d / b }, 1))
            - G({ 0, 0, a / b, c / b, 0, d / b }, 1)
            - G({ 0, 0, c / b, 0, a / b, d / b }, 1)
            - G({ 0, 0, c / b, 0, d / b, a / b }, 1)
            - G({ 0, 0, c / b, a / b, 0, d / b }, 1)
            - 2. * G({ 0, c / b, 0, 0, a / b, d / b }, 1)
            - 2. * G({ 0, c / b, 0, 0, d / b, a / b }, 1)
            - G({ 0, c / b, 0, a / b, 0, d / b }, 1)
            - G({ 0, c / b, 0, d / b, 0, a / b }, 1)
            + G({ a / b, 0, 0, c / b, 0, d / b }, 1)
            + 2. * G({ a / b, 0, c / b, 0, 0, d / b }, 1)
            + G({ a / b, 0, c / b, 0, d / b, x / b }, 1)
            + 3. * G({ a / b, c / b, 0, 0, 0, d / b }, 1)
            + 2. * G({ a / b, c / b, 0, 0, d / b, x / b }, 1)
            + G({ a / b, c / b, 0, d / b, 0, x / b }, 1)
            - 3. * G({ c / b, 0, 0, 0, a / b, d / b }, 1)
            - 3. * G({ c / b, 0, 0, 0, d / b, a / b }, 1)
            - G({ c / b, 0, 0, a / b, 0, d / b }, 1)
            - 2. * G({ c / b, 0, 0, d / b, 0, a / b }, 1)
            + G({ c / b, 0, a / b, 0, 0, d / b }, 1)
            + G({ c / b, 0, a / b, 0, d / b, x / b }, 1)
            + G({ c / b, 0, a / b, d / b, 0, x / b }, 1)
            - G({ c / b, 0, d / b, 0, 0, a / b }, 1)
            + G({ c / b, 0, d / b, a / b, 0, x / b }, 1)
            + 3. * G({ c / b, a / b, 0, 0, 0, d / b }, 1)
            + 2. * G({ c / b, a / b, 0, 0, d / b, x / b }, 1)
            + G({ c / b, a / b, 0, d / b, 0, x / b }, 1)
            + (sy[9] + sy[10] + sy[11] + G({ a / b, 0, c / b, x / b }, 1)
                  + G({ a / b, c / b, 0, x / b }, 1) + G({ c / b, a / b, 0, x / b }, 1))
                * G({ 0, d }, { 1, sd }, x)
            - 2. * sy[3] * G({ 0, 0, a }, { 1, 1, sa }, x)
            + (sy[0] - sy[1]) * G({ 0, 0, d }, { 1, 1, sd }, x)
            + sy[4] * G({ 0, a, 0, d }, { 1, sa, 1, sd }, x)
            + G({ a / b, c / b }, 1) * (-sy[5] + G({ 0, c, 0, d }, { 1, sc, 1, sd }, x))
            + G({ a / b }, 1) * (sy[13] - G({ 0, 0, c, 0, d }, { 1, 1, sc, 1, sd }, x))
            + G({ c / b }, 1) * (-sy[13] + G({ 0, a, 0, 0, d }, { 1, sa, 1, 1, sd }, x))
            + G({ 0, a, 0, c, 0, d }, { 1, sa, 1, sc, 1, sd }, x)
            + (-sy[14] - sy[15] - sy[16] - sy[17] - 2. * sy[18] - 2. * sy[19] - sy[20]
                  - sy[21] + sy[3] * sy[6] + (sy[9] + sy[10] + sy[11] + sy[7]) * sy[8])
                * Log(-x, sb)
            + (-sy[9] / 2. - sy[10] / 2. - sy[11] / 2. - sy[7] / 2.) * pow(sy[8], 2.)
            + 2. * sy[10] * Zeta(2) + 2. * sy[11] * Zeta(2) + 2. * sy[7] * Zeta(2)
            + sy[9] * (sy[12] + 2. * Zeta(2)) };
        if (c != x) {
            res += (-sy[0] + G({ a / b, 0, x / b }, 1))
                * G({ c, 0, d }, { sc, 1, sd }, x);
        }
        if (d != x) {
            res += (-sy[20] - G({ a / b, 0, c / b, 0, d / b }, 1)
                       + G({ a / b, 0, c / b, 0, x / b }, 1)
                       - 2. * G({ a / b, c / b, 0, 0, d / b }, 1)
                       + 2. * G({ a / b, c / b, 0, 0, x / b }, 1)
                       + G({ c / b, 0, a / b, 0, x / b }, 1)
                       - 2. * G({ c / b, a / b, 0, 0, d / b }, 1)
                       + 2. * G({ c / b, a / b, 0, 0, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_0abc0d_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (a == b) {
        if (a == c) { // aaad
            const vector<complex<double>> sy = { G({ a, 0, d }, { sa, 1, sd }, x),
                G({ 0, 0 }, { 1, 1 }, x), G({ 0, 1, d / a, 1 }, 1),
                G({ 0, d / a, 1, 1 }, 1), Log(a, sa), G({ 0, 1, 1, d / a }, 1),
                G({ 0, d }, { 1, sd }, x), G({ a, a, 0, d }, { sa, sa, 1, sd }, x) };
            complex<double> res { sy[1] * (-sy[2] - sy[3]) - sy[7] * G({ 0, x / a }, 1)
                + sy[0] * G({ 0, x / a, 1 }, 1)
                + sy[6] * (-sy[2] - sy[3] - G({ 0, x / a, 1, 1 }, 1))
                - 3. * G({ 0, 0, 0, 1, 1, d / a }, 1)
                - 3. * G({ 0, 0, 0, 1, d / a, 1 }, 1)
                - 3. * G({ 0, 0, 0, d / a, 1, 1 }, 1)
                - 2. * G({ 0, 0, 1, 0, 1, d / a }, 1)
                - 2. * G({ 0, 0, 1, 0, d / a, 1 }, 1)
                - 2. * G({ 0, 0, 1, 1, 0, d / a }, 1)
                - 2. * G({ 0, 0, 1, 1, d / a, x / a }, 1)
                - 2. * G({ 0, 0, 1, d / a, 1, x / a }, 1)
                - 2. * G({ 0, 0, 1, d / a, x / a, 1 }, 1)
                - 2. * G({ 0, 0, d / a, 1, 1, x / a }, 1)
                - 2. * G({ 0, 0, d / a, 1, x / a, 1 }, 1)
                - 2. * G({ 0, 0, d / a, x / a, 1, 1 }, 1) - G({ 0, 1, 0, 0, 1, d / a }, 1)
                - G({ 0, 1, 0, 0, d / a, 1 }, 1) - G({ 0, 1, 0, 1, 0, d / a }, 1)
                - G({ 0, 1, 0, 1, d / a, x / a }, 1) - G({ 0, 1, 0, d / a, 1, x / a }, 1)
                - G({ 0, 1, 0, d / a, x / a, 1 }, 1) - G({ 0, 1, 1, 0, 0, d / a }, 1)
                - G({ 0, 1, 1, 0, d / a, x / a }, 1) - G({ 0, 1, 1, d / a, 0, x / a }, 1)
                - G({ 0, 1, d / a, 0, 1, x / a }, 1) - G({ 0, 1, d / a, 0, x / a, 1 }, 1)
                - G({ 0, 1, d / a, 1, 0, x / a }, 1) - G({ 0, d / a, 0, 1, 1, x / a }, 1)
                - G({ 0, d / a, 0, 1, x / a, 1 }, 1) - G({ 0, d / a, 0, x / a, 1, 1 }, 1)
                - G({ 0, d / a, 1, 0, 1, x / a }, 1) - G({ 0, d / a, 1, 0, x / a, 1 }, 1)
                - G({ 0, d / a, 1, 1, 0, x / a }, 1)
                + G({ 0, 0, a, a, 0, d }, { 1, 1, sa, sa, 1, sd }, x)
                + sy[4] * (-sy[2] - sy[3] - sy[5]) * Log(-x, sa)
                + (sy[2] / 2. + sy[3] / 2. + sy[5] / 2.) * pow(sy[4], 2.)
                + sy[5] * (-sy[1] - sy[6] - 2. * Zeta(2)) - 2. * sy[2] * Zeta(2)
                - 2. * sy[3] * Zeta(2) - sy[7] * Zeta(2)
                + G({ 0, a, 0, d }, { 1, sa, 1, sd }, x) * Zeta(2) - sy[0] * Zeta(3)
                + G({ 0, 0, d }, { 1, 1, sd }, x) * Zeta(3) };
            if (d != x) {
                res += (2. * G({ 0, 0, 1, 1, d / a }, 1)
                           - 2. * G({ 0, 0, 1, 1, x / a }, 1)
                           + 2. * G({ 0, 0, 1, d / a, 1 }, 1)
                           - 2. * G({ 0, 0, 1, x / a, 1 }, 1)
                           + 2. * G({ 0, 0, d / a, 1, 1 }, 1)
                           - 2. * G({ 0, 0, x / a, 1, 1 }, 1)
                           + G({ 0, 1, 0, 1, d / a }, 1) - G({ 0, 1, 0, 1, x / a }, 1)
                           + G({ 0, 1, 0, d / a, 1 }, 1) - G({ 0, 1, 0, x / a, 1 }, 1)
                           + G({ 0, 1, 1, 0, d / a }, 1) - G({ 0, 1, 1, 0, x / a }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        } else { // aacd

            const vector<complex<double>> sy
                = { G({ 0, c / a, 1 }, 1), G({ 0, c, 0, d }, { 1, sc, 1, sd }, x),
                      G({ a, c, 0, d }, { sa, sc, 1, sd }, x), Log(a, sa),
                      G({ c / a, 0, 1, d / a }, 1), G({ c / a, 0, d / a, 1 }, 1),
                      G({ c / a, 1, 0, d / a }, 1), G({ 0, 0 }, { 1, 1 }, x) };
            complex<double> res { (sy[5] + sy[6]) * sy[7] - sy[2] * G({ 0, x / a }, 1)
                + G({ 0, 0, c / a, 0, 1, d / a }, 1) + G({ 0, 0, c / a, 0, d / a, 1 }, 1)
                + G({ 0, 0, c / a, 1, 0, d / a }, 1)
                + 2. * G({ 0, c / a, 0, 0, 1, d / a }, 1)
                + 2. * G({ 0, c / a, 0, 0, d / a, 1 }, 1)
                + 2. * G({ 0, c / a, 0, 1, 0, d / a }, 1)
                + G({ 0, c / a, 0, 1, d / a, x / a }, 1)
                + G({ 0, c / a, 0, d / a, 1, x / a }, 1)
                + G({ 0, c / a, 0, d / a, x / a, 1 }, 1)
                + 2. * G({ 0, c / a, 1, 0, 0, d / a }, 1)
                + G({ 0, c / a, 1, 0, d / a, x / a }, 1)
                + 3. * G({ c / a, 0, 0, 0, 1, d / a }, 1)
                + 3. * G({ c / a, 0, 0, 0, d / a, 1 }, 1)
                + 3. * G({ c / a, 0, 0, 1, 0, d / a }, 1)
                + 2. * G({ c / a, 0, 0, 1, d / a, x / a }, 1)
                + 2. * G({ c / a, 0, 0, d / a, 1, x / a }, 1)
                + 2. * G({ c / a, 0, 0, d / a, x / a, 1 }, 1)
                + 3. * G({ c / a, 0, 1, 0, 0, d / a }, 1)
                + 2. * G({ c / a, 0, 1, 0, d / a, x / a }, 1)
                + G({ c / a, 0, 1, d / a, 0, x / a }, 1)
                + G({ c / a, 0, d / a, 0, 1, x / a }, 1)
                + G({ c / a, 0, d / a, 0, x / a, 1 }, 1)
                + G({ c / a, 0, d / a, 1, 0, x / a }, 1)
                + 3. * G({ c / a, 1, 0, 0, 0, d / a }, 1)
                + 2. * G({ c / a, 1, 0, 0, d / a, x / a }, 1)
                + G({ c / a, 1, 0, d / a, 0, x / a }, 1)
                + (sy[4] + sy[5] + sy[6] + G({ 0, c / a, 1, x / a }, 1)
                      + G({ 0, c / a, x / a, 1 }, 1) + G({ c / a, 0, 1, x / a }, 1)
                      + G({ c / a, 0, x / a, 1 }, 1) + G({ c / a, 1, 0, x / a }, 1))
                    * G({ 0, d }, { 1, sd }, x)
                + sy[0] * G({ 0, 0, d }, { 1, 1, sd }, x)
                + G({ c / a, 1 }, 1) * (sy[1] - G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x))
                + G({ 0, 0, a, c, 0, d }, { 1, 1, sa, sc, 1, sd }, x)
                + sy[3] * (sy[4] + sy[5] + sy[6]) * Log(-x, sa)
                + (-sy[4] / 2. - sy[5] / 2. - sy[6] / 2.) * pow(sy[3], 2.)
                + sy[1] * Zeta(2) - sy[2] * Zeta(2) + 2. * sy[5] * Zeta(2)
                + 2. * sy[6] * Zeta(2) + sy[4] * (sy[7] + 2. * Zeta(2)) };
            if (c != x) {
                res += (-sy[0] + G({ 0, x / a, 1 }, 1))
                    * G({ c, 0, d }, { sc, 1, sd }, x);
            }
            if (d != x) {
                res += (-G({ 0, c / a, 0, 1, d / a }, 1) + G({ 0, c / a, 0, 1, x / a }, 1)
                           - G({ 0, c / a, 0, d / a, 1 }, 1)
                           + G({ 0, c / a, 0, x / a, 1 }, 1)
                           - G({ 0, c / a, 1, 0, d / a }, 1)
                           + G({ 0, c / a, 1, 0, x / a }, 1)
                           - 2. * G({ c / a, 0, 0, 1, d / a }, 1)
                           + 2. * G({ c / a, 0, 0, 1, x / a }, 1)
                           - 2. * G({ c / a, 0, 0, d / a, 1 }, 1)
                           + 2. * G({ c / a, 0, 0, x / a, 1 }, 1)
                           - 2. * G({ c / a, 0, 1, 0, d / a }, 1)
                           + 2. * G({ c / a, 0, 1, 0, x / a }, 1)
                           - 2. * G({ c / a, 1, 0, 0, d / a }, 1)
                           + 2. * G({ c / a, 1, 0, 0, x / a }, 1))
                    * G({ d }, { sd }, x);
            }
            return res;
        }
    } else { // abcd
        const vector<complex<double>> sy = { G({ 0, b / a, c / a }, 1),
            G({ b / a, 0, c / a }, 1), G({ 0, c, 0, d }, { 1, sc, 1, sd }, x),
            G({ 0, b / a }, 1), Log(a, sa), G({ b / a, c / a, 0, d / a }, 1) };
        complex<double> res { -(sy[2] * sy[3]) - G({ 0, 0, b / a, c / a, 0, d / a }, 1)
            - G({ 0, b / a, 0, c / a, 0, d / a }, 1)
            - 2. * G({ 0, b / a, c / a, 0, 0, d / a }, 1)
            - G({ 0, b / a, c / a, 0, d / a, x / a }, 1)
            - G({ b / a, 0, 0, c / a, 0, d / a }, 1)
            - 2. * G({ b / a, 0, c / a, 0, 0, d / a }, 1)
            - G({ b / a, 0, c / a, 0, d / a, x / a }, 1)
            - 3. * G({ b / a, c / a, 0, 0, 0, d / a }, 1)
            - 2. * G({ b / a, c / a, 0, 0, d / a, x / a }, 1)
            - G({ b / a, c / a, 0, d / a, 0, x / a }, 1)
            + (-sy[5] - G({ 0, b / a, c / a, x / a }, 1)
                  - G({ b / a, 0, c / a, x / a }, 1) - G({ b / a, c / a, 0, x / a }, 1))
                * G({ 0, d }, { 1, sd }, x)
            + (-sy[0] - sy[1]) * G({ 0, 0, d }, { 1, 1, sd }, x)
            + G({ b / a, c / a }, 1) * (-sy[2] + G({ 0, 0, 0, d }, { 1, 1, 1, sd }, x))
            + G({ b / a }, 1)
                * (G({ 0, 0, c, 0, d }, { 1, 1, sc, 1, sd }, x)
                    - G({ 0, b, c, 0, d }, { 1, sb, sc, 1, sd }, x))
            + G({ 0, 0, b, c, 0, d }, { 1, 1, sb, sc, 1, sd }, x)
            - sy[4] * sy[5] * Log(-x, sa) + (sy[5] * pow(sy[4], 2.)) / 2.
            + sy[5] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
        if (b != x) {
            res += (sy[3] - G({ 0, x / a }, 1)) * G({ b, c, 0, d }, { sb, sc, 1, sd }, x);
        }
        if (c != x) {
            res += (sy[0] + sy[1] - G({ 0, b / a, x / a }, 1) - G({ b / a, 0, x / a }, 1))
                * G({ c, 0, d }, { sc, 1, sd }, x);
        }
        if (d != x) {
            res += (G({ 0, b / a, c / a, 0, d / a }, 1)
                       - G({ 0, b / a, c / a, 0, x / a }, 1)
                       + G({ b / a, 0, c / a, 0, d / a }, 1)
                       - G({ b / a, 0, c / a, 0, x / a }, 1)
                       + 2. * G({ b / a, c / a, 0, 0, d / a }, 1)
                       - 2. * G({ b / a, c / a, 0, 0, x / a }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G6_a0b0cd_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy
        = { G({ c / d, 0, b / d }, 1), G({ a }, { sa }, x), G({ 0, b }, { 1, sb }, x),
              G({ c / d, 0, x / d }, 1), G({ a, 0, b }, { sa, 1, sb }, x),
              G({ 0, c / d }, 1), G({ 0, a, 0, b }, { 1, sa, 1, sb }, x),
              G({ a, 0, 0, b }, { sa, 1, 1, sb }, x), G({ c / d, x / d }, 1),
              G({ 0, 0, b }, { 1, 1, sb }, x), G({ 0, a }, { 1, sa }, x), G({ c / d }, 1),
              Log(d, sd), G({ 0, c / d, 0, b / d }, 1), G({ c / d, 0, 0, b / d }, 1),
              G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x),
              G({ c / d, 0, b / d, 0, a / d }, 1), G({ 0, 0 }, { 1, 1 }, x) };
    complex<double> res { sy[1] * sy[2] * (sy[0] - sy[3])
        + sy[12] * (sy[16] + sy[4] * sy[5]) + sy[5] * (sy[6] + 2. * sy[7])
        + sy[10]
            * (sy[9] * sy[11] - sy[0] * sy[12] + sy[13] + 2. * sy[14] + sy[2] * sy[8])
        + sy[4] * (sy[3] - G({ 0, 0, c / d }, 1)) - G({ 0, c / d, 0, b / d, 0, a / d }, 1)
        - 2. * G({ c / d, 0, 0, b / d, 0, a / d }, 1)
        - 2. * G({ c / d, 0, b / d, 0, 0, a / d }, 1)
        - G({ c / d, 0, b / d, 0, a / d, x / d }, 1)
        + sy[15] * (-G({ x / d }, 1) + G({ d }, { sd }, x))
        - sy[0] * G({ 0, 0, a }, { 1, 1, sa }, x)
        + sy[8]
            * (-2. * sy[9] * sy[1] + sy[6] + 2. * sy[7]
                - G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
        - G({ 0, a, 0, b, 0, c }, { 1, sa, 1, sb, 1, sc }, x)
        - 2. * G({ a, 0, 0, b, 0, c }, { sa, 1, 1, sb, 1, sc }, x)
        - 2. * G({ a, 0, b, 0, 0, c }, { sa, 1, sb, 1, 1, sc }, x)
        + (sy[0] * sy[10] - sy[16] + (sy[0] * sy[12] - sy[13] - 2. * sy[14]) * sy[1]
              - sy[4] * sy[5] + sy[11] * (sy[12] * sy[4] + sy[6] + 2. * sy[7]))
            * Log(-x, sd)
        + sy[1]
            * (sy[12] * (sy[13] + 2. * sy[14]) + sy[16] + sy[0] * sy[17]
                - G({ 0, 0, c / d, 0, b / d }, 1) - 2. * G({ 0, c / d, 0, 0, b / d }, 1)
                - 3. * G({ c / d, 0, 0, 0, b / d }, 1) - (sy[0] * pow(sy[12], 2.)) / 2.
                + 2. * sy[0] * Zeta(2))
        + sy[11]
            * (sy[15] + sy[12] * (-sy[6] - 2. * sy[7])
                - 3. * sy[1] * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
                - G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x)
                - (sy[4] * pow(sy[12], 2.)) / 2. + sy[4] * (sy[17] + 2. * Zeta(2))) };
    if (b != x) {
        res += (-(sy[0] * sy[10]) + sy[10] * sy[3]) * G({ b }, { sb }, x)
            + (sy[0] - sy[3]) * G({ b, 0, a }, { sb, 1, sa }, x);
    }
    return res;
}
complex<double> G6_a0b0cd_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (c == d) { // abcc
        const vector<complex<double>> sy = { G({ 0, a }, { 1, sa }, x),
            G({ 0, b }, { 1, sb }, x), G({ x / c, 1 }, 1), G({ a }, { sa }, x),
            G({ 0, 1, b / c }, 1), G({ 0, 1, x / c }, 1), G({ 0, b / c, 1 }, 1),
            G({ 0, x / c, 1 }, 1), G({ a, 0, b }, { sa, 1, sb }, x),
            G({ 0, b / c, 0, a / c }, 1), G({ a, 0, b, c }, { sa, 1, sb, sc }, x),
            G({ c }, { sc }, x), G({ 0, b / c, 0, a / c, 1 }, 1) };
        complex<double> res { -(sy[11] * sy[12]) + sy[0] * sy[1] * sy[2]
            + sy[1] * sy[3] * (sy[4] - sy[5] + sy[6] - sy[7]) - sy[4] * sy[8]
            + sy[5] * sy[8] + sy[7] * sy[8]
            + sy[3] * (sy[12] - G({ 0, b / c, 0, 0, 1 }, 1))
            + sy[11] * G({ 0, b / c, 0, a / c, x / c }, 1)
            + G({ 0, b / c, 0, a / c, 0, 1 }, 1) - G({ 0, b / c, 0, a / c, x / c, 1 }, 1)
            + sy[9] * G({ 0, c }, { 1, sc }, x)
            + (-sy[9] + G({ 0, b / c, 0, 1 }, 1)) * G({ a, c }, { sa, sc }, x)
            - sy[6] * G({ a, 0, c }, { sa, 1, sc }, x)
            + sy[2]
                * (-2. * sy[3] * G({ 0, 0, b }, { 1, 1, sb }, x)
                    + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                    - G({ 0, b, 0, a }, { 1, sb, 1, sa }, x)
                    + 2. * G({ a, 0, 0, b }, { sa, 1, 1, sb }, x))
            + G({ 0, b / c }, 1) * (-sy[10] + G({ a, 0, 0, c }, { sa, 1, 1, sc }, x))
            + G({ a, 0, b, 0, 0, c }, { sa, 1, sb, 1, 1, sc }, x) - sy[10] * Zeta(2)
            + sy[8] * Zeta(3) };
        if (b != x) {
            res += (-(sy[0] * sy[4]) + sy[0] * sy[5] - sy[0] * sy[6] + sy[0] * sy[7])
                    * G({ b }, { sb }, x)
                + (sy[4] - sy[5] + sy[6] - sy[7]) * G({ b, 0, a }, { sb, 1, sa }, x);
        }
        return res;
    } else { // abcd
        const vector<complex<double>> sy = { G({ 0, d / c, b / c }, 1),
            G({ a, 0, b }, { sa, 1, sb }, x), G({ 0, b / c, d / c }, 1),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ d / c, 0, b / c }, 1),
            G({ d / c, 0, x / c }, 1), G({ a }, { sa }, x), G({ 0, b }, { 1, sb }, x),
            G({ 0, b / c, 0, a / c }, 1), G({ a, d }, { sa, sd }, x),
            G({ d / c, x / c }, 1), G({ 0, 0, b }, { 1, 1, sb }, x),
            G({ 0, a, 0, b }, { 1, sa, 1, sb }, x),
            G({ a, 0, 0, b }, { sa, 1, 1, sb }, x),
            G({ a, 0, b, d }, { sa, 1, sb, sd }, x), G({ 0, a }, { 1, sa }, x),
            G({ d / c }, 1), Log(c, sc), G({ 0, 0, b / c, d / c }, 1),
            G({ 0, 0, d / c, b / c }, 1), G({ 0, d / c, 0, b / c }, 1),
            G({ d / c, 0, 0, b / c }, 1), G({ 0, b / c, 0, a / c, d / c }, 1),
            G({ 0, b / c, 0, d / c, a / c }, 1), G({ 0, b / c, d / c, 0, a / c }, 1),
            G({ 0, d / c, b / c, 0, a / c }, 1), G({ d / c, 0, b / c, 0, a / c }, 1),
            G({ 0, 0 }, { 1, 1 }, x) };
        complex<double> res { -(sy[0] * sy[1])
            + sy[17] * (-sy[22] - sy[23] - sy[24] - sy[25] - sy[26])
            + sy[3] * (sy[0] + sy[4]) + (-sy[4] + sy[5]) * sy[6] * sy[7]
            + sy[15]
                * (-(sy[11] * sy[16]) - 2. * sy[18] - 2. * sy[19] - 2. * sy[20]
                    - 2. * sy[21] + sy[17] * (sy[0] + sy[2] + sy[4]) - sy[10] * sy[7])
            - sy[9] * sy[8] + sy[14] * G({ 0, d / c }, 1)
            + sy[1] * (-sy[5] - G({ 0, 0, d / c }, 1))
            + sy[9] * G({ 0, b / c, 0, d / c }, 1)
            + 2. * G({ 0, 0, b / c, 0, a / c, d / c }, 1)
            + 2. * G({ 0, 0, b / c, 0, d / c, a / c }, 1)
            + 2. * G({ 0, 0, b / c, d / c, 0, a / c }, 1)
            + 2. * G({ 0, 0, d / c, b / c, 0, a / c }, 1)
            + 2. * G({ 0, b / c, 0, 0, a / c, d / c }, 1)
            + 2. * G({ 0, b / c, 0, 0, d / c, a / c }, 1)
            + G({ 0, b / c, 0, a / c, 0, d / c }, 1)
            + G({ 0, b / c, 0, a / c, d / c, x / c }, 1)
            + 2. * G({ 0, b / c, 0, d / c, 0, a / c }, 1)
            + G({ 0, b / c, 0, d / c, a / c, x / c }, 1)
            + 2. * G({ 0, b / c, d / c, 0, 0, a / c }, 1)
            + G({ 0, b / c, d / c, 0, a / c, x / c }, 1)
            + 2. * G({ 0, d / c, 0, b / c, 0, a / c }, 1)
            + 2. * G({ 0, d / c, b / c, 0, 0, a / c }, 1)
            + G({ 0, d / c, b / c, 0, a / c, x / c }, 1)
            + 2. * G({ d / c, 0, 0, b / c, 0, a / c }, 1)
            + 2. * G({ d / c, 0, b / c, 0, 0, a / c }, 1)
            + G({ d / c, 0, b / c, 0, a / c, x / c }, 1)
            + sy[8] * G({ 0, d }, { 1, sd }, x)
            + sy[2] * (sy[3] - G({ a, 0, d }, { sa, 1, sd }, x))
            + sy[10]
                * (-sy[12] - 2. * sy[13] + 2. * sy[11] * sy[6]
                    + G({ 0, b, 0, a }, { 1, sb, 1, sa }, x))
            + G({ 0, b / c }, 1) * (-sy[14] + G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
            + G({ a, 0, b, 0, 0, d }, { sa, 1, sb, 1, 1, sd }, x)
            + (sy[16] * (-sy[12] - 2. * sy[13] - sy[17] * sy[1]) + sy[22] + sy[23]
                  + sy[24] + sy[25] + sy[26] + sy[15] * (-sy[0] - sy[2] - sy[4])
                  + (2. * sy[18] + 2. * sy[19] + 2. * sy[20] + 2. * sy[21]
                        + sy[17] * (-sy[0] - sy[2] - sy[4]))
                      * sy[6])
                * Log(-x, sc)
            + sy[16]
                * ((sy[12] + 2. * sy[13]) * sy[17] - sy[1] * sy[27]
                    + 3. * sy[6] * G({ 0, 0, 0, b }, { 1, 1, 1, sb }, x)
                    + G({ 0, 0, b, 0, a }, { 1, 1, sb, 1, sa }, x)
                    - G({ a, 0, b, 0, d }, { sa, 1, sb, 1, sd }, x)
                    + (sy[1] * pow(sy[17], 2.)) / 2. - 2. * sy[1] * Zeta(2))
            + sy[6]
                * (sy[17] * (-2. * sy[18] - 2. * sy[19] - 2. * sy[20] - 2. * sy[21])
                    - sy[23] - sy[24] - sy[25] - sy[26] + sy[27] * (-sy[0] - sy[4])
                    + 3. * G({ 0, 0, 0, b / c, d / c }, 1)
                    + 3. * G({ 0, 0, 0, d / c, b / c }, 1)
                    + 3. * G({ 0, 0, d / c, 0, b / c }, 1)
                    - G({ 0, b / c, 0, 0, d / c }, 1)
                    + 3. * G({ 0, d / c, 0, 0, b / c }, 1)
                    + 3. * G({ d / c, 0, 0, 0, b / c }, 1)
                    + (sy[0] / 2. + sy[2] / 2. + sy[4] / 2.) * pow(sy[17], 2.)
                    + sy[2] * (-sy[27] - 2. * Zeta(2)) - 2. * sy[0] * Zeta(2)
                    - 2. * sy[4] * Zeta(2)) };
        if (b != x) {
            res += (sy[15] * sy[4] - sy[15] * sy[5]) * G({ b }, { sb }, x)
                + (-sy[4] + sy[5]) * G({ b, 0, a }, { sb, 1, sa }, x);
        }
        if (d != x) {
            res += (-sy[22] + G({ 0, b / c, 0, a / c, x / c }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_a0b0cd_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ 0, c, d }, { 1, sc, sd }, x), Log(b, sb),
        G({ a }, { sa }, x), G({ 0, c / b, d / b }, 1), G({ a, 0, d }, { sa, 1, sd }, x),
        G({ a, c, d }, { sa, sc, sd }, x), G({ a, d }, { sa, sd }, x),
        G({ 0, 0, c / b, a / b }, 1), G({ 0, 0, a / b, c / b }, 1),
        G({ 0, a / b, 0, c / b }, 1), G({ 0, c / b, 0, a / b }, 1),
        G({ a, 0, c, d }, { sa, 1, sc, sd }, x), G({ 0, a }, { 1, sa }, x),
        G({ 0, 0, a / b, c / b, d / b }, 1), G({ 0, 0, c / b, a / b, d / b }, 1),
        G({ 0, 0, c / b, d / b, a / b }, 1), G({ 0, a / b, 0, c / b, d / b }, 1),
        G({ 0, c / b, 0, a / b, d / b }, 1), G({ 0, c / b, 0, d / b, a / b }, 1),
        G({ 0, c / b, d / b, 0, a / b }, 1) };
    complex<double> res { sy[1]
            * (-2. * sy[13] - 2. * sy[14] - 2. * sy[15] - sy[16] - sy[17] - sy[18]
                - sy[19] + sy[12] * sy[3])
        - 2. * sy[6] * sy[7] + (2. * sy[0] - 2. * sy[5]) * G({ 0, 0, a / b }, 1)
        + (-2. * sy[4] + 2. * sy[5]) * G({ 0, 0, c / b }, 1)
        + sy[0] * G({ 0, a / b, x / b }, 1)
        + sy[6]
            * (-sy[10] + 2. * G({ 0, 0, c / b, d / b }, 1) + G({ 0, c / b, 0, d / b }, 1))
        + 6. * G({ 0, 0, 0, a / b, c / b, d / b }, 1)
        + 6. * G({ 0, 0, 0, c / b, a / b, d / b }, 1)
        + 6. * G({ 0, 0, 0, c / b, d / b, a / b }, 1)
        + 4. * G({ 0, 0, a / b, 0, c / b, d / b }, 1)
        + 2. * G({ 0, 0, a / b, c / b, 0, d / b }, 1)
        + 2. * G({ 0, 0, a / b, c / b, d / b, x / b }, 1)
        + 4. * G({ 0, 0, c / b, 0, a / b, d / b }, 1)
        + 4. * G({ 0, 0, c / b, 0, d / b, a / b }, 1)
        + 2. * G({ 0, 0, c / b, a / b, 0, d / b }, 1)
        + 2. * G({ 0, 0, c / b, a / b, d / b, x / b }, 1)
        + 4. * G({ 0, 0, c / b, d / b, 0, a / b }, 1)
        + 2. * G({ 0, 0, c / b, d / b, a / b, x / b }, 1)
        + 2. * G({ 0, a / b, 0, 0, c / b, d / b }, 1)
        + G({ 0, a / b, 0, c / b, 0, d / b }, 1)
        + G({ 0, a / b, 0, c / b, d / b, x / b }, 1)
        + 2. * G({ 0, c / b, 0, 0, a / b, d / b }, 1)
        + 2. * G({ 0, c / b, 0, 0, d / b, a / b }, 1)
        + G({ 0, c / b, 0, a / b, 0, d / b }, 1)
        + G({ 0, c / b, 0, a / b, d / b, x / b }, 1)
        + 2. * G({ 0, c / b, 0, d / b, 0, a / b }, 1)
        + G({ 0, c / b, 0, d / b, a / b, x / b }, 1)
        + 2. * G({ 0, c / b, d / b, 0, 0, a / b }, 1)
        + G({ 0, c / b, d / b, 0, a / b, x / b }, 1)
        + (sy[9] + sy[10] + 2. * sy[7] + 2. * sy[8]) * G({ 0, d }, { 1, sd }, x)
        + sy[3] * (-sy[4] + G({ 0, 0, a }, { 1, 1, sa }, x))
        + G({ 0, a / b }, 1) * (-sy[11] + G({ 0, 0, c, d }, { 1, 1, sc, sd }, x))
        + G({ 0, c / b }, 1) * (-sy[11] + G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
        + G({ a, 0, 0, 0, c, d }, { sa, 1, 1, 1, sc, sd }, x)
        + (2. * sy[13] + 2. * sy[14] + 2. * sy[15] + sy[16] + sy[17] + sy[18] + sy[19]
              - sy[12] * sy[3] - sy[1] * sy[2] * sy[3])
            * Log(-x, sb)
        + (sy[2] * sy[3] * pow(sy[1], 2.)) / 2.
        + sy[2]
            * (-2. * sy[15] - sy[18] - sy[19] - 3. * G({ 0, 0, 0, c / b, d / b }, 1)
                - 2. * G({ 0, 0, c / b, 0, d / b }, 1) - G({ 0, c / b, 0, 0, d / b }, 1)
                + sy[3] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
    if (c != x) {
        res += (-sy[9] - 2. * sy[8] + 2. * G({ 0, 0, a / b, x / b }, 1)
                   + G({ 0, a / b, 0, x / b }, 1))
            * G({ c, d }, { sc, sd }, x);
    }
    if (d != x) {
        res += (-2. * sy[13] - 2. * sy[14] - sy[16] - sy[17]
                   + 2. * G({ 0, 0, a / b, c / b, x / b }, 1)
                   + 2. * G({ 0, 0, c / b, a / b, x / b }, 1)
                   + G({ 0, a / b, 0, c / b, x / b }, 1)
                   + G({ 0, c / b, 0, a / b, x / b }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}
complex<double> G6_a0b0cd_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ 0, b / a }, 1), G({ 0, b / a, 0, c / a }, 1),
        G({ 0, b / a, 0, c / a, d / a }, 1) };
    complex<double> res { 2. * G({ 0, 0, b / a, 0, c / a, d / a }, 1)
        + 2. * G({ 0, b / a, 0, 0, c / a, d / a }, 1)
        + G({ 0, b / a, 0, c / a, 0, d / a }, 1)
        + G({ 0, b / a, 0, c / a, d / a, x / a }, 1) + sy[1] * G({ 0, d }, { 1, sd }, x)
        + G({ 0, b / a, x / a }, 1) * G({ 0, c, d }, { 1, sc, sd }, x)
        + sy[0] * G({ 0, 0, c, d }, { 1, 1, sc, sd }, x)
        + G({ x / a }, 1) * G({ 0, b, 0, c, d }, { 1, sb, 1, sc, sd }, x)
        + G({ 0, 0, b, 0, c, d }, { 1, 1, sb, 1, sc, sd }, x) - sy[2] * Log(a, sa)
        + sy[2] * Log(-x, sa) };
    if (b != x) {
        res += (-sy[0] + G({ 0, x / a }, 1)) * G({ b, 0, c, d }, { sb, 1, sc, sd }, x);
    }
    if (c != x) {
        res += (-sy[1] + G({ 0, b / a, 0, x / a }, 1)) * G({ c, d }, { sc, sd }, x);
    }
    if (d != x) {
        res += (-sy[2] + G({ 0, b / a, 0, c / a, x / a }, 1)) * G({ d }, { sd }, x);
    }
    return res;
}

complex<double> G6_a0bc0d_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy
        = { Log(d, sd), G({ a }, { sa }, x), G({ 0, c / d, b / d }, 1),
              G({ a, 0, b }, { sa, 1, sb }, x), G({ 0, a }, { 1, sa }, x),
              G({ 0, 0, c / d, b / d }, 1), G({ 0, c / d, 0, b / d }, 1),
              G({ 0, c / d }, 1), G({ a, 0, b, c }, { sa, 1, sb, sc }, x),
              G({ x / d }, 1), G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x),
              G({ 0, c / d, b / d, 0, a / d }, 1),
              G({ a, 0, 0, b, c }, { sa, 1, 1, sb, sc }, x),
              G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x) };
    complex<double> res { sy[9] * sy[10] + sy[9] * (2. * sy[12] + sy[13])
        + sy[4] * (2. * sy[5] + sy[6])
        + sy[0] * (sy[11] - sy[2] * sy[4] + sy[1] * (2. * sy[5] + sy[6]) - sy[3] * sy[7])
        + 2. * sy[3] * G({ 0, 0, c / d }, 1) - 2. * G({ 0, 0, c / d, b / d, 0, a / d }, 1)
        - G({ 0, c / d, 0, b / d, 0, a / d }, 1)
        - 2. * G({ 0, c / d, b / d, 0, 0, a / d }, 1)
        - G({ 0, c / d, b / d, 0, a / d, x / d }, 1)
        + (-sy[10] - 2. * sy[12] - sy[13]) * G({ d }, { sd }, x)
        + sy[8] * (G({ 0, x / d }, 1) + G({ 0, d }, { 1, sd }, x))
        + sy[2] * (sy[3] - G({ 0, 0, a }, { 1, 1, sa }, x))
        + sy[7]
            * (-sy[8] - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                - 2. * G({ a, 0, 0, b }, { sa, 1, 1, sb }, x))
        + G({ 0, 0, a, 0, b, c }, { 1, 1, sa, 1, sb, sc }, x)
        + 2. * G({ 0, a, 0, 0, b, c }, { 1, sa, 1, 1, sb, sc }, x)
        + G({ 0, a, 0, b, 0, c }, { 1, sa, 1, sb, 1, sc }, x)
        + 3. * G({ a, 0, 0, 0, b, c }, { sa, 1, 1, 1, sb, sc }, x)
        + 2. * G({ a, 0, 0, b, 0, c }, { sa, 1, 1, sb, 1, sc }, x)
        + G({ a, 0, b, 0, 0, c }, { sa, 1, sb, 1, 1, sc }, x)
        + (-sy[11] + sy[0] * sy[1] * sy[2] + sy[2] * sy[4] + sy[1] * (-2. * sy[5] - sy[6])
              + sy[3] * sy[7])
            * Log(-x, sd)
        - (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.
        + sy[1]
            * (sy[11] - 3. * G({ 0, 0, 0, c / d, b / d }, 1)
                - 2. * G({ 0, 0, c / d, 0, b / d }, 1) - G({ 0, c / d, 0, 0, b / d }, 1)
                + sy[2] * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2))) };
    return res;
}
complex<double> G6_a0bc0d_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ a, 0, b }, { sa, 1, sb }, x),
        G({ 0, d / c, b / c }, 1), G({ 0, b / c, d / c }, 1),
        G({ 0, 0, a }, { 1, 1, sa }, x), G({ a, 0, d }, { sa, 1, sd }, x),
        G({ b / c, 0, a / c }, 1), G({ b / c, 0, d / c }, 1), Log(c, sc),
        G({ a }, { sa }, x), G({ a, d }, { sa, sd }, x), G({ 0, b / c, 0, a / c }, 1),
        G({ 0, a }, { 1, sa }, x), G({ 0, 0, b / c, d / c }, 1),
        G({ 0, 0, d / c, b / c }, 1), G({ 0, b / c, 0, d / c }, 1),
        G({ 0, d / c, 0, b / c }, 1), G({ 0, d / c }, 1),
        G({ a, 0, b, d }, { sa, 1, sb, sd }, x), G({ b / c, 0, 0, a / c }, 1),
        G({ 0, b / c, 0, a / c, d / c }, 1), G({ 0, b / c, 0, d / c, a / c }, 1),
        G({ 0, b / c, d / c, 0, a / c }, 1), G({ 0, d / c, b / c, 0, a / c }, 1),
        G({ b / c, 0, 0, a / c, d / c }, 1), G({ b / c, 0, 0, d / c, a / c }, 1),
        G({ b / c, 0, a / c, 0, d / c }, 1), G({ b / c, 0, d / c, 0, a / c }, 1),
        G({ 0, 0 }, { 1, 1 }, x) };
    complex<double> res { sy[9] * sy[10]
        + sy[11] * (2. * sy[12] + 2. * sy[13] + sy[14] + sy[15]) + sy[0] * sy[1]
        + sy[2] * (-sy[3] + sy[4]) + sy[3] * (-sy[1] - sy[6]) + sy[4] * (sy[5] + sy[6])
        + sy[7]
            * (-(sy[0] * sy[16]) + sy[19] + sy[20] + sy[21] + sy[22] + 2. * sy[23]
                + 2. * sy[24] + sy[25] + sy[26] - sy[11] * sy[2]
                + sy[11] * (-sy[1] - sy[6])
                + (2. * sy[12] + 2. * sy[13] + sy[14] + sy[15]) * sy[8])
        + 2. * sy[0] * G({ 0, 0, d / c }, 1)
        + sy[9] * (-sy[14] + 2. * sy[18] - 2. * G({ b / c, 0, 0, d / c }, 1))
        - 2. * G({ 0, 0, b / c, 0, a / c, d / c }, 1)
        - 2. * G({ 0, 0, b / c, 0, d / c, a / c }, 1)
        - 2. * G({ 0, 0, b / c, d / c, 0, a / c }, 1)
        - 2. * G({ 0, 0, d / c, b / c, 0, a / c }, 1)
        - 4. * G({ 0, b / c, 0, 0, a / c, d / c }, 1)
        - 4. * G({ 0, b / c, 0, 0, d / c, a / c }, 1)
        - 2. * G({ 0, b / c, 0, a / c, 0, d / c }, 1)
        - G({ 0, b / c, 0, a / c, d / c, x / c }, 1)
        - 3. * G({ 0, b / c, 0, d / c, 0, a / c }, 1)
        - G({ 0, b / c, 0, d / c, a / c, x / c }, 1)
        - 2. * G({ 0, b / c, d / c, 0, 0, a / c }, 1)
        - G({ 0, b / c, d / c, 0, a / c, x / c }, 1)
        - G({ 0, d / c, 0, b / c, 0, a / c }, 1)
        - 2. * G({ 0, d / c, b / c, 0, 0, a / c }, 1)
        - G({ 0, d / c, b / c, 0, a / c, x / c }, 1)
        - 6. * G({ b / c, 0, 0, 0, a / c, d / c }, 1)
        - 6. * G({ b / c, 0, 0, 0, d / c, a / c }, 1)
        - 4. * G({ b / c, 0, 0, a / c, 0, d / c }, 1)
        - 2. * G({ b / c, 0, 0, a / c, d / c, x / c }, 1)
        - 4. * G({ b / c, 0, 0, d / c, 0, a / c }, 1)
        - 2. * G({ b / c, 0, 0, d / c, a / c, x / c }, 1)
        - 2. * G({ b / c, 0, a / c, 0, 0, d / c }, 1)
        - G({ b / c, 0, a / c, 0, d / c, x / c }, 1)
        - 2. * G({ b / c, 0, d / c, 0, 0, a / c }, 1)
        - G({ b / c, 0, d / c, 0, a / c, x / c }, 1)
        + (-sy[10] - 2. * sy[18] - G({ b / c, 0, a / c, x / c }, 1))
            * G({ 0, d }, { 1, sd }, x)
        - sy[5] * G({ 0, 0, d }, { 1, 1, sd }, x)
        + sy[16]
            * (-sy[17] - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                - 2. * G({ a, 0, 0, b }, { sa, 1, 1, sb }, x))
        + G({ 0, b / c }, 1) * (sy[17] - G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
        + G({ b / c }, 1)
            * (-G({ a, 0, 0, 0, d }, { sa, 1, 1, 1, sd }, x)
                + G({ a, 0, b, 0, d }, { sa, 1, sb, 1, sd }, x))
        + G({ a, 0, b, 0, 0, d }, { sa, 1, sb, 1, 1, sd }, x)
        + (sy[0] * sy[16] - sy[19] - sy[20] - sy[21] - sy[22] - 2. * sy[23] - 2. * sy[24]
              - sy[25] - sy[26] + sy[11] * sy[2] + sy[11] * (sy[1] + sy[6])
              + (-2. * sy[12] - 2. * sy[13] - sy[14] - sy[15]) * sy[8]
              + (sy[1] + sy[2] + sy[6]) * sy[7] * sy[8])
            * Log(-x, sc)
        + (-sy[1] / 2. - sy[2] / 2. - sy[6] / 2.) * sy[8] * pow(sy[7], 2.)
        + sy[8]
            * (sy[20] + sy[21] + sy[22] + 2. * sy[24] + sy[26] + sy[27] * (sy[1] + sy[6])
                - 3. * G({ 0, 0, 0, b / c, d / c }, 1)
                - 3. * G({ 0, 0, 0, d / c, b / c }, 1) - G({ 0, 0, b / c, 0, d / c }, 1)
                - 2. * G({ 0, 0, d / c, 0, b / c }, 1) + G({ 0, b / c, 0, 0, d / c }, 1)
                - G({ 0, d / c, 0, 0, b / c }, 1) + 3. * G({ b / c, 0, 0, 0, d / c }, 1)
                + 2. * sy[1] * Zeta(2) + 2. * sy[6] * Zeta(2)
                + sy[2] * (sy[27] + 2. * Zeta(2))) };
    if (d != x) {
        res += (sy[19] + 2. * sy[23] + sy[25] - G({ 0, b / c, 0, a / c, x / c }, 1)
                   - 2. * G({ b / c, 0, 0, a / c, x / c }, 1)
                   - G({ b / c, 0, a / c, 0, x / c }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}
complex<double> G6_a0bc0d_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (b == c) {
        // abbd
        const vector<complex<double>> sy = { Log(b, sb), G({ a }, { sa }, x),
            G({ 0, 1, d / b }, 1), G({ 0, d / b, 1 }, 1), G({ 0, 0, a }, { 1, 1, sa }, x),
            G({ 0, a / b, 1 }, 1), G({ a, 0, d }, { sa, 1, sd }, x),
            G({ b, 0, d }, { sb, 1, sd }, x), G({ a, d }, { sa, sd }, x),
            G({ 0, 0, 1, a / b }, 1), G({ 0, 0, a / b, 1 }, 1), G({ 0, 1, 0, a / b }, 1),
            G({ a, b, 0, d }, { sa, sb, 1, sd }, x), G({ 0, a }, { 1, sa }, x),
            G({ 0, 0, 1, a / b, d / b }, 1), G({ 0, 0, 1, d / b, a / b }, 1),
            G({ 0, 0, a / b, 1, d / b }, 1), G({ 0, 0, a / b, d / b, 1 }, 1),
            G({ 0, 0, d / b, 1, a / b }, 1), G({ 0, 0, d / b, a / b, 1 }, 1),
            G({ 0, 1, 0, a / b, d / b }, 1), G({ 0, 1, 0, d / b, a / b }, 1),
            G({ 0, 1, d / b, 0, a / b }, 1), G({ 0, a / b, 0, 1, d / b }, 1),
            G({ 0, a / b, 0, d / b, 1 }, 1), G({ 0, d / b, 0, 1, a / b }, 1),
            G({ 0, d / b, 0, a / b, 1 }, 1), G({ 0, d / b, 1, 0, a / b }, 1),
            G({ 0, 0 }, { 1, 1 }, x) };
        complex<double> res { sy[0]
                * (2. * sy[14] + 2. * sy[15] + 2. * sy[16] + 2. * sy[17] + 2. * sy[18]
                    + 2. * sy[19] + sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25]
                    + sy[26] + sy[27] - sy[13] * sy[2] - sy[13] * sy[3])
            - sy[3] * sy[4] + sy[3] * sy[6] + sy[5] * sy[6] + sy[2] * (-sy[4] + sy[6])
            - sy[5] * sy[7] + 2. * sy[9] * sy[8] + sy[7] * G({ 0, a / b, x / b }, 1)
            + sy[8]
                * (2. * sy[10] + sy[11] - 2. * G({ 0, 0, 1, d / b }, 1)
                    - 2. * G({ 0, 0, d / b, 1 }, 1) - G({ 0, 1, 0, d / b }, 1))
            - 6. * G({ 0, 0, 0, 1, a / b, d / b }, 1)
            - 6. * G({ 0, 0, 0, 1, d / b, a / b }, 1)
            - 6. * G({ 0, 0, 0, a / b, 1, d / b }, 1)
            - 6. * G({ 0, 0, 0, a / b, d / b, 1 }, 1)
            - 6. * G({ 0, 0, 0, d / b, 1, a / b }, 1)
            - 6. * G({ 0, 0, 0, d / b, a / b, 1 }, 1)
            - 4. * G({ 0, 0, 1, 0, a / b, d / b }, 1)
            - 4. * G({ 0, 0, 1, 0, d / b, a / b }, 1)
            - 2. * G({ 0, 0, 1, a / b, 0, d / b }, 1)
            - 2. * G({ 0, 0, 1, a / b, d / b, x / b }, 1)
            - 4. * G({ 0, 0, 1, d / b, 0, a / b }, 1)
            - 2. * G({ 0, 0, 1, d / b, a / b, x / b }, 1)
            - 4. * G({ 0, 0, a / b, 0, 1, d / b }, 1)
            - 4. * G({ 0, 0, a / b, 0, d / b, 1 }, 1)
            - 2. * G({ 0, 0, a / b, 1, 0, d / b }, 1)
            - 2. * G({ 0, 0, a / b, 1, d / b, x / b }, 1)
            - 2. * G({ 0, 0, a / b, d / b, 1, x / b }, 1)
            - 2. * G({ 0, 0, a / b, d / b, x / b, 1 }, 1)
            - 4. * G({ 0, 0, d / b, 0, 1, a / b }, 1)
            - 4. * G({ 0, 0, d / b, 0, a / b, 1 }, 1)
            - 4. * G({ 0, 0, d / b, 1, 0, a / b }, 1)
            - 2. * G({ 0, 0, d / b, 1, a / b, x / b }, 1)
            - 2. * G({ 0, 0, d / b, a / b, 1, x / b }, 1)
            - 2. * G({ 0, 0, d / b, a / b, x / b, 1 }, 1)
            - 2. * G({ 0, 1, 0, 0, a / b, d / b }, 1)
            - 2. * G({ 0, 1, 0, 0, d / b, a / b }, 1) - G({ 0, 1, 0, a / b, 0, d / b }, 1)
            - G({ 0, 1, 0, a / b, d / b, x / b }, 1)
            - 2. * G({ 0, 1, 0, d / b, 0, a / b }, 1)
            - G({ 0, 1, 0, d / b, a / b, x / b }, 1)
            - 2. * G({ 0, 1, d / b, 0, 0, a / b }, 1)
            - G({ 0, 1, d / b, 0, a / b, x / b }, 1)
            - 2. * G({ 0, a / b, 0, 0, 1, d / b }, 1)
            - 2. * G({ 0, a / b, 0, 0, d / b, 1 }, 1) - G({ 0, a / b, 0, 1, 0, d / b }, 1)
            - G({ 0, a / b, 0, 1, d / b, x / b }, 1)
            - G({ 0, a / b, 0, d / b, 1, x / b }, 1)
            - G({ 0, a / b, 0, d / b, x / b, 1 }, 1)
            - 2. * G({ 0, d / b, 0, 0, 1, a / b }, 1)
            - 2. * G({ 0, d / b, 0, 0, a / b, 1 }, 1)
            - 2. * G({ 0, d / b, 0, 1, 0, a / b }, 1)
            - G({ 0, d / b, 0, 1, a / b, x / b }, 1)
            - G({ 0, d / b, 0, a / b, 1, x / b }, 1)
            - G({ 0, d / b, 0, a / b, x / b, 1 }, 1)
            - 2. * G({ 0, d / b, 1, 0, 0, a / b }, 1)
            - G({ 0, d / b, 1, 0, a / b, x / b }, 1)
            + (-2. * sy[9] - 2. * sy[10] - sy[11] - G({ 0, a / b, x / b, 1 }, 1))
                * G({ 0, d }, { 1, sd }, x)
            + G({ 0, a / b }, 1) * (-sy[12] + G({ 0, b, 0, d }, { 1, sb, 1, sd }, x))
            + G({ a, 0, 0, b, 0, d }, { sa, 1, 1, sb, 1, sd }, x)
            + (-2. * sy[14] - 2. * sy[15] - 2. * sy[16] - 2. * sy[17] - 2. * sy[18]
                  - 2. * sy[19] - sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25]
                  - sy[26] - sy[27] + sy[13] * sy[2] + sy[13] * sy[3]
                  + sy[0] * sy[1] * (sy[2] + sy[3]))
                * Log(-x, sb)
            + sy[1] * (-sy[2] / 2. - sy[3] / 2.) * pow(sy[0], 2.) - sy[12] * Zeta(2)
            + G({ a, 0, 0, d }, { sa, 1, 1, sd }, x) * Zeta(2)
            + sy[1]
                * (2. * sy[15] + 2. * sy[18] + 2. * sy[19] + sy[21] + sy[22] + sy[25]
                    + sy[26] + sy[27] + sy[28] * sy[3] + 3. * G({ 0, 0, 0, 1, d / b }, 1)
                    + 3. * G({ 0, 0, 0, d / b, 1 }, 1) + 2. * G({ 0, 0, 1, 0, d / b }, 1)
                    + G({ 0, 1, 0, 0, d / b }, 1) + 2. * sy[3] * Zeta(2)
                    + sy[2] * (sy[28] + 2. * Zeta(2))) };
        if (d != x) {
            res += (2. * sy[14] + 2. * sy[16] + 2. * sy[17] + sy[20] + sy[23] + sy[24]
                       - 2. * G({ 0, 0, 1, a / b, x / b }, 1)
                       - 2. * G({ 0, 0, a / b, 1, x / b }, 1)
                       - 2. * G({ 0, 0, a / b, x / b, 1 }, 1)
                       - G({ 0, 1, 0, a / b, x / b }, 1) - G({ 0, a / b, 0, 1, x / b }, 1)
                       - G({ 0, a / b, 0, x / b, 1 }, 1))
                * G({ d }, { sd }, x);
        }
        return res;

    } else { // abcd
        const vector<complex<double>> sy = { G({ 0, c / b, a / b }, 1),
            G({ a, 0, d }, { sa, 1, sd }, x), G({ 0, a / b, c / b }, 1),
            G({ c / b, 0, a / b }, 1), G({ c / b, 0, d / b }, 1), Log(b, sb),
            G({ a }, { sa }, x), G({ a, d }, { sa, sd }, x), G({ 0, c / b, 0, a / b }, 1),
            G({ a, c, 0, d }, { sa, sc, 1, sd }, x), G({ c / b, 0, 0, a / b }, 1),
            G({ 0, a }, { 1, sa }, x), G({ 0, a / b, c / b, 0, d / b }, 1),
            G({ 0, c / b, 0, a / b, d / b }, 1), G({ 0, c / b, 0, d / b, a / b }, 1),
            G({ 0, c / b, a / b, 0, d / b }, 1), G({ c / b, 0, 0, a / b, d / b }, 1),
            G({ c / b, 0, 0, d / b, a / b }, 1), G({ c / b, 0, a / b, 0, d / b }, 1),
            G({ c / b, 0, d / b, 0, a / b }, 1) };
        complex<double> res { -(sy[0] * sy[1]) + sy[1] * (-sy[3] - sy[4])
            + (-sy[12] - sy[13] - sy[14] - sy[15] - 2. * sy[16] - 2. * sy[17] - sy[18]
                  - sy[19] + sy[11] * sy[4])
                * sy[5]
            - sy[7] * sy[8]
            + sy[7]
                * (-2. * sy[10] + G({ 0, c / b, 0, d / b }, 1)
                    + 2. * G({ c / b, 0, 0, d / b }, 1))
            + 2. * G({ 0, 0, a / b, c / b, 0, d / b }, 1)
            + 2. * G({ 0, 0, c / b, 0, a / b, d / b }, 1)
            + 2. * G({ 0, 0, c / b, 0, d / b, a / b }, 1)
            + 2. * G({ 0, 0, c / b, a / b, 0, d / b }, 1)
            + G({ 0, a / b, 0, c / b, 0, d / b }, 1)
            + 2. * G({ 0, a / b, c / b, 0, 0, d / b }, 1)
            + G({ 0, a / b, c / b, 0, d / b, x / b }, 1)
            + 4. * G({ 0, c / b, 0, 0, a / b, d / b }, 1)
            + 4. * G({ 0, c / b, 0, 0, d / b, a / b }, 1)
            + 3. * G({ 0, c / b, 0, a / b, 0, d / b }, 1)
            + G({ 0, c / b, 0, a / b, d / b, x / b }, 1)
            + 2. * G({ 0, c / b, 0, d / b, 0, a / b }, 1)
            + G({ 0, c / b, 0, d / b, a / b, x / b }, 1)
            + 2. * G({ 0, c / b, a / b, 0, 0, d / b }, 1)
            + G({ 0, c / b, a / b, 0, d / b, x / b }, 1)
            + 6. * G({ c / b, 0, 0, 0, a / b, d / b }, 1)
            + 6. * G({ c / b, 0, 0, 0, d / b, a / b }, 1)
            + 4. * G({ c / b, 0, 0, a / b, 0, d / b }, 1)
            + 2. * G({ c / b, 0, 0, a / b, d / b, x / b }, 1)
            + 4. * G({ c / b, 0, 0, d / b, 0, a / b }, 1)
            + 2. * G({ c / b, 0, 0, d / b, a / b, x / b }, 1)
            + 2. * G({ c / b, 0, a / b, 0, 0, d / b }, 1)
            + G({ c / b, 0, a / b, 0, d / b, x / b }, 1)
            + 2. * G({ c / b, 0, d / b, 0, 0, a / b }, 1)
            + G({ c / b, 0, d / b, 0, a / b, x / b }, 1)
            + (2. * sy[10] + sy[8] + G({ 0, a / b, c / b, x / b }, 1)
                  + G({ 0, c / b, a / b, x / b }, 1) + G({ c / b, 0, a / b, x / b }, 1))
                * G({ 0, d }, { 1, sd }, x)
            + sy[4] * G({ 0, 0, a }, { 1, 1, sa }, x)
            + (sy[0] + sy[2] + sy[3]) * G({ 0, 0, d }, { 1, 1, sd }, x)
            + G({ 0, a / b }, 1) * (-sy[9] + G({ 0, c, 0, d }, { 1, sc, 1, sd }, x))
            + G({ 0, c / b }, 1) * (sy[9] - G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
            + G({ c / b }, 1)
                * (G({ a, 0, 0, 0, d }, { sa, 1, 1, 1, sd }, x)
                    - G({ a, 0, c, 0, d }, { sa, 1, sc, 1, sd }, x))
            + G({ a, 0, 0, c, 0, d }, { sa, 1, 1, sc, 1, sd }, x)
            + (sy[12] + sy[13] + sy[14] + sy[15] + 2. * sy[16] + 2. * sy[17] + sy[18]
                  + sy[19] - sy[11] * sy[4] - sy[4] * sy[5] * sy[6])
                * Log(-x, sb)
            + (sy[4] * sy[6] * pow(sy[5], 2.)) / 2.
            + sy[6]
                * (-sy[14] - 2. * sy[17] - sy[19] - G({ 0, 0, c / b, 0, d / b }, 1)
                    - 2. * G({ 0, c / b, 0, 0, d / b }, 1)
                    - 3. * G({ c / b, 0, 0, 0, d / b }, 1)
                    + sy[4] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
        if (c != x) {
            res += (-sy[2] + G({ 0, a / b, x / b }, 1))
                * G({ c, 0, d }, { sc, 1, sd }, x);
        }
        if (d != x) {
            res += (-sy[12] - sy[13] - sy[15] - 2. * sy[16] - sy[18]
                       + G({ 0, a / b, c / b, 0, x / b }, 1)
                       + G({ 0, c / b, 0, a / b, x / b }, 1)
                       + G({ 0, c / b, a / b, 0, x / b }, 1)
                       + 2. * G({ c / b, 0, 0, a / b, x / b }, 1)
                       + G({ c / b, 0, a / b, 0, x / b }, 1))
                * G({ d }, { sd }, x);
        }
        return res;
    }
}
complex<double> G6_a0bc0d_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy = { G({ 0, b / a, c / a }, 1), G({ 0, b / a }, 1),
        G({ 0, b / a, c / a, 0, d / a }, 1) };
    complex<double> res { 2. * G({ 0, 0, b / a, c / a, 0, d / a }, 1)
        + G({ 0, b / a, 0, c / a, 0, d / a }, 1)
        + 2. * G({ 0, b / a, c / a, 0, 0, d / a }, 1)
        + G({ 0, b / a, c / a, 0, d / a, x / a }, 1)
        + G({ 0, b / a, c / a, x / a }, 1) * G({ 0, d }, { 1, sd }, x)
        + sy[0] * G({ 0, 0, d }, { 1, 1, sd }, x)
        + sy[1] * G({ 0, c, 0, d }, { 1, sc, 1, sd }, x)
        + G({ x / a }, 1) * G({ 0, b, c, 0, d }, { 1, sb, sc, 1, sd }, x)
        + G({ 0, 0, b, c, 0, d }, { 1, 1, sb, sc, 1, sd }, x) - sy[2] * Log(a, sa)
        + sy[2] * Log(-x, sa) };
    if (b != x) {
        res += (-sy[1] + G({ 0, x / a }, 1)) * G({ b, c, 0, d }, { sb, sc, 1, sd }, x);
    }
    if (c != x) {
        res += (-sy[0] + G({ 0, b / a, x / a }, 1)) * G({ c, 0, d }, { sc, 1, sd }, x);
    }
    if (d != x) {
        res += (-sy[2] + G({ 0, b / a, c / a, 0, x / a }, 1)) * G({ d }, { sd }, x);
    }
    return res;
}

complex<double> G6_ab0c0d_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { Log(d, sd), G({ 0, c / d }, 1),
        G({ a, b }, { sa, sb }, x), G({ 0, 0, c / d }, 1),
        G({ 0, a, b }, { 1, sa, sb }, x), G({ a, 0, b }, { sa, 1, sb }, x),
        G({ 0, c / d, 0, b / d }, 1), G({ a, b, 0, c }, { sa, sb, 1, sc }, x),
        G({ x / d }, 1), G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x),
        G({ a }, { sa }, x), G({ 0, c / d, 0, b / d, a / d }, 1),
        G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x),
        G({ a, b, 0, 0, c }, { sa, sb, 1, 1, sc }, x) };
    complex<double> res { sy[3] * (-2. * sy[4] - 2. * sy[5])
        + sy[0] * (sy[11] - 2. * sy[2] * sy[3] + sy[1] * (sy[4] + sy[5]) - sy[10] * sy[6])
        + sy[9] * sy[8] + (sy[12] + 2. * sy[13]) * sy[8]
        + sy[2] * (-sy[6] + 3. * G({ 0, 0, 0, c / d }, 1))
        + sy[10]
            * (sy[11] + 2. * G({ 0, 0, c / d, 0, b / d }, 1)
                + 2. * G({ 0, c / d, 0, 0, b / d }, 1))
        - 2. * G({ 0, 0, c / d, 0, b / d, a / d }, 1)
        - 2. * G({ 0, c / d, 0, 0, b / d, a / d }, 1)
        - G({ 0, c / d, 0, b / d, 0, a / d }, 1)
        - G({ 0, c / d, 0, b / d, a / d, x / d }, 1)
        + (-sy[9] - sy[12] - 2. * sy[13]) * G({ d }, { sd }, x)
        - sy[6] * G({ 0, a }, { 1, sa }, x)
        + sy[7] * (G({ 0, x / d }, 1) + G({ 0, d }, { 1, sd }, x))
        + G({ 0, 0, a, b, 0, c }, { 1, 1, sa, sb, 1, sc }, x)
        + G({ 0, a, 0, b, 0, c }, { 1, sa, 1, sb, 1, sc }, x)
        + 2. * G({ 0, a, b, 0, 0, c }, { 1, sa, sb, 1, 1, sc }, x)
        + G({ a, 0, 0, b, 0, c }, { sa, 1, 1, sb, 1, sc }, x)
        + 2. * G({ a, 0, b, 0, 0, c }, { sa, 1, sb, 1, 1, sc }, x)
        + 3. * G({ a, b, 0, 0, 0, c }, { sa, sb, 1, 1, 1, sc }, x)
        + (-sy[11] - sy[0] * sy[1] * sy[2] + 2. * sy[2] * sy[3] + sy[1] * (-sy[4] - sy[5])
              + sy[10] * sy[6])
            * Log(-x, sd)
        + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.
        + sy[1]
            * (-sy[7] + G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                + G({ a, 0, 0, b }, { sa, 1, 1, sb }, x)
                + sy[2] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
    return res;
}
complex<double> G6_ab0c0d_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy
        = { Log(c, sc), G({ 0, d / c }, 1), G({ a, b }, { sa, sb }, x),
              G({ 0, b / c, a / c }, 1), G({ 0, a, b }, { 1, sa, sb }, x),
              G({ 0, d / c, x / c }, 1), G({ a, 0, b }, { sa, 1, sb }, x),
              G({ a, 0, d }, { sa, 1, sd }, x), G({ a, b, d }, { sa, sb, sd }, x),
              G({ a, d }, { sa, sd }, x), G({ 0, 0, b / c, a / c }, 1),
              G({ 0, a }, { 1, sa }, x), G({ 0, 0, b / c, d / c }, 1),
              G({ 0, b / c, 0, a / c }, 1), G({ 0, b / c, 0, d / c }, 1),
              G({ 0, 0, d / c, b / c }, 1), G({ 0, d / c, 0, b / c }, 1),
              G({ 0, d / c, 0, x / c }, 1), G({ a, b, 0, d }, { sa, sb, 1, sd }, x),
              G({ a }, { sa }, x), G({ 0, 0, b / c, a / c, d / c }, 1),
              G({ 0, 0, b / c, d / c, a / c }, 1), G({ 0, 0, d / c, b / c, a / c }, 1),
              G({ 0, b / c, 0, a / c, d / c }, 1), G({ 0, b / c, 0, d / c, a / c }, 1),
              G({ 0, b / c, a / c, 0, d / c }, 1), G({ 0, d / c, 0, b / c, a / c }, 1) };
    complex<double> res { 2. * sy[9] * sy[10] - 2. * sy[11] * sy[12]
        + sy[9] * (-2. * sy[12] + sy[13] - sy[14])
        + sy[11] * (-sy[14] - 2. * sy[15] - sy[16]) - sy[4] * sy[5]
        + sy[0]
            * ((-2. * sy[12] - sy[14] - 2. * sy[15] - sy[16]) * sy[19] + 2. * sy[20]
                + 2. * sy[21] + 2. * sy[22] + sy[23] + sy[24] + sy[25] + sy[26]
                + sy[1] * (sy[4] + sy[6]))
        + sy[3] * sy[7] + (2. * sy[7] - 2. * sy[8]) * G({ 0, 0, b / c }, 1)
        + 2. * sy[8] * G({ 0, 0, d / c }, 1)
        + sy[2] * (-2. * sy[15] - sy[17] - 3. * G({ 0, 0, 0, d / c }, 1))
        - 6. * G({ 0, 0, 0, b / c, a / c, d / c }, 1)
        - 6. * G({ 0, 0, 0, b / c, d / c, a / c }, 1)
        - 6. * G({ 0, 0, 0, d / c, b / c, a / c }, 1)
        - 4. * G({ 0, 0, b / c, 0, a / c, d / c }, 1)
        - 4. * G({ 0, 0, b / c, 0, d / c, a / c }, 1)
        - 4. * G({ 0, 0, b / c, a / c, 0, d / c }, 1)
        - 2. * G({ 0, 0, b / c, a / c, d / c, x / c }, 1)
        - 2. * G({ 0, 0, b / c, d / c, 0, a / c }, 1)
        - 2. * G({ 0, 0, b / c, d / c, a / c, x / c }, 1)
        - 4. * G({ 0, 0, d / c, 0, b / c, a / c }, 1)
        - 2. * G({ 0, 0, d / c, b / c, 0, a / c }, 1)
        - 2. * G({ 0, 0, d / c, b / c, a / c, x / c }, 1)
        - 2. * G({ 0, b / c, 0, 0, a / c, d / c }, 1)
        - 2. * G({ 0, b / c, 0, 0, d / c, a / c }, 1)
        - 2. * G({ 0, b / c, 0, a / c, 0, d / c }, 1)
        - G({ 0, b / c, 0, a / c, d / c, x / c }, 1)
        - G({ 0, b / c, 0, d / c, 0, a / c }, 1)
        - G({ 0, b / c, 0, d / c, a / c, x / c }, 1)
        - 2. * G({ 0, b / c, a / c, 0, 0, d / c }, 1)
        - G({ 0, b / c, a / c, 0, d / c, x / c }, 1)
        - 2. * G({ 0, d / c, 0, 0, b / c, a / c }, 1)
        - G({ 0, d / c, 0, b / c, 0, a / c }, 1)
        - G({ 0, d / c, 0, b / c, a / c, x / c }, 1)
        + sy[19]
            * (2. * sy[21] + 2. * sy[22] + sy[24] + sy[26]
                + 6. * G({ 0, 0, 0, b / c, d / c }, 1)
                + 6. * G({ 0, 0, 0, d / c, b / c }, 1)
                + 4. * G({ 0, 0, b / c, 0, d / c }, 1)
                + 4. * G({ 0, 0, d / c, 0, b / c }, 1)
                + 2. * G({ 0, b / c, 0, 0, d / c }, 1)
                + 2. * G({ 0, d / c, 0, 0, b / c }, 1)
                + sy[5] * G({ 0, b }, { 1, sb }, x))
        + (-2. * sy[10] - sy[13] - G({ 0, b / c, a / c, x / c }, 1))
            * G({ 0, d }, { 1, sd }, x)
        - sy[3] * G({ 0, 0, d }, { 1, 1, sd }, x)
        + sy[5] * (-sy[6] - G({ 0, b, a }, { 1, sb, sa }, x))
        + G({ 0, b / c }, 1) * (-sy[18] + G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
        + G({ a, b, 0, 0, 0, d }, { sa, sb, 1, 1, 1, sd }, x)
        + ((2. * sy[12] + sy[14] + 2. * sy[15] + sy[16]) * sy[19] - 2. * sy[20]
              - 2. * sy[21] - 2. * sy[22] - sy[23] - sy[24] - sy[25] - sy[26]
              - sy[0] * sy[1] * sy[2] + sy[1] * (-sy[4] - sy[6]))
            * Log(-x, sc)
        + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.
        + sy[1]
            * (-sy[18] + sy[19] * G({ 0, 0, b }, { 1, 1, sb }, x)
                - G({ 0, 0, b, a }, { 1, 1, sb, sa }, x)
                + sy[2] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
    if (b != x) {
        res += (-(sy[16] * sy[19]) + sy[17] * sy[19]) * G({ b }, { sb }, x)
            + (sy[16] - sy[17]) * G({ b, a }, { sb, sa }, x);
    }
    if (d != x) {
        res += (2. * sy[20] + sy[23] + sy[25] - 2. * G({ 0, 0, b / c, a / c, x / c }, 1)
                   - G({ 0, b / c, 0, a / c, x / c }, 1)
                   - G({ 0, b / c, a / c, 0, x / c }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}
complex<double> G6_ab0c0d_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    // abcd
    const vector<complex<double>> sy
        = { G({ 0, c / b, a / b }, 1), G({ 0, a / b, c / b }, 1),
              G({ a / b, 0, c / b }, 1), G({ 0, c, 0, d }, { 1, sc, 1, sd }, x),
              G({ a, d }, { sa, sd }, x), G({ 0, c / b, 0, a / b }, 1),
              G({ 0, c / b, 0, d / b }, 1), G({ a, c, 0, d }, { sa, sc, 1, sd }, x),
              G({ a }, { sa }, x), G({ 0, c / b, 0, d / b, a / b }, 1),
              G({ 0, a / b, c / b, 0, d / b }, 1), G({ 0, c / b, 0, a / b, d / b }, 1),
              G({ 0, c / b, a / b, 0, d / b }, 1), G({ a / b, 0, c / b, 0, d / b }, 1) };
    complex<double> res { sy[4] * sy[5] - sy[4] * sy[6]
        + (-sy[3] + sy[7]) * G({ 0, a / b }, 1) - sy[3] * G({ a / b, x / b }, 1)
        + sy[8]
            * (sy[9] + 2. * G({ 0, 0, c / b, 0, d / b }, 1)
                + 2. * G({ 0, c / b, 0, 0, d / b }, 1))
        - 2. * G({ 0, 0, a / b, c / b, 0, d / b }, 1)
        - 2. * G({ 0, 0, c / b, 0, a / b, d / b }, 1)
        - 2. * G({ 0, 0, c / b, 0, d / b, a / b }, 1)
        - 2. * G({ 0, 0, c / b, a / b, 0, d / b }, 1)
        - 2. * G({ 0, a / b, 0, c / b, 0, d / b }, 1)
        - 2. * G({ 0, a / b, c / b, 0, 0, d / b }, 1)
        - G({ 0, a / b, c / b, 0, d / b, x / b }, 1)
        - 2. * G({ 0, c / b, 0, 0, a / b, d / b }, 1)
        - 2. * G({ 0, c / b, 0, 0, d / b, a / b }, 1)
        - 2. * G({ 0, c / b, 0, a / b, 0, d / b }, 1)
        - G({ 0, c / b, 0, a / b, d / b, x / b }, 1)
        - G({ 0, c / b, 0, d / b, 0, a / b }, 1)
        - G({ 0, c / b, 0, d / b, a / b, x / b }, 1)
        - 2. * G({ 0, c / b, a / b, 0, 0, d / b }, 1)
        - G({ 0, c / b, a / b, 0, d / b, x / b }, 1)
        - 2. * G({ a / b, 0, 0, c / b, 0, d / b }, 1)
        - 2. * G({ a / b, 0, c / b, 0, 0, d / b }, 1)
        - G({ a / b, 0, c / b, 0, d / b, x / b }, 1) - sy[6] * G({ 0, a }, { 1, sa }, x)
        + (-sy[5] - G({ 0, a / b, c / b, x / b }, 1) - G({ 0, c / b, a / b, x / b }, 1)
              - G({ a / b, 0, c / b, x / b }, 1))
            * G({ 0, d }, { 1, sd }, x)
        + (-sy[0] - sy[1] - sy[2]) * G({ 0, 0, d }, { 1, 1, sd }, x)
        + sy[0] * G({ a, 0, d }, { sa, 1, sd }, x)
        + G({ 0, c / b }, 1) * (-sy[7] + G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
        + G({ a / b }, 1)
            * (-G({ 0, 0, c, 0, d }, { 1, 1, sc, 1, sd }, x)
                + G({ a, 0, c, 0, d }, { sa, 1, sc, 1, sd }, x))
        + G({ a, 0, 0, c, 0, d }, { sa, 1, 1, sc, 1, sd }, x)
        + (sy[9] + sy[10] + sy[11] + sy[12] + sy[13] - sy[6] * sy[8]) * Log(b, sb)
        + (-sy[9] - sy[10] - sy[11] - sy[12] - sy[13] + sy[6] * sy[8]) * Log(-x, sb) };
    if (c != x) {
        res += (sy[1] + sy[2] - G({ 0, a / b, x / b }, 1) - G({ a / b, 0, x / b }, 1))
            * G({ c, 0, d }, { sc, 1, sd }, x);
    }
    if (d != x) {
        res += (sy[10] + sy[11] + sy[12] + sy[13] - G({ 0, a / b, c / b, 0, x / b }, 1)
                   - G({ 0, c / b, 0, a / b, x / b }, 1)
                   - G({ 0, c / b, a / b, 0, x / b }, 1)
                   - G({ a / b, 0, c / b, 0, x / b }, 1))
            * G({ d }, { sd }, x);
    }
    return res;
}
complex<double> G6_ab0c0d_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, int sa, int sb, int sc, int sd, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    double x = 1.;
    if (a == b) {
        // aacd

        const vector<complex<double>> sy = { G({ 0, 1, c / a }, 1), G({ 0, c / a, 1 }, 1),
            G({ 0, 1, c / a, 0, d / a }, 1), G({ 0, c / a, 0, 1, d / a }, 1),
            G({ 0, c / a, 0, d / a, 1 }, 1), G({ 0, c / a, 1, 0, d / a }, 1) };
        complex<double> res { -2. * G({ 0, 0, 1, c / a, 0, d / a }, 1)
            - 2. * G({ 0, 0, c / a, 0, 1, d / a }, 1)
            - 2. * G({ 0, 0, c / a, 0, d / a, 1 }, 1)
            - 2. * G({ 0, 0, c / a, 1, 0, d / a }, 1) - G({ 0, 1, 0, c / a, 0, d / a }, 1)
            - 2. * G({ 0, 1, c / a, 0, 0, d / a }, 1)
            - G({ 0, 1, c / a, 0, d / a, x / a }, 1)
            - 2. * G({ 0, c / a, 0, 0, 1, d / a }, 1)
            - 2. * G({ 0, c / a, 0, 0, d / a, 1 }, 1)
            - 2. * G({ 0, c / a, 0, 1, 0, d / a }, 1)
            - G({ 0, c / a, 0, 1, d / a, x / a }, 1)
            - G({ 0, c / a, 0, d / a, 1, x / a }, 1)
            - G({ 0, c / a, 0, d / a, x / a, 1 }, 1)
            - 2. * G({ 0, c / a, 1, 0, 0, d / a }, 1)
            - G({ 0, c / a, 1, 0, d / a, x / a }, 1)
            + (-G({ 0, 1, c / a, x / a }, 1) - G({ 0, c / a, 1, x / a }, 1)
                  - G({ 0, c / a, x / a, 1 }, 1))
                * G({ 0, d }, { 1, sd }, x)
            + (-sy[0] - sy[1]) * G({ 0, 0, d }, { 1, 1, sd }, x)
            - G({ x / a, 1 }, 1) * G({ 0, c, 0, d }, { 1, sc, 1, sd }, x)
            + G({ x / a }, 1) * G({ a, 0, c, 0, d }, { sa, 1, sc, 1, sd }, x)
            + G({ 0, a, 0, c, 0, d }, { 1, sa, 1, sc, 1, sd }, x)
            + (sy[2] + sy[3] + sy[4] + sy[5]) * Log(a, sa)
            + (-sy[2] - sy[3] - sy[4] - sy[5]) * Log(-x, sa) };
        if (c != x) {
            res += (sy[0] + sy[1] - G({ 0, 1, x / a }, 1) - G({ 0, x / a, 1 }, 1))
                * G({ c, 0, d }, { sc, 1, sd }, x);
        }
        if (d != x) {
            res += (sy[2] + sy[3] + sy[4] + sy[5] - G({ 0, 1, c / a, 0, x / a }, 1)
                       - G({ 0, c / a, 0, 1, x / a }, 1) - G({ 0, c / a, 0, x / a, 1 }, 1)
                       - G({ 0, c / a, 1, 0, x / a }, 1))
                * G({ d }, { sd }, x);
        }
        return res;

    } else { // abcd
        const vector<complex<double>> sy = { G({ b / a, 0, c / a }, 1), G({ b / a }, 1),
            G({ b / a, 0, c / a, 0, d / a }, 1) };
        complex<double> res { G({ 0, b / a, 0, c / a, 0, d / a }, 1)
            + 2. * G({ b / a, 0, 0, c / a, 0, d / a }, 1)
            + 2. * G({ b / a, 0, c / a, 0, 0, d / a }, 1)
            + G({ b / a, 0, c / a, 0, d / a, x / a }, 1)
            + G({ b / a, 0, c / a, x / a }, 1) * G({ 0, d }, { 1, sd }, x)
            + sy[0] * G({ 0, 0, d }, { 1, 1, sd }, x)
            + G({ b / a, x / a }, 1) * G({ 0, c, 0, d }, { 1, sc, 1, sd }, x)
            + sy[1] * G({ 0, 0, c, 0, d }, { 1, 1, sc, 1, sd }, x)
            + G({ 0, b, 0, c, 0, d }, { 1, sb, 1, sc, 1, sd }, x) - sy[2] * Log(a, sa)
            + sy[2] * Log(-x, sa) };
        if (b != x) {
            res += (-sy[1] + G({ x / a }, 1))
                * G({ b, 0, c, 0, d }, { sb, 1, sc, 1, sd }, x);
        }
        if (c != x) {
            res += (-sy[0] + G({ b / a, 0, x / a }, 1))
                * G({ c, 0, d }, { sc, 1, sd }, x);
        }
        if (d != x) {
            res += (-sy[2] + G({ b / a, 0, c / a, 0, x / a }, 1)) * G({ d }, { sd }, x);
        }
        return res;
    }
}

complex<double> G6_00abcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_00abcd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_00abcd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_00abcd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_00abcd_b(a, b, c, d, sa, sb, sc, sd, x);
}
complex<double> G6_a00bcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_a00bcd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_a00bcd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_a00bcd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_a00bcd_b(a, b, c, d, sa, sb, sc, sd, x);
}
complex<double> G6_ab00cd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_ab00cd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_ab00cd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_ab00cd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_ab00cd_b(a, b, c, d, sa, sb, sc, sd, x);
}
complex<double> G6_abc00d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_abc00d_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_abc00d_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_abc00d_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_abc00d_b(a, b, c, d, sa, sb, sc, sd, x);
}
complex<double> G6_0a0bcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_0a0bcd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_0a0bcd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_0a0bcd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_0a0bcd_b(a, b, c, d, sa, sb, sc, sd, x);
}
complex<double> G6_0ab0cd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_0ab0cd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_0ab0cd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_0ab0cd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_0ab0cd_b(a, b, c, d, sa, sb, sc, sd, x);
}
complex<double> G6_0abc0d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_0abc0d_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_0abc0d_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_0abc0d_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_0abc0d_b(a, b, c, d, sa, sb, sc, sd, x);
}
complex<double> G6_a0b0cd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_a0b0cd_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_a0b0cd_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_a0b0cd_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_a0b0cd_b(a, b, c, d, sa, sb, sc, sd, x);
}
complex<double> G6_a0bc0d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_a0bc0d_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_a0bc0d_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_a0bc0d_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_a0bc0d_b(a, b, c, d, sa, sb, sc, sd, x);
}
complex<double> G6_ab0c0d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x)
{
    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_ab0c0d_d(a, b, c, d, sa, sb, sc, sd, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_ab0c0d_c(a, b, c, d, sa, sb, sc, sd, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_ab0c0d_a(a, b, c, d, sa, sb, sc, sd, x);

    // b is smallest
    return G6_ab0c0d_b(a, b, c, d, sa, sb, sc, sd, x);
}
