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

complex<double> G6_0abcde_e(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { G({ d / e, c / e, b / e }, 1),
        G({ 0, a, b }, { 1, sa, sb }, x), G({ d / e, c / e }, 1),
        G({ 0, a, b, c }, { 1, sa, sb, sc }, x), G({ 0, a }, { 1, sa }, x),
        G({ d / e, c / e, b / e, a / e }, 1), Log(e, se),
        G({ 0, a, b, c, d }, { 1, sa, sb, sc, sd }, x), G({ d / e }, 1),
        G({ 0, d / e, c / e, b / e, a / e }, 1), G({ d / e, 0, c / e, b / e, a / e }, 1),
        G({ d / e, c / e, 0, b / e, a / e }, 1),
        G({ d / e, c / e, b / e, 0, a / e }, 1) };
    complex<double> res { sy[6]
            * (-sy[9] - sy[10] - sy[11] - sy[12] - sy[1] * sy[2] + sy[0] * sy[4]
                + sy[3] * sy[8])
        - sy[3] * G({ 0, d / e }, 1)
        + sy[1] * (sy[0] + G({ 0, d / e, c / e }, 1) + G({ d / e, 0, c / e }, 1))
        + sy[4]
            * (-sy[5] - G({ 0, d / e, c / e, b / e }, 1)
                - G({ d / e, 0, c / e, b / e }, 1) - G({ d / e, c / e, 0, b / e }, 1))
        + G({ 0, 0, d / e, c / e, b / e, a / e }, 1)
        + G({ 0, d / e, 0, c / e, b / e, a / e }, 1)
        + G({ 0, d / e, c / e, 0, b / e, a / e }, 1)
        + G({ 0, d / e, c / e, b / e, 0, a / e }, 1)
        + G({ d / e, 0, 0, c / e, b / e, a / e }, 1)
        + G({ d / e, 0, c / e, 0, b / e, a / e }, 1)
        + G({ d / e, 0, c / e, b / e, 0, a / e }, 1)
        + G({ d / e, c / e, 0, 0, b / e, a / e }, 1)
        + G({ d / e, c / e, 0, b / e, 0, a / e }, 1)
        + G({ d / e, c / e, b / e, 0, 0, a / e }, 1)
        - G({ d / e, c / e, b / e, a / e, 0, x / e }, 1)
        + sy[7] * (-G({ x / e }, 1) + G({ e }, { se }, x))
        + 2. * sy[0] * G({ 0, 0, a }, { 1, 1, sa }, x)
        + sy[2]
            * (-sy[3] - 2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
        + sy[8]
            * (sy[7] + 2. * G({ 0, 0, a, b, c }, { 1, 1, sa, sb, sc }, x)
                + G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x)
                + G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x))
        - 2. * G({ 0, 0, a, b, c, d }, { 1, 1, sa, sb, sc, sd }, x)
        - G({ 0, a, 0, b, c, d }, { 1, sa, 1, sb, sc, sd }, x)
        - G({ 0, a, b, 0, c, d }, { 1, sa, sb, 1, sc, sd }, x)
        - G({ 0, a, b, c, 0, d }, { 1, sa, sb, sc, 1, sd }, x)
        + (sy[9] + sy[10] + sy[11] + sy[12] + sy[1] * sy[2] - sy[0] * sy[4]
              - sy[5] * sy[6] - sy[3] * sy[8])
            * Log(-x, se)
        + (sy[5] * pow(sy[6], 2.)) / 2.
        + sy[5] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
    return res;
}
complex<double> G6_0abcde_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (d == e) {
        const vector<complex<double>> sy = { G({ 0, a, d }, { 1, sa, sd }, x),
            G({ c / d, b / d, 1 }, 1), G({ c / d, b / d, a / d }, 1),
            G({ 0, a, b, d }, { 1, sa, sb, sd }, x),
            G({ 0, a, b, c }, { 1, sa, sb, sc }, x), G({ c / d, b / d, a / d, 1 }, 1) };
        complex<double> res { -(sy[0] * sy[1]) + sy[0] * sy[2]
            + (sy[3] - sy[4]) * G({ c / d, 1 }, 1)
            + G({ c / d, b / d, a / d, 0, 0, 1 }, 1)
            - G({ c / d, b / d, a / d, 0, x / d, 1 }, 1)
            + (-G({ c / d, b / d, a / d, 0, 1 }, 1)
                  + G({ c / d, b / d, a / d, 0, x / d }, 1))
                * G({ d }, { sd }, x)
            + (-sy[5] + G({ c / d, b / d, 0, 1 }, 1)) * G({ 0, a }, { 1, sa }, x)
            + sy[5] * G({ 0, d }, { 1, sd }, x) - sy[2] * G({ 0, 0, d }, { 1, 1, sd }, x)
            + (sy[1] - G({ c / d, 0, 1 }, 1)) * G({ 0, a, b }, { 1, sa, sb }, x)
            + G({ c / d, b / d }, 1) * (-sy[3] + G({ 0, a, 0, d }, { 1, sa, 1, sd }, x))
            + G({ c / d }, 1)
                * (-G({ 0, a, b, 0, d }, { 1, sa, sb, 1, sd }, x)
                    + G({ 0, a, b, c, d }, { 1, sa, sb, sc, sd }, x))
            + G({ 0, a, b, c, 0, d }, { 1, sa, sb, sc, 1, sd }, x) - sy[4] * Zeta(2) };
        return res;
    } else {
        const vector<complex<double>> sy = { G({ c / d, b / d, a / d }, 1),
            G({ 0, a, e }, { 1, sa, se }, x), G({ 0, 0, a }, { 1, 1, sa }, x),
            G({ c / d, b / d, e / d }, 1), G({ c / d, e / d, b / d }, 1),
            G({ e / d, c / d, b / d }, 1), G({ 0, a, b }, { 1, sa, sb }, x),
            G({ 0, a, b, c }, { 1, sa, sb, sc }, x), G({ e / d, c / d }, 1),
            G({ 0, 0, a, b }, { 1, 1, sa, sb }, x),
            G({ 0, a, 0, b }, { 1, sa, 1, sb }, x),
            G({ 0, a, b, e }, { 1, sa, sb, se }, x), G({ c / d, e / d }, 1), Log(d, sd),
            G({ c / d, b / d, a / d, e / d }, 1), G({ c / d, b / d, e / d, a / d }, 1),
            G({ c / d, e / d, b / d, a / d }, 1), G({ e / d, c / d, b / d, a / d }, 1),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0, a }, { 1, sa }, x), G({ e / d }, 1),
            G({ 0, a, b, c, e }, { 1, sa, sb, sc, se }, x),
            G({ 0, c / d, b / d, a / d, e / d }, 1),
            G({ 0, c / d, b / d, e / d, a / d }, 1),
            G({ 0, c / d, e / d, b / d, a / d }, 1),
            G({ 0, e / d, c / d, b / d, a / d }, 1),
            G({ c / d, 0, b / d, a / d, e / d }, 1),
            G({ c / d, 0, b / d, e / d, a / d }, 1),
            G({ c / d, 0, e / d, b / d, a / d }, 1),
            G({ c / d, b / d, 0, a / d, e / d }, 1),
            G({ c / d, b / d, 0, e / d, a / d }, 1),
            G({ c / d, b / d, e / d, 0, a / d }, 1),
            G({ c / d, e / d, 0, b / d, a / d }, 1),
            G({ c / d, e / d, b / d, 0, a / d }, 1),
            G({ e / d, 0, c / d, b / d, a / d }, 1),
            G({ e / d, c / d, 0, b / d, a / d }, 1),
            G({ e / d, c / d, b / d, 0, a / d }, 1) };
        complex<double> res { (2. * sy[9] + sy[10] + sy[11]) * sy[12]
            + (sy[15] + sy[16] + sy[17]) * sy[18] + sy[0] * sy[1] - sy[1] * sy[3]
            - 2. * sy[2] * sy[3] + sy[2] * (-2. * sy[4] - 2. * sy[5])
            + (2. * sy[9] + sy[10] + sy[7]) * sy[8]
            + sy[13]
                * (sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27] + sy[28] + sy[29]
                    + sy[30] + sy[31] + sy[32] + sy[33] + sy[34] + sy[35] + sy[36]
                    + sy[19] * (-sy[3] - sy[4] - sy[5]) + sy[12] * sy[6] - sy[20] * sy[7]
                    + sy[6] * sy[8])
            + sy[7] * G({ 0, e / d }, 1)
            + sy[6]
                * (-sy[4] - sy[5] - G({ 0, c / d, e / d }, 1) - G({ 0, e / d, c / d }, 1)
                    - G({ c / d, 0, e / d }, 1) - G({ e / d, 0, c / d }, 1))
            + sy[19]
                * (sy[15] + sy[16] + sy[17] + G({ 0, c / d, b / d, e / d }, 1)
                    + G({ 0, c / d, e / d, b / d }, 1) + G({ 0, e / d, c / d, b / d }, 1)
                    + G({ c / d, 0, b / d, e / d }, 1) + G({ c / d, 0, e / d, b / d }, 1)
                    + G({ c / d, b / d, 0, e / d }, 1) + G({ c / d, e / d, 0, b / d }, 1)
                    + G({ e / d, 0, c / d, b / d }, 1) + G({ e / d, c / d, 0, b / d }, 1))
            - G({ 0, 0, c / d, b / d, a / d, e / d }, 1)
            - G({ 0, 0, c / d, b / d, e / d, a / d }, 1)
            - G({ 0, 0, c / d, e / d, b / d, a / d }, 1)
            - G({ 0, 0, e / d, c / d, b / d, a / d }, 1)
            - G({ 0, c / d, 0, b / d, a / d, e / d }, 1)
            - G({ 0, c / d, 0, b / d, e / d, a / d }, 1)
            - G({ 0, c / d, 0, e / d, b / d, a / d }, 1)
            - G({ 0, c / d, b / d, 0, a / d, e / d }, 1)
            - G({ 0, c / d, b / d, 0, e / d, a / d }, 1)
            - G({ 0, c / d, b / d, e / d, 0, a / d }, 1)
            - G({ 0, c / d, e / d, 0, b / d, a / d }, 1)
            - G({ 0, c / d, e / d, b / d, 0, a / d }, 1)
            - G({ 0, e / d, 0, c / d, b / d, a / d }, 1)
            - G({ 0, e / d, c / d, 0, b / d, a / d }, 1)
            - G({ 0, e / d, c / d, b / d, 0, a / d }, 1)
            - G({ c / d, 0, 0, b / d, a / d, e / d }, 1)
            - G({ c / d, 0, 0, b / d, e / d, a / d }, 1)
            - G({ c / d, 0, 0, e / d, b / d, a / d }, 1)
            - G({ c / d, 0, b / d, 0, a / d, e / d }, 1)
            - G({ c / d, 0, b / d, 0, e / d, a / d }, 1)
            - G({ c / d, 0, b / d, e / d, 0, a / d }, 1)
            - G({ c / d, 0, e / d, 0, b / d, a / d }, 1)
            - G({ c / d, 0, e / d, b / d, 0, a / d }, 1)
            - G({ c / d, b / d, 0, 0, a / d, e / d }, 1)
            - G({ c / d, b / d, 0, 0, e / d, a / d }, 1)
            - G({ c / d, b / d, 0, e / d, 0, a / d }, 1)
            + G({ c / d, b / d, a / d, 0, 0, e / d }, 1)
            + G({ c / d, b / d, a / d, 0, e / d, x / d }, 1)
            + G({ c / d, b / d, a / d, e / d, 0, x / d }, 1)
            - G({ c / d, b / d, e / d, 0, 0, a / d }, 1)
            + G({ c / d, b / d, e / d, a / d, 0, x / d }, 1)
            - G({ c / d, e / d, 0, 0, b / d, a / d }, 1)
            - G({ c / d, e / d, 0, b / d, 0, a / d }, 1)
            - G({ c / d, e / d, b / d, 0, 0, a / d }, 1)
            + G({ c / d, e / d, b / d, a / d, 0, x / d }, 1)
            - G({ e / d, 0, 0, c / d, b / d, a / d }, 1)
            - G({ e / d, 0, c / d, 0, b / d, a / d }, 1)
            - G({ e / d, 0, c / d, b / d, 0, a / d }, 1)
            - G({ e / d, c / d, 0, 0, b / d, a / d }, 1)
            - G({ e / d, c / d, 0, b / d, 0, a / d }, 1)
            - G({ e / d, c / d, b / d, 0, 0, a / d }, 1)
            + G({ e / d, c / d, b / d, a / d, 0, x / d }, 1)
            - sy[0] * G({ 0, 0, e }, { 1, 1, se }, x)
            + G({ c / d, b / d }, 1) * (-sy[11] + G({ 0, a, 0, e }, { 1, sa, 1, se }, x))
            + sy[20]
                * (-sy[21] - 2. * G({ 0, 0, a, b, c }, { 1, 1, sa, sb, sc }, x)
                    - G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x)
                    - G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x))
            + G({ c / d }, 1) * (sy[21] - G({ 0, a, b, 0, e }, { 1, sa, sb, 1, se }, x))
            + G({ 0, a, b, c, 0, e }, { 1, sa, sb, sc, 1, se }, x)
            + (sy[13] * (sy[14] + sy[15] + sy[16] + sy[17]) - sy[22] - sy[23] - sy[24]
                  - sy[25] - sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31] - sy[32]
                  - sy[33] - sy[34] - sy[35] - sy[36] + sy[19] * (sy[3] + sy[4] + sy[5])
                  - sy[12] * sy[6] + sy[20] * sy[7] - sy[6] * sy[8])
                * Log(-x, sd)
            + (-sy[14] / 2. - sy[15] / 2. - sy[16] / 2. - sy[17] / 2.) * pow(sy[13], 2.)
            + 2. * sy[15] * Zeta(2) + 2. * sy[16] * Zeta(2) + 2. * sy[17] * Zeta(2)
            + sy[14] * (sy[18] + G({ 0, e }, { 1, se }, x) + 2. * Zeta(2)) };
        if (e != x) {
            res += (-G({ c / d, b / d, a / d, 0, e / d }, 1)
                       + G({ c / d, b / d, a / d, 0, x / d }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_0abcde_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (c == d) {
        if (c == e) {
            const vector<complex<double>> sy = { G({ 0, a, c }, { 1, sa, sc }, x),
                G({ b / c, 1, 1 }, 1), G({ b / c, a / c, 1 }, 1),
                G({ 0, a, c, c }, { 1, sa, sc, sc }, x), G({ c, c }, { sc, sc }, x),
                G({ b / c, a / c, 0, 1 }, 1), G({ b / c, a / c, 1, 1 }, 1) };
            complex<double> res { -(sy[0] * sy[1]) + sy[0] * sy[2] + sy[4] * sy[5]
                - sy[4] * G({ b / c, a / c, 0, x / c }, 1)
                + G({ b / c, a / c, 0, 0, 1, 1 }, 1)
                - G({ b / c, a / c, 0, x / c, 1, 1 }, 1)
                + (-G({ b / c, a / c, 0, 1, 1 }, 1) + G({ b / c, a / c, 0, x / c, 1 }, 1))
                    * G({ c }, { sc }, x)
                + (-sy[6] + G({ b / c, 0, 1, 1 }, 1)) * G({ 0, a }, { 1, sa }, x)
                + (-sy[5] + sy[6]) * G({ 0, c }, { 1, sc }, x)
                - sy[2] * G({ 0, c, c }, { 1, sc, sc }, x)
                + G({ b / c, a / c }, 1)
                    * (-sy[3] + G({ 0, 0, c, c }, { 1, 1, sc, sc }, x))
                + G({ b / c, 1 }, 1) * (sy[3] - G({ 0, a, b, c }, { 1, sa, sb, sc }, x))
                + G({ b / c }, 1)
                    * (-G({ 0, a, 0, c, c }, { 1, sa, 1, sc, sc }, x)
                        + G({ 0, a, b, c, c }, { 1, sa, sb, sc, sc }, x))
                + G({ 0, a, b, 0, c, c }, { 1, sa, sb, 1, sc, sc }, x)
                + G({ 0, a, b }, { 1, sa, sb }, x) * (sy[1] - Zeta(3)) };
            return res;
        } else {
            const vector<complex<double>> sy = { G({ b / c, a / c, 1 }, 1),
                G({ b / c, e / c, 1 }, 1), G({ 0, 0, a }, { 1, 1, sa }, x),
                G({ e / c, 1, b / c }, 1), G({ e / c, b / c, 1 }, 1),
                G({ 0, a, b }, { 1, sa, sb }, x), G({ e / c, 1 }, 1),
                G({ 0, a, b, e }, { 1, sa, sb, se }, x),
                G({ 0, a, c, e }, { 1, sa, sc, se }, x), G({ c, e }, { sc, se }, x),
                G({ b / c, a / c, 0, 1 }, 1), G({ b / c, a / c, e / c, 1 }, 1),
                Log(c, sc), G({ b / c, e / c, 1, a / c }, 1),
                G({ b / c, e / c, a / c, 1 }, 1), G({ e / c, 1, b / c, a / c }, 1),
                G({ e / c, b / c, 1, a / c }, 1), G({ e / c, b / c, a / c, 1 }, 1),
                G({ 0, 0 }, { 1, 1 }, x), G({ 0, a }, { 1, sa }, x),
                G({ 0, b / c, a / c, e / c, 1 }, 1), G({ 0, b / c, e / c, 1, a / c }, 1),
                G({ 0, b / c, e / c, a / c, 1 }, 1), G({ 0, e / c, 1, b / c, a / c }, 1),
                G({ 0, e / c, b / c, 1, a / c }, 1), G({ 0, e / c, b / c, a / c, 1 }, 1),
                G({ b / c, 0, a / c, e / c, 1 }, 1), G({ b / c, 0, e / c, 1, a / c }, 1),
                G({ b / c, 0, e / c, a / c, 1 }, 1), G({ b / c, e / c, 0, 1, a / c }, 1),
                G({ b / c, e / c, 0, a / c, 1 }, 1), G({ b / c, e / c, 1, 0, a / c }, 1),
                G({ e / c, 0, 1, b / c, a / c }, 1), G({ e / c, 0, b / c, 1, a / c }, 1),
                G({ e / c, 0, b / c, a / c, 1 }, 1), G({ e / c, 1, 0, b / c, a / c }, 1),
                G({ e / c, 1, b / c, 0, a / c }, 1), G({ e / c, b / c, 0, 1, a / c }, 1),
                G({ e / c, b / c, 0, a / c, 1 }, 1),
                G({ e / c, b / c, 1, 0, a / c }, 1) };
            complex<double> res { sy[9] * sy[10]
                + (sy[13] + sy[14] + sy[15] + sy[16] + sy[17]) * sy[18]
                - 2. * sy[1] * sy[2] + sy[2] * (-2. * sy[3] - 2. * sy[4])
                + sy[12]
                    * (sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25] + sy[26]
                        + sy[27] + sy[28] + sy[29] + sy[30] + sy[31] + sy[32] + sy[33]
                        + sy[34] + sy[35] + sy[36] + sy[37] + sy[38] + sy[39]
                        + sy[19] * (-sy[1] - sy[3] - sy[4]) + sy[5] * sy[6])
                + (-sy[7] + sy[8]) * G({ b / c, 1 }, 1)
                + sy[5] * (-sy[3] - sy[4] - G({ 0, e / c, 1 }, 1))
                - sy[9] * G({ b / c, a / c, 0, x / c }, 1)
                + sy[19]
                    * (sy[13] + sy[14] + sy[15] + sy[16] + sy[17]
                        + G({ 0, b / c, e / c, 1 }, 1) + G({ 0, e / c, 1, b / c }, 1)
                        + G({ 0, e / c, b / c, 1 }, 1) + G({ b / c, 0, e / c, 1 }, 1)
                        + G({ e / c, 0, 1, b / c }, 1) + G({ e / c, 0, b / c, 1 }, 1)
                        + G({ e / c, 1, 0, b / c }, 1))
                - G({ 0, 0, b / c, a / c, e / c, 1 }, 1)
                - G({ 0, 0, b / c, e / c, 1, a / c }, 1)
                - G({ 0, 0, b / c, e / c, a / c, 1 }, 1)
                - G({ 0, 0, e / c, 1, b / c, a / c }, 1)
                - G({ 0, 0, e / c, b / c, 1, a / c }, 1)
                - G({ 0, 0, e / c, b / c, a / c, 1 }, 1)
                - G({ 0, b / c, 0, a / c, e / c, 1 }, 1)
                - G({ 0, b / c, 0, e / c, 1, a / c }, 1)
                - G({ 0, b / c, 0, e / c, a / c, 1 }, 1)
                - G({ 0, b / c, e / c, 0, 1, a / c }, 1)
                - G({ 0, b / c, e / c, 0, a / c, 1 }, 1)
                - G({ 0, b / c, e / c, 1, 0, a / c }, 1)
                - G({ 0, e / c, 0, 1, b / c, a / c }, 1)
                - G({ 0, e / c, 0, b / c, 1, a / c }, 1)
                - G({ 0, e / c, 0, b / c, a / c, 1 }, 1)
                - G({ 0, e / c, 1, 0, b / c, a / c }, 1)
                - G({ 0, e / c, 1, b / c, 0, a / c }, 1)
                - G({ 0, e / c, b / c, 0, 1, a / c }, 1)
                - G({ 0, e / c, b / c, 0, a / c, 1 }, 1)
                - G({ 0, e / c, b / c, 1, 0, a / c }, 1)
                - G({ b / c, 0, 0, a / c, e / c, 1 }, 1)
                - G({ b / c, 0, 0, e / c, 1, a / c }, 1)
                - G({ b / c, 0, 0, e / c, a / c, 1 }, 1)
                - G({ b / c, 0, e / c, 0, 1, a / c }, 1)
                - G({ b / c, 0, e / c, 0, a / c, 1 }, 1)
                - G({ b / c, 0, e / c, 1, 0, a / c }, 1)
                + G({ b / c, a / c, 0, 0, e / c, 1 }, 1)
                + G({ b / c, a / c, 0, e / c, 1, x / c }, 1)
                + G({ b / c, a / c, 0, e / c, x / c, 1 }, 1)
                + G({ b / c, a / c, e / c, 0, 1, x / c }, 1)
                + G({ b / c, a / c, e / c, 0, x / c, 1 }, 1)
                + G({ b / c, a / c, e / c, 1, 0, x / c }, 1)
                - G({ b / c, e / c, 0, 0, 1, a / c }, 1)
                - G({ b / c, e / c, 0, 0, a / c, 1 }, 1)
                - G({ b / c, e / c, 0, 1, 0, a / c }, 1)
                - G({ b / c, e / c, 1, 0, 0, a / c }, 1)
                + G({ b / c, e / c, 1, a / c, 0, x / c }, 1)
                + G({ b / c, e / c, a / c, 0, 1, x / c }, 1)
                + G({ b / c, e / c, a / c, 0, x / c, 1 }, 1)
                + G({ b / c, e / c, a / c, 1, 0, x / c }, 1)
                - G({ e / c, 0, 0, 1, b / c, a / c }, 1)
                - G({ e / c, 0, 0, b / c, 1, a / c }, 1)
                - G({ e / c, 0, 0, b / c, a / c, 1 }, 1)
                - G({ e / c, 0, 1, 0, b / c, a / c }, 1)
                - G({ e / c, 0, 1, b / c, 0, a / c }, 1)
                - G({ e / c, 0, b / c, 0, 1, a / c }, 1)
                - G({ e / c, 0, b / c, 0, a / c, 1 }, 1)
                - G({ e / c, 0, b / c, 1, 0, a / c }, 1)
                - G({ e / c, 1, 0, 0, b / c, a / c }, 1)
                - G({ e / c, 1, 0, b / c, 0, a / c }, 1)
                - G({ e / c, 1, b / c, 0, 0, a / c }, 1)
                + G({ e / c, 1, b / c, a / c, 0, x / c }, 1)
                - G({ e / c, b / c, 0, 0, 1, a / c }, 1)
                - G({ e / c, b / c, 0, 0, a / c, 1 }, 1)
                - G({ e / c, b / c, 0, 1, 0, a / c }, 1)
                - G({ e / c, b / c, 1, 0, 0, a / c }, 1)
                + G({ e / c, b / c, 1, a / c, 0, x / c }, 1)
                + G({ e / c, b / c, a / c, 0, 1, x / c }, 1)
                + G({ e / c, b / c, a / c, 0, x / c, 1 }, 1)
                + G({ e / c, b / c, a / c, 1, 0, x / c }, 1)
                + (-sy[10] + sy[11]) * G({ 0, e }, { 1, se }, x)
                + (sy[0] - sy[1]) * G({ 0, a, e }, { 1, sa, se }, x)
                - sy[0] * G({ 0, c, e }, { 1, sc, se }, x)
                + G({ b / c, a / c }, 1)
                    * (-sy[8] + G({ 0, 0, c, e }, { 1, 1, sc, se }, x))
                + sy[6]
                    * (sy[7] + 2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                        + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
                + G({ b / c }, 1)
                    * (-G({ 0, a, 0, c, e }, { 1, sa, 1, sc, se }, x)
                        + G({ 0, a, b, c, e }, { 1, sa, sb, sc, se }, x))
                + G({ 0, a, b, 0, c, e }, { 1, sa, sb, 1, sc, se }, x)
                + (sy[12] * (sy[11] + sy[13] + sy[14] + sy[15] + sy[16] + sy[17]) - sy[20]
                      - sy[21] - sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27]
                      - sy[28] - sy[29] - sy[30] - sy[31] - sy[32] - sy[33] - sy[34]
                      - sy[35] - sy[36] - sy[37] - sy[38] - sy[39]
                      + sy[19] * (sy[1] + sy[3] + sy[4]) - sy[5] * sy[6])
                    * Log(-x, sc)
                + (-sy[11] / 2. - sy[13] / 2. - sy[14] / 2. - sy[15] / 2. - sy[16] / 2.
                      - sy[17] / 2.)
                    * pow(sy[12], 2.)
                + 2. * sy[13] * Zeta(2) + 2. * sy[14] * Zeta(2) + 2. * sy[15] * Zeta(2)
                + 2. * sy[16] * Zeta(2) + 2. * sy[17] * Zeta(2)
                + sy[11] * (sy[18] + 2. * Zeta(2)) };
            if (e != x) {
                res += (-G({ b / c, a / c, 0, e / c, 1 }, 1)
                           + G({ b / c, a / c, 0, x / c, 1 }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else {
        const vector<complex<double>> sy = { G({ b / c, a / c, d / c }, 1),
            G({ 0, a, e }, { 1, sa, se }, x), G({ b / c, d / c, a / c }, 1),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ b / c, d / c, e / c }, 1),
            G({ d / c, b / c, a / c }, 1), G({ d / c, b / c, e / c }, 1),
            G({ 0, a, b }, { 1, sa, sb }, x), G({ d / c, e / c, b / c }, 1),
            G({ d / c, b / c }, 1), G({ 0, a, 0, e }, { 1, sa, 1, se }, x),
            G({ d / c, e / c }, 1), G({ 0, a, b, e }, { 1, sa, sb, se }, x),
            G({ 0, a, d, e }, { 1, sa, sd, se }, x), G({ b / c, a / c, 0, d / c }, 1),
            G({ b / c, a / c, d / c, e / c }, 1), G({ b / c, d / c, a / c, e / c }, 1),
            G({ d / c, b / c, a / c, e / c }, 1), G({ 0, 0 }, { 1, 1 }, x),
            G({ b / c, d / c, e / c, a / c }, 1), G({ d / c, b / c, e / c, a / c }, 1),
            G({ d / c, e / c, b / c, a / c }, 1), G({ 0, a }, { 1, sa }, x), Log(c, sc),
            G({ 0, a, b, d, e }, { 1, sa, sb, sd, se }, x),
            G({ 0, b / c, a / c, d / c, e / c }, 1),
            G({ 0, b / c, d / c, a / c, e / c }, 1),
            G({ 0, b / c, d / c, e / c, a / c }, 1),
            G({ 0, d / c, b / c, a / c, e / c }, 1),
            G({ 0, d / c, b / c, e / c, a / c }, 1),
            G({ 0, d / c, e / c, b / c, a / c }, 1),
            G({ b / c, 0, a / c, d / c, e / c }, 1),
            G({ b / c, 0, d / c, a / c, e / c }, 1),
            G({ b / c, 0, d / c, e / c, a / c }, 1),
            G({ b / c, d / c, 0, a / c, e / c }, 1),
            G({ b / c, d / c, 0, e / c, a / c }, 1),
            G({ b / c, d / c, e / c, 0, a / c }, 1),
            G({ d / c, 0, b / c, a / c, e / c }, 1),
            G({ d / c, 0, b / c, e / c, a / c }, 1),
            G({ d / c, 0, e / c, b / c, a / c }, 1),
            G({ d / c, b / c, 0, a / c, e / c }, 1),
            G({ d / c, b / c, 0, e / c, a / c }, 1),
            G({ d / c, b / c, e / c, 0, a / c }, 1),
            G({ d / c, e / c, 0, b / c, a / c }, 1),
            G({ d / c, e / c, b / c, 0, a / c }, 1) };
        complex<double> res { -(sy[9] * sy[10]) + sy[9] * sy[12]
            + sy[18] * (-sy[16] - sy[17] - sy[19] - sy[20] - sy[21]) - sy[1] * sy[2]
            + 2. * sy[3] * sy[4] + sy[1] * (sy[4] - sy[5] + sy[6])
            + sy[3] * (2. * sy[6] + 2. * sy[8])
            + sy[23]
                * (-sy[25] - sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31] - sy[32]
                    - sy[33] - sy[34] - sy[35] - sy[36] - sy[37] - sy[38] - sy[39]
                    - sy[40] - sy[41] - sy[42] - sy[43] - sy[44] - sy[11] * sy[7]
                    + sy[22] * (sy[4] + sy[6] + sy[8]))
            + (-sy[10] + sy[13]) * G({ b / c, d / c }, 1)
            + sy[7] * (sy[8] + G({ 0, d / c, e / c }, 1) + G({ d / c, 0, e / c }, 1))
            + sy[22]
                * (-sy[19] - sy[20] - sy[21] - G({ 0, b / c, d / c, e / c }, 1)
                    - G({ 0, d / c, b / c, e / c }, 1) - G({ 0, d / c, e / c, b / c }, 1)
                    - G({ b / c, 0, d / c, e / c }, 1) - G({ b / c, d / c, 0, e / c }, 1)
                    - G({ d / c, 0, b / c, e / c }, 1) - G({ d / c, 0, e / c, b / c }, 1)
                    - G({ d / c, b / c, 0, e / c }, 1) - G({ d / c, e / c, 0, b / c }, 1))
            + G({ 0, 0, b / c, a / c, d / c, e / c }, 1)
            + G({ 0, 0, b / c, d / c, a / c, e / c }, 1)
            + G({ 0, 0, b / c, d / c, e / c, a / c }, 1)
            + G({ 0, 0, d / c, b / c, a / c, e / c }, 1)
            + G({ 0, 0, d / c, b / c, e / c, a / c }, 1)
            + G({ 0, 0, d / c, e / c, b / c, a / c }, 1)
            + G({ 0, b / c, 0, a / c, d / c, e / c }, 1)
            + G({ 0, b / c, 0, d / c, a / c, e / c }, 1)
            + G({ 0, b / c, 0, d / c, e / c, a / c }, 1)
            + G({ 0, b / c, d / c, 0, a / c, e / c }, 1)
            + G({ 0, b / c, d / c, 0, e / c, a / c }, 1)
            + G({ 0, b / c, d / c, e / c, 0, a / c }, 1)
            + G({ 0, d / c, 0, b / c, a / c, e / c }, 1)
            + G({ 0, d / c, 0, b / c, e / c, a / c }, 1)
            + G({ 0, d / c, 0, e / c, b / c, a / c }, 1)
            + G({ 0, d / c, b / c, 0, a / c, e / c }, 1)
            + G({ 0, d / c, b / c, 0, e / c, a / c }, 1)
            + G({ 0, d / c, b / c, e / c, 0, a / c }, 1)
            + G({ 0, d / c, e / c, 0, b / c, a / c }, 1)
            + G({ 0, d / c, e / c, b / c, 0, a / c }, 1)
            + G({ b / c, 0, 0, a / c, d / c, e / c }, 1)
            + G({ b / c, 0, 0, d / c, a / c, e / c }, 1)
            + G({ b / c, 0, 0, d / c, e / c, a / c }, 1)
            + G({ b / c, 0, d / c, 0, a / c, e / c }, 1)
            + G({ b / c, 0, d / c, 0, e / c, a / c }, 1)
            + G({ b / c, 0, d / c, e / c, 0, a / c }, 1)
            - G({ b / c, a / c, 0, 0, d / c, e / c }, 1)
            - G({ b / c, a / c, 0, d / c, 0, e / c }, 1)
            - G({ b / c, a / c, 0, d / c, e / c, x / c }, 1)
            - G({ b / c, a / c, d / c, 0, 0, e / c }, 1)
            - G({ b / c, a / c, d / c, 0, e / c, x / c }, 1)
            - G({ b / c, a / c, d / c, e / c, 0, x / c }, 1)
            + G({ b / c, d / c, 0, 0, a / c, e / c }, 1)
            + G({ b / c, d / c, 0, 0, e / c, a / c }, 1)
            + G({ b / c, d / c, 0, e / c, 0, a / c }, 1)
            - G({ b / c, d / c, a / c, 0, 0, e / c }, 1)
            - G({ b / c, d / c, a / c, 0, e / c, x / c }, 1)
            - G({ b / c, d / c, a / c, e / c, 0, x / c }, 1)
            + G({ b / c, d / c, e / c, 0, 0, a / c }, 1)
            - G({ b / c, d / c, e / c, a / c, 0, x / c }, 1)
            + G({ d / c, 0, 0, b / c, a / c, e / c }, 1)
            + G({ d / c, 0, 0, b / c, e / c, a / c }, 1)
            + G({ d / c, 0, 0, e / c, b / c, a / c }, 1)
            + G({ d / c, 0, b / c, 0, a / c, e / c }, 1)
            + G({ d / c, 0, b / c, 0, e / c, a / c }, 1)
            + G({ d / c, 0, b / c, e / c, 0, a / c }, 1)
            + G({ d / c, 0, e / c, 0, b / c, a / c }, 1)
            + G({ d / c, 0, e / c, b / c, 0, a / c }, 1)
            + G({ d / c, b / c, 0, 0, a / c, e / c }, 1)
            + G({ d / c, b / c, 0, 0, e / c, a / c }, 1)
            + G({ d / c, b / c, 0, e / c, 0, a / c }, 1)
            - G({ d / c, b / c, a / c, 0, 0, e / c }, 1)
            - G({ d / c, b / c, a / c, 0, e / c, x / c }, 1)
            - G({ d / c, b / c, a / c, e / c, 0, x / c }, 1)
            + G({ d / c, b / c, e / c, 0, 0, a / c }, 1)
            - G({ d / c, b / c, e / c, a / c, 0, x / c }, 1)
            + G({ d / c, e / c, 0, 0, b / c, a / c }, 1)
            + G({ d / c, e / c, 0, b / c, 0, a / c }, 1)
            + G({ d / c, e / c, b / c, 0, 0, a / c }, 1)
            - G({ d / c, e / c, b / c, a / c, 0, x / c }, 1)
            + (-sy[14] - sy[15] - sy[16] - sy[17]) * G({ 0, e }, { 1, se }, x)
            + (sy[0] + sy[2] + sy[5]) * G({ 0, 0, e }, { 1, 1, se }, x)
            - sy[0] * G({ 0, d, e }, { 1, sd, se }, x)
            + G({ b / c, a / c }, 1) * (-sy[13] + G({ 0, 0, d, e }, { 1, 1, sd, se }, x))
            + sy[11]
                * (-sy[12] - 2. * G({ 0, 0, a, b }, { 1, 1, sa, sb }, x)
                    - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x))
            + G({ b / c }, 1) * (sy[24] - G({ 0, a, 0, d, e }, { 1, sa, 1, sd, se }, x))
            + G({ d / c }, 1) * (-sy[24] + G({ 0, a, b, 0, e }, { 1, sa, sb, 1, se }, x))
            + G({ 0, a, b, 0, d, e }, { 1, sa, sb, 1, sd, se }, x)
            + ((-sy[15] - sy[16] - sy[17] - sy[19] - sy[20] - sy[21]) * sy[23] + sy[25]
                  + sy[26] + sy[27] + sy[28] + sy[29] + sy[30] + sy[31] + sy[32] + sy[33]
                  + sy[34] + sy[35] + sy[36] + sy[37] + sy[38] + sy[39] + sy[40] + sy[41]
                  + sy[42] + sy[43] + sy[44] + sy[11] * sy[7]
                  + sy[22] * (-sy[4] - sy[6] - sy[8]))
                * Log(-x, sc)
            + (sy[15] / 2. + sy[16] / 2. + sy[17] / 2. + sy[19] / 2. + sy[20] / 2.
                  + sy[21] / 2.)
                * pow(sy[23], 2.)
            + sy[15] * (-sy[18] - 2. * Zeta(2)) - 2. * sy[16] * Zeta(2)
            - 2. * sy[17] * Zeta(2) - 2. * sy[19] * Zeta(2) - 2. * sy[20] * Zeta(2)
            - 2. * sy[21] * Zeta(2) };
        if (d != x) {
            res += (sy[14] - G({ b / c, a / c, 0, x / c }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (G({ b / c, a / c, 0, d / c, e / c }, 1)
                       - G({ b / c, a / c, 0, d / c, x / c }, 1)
                       + G({ b / c, a / c, d / c, 0, e / c }, 1)
                       - G({ b / c, a / c, d / c, 0, x / c }, 1)
                       + G({ b / c, d / c, a / c, 0, e / c }, 1)
                       - G({ b / c, d / c, a / c, 0, x / c }, 1)
                       + G({ d / c, b / c, a / c, 0, e / c }, 1)
                       - G({ d / c, b / c, a / c, 0, x / c }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_0abcde_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (b == c) {
        if (b == d) {
            if (b == e) { // abbbb
                const vector<complex<double>> sy
                    = { G({ a / b, 0, 1 }, 1), G({ a / b, 1, 1 }, 1),
                          G({ b, b, b }, { sb, sb, sb }, x), G({ b, b }, { sb, sb }, x),
                          G({ a / b, 0, 1, 1 }, 1), G({ a / b, 1, 1, 1 }, 1) };
                complex<double> res { -(sy[0] * sy[2]) + sy[3] * sy[4]
                    + sy[2] * G({ a / b, 0, x / b }, 1)
                    - sy[3] * G({ a / b, 0, x / b, 1 }, 1)
                    + G({ a / b, 0, 0, 1, 1, 1 }, 1) - G({ a / b, 0, x / b, 1, 1, 1 }, 1)
                    + (-G({ a / b, 0, 1, 1, 1 }, 1) + G({ a / b, 0, x / b, 1, 1 }, 1))
                        * G({ b }, { sb }, x)
                    + (-sy[4] + sy[5]) * G({ 0, b }, { 1, sb }, x)
                    + sy[1] * G({ 0, a, b }, { 1, sa, sb }, x)
                    + (sy[0] - sy[1]) * G({ 0, b, b }, { 1, sb, sb }, x)
                    + G({ a / b, 1 }, 1)
                        * (-G({ 0, a, b, b }, { 1, sa, sb, sb }, x)
                            + G({ 0, b, b, b }, { 1, sb, sb, sb }, x))
                    + G({ a / b }, 1)
                        * (-G({ 0, 0, b, b, b }, { 1, 1, sb, sb, sb }, x)
                            + G({ 0, a, b, b, b }, { 1, sa, sb, sb, sb }, x))
                    + G({ 0, a, 0, b, b, b }, { 1, sa, 1, sb, sb, sb }, x)
                    + G({ 0, a }, { 1, sa }, x) * (-sy[5] - Zeta(4)) };
                return res;
            } else { // abbbe
                const vector<complex<double>> sy = { G({ a / b, 0, 1 }, 1),
                    G({ a / b, 1, 1 }, 1), G({ b, b, e }, { sb, sb, se }, x),
                    G({ e / b, 1, 1 }, 1), G({ b, e }, { sb, se }, x),
                    G({ a / b, 0, 1, 1 }, 1), G({ a / b, e / b, 1, 1 }, 1), Log(b, sb),
                    G({ e / b, 1, 1, a / b }, 1), G({ e / b, 1, a / b, 1 }, 1),
                    G({ e / b, a / b, 1, 1 }, 1), G({ 0, 0 }, { 1, 1 }, x),
                    G({ 0, a }, { 1, sa }, x), G({ 0, a / b, e / b, 1, 1 }, 1),
                    G({ 0, e / b, 1, 1, a / b }, 1), G({ 0, e / b, 1, a / b, 1 }, 1),
                    G({ 0, e / b, a / b, 1, 1 }, 1), G({ e / b, 0, 1, 1, a / b }, 1),
                    G({ e / b, 0, 1, a / b, 1 }, 1), G({ e / b, 0, a / b, 1, 1 }, 1),
                    G({ e / b, 1, 0, 1, a / b }, 1), G({ e / b, 1, 0, a / b, 1 }, 1),
                    G({ e / b, 1, 1, 0, a / b }, 1) };
                complex<double> res { -(sy[0] * sy[2]) + sy[4] * sy[5]
                    + (sy[13] + sy[14] + sy[15] + sy[16] + sy[17] + sy[18] + sy[19]
                          + sy[20] + sy[21] + sy[22] - sy[12] * sy[3])
                        * sy[7]
                    + sy[11] * (sy[9] + sy[10] + sy[8])
                    + sy[2] * G({ a / b, 0, x / b }, 1)
                    + sy[12] * (sy[9] + sy[10] + sy[8] + G({ 0, e / b, 1, 1 }, 1))
                    - sy[4] * G({ a / b, 0, x / b, 1 }, 1)
                    - G({ 0, 0, a / b, e / b, 1, 1 }, 1)
                    - G({ 0, 0, e / b, 1, 1, a / b }, 1)
                    - G({ 0, 0, e / b, 1, a / b, 1 }, 1)
                    - G({ 0, 0, e / b, a / b, 1, 1 }, 1)
                    - G({ 0, e / b, 0, 1, 1, a / b }, 1)
                    - G({ 0, e / b, 0, 1, a / b, 1 }, 1)
                    - G({ 0, e / b, 0, a / b, 1, 1 }, 1)
                    - G({ 0, e / b, 1, 0, 1, a / b }, 1)
                    - G({ 0, e / b, 1, 0, a / b, 1 }, 1)
                    - G({ 0, e / b, 1, 1, 0, a / b }, 1)
                    + G({ a / b, 0, 0, e / b, 1, 1 }, 1)
                    + G({ a / b, 0, e / b, 1, 1, x / b }, 1)
                    + G({ a / b, 0, e / b, 1, x / b, 1 }, 1)
                    + G({ a / b, 0, e / b, x / b, 1, 1 }, 1)
                    + G({ a / b, e / b, 0, 1, 1, x / b }, 1)
                    + G({ a / b, e / b, 0, 1, x / b, 1 }, 1)
                    + G({ a / b, e / b, 0, x / b, 1, 1 }, 1)
                    + G({ a / b, e / b, 1, 0, 1, x / b }, 1)
                    + G({ a / b, e / b, 1, 0, x / b, 1 }, 1)
                    + G({ a / b, e / b, 1, 1, 0, x / b }, 1)
                    - G({ e / b, 0, 0, 1, 1, a / b }, 1)
                    - G({ e / b, 0, 0, 1, a / b, 1 }, 1)
                    - G({ e / b, 0, 0, a / b, 1, 1 }, 1)
                    - G({ e / b, 0, 1, 0, 1, a / b }, 1)
                    - G({ e / b, 0, 1, 0, a / b, 1 }, 1)
                    - G({ e / b, 0, 1, 1, 0, a / b }, 1)
                    - G({ e / b, 1, 0, 0, 1, a / b }, 1)
                    - G({ e / b, 1, 0, 0, a / b, 1 }, 1)
                    - G({ e / b, 1, 0, 1, 0, a / b }, 1)
                    - G({ e / b, 1, 1, 0, 0, a / b }, 1)
                    + G({ e / b, 1, 1, a / b, 0, x / b }, 1)
                    + G({ e / b, 1, a / b, 0, 1, x / b }, 1)
                    + G({ e / b, 1, a / b, 0, x / b, 1 }, 1)
                    + G({ e / b, 1, a / b, 1, 0, x / b }, 1)
                    + G({ e / b, a / b, 0, 1, 1, x / b }, 1)
                    + G({ e / b, a / b, 0, 1, x / b, 1 }, 1)
                    + G({ e / b, a / b, 0, x / b, 1, 1 }, 1)
                    + G({ e / b, a / b, 1, 0, 1, x / b }, 1)
                    + G({ e / b, a / b, 1, 0, x / b, 1 }, 1)
                    + G({ e / b, a / b, 1, 1, 0, x / b }, 1)
                    + (-sy[5] + sy[6]) * G({ 0, e }, { 1, se }, x)
                    - 2. * sy[3] * G({ 0, 0, a }, { 1, 1, sa }, x)
                    + (sy[1] - sy[3]) * G({ 0, a, e }, { 1, sa, se }, x)
                    + (sy[0] - sy[1]) * G({ 0, b, e }, { 1, sb, se }, x)
                    + G({ a / b, 1 }, 1)
                        * (-G({ 0, a, b, e }, { 1, sa, sb, se }, x)
                            + G({ 0, b, b, e }, { 1, sb, sb, se }, x))
                    + G({ a / b }, 1)
                        * (-G({ 0, 0, b, b, e }, { 1, 1, sb, sb, se }, x)
                            + G({ 0, a, b, b, e }, { 1, sa, sb, sb, se }, x))
                    + G({ 0, a, 0, b, b, e }, { 1, sa, 1, sb, sb, se }, x)
                    + (-sy[13] - sy[14] - sy[15] - sy[16] - sy[17] - sy[18] - sy[19]
                          - sy[20] - sy[21] - sy[22] + sy[12] * sy[3]
                          + sy[7] * (sy[9] + sy[10] + sy[6] + sy[8]))
                        * Log(-x, sb)
                    + (-sy[9] / 2. - sy[10] / 2. - sy[6] / 2. - sy[8] / 2.)
                        * pow(sy[7], 2.)
                    + 2. * sy[9] * Zeta(2) + 2. * sy[10] * Zeta(2) + 2. * sy[8] * Zeta(2)
                    + sy[6] * (sy[11] + 2. * Zeta(2)) };

                if (e != x) {
                    res += (-G({ a / b, 0, e / b, 1, 1 }, 1)
                               + G({ a / b, 0, x / b, 1, 1 }, 1))
                        * G({ e }, { se }, x);
                }
                return res;
            }
        } else { // abbde
            const vector<complex<double>> sy = { G({ a / b, 0, 1 }, 1),
                G({ a / b, d / b, 1 }, 1), G({ b, d, e }, { sb, sd, se }, x),
                G({ 0, a, e }, { 1, sa, se }, x), G({ d / b, 1, a / b }, 1),
                G({ 0, 0, a }, { 1, 1, sa }, x), G({ d / b, 1, e / b }, 1),
                G({ d / b, a / b, 1 }, 1), G({ d / b, e / b, 1 }, 1),
                G({ 0, a, d, e }, { 1, sa, sd, se }, x), G({ a / b, 0, d / b, 1 }, 1),
                G({ a / b, d / b, 1, e / b }, 1), G({ a / b, d / b, e / b, 1 }, 1),
                G({ d / b, 1, a / b, e / b }, 1), G({ d / b, a / b, 1, e / b }, 1),
                G({ d / b, a / b, e / b, 1 }, 1), G({ 0, a }, { 1, sa }, x),
                G({ d / b, 1, e / b, a / b }, 1), G({ d / b, e / b, 1, a / b }, 1),
                G({ d / b, e / b, a / b, 1 }, 1), G({ 0, 0 }, { 1, 1 }, x), Log(b, sb),
                G({ 0, a / b, d / b, 1, e / b }, 1), G({ 0, a / b, d / b, e / b, 1 }, 1),
                G({ 0, d / b, 1, a / b, e / b }, 1), G({ 0, d / b, 1, e / b, a / b }, 1),
                G({ 0, d / b, a / b, 1, e / b }, 1), G({ 0, d / b, a / b, e / b, 1 }, 1),
                G({ 0, d / b, e / b, 1, a / b }, 1), G({ 0, d / b, e / b, a / b, 1 }, 1),
                G({ d / b, 0, 1, a / b, e / b }, 1), G({ d / b, 0, 1, e / b, a / b }, 1),
                G({ d / b, 0, a / b, 1, e / b }, 1), G({ d / b, 0, a / b, e / b, 1 }, 1),
                G({ d / b, 0, e / b, 1, a / b }, 1), G({ d / b, 0, e / b, a / b, 1 }, 1),
                G({ d / b, 1, 0, a / b, e / b }, 1), G({ d / b, 1, 0, e / b, a / b }, 1),
                G({ d / b, 1, e / b, 0, a / b }, 1), G({ d / b, e / b, 0, 1, a / b }, 1),
                G({ d / b, e / b, 0, a / b, 1 }, 1),
                G({ d / b, e / b, 1, 0, a / b }, 1) };
            complex<double> res {
                (-sy[12] - sy[13] - sy[14] - sy[15] - sy[17] - sy[18] - sy[19]) * sy[20]
                - sy[0] * sy[2] - sy[3] * sy[4] + 2. * sy[5] * sy[6] + 2. * sy[5] * sy[8]
                + sy[3] * (sy[6] - sy[7] + sy[8])
                + sy[21]
                    * (-sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27] - sy[28]
                        - sy[29] - sy[30] - sy[31] - sy[32] - sy[33] - sy[34] - sy[35]
                        - sy[36] - sy[37] - sy[38] - sy[39] - sy[40] - sy[41]
                        + sy[16] * (sy[6] + sy[8]))
                + sy[2] * G({ a / b, 0, x / b }, 1)
                + sy[16]
                    * (-sy[17] - sy[18] - sy[19] - G({ 0, d / b, 1, e / b }, 1)
                        - G({ 0, d / b, e / b, 1 }, 1) - G({ d / b, 0, 1, e / b }, 1)
                        - G({ d / b, 0, e / b, 1 }, 1) - G({ d / b, 1, 0, e / b }, 1))
                + G({ 0, 0, a / b, d / b, 1, e / b }, 1)
                + G({ 0, 0, a / b, d / b, e / b, 1 }, 1)
                + G({ 0, 0, d / b, 1, a / b, e / b }, 1)
                + G({ 0, 0, d / b, 1, e / b, a / b }, 1)
                + G({ 0, 0, d / b, a / b, 1, e / b }, 1)
                + G({ 0, 0, d / b, a / b, e / b, 1 }, 1)
                + G({ 0, 0, d / b, e / b, 1, a / b }, 1)
                + G({ 0, 0, d / b, e / b, a / b, 1 }, 1)
                + G({ 0, d / b, 0, 1, a / b, e / b }, 1)
                + G({ 0, d / b, 0, 1, e / b, a / b }, 1)
                + G({ 0, d / b, 0, a / b, 1, e / b }, 1)
                + G({ 0, d / b, 0, a / b, e / b, 1 }, 1)
                + G({ 0, d / b, 0, e / b, 1, a / b }, 1)
                + G({ 0, d / b, 0, e / b, a / b, 1 }, 1)
                + G({ 0, d / b, 1, 0, a / b, e / b }, 1)
                + G({ 0, d / b, 1, 0, e / b, a / b }, 1)
                + G({ 0, d / b, 1, e / b, 0, a / b }, 1)
                + G({ 0, d / b, e / b, 0, 1, a / b }, 1)
                + G({ 0, d / b, e / b, 0, a / b, 1 }, 1)
                + G({ 0, d / b, e / b, 1, 0, a / b }, 1)
                - G({ a / b, 0, 0, d / b, 1, e / b }, 1)
                - G({ a / b, 0, 0, d / b, e / b, 1 }, 1)
                - G({ a / b, 0, d / b, 0, 1, e / b }, 1)
                - G({ a / b, 0, d / b, 0, e / b, 1 }, 1)
                - G({ a / b, 0, d / b, 1, 0, e / b }, 1)
                - G({ a / b, 0, d / b, 1, e / b, x / b }, 1)
                - G({ a / b, 0, d / b, e / b, 1, x / b }, 1)
                - G({ a / b, 0, d / b, e / b, x / b, 1 }, 1)
                - G({ a / b, d / b, 0, 0, 1, e / b }, 1)
                - G({ a / b, d / b, 0, 0, e / b, 1 }, 1)
                - G({ a / b, d / b, 0, 1, 0, e / b }, 1)
                - G({ a / b, d / b, 0, 1, e / b, x / b }, 1)
                - G({ a / b, d / b, 0, e / b, 1, x / b }, 1)
                - G({ a / b, d / b, 0, e / b, x / b, 1 }, 1)
                - G({ a / b, d / b, 1, 0, 0, e / b }, 1)
                - G({ a / b, d / b, 1, 0, e / b, x / b }, 1)
                - G({ a / b, d / b, 1, e / b, 0, x / b }, 1)
                - G({ a / b, d / b, e / b, 0, 1, x / b }, 1)
                - G({ a / b, d / b, e / b, 0, x / b, 1 }, 1)
                - G({ a / b, d / b, e / b, 1, 0, x / b }, 1)
                + G({ d / b, 0, 0, 1, a / b, e / b }, 1)
                + G({ d / b, 0, 0, 1, e / b, a / b }, 1)
                + G({ d / b, 0, 0, a / b, 1, e / b }, 1)
                + G({ d / b, 0, 0, a / b, e / b, 1 }, 1)
                + G({ d / b, 0, 0, e / b, 1, a / b }, 1)
                + G({ d / b, 0, 0, e / b, a / b, 1 }, 1)
                + G({ d / b, 0, 1, 0, a / b, e / b }, 1)
                + G({ d / b, 0, 1, 0, e / b, a / b }, 1)
                + G({ d / b, 0, 1, e / b, 0, a / b }, 1)
                + G({ d / b, 0, e / b, 0, 1, a / b }, 1)
                + G({ d / b, 0, e / b, 0, a / b, 1 }, 1)
                + G({ d / b, 0, e / b, 1, 0, a / b }, 1)
                + G({ d / b, 1, 0, 0, a / b, e / b }, 1)
                + G({ d / b, 1, 0, 0, e / b, a / b }, 1)
                + G({ d / b, 1, 0, e / b, 0, a / b }, 1)
                - G({ d / b, 1, a / b, 0, 0, e / b }, 1)
                - G({ d / b, 1, a / b, 0, e / b, x / b }, 1)
                - G({ d / b, 1, a / b, e / b, 0, x / b }, 1)
                + G({ d / b, 1, e / b, 0, 0, a / b }, 1)
                - G({ d / b, 1, e / b, a / b, 0, x / b }, 1)
                - G({ d / b, a / b, 0, 0, 1, e / b }, 1)
                - G({ d / b, a / b, 0, 0, e / b, 1 }, 1)
                - G({ d / b, a / b, 0, 1, 0, e / b }, 1)
                - G({ d / b, a / b, 0, 1, e / b, x / b }, 1)
                - G({ d / b, a / b, 0, e / b, 1, x / b }, 1)
                - G({ d / b, a / b, 0, e / b, x / b, 1 }, 1)
                - G({ d / b, a / b, 1, 0, 0, e / b }, 1)
                - G({ d / b, a / b, 1, 0, e / b, x / b }, 1)
                - G({ d / b, a / b, 1, e / b, 0, x / b }, 1)
                - G({ d / b, a / b, e / b, 0, 1, x / b }, 1)
                - G({ d / b, a / b, e / b, 0, x / b, 1 }, 1)
                - G({ d / b, a / b, e / b, 1, 0, x / b }, 1)
                + G({ d / b, e / b, 0, 0, 1, a / b }, 1)
                + G({ d / b, e / b, 0, 0, a / b, 1 }, 1)
                + G({ d / b, e / b, 0, 1, 0, a / b }, 1)
                + G({ d / b, e / b, 1, 0, 0, a / b }, 1)
                - G({ d / b, e / b, 1, a / b, 0, x / b }, 1)
                - G({ d / b, e / b, a / b, 0, 1, x / b }, 1)
                - G({ d / b, e / b, a / b, 0, x / b, 1 }, 1)
                - G({ d / b, e / b, a / b, 1, 0, x / b }, 1)
                + (-sy[10] - sy[11] - sy[12] - sy[13] - sy[14] - sy[15])
                    * G({ 0, e }, { 1, se }, x)
                + (sy[1] + sy[4] + sy[7]) * G({ 0, 0, e }, { 1, 1, se }, x)
                + (sy[0] - sy[1]) * G({ 0, d, e }, { 1, sd, se }, x)
                + G({ d / b, 1 }, 1) * (sy[9] - G({ 0, a, 0, e }, { 1, sa, 1, se }, x))
                + G({ a / b, 1 }, 1) * (-sy[9] + G({ 0, b, d, e }, { 1, sb, sd, se }, x))
                + G({ a / b }, 1)
                    * (-G({ 0, 0, b, d, e }, { 1, 1, sb, sd, se }, x)
                        + G({ 0, a, b, d, e }, { 1, sa, sb, sd, se }, x))
                + G({ 0, a, 0, b, d, e }, { 1, sa, 1, sb, sd, se }, x)
                + ((-sy[11] - sy[12] - sy[13] - sy[14] - sy[15] - sy[17] - sy[18]
                       - sy[19])
                          * sy[21]
                      + sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27] + sy[28]
                      + sy[29] + sy[30] + sy[31] + sy[32] + sy[33] + sy[34] + sy[35]
                      + sy[36] + sy[37] + sy[38] + sy[39] + sy[40] + sy[41]
                      + sy[16] * (-sy[6] - sy[8]))
                    * Log(-x, sb)
                + (sy[11] / 2. + sy[12] / 2. + sy[13] / 2. + sy[14] / 2. + sy[15] / 2.
                      + sy[17] / 2. + sy[18] / 2. + sy[19] / 2.)
                    * pow(sy[21], 2.)
                + sy[11] * (-sy[20] - 2. * Zeta(2)) - 2. * sy[12] * Zeta(2)
                - 2. * sy[13] * Zeta(2) - 2. * sy[14] * Zeta(2) - 2. * sy[15] * Zeta(2)
                - 2. * sy[17] * Zeta(2) - 2. * sy[18] * Zeta(2) - 2. * sy[19] * Zeta(2)
            };
            if (d != x) {
                res += (sy[10] - G({ a / b, 0, x / b, 1 }, 1))
                    * G({ d, e }, { sd, se }, x);
            }
            if (e != x) {
                res += (G({ a / b, 0, d / b, 1, e / b }, 1)
                           - G({ a / b, 0, d / b, 1, x / b }, 1)
                           + G({ a / b, 0, d / b, e / b, 1 }, 1)
                           - G({ a / b, 0, d / b, x / b, 1 }, 1)
                           + G({ a / b, d / b, 0, 1, e / b }, 1)
                           - G({ a / b, d / b, 0, 1, x / b }, 1)
                           + G({ a / b, d / b, 0, e / b, 1 }, 1)
                           - G({ a / b, d / b, 0, x / b, 1 }, 1)
                           + G({ a / b, d / b, 1, 0, e / b }, 1)
                           - G({ a / b, d / b, 1, 0, x / b }, 1)
                           + G({ d / b, 1, a / b, 0, e / b }, 1)
                           - G({ d / b, 1, a / b, 0, x / b }, 1)
                           + G({ d / b, a / b, 0, 1, e / b }, 1)
                           - G({ d / b, a / b, 0, 1, x / b }, 1)
                           + G({ d / b, a / b, 0, e / b, 1 }, 1)
                           - G({ d / b, a / b, 0, x / b, 1 }, 1)
                           + G({ d / b, a / b, 1, 0, e / b }, 1)
                           - G({ d / b, a / b, 1, 0, x / b }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else { // abcde
        const vector<complex<double>> sy = { G({ a / b, 0, c / b }, 1),
            G({ a / b, c / b, d / b }, 1), G({ c / b, a / b, d / b }, 1),
            G({ c / b, d / b, a / b }, 1), G({ 0, a, e }, { 1, sa, se }, x),
            G({ c / b, d / b, e / b }, 1), G({ c / b, a / b }, 1),
            G({ 0, 0, d, e }, { 1, 1, sd, se }, x),
            G({ 0, a, d, e }, { 1, sa, sd, se }, x), G({ a / b, 0, c / b, d / b }, 1),
            G({ a / b, c / b, 0, d / b }, 1), G({ a / b, c / b, d / b, e / b }, 1),
            G({ c / b, a / b, 0, d / b }, 1), G({ c / b, a / b, d / b, e / b }, 1),
            G({ c / b, d / b, a / b, e / b }, 1), Log(b, sb),
            G({ c / b, d / b, e / b, a / b }, 1), G({ 0, a }, { 1, sa }, x),
            G({ 0, 0 }, { 1, 1 }, x), G({ 0, a, c, d, e }, { 1, sa, sc, sd, se }, x),
            G({ 0, a / b, c / b, d / b, e / b }, 1),
            G({ 0, c / b, a / b, d / b, e / b }, 1),
            G({ 0, c / b, d / b, a / b, e / b }, 1),
            G({ 0, c / b, d / b, e / b, a / b }, 1),
            G({ c / b, 0, a / b, d / b, e / b }, 1),
            G({ c / b, 0, d / b, a / b, e / b }, 1),
            G({ c / b, 0, d / b, e / b, a / b }, 1),
            G({ c / b, d / b, 0, a / b, e / b }, 1),
            G({ c / b, d / b, 0, e / b, a / b }, 1),
            G({ c / b, d / b, e / b, 0, a / b }, 1) };
        complex<double> res { (sy[13] + sy[14] + sy[16]) * sy[18] + sy[3] * sy[4]
            - sy[4] * sy[5]
            + sy[15]
                * (sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27]
                    + sy[28] + sy[29] - sy[17] * sy[5])
            - sy[6] * sy[7] + sy[6] * sy[8]
            + sy[17]
                * (sy[16] + G({ 0, c / b, d / b, e / b }, 1)
                    + G({ c / b, 0, d / b, e / b }, 1) + G({ c / b, d / b, 0, e / b }, 1))
            - G({ 0, 0, a / b, c / b, d / b, e / b }, 1)
            - G({ 0, 0, c / b, a / b, d / b, e / b }, 1)
            - G({ 0, 0, c / b, d / b, a / b, e / b }, 1)
            - G({ 0, 0, c / b, d / b, e / b, a / b }, 1)
            - G({ 0, c / b, 0, a / b, d / b, e / b }, 1)
            - G({ 0, c / b, 0, d / b, a / b, e / b }, 1)
            - G({ 0, c / b, 0, d / b, e / b, a / b }, 1)
            - G({ 0, c / b, d / b, 0, a / b, e / b }, 1)
            - G({ 0, c / b, d / b, 0, e / b, a / b }, 1)
            - G({ 0, c / b, d / b, e / b, 0, a / b }, 1)
            + G({ a / b, 0, 0, c / b, d / b, e / b }, 1)
            + G({ a / b, 0, c / b, 0, d / b, e / b }, 1)
            + G({ a / b, 0, c / b, d / b, 0, e / b }, 1)
            + G({ a / b, 0, c / b, d / b, e / b, x / b }, 1)
            + G({ a / b, c / b, 0, 0, d / b, e / b }, 1)
            + G({ a / b, c / b, 0, d / b, 0, e / b }, 1)
            + G({ a / b, c / b, 0, d / b, e / b, x / b }, 1)
            + G({ a / b, c / b, d / b, 0, 0, e / b }, 1)
            + G({ a / b, c / b, d / b, 0, e / b, x / b }, 1)
            + G({ a / b, c / b, d / b, e / b, 0, x / b }, 1)
            - G({ c / b, 0, 0, a / b, d / b, e / b }, 1)
            - G({ c / b, 0, 0, d / b, a / b, e / b }, 1)
            - G({ c / b, 0, 0, d / b, e / b, a / b }, 1)
            - G({ c / b, 0, d / b, 0, a / b, e / b }, 1)
            - G({ c / b, 0, d / b, 0, e / b, a / b }, 1)
            - G({ c / b, 0, d / b, e / b, 0, a / b }, 1)
            + G({ c / b, a / b, 0, 0, d / b, e / b }, 1)
            + G({ c / b, a / b, 0, d / b, 0, e / b }, 1)
            + G({ c / b, a / b, 0, d / b, e / b, x / b }, 1)
            + G({ c / b, a / b, d / b, 0, 0, e / b }, 1)
            + G({ c / b, a / b, d / b, 0, e / b, x / b }, 1)
            + G({ c / b, a / b, d / b, e / b, 0, x / b }, 1)
            - G({ c / b, d / b, 0, 0, a / b, e / b }, 1)
            - G({ c / b, d / b, 0, 0, e / b, a / b }, 1)
            - G({ c / b, d / b, 0, e / b, 0, a / b }, 1)
            + G({ c / b, d / b, a / b, 0, 0, e / b }, 1)
            + G({ c / b, d / b, a / b, 0, e / b, x / b }, 1)
            + G({ c / b, d / b, a / b, e / b, 0, x / b }, 1)
            - G({ c / b, d / b, e / b, 0, 0, a / b }, 1)
            + G({ c / b, d / b, e / b, a / b, 0, x / b }, 1)
            + (sy[9] + sy[10] + sy[11] + sy[12] + sy[13] + sy[14])
                * G({ 0, e }, { 1, se }, x)
            - 2. * sy[5] * G({ 0, 0, a }, { 1, 1, sa }, x)
            + (-sy[1] - sy[2] - sy[3]) * G({ 0, 0, e }, { 1, 1, se }, x)
            + (sy[0] + sy[1] + sy[2]) * G({ 0, d, e }, { 1, sd, se }, x)
            + G({ c / b, d / b }, 1) * (-sy[8] + G({ 0, a, 0, e }, { 1, sa, 1, se }, x))
            + G({ a / b, c / b }, 1) * (-sy[7] + G({ 0, c, d, e }, { 1, sc, sd, se }, x))
            + G({ a / b }, 1) * (sy[19] - G({ 0, 0, c, d, e }, { 1, 1, sc, sd, se }, x))
            + G({ c / b }, 1) * (-sy[19] + G({ 0, a, 0, d, e }, { 1, sa, 1, sd, se }, x))
            + G({ 0, a, 0, c, d, e }, { 1, sa, 1, sc, sd, se }, x)
            + (sy[15] * (sy[11] + sy[13] + sy[14] + sy[16]) - sy[20] - sy[21] - sy[22]
                  - sy[23] - sy[24] - sy[25] - sy[26] - sy[27] - sy[28] - sy[29]
                  + sy[17] * sy[5])
                * Log(-x, sb)
            + (-sy[11] / 2. - sy[13] / 2. - sy[14] / 2. - sy[16] / 2.) * pow(sy[15], 2.)
            + 2. * sy[13] * Zeta(2) + 2. * sy[14] * Zeta(2) + 2. * sy[16] * Zeta(2)
            + sy[11] * (sy[18] + 2. * Zeta(2)) };
        if (c != x) {
            res += (-sy[0] + G({ a / b, 0, x / b }, 1))
                * G({ c, d, e }, { sc, sd, se }, x);
        }
        if (d != x) {
            res += (-sy[9] - sy[10] - sy[12] + G({ a / b, 0, c / b, x / b }, 1)
                       + G({ a / b, c / b, 0, x / b }, 1)
                       + G({ c / b, a / b, 0, x / b }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (-G({ a / b, 0, c / b, d / b, e / b }, 1)
                       + G({ a / b, 0, c / b, d / b, x / b }, 1)
                       - G({ a / b, c / b, 0, d / b, e / b }, 1)
                       + G({ a / b, c / b, 0, d / b, x / b }, 1)
                       - G({ a / b, c / b, d / b, 0, e / b }, 1)
                       + G({ a / b, c / b, d / b, 0, x / b }, 1)
                       - G({ c / b, a / b, 0, d / b, e / b }, 1)
                       + G({ c / b, a / b, 0, d / b, x / b }, 1)
                       - G({ c / b, a / b, d / b, 0, e / b }, 1)
                       + G({ c / b, a / b, d / b, 0, x / b }, 1)
                       - G({ c / b, d / b, a / b, 0, e / b }, 1)
                       + G({ c / b, d / b, a / b, 0, x / b }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_0abcde_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (a == b) {
        if (a == c) {
            if (a == d) { // aaaae
                const vector<complex<double>> sy = { G({ a, a, e }, { sa, sa, se }, x),
                    G({ a, a, a, e }, { sa, sa, sa, se }, x), Log(a, sa),
                    G({ e / a, 1, 1, 1 }, 1), G({ 0, e }, { 1, se }, x) };
                complex<double> res { -(sy[1] * G({ 0, x / a }, 1))
                    + sy[0] * G({ 0, x / a, 1 }, 1) + G({ 0, 0, e / a, 1, 1, 1 }, 1)
                    + G({ 0, e / a, 1, 1, 1, x / a }, 1)
                    + G({ 0, e / a, 1, 1, x / a, 1 }, 1)
                    + G({ 0, e / a, 1, x / a, 1, 1 }, 1)
                    + G({ 0, e / a, x / a, 1, 1, 1 }, 1)
                    + G({ e / a, 0, 1, 1, 1, x / a }, 1)
                    + G({ e / a, 0, 1, 1, x / a, 1 }, 1)
                    + G({ e / a, 0, 1, x / a, 1, 1 }, 1)
                    + G({ e / a, 0, x / a, 1, 1, 1 }, 1)
                    + G({ e / a, 1, 0, 1, 1, x / a }, 1)
                    + G({ e / a, 1, 0, 1, x / a, 1 }, 1)
                    + G({ e / a, 1, 0, x / a, 1, 1 }, 1)
                    + G({ e / a, 1, 1, 0, 1, x / a }, 1)
                    + G({ e / a, 1, 1, 0, x / a, 1 }, 1)
                    + G({ e / a, 1, 1, 1, 0, x / a }, 1)
                    + G({ 0, 0, a, a, a, e }, { 1, 1, sa, sa, sa, se }, x)
                    + sy[2] * sy[3] * Log(-x, sa) - (sy[3] * pow(sy[2], 2.)) / 2.
                    - sy[1] * Zeta(2) + G({ 0, a, a, e }, { 1, sa, sa, se }, x) * Zeta(2)
                    + sy[3] * (sy[4] + G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2))
                    + G({ a, e }, { sa, se }, x) * (-G({ 0, x / a, 1, 1 }, 1) - Zeta(4))
                    + sy[4] * Zeta(4) - sy[0] * Zeta(3)
                    + G({ 0, a, e }, { 1, sa, se }, x) * Zeta(3) };
                if (e != x) {
                    res += (-G({ 0, e / a, 1, 1, 1 }, 1) + G({ 0, x / a, 1, 1, 1 }, 1))
                        * G({ e }, { se }, x);
                }
                return res;
            } else { // aaade
                const vector<complex<double>> sy = { G({ a, d, e }, { sa, sd, se }, x),
                    G({ d / a, 1, 1 }, 1), G({ 0, d, e }, { 1, sd, se }, x),
                    G({ a, a, d, e }, { sa, sa, sd, se }, x), G({ 0, 0 }, { 1, 1 }, x),
                    G({ d / a, 1, e / a, 1 }, 1), G({ d / a, e / a, 1, 1 }, 1),
                    G({ 0, d / a, 1, 1 }, 1), G({ d / a, 1, 1, e / a }, 1), Log(a, sa) };
                complex<double> res { -(sy[1] * sy[2]) + sy[4] * (-sy[5] - sy[6])
                    - sy[3] * G({ 0, x / a }, 1) + sy[0] * G({ 0, x / a, 1 }, 1)
                    - G({ 0, 0, d / a, 1, 1, e / a }, 1)
                    - G({ 0, 0, d / a, 1, e / a, 1 }, 1)
                    - G({ 0, 0, d / a, e / a, 1, 1 }, 1)
                    - G({ 0, d / a, 0, 1, 1, e / a }, 1)
                    - G({ 0, d / a, 0, 1, e / a, 1 }, 1)
                    - G({ 0, d / a, 0, e / a, 1, 1 }, 1)
                    - G({ 0, d / a, 1, 0, 1, e / a }, 1)
                    - G({ 0, d / a, 1, 0, e / a, 1 }, 1)
                    - G({ 0, d / a, 1, 1, 0, e / a }, 1)
                    - G({ 0, d / a, 1, 1, e / a, x / a }, 1)
                    - G({ 0, d / a, 1, e / a, 1, x / a }, 1)
                    - G({ 0, d / a, 1, e / a, x / a, 1 }, 1)
                    - G({ 0, d / a, e / a, 1, 1, x / a }, 1)
                    - G({ 0, d / a, e / a, 1, x / a, 1 }, 1)
                    - G({ 0, d / a, e / a, x / a, 1, 1 }, 1)
                    - G({ d / a, 0, 0, 1, 1, e / a }, 1)
                    - G({ d / a, 0, 0, 1, e / a, 1 }, 1)
                    - G({ d / a, 0, 0, e / a, 1, 1 }, 1)
                    - G({ d / a, 0, 1, 0, 1, e / a }, 1)
                    - G({ d / a, 0, 1, 0, e / a, 1 }, 1)
                    - G({ d / a, 0, 1, 1, 0, e / a }, 1)
                    - G({ d / a, 0, 1, 1, e / a, x / a }, 1)
                    - G({ d / a, 0, 1, e / a, 1, x / a }, 1)
                    - G({ d / a, 0, 1, e / a, x / a, 1 }, 1)
                    - G({ d / a, 0, e / a, 1, 1, x / a }, 1)
                    - G({ d / a, 0, e / a, 1, x / a, 1 }, 1)
                    - G({ d / a, 0, e / a, x / a, 1, 1 }, 1)
                    - G({ d / a, 1, 0, 0, 1, e / a }, 1)
                    - G({ d / a, 1, 0, 0, e / a, 1 }, 1)
                    - G({ d / a, 1, 0, 1, 0, e / a }, 1)
                    - G({ d / a, 1, 0, 1, e / a, x / a }, 1)
                    - G({ d / a, 1, 0, e / a, 1, x / a }, 1)
                    - G({ d / a, 1, 0, e / a, x / a, 1 }, 1)
                    - G({ d / a, 1, 1, 0, 0, e / a }, 1)
                    - G({ d / a, 1, 1, 0, e / a, x / a }, 1)
                    - G({ d / a, 1, 1, e / a, 0, x / a }, 1)
                    - G({ d / a, 1, e / a, 0, 1, x / a }, 1)
                    - G({ d / a, 1, e / a, 0, x / a, 1 }, 1)
                    - G({ d / a, 1, e / a, 1, 0, x / a }, 1)
                    - G({ d / a, e / a, 0, 1, 1, x / a }, 1)
                    - G({ d / a, e / a, 0, 1, x / a, 1 }, 1)
                    - G({ d / a, e / a, 0, x / a, 1, 1 }, 1)
                    - G({ d / a, e / a, 1, 0, 1, x / a }, 1)
                    - G({ d / a, e / a, 1, 0, x / a, 1 }, 1)
                    - G({ d / a, e / a, 1, 1, 0, x / a }, 1)
                    + (-sy[5] - sy[6] - sy[7] - sy[8]) * G({ 0, e }, { 1, se }, x)
                    + sy[1] * G({ 0, 0, e }, { 1, 1, se }, x)
                    + G({ 0, 0, a, a, d, e }, { 1, 1, sa, sa, sd, se }, x)
                    + sy[9] * (-sy[5] - sy[6] - sy[8]) * Log(-x, sa)
                    + (sy[5] / 2. + sy[6] / 2. + sy[8] / 2.) * pow(sy[9], 2.)
                    + sy[8] * (-sy[4] - 2. * Zeta(2)) - sy[3] * Zeta(2)
                    - 2. * sy[5] * Zeta(2) - 2. * sy[6] * Zeta(2)
                    + G({ 0, a, d, e }, { 1, sa, sd, se }, x) * Zeta(2) - sy[0] * Zeta(3)
                    + sy[2] * Zeta(3) };
                if (d != x) {
                    res += (sy[7] - G({ 0, x / a, 1, 1 }, 1))
                        * G({ d, e }, { sd, se }, x);
                }
                if (e != x) {
                    res += (G({ 0, d / a, 1, 1, e / a }, 1)
                               - G({ 0, d / a, 1, 1, x / a }, 1)
                               + G({ 0, d / a, 1, e / a, 1 }, 1)
                               - G({ 0, d / a, 1, x / a, 1 }, 1)
                               + G({ 0, d / a, e / a, 1, 1 }, 1)
                               - G({ 0, d / a, x / a, 1, 1 }, 1)
                               + G({ d / a, 0, 1, 1, e / a }, 1)
                               - G({ d / a, 0, 1, 1, x / a }, 1)
                               + G({ d / a, 0, 1, e / a, 1 }, 1)
                               - G({ d / a, 0, 1, x / a, 1 }, 1)
                               + G({ d / a, 0, e / a, 1, 1 }, 1)
                               - G({ d / a, 0, x / a, 1, 1 }, 1)
                               + G({ d / a, 1, 0, 1, e / a }, 1)
                               - G({ d / a, 1, 0, 1, x / a }, 1)
                               + G({ d / a, 1, 0, e / a, 1 }, 1)
                               - G({ d / a, 1, 0, x / a, 1 }, 1)
                               + G({ d / a, 1, 1, 0, e / a }, 1)
                               - G({ d / a, 1, 1, 0, x / a }, 1))
                        * G({ e }, { se }, x);
                }
                return res;
            }
        } else { // aacde
            const vector<complex<double>> sy = { G({ 0, c / a, 1 }, 1),
                G({ 0, d, e }, { 1, sd, se }, x), G({ c / a, 1, d / a }, 1),
                G({ c / a, d / a, 1 }, 1), G({ 0, c, d, e }, { 1, sc, sd, se }, x),
                G({ a, c, d, e }, { sa, sc, sd, se }, x), Log(a, sa),
                G({ c / a, 1, d / a, e / a }, 1), G({ c / a, d / a, 1, e / a }, 1),
                G({ c / a, d / a, e / a, 1 }, 1), G({ 0, 0 }, { 1, 1 }, x),
                G({ 0, c / a, 1, d / a }, 1), G({ 0, c / a, d / a, 1 }, 1),
                G({ c / a, 0, 1, d / a }, 1), G({ c / a, 0, d / a, 1 }, 1),
                G({ c / a, 1, 0, d / a }, 1) };
            complex<double> res { sy[0] * sy[1] + sy[1] * (sy[2] + sy[3])
                + sy[10] * (sy[9] + sy[8]) - sy[5] * G({ 0, x / a }, 1)
                + G({ 0, 0, c / a, 1, d / a, e / a }, 1)
                + G({ 0, 0, c / a, d / a, 1, e / a }, 1)
                + G({ 0, 0, c / a, d / a, e / a, 1 }, 1)
                + G({ 0, c / a, 0, 1, d / a, e / a }, 1)
                + G({ 0, c / a, 0, d / a, 1, e / a }, 1)
                + G({ 0, c / a, 0, d / a, e / a, 1 }, 1)
                + G({ 0, c / a, 1, 0, d / a, e / a }, 1)
                + G({ 0, c / a, 1, d / a, 0, e / a }, 1)
                + G({ 0, c / a, 1, d / a, e / a, x / a }, 1)
                + G({ 0, c / a, d / a, 0, 1, e / a }, 1)
                + G({ 0, c / a, d / a, 0, e / a, 1 }, 1)
                + G({ 0, c / a, d / a, 1, 0, e / a }, 1)
                + G({ 0, c / a, d / a, 1, e / a, x / a }, 1)
                + G({ 0, c / a, d / a, e / a, 1, x / a }, 1)
                + G({ 0, c / a, d / a, e / a, x / a, 1 }, 1)
                + G({ c / a, 0, 0, 1, d / a, e / a }, 1)
                + G({ c / a, 0, 0, d / a, 1, e / a }, 1)
                + G({ c / a, 0, 0, d / a, e / a, 1 }, 1)
                + G({ c / a, 0, 1, 0, d / a, e / a }, 1)
                + G({ c / a, 0, 1, d / a, 0, e / a }, 1)
                + G({ c / a, 0, 1, d / a, e / a, x / a }, 1)
                + G({ c / a, 0, d / a, 0, 1, e / a }, 1)
                + G({ c / a, 0, d / a, 0, e / a, 1 }, 1)
                + G({ c / a, 0, d / a, 1, 0, e / a }, 1)
                + G({ c / a, 0, d / a, 1, e / a, x / a }, 1)
                + G({ c / a, 0, d / a, e / a, 1, x / a }, 1)
                + G({ c / a, 0, d / a, e / a, x / a, 1 }, 1)
                + G({ c / a, 1, 0, 0, d / a, e / a }, 1)
                + G({ c / a, 1, 0, d / a, 0, e / a }, 1)
                + G({ c / a, 1, 0, d / a, e / a, x / a }, 1)
                + G({ c / a, 1, d / a, 0, 0, e / a }, 1)
                + G({ c / a, 1, d / a, 0, e / a, x / a }, 1)
                + G({ c / a, 1, d / a, e / a, 0, x / a }, 1)
                + G({ c / a, d / a, 0, 0, 1, e / a }, 1)
                + G({ c / a, d / a, 0, 0, e / a, 1 }, 1)
                + G({ c / a, d / a, 0, 1, 0, e / a }, 1)
                + G({ c / a, d / a, 0, 1, e / a, x / a }, 1)
                + G({ c / a, d / a, 0, e / a, 1, x / a }, 1)
                + G({ c / a, d / a, 0, e / a, x / a, 1 }, 1)
                + G({ c / a, d / a, 1, 0, 0, e / a }, 1)
                + G({ c / a, d / a, 1, 0, e / a, x / a }, 1)
                + G({ c / a, d / a, 1, e / a, 0, x / a }, 1)
                + G({ c / a, d / a, e / a, 0, 1, x / a }, 1)
                + G({ c / a, d / a, e / a, 0, x / a, 1 }, 1)
                + G({ c / a, d / a, e / a, 1, 0, x / a }, 1)
                + (sy[9] + sy[11] + sy[12] + sy[13] + sy[14] + sy[15] + sy[7] + sy[8])
                    * G({ 0, e }, { 1, se }, x)
                + (-sy[2] - sy[3]) * G({ 0, 0, e }, { 1, 1, se }, x)
                + G({ c / a, 1 }, 1) * (sy[4] - G({ 0, 0, d, e }, { 1, 1, sd, se }, x))
                + G({ 0, 0, a, c, d, e }, { 1, 1, sa, sc, sd, se }, x)
                + sy[6] * (sy[9] + sy[7] + sy[8]) * Log(-x, sa)
                + (-sy[9] / 2. - sy[7] / 2. - sy[8] / 2.) * pow(sy[6], 2.)
                + 2. * sy[9] * Zeta(2) + sy[4] * Zeta(2) - sy[5] * Zeta(2)
                + 2. * sy[8] * Zeta(2) + sy[7] * (sy[10] + 2. * Zeta(2)) };
            if (c != x) {
                res += (-sy[0] + G({ 0, x / a, 1 }, 1))
                    * G({ c, d, e }, { sc, sd, se }, x);
            }
            if (d != x) {
                res += (-sy[11] - sy[12] - sy[13] - sy[14] - sy[15]
                           + G({ 0, c / a, 1, x / a }, 1) + G({ 0, c / a, x / a, 1 }, 1)
                           + G({ c / a, 0, 1, x / a }, 1) + G({ c / a, 0, x / a, 1 }, 1)
                           + G({ c / a, 1, 0, x / a }, 1))
                    * G({ d, e }, { sd, se }, x);
            }
            if (e != x) {
                res += (-G({ 0, c / a, 1, d / a, e / a }, 1)
                           + G({ 0, c / a, 1, d / a, x / a }, 1)
                           - G({ 0, c / a, d / a, 1, e / a }, 1)
                           + G({ 0, c / a, d / a, 1, x / a }, 1)
                           - G({ 0, c / a, d / a, e / a, 1 }, 1)
                           + G({ 0, c / a, d / a, x / a, 1 }, 1)
                           - G({ c / a, 0, 1, d / a, e / a }, 1)
                           + G({ c / a, 0, 1, d / a, x / a }, 1)
                           - G({ c / a, 0, d / a, 1, e / a }, 1)
                           + G({ c / a, 0, d / a, 1, x / a }, 1)
                           - G({ c / a, 0, d / a, e / a, 1 }, 1)
                           + G({ c / a, 0, d / a, x / a, 1 }, 1)
                           - G({ c / a, 1, 0, d / a, e / a }, 1)
                           + G({ c / a, 1, 0, d / a, x / a }, 1)
                           - G({ c / a, 1, d / a, 0, e / a }, 1)
                           + G({ c / a, 1, d / a, 0, x / a }, 1)
                           - G({ c / a, d / a, 0, 1, e / a }, 1)
                           + G({ c / a, d / a, 0, 1, x / a }, 1)
                           - G({ c / a, d / a, 0, e / a, 1 }, 1)
                           + G({ c / a, d / a, 0, x / a, 1 }, 1)
                           - G({ c / a, d / a, 1, 0, e / a }, 1)
                           + G({ c / a, d / a, 1, 0, x / a }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else { // abcde
        const vector<complex<double>> sy = { G({ 0, b / a, c / a }, 1),
            G({ 0, d, e }, { 1, sd, se }, x), G({ b / a, 0, c / a }, 1),
            G({ b / a, c / a, d / a }, 1), G({ 0, c, d, e }, { 1, sc, sd, se }, x),
            G({ 0, b / a }, 1), G({ 0, b / a, c / a, d / a }, 1),
            G({ b / a, 0, c / a, d / a }, 1), G({ b / a, c / a, 0, d / a }, 1),
            G({ b / a, c / a, d / a, e / a }, 1), Log(a, sa) };
        complex<double> res { -(sy[0] * sy[1]) + sy[1] * (-sy[2] - sy[3]) - sy[4] * sy[5]
            - G({ 0, 0, b / a, c / a, d / a, e / a }, 1)
            - G({ 0, b / a, 0, c / a, d / a, e / a }, 1)
            - G({ 0, b / a, c / a, 0, d / a, e / a }, 1)
            - G({ 0, b / a, c / a, d / a, 0, e / a }, 1)
            - G({ 0, b / a, c / a, d / a, e / a, x / a }, 1)
            - G({ b / a, 0, 0, c / a, d / a, e / a }, 1)
            - G({ b / a, 0, c / a, 0, d / a, e / a }, 1)
            - G({ b / a, 0, c / a, d / a, 0, e / a }, 1)
            - G({ b / a, 0, c / a, d / a, e / a, x / a }, 1)
            - G({ b / a, c / a, 0, 0, d / a, e / a }, 1)
            - G({ b / a, c / a, 0, d / a, 0, e / a }, 1)
            - G({ b / a, c / a, 0, d / a, e / a, x / a }, 1)
            - G({ b / a, c / a, d / a, 0, 0, e / a }, 1)
            - G({ b / a, c / a, d / a, 0, e / a, x / a }, 1)
            - G({ b / a, c / a, d / a, e / a, 0, x / a }, 1)
            + (-sy[9] - sy[6] - sy[7] - sy[8]) * G({ 0, e }, { 1, se }, x)
            + sy[3] * G({ 0, 0, e }, { 1, 1, se }, x)
            + G({ b / a, c / a }, 1) * (-sy[4] + G({ 0, 0, d, e }, { 1, 1, sd, se }, x))
            + G({ b / a }, 1)
                * (G({ 0, 0, c, d, e }, { 1, 1, sc, sd, se }, x)
                    - G({ 0, b, c, d, e }, { 1, sb, sc, sd, se }, x))
            + G({ 0, 0, b, c, d, e }, { 1, 1, sb, sc, sd, se }, x)
            - sy[9] * sy[10] * Log(-x, sa) + (sy[9] * pow(sy[10], 2.)) / 2.
            + sy[9] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)) };
        if (b != x) {
            res += (sy[5] - G({ 0, x / a }, 1))
                * G({ b, c, d, e }, { sb, sc, sd, se }, x);
        }
        if (c != x) {
            res += (sy[0] + sy[2] - G({ 0, b / a, x / a }, 1) - G({ b / a, 0, x / a }, 1))
                * G({ c, d, e }, { sc, sd, se }, x);
        }
        if (d != x) {
            res += (sy[6] + sy[7] + sy[8] - G({ 0, b / a, c / a, x / a }, 1)
                       - G({ b / a, 0, c / a, x / a }, 1)
                       - G({ b / a, c / a, 0, x / a }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (G({ 0, b / a, c / a, d / a, e / a }, 1)
                       - G({ 0, b / a, c / a, d / a, x / a }, 1)
                       + G({ b / a, 0, c / a, d / a, e / a }, 1)
                       - G({ b / a, 0, c / a, d / a, x / a }, 1)
                       + G({ b / a, c / a, 0, d / a, e / a }, 1)
                       - G({ b / a, c / a, 0, d / a, x / a }, 1)
                       + G({ b / a, c / a, d / a, 0, e / a }, 1)
                       - G({ b / a, c / a, d / a, 0, x / a }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}

complex<double> G6_a0bcde_e(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { Log(e, se), G({ a }, { sa }, x),
        G({ d / e, c / e, b / e }, 1), G({ a, 0, b }, { sa, 1, sb }, x),
        G({ d / e, c / e }, 1), G({ a, 0, b, c }, { sa, 1, sb, sc }, x),
        G({ 0, a }, { 1, sa }, x), G({ 0, d / e, c / e, b / e }, 1),
        G({ d / e, 0, c / e, b / e }, 1), G({ d / e, c / e, 0, b / e }, 1),
        G({ a, 0, b, c, d }, { sa, 1, sb, sc, sd }, x), G({ d / e }, 1),
        G({ d / e, c / e, b / e, 0, a / e }, 1) };
    complex<double> res { sy[6] * (sy[9] + sy[7] + sy[8])
        + sy[0]
            * (sy[12] - sy[3] * sy[4] + sy[11] * sy[5] - sy[2] * sy[6]
                + sy[1] * (sy[9] + sy[7] + sy[8]))
        - sy[5] * G({ 0, d / e }, 1)
        + sy[3] * (sy[2] + G({ 0, d / e, c / e }, 1) + G({ d / e, 0, c / e }, 1))
        - G({ 0, d / e, c / e, b / e, 0, a / e }, 1)
        - G({ d / e, 0, c / e, b / e, 0, a / e }, 1)
        - G({ d / e, c / e, 0, b / e, 0, a / e }, 1)
        - 2. * G({ d / e, c / e, b / e, 0, 0, a / e }, 1)
        - G({ d / e, c / e, b / e, 0, a / e, x / e }, 1)
        + sy[10] * (-G({ x / e }, 1) + G({ e }, { se }, x))
        - sy[2] * G({ 0, 0, a }, { 1, 1, sa }, x)
        + sy[4]
            * (-sy[5] - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                - 2. * G({ a, 0, 0, b }, { sa, 1, 1, sb }, x))
        + sy[11]
            * (sy[10] + G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x)
                + 2. * G({ a, 0, 0, b, c }, { sa, 1, 1, sb, sc }, x)
                + G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x))
        - G({ 0, a, 0, b, c, d }, { 1, sa, 1, sb, sc, sd }, x)
        - 2. * G({ a, 0, 0, b, c, d }, { sa, 1, 1, sb, sc, sd }, x)
        - G({ a, 0, b, 0, c, d }, { sa, 1, sb, 1, sc, sd }, x)
        - G({ a, 0, b, c, 0, d }, { sa, 1, sb, sc, 1, sd }, x)
        + (-sy[12] + sy[0] * sy[1] * sy[2] + sy[3] * sy[4] - sy[11] * sy[5]
              + sy[2] * sy[6] + sy[1] * (-sy[9] - sy[7] - sy[8]))
            * Log(-x, se)
        - (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.
        + sy[1]
            * (sy[12] - G({ 0, 0, d / e, c / e, b / e }, 1)
                - G({ 0, d / e, 0, c / e, b / e }, 1)
                - G({ 0, d / e, c / e, 0, b / e }, 1)
                - G({ d / e, 0, 0, c / e, b / e }, 1)
                - G({ d / e, 0, c / e, 0, b / e }, 1)
                - G({ d / e, c / e, 0, 0, b / e }, 1)
                + sy[2] * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2))) };
    return res;
}
complex<double> G6_a0bcde_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (d == e) {
        const vector<complex<double>> sy = { G({ c / d, b / d, 1 }, 1),
            G({ a, 0, b, d }, { sa, 1, sb, sd }, x),
            G({ a, 0, b, c }, { sa, 1, sb, sc }, x), G({ c / d, b / d, 0, a / d }, 1),
            G({ d }, { sd }, x), G({ c / d, b / d, 0, a / d, 1 }, 1) };
        complex<double> res { -(sy[4] * sy[5]) + (sy[1] - sy[2]) * G({ c / d, 1 }, 1)
            + sy[4] * G({ c / d, b / d, 0, a / d, x / d }, 1)
            + G({ c / d, b / d, 0, a / d, 0, 1 }, 1)
            - G({ c / d, b / d, 0, a / d, x / d, 1 }, 1)
            + (sy[5] - G({ c / d, b / d, 0, 0, 1 }, 1)) * G({ a }, { sa }, x)
            + sy[3] * G({ 0, d }, { 1, sd }, x)
            + (-sy[3] + G({ c / d, b / d, 0, 1 }, 1)) * G({ a, d }, { sa, sd }, x)
            + (sy[0] - G({ c / d, 0, 1 }, 1)) * G({ a, 0, b }, { sa, 1, sb }, x)
            - sy[0] * G({ a, 0, d }, { sa, 1, sd }, x)
            + G({ c / d, b / d }, 1) * (-sy[1] + G({ a, 0, 0, d }, { sa, 1, 1, sd }, x))
            + G({ c / d }, 1)
                * (-G({ a, 0, b, 0, d }, { sa, 1, sb, 1, sd }, x)
                    + G({ a, 0, b, c, d }, { sa, 1, sb, sc, sd }, x))
            + G({ a, 0, b, c, 0, d }, { sa, 1, sb, sc, 1, sd }, x) - sy[2] * Zeta(2) };
        return res;
    } else {
        const vector<complex<double>> sy = { G({ 0, 0, a }, { 1, 1, sa }, x),
            G({ c / d, b / d, e / d }, 1), G({ a, 0, b }, { sa, 1, sb }, x),
            G({ c / d, e / d, b / d }, 1), G({ e / d, c / d, b / d }, 1), Log(d, sd),
            G({ a }, { sa }, x), G({ a, 0, b, c }, { sa, 1, sb, sc }, x),
            G({ e / d, c / d }, 1), G({ 0, a, 0, b }, { 1, sa, 1, sb }, x),
            G({ a, 0, 0, b }, { sa, 1, 1, sb }, x),
            G({ a, 0, b, e }, { sa, 1, sb, se }, x), G({ c / d, e / d }, 1),
            G({ c / d, b / d, 0, a / d }, 1), G({ a, e }, { sa, se }, x),
            G({ 0, a }, { 1, sa }, x), G({ 0, c / d, b / d, e / d }, 1),
            G({ 0, c / d, e / d, b / d }, 1), G({ 0, e / d, c / d, b / d }, 1),
            G({ c / d, 0, b / d, e / d }, 1), G({ c / d, 0, e / d, b / d }, 1),
            G({ c / d, e / d, 0, b / d }, 1), G({ e / d, 0, c / d, b / d }, 1),
            G({ e / d, c / d, 0, b / d }, 1), G({ e / d }, 1),
            G({ a, 0, b, c, e }, { sa, 1, sb, sc, se }, x),
            G({ c / d, b / d, 0, a / d, e / d }, 1),
            G({ c / d, b / d, 0, e / d, a / d }, 1),
            G({ c / d, b / d, e / d, 0, a / d }, 1),
            G({ c / d, e / d, b / d, 0, a / d }, 1),
            G({ e / d, c / d, b / d, 0, a / d }, 1), G({ 0, 0 }, { 1, 1 }, x) };
        complex<double> res { (sy[9] + 2. * sy[10] + sy[11]) * sy[12] - sy[13] * sy[14]
            + sy[15]
                * (-sy[16] - sy[17] - sy[18] - sy[19] - sy[20] - sy[21] - sy[22] - sy[23])
            + sy[0] * (sy[3] + sy[4]) + (sy[9] + 2. * sy[10] + sy[7]) * sy[8]
            + sy[5]
                * (sy[15] * sy[1] - sy[26] - sy[27] - sy[28] + sy[12] * sy[2] - sy[29]
                    - sy[30] + sy[15] * (sy[3] + sy[4])
                    + (-sy[16] - sy[17] - sy[18] - sy[19] - sy[20] - sy[21] - sy[22]
                          - sy[23])
                        * sy[6]
                    - sy[24] * sy[7] + sy[2] * sy[8])
            + sy[7] * G({ 0, e / d }, 1)
            + sy[2]
                * (-sy[3] - sy[4] - G({ 0, c / d, e / d }, 1) - G({ 0, e / d, c / d }, 1)
                    - G({ c / d, 0, e / d }, 1) - G({ e / d, 0, c / d }, 1))
            + sy[14] * G({ c / d, b / d, 0, e / d }, 1)
            + G({ 0, c / d, b / d, 0, a / d, e / d }, 1)
            + G({ 0, c / d, b / d, 0, e / d, a / d }, 1)
            + G({ 0, c / d, b / d, e / d, 0, a / d }, 1)
            + G({ 0, c / d, e / d, b / d, 0, a / d }, 1)
            + G({ 0, e / d, c / d, b / d, 0, a / d }, 1)
            + G({ c / d, 0, b / d, 0, a / d, e / d }, 1)
            + G({ c / d, 0, b / d, 0, e / d, a / d }, 1)
            + G({ c / d, 0, b / d, e / d, 0, a / d }, 1)
            + G({ c / d, 0, e / d, b / d, 0, a / d }, 1)
            + 2. * G({ c / d, b / d, 0, 0, a / d, e / d }, 1)
            + 2. * G({ c / d, b / d, 0, 0, e / d, a / d }, 1)
            + G({ c / d, b / d, 0, a / d, 0, e / d }, 1)
            + G({ c / d, b / d, 0, a / d, e / d, x / d }, 1)
            + 2. * G({ c / d, b / d, 0, e / d, 0, a / d }, 1)
            + G({ c / d, b / d, 0, e / d, a / d, x / d }, 1)
            + 2. * G({ c / d, b / d, e / d, 0, 0, a / d }, 1)
            + G({ c / d, b / d, e / d, 0, a / d, x / d }, 1)
            + G({ c / d, e / d, 0, b / d, 0, a / d }, 1)
            + 2. * G({ c / d, e / d, b / d, 0, 0, a / d }, 1)
            + G({ c / d, e / d, b / d, 0, a / d, x / d }, 1)
            + G({ e / d, 0, c / d, b / d, 0, a / d }, 1)
            + G({ e / d, c / d, 0, b / d, 0, a / d }, 1)
            + 2. * G({ e / d, c / d, b / d, 0, 0, a / d }, 1)
            + G({ e / d, c / d, b / d, 0, a / d, x / d }, 1)
            + sy[13] * G({ 0, e }, { 1, se }, x)
            + sy[1] * (sy[0] - G({ a, 0, e }, { sa, 1, se }, x))
            + G({ c / d, b / d }, 1) * (-sy[11] + G({ a, 0, 0, e }, { sa, 1, 1, se }, x))
            + sy[24]
                * (-sy[25] - G({ 0, a, 0, b, c }, { 1, sa, 1, sb, sc }, x)
                    - 2. * G({ a, 0, 0, b, c }, { sa, 1, 1, sb, sc }, x)
                    - G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x))
            + G({ c / d }, 1) * (sy[25] - G({ a, 0, b, 0, e }, { sa, 1, sb, 1, se }, x))
            + G({ a, 0, b, c, 0, e }, { sa, 1, sb, sc, 1, se }, x)
            + (-(sy[15] * sy[1]) + sy[26] + sy[27] + sy[28] - sy[12] * sy[2] + sy[29]
                  + sy[30] + sy[15] * (-sy[3] - sy[4])
                  + (sy[16] + sy[17] + sy[18] + sy[19] + sy[20] + sy[21] + sy[22]
                        + sy[23])
                      * sy[6]
                  + (-sy[1] - sy[3] - sy[4]) * sy[5] * sy[6] + sy[24] * sy[7]
                  - sy[2] * sy[8])
                * Log(-x, sd)
            + (sy[1] / 2. + sy[3] / 2. + sy[4] / 2.) * sy[6] * pow(sy[5], 2.)
            + sy[6]
                * (-sy[27] - sy[28] - sy[29] - sy[30] + sy[31] * (-sy[3] - sy[4])
                    + G({ 0, 0, c / d, b / d, e / d }, 1)
                    + G({ 0, 0, c / d, e / d, b / d }, 1)
                    + G({ 0, 0, e / d, c / d, b / d }, 1)
                    + G({ 0, c / d, 0, b / d, e / d }, 1)
                    + G({ 0, c / d, 0, e / d, b / d }, 1)
                    + G({ 0, c / d, e / d, 0, b / d }, 1)
                    + G({ 0, e / d, 0, c / d, b / d }, 1)
                    + G({ 0, e / d, c / d, 0, b / d }, 1)
                    + G({ c / d, 0, 0, b / d, e / d }, 1)
                    + G({ c / d, 0, 0, e / d, b / d }, 1)
                    + G({ c / d, 0, e / d, 0, b / d }, 1)
                    - G({ c / d, b / d, 0, 0, e / d }, 1)
                    + G({ c / d, e / d, 0, 0, b / d }, 1)
                    + G({ e / d, 0, 0, c / d, b / d }, 1)
                    + G({ e / d, 0, c / d, 0, b / d }, 1)
                    + G({ e / d, c / d, 0, 0, b / d }, 1)
                    + sy[1] * (-sy[31] - 2. * Zeta(2)) - 2. * sy[3] * Zeta(2)
                    - 2. * sy[4] * Zeta(2)) };
        if (e != x) {
            res += (-sy[26] + G({ c / d, b / d, 0, a / d, x / d }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_a0bcde_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (c == d) {
        if (c == e) {
            const vector<complex<double>> sy = { G({ a, c, c }, { sa, sc, sc }, x),
                G({ b / c, 0, 1 }, 1), G({ b / c, 0, a / c }, 1), G({ b / c, 1, 1 }, 1),
                G({ b / c, 0, a / c, 1 }, 1), G({ c, c }, { sc, sc }, x),
                G({ c }, { sc }, x), G({ b / c, 0, a / c, 1, 1 }, 1) };
            complex<double> res { -(sy[0] * sy[1]) + sy[0] * sy[2] + sy[4] * sy[5]
                - sy[6] * sy[7] - sy[5] * G({ b / c, 0, a / c, x / c }, 1)
                + sy[6] * G({ b / c, 0, a / c, x / c, 1 }, 1)
                + G({ b / c, 0, a / c, 0, 1, 1 }, 1)
                - G({ b / c, 0, a / c, x / c, 1, 1 }, 1)
                + (sy[7] - G({ b / c, 0, 0, 1, 1 }, 1)) * G({ a }, { sa }, x)
                + (-sy[4] + G({ b / c, 0, 1, 1 }, 1)) * G({ a, c }, { sa, sc }, x)
                - sy[2] * G({ 0, c, c }, { 1, sc, sc }, x)
                + (sy[1] - sy[3]) * G({ a, 0, c }, { sa, 1, sc }, x)
                + G({ b / c, 1 }, 1)
                    * (-G({ a, 0, b, c }, { sa, 1, sb, sc }, x)
                        + G({ a, 0, c, c }, { sa, 1, sc, sc }, x))
                + G({ b / c }, 1)
                    * (-G({ a, 0, 0, c, c }, { sa, 1, 1, sc, sc }, x)
                        + G({ a, 0, b, c, c }, { sa, 1, sb, sc, sc }, x))
                + G({ a, 0, b, 0, c, c }, { sa, 1, sb, 1, sc, sc }, x)
                + G({ a, 0, b }, { sa, 1, sb }, x) * (sy[3] - Zeta(3)) };
            return res;
        } else {
            const vector<complex<double>> sy = { G({ a, c, e }, { sa, sc, se }, x),
                G({ b / c, 0, 1 }, 1), G({ b / c, 0, a / c }, 1),
                G({ b / c, e / c, 1 }, 1), G({ 0, 0, a }, { 1, 1, sa }, x),
                G({ a, 0, b }, { sa, 1, sb }, x), G({ e / c, 1, b / c }, 1),
                G({ e / c, b / c, 1 }, 1), Log(c, sc), G({ a }, { sa }, x),
                G({ e / c, 1 }, 1), G({ a, 0, b, e }, { sa, 1, sb, se }, x),
                G({ c, e }, { sc, se }, x), G({ b / c, 0, a / c, 1 }, 1),
                G({ 0, a }, { 1, sa }, x), G({ 0, b / c, e / c, 1 }, 1),
                G({ 0, e / c, 1, b / c }, 1), G({ 0, e / c, b / c, 1 }, 1),
                G({ e / c, 0, 1, b / c }, 1), G({ e / c, 0, b / c, 1 }, 1),
                G({ e / c, 1, 0, b / c }, 1), G({ b / c, 0, a / c, e / c, 1 }, 1),
                G({ b / c, 0, e / c, 1, a / c }, 1), G({ b / c, 0, e / c, a / c, 1 }, 1),
                G({ b / c, e / c, 0, 1, a / c }, 1), G({ b / c, e / c, 0, a / c, 1 }, 1),
                G({ b / c, e / c, 1, 0, a / c }, 1), G({ e / c, 1, b / c, 0, a / c }, 1),
                G({ e / c, b / c, 0, 1, a / c }, 1), G({ e / c, b / c, 0, a / c, 1 }, 1),
                G({ e / c, b / c, 1, 0, a / c }, 1), G({ 0, 0 }, { 1, 1 }, x) };
            complex<double> res { sy[12] * sy[13] - sy[0] * sy[1]
                + sy[14] * (-sy[15] - sy[16] - sy[17] - sy[18] - sy[19] - sy[20])
                + sy[0] * sy[2] + sy[3] * sy[4] + sy[4] * (sy[6] + sy[7])
                + (sy[9] * (-sy[15] - sy[16] - sy[17] - sy[18] - sy[19] - sy[20]) - sy[21]
                      - sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27] - sy[28]
                      - sy[29] - sy[30] + sy[14] * sy[3] + sy[10] * sy[5]
                      + sy[14] * (sy[6] + sy[7]))
                    * sy[8]
                + sy[5] * (-sy[6] - sy[7] - G({ 0, e / c, 1 }, 1))
                - sy[12] * G({ b / c, 0, a / c, x / c }, 1)
                + G({ 0, b / c, 0, a / c, e / c, 1 }, 1)
                + G({ 0, b / c, 0, e / c, 1, a / c }, 1)
                + G({ 0, b / c, 0, e / c, a / c, 1 }, 1)
                + G({ 0, b / c, e / c, 0, 1, a / c }, 1)
                + G({ 0, b / c, e / c, 0, a / c, 1 }, 1)
                + G({ 0, b / c, e / c, 1, 0, a / c }, 1)
                + G({ 0, e / c, 1, b / c, 0, a / c }, 1)
                + G({ 0, e / c, b / c, 0, 1, a / c }, 1)
                + G({ 0, e / c, b / c, 0, a / c, 1 }, 1)
                + G({ 0, e / c, b / c, 1, 0, a / c }, 1)
                + 2. * G({ b / c, 0, 0, a / c, e / c, 1 }, 1)
                + 2. * G({ b / c, 0, 0, e / c, 1, a / c }, 1)
                + 2. * G({ b / c, 0, 0, e / c, a / c, 1 }, 1)
                + G({ b / c, 0, a / c, 0, e / c, 1 }, 1)
                + G({ b / c, 0, a / c, e / c, 1, x / c }, 1)
                + G({ b / c, 0, a / c, e / c, x / c, 1 }, 1)
                + 2. * G({ b / c, 0, e / c, 0, 1, a / c }, 1)
                + 2. * G({ b / c, 0, e / c, 0, a / c, 1 }, 1)
                + 2. * G({ b / c, 0, e / c, 1, 0, a / c }, 1)
                + G({ b / c, 0, e / c, 1, a / c, x / c }, 1)
                + G({ b / c, 0, e / c, a / c, 1, x / c }, 1)
                + G({ b / c, 0, e / c, a / c, x / c, 1 }, 1)
                + 2. * G({ b / c, e / c, 0, 0, 1, a / c }, 1)
                + 2. * G({ b / c, e / c, 0, 0, a / c, 1 }, 1)
                + 2. * G({ b / c, e / c, 0, 1, 0, a / c }, 1)
                + G({ b / c, e / c, 0, 1, a / c, x / c }, 1)
                + G({ b / c, e / c, 0, a / c, 1, x / c }, 1)
                + G({ b / c, e / c, 0, a / c, x / c, 1 }, 1)
                + 2. * G({ b / c, e / c, 1, 0, 0, a / c }, 1)
                + G({ b / c, e / c, 1, 0, a / c, x / c }, 1)
                + G({ e / c, 0, 1, b / c, 0, a / c }, 1)
                + G({ e / c, 0, b / c, 0, 1, a / c }, 1)
                + G({ e / c, 0, b / c, 0, a / c, 1 }, 1)
                + G({ e / c, 0, b / c, 1, 0, a / c }, 1)
                + G({ e / c, 1, 0, b / c, 0, a / c }, 1)
                + 2. * G({ e / c, 1, b / c, 0, 0, a / c }, 1)
                + G({ e / c, 1, b / c, 0, a / c, x / c }, 1)
                + 2. * G({ e / c, b / c, 0, 0, 1, a / c }, 1)
                + 2. * G({ e / c, b / c, 0, 0, a / c, 1 }, 1)
                + 2. * G({ e / c, b / c, 0, 1, 0, a / c }, 1)
                + G({ e / c, b / c, 0, 1, a / c, x / c }, 1)
                + G({ e / c, b / c, 0, a / c, 1, x / c }, 1)
                + G({ e / c, b / c, 0, a / c, x / c, 1 }, 1)
                + 2. * G({ e / c, b / c, 1, 0, 0, a / c }, 1)
                + G({ e / c, b / c, 1, 0, a / c, x / c }, 1)
                + (-sy[13] + G({ b / c, 0, e / c, 1 }, 1)) * G({ a, e }, { sa, se }, x)
                - sy[2] * G({ 0, c, e }, { 1, sc, se }, x)
                + (sy[1] - sy[3]) * G({ a, 0, e }, { sa, 1, se }, x)
                + sy[10]
                    * (sy[11] + G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                        + 2. * G({ a, 0, 0, b }, { sa, 1, 1, sb }, x))
                + G({ b / c, 1 }, 1) * (-sy[11] + G({ a, 0, c, e }, { sa, 1, sc, se }, x))
                + G({ b / c }, 1)
                    * (-G({ a, 0, 0, c, e }, { sa, 1, 1, sc, se }, x)
                        + G({ a, 0, b, c, e }, { sa, 1, sb, sc, se }, x))
                + G({ a, 0, b, 0, c, e }, { sa, 1, sb, 1, sc, se }, x)
                + (sy[9] * (sy[15] + sy[16] + sy[17] + sy[18] + sy[19] + sy[20]) + sy[21]
                      + sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27] + sy[28]
                      + sy[29] + sy[30] - sy[14] * sy[3] - sy[10] * sy[5]
                      + sy[14] * (-sy[6] - sy[7])
                      + sy[9] * (-sy[3] - sy[6] - sy[7]) * sy[8])
                    * Log(-x, sc)
                + sy[9] * (sy[3] / 2. + sy[6] / 2. + sy[7] / 2.) * pow(sy[8], 2.)
                + sy[9]
                    * (-sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27] - sy[28]
                        - sy[29] - sy[30] + sy[31] * (-sy[6] - sy[7])
                        + G({ 0, 0, b / c, e / c, 1 }, 1)
                        + G({ 0, 0, e / c, 1, b / c }, 1)
                        + G({ 0, 0, e / c, b / c, 1 }, 1)
                        + G({ 0, e / c, 0, 1, b / c }, 1)
                        + G({ 0, e / c, 0, b / c, 1 }, 1)
                        + G({ 0, e / c, 1, 0, b / c }, 1)
                        - G({ b / c, 0, 0, e / c, 1 }, 1)
                        + G({ e / c, 0, 0, 1, b / c }, 1)
                        + G({ e / c, 0, 0, b / c, 1 }, 1)
                        + G({ e / c, 0, 1, 0, b / c }, 1)
                        + G({ e / c, 1, 0, 0, b / c }, 1)
                        + sy[3] * (-sy[31] - 2. * Zeta(2)) - 2. * sy[6] * Zeta(2)
                        - 2. * sy[7] * Zeta(2)) };
            if (e != x) {
                res += (-sy[21] + G({ b / c, 0, a / c, x / c, 1 }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else {
        const vector<complex<double>> sy = { G({ b / c, 0, a / c }, 1),
            G({ a, d, e }, { sa, sd, se }, x), G({ b / c, 0, d / c }, 1),
            G({ 0, 0, a }, { 1, 1, sa }, x), G({ b / c, d / c, e / c }, 1),
            G({ d / c, b / c, e / c }, 1), G({ d / c, e / c, b / c }, 1), Log(c, sc),
            G({ a }, { sa }, x), G({ a, 0, b }, { sa, 1, sb }, x), G({ d / c, b / c }, 1),
            G({ a, 0, 0, e }, { sa, 1, 1, se }, x), G({ d / c, e / c }, 1),
            G({ a, 0, b, e }, { sa, 1, sb, se }, x), G({ a, e }, { sa, se }, x),
            G({ b / c, 0, d / c, a / c }, 1), G({ b / c, 0, a / c, d / c }, 1),
            G({ b / c, d / c, 0, a / c }, 1), G({ d / c, b / c, 0, a / c }, 1),
            G({ 0, a }, { 1, sa }, x), G({ 0, b / c, d / c, e / c }, 1),
            G({ 0, d / c, b / c, e / c }, 1), G({ 0, d / c, e / c, b / c }, 1),
            G({ d / c, 0, b / c, e / c }, 1), G({ d / c, 0, e / c, b / c }, 1),
            G({ d / c, e / c, 0, b / c }, 1),
            G({ a, 0, b, d, e }, { sa, 1, sb, sd, se }, x),
            G({ b / c, 0, a / c, d / c, e / c }, 1),
            G({ b / c, 0, d / c, a / c, e / c }, 1),
            G({ b / c, 0, d / c, e / c, a / c }, 1),
            G({ b / c, d / c, 0, a / c, e / c }, 1),
            G({ b / c, d / c, 0, e / c, a / c }, 1),
            G({ b / c, d / c, e / c, 0, a / c }, 1),
            G({ d / c, b / c, 0, a / c, e / c }, 1),
            G({ d / c, b / c, 0, e / c, a / c }, 1),
            G({ d / c, b / c, e / c, 0, a / c }, 1),
            G({ d / c, e / c, b / c, 0, a / c }, 1), G({ 0, 0 }, { 1, 1 }, x) };
        complex<double> res { -(sy[10] * sy[11]) + sy[10] * sy[13] + sy[14] * sy[15]
            + sy[0] * sy[1]
            + sy[19] * (sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25])
            - sy[1] * sy[2] - sy[3] * sy[4] + sy[3] * (-sy[5] - sy[6])
            + sy[7]
                * (-(sy[9] * sy[12]) + sy[27] + sy[28] + sy[29] + sy[30] + sy[31] + sy[32]
                    + sy[33] + sy[34] + sy[35] + sy[36] - sy[19] * sy[4]
                    + sy[19] * (-sy[5] - sy[6])
                    + (sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25]) * sy[8])
            + sy[9] * (sy[6] + G({ 0, d / c, e / c }, 1) + G({ d / c, 0, e / c }, 1))
            + sy[14]
                * (sy[17] + sy[18] - G({ b / c, 0, d / c, e / c }, 1)
                    - G({ b / c, d / c, 0, e / c }, 1) - G({ d / c, b / c, 0, e / c }, 1))
            - G({ 0, b / c, 0, a / c, d / c, e / c }, 1)
            - G({ 0, b / c, 0, d / c, a / c, e / c }, 1)
            - G({ 0, b / c, 0, d / c, e / c, a / c }, 1)
            - G({ 0, b / c, d / c, 0, a / c, e / c }, 1)
            - G({ 0, b / c, d / c, 0, e / c, a / c }, 1)
            - G({ 0, b / c, d / c, e / c, 0, a / c }, 1)
            - G({ 0, d / c, b / c, 0, a / c, e / c }, 1)
            - G({ 0, d / c, b / c, 0, e / c, a / c }, 1)
            - G({ 0, d / c, b / c, e / c, 0, a / c }, 1)
            - G({ 0, d / c, e / c, b / c, 0, a / c }, 1)
            - 2. * G({ b / c, 0, 0, a / c, d / c, e / c }, 1)
            - 2. * G({ b / c, 0, 0, d / c, a / c, e / c }, 1)
            - 2. * G({ b / c, 0, 0, d / c, e / c, a / c }, 1)
            - G({ b / c, 0, a / c, 0, d / c, e / c }, 1)
            - G({ b / c, 0, a / c, d / c, 0, e / c }, 1)
            - G({ b / c, 0, a / c, d / c, e / c, x / c }, 1)
            - 2. * G({ b / c, 0, d / c, 0, a / c, e / c }, 1)
            - 2. * G({ b / c, 0, d / c, 0, e / c, a / c }, 1)
            - G({ b / c, 0, d / c, a / c, 0, e / c }, 1)
            - G({ b / c, 0, d / c, a / c, e / c, x / c }, 1)
            - 2. * G({ b / c, 0, d / c, e / c, 0, a / c }, 1)
            - G({ b / c, 0, d / c, e / c, a / c, x / c }, 1)
            - 2. * G({ b / c, d / c, 0, 0, a / c, e / c }, 1)
            - 2. * G({ b / c, d / c, 0, 0, e / c, a / c }, 1)
            - G({ b / c, d / c, 0, a / c, 0, e / c }, 1)
            - G({ b / c, d / c, 0, a / c, e / c, x / c }, 1)
            - 2. * G({ b / c, d / c, 0, e / c, 0, a / c }, 1)
            - G({ b / c, d / c, 0, e / c, a / c, x / c }, 1)
            - 2. * G({ b / c, d / c, e / c, 0, 0, a / c }, 1)
            - G({ b / c, d / c, e / c, 0, a / c, x / c }, 1)
            - G({ d / c, 0, b / c, 0, a / c, e / c }, 1)
            - G({ d / c, 0, b / c, 0, e / c, a / c }, 1)
            - G({ d / c, 0, b / c, e / c, 0, a / c }, 1)
            - G({ d / c, 0, e / c, b / c, 0, a / c }, 1)
            - 2. * G({ d / c, b / c, 0, 0, a / c, e / c }, 1)
            - 2. * G({ d / c, b / c, 0, 0, e / c, a / c }, 1)
            - G({ d / c, b / c, 0, a / c, 0, e / c }, 1)
            - G({ d / c, b / c, 0, a / c, e / c, x / c }, 1)
            - 2. * G({ d / c, b / c, 0, e / c, 0, a / c }, 1)
            - G({ d / c, b / c, 0, e / c, a / c, x / c }, 1)
            - 2. * G({ d / c, b / c, e / c, 0, 0, a / c }, 1)
            - G({ d / c, b / c, e / c, 0, a / c, x / c }, 1)
            - G({ d / c, e / c, 0, b / c, 0, a / c }, 1)
            - 2. * G({ d / c, e / c, b / c, 0, 0, a / c }, 1)
            - G({ d / c, e / c, b / c, 0, a / c, x / c }, 1)
            + (-sy[15] - sy[16] - sy[17] - sy[18]) * G({ 0, e }, { 1, se }, x)
            - sy[0] * G({ 0, d, e }, { 1, sd, se }, x)
            + (sy[2] + sy[4] + sy[5]) * G({ a, 0, e }, { sa, 1, se }, x)
            + sy[12]
                * (-sy[13] - G({ 0, a, 0, b }, { 1, sa, 1, sb }, x)
                    - 2. * G({ a, 0, 0, b }, { sa, 1, 1, sb }, x))
            + G({ b / c, d / c }, 1) * (-sy[11] + G({ a, 0, d, e }, { sa, 1, sd, se }, x))
            + G({ b / c }, 1) * (sy[26] - G({ a, 0, 0, d, e }, { sa, 1, 1, sd, se }, x))
            + G({ d / c }, 1) * (-sy[26] + G({ a, 0, b, 0, e }, { sa, 1, sb, 1, se }, x))
            + G({ a, 0, b, 0, d, e }, { sa, 1, sb, 1, sd, se }, x)
            + (sy[9] * sy[12] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31] - sy[32]
                  - sy[33] - sy[34] - sy[35] - sy[36] + sy[19] * sy[4]
                  + sy[19] * (sy[5] + sy[6])
                  + (-sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25]) * sy[8]
                  + (sy[4] + sy[5] + sy[6]) * sy[7] * sy[8])
                * Log(-x, sc)
            + (-sy[4] / 2. - sy[5] / 2. - sy[6] / 2.) * sy[8] * pow(sy[7], 2.)
            + sy[8]
                * (sy[29] + sy[31] + sy[32] + sy[34] + sy[35] + sy[36]
                    + sy[37] * (sy[5] + sy[6]) - G({ 0, 0, b / c, d / c, e / c }, 1)
                    - G({ 0, 0, d / c, b / c, e / c }, 1)
                    - G({ 0, 0, d / c, e / c, b / c }, 1)
                    - G({ 0, d / c, 0, b / c, e / c }, 1)
                    - G({ 0, d / c, 0, e / c, b / c }, 1)
                    - G({ 0, d / c, e / c, 0, b / c }, 1)
                    + G({ b / c, 0, 0, d / c, e / c }, 1)
                    + G({ b / c, 0, d / c, 0, e / c }, 1)
                    + G({ b / c, d / c, 0, 0, e / c }, 1)
                    - G({ d / c, 0, 0, b / c, e / c }, 1)
                    - G({ d / c, 0, 0, e / c, b / c }, 1)
                    - G({ d / c, 0, e / c, 0, b / c }, 1)
                    + G({ d / c, b / c, 0, 0, e / c }, 1)
                    - G({ d / c, e / c, 0, 0, b / c }, 1) + 2. * sy[5] * Zeta(2)
                    + 2. * sy[6] * Zeta(2) + sy[4] * (sy[37] + 2. * Zeta(2))) };
        if (d != x) {
            res += (sy[16] - G({ b / c, 0, a / c, x / c }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (sy[27] + sy[28] + sy[30] + sy[33]
                       - G({ b / c, 0, a / c, d / c, x / c }, 1)
                       - G({ b / c, 0, d / c, a / c, x / c }, 1)
                       - G({ b / c, d / c, 0, a / c, x / c }, 1)
                       - G({ d / c, b / c, 0, a / c, x / c }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_a0bcde_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (b == c) {
        if (b == d) {
            if (b == e) { // abbbb
                const vector<complex<double>> sy = { G({ a, b, b }, { sa, sb, sb }, x),
                    G({ b, b, b }, { sb, sb, sb }, x), G({ a, b }, { sa, sb }, x),
                    G({ 0, a / b, 1, 1 }, 1), G({ b, b }, { sb, sb }, x),
                    G({ 0, a / b }, 1), G({ a, b, b, b }, { sa, sb, sb, sb }, x),
                    G({ a }, { sa }, x), G({ 0, a / b, 1, 1, 1 }, 1),
                    G({ b }, { sb }, x) };
                complex<double> res { -(sy[2] * sy[3]) + sy[3] * sy[4] - sy[5] * sy[6]
                    - sy[9] * sy[8] + sy[7] * sy[8]
                    + (sy[0] - sy[1]) * G({ 0, a / b, 1 }, 1)
                    + sy[1] * G({ 0, a / b, x / b }, 1)
                    - sy[4] * G({ 0, a / b, x / b, 1 }, 1)
                    + sy[9] * G({ 0, a / b, x / b, 1, 1 }, 1)
                    + G({ 0, a / b, 0, 1, 1, 1 }, 1) - G({ 0, a / b, x / b, 1, 1, 1 }, 1)
                    + sy[5] * G({ 0, b, b, b }, { 1, sb, sb, sb }, x)
                    + G({ a, 0, 0, b, b, b }, { sa, 1, 1, sb, sb, sb }, x)
                    - sy[6] * Zeta(2) + G({ a, 0, b, b }, { sa, 1, sb, sb }, x) * Zeta(2)
                    - sy[2] * Zeta(4)
                    + (-sy[0] + G({ a, 0, b }, { sa, 1, sb }, x)) * Zeta(3)
                    - sy[7] * (Zeta(2) * Zeta(3) - 2. * Zeta(5)) };
                return res;
            } else { // abbbe
                const vector<complex<double>> sy = { G({ a, b, e }, { sa, sb, se }, x),
                    G({ b, b, e }, { sb, sb, se }, x), Log(b, sb), G({ a }, { sa }, x),
                    G({ e / b, 1, 1 }, 1), G({ a, 0, e }, { sa, 1, se }, x),
                    G({ b, e }, { sb, se }, x), G({ 0, a / b, 1, 1 }, 1),
                    G({ a, b, b, e }, { sa, sb, sb, se }, x), G({ 0, a }, { 1, sa }, x),
                    G({ 0, a / b, e / b, 1, 1 }, 1), G({ 0, e / b, 1, 1, a / b }, 1),
                    G({ 0, e / b, 1, a / b, 1 }, 1), G({ 0, e / b, a / b, 1, 1 }, 1),
                    G({ e / b, 0, 1, 1, a / b }, 1), G({ e / b, 0, 1, a / b, 1 }, 1),
                    G({ e / b, 0, a / b, 1, 1 }, 1), G({ e / b, 1, 0, 1, a / b }, 1),
                    G({ e / b, 1, 0, a / b, 1 }, 1), G({ e / b, 1, 1, 0, a / b }, 1) };
                complex<double> res { sy[2]
                        * (-sy[10] - sy[11] - sy[12] - sy[13] - sy[14] - sy[15] - sy[16]
                            - sy[17] - sy[18] - sy[19] + sy[9] * sy[4])
                    + sy[6] * sy[7] + (sy[0] - sy[1]) * G({ 0, a / b, 1 }, 1)
                    + sy[1] * G({ 0, a / b, x / b }, 1)
                    - sy[6] * G({ 0, a / b, x / b, 1 }, 1)
                    + 2. * G({ 0, 0, a / b, e / b, 1, 1 }, 1)
                    + 2. * G({ 0, 0, e / b, 1, 1, a / b }, 1)
                    + 2. * G({ 0, 0, e / b, 1, a / b, 1 }, 1)
                    + 2. * G({ 0, 0, e / b, a / b, 1, 1 }, 1)
                    + G({ 0, a / b, 0, e / b, 1, 1 }, 1)
                    + G({ 0, a / b, e / b, 1, 1, x / b }, 1)
                    + G({ 0, a / b, e / b, 1, x / b, 1 }, 1)
                    + G({ 0, a / b, e / b, x / b, 1, 1 }, 1)
                    + 2. * G({ 0, e / b, 0, 1, 1, a / b }, 1)
                    + 2. * G({ 0, e / b, 0, 1, a / b, 1 }, 1)
                    + 2. * G({ 0, e / b, 0, a / b, 1, 1 }, 1)
                    + 2. * G({ 0, e / b, 1, 0, 1, a / b }, 1)
                    + 2. * G({ 0, e / b, 1, 0, a / b, 1 }, 1)
                    + 2. * G({ 0, e / b, 1, 1, 0, a / b }, 1)
                    + G({ 0, e / b, 1, 1, a / b, x / b }, 1)
                    + G({ 0, e / b, 1, a / b, 1, x / b }, 1)
                    + G({ 0, e / b, 1, a / b, x / b, 1 }, 1)
                    + G({ 0, e / b, a / b, 1, 1, x / b }, 1)
                    + G({ 0, e / b, a / b, 1, x / b, 1 }, 1)
                    + G({ 0, e / b, a / b, x / b, 1, 1 }, 1)
                    + 2. * G({ e / b, 0, 0, 1, 1, a / b }, 1)
                    + 2. * G({ e / b, 0, 0, 1, a / b, 1 }, 1)
                    + 2. * G({ e / b, 0, 0, a / b, 1, 1 }, 1)
                    + 2. * G({ e / b, 0, 1, 0, 1, a / b }, 1)
                    + 2. * G({ e / b, 0, 1, 0, a / b, 1 }, 1)
                    + 2. * G({ e / b, 0, 1, 1, 0, a / b }, 1)
                    + G({ e / b, 0, 1, 1, a / b, x / b }, 1)
                    + G({ e / b, 0, 1, a / b, 1, x / b }, 1)
                    + G({ e / b, 0, 1, a / b, x / b, 1 }, 1)
                    + G({ e / b, 0, a / b, 1, 1, x / b }, 1)
                    + G({ e / b, 0, a / b, 1, x / b, 1 }, 1)
                    + G({ e / b, 0, a / b, x / b, 1, 1 }, 1)
                    + 2. * G({ e / b, 1, 0, 0, 1, a / b }, 1)
                    + 2. * G({ e / b, 1, 0, 0, a / b, 1 }, 1)
                    + 2. * G({ e / b, 1, 0, 1, 0, a / b }, 1)
                    + G({ e / b, 1, 0, 1, a / b, x / b }, 1)
                    + G({ e / b, 1, 0, a / b, 1, x / b }, 1)
                    + G({ e / b, 1, 0, a / b, x / b, 1 }, 1)
                    + 2. * G({ e / b, 1, 1, 0, 0, a / b }, 1)
                    + G({ e / b, 1, 1, 0, a / b, x / b }, 1)
                    + (-sy[7] + G({ 0, e / b, 1, 1 }, 1)) * G({ a, e }, { sa, se }, x)
                    + sy[4] * (-sy[5] + G({ 0, 0, a }, { 1, 1, sa }, x))
                    + G({ 0, a / b }, 1)
                        * (-sy[8] + G({ 0, b, b, e }, { 1, sb, sb, se }, x))
                    + G({ a, 0, 0, b, b, e }, { sa, 1, 1, sb, sb, se }, x)
                    + (sy[10] + sy[11] + sy[12] + sy[13] + sy[14] + sy[15] + sy[16]
                          + sy[17] + sy[18] + sy[19] - sy[9] * sy[4]
                          - sy[2] * sy[3] * sy[4])
                        * Log(-x, sb)
                    + (sy[3] * sy[4] * pow(sy[2], 2.)) / 2.
                    + sy[3]
                        * (-sy[11] - sy[12] - sy[13] - sy[14] - sy[15] - sy[16] - sy[17]
                            - sy[18] - sy[19] - G({ 0, 0, e / b, 1, 1 }, 1)
                            + sy[4] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2)))
                    - sy[8] * Zeta(2) + G({ a, 0, b, e }, { sa, 1, sb, se }, x) * Zeta(2)
                    - sy[0] * Zeta(3) + sy[5] * Zeta(3) };

                if (e != x) {
                    res += (-sy[10] + G({ 0, a / b, x / b, 1, 1 }, 1))
                        * G({ e }, { se }, x);
                }
                return res;
            }
        } else { // abbde
            const vector<complex<double>> sy = { G({ a, 0, e }, { sa, 1, se }, x),
                G({ a, d, e }, { sa, sd, se }, x), G({ b, d, e }, { sb, sd, se }, x),
                G({ 0, 0, a }, { 1, 1, sa }, x), G({ d / b, 1, e / b }, 1), Log(b, sb),
                G({ a }, { sa }, x), G({ d / b, e / b, 1 }, 1),
                G({ a, e }, { sa, se }, x), G({ 0, d / b, 1, a / b }, 1),
                G({ a, 0, d, e }, { sa, 1, sd, se }, x),
                G({ a, b, d, e }, { sa, sb, sd, se }, x), G({ 0, a / b, d / b, 1 }, 1),
                G({ 0, d / b, a / b, 1 }, 1), G({ d / b, 0, 1, a / b }, 1),
                G({ d / b, 0, a / b, 1 }, 1), G({ d / b, 1, 0, a / b }, 1),
                G({ 0, a }, { 1, sa }, x), G({ 0, a / b, d / b, 1, e / b }, 1),
                G({ 0, a / b, d / b, e / b, 1 }, 1), G({ 0, d / b, 1, a / b, e / b }, 1),
                G({ 0, d / b, 1, e / b, a / b }, 1), G({ 0, d / b, a / b, 1, e / b }, 1),
                G({ 0, d / b, a / b, e / b, 1 }, 1), G({ 0, d / b, e / b, 1, a / b }, 1),
                G({ 0, d / b, e / b, a / b, 1 }, 1), G({ d / b, 0, 1, a / b, e / b }, 1),
                G({ d / b, 0, 1, e / b, a / b }, 1), G({ d / b, 0, a / b, 1, e / b }, 1),
                G({ d / b, 0, a / b, e / b, 1 }, 1), G({ d / b, 0, e / b, 1, a / b }, 1),
                G({ d / b, 0, e / b, a / b, 1 }, 1), G({ d / b, 1, 0, a / b, e / b }, 1),
                G({ d / b, 1, 0, e / b, a / b }, 1), G({ d / b, 1, e / b, 0, a / b }, 1),
                G({ d / b, e / b, 0, 1, a / b }, 1), G({ d / b, e / b, 0, a / b, 1 }, 1),
                G({ d / b, e / b, 1, 0, a / b }, 1), G({ 0, 0 }, { 1, 1 }, x) };
            complex<double> res { -(sy[3] * sy[4]) - sy[3] * sy[7]
                + sy[0] * (sy[4] + sy[7])
                + sy[5]
                    * (sy[18] + sy[19] + sy[20] + sy[21] + sy[22] + sy[23] + sy[24]
                        + sy[25] + sy[26] + sy[27] + sy[28] + sy[29] + sy[30] + sy[31]
                        + sy[32] + sy[33] + sy[34] + sy[35] + sy[36] + sy[37]
                        - sy[17] * sy[4] - sy[17] * sy[7])
                + sy[9] * sy[8] + (sy[1] - sy[2]) * G({ 0, a / b, 1 }, 1)
                + sy[2] * G({ 0, a / b, x / b }, 1)
                + (sy[0] - sy[1]) * G({ 0, d / b, 1 }, 1)
                + sy[8]
                    * (sy[13] + sy[14] + sy[15] + sy[16] - G({ 0, d / b, 1, e / b }, 1)
                        - G({ 0, d / b, e / b, 1 }, 1) - G({ d / b, 0, 1, e / b }, 1)
                        - G({ d / b, 0, e / b, 1 }, 1) - G({ d / b, 1, 0, e / b }, 1))
                - 2. * G({ 0, 0, a / b, d / b, 1, e / b }, 1)
                - 2. * G({ 0, 0, a / b, d / b, e / b, 1 }, 1)
                - 2. * G({ 0, 0, d / b, 1, a / b, e / b }, 1)
                - 2. * G({ 0, 0, d / b, 1, e / b, a / b }, 1)
                - 2. * G({ 0, 0, d / b, a / b, 1, e / b }, 1)
                - 2. * G({ 0, 0, d / b, a / b, e / b, 1 }, 1)
                - 2. * G({ 0, 0, d / b, e / b, 1, a / b }, 1)
                - 2. * G({ 0, 0, d / b, e / b, a / b, 1 }, 1)
                - G({ 0, a / b, 0, d / b, 1, e / b }, 1)
                - G({ 0, a / b, 0, d / b, e / b, 1 }, 1)
                - G({ 0, a / b, d / b, 0, 1, e / b }, 1)
                - G({ 0, a / b, d / b, 0, e / b, 1 }, 1)
                - G({ 0, a / b, d / b, 1, 0, e / b }, 1)
                - G({ 0, a / b, d / b, 1, e / b, x / b }, 1)
                - G({ 0, a / b, d / b, e / b, 1, x / b }, 1)
                - G({ 0, a / b, d / b, e / b, x / b, 1 }, 1)
                - 2. * G({ 0, d / b, 0, 1, a / b, e / b }, 1)
                - 2. * G({ 0, d / b, 0, 1, e / b, a / b }, 1)
                - 2. * G({ 0, d / b, 0, a / b, 1, e / b }, 1)
                - 2. * G({ 0, d / b, 0, a / b, e / b, 1 }, 1)
                - 2. * G({ 0, d / b, 0, e / b, 1, a / b }, 1)
                - 2. * G({ 0, d / b, 0, e / b, a / b, 1 }, 1)
                - 2. * G({ 0, d / b, 1, 0, a / b, e / b }, 1)
                - 2. * G({ 0, d / b, 1, 0, e / b, a / b }, 1)
                - G({ 0, d / b, 1, a / b, 0, e / b }, 1)
                - G({ 0, d / b, 1, a / b, e / b, x / b }, 1)
                - 2. * G({ 0, d / b, 1, e / b, 0, a / b }, 1)
                - G({ 0, d / b, 1, e / b, a / b, x / b }, 1)
                - G({ 0, d / b, a / b, 0, 1, e / b }, 1)
                - G({ 0, d / b, a / b, 0, e / b, 1 }, 1)
                - G({ 0, d / b, a / b, 1, 0, e / b }, 1)
                - G({ 0, d / b, a / b, 1, e / b, x / b }, 1)
                - G({ 0, d / b, a / b, e / b, 1, x / b }, 1)
                - G({ 0, d / b, a / b, e / b, x / b, 1 }, 1)
                - 2. * G({ 0, d / b, e / b, 0, 1, a / b }, 1)
                - 2. * G({ 0, d / b, e / b, 0, a / b, 1 }, 1)
                - 2. * G({ 0, d / b, e / b, 1, 0, a / b }, 1)
                - G({ 0, d / b, e / b, 1, a / b, x / b }, 1)
                - G({ 0, d / b, e / b, a / b, 1, x / b }, 1)
                - G({ 0, d / b, e / b, a / b, x / b, 1 }, 1)
                - 2. * G({ d / b, 0, 0, 1, a / b, e / b }, 1)
                - 2. * G({ d / b, 0, 0, 1, e / b, a / b }, 1)
                - 2. * G({ d / b, 0, 0, a / b, 1, e / b }, 1)
                - 2. * G({ d / b, 0, 0, a / b, e / b, 1 }, 1)
                - 2. * G({ d / b, 0, 0, e / b, 1, a / b }, 1)
                - 2. * G({ d / b, 0, 0, e / b, a / b, 1 }, 1)
                - 2. * G({ d / b, 0, 1, 0, a / b, e / b }, 1)
                - 2. * G({ d / b, 0, 1, 0, e / b, a / b }, 1)
                - G({ d / b, 0, 1, a / b, 0, e / b }, 1)
                - G({ d / b, 0, 1, a / b, e / b, x / b }, 1)
                - 2. * G({ d / b, 0, 1, e / b, 0, a / b }, 1)
                - G({ d / b, 0, 1, e / b, a / b, x / b }, 1)
                - G({ d / b, 0, a / b, 0, 1, e / b }, 1)
                - G({ d / b, 0, a / b, 0, e / b, 1 }, 1)
                - G({ d / b, 0, a / b, 1, 0, e / b }, 1)
                - G({ d / b, 0, a / b, 1, e / b, x / b }, 1)
                - G({ d / b, 0, a / b, e / b, 1, x / b }, 1)
                - G({ d / b, 0, a / b, e / b, x / b, 1 }, 1)
                - 2. * G({ d / b, 0, e / b, 0, 1, a / b }, 1)
                - 2. * G({ d / b, 0, e / b, 0, a / b, 1 }, 1)
                - 2. * G({ d / b, 0, e / b, 1, 0, a / b }, 1)
                - G({ d / b, 0, e / b, 1, a / b, x / b }, 1)
                - G({ d / b, 0, e / b, a / b, 1, x / b }, 1)
                - G({ d / b, 0, e / b, a / b, x / b, 1 }, 1)
                - 2. * G({ d / b, 1, 0, 0, a / b, e / b }, 1)
                - 2. * G({ d / b, 1, 0, 0, e / b, a / b }, 1)
                - G({ d / b, 1, 0, a / b, 0, e / b }, 1)
                - G({ d / b, 1, 0, a / b, e / b, x / b }, 1)
                - 2. * G({ d / b, 1, 0, e / b, 0, a / b }, 1)
                - G({ d / b, 1, 0, e / b, a / b, x / b }, 1)
                - 2. * G({ d / b, 1, e / b, 0, 0, a / b }, 1)
                - G({ d / b, 1, e / b, 0, a / b, x / b }, 1)
                - 2. * G({ d / b, e / b, 0, 0, 1, a / b }, 1)
                - 2. * G({ d / b, e / b, 0, 0, a / b, 1 }, 1)
                - 2. * G({ d / b, e / b, 0, 1, 0, a / b }, 1)
                - G({ d / b, e / b, 0, 1, a / b, x / b }, 1)
                - G({ d / b, e / b, 0, a / b, 1, x / b }, 1)
                - G({ d / b, e / b, 0, a / b, x / b, 1 }, 1)
                - 2. * G({ d / b, e / b, 1, 0, 0, a / b }, 1)
                - G({ d / b, e / b, 1, 0, a / b, x / b }, 1)
                + (-sy[9] - sy[12] - sy[13] - sy[14] - sy[15] - sy[16])
                    * G({ 0, e }, { 1, se }, x)
                + G({ 0, a / b }, 1) * (-sy[11] + G({ 0, b, d, e }, { 1, sb, sd, se }, x))
                + G({ d / b, 1 }, 1) * (sy[10] - G({ a, 0, 0, e }, { sa, 1, 1, se }, x))
                + G({ a, 0, 0, b, d, e }, { sa, 1, 1, sb, sd, se }, x)
                + (-sy[18] - sy[19] - sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25]
                      - sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31] - sy[32]
                      - sy[33] - sy[34] - sy[35] - sy[36] - sy[37] + sy[17] * sy[4]
                      + sy[17] * sy[7] + sy[5] * sy[6] * (sy[4] + sy[7]))
                    * Log(-x, sb)
                + sy[6] * (-sy[4] / 2. - sy[7] / 2.) * pow(sy[5], 2.) + sy[10] * Zeta(2)
                - sy[11] * Zeta(2)
                + sy[6]
                    * (sy[21] + sy[24] + sy[25] + sy[27] + sy[30] + sy[31] + sy[33]
                        + sy[34] + sy[35] + sy[36] + sy[37] + sy[38] * sy[7]
                        + G({ 0, 0, d / b, 1, e / b }, 1)
                        + G({ 0, 0, d / b, e / b, 1 }, 1)
                        + G({ 0, d / b, 0, 1, e / b }, 1)
                        + G({ 0, d / b, 0, e / b, 1 }, 1)
                        + G({ 0, d / b, 1, 0, e / b }, 1)
                        + G({ d / b, 0, 0, 1, e / b }, 1)
                        + G({ d / b, 0, 0, e / b, 1 }, 1)
                        + G({ d / b, 0, 1, 0, e / b }, 1)
                        + G({ d / b, 1, 0, 0, e / b }, 1) + 2. * sy[7] * Zeta(2)
                        + sy[4] * (sy[38] + 2. * Zeta(2))) };
            if (d != x) {
                res += (sy[12] - G({ 0, a / b, x / b, 1 }, 1))
                    * G({ d, e }, { sd, se }, x);
            }
            if (e != x) {
                res += (sy[18] + sy[19] + sy[20] + sy[22] + sy[23] + sy[26] + sy[28]
                           + sy[29] + sy[32] - G({ 0, a / b, d / b, 1, x / b }, 1)
                           - G({ 0, a / b, d / b, x / b, 1 }, 1)
                           - G({ 0, d / b, 1, a / b, x / b }, 1)
                           - G({ 0, d / b, a / b, 1, x / b }, 1)
                           - G({ 0, d / b, a / b, x / b, 1 }, 1)
                           - G({ d / b, 0, 1, a / b, x / b }, 1)
                           - G({ d / b, 0, a / b, 1, x / b }, 1)
                           - G({ d / b, 0, a / b, x / b, 1 }, 1)
                           - G({ d / b, 1, 0, a / b, x / b }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else { // abcde
        const vector<complex<double>> sy = { G({ 0, a / b, c / b }, 1),
            G({ 0, d, e }, { 1, sd, se }, x), G({ 0, c / b, a / b }, 1),
            G({ a, d, e }, { sa, sd, se }, x), G({ a, 0, e }, { sa, 1, se }, x),
            G({ c / b, 0, a / b }, 1), G({ c / b, 0, d / b }, 1),
            G({ c / b, d / b, e / b }, 1), Log(b, sb), G({ a }, { sa }, x),
            G({ a, e }, { sa, se }, x), G({ 0, c / b, d / b, a / b }, 1),
            G({ a, 0, d, e }, { sa, 1, sd, se }, x),
            G({ a, c, d, e }, { sa, sc, sd, se }, x), G({ 0, a / b, c / b, d / b }, 1),
            G({ 0, c / b, a / b, d / b }, 1), G({ c / b, 0, a / b, d / b }, 1),
            G({ c / b, 0, d / b, a / b }, 1), G({ c / b, d / b, 0, a / b }, 1),
            G({ 0, a }, { 1, sa }, x), G({ 0, a / b, c / b, d / b, e / b }, 1),
            G({ 0, c / b, a / b, d / b, e / b }, 1),
            G({ 0, c / b, d / b, a / b, e / b }, 1),
            G({ 0, c / b, d / b, e / b, a / b }, 1),
            G({ c / b, 0, a / b, d / b, e / b }, 1),
            G({ c / b, 0, d / b, a / b, e / b }, 1),
            G({ c / b, 0, d / b, e / b, a / b }, 1),
            G({ c / b, d / b, 0, a / b, e / b }, 1),
            G({ c / b, d / b, 0, e / b, a / b }, 1),
            G({ c / b, d / b, e / b, 0, a / b }, 1) };
        complex<double> res { -(sy[10] * sy[11]) + sy[0] * sy[1] - sy[2] * sy[3]
            + sy[1] * (sy[2] + sy[5]) + sy[3] * (-sy[5] + sy[6])
            + sy[4] * (-sy[6] - sy[7])
            + (-sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27]
                  - sy[28] - sy[29] + sy[19] * sy[7])
                * sy[8]
            + (-sy[12] + sy[13]) * G({ 0, c / b }, 1)
            + (sy[3] - sy[4]) * G({ 0, c / b, d / b }, 1)
            + sy[10]
                * (-sy[17] - sy[18] + G({ 0, c / b, d / b, e / b }, 1)
                    + G({ c / b, 0, d / b, e / b }, 1) + G({ c / b, d / b, 0, e / b }, 1))
            + 2. * G({ 0, 0, a / b, c / b, d / b, e / b }, 1)
            + 2. * G({ 0, 0, c / b, a / b, d / b, e / b }, 1)
            + 2. * G({ 0, 0, c / b, d / b, a / b, e / b }, 1)
            + 2. * G({ 0, 0, c / b, d / b, e / b, a / b }, 1)
            + G({ 0, a / b, 0, c / b, d / b, e / b }, 1)
            + G({ 0, a / b, c / b, 0, d / b, e / b }, 1)
            + G({ 0, a / b, c / b, d / b, 0, e / b }, 1)
            + G({ 0, a / b, c / b, d / b, e / b, x / b }, 1)
            + 2. * G({ 0, c / b, 0, a / b, d / b, e / b }, 1)
            + 2. * G({ 0, c / b, 0, d / b, a / b, e / b }, 1)
            + 2. * G({ 0, c / b, 0, d / b, e / b, a / b }, 1)
            + G({ 0, c / b, a / b, 0, d / b, e / b }, 1)
            + G({ 0, c / b, a / b, d / b, 0, e / b }, 1)
            + G({ 0, c / b, a / b, d / b, e / b, x / b }, 1)
            + 2. * G({ 0, c / b, d / b, 0, a / b, e / b }, 1)
            + 2. * G({ 0, c / b, d / b, 0, e / b, a / b }, 1)
            + G({ 0, c / b, d / b, a / b, 0, e / b }, 1)
            + G({ 0, c / b, d / b, a / b, e / b, x / b }, 1)
            + 2. * G({ 0, c / b, d / b, e / b, 0, a / b }, 1)
            + G({ 0, c / b, d / b, e / b, a / b, x / b }, 1)
            + 2. * G({ c / b, 0, 0, a / b, d / b, e / b }, 1)
            + 2. * G({ c / b, 0, 0, d / b, a / b, e / b }, 1)
            + 2. * G({ c / b, 0, 0, d / b, e / b, a / b }, 1)
            + G({ c / b, 0, a / b, 0, d / b, e / b }, 1)
            + G({ c / b, 0, a / b, d / b, 0, e / b }, 1)
            + G({ c / b, 0, a / b, d / b, e / b, x / b }, 1)
            + 2. * G({ c / b, 0, d / b, 0, a / b, e / b }, 1)
            + 2. * G({ c / b, 0, d / b, 0, e / b, a / b }, 1)
            + G({ c / b, 0, d / b, a / b, 0, e / b }, 1)
            + G({ c / b, 0, d / b, a / b, e / b, x / b }, 1)
            + 2. * G({ c / b, 0, d / b, e / b, 0, a / b }, 1)
            + G({ c / b, 0, d / b, e / b, a / b, x / b }, 1)
            + 2. * G({ c / b, d / b, 0, 0, a / b, e / b }, 1)
            + 2. * G({ c / b, d / b, 0, 0, e / b, a / b }, 1)
            + G({ c / b, d / b, 0, a / b, 0, e / b }, 1)
            + G({ c / b, d / b, 0, a / b, e / b, x / b }, 1)
            + 2. * G({ c / b, d / b, 0, e / b, 0, a / b }, 1)
            + G({ c / b, d / b, 0, e / b, a / b, x / b }, 1)
            + 2. * G({ c / b, d / b, e / b, 0, 0, a / b }, 1)
            + G({ c / b, d / b, e / b, 0, a / b, x / b }, 1)
            + (sy[11] + sy[14] + sy[15] + sy[16] + sy[17] + sy[18])
                * G({ 0, e }, { 1, se }, x)
            + sy[7] * G({ 0, 0, a }, { 1, 1, sa }, x)
            + G({ 0, a / b }, 1) * (-sy[13] + G({ 0, c, d, e }, { 1, sc, sd, se }, x))
            + G({ c / b, d / b }, 1) * (-sy[12] + G({ a, 0, 0, e }, { sa, 1, 1, se }, x))
            + G({ c / b }, 1)
                * (G({ a, 0, 0, d, e }, { sa, 1, 1, sd, se }, x)
                    - G({ a, 0, c, d, e }, { sa, 1, sc, sd, se }, x))
            + G({ a, 0, 0, c, d, e }, { sa, 1, 1, sc, sd, se }, x)
            + (sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27]
                  + sy[28] + sy[29] - sy[19] * sy[7] - sy[9] * sy[7] * sy[8])
                * Log(-x, sb)
            + (sy[9] * sy[7] * pow(sy[8], 2.)) / 2.
            + sy[9]
                * (-sy[23] - sy[26] - sy[28] - sy[29]
                    - G({ 0, 0, c / b, d / b, e / b }, 1)
                    - G({ 0, c / b, 0, d / b, e / b }, 1)
                    - G({ 0, c / b, d / b, 0, e / b }, 1)
                    - G({ c / b, 0, 0, d / b, e / b }, 1)
                    - G({ c / b, 0, d / b, 0, e / b }, 1)
                    - G({ c / b, d / b, 0, 0, e / b }, 1)
                    + sy[7] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
        if (c != x) {
            res += (-sy[0] + G({ 0, a / b, x / b }, 1))
                * G({ c, d, e }, { sc, sd, se }, x);
        }
        if (d != x) {
            res += (-sy[14] - sy[15] - sy[16] + G({ 0, a / b, c / b, x / b }, 1)
                       + G({ 0, c / b, a / b, x / b }, 1)
                       + G({ c / b, 0, a / b, x / b }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (-sy[20] - sy[21] - sy[22] - sy[24] - sy[25] - sy[27]
                       + G({ 0, a / b, c / b, d / b, x / b }, 1)
                       + G({ 0, c / b, a / b, d / b, x / b }, 1)
                       + G({ 0, c / b, d / b, a / b, x / b }, 1)
                       + G({ c / b, 0, a / b, d / b, x / b }, 1)
                       + G({ c / b, 0, d / b, a / b, x / b }, 1)
                       + G({ c / b, d / b, 0, a / b, x / b }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_a0bcde_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;

    // abcde
    const vector<complex<double>> sy
        = { G({ 0, b / a, c / a }, 1), G({ 0, b / a, c / a, d / a }, 1),
              G({ 0, b / a }, 1), G({ 0, b / a, c / a, d / a, e / a }, 1) };
    complex<double> res { 2. * G({ 0, 0, b / a, c / a, d / a, e / a }, 1)
        + G({ 0, b / a, 0, c / a, d / a, e / a }, 1)
        + G({ 0, b / a, c / a, 0, d / a, e / a }, 1)
        + G({ 0, b / a, c / a, d / a, 0, e / a }, 1)
        + G({ 0, b / a, c / a, d / a, e / a, x / a }, 1)
        + sy[1] * G({ 0, e }, { 1, se }, x) + sy[0] * G({ 0, d, e }, { 1, sd, se }, x)
        + sy[2] * G({ 0, c, d, e }, { 1, sc, sd, se }, x)
        + G({ x / a }, 1) * G({ 0, b, c, d, e }, { 1, sb, sc, sd, se }, x)
        + G({ 0, 0, b, c, d, e }, { 1, 1, sb, sc, sd, se }, x) - sy[3] * Log(a, sa)
        + sy[3] * Log(-x, sa) };
    if (b != x) {
        res += (-sy[2] + G({ 0, x / a }, 1)) * G({ b, c, d, e }, { sb, sc, sd, se }, x);
    }
    if (c != x) {
        res += (-sy[0] + G({ 0, b / a, x / a }, 1)) * G({ c, d, e }, { sc, sd, se }, x);
    }
    if (d != x) {
        res += (-sy[1] + G({ 0, b / a, c / a, x / a }, 1)) * G({ d, e }, { sd, se }, x);
    }
    if (e != x) {
        res += (-sy[3] + G({ 0, b / a, c / a, d / a, x / a }, 1)) * G({ e }, { se }, x);
    }
    return res;
}

complex<double> G6_ab0cde_e(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    const vector<complex<double>> sy
        = { Log(e, se), G({ a, b }, { sa, sb }, x), G({ d / e, c / e }, 1),
              G({ 0, d / e, c / e }, 1), G({ a, 0, b }, { sa, 1, sb }, x),
              G({ d / e, 0, c / e }, 1), G({ d / e, c / e, x / e }, 1),
              G({ 0, a, b }, { 1, sa, sb }, x), G({ a }, { sa }, x),
              G({ a, b, 0, c }, { sa, sb, 1, sc }, x), G({ d / e, c / e, 0, b / e }, 1),
              G({ a, b, 0, c, d }, { sa, sb, 1, sc, sd }, x), G({ d / e }, 1),
              G({ d / e, c / e, 0, b / e, a / e }, 1), G({ d / e, c / e, 0, x / e }, 1) };
    complex<double> res { -(sy[3] * sy[4]) + sy[4] * (-sy[5] - sy[6])
        + (-sy[3] - sy[5] - sy[6]) * sy[7]
        + sy[0]
            * (sy[9] * sy[12] + sy[13] + sy[1] * (-sy[3] - sy[5])
                + sy[2] * (sy[4] + sy[7]) - sy[10] * sy[8])
        - sy[9] * G({ 0, d / e }, 1) - G({ 0, d / e, c / e, 0, b / e, a / e }, 1)
        - G({ d / e, 0, c / e, 0, b / e, a / e }, 1)
        - 2. * G({ d / e, c / e, 0, 0, b / e, a / e }, 1)
        - G({ d / e, c / e, 0, b / e, 0, a / e }, 1)
        - G({ d / e, c / e, 0, b / e, a / e, x / e }, 1)
        + sy[11] * (-G({ x / e }, 1) + G({ e }, { se }, x))
        - sy[10] * G({ 0, a }, { 1, sa }, x)
        + sy[8]
            * (sy[13] + G({ 0, d / e, c / e, 0, b / e }, 1)
                + G({ d / e, 0, c / e, 0, b / e }, 1)
                + 2. * G({ d / e, c / e, 0, 0, b / e }, 1)
                + sy[6] * G({ 0, b }, { 1, sb }, x))
        - sy[6] * G({ 0, b, a }, { 1, sb, sa }, x)
        + sy[2]
            * (-sy[9] + sy[8] * G({ 0, 0, b }, { 1, 1, sb }, x)
                - G({ 0, 0, b, a }, { 1, 1, sb, sa }, x))
        + sy[12]
            * (sy[11] + G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x)
                + G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x)
                + 2. * G({ a, b, 0, 0, c }, { sa, sb, 1, 1, sc }, x))
        - G({ 0, a, b, 0, c, d }, { 1, sa, sb, 1, sc, sd }, x)
        - G({ a, 0, b, 0, c, d }, { sa, 1, sb, 1, sc, sd }, x)
        - 2. * G({ a, b, 0, 0, c, d }, { sa, sb, 1, 1, sc, sd }, x)
        - G({ a, b, 0, c, 0, d }, { sa, sb, 1, sc, 1, sd }, x)
        + (-(sy[9] * sy[12]) - sy[13] - sy[0] * sy[1] * sy[2] + sy[1] * (sy[3] + sy[5])
              + sy[2] * (-sy[4] - sy[7]) + sy[10] * sy[8])
            * Log(-x, se)
        + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.
        + sy[1]
            * (-sy[14] + G({ 0, 0, d / e, c / e }, 1) + G({ 0, d / e, 0, c / e }, 1)
                + G({ d / e, 0, 0, c / e }, 1)
                + sy[2] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
    if (b != x) {
        res += (-(sy[10] * sy[8]) + sy[14] * sy[8]) * G({ b }, { sb }, x)
            + (sy[10] - sy[14]) * G({ b, a }, { sb, sa }, x);
    }
    return res;
}
complex<double> G6_ab0cde_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (d == e) {
        const vector<complex<double>> sy = { G({ c / d, 0, b / d }, 1),
            G({ c / d, x / d, 1 }, 1), G({ a, b, 0, c }, { sa, sb, 1, sc }, x),
            G({ c / d, 0, b / d, 1 }, 1), G({ c / d, 0, b / d, a / d }, 1),
            G({ c / d, 0, 1, b / d }, 1), G({ c / d, 0, 1, x / d }, 1),
            G({ c / d, 0, x / d, 1 }, 1), G({ d }, { sd }, x),
            G({ c / d, 0, b / d, a / d, 1 }, 1), G({ a }, { sa }, x) };
        complex<double> res { -(sy[9] * sy[8])
            + sy[8] * G({ c / d, 0, b / d, a / d, x / d }, 1)
            + G({ c / d, 0, b / d, a / d, 0, 1 }, 1)
            - G({ c / d, 0, b / d, a / d, x / d, 1 }, 1)
            + sy[10]
                * (sy[9] - G({ c / d, 0, b / d, 0, 1 }, 1)
                    + sy[1] * G({ 0, b }, { 1, sb }, x))
            + sy[4] * G({ 0, d }, { 1, sd }, x)
            + (sy[5] - sy[6] - sy[7] + G({ c / d, 0, 0, 1 }, 1))
                * G({ a, b }, { sa, sb }, x)
            + (sy[3] - sy[4]) * G({ a, d }, { sa, sd }, x)
            + sy[1]
                * (-G({ 0, a, b }, { 1, sa, sb }, x) - G({ 0, b, a }, { 1, sb, sa }, x)
                    - G({ a, 0, b }, { sa, 1, sb }, x))
            - sy[0] * G({ a, 0, d }, { sa, 1, sd }, x)
            + (sy[0] - G({ c / d, 0, 1 }, 1)) * G({ a, b, d }, { sa, sb, sd }, x)
            + G({ c / d, 1 }, 1) * (-sy[2] + G({ a, b, 0, d }, { sa, sb, 1, sd }, x))
            + G({ c / d }, 1)
                * (-G({ a, b, 0, 0, d }, { sa, sb, 1, 1, sd }, x)
                    + G({ a, b, 0, c, d }, { sa, sb, 1, sc, sd }, x))
            + G({ a, b, 0, c, 0, d }, { sa, sb, 1, sc, 1, sd }, x) - sy[2] * Zeta(2) };
        if (b != x) {
            res += (-(sy[10] * sy[3]) - sy[10] * sy[5] + sy[10] * sy[6] + sy[10] * sy[7])
                    * G({ b }, { sb }, x)
                + (sy[3] + sy[5] - sy[6] - sy[7]) * G({ b, a }, { sb, sa }, x);
        }
        return res;
    } else {
        const vector<complex<double>> sy = { Log(d, sd), G({ a, b }, { sa, sb }, x),
            G({ c / d, e / d }, 1), G({ e / d, c / d }, 1), G({ 0, c / d, e / d }, 1),
            G({ a, 0, b }, { sa, 1, sb }, x), G({ 0, e / d, c / d }, 1),
            G({ c / d, 0, b / d }, 1), G({ a, b, e }, { sa, sb, se }, x),
            G({ 0, b, a }, { 1, sb, sa }, x), G({ c / d, e / d, x / d }, 1),
            G({ e / d, c / d, x / d }, 1), G({ e / d, 0, c / d }, 1),
            G({ 0, a, b }, { 1, sa, sb }, x), G({ a, b, 0, c }, { sa, sb, 1, sc }, x),
            G({ a }, { sa }, x), G({ 0, 0, b }, { 1, 1, sb }, x),
            G({ 0, 0, b, a }, { 1, 1, sb, sa }, x), G({ c / d, 0, b / d, a / d }, 1),
            G({ a, e }, { sa, se }, x), G({ 0, a }, { 1, sa }, x),
            G({ c / d, 0, b / d, e / d }, 1), G({ c / d, 0, e / d, b / d }, 1),
            G({ c / d, e / d, 0, b / d }, 1), G({ e / d, c / d, 0, b / d }, 1),
            G({ e / d }, 1), G({ a, b, 0, c, e }, { sa, sb, 1, sc, se }, x),
            G({ c / d, 0, b / d, a / d, e / d }, 1),
            G({ c / d, 0, b / d, e / d, a / d }, 1),
            G({ c / d, 0, e / d, b / d, a / d }, 1),
            G({ c / d, e / d, 0, b / d, a / d }, 1),
            G({ e / d, c / d, 0, b / d, a / d }, 1), G({ 0, 0 }, { 1, 1 }, x),
            G({ c / d, e / d, 0, x / d }, 1), G({ e / d, c / d, 0, x / d }, 1) };
        complex<double> res { sy[9] * sy[10] + sy[9] * sy[11] - sy[18] * sy[19]
            + sy[19] * sy[21] + sy[20] * sy[21] + sy[20] * (sy[22] + sy[23] + sy[24])
            + (sy[14] - sy[15] * sy[16] + sy[17]) * sy[3]
            + (sy[10] + sy[11] + sy[12]) * sy[5] + sy[4] * sy[5] + sy[5] * sy[6]
            + sy[13] * (sy[10] + sy[11] + sy[12] + sy[4] + sy[6])
            + sy[0]
                * (sy[15] * (sy[21] + sy[22] + sy[23] + sy[24]) - sy[14] * sy[25] - sy[27]
                    - sy[28] - sy[29] - sy[30] - sy[31] + sy[2] * (-sy[13] - sy[5])
                    + sy[3] * (-sy[13] - sy[5]) + sy[1] * (sy[12] + sy[4] + sy[6]))
            + sy[7] * sy[8] + sy[14] * G({ 0, e / d }, 1)
            - sy[8] * G({ c / d, 0, e / d }, 1)
            + G({ 0, c / d, 0, b / d, a / d, e / d }, 1)
            + G({ 0, c / d, 0, b / d, e / d, a / d }, 1)
            + G({ 0, c / d, 0, e / d, b / d, a / d }, 1)
            + G({ 0, c / d, e / d, 0, b / d, a / d }, 1)
            + G({ 0, e / d, c / d, 0, b / d, a / d }, 1)
            + 2. * G({ c / d, 0, 0, b / d, a / d, e / d }, 1)
            + 2. * G({ c / d, 0, 0, b / d, e / d, a / d }, 1)
            + 2. * G({ c / d, 0, 0, e / d, b / d, a / d }, 1)
            + G({ c / d, 0, b / d, 0, a / d, e / d }, 1)
            + G({ c / d, 0, b / d, 0, e / d, a / d }, 1)
            + G({ c / d, 0, b / d, a / d, 0, e / d }, 1)
            + G({ c / d, 0, b / d, a / d, e / d, x / d }, 1)
            + G({ c / d, 0, b / d, e / d, 0, a / d }, 1)
            + G({ c / d, 0, b / d, e / d, a / d, x / d }, 1)
            + 2. * G({ c / d, 0, e / d, 0, b / d, a / d }, 1)
            + G({ c / d, 0, e / d, b / d, 0, a / d }, 1)
            + G({ c / d, 0, e / d, b / d, a / d, x / d }, 1)
            + 2. * G({ c / d, e / d, 0, 0, b / d, a / d }, 1)
            + G({ c / d, e / d, 0, b / d, 0, a / d }, 1)
            + G({ c / d, e / d, 0, b / d, a / d, x / d }, 1)
            + G({ e / d, 0, c / d, 0, b / d, a / d }, 1)
            + 2. * G({ e / d, c / d, 0, 0, b / d, a / d }, 1)
            + G({ e / d, c / d, 0, b / d, 0, a / d }, 1)
            + G({ e / d, c / d, 0, b / d, a / d, x / d }, 1)
            + sy[15]
                * (-sy[28] - sy[29] - sy[30] - sy[31]
                    - G({ 0, c / d, 0, b / d, e / d }, 1)
                    - G({ 0, c / d, 0, e / d, b / d }, 1)
                    - G({ 0, c / d, e / d, 0, b / d }, 1)
                    - G({ 0, e / d, c / d, 0, b / d }, 1)
                    - 2. * G({ c / d, 0, 0, b / d, e / d }, 1)
                    - 2. * G({ c / d, 0, 0, e / d, b / d }, 1)
                    - G({ c / d, 0, b / d, 0, e / d }, 1)
                    - 2. * G({ c / d, 0, e / d, 0, b / d }, 1)
                    - 2. * G({ c / d, e / d, 0, 0, b / d }, 1)
                    - G({ e / d, 0, c / d, 0, b / d }, 1)
                    - 2. * G({ e / d, c / d, 0, 0, b / d }, 1)
                    + (-sy[10] - sy[11]) * G({ 0, b }, { 1, sb }, x))
            + sy[18] * G({ 0, e }, { 1, se }, x)
            - sy[7] * G({ a, 0, e }, { sa, 1, se }, x)
            + sy[2]
                * (-(sy[15] * sy[16]) + sy[17] + G({ a, b, 0, e }, { sa, sb, 1, se }, x))
            + sy[25]
                * (-sy[26] - G({ 0, a, b, 0, c }, { 1, sa, sb, 1, sc }, x)
                    - G({ a, 0, b, 0, c }, { sa, 1, sb, 1, sc }, x)
                    - 2. * G({ a, b, 0, 0, c }, { sa, sb, 1, 1, sc }, x))
            + G({ c / d }, 1) * (sy[26] - G({ a, b, 0, 0, e }, { sa, sb, 1, 1, se }, x))
            + G({ a, b, 0, c, 0, e }, { sa, sb, 1, sc, 1, se }, x)
            + (sy[15] * (-sy[21] - sy[22] - sy[23] - sy[24]) + sy[14] * sy[25] + sy[27]
                  + sy[28] + sy[29] + sy[30] + sy[31] + sy[0] * sy[1] * (sy[2] + sy[3])
                  + sy[2] * (sy[13] + sy[5]) + sy[3] * (sy[13] + sy[5])
                  + sy[1] * (-sy[12] - sy[4] - sy[6]))
                * Log(-x, sd)
            + sy[1] * (-sy[2] / 2. - sy[3] / 2.) * pow(sy[0], 2.)
            + sy[1]
                * (sy[22] + sy[33] + sy[34] + sy[32] * sy[3]
                    - G({ 0, 0, c / d, e / d }, 1) - G({ 0, 0, e / d, c / d }, 1)
                    - G({ 0, e / d, 0, c / d }, 1) + G({ c / d, 0, 0, e / d }, 1)
                    - G({ e / d, 0, 0, c / d }, 1) + 2. * sy[3] * Zeta(2)
                    + sy[2] * (sy[32] + 2. * Zeta(2))) };
        if (b != x) {
            res += (sy[15] * sy[23] + sy[15] * sy[24] - sy[15] * sy[33] - sy[15] * sy[34])
                    * G({ b }, { sb }, x)
                + (-sy[23] - sy[24] + sy[33] + sy[34]) * G({ b, a }, { sb, sa }, x);
        }
        if (e != x) {
            res += (-sy[27] + G({ c / d, 0, b / d, a / d, x / d }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_ab0cde_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (c == d) {
        if (c == e) {
            const vector<complex<double>> sy = { G({ a, b, c }, { sa, sb, sc }, x),
                G({ a, c, c }, { sa, sc, sc }, x), G({ x / c, 1, 1 }, 1),
                G({ 0, b / c, 1, 1 }, 1), G({ 0, b / c, a / c, 1 }, 1),
                G({ c, c }, { sc, sc }, x), G({ a, b, c, c }, { sa, sb, sc, sc }, x),
                G({ c }, { sc }, x), G({ 0, b / c, a / c, 1, 1 }, 1), G({ a }, { sa }, x),
                G({ 0, 1, 1, b / c }, 1), G({ 0, 1, 1, x / c }, 1),
                G({ 0, 1, b / c, 1 }, 1), G({ 0, 1, x / c, 1 }, 1),
                G({ 0, x / c, 1, 1 }, 1) };
            complex<double> res { sy[4] * sy[5] - sy[7] * sy[8]
                + (sy[0] - sy[1]) * G({ 0, b / c, 1 }, 1)
                - sy[5] * G({ 0, b / c, a / c, x / c }, 1)
                + sy[7] * G({ 0, b / c, a / c, x / c, 1 }, 1)
                + G({ 0, b / c, a / c, 0, 1, 1 }, 1)
                - G({ 0, b / c, a / c, x / c, 1, 1 }, 1)
                + sy[9]
                    * (sy[8] - G({ 0, b / c, 0, 1, 1 }, 1)
                        + sy[2] * G({ 0, b }, { 1, sb }, x))
                + (sy[3] - sy[4]) * G({ a, c }, { sa, sc }, x)
                + G({ 0, b / c, a / c }, 1) * (sy[1] - G({ 0, c, c }, { 1, sc, sc }, x))
                + sy[2]
                    * (-G({ 0, a, b }, { 1, sa, sb }, x)
                        - G({ 0, b, a }, { 1, sb, sa }, x)
                        - G({ a, 0, b }, { sa, 1, sb }, x))
                + G({ 0, b / c }, 1) * (-sy[6] + G({ a, 0, c, c }, { sa, 1, sc, sc }, x))
                + G({ a, b, 0, 0, c, c }, { sa, sb, 1, 1, sc, sc }, x) - sy[6] * Zeta(2)
                + G({ a, b, 0, c }, { sa, sb, 1, sc }, x) * Zeta(2)
                + G({ a, b }, { sa, sb }, x)
                    * (sy[10] - sy[11] + sy[12] - sy[13] - sy[14] + Zeta(4) / 4.)
                - sy[0] * Zeta(3) };
            if (b != x) {
                res += (-(sy[9] * sy[10]) + sy[9] * sy[11] - sy[9] * sy[12]
                           + sy[9] * sy[13] + sy[9] * sy[14] - sy[9] * sy[3])
                        * G({ b }, { sb }, x)
                    + (sy[10] - sy[11] + sy[12] - sy[13] - sy[14] + sy[3])
                        * G({ b, a }, { sb, sa }, x);
            }
            return res;
        } else {
            const vector<complex<double>> sy = { Log(c, sc), G({ a, b }, { sa, sb }, x),
                G({ e / c, 1 }, 1), G({ a, b, e }, { sa, sb, se }, x),
                G({ a, c, e }, { sa, sc, se }, x), G({ 0, b, a }, { 1, sb, sa }, x),
                G({ e / c, 1, x / c }, 1), G({ e / c, x / c, 1 }, 1),
                G({ 0, a, b }, { 1, sa, sb }, x), G({ a, 0, b }, { sa, 1, sb }, x),
                G({ c, e }, { sc, se }, x), G({ 0, b / c, a / c, 1 }, 1),
                G({ 0, a }, { 1, sa }, x), G({ 0, b / c, e / c, 1 }, 1),
                G({ a }, { sa }, x), G({ a, b, 0, e }, { sa, sb, 1, se }, x),
                G({ a, b, c, e }, { sa, sb, sc, se }, x), G({ 0, e / c, 1, b / c }, 1),
                G({ 0, e / c, b / c, 1 }, 1), G({ e / c, 0, 1, b / c }, 1),
                G({ e / c, 0, b / c, 1 }, 1), G({ e / c, 1, 0, b / c }, 1),
                G({ 0, b / c, a / c, e / c, 1 }, 1), G({ 0, b / c, e / c, 1, a / c }, 1),
                G({ 0, b / c, e / c, a / c, 1 }, 1), G({ 0, e / c, 1, b / c, a / c }, 1),
                G({ 0, e / c, b / c, 1, a / c }, 1), G({ 0, e / c, b / c, a / c, 1 }, 1),
                G({ e / c, 0, 1, b / c, a / c }, 1), G({ e / c, 0, b / c, 1, a / c }, 1),
                G({ e / c, 0, b / c, a / c, 1 }, 1), G({ e / c, 1, 0, b / c, a / c }, 1),
                G({ e / c, 0, 1, x / c }, 1), G({ e / c, 0, x / c, 1 }, 1),
                G({ e / c, 1, 0, x / c }, 1) };
            complex<double> res { sy[10] * sy[11] + sy[12] * sy[13]
                + sy[12] * (sy[17] + sy[18] + sy[19] + sy[20] + sy[21]) + sy[5] * sy[6]
                + sy[5] * sy[7] + sy[9] * (sy[6] + sy[7])
                + sy[0]
                    * (sy[14] * (sy[13] + sy[17] + sy[18] + sy[19] + sy[20] + sy[21])
                        - sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27] - sy[28]
                        - sy[29] - sy[30] - sy[31] + sy[2] * (-sy[9] - sy[8]))
                + (sy[6] + sy[7]) * sy[8] + (sy[3] - sy[4]) * G({ 0, b / c, 1 }, 1)
                - sy[3] * G({ 0, e / c, 1 }, 1)
                - sy[10] * G({ 0, b / c, a / c, x / c }, 1)
                + 2. * G({ 0, 0, b / c, a / c, e / c, 1 }, 1)
                + 2. * G({ 0, 0, b / c, e / c, 1, a / c }, 1)
                + 2. * G({ 0, 0, b / c, e / c, a / c, 1 }, 1)
                + 2. * G({ 0, 0, e / c, 1, b / c, a / c }, 1)
                + 2. * G({ 0, 0, e / c, b / c, 1, a / c }, 1)
                + 2. * G({ 0, 0, e / c, b / c, a / c, 1 }, 1)
                + G({ 0, b / c, 0, a / c, e / c, 1 }, 1)
                + G({ 0, b / c, 0, e / c, 1, a / c }, 1)
                + G({ 0, b / c, 0, e / c, a / c, 1 }, 1)
                + G({ 0, b / c, a / c, 0, e / c, 1 }, 1)
                + G({ 0, b / c, a / c, e / c, 1, x / c }, 1)
                + G({ 0, b / c, a / c, e / c, x / c, 1 }, 1)
                + G({ 0, b / c, e / c, 0, 1, a / c }, 1)
                + G({ 0, b / c, e / c, 0, a / c, 1 }, 1)
                + G({ 0, b / c, e / c, 1, 0, a / c }, 1)
                + G({ 0, b / c, e / c, 1, a / c, x / c }, 1)
                + G({ 0, b / c, e / c, a / c, 1, x / c }, 1)
                + G({ 0, b / c, e / c, a / c, x / c, 1 }, 1)
                + 2. * G({ 0, e / c, 0, 1, b / c, a / c }, 1)
                + 2. * G({ 0, e / c, 0, b / c, 1, a / c }, 1)
                + 2. * G({ 0, e / c, 0, b / c, a / c, 1 }, 1)
                + 2. * G({ 0, e / c, 1, 0, b / c, a / c }, 1)
                + G({ 0, e / c, 1, b / c, 0, a / c }, 1)
                + G({ 0, e / c, 1, b / c, a / c, x / c }, 1)
                + G({ 0, e / c, b / c, 0, 1, a / c }, 1)
                + G({ 0, e / c, b / c, 0, a / c, 1 }, 1)
                + G({ 0, e / c, b / c, 1, 0, a / c }, 1)
                + G({ 0, e / c, b / c, 1, a / c, x / c }, 1)
                + G({ 0, e / c, b / c, a / c, 1, x / c }, 1)
                + G({ 0, e / c, b / c, a / c, x / c, 1 }, 1)
                + 2. * G({ e / c, 0, 0, 1, b / c, a / c }, 1)
                + 2. * G({ e / c, 0, 0, b / c, 1, a / c }, 1)
                + 2. * G({ e / c, 0, 0, b / c, a / c, 1 }, 1)
                + 2. * G({ e / c, 0, 1, 0, b / c, a / c }, 1)
                + G({ e / c, 0, 1, b / c, 0, a / c }, 1)
                + G({ e / c, 0, 1, b / c, a / c, x / c }, 1)
                + G({ e / c, 0, b / c, 0, 1, a / c }, 1)
                + G({ e / c, 0, b / c, 0, a / c, 1 }, 1)
                + G({ e / c, 0, b / c, 1, 0, a / c }, 1)
                + G({ e / c, 0, b / c, 1, a / c, x / c }, 1)
                + G({ e / c, 0, b / c, a / c, 1, x / c }, 1)
                + G({ e / c, 0, b / c, a / c, x / c, 1 }, 1)
                + 2. * G({ e / c, 1, 0, 0, b / c, a / c }, 1)
                + G({ e / c, 1, 0, b / c, 0, a / c }, 1)
                + G({ e / c, 1, 0, b / c, a / c, x / c }, 1)
                + sy[14]
                    * (-sy[23] - sy[24] - sy[25] - sy[26] - sy[27] - sy[28] - sy[29]
                        - sy[30] - sy[31] - 2. * G({ 0, 0, b / c, e / c, 1 }, 1)
                        - 2. * G({ 0, 0, e / c, 1, b / c }, 1)
                        - 2. * G({ 0, 0, e / c, b / c, 1 }, 1)
                        - G({ 0, b / c, 0, e / c, 1 }, 1)
                        - 2. * G({ 0, e / c, 0, 1, b / c }, 1)
                        - 2. * G({ 0, e / c, 0, b / c, 1 }, 1)
                        - 2. * G({ 0, e / c, 1, 0, b / c }, 1)
                        - 2. * G({ e / c, 0, 0, 1, b / c }, 1)
                        - 2. * G({ e / c, 0, 0, b / c, 1 }, 1)
                        - 2. * G({ e / c, 0, 1, 0, b / c }, 1)
                        - 2. * G({ e / c, 1, 0, 0, b / c }, 1)
                        + (-sy[6] - sy[7]) * G({ 0, b }, { 1, sb }, x))
                + (-sy[11] + sy[13]) * G({ a, e }, { sa, se }, x)
                + G({ 0, b / c, a / c }, 1) * (sy[4] - G({ 0, c, e }, { 1, sc, se }, x))
                + sy[2]
                    * (sy[15] - sy[14] * G({ 0, 0, b }, { 1, 1, sb }, x)
                        + G({ 0, 0, b, a }, { 1, 1, sb, sa }, x))
                + G({ 0, b / c }, 1) * (-sy[16] + G({ a, 0, c, e }, { sa, 1, sc, se }, x))
                + G({ a, b, 0, 0, c, e }, { sa, sb, 1, 1, sc, se }, x)
                + (sy[14] * (-sy[13] - sy[17] - sy[18] - sy[19] - sy[20] - sy[21])
                      + sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27] + sy[28]
                      + sy[0] * sy[1] * sy[2] + sy[29] + sy[30] + sy[31]
                      + sy[2] * (sy[9] + sy[8]))
                    * Log(-x, sc)
                - (sy[1] * sy[2] * pow(sy[0], 2.)) / 2. + sy[15] * Zeta(2)
                - sy[16] * Zeta(2)
                + sy[1]
                    * (sy[17] + sy[18] + sy[32] + sy[33] + sy[34]
                        + G({ 0, 0, e / c, 1 }, 1)
                        + sy[2] * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2))) };
            if (b != x) {
                res += (sy[14] * sy[19] + sy[14] * sy[20] + sy[14] * sy[21]
                           - sy[14] * sy[32] - sy[14] * sy[33] - sy[14] * sy[34])
                        * G({ b }, { sb }, x)
                    + (-sy[19] - sy[20] - sy[21] + sy[32] + sy[33] + sy[34])
                        * G({ b, a }, { sb, sa }, x);
            }
            if (e != x) {
                res += (-sy[22] + G({ 0, b / c, a / c, x / c, 1 }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else {
        const vector<complex<double>> sy = { Log(c, sc), G({ a, b }, { sa, sb }, x),
            G({ d / c, e / c }, 1), G({ 0, d / c, b / c }, 1),
            G({ a, b, e }, { sa, sb, se }, x), G({ a, 0, e }, { sa, 1, se }, x),
            G({ a, d, e }, { sa, sd, se }, x), G({ d / c, 0, b / c }, 1),
            G({ 0, a, b }, { 1, sa, sb }, x), G({ d / c, e / c, x / c }, 1),
            G({ a, 0, b }, { sa, 1, sb }, x), G({ a, e }, { sa, se }, x),
            G({ 0, b / c, d / c, a / c }, 1), G({ 0, a }, { 1, sa }, x),
            G({ 0, b / c, d / c, e / c }, 1), G({ a }, { sa }, x),
            G({ a, b, 0, e }, { sa, sb, 1, se }, x),
            G({ a, b, d, e }, { sa, sb, sd, se }, x), G({ 0, b / c, a / c, d / c }, 1),
            G({ 0, d / c, b / c, a / c }, 1), G({ d / c, 0, b / c, a / c }, 1),
            G({ 0, d / c, b / c, e / c }, 1), G({ d / c, 0, b / c, e / c }, 1),
            G({ 0, d / c, e / c, b / c }, 1), G({ d / c, 0, e / c, b / c }, 1),
            G({ d / c, e / c, 0, b / c }, 1), G({ 0, b / c, a / c, d / c, e / c }, 1),
            G({ 0, b / c, d / c, a / c, e / c }, 1),
            G({ 0, b / c, d / c, e / c, a / c }, 1),
            G({ 0, d / c, b / c, a / c, e / c }, 1),
            G({ 0, d / c, b / c, e / c, a / c }, 1),
            G({ 0, d / c, e / c, b / c, a / c }, 1),
            G({ d / c, 0, b / c, a / c, e / c }, 1),
            G({ d / c, 0, b / c, e / c, a / c }, 1),
            G({ d / c, 0, e / c, b / c, a / c }, 1),
            G({ d / c, e / c, 0, b / c, a / c }, 1), G({ d / c, e / c, 0, x / c }, 1) };
        complex<double> res { -(sy[9] * sy[10]) + sy[11] * sy[12] - sy[13] * sy[14]
            + sy[11] * (-sy[14] + sy[19] + sy[20] - sy[21] - sy[22])
            + sy[13] * (-sy[21] - sy[22] - sy[23] - sy[24] - sy[25]) - sy[3] * sy[4]
            + sy[5] * (sy[3] + sy[7]) - sy[9] * sy[8]
            + sy[0]
                * (sy[15] * (-sy[14] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25])
                    + sy[26] + sy[27] + sy[28] + sy[29] + sy[30] + sy[31] + sy[32]
                    + sy[33] + sy[34] + sy[35] + sy[2] * (sy[10] + sy[8]))
            + (-sy[16] + sy[17]) * G({ 0, d / c }, 1)
            + (sy[5] - sy[6]) * G({ 0, b / c, d / c }, 1)
            + sy[4] * (-sy[7] + G({ 0, d / c, e / c }, 1) + G({ d / c, 0, e / c }, 1))
            - 2. * G({ 0, 0, b / c, a / c, d / c, e / c }, 1)
            - 2. * G({ 0, 0, b / c, d / c, a / c, e / c }, 1)
            - 2. * G({ 0, 0, b / c, d / c, e / c, a / c }, 1)
            - 2. * G({ 0, 0, d / c, b / c, a / c, e / c }, 1)
            - 2. * G({ 0, 0, d / c, b / c, e / c, a / c }, 1)
            - 2. * G({ 0, 0, d / c, e / c, b / c, a / c }, 1)
            - G({ 0, b / c, 0, a / c, d / c, e / c }, 1)
            - G({ 0, b / c, 0, d / c, a / c, e / c }, 1)
            - G({ 0, b / c, 0, d / c, e / c, a / c }, 1)
            - G({ 0, b / c, a / c, 0, d / c, e / c }, 1)
            - G({ 0, b / c, a / c, d / c, 0, e / c }, 1)
            - G({ 0, b / c, a / c, d / c, e / c, x / c }, 1)
            - G({ 0, b / c, d / c, 0, a / c, e / c }, 1)
            - G({ 0, b / c, d / c, 0, e / c, a / c }, 1)
            - G({ 0, b / c, d / c, a / c, 0, e / c }, 1)
            - G({ 0, b / c, d / c, a / c, e / c, x / c }, 1)
            - G({ 0, b / c, d / c, e / c, 0, a / c }, 1)
            - G({ 0, b / c, d / c, e / c, a / c, x / c }, 1)
            - 2. * G({ 0, d / c, 0, b / c, a / c, e / c }, 1)
            - 2. * G({ 0, d / c, 0, b / c, e / c, a / c }, 1)
            - 2. * G({ 0, d / c, 0, e / c, b / c, a / c }, 1)
            - G({ 0, d / c, b / c, 0, a / c, e / c }, 1)
            - G({ 0, d / c, b / c, 0, e / c, a / c }, 1)
            - G({ 0, d / c, b / c, a / c, 0, e / c }, 1)
            - G({ 0, d / c, b / c, a / c, e / c, x / c }, 1)
            - G({ 0, d / c, b / c, e / c, 0, a / c }, 1)
            - G({ 0, d / c, b / c, e / c, a / c, x / c }, 1)
            - 2. * G({ 0, d / c, e / c, 0, b / c, a / c }, 1)
            - G({ 0, d / c, e / c, b / c, 0, a / c }, 1)
            - G({ 0, d / c, e / c, b / c, a / c, x / c }, 1)
            - 2. * G({ d / c, 0, 0, b / c, a / c, e / c }, 1)
            - 2. * G({ d / c, 0, 0, b / c, e / c, a / c }, 1)
            - 2. * G({ d / c, 0, 0, e / c, b / c, a / c }, 1)
            - G({ d / c, 0, b / c, 0, a / c, e / c }, 1)
            - G({ d / c, 0, b / c, 0, e / c, a / c }, 1)
            - G({ d / c, 0, b / c, a / c, 0, e / c }, 1)
            - G({ d / c, 0, b / c, a / c, e / c, x / c }, 1)
            - G({ d / c, 0, b / c, e / c, 0, a / c }, 1)
            - G({ d / c, 0, b / c, e / c, a / c, x / c }, 1)
            - 2. * G({ d / c, 0, e / c, 0, b / c, a / c }, 1)
            - G({ d / c, 0, e / c, b / c, 0, a / c }, 1)
            - G({ d / c, 0, e / c, b / c, a / c, x / c }, 1)
            - 2. * G({ d / c, e / c, 0, 0, b / c, a / c }, 1)
            - G({ d / c, e / c, 0, b / c, 0, a / c }, 1)
            - G({ d / c, e / c, 0, b / c, a / c, x / c }, 1)
            + sy[15]
                * (sy[28] + sy[30] + sy[31] + sy[33] + sy[34] + sy[35]
                    + 2. * G({ 0, 0, b / c, d / c, e / c }, 1)
                    + 2. * G({ 0, 0, d / c, b / c, e / c }, 1)
                    + 2. * G({ 0, 0, d / c, e / c, b / c }, 1)
                    + G({ 0, b / c, 0, d / c, e / c }, 1)
                    + G({ 0, b / c, d / c, 0, e / c }, 1)
                    + 2. * G({ 0, d / c, 0, b / c, e / c }, 1)
                    + 2. * G({ 0, d / c, 0, e / c, b / c }, 1)
                    + G({ 0, d / c, b / c, 0, e / c }, 1)
                    + 2. * G({ 0, d / c, e / c, 0, b / c }, 1)
                    + 2. * G({ d / c, 0, 0, b / c, e / c }, 1)
                    + 2. * G({ d / c, 0, 0, e / c, b / c }, 1)
                    + G({ d / c, 0, b / c, 0, e / c }, 1)
                    + 2. * G({ d / c, 0, e / c, 0, b / c }, 1)
                    + 2. * G({ d / c, e / c, 0, 0, b / c }, 1)
                    + sy[9] * G({ 0, b }, { 1, sb }, x))
            + (-sy[12] - sy[18] - sy[19] - sy[20]) * G({ 0, e }, { 1, se }, x)
            - sy[9] * G({ 0, b, a }, { 1, sb, sa }, x)
            + G({ 0, b / c, a / c }, 1) * (sy[6] - G({ 0, d, e }, { 1, sd, se }, x))
            + sy[2]
                * (-sy[16] + sy[15] * G({ 0, 0, b }, { 1, 1, sb }, x)
                    - G({ 0, 0, b, a }, { 1, 1, sb, sa }, x))
            + G({ 0, b / c }, 1) * (-sy[17] + G({ a, 0, d, e }, { sa, 1, sd, se }, x))
            + G({ d / c }, 1)
                * (G({ a, b, 0, 0, e }, { sa, sb, 1, 1, se }, x)
                    - G({ a, b, 0, d, e }, { sa, sb, 1, sd, se }, x))
            + G({ a, b, 0, 0, d, e }, { sa, sb, 1, 1, sd, se }, x)
            + (sy[15] * (sy[14] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25]) - sy[26]
                  - sy[27] - sy[28] - sy[0] * sy[1] * sy[2] - sy[29] - sy[30] - sy[31]
                  - sy[32] - sy[33] - sy[34] - sy[35] + sy[2] * (-sy[10] - sy[8]))
                * Log(-x, sc)
            + (sy[1] * sy[2] * pow(sy[0], 2.)) / 2.
            + sy[1]
                * (-sy[23] - sy[24] - sy[36] - G({ 0, 0, d / c, e / c }, 1)
                    - G({ 0, d / c, 0, e / c }, 1) - G({ d / c, 0, 0, e / c }, 1)
                    + sy[2] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
        if (b != x) {
            res += (-(sy[15] * sy[25]) + sy[15] * sy[36]) * G({ b }, { sb }, x)
                + (sy[25] - sy[36]) * G({ b, a }, { sb, sa }, x);
        }
        if (d != x) {
            res += (sy[18] - G({ 0, b / c, a / c, x / c }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (sy[26] + sy[27] + sy[29] + sy[32]
                       - G({ 0, b / c, a / c, d / c, x / c }, 1)
                       - G({ 0, b / c, d / c, a / c, x / c }, 1)
                       - G({ 0, d / c, b / c, a / c, x / c }, 1)
                       - G({ d / c, 0, b / c, a / c, x / c }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_ab0cde_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;

    // abcde
    const vector<complex<double>> sy = { G({ 0, a / b, c / b }, 1),
        G({ 0, d, e }, { 1, sd, se }, x), G({ a, d, e }, { sa, sd, se }, x),
        G({ 0, c / b, a / b }, 1), G({ a / b, 0, c / b }, 1),
        G({ 0, c, d, e }, { 1, sc, sd, se }, x), G({ a, e }, { sa, se }, x),
        G({ 0, c / b, d / b, a / b }, 1), G({ 0, c / b, d / b, e / b }, 1),
        G({ a, c, d, e }, { sa, sc, sd, se }, x), G({ 0, a / b, c / b, d / b }, 1),
        G({ 0, c / b, a / b, d / b }, 1), G({ a / b, 0, c / b, d / b }, 1),
        G({ a }, { sa }, x), G({ 0, c / b, d / b, e / b, a / b }, 1),
        G({ 0, a / b, c / b, d / b, e / b }, 1), G({ 0, c / b, a / b, d / b, e / b }, 1),
        G({ 0, c / b, d / b, a / b, e / b }, 1),
        G({ a / b, 0, c / b, d / b, e / b }, 1) };
    complex<double> res { -(sy[0] * sy[1]) + sy[2] * sy[3] + sy[1] * (-sy[3] - sy[4])
        + sy[6] * sy[7] - sy[6] * sy[8] + (sy[9] - sy[5]) * G({ 0, a / b }, 1)
        - sy[5] * G({ a / b, x / b }, 1)
        + sy[13]
            * (sy[14] + 2. * G({ 0, 0, c / b, d / b, e / b }, 1)
                + G({ 0, c / b, 0, d / b, e / b }, 1)
                + G({ 0, c / b, d / b, 0, e / b }, 1))
        - 2. * G({ 0, 0, a / b, c / b, d / b, e / b }, 1)
        - 2. * G({ 0, 0, c / b, a / b, d / b, e / b }, 1)
        - 2. * G({ 0, 0, c / b, d / b, a / b, e / b }, 1)
        - 2. * G({ 0, 0, c / b, d / b, e / b, a / b }, 1)
        - 2. * G({ 0, a / b, 0, c / b, d / b, e / b }, 1)
        - G({ 0, a / b, c / b, 0, d / b, e / b }, 1)
        - G({ 0, a / b, c / b, d / b, 0, e / b }, 1)
        - G({ 0, a / b, c / b, d / b, e / b, x / b }, 1)
        - G({ 0, c / b, 0, a / b, d / b, e / b }, 1)
        - G({ 0, c / b, 0, d / b, a / b, e / b }, 1)
        - G({ 0, c / b, 0, d / b, e / b, a / b }, 1)
        - G({ 0, c / b, a / b, 0, d / b, e / b }, 1)
        - G({ 0, c / b, a / b, d / b, 0, e / b }, 1)
        - G({ 0, c / b, a / b, d / b, e / b, x / b }, 1)
        - G({ 0, c / b, d / b, 0, a / b, e / b }, 1)
        - G({ 0, c / b, d / b, 0, e / b, a / b }, 1)
        - G({ 0, c / b, d / b, a / b, 0, e / b }, 1)
        - G({ 0, c / b, d / b, a / b, e / b, x / b }, 1)
        - G({ 0, c / b, d / b, e / b, 0, a / b }, 1)
        - G({ 0, c / b, d / b, e / b, a / b, x / b }, 1)
        - 2. * G({ a / b, 0, 0, c / b, d / b, e / b }, 1)
        - G({ a / b, 0, c / b, 0, d / b, e / b }, 1)
        - G({ a / b, 0, c / b, d / b, 0, e / b }, 1)
        - G({ a / b, 0, c / b, d / b, e / b, x / b }, 1)
        - sy[8] * G({ 0, a }, { 1, sa }, x)
        + (-sy[10] - sy[11] - sy[12] - sy[7]) * G({ 0, e }, { 1, se }, x)
        + G({ 0, c / b, d / b }, 1) * (-sy[2] + G({ a, 0, e }, { sa, 1, se }, x))
        + G({ 0, c / b }, 1) * (-sy[9] + G({ a, 0, d, e }, { sa, 1, sd, se }, x))
        + G({ a / b }, 1)
            * (-G({ 0, 0, c, d, e }, { 1, 1, sc, sd, se }, x)
                + G({ a, 0, c, d, e }, { sa, 1, sc, sd, se }, x))
        + G({ a, 0, 0, c, d, e }, { sa, 1, 1, sc, sd, se }, x)
        + (sy[14] + sy[15] + sy[16] + sy[17] + sy[18] - sy[13] * sy[8]) * Log(b, sb)
        + (-sy[14] - sy[15] - sy[16] - sy[17] - sy[18] + sy[13] * sy[8]) * Log(-x, sb) };
    if (c != x) {
        res += (sy[0] + sy[4] - G({ 0, a / b, x / b }, 1) - G({ a / b, 0, x / b }, 1))
            * G({ c, d, e }, { sc, sd, se }, x);
    }
    if (d != x) {
        res += (sy[10] + sy[11] + sy[12] - G({ 0, a / b, c / b, x / b }, 1)
                   - G({ 0, c / b, a / b, x / b }, 1) - G({ a / b, 0, c / b, x / b }, 1))
            * G({ d, e }, { sd, se }, x);
    }
    if (e != x) {
        res += (sy[15] + sy[16] + sy[17] + sy[18]
                   - G({ 0, a / b, c / b, d / b, x / b }, 1)
                   - G({ 0, c / b, a / b, d / b, x / b }, 1)
                   - G({ 0, c / b, d / b, a / b, x / b }, 1)
                   - G({ a / b, 0, c / b, d / b, x / b }, 1))
            * G({ e }, { se }, x);
    }
    return res;
}
complex<double> G6_ab0cde_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (a == b) {
        // aacde
        const vector<complex<double>> sy = { G({ 0, 1, c / a }, 1),
            G({ 0, d, e }, { 1, sd, se }, x), G({ 0, c / a, 1 }, 1),
            G({ 0, 1, c / a, d / a }, 1), G({ 0, c / a, 1, d / a }, 1),
            G({ 0, c / a, d / a, 1 }, 1), G({ 0, 1, c / a, d / a, e / a }, 1),
            G({ 0, c / a, 1, d / a, e / a }, 1), G({ 0, c / a, d / a, 1, e / a }, 1),
            G({ 0, c / a, d / a, e / a, 1 }, 1) };
        complex<double> res { -(sy[0] * sy[1]) - sy[1] * sy[2]
            - 2. * G({ 0, 0, 1, c / a, d / a, e / a }, 1)
            - 2. * G({ 0, 0, c / a, 1, d / a, e / a }, 1)
            - 2. * G({ 0, 0, c / a, d / a, 1, e / a }, 1)
            - 2. * G({ 0, 0, c / a, d / a, e / a, 1 }, 1)
            - G({ 0, 1, 0, c / a, d / a, e / a }, 1)
            - G({ 0, 1, c / a, 0, d / a, e / a }, 1)
            - G({ 0, 1, c / a, d / a, 0, e / a }, 1)
            - G({ 0, 1, c / a, d / a, e / a, x / a }, 1)
            - G({ 0, c / a, 0, 1, d / a, e / a }, 1)
            - G({ 0, c / a, 0, d / a, 1, e / a }, 1)
            - G({ 0, c / a, 0, d / a, e / a, 1 }, 1)
            - G({ 0, c / a, 1, 0, d / a, e / a }, 1)
            - G({ 0, c / a, 1, d / a, 0, e / a }, 1)
            - G({ 0, c / a, 1, d / a, e / a, x / a }, 1)
            - G({ 0, c / a, d / a, 0, 1, e / a }, 1)
            - G({ 0, c / a, d / a, 0, e / a, 1 }, 1)
            - G({ 0, c / a, d / a, 1, 0, e / a }, 1)
            - G({ 0, c / a, d / a, 1, e / a, x / a }, 1)
            - G({ 0, c / a, d / a, e / a, 1, x / a }, 1)
            - G({ 0, c / a, d / a, e / a, x / a, 1 }, 1)
            + (-sy[3] - sy[4] - sy[5]) * G({ 0, e }, { 1, se }, x)
            - G({ x / a, 1 }, 1) * G({ 0, c, d, e }, { 1, sc, sd, se }, x)
            + G({ x / a }, 1) * G({ a, 0, c, d, e }, { sa, 1, sc, sd, se }, x)
            + G({ 0, a, 0, c, d, e }, { 1, sa, 1, sc, sd, se }, x)
            + (sy[9] + sy[6] + sy[7] + sy[8]) * Log(a, sa)
            + (-sy[9] - sy[6] - sy[7] - sy[8]) * Log(-x, sa) };
        if (c != x) {
            res += (sy[0] + sy[2] - G({ 0, 1, x / a }, 1) - G({ 0, x / a, 1 }, 1))
                * G({ c, d, e }, { sc, sd, se }, x);
        }
        if (d != x) {
            res += (sy[3] + sy[4] + sy[5] - G({ 0, 1, c / a, x / a }, 1)
                       - G({ 0, c / a, 1, x / a }, 1) - G({ 0, c / a, x / a, 1 }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (sy[9] + sy[6] + sy[7] + sy[8] - G({ 0, 1, c / a, d / a, x / a }, 1)
                       - G({ 0, c / a, 1, d / a, x / a }, 1)
                       - G({ 0, c / a, d / a, 1, x / a }, 1)
                       - G({ 0, c / a, d / a, x / a, 1 }, 1))
                * G({ e }, { se }, x);
        }
        return res;

    } else { // abcde
        const vector<complex<double>> sy
            = { G({ b / a, 0, c / a }, 1), G({ b / a, 0, c / a, d / a }, 1),
                  G({ b / a }, 1), G({ b / a, 0, c / a, d / a, e / a }, 1) };
        complex<double> res { G({ 0, b / a, 0, c / a, d / a, e / a }, 1)
            + 2. * G({ b / a, 0, 0, c / a, d / a, e / a }, 1)
            + G({ b / a, 0, c / a, 0, d / a, e / a }, 1)
            + G({ b / a, 0, c / a, d / a, 0, e / a }, 1)
            + G({ b / a, 0, c / a, d / a, e / a, x / a }, 1)
            + sy[1] * G({ 0, e }, { 1, se }, x) + sy[0] * G({ 0, d, e }, { 1, sd, se }, x)
            + G({ b / a, x / a }, 1) * G({ 0, c, d, e }, { 1, sc, sd, se }, x)
            + sy[2] * G({ 0, 0, c, d, e }, { 1, 1, sc, sd, se }, x)
            + G({ 0, b, 0, c, d, e }, { 1, sb, 1, sc, sd, se }, x) - sy[3] * Log(a, sa)
            + sy[3] * Log(-x, sa) };
        if (b != x) {
            res += (-sy[2] + G({ x / a }, 1))
                * G({ b, 0, c, d, e }, { sb, 1, sc, sd, se }, x);
        }
        if (c != x) {
            res += (-sy[0] + G({ b / a, 0, x / a }, 1))
                * G({ c, d, e }, { sc, sd, se }, x);
        }
        if (d != x) {
            res += (-sy[1] + G({ b / a, 0, c / a, x / a }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (-sy[3] + G({ b / a, 0, c / a, d / a, x / a }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}

complex<double> G6_abc0de_e(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { G({ a, b }, { sa, sb }, x),
        G({ d / e, x / e }, 1), G({ d / e, 0, c / e }, 1),
        G({ a, b, c }, { sa, sb, sc }, x), G({ d / e, 0, x / e }, 1), G({ 0, d / e }, 1),
        G({ 0, a, b, c }, { 1, sa, sb, sc }, x), G({ a, 0, b, c }, { sa, 1, sb, sc }, x),
        G({ a, b, 0, c }, { sa, sb, 1, sc }, x), G({ a }, { sa }, x),
        G({ d / e, 0, c / e, b / e }, 1), G({ d / e }, 1), Log(e, se),
        G({ d / e, 0, c / e, x / e }, 1), G({ a, b, c, 0, d }, { sa, sb, sc, 1, sd }, x),
        G({ d / e, 0, c / e, b / e, a / e }, 1) };
    complex<double> res { sy[12] * (sy[15] + sy[3] * sy[5])
        + sy[5] * (sy[6] + sy[7] + sy[8]) + sy[3] * (sy[4] - G({ 0, 0, d / e }, 1))
        - G({ 0, d / e, 0, c / e, b / e, a / e }, 1)
        - 2. * G({ d / e, 0, 0, c / e, b / e, a / e }, 1)
        - G({ d / e, 0, c / e, 0, b / e, a / e }, 1)
        - G({ d / e, 0, c / e, b / e, 0, a / e }, 1)
        - G({ d / e, 0, c / e, b / e, a / e, x / e }, 1)
        + sy[14] * (-G({ x / e }, 1) + G({ e }, { se }, x))
        - sy[10] * G({ 0, a }, { 1, sa }, x)
        + sy[9]
            * (-(sy[10] * sy[12]) + sy[15] + G({ 0, d / e, 0, c / e, b / e }, 1)
                + 2. * G({ d / e, 0, 0, c / e, b / e }, 1)
                + G({ d / e, 0, c / e, 0, b / e }, 1) + sy[2] * G({ 0, b }, { 1, sb }, x))
        - sy[0] * sy[1] * G({ 0, c }, { 1, sc }, x)
        + sy[0]
            * (-sy[13] + sy[12] * sy[2] - G({ 0, d / e, 0, c / e }, 1)
                - 2. * G({ d / e, 0, 0, c / e }, 1)
                - sy[11] * G({ 0, 0, c }, { 1, 1, sc }, x))
        - sy[2] * G({ 0, b, a }, { 1, sb, sa }, x)
        + sy[1]
            * (sy[6] + sy[7] + sy[8] + sy[9] * G({ 0, c, b }, { 1, sc, sb }, x)
                - G({ 0, c, b, a }, { 1, sc, sb, sa }, x))
        - G({ 0, a, b, c, 0, d }, { 1, sa, sb, sc, 1, sd }, x)
        - G({ a, 0, b, c, 0, d }, { sa, 1, sb, sc, 1, sd }, x)
        - G({ a, b, 0, c, 0, d }, { sa, sb, 1, sc, 1, sd }, x)
        - 2. * G({ a, b, c, 0, 0, d }, { sa, sb, sc, 1, 1, sd }, x)
        + (sy[9] * sy[10] - sy[15] - sy[0] * sy[2] - sy[3] * sy[5]
              + sy[11] * (sy[12] * sy[3] + sy[6] + sy[7] + sy[8]))
            * Log(-x, se)
        + sy[11]
            * (sy[14] + sy[12] * (-sy[6] - sy[7] - sy[8])
                + sy[9] * G({ 0, 0, c, b }, { 1, 1, sc, sb }, x)
                - G({ 0, 0, c, b, a }, { 1, 1, sc, sb, sa }, x)
                - (sy[3] * pow(sy[12], 2.)) / 2.
                + sy[3] * (G({ 0, 0 }, { 1, 1 }, x) + 2. * Zeta(2))) };
    if (b != x) {
        res += (-(sy[9] * sy[10]) + sy[9] * sy[13]) * G({ b }, { sb }, x)
            + (sy[10] - sy[13]) * G({ b, a }, { sb, sa }, x);
    }
    if (c != x) {
        res += (sy[0] * sy[2] - sy[0] * sy[4]) * G({ c }, { sc }, x)
            + (-(sy[9] * sy[2]) + sy[9] * sy[4]) * G({ c, b }, { sc, sb }, x)
            + (sy[2] - sy[4]) * G({ c, b, a }, { sc, sb, sa }, x);
    }
    return res;
}
complex<double> G6_abc0de_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (d == e) {
        const vector<complex<double>> sy = { G({ a, b }, { sa, sb }, x),
            G({ x / d, 1 }, 1), G({ 0, a, b }, { 1, sa, sb }, x), G({ 0, c / d, 1 }, 1),
            G({ 0, b, a }, { 1, sb, sa }, x), G({ 0, 1, c / d }, 1),
            G({ a, 0, b }, { sa, 1, sb }, x), G({ a, b, c }, { sa, sb, sc }, x),
            G({ a, b, d }, { sa, sb, sd }, x), G({ 0, c / d, b / d, 1 }, 1),
            G({ 0, c / d, b / d, a / d }, 1), G({ 0, 1, c / d, b / d }, 1),
            G({ 0, 1, c / d, x / d }, 1), G({ 0, c / d, 1, b / d }, 1),
            G({ 0, c / d, 1, x / d }, 1), G({ 0, c / d, x / d, 1 }, 1),
            G({ a }, { sa }, x), G({ a, b, c, d }, { sa, sb, sc, sd }, x),
            G({ d }, { sd }, x), G({ 0, c / d, b / d, a / d, 1 }, 1),
            G({ 0, 1, x / d }, 1), G({ 0, x / d, 1 }, 1) };
        complex<double> res { -(sy[18] * sy[19]) - sy[2] * sy[3] - sy[3] * sy[4]
            + sy[5] * (-sy[2] - sy[4] - sy[6] - sy[7]) + sy[3] * (-sy[6] - sy[8])
            + sy[0]
                * (sy[11] - sy[12] + sy[13] - sy[14] - sy[15] + G({ 0, c / d, 0, 1 }, 1))
            + sy[18] * G({ 0, c / d, b / d, a / d, x / d }, 1)
            + G({ 0, c / d, b / d, a / d, 0, 1 }, 1)
            - G({ 0, c / d, b / d, a / d, x / d, 1 }, 1)
            + sy[16]
                * (sy[19] - G({ 0, c / d, b / d, 0, 1 }, 1)
                    + (sy[3] + sy[5]) * G({ 0, b }, { 1, sb }, x))
            - sy[0] * sy[1] * G({ 0, c }, { 1, sc }, x)
            + sy[10] * G({ 0, d }, { 1, sd }, x)
            + (sy[9] - sy[10]) * G({ a, d }, { sa, sd }, x)
            + G({ 0, c / d, b / d }, 1) * (sy[8] - G({ a, 0, d }, { sa, 1, sd }, x))
            + sy[1]
                * (sy[16] * G({ 0, c, b }, { 1, sc, sb }, x)
                    + G({ 0, a, b, c }, { 1, sa, sb, sc }, x)
                    - G({ 0, c, b, a }, { 1, sc, sb, sa }, x)
                    + G({ a, 0, b, c }, { sa, 1, sb, sc }, x)
                    + G({ a, b, 0, c }, { sa, sb, 1, sc }, x))
            + G({ 0, c / d }, 1) * (-sy[17] + G({ a, b, 0, d }, { sa, sb, 1, sd }, x))
            + G({ a, b, c, 0, 0, d }, { sa, sb, sc, 1, 1, sd }, x) - sy[17] * Zeta(2)
            + sy[7] * (sy[20] + sy[21] + Zeta(3)) };
        if (b != x) {
            res += (-(sy[9] * sy[16]) - sy[11] * sy[16] + sy[12] * sy[16]
                       - sy[13] * sy[16] + sy[14] * sy[16] + sy[15] * sy[16])
                    * G({ b }, { sb }, x)
                + (sy[9] + sy[11] - sy[12] + sy[13] - sy[14] - sy[15])
                    * G({ b, a }, { sb, sa }, x);
        }
        if (c != x) {
            res += (-(sy[0] * sy[20]) - sy[0] * sy[21] + sy[0] * sy[3] + sy[0] * sy[5])
                    * G({ c }, { sc }, x)
                + (sy[16] * sy[20] + sy[16] * sy[21] - sy[16] * sy[3] - sy[16] * sy[5])
                    * G({ c, b }, { sc, sb }, x)
                + (-sy[20] - sy[21] + sy[3] + sy[5]) * G({ c, b, a }, { sc, sb, sa }, x);
        }
        return res;
    } else {
        const vector<complex<double>> sy = { G({ a, b }, { sa, sb }, x),
            G({ e / d, x / d }, 1), G({ 0, a, b }, { 1, sa, sb }, x),
            G({ 0, e / d, c / d }, 1), G({ a, 0, b }, { sa, 1, sb }, x),
            G({ a, b, c }, { sa, sb, sc }, x), G({ 0, c / d, e / d }, 1),
            G({ a, b, e }, { sa, sb, se }, x), G({ e / d, 0, c / d }, 1),
            G({ e / d, 0, x / d }, 1), G({ 0, c / d, b / d, a / d }, 1),
            G({ a, e }, { sa, se }, x), G({ 0, a }, { 1, sa }, x),
            G({ 0, c / d, b / d, e / d }, 1), G({ a }, { sa }, x),
            G({ 0, a, b, c }, { 1, sa, sb, sc }, x),
            G({ a, 0, b, c }, { sa, 1, sb, sc }, x),
            G({ a, b, 0, c }, { sa, sb, 1, sc }, x),
            G({ a, b, c, e }, { sa, sb, sc, se }, x), G({ 0, c / d, e / d, b / d }, 1),
            G({ 0, e / d, c / d, b / d }, 1), G({ e / d, 0, c / d, b / d }, 1),
            G({ e / d }, 1), Log(d, sd), G({ e / d, 0, c / d, x / d }, 1),
            G({ 0, c / d, b / d, a / d, e / d }, 1),
            G({ 0, c / d, b / d, e / d, a / d }, 1),
            G({ 0, c / d, e / d, b / d, a / d }, 1),
            G({ 0, e / d, c / d, b / d, a / d }, 1),
            G({ e / d, 0, c / d, b / d, a / d }, 1) };
        complex<double> res { -(sy[10] * sy[11]) + sy[11] * sy[13] + sy[12] * sy[13]
            + sy[12] * (sy[19] + sy[20] + sy[21])
            + sy[23] * (-sy[25] - sy[26] - sy[27] - sy[28] - sy[29]) - sy[2] * sy[3]
            + sy[3] * (-sy[4] - sy[5]) + sy[6] * (-sy[2] - sy[4] - sy[7])
            + sy[18] * G({ 0, e / d }, 1) + sy[5] * (-sy[9] - G({ 0, 0, e / d }, 1))
            + 2. * G({ 0, 0, c / d, b / d, a / d, e / d }, 1)
            + 2. * G({ 0, 0, c / d, b / d, e / d, a / d }, 1)
            + 2. * G({ 0, 0, c / d, e / d, b / d, a / d }, 1)
            + 2. * G({ 0, 0, e / d, c / d, b / d, a / d }, 1)
            + G({ 0, c / d, 0, b / d, a / d, e / d }, 1)
            + G({ 0, c / d, 0, b / d, e / d, a / d }, 1)
            + G({ 0, c / d, 0, e / d, b / d, a / d }, 1)
            + G({ 0, c / d, b / d, 0, a / d, e / d }, 1)
            + G({ 0, c / d, b / d, 0, e / d, a / d }, 1)
            + G({ 0, c / d, b / d, a / d, 0, e / d }, 1)
            + G({ 0, c / d, b / d, a / d, e / d, x / d }, 1)
            + G({ 0, c / d, b / d, e / d, 0, a / d }, 1)
            + G({ 0, c / d, b / d, e / d, a / d, x / d }, 1)
            + G({ 0, c / d, e / d, 0, b / d, a / d }, 1)
            + G({ 0, c / d, e / d, b / d, 0, a / d }, 1)
            + G({ 0, c / d, e / d, b / d, a / d, x / d }, 1)
            + 2. * G({ 0, e / d, 0, c / d, b / d, a / d }, 1)
            + G({ 0, e / d, c / d, 0, b / d, a / d }, 1)
            + G({ 0, e / d, c / d, b / d, 0, a / d }, 1)
            + G({ 0, e / d, c / d, b / d, a / d, x / d }, 1)
            + 2. * G({ e / d, 0, 0, c / d, b / d, a / d }, 1)
            + G({ e / d, 0, c / d, 0, b / d, a / d }, 1)
            + G({ e / d, 0, c / d, b / d, 0, a / d }, 1)
            + G({ e / d, 0, c / d, b / d, a / d, x / d }, 1)
            + sy[14]
                * ((sy[13] + sy[19] + sy[20] + sy[21]) * sy[23] - sy[26] - sy[27] - sy[28]
                    - sy[29] - 2. * G({ 0, 0, c / d, b / d, e / d }, 1)
                    - 2. * G({ 0, 0, c / d, e / d, b / d }, 1)
                    - 2. * G({ 0, 0, e / d, c / d, b / d }, 1)
                    - G({ 0, c / d, 0, b / d, e / d }, 1)
                    - G({ 0, c / d, 0, e / d, b / d }, 1)
                    - G({ 0, c / d, b / d, 0, e / d }, 1)
                    - G({ 0, c / d, e / d, 0, b / d }, 1)
                    - 2. * G({ 0, e / d, 0, c / d, b / d }, 1)
                    - G({ 0, e / d, c / d, 0, b / d }, 1)
                    - 2. * G({ e / d, 0, 0, c / d, b / d }, 1)
                    - G({ e / d, 0, c / d, 0, b / d }, 1)
                    - sy[8] * G({ 0, b }, { 1, sb }, x))
            + sy[0] * sy[1] * G({ 0, c }, { 1, sc }, x)
            + sy[10] * G({ 0, e }, { 1, se }, x)
            + sy[0]
                * (sy[19] + sy[20] + sy[24] + sy[23] * (-sy[3] - sy[6] - sy[8])
                    + 2. * G({ 0, 0, c / d, e / d }, 1)
                    + 2. * G({ 0, 0, e / d, c / d }, 1) + G({ 0, c / d, 0, e / d }, 1)
                    + 2. * G({ 0, e / d, 0, c / d }, 1)
                    + 2. * G({ e / d, 0, 0, c / d }, 1)
                    + sy[22] * G({ 0, 0, c }, { 1, 1, sc }, x))
            + sy[8] * G({ 0, b, a }, { 1, sb, sa }, x)
            + G({ 0, c / d, b / d }, 1) * (sy[7] - G({ a, 0, e }, { sa, 1, se }, x))
            + sy[1]
                * (-sy[15] - sy[16] - sy[17] - sy[14] * G({ 0, c, b }, { 1, sc, sb }, x)
                    + G({ 0, c, b, a }, { 1, sc, sb, sa }, x))
            + G({ 0, c / d }, 1) * (-sy[18] + G({ a, b, 0, e }, { sa, sb, 1, se }, x))
            + G({ a, b, c, 0, 0, e }, { sa, sb, sc, 1, 1, se }, x)
            + (sy[14] * (-sy[13] - sy[19] - sy[20] - sy[21]) + sy[25] + sy[26] + sy[27]
                  + sy[28] + sy[29]
                  + sy[22] * (-sy[15] - sy[16] - sy[17] - sy[23] * sy[5])
                  + sy[0] * (sy[3] + sy[6] + sy[8]))
                * Log(-x, sd)
            + sy[22]
                * ((sy[15] + sy[16] + sy[17]) * sy[23]
                    - sy[14] * G({ 0, 0, c, b }, { 1, 1, sc, sb }, x)
                    + G({ 0, 0, c, b, a }, { 1, 1, sc, sb, sa }, x)
                    - G({ a, b, c, 0, e }, { sa, sb, sc, 1, se }, x)
                    + (sy[5] * pow(sy[23], 2.)) / 2.
                    + sy[5] * (-G({ 0, 0 }, { 1, 1 }, x) - 2. * Zeta(2))) };
        if (b != x) {
            res += (sy[14] * sy[21] - sy[14] * sy[24]) * G({ b }, { sb }, x)
                + (-sy[21] + sy[24]) * G({ b, a }, { sb, sa }, x);
        }
        if (c != x) {
            res += (sy[0] * sy[9] - sy[0] * sy[8]) * G({ c }, { sc }, x)
                + (-(sy[9] * sy[14]) + sy[14] * sy[8]) * G({ c, b }, { sc, sb }, x)
                + (sy[9] - sy[8]) * G({ c, b, a }, { sc, sb, sa }, x);
        }
        if (e != x) {
            res += (-sy[25] + G({ 0, c / d, b / d, a / d, x / d }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_abc0de_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    const vector<complex<double>> sy = { G({ 0, d / c, e / c }, 1),
        G({ a, b, e }, { sa, sb, se }, x), G({ 0, d / c, b / c }, 1),
        G({ 0, d, e }, { 1, sd, se }, x), G({ a, d, e }, { sa, sd, se }, x),
        G({ a, 0, e }, { sa, 1, se }, x), G({ b / c, 0, d / c }, 1),
        G({ b / c, 0, a / c }, 1), G({ a, e }, { sa, se }, x),
        G({ 0, b / c, d / c, a / c }, 1), G({ 0, a }, { 1, sa }, x),
        G({ 0, b / c, d / c, e / c }, 1), G({ a, b }, { sa, sb }, x),
        G({ 0, d / c, e / c, b / c }, 1), G({ a, 0, d, e }, { sa, 1, sd, se }, x),
        G({ a, b, d, e }, { sa, sb, sd, se }, x), G({ 0, d / c, b / c, e / c }, 1),
        G({ b / c, 0, d / c, e / c }, 1), G({ 0, d / c, b / c, a / c }, 1),
        G({ b / c, 0, d / c, a / c }, 1), G({ 0, b / c, a / c, d / c }, 1),
        G({ b / c, 0, a / c, d / c }, 1), G({ b / c, a / c, 0, d / c }, 1),
        G({ a }, { sa }, x), G({ 0, b / c, d / c, e / c, a / c }, 1),
        G({ 0, d / c, b / c, e / c, a / c }, 1), G({ 0, d / c, e / c, b / c, a / c }, 1),
        G({ b / c, 0, d / c, e / c, a / c }, 1), G({ 0, b / c, a / c, d / c, e / c }, 1),
        G({ 0, b / c, d / c, a / c, e / c }, 1), G({ 0, d / c, b / c, a / c, e / c }, 1),
        G({ b / c, 0, a / c, d / c, e / c }, 1), G({ b / c, 0, d / c, a / c, e / c }, 1),
        G({ b / c, a / c, 0, d / c, e / c }, 1) };
    complex<double> res { sy[10] * sy[11] + sy[10] * (sy[13] + sy[16] + sy[17])
        + sy[1] * sy[2] + sy[5] * (-sy[2] - sy[6]) + sy[4] * (sy[6] - sy[7])
        - sy[9] * sy[8] + (sy[11] + sy[16] + sy[17] - sy[18] - sy[19]) * sy[8]
        + (-sy[14] + sy[15]) * G({ 0, b / c }, 1)
        + (sy[3] - sy[4]) * G({ 0, b / c, a / c }, 1)
        + (sy[4] - sy[5]) * G({ 0, b / c, d / c }, 1)
        + sy[3] * (sy[7] + G({ b / c, a / c, x / c }, 1))
        + sy[12]
            * (sy[13] + 2. * G({ 0, 0, d / c, e / c }, 1) + G({ 0, d / c, 0, e / c }, 1))
        + sy[23]
            * (-sy[24] - sy[25] - sy[26] - sy[27]
                - 2. * G({ 0, 0, b / c, d / c, e / c }, 1)
                - 2. * G({ 0, 0, d / c, b / c, e / c }, 1)
                - 2. * G({ 0, 0, d / c, e / c, b / c }, 1)
                - 2. * G({ 0, b / c, 0, d / c, e / c }, 1)
                - G({ 0, b / c, d / c, 0, e / c }, 1)
                - G({ 0, d / c, 0, b / c, e / c }, 1)
                - G({ 0, d / c, 0, e / c, b / c }, 1)
                - G({ 0, d / c, b / c, 0, e / c }, 1)
                - G({ 0, d / c, e / c, 0, b / c }, 1)
                - 2. * G({ b / c, 0, 0, d / c, e / c }, 1)
                - G({ b / c, 0, d / c, 0, e / c }, 1))
        + 2. * G({ 0, 0, b / c, a / c, d / c, e / c }, 1)
        + 2. * G({ 0, 0, b / c, d / c, a / c, e / c }, 1)
        + 2. * G({ 0, 0, b / c, d / c, e / c, a / c }, 1)
        + 2. * G({ 0, 0, d / c, b / c, a / c, e / c }, 1)
        + 2. * G({ 0, 0, d / c, b / c, e / c, a / c }, 1)
        + 2. * G({ 0, 0, d / c, e / c, b / c, a / c }, 1)
        + 2. * G({ 0, b / c, 0, a / c, d / c, e / c }, 1)
        + 2. * G({ 0, b / c, 0, d / c, a / c, e / c }, 1)
        + 2. * G({ 0, b / c, 0, d / c, e / c, a / c }, 1)
        + 2. * G({ 0, b / c, a / c, 0, d / c, e / c }, 1)
        + G({ 0, b / c, a / c, d / c, 0, e / c }, 1)
        + G({ 0, b / c, a / c, d / c, e / c, x / c }, 1)
        + G({ 0, b / c, d / c, 0, a / c, e / c }, 1)
        + G({ 0, b / c, d / c, 0, e / c, a / c }, 1)
        + G({ 0, b / c, d / c, a / c, 0, e / c }, 1)
        + G({ 0, b / c, d / c, a / c, e / c, x / c }, 1)
        + G({ 0, b / c, d / c, e / c, 0, a / c }, 1)
        + G({ 0, b / c, d / c, e / c, a / c, x / c }, 1)
        + G({ 0, d / c, 0, b / c, a / c, e / c }, 1)
        + G({ 0, d / c, 0, b / c, e / c, a / c }, 1)
        + G({ 0, d / c, 0, e / c, b / c, a / c }, 1)
        + G({ 0, d / c, b / c, 0, a / c, e / c }, 1)
        + G({ 0, d / c, b / c, 0, e / c, a / c }, 1)
        + G({ 0, d / c, b / c, a / c, 0, e / c }, 1)
        + G({ 0, d / c, b / c, a / c, e / c, x / c }, 1)
        + G({ 0, d / c, b / c, e / c, 0, a / c }, 1)
        + G({ 0, d / c, b / c, e / c, a / c, x / c }, 1)
        + G({ 0, d / c, e / c, 0, b / c, a / c }, 1)
        + G({ 0, d / c, e / c, b / c, 0, a / c }, 1)
        + G({ 0, d / c, e / c, b / c, a / c, x / c }, 1)
        + 2. * G({ b / c, 0, 0, a / c, d / c, e / c }, 1)
        + 2. * G({ b / c, 0, 0, d / c, a / c, e / c }, 1)
        + 2. * G({ b / c, 0, 0, d / c, e / c, a / c }, 1)
        + 2. * G({ b / c, 0, a / c, 0, d / c, e / c }, 1)
        + G({ b / c, 0, a / c, d / c, 0, e / c }, 1)
        + G({ b / c, 0, a / c, d / c, e / c, x / c }, 1)
        + G({ b / c, 0, d / c, 0, a / c, e / c }, 1)
        + G({ b / c, 0, d / c, 0, e / c, a / c }, 1)
        + G({ b / c, 0, d / c, a / c, 0, e / c }, 1)
        + G({ b / c, 0, d / c, a / c, e / c, x / c }, 1)
        + G({ b / c, 0, d / c, e / c, 0, a / c }, 1)
        + G({ b / c, 0, d / c, e / c, a / c, x / c }, 1)
        + 2. * G({ b / c, a / c, 0, 0, d / c, e / c }, 1)
        + G({ b / c, a / c, 0, d / c, 0, e / c }, 1)
        + G({ b / c, a / c, 0, d / c, e / c, x / c }, 1)
        + (sy[9] + sy[18] + sy[19] + sy[20] + sy[21] + sy[22]) * G({ 0, e }, { 1, se }, x)
        + sy[0]
            * (-sy[1] - G({ 0, a, b }, { 1, sa, sb }, x)
                - G({ a, 0, b }, { sa, 1, sb }, x))
        + G({ b / c, a / c }, 1) * (-sy[14] + G({ 0, 0, d, e }, { 1, 1, sd, se }, x))
        + G({ 0, d / c }, 1) * (-sy[15] + G({ a, b, 0, e }, { sa, sb, 1, se }, x))
        + G({ b / c }, 1)
            * (-G({ a, 0, 0, d, e }, { sa, 1, 1, sd, se }, x)
                + G({ a, b, 0, d, e }, { sa, sb, 1, sd, se }, x))
        + G({ a, b, 0, 0, d, e }, { sa, sb, 1, 1, sd, se }, x)
        + (-(sy[0] * sy[12]) + (sy[11] + sy[13] + sy[16] + sy[17]) * sy[23] - sy[24]
              - sy[25] - sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31] - sy[32]
              - sy[33])
            * Log(c, sc)
        + (sy[0] * sy[12] + (-sy[11] - sy[13] - sy[16] - sy[17]) * sy[23] + sy[24]
              + sy[25] + sy[26] + sy[27] + sy[28] + sy[29] + sy[30] + sy[31] + sy[32]
              + sy[33])
            * Log(-x, sc) };
    if (d != x) {
        res += (-sy[20] - sy[21] - sy[22] + G({ 0, b / c, a / c, x / c }, 1)
                   + G({ b / c, 0, a / c, x / c }, 1) + G({ b / c, a / c, 0, x / c }, 1))
            * G({ d, e }, { sd, se }, x);
    }
    if (e != x) {
        res += (-sy[28] - sy[29] - sy[30] - sy[31] - sy[32] - sy[33]
                   + G({ 0, b / c, a / c, d / c, x / c }, 1)
                   + G({ 0, b / c, d / c, a / c, x / c }, 1)
                   + G({ 0, d / c, b / c, a / c, x / c }, 1)
                   + G({ b / c, 0, a / c, d / c, x / c }, 1)
                   + G({ b / c, 0, d / c, a / c, x / c }, 1)
                   + G({ b / c, a / c, 0, d / c, x / c }, 1))
            * G({ e }, { se }, x);
    }
    return res;
}
complex<double> G6_abc0de_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (b == c) {
        // abbde
        const vector<complex<double>> sy = { G({ 0, d / b, 1 }, 1),
            G({ a, 0, e }, { sa, 1, se }, x), G({ 0, d, e }, { 1, sd, se }, x),
            G({ a, d, e }, { sa, sd, se }, x), G({ 0, a / b, 1 }, 1),
            G({ a, e }, { sa, se }, x), G({ 0, 1, d / b, a / b }, 1),
            G({ 0, a }, { 1, sa }, x), G({ 0, 1, d / b, e / b }, 1),
            G({ 0, d / b, 1, e / b }, 1), G({ 0, d / b, e / b, 1 }, 1),
            G({ 0, d / b, 1, a / b }, 1), G({ 0, d / b, a / b, 1 }, 1),
            G({ 0, 1, a / b, d / b }, 1), G({ 0, a / b, 1, d / b }, 1),
            G({ 0, a / b, d / b, 1 }, 1), G({ a / b, 0, 1, d / b }, 1),
            G({ a / b, 0, d / b, 1 }, 1), G({ b, 0, d, e }, { sb, 1, sd, se }, x),
            G({ a }, { sa }, x), G({ 0, 1, d / b, e / b, a / b }, 1),
            G({ 0, d / b, 1, e / b, a / b }, 1), G({ 0, d / b, e / b, 1, a / b }, 1),
            G({ 0, d / b, e / b, a / b, 1 }, 1), G({ 0, 1, a / b, d / b, e / b }, 1),
            G({ 0, 1, d / b, a / b, e / b }, 1), G({ 0, a / b, 1, d / b, e / b }, 1),
            G({ 0, a / b, d / b, 1, e / b }, 1), G({ 0, a / b, d / b, e / b, 1 }, 1),
            G({ 0, d / b, 1, a / b, e / b }, 1), G({ 0, d / b, a / b, 1, e / b }, 1),
            G({ 0, d / b, a / b, e / b, 1 }, 1), G({ a / b, 0, 1, d / b, e / b }, 1),
            G({ a / b, 0, d / b, 1, e / b }, 1), G({ a / b, 0, d / b, e / b, 1 }, 1) };
        complex<double> res { -(sy[0] * sy[1]) + sy[0] * sy[3] - sy[3] * sy[4]
            - sy[5] * sy[6] + (sy[9] + sy[10]) * sy[7] + sy[7] * sy[8]
            + sy[5] * (sy[9] + sy[10] - sy[11] - sy[12] + sy[8])
            - sy[18] * G({ a / b, x / b }, 1) + (sy[2] - sy[3]) * G({ 0, 1, a / b }, 1)
            + (-sy[1] + sy[3]) * G({ 0, 1, d / b }, 1)
            + sy[2] * (sy[4] + G({ a / b, x / b, 1 }, 1))
            + sy[19]
                * (-sy[20] - sy[21] - sy[22] - sy[23]
                    - 2. * G({ 0, 0, 1, d / b, e / b }, 1)
                    - 2. * G({ 0, 0, d / b, 1, e / b }, 1)
                    - 2. * G({ 0, 0, d / b, e / b, 1 }, 1)
                    - G({ 0, 1, 0, d / b, e / b }, 1) - G({ 0, 1, d / b, 0, e / b }, 1)
                    - G({ 0, d / b, 0, 1, e / b }, 1) - G({ 0, d / b, 0, e / b, 1 }, 1)
                    - G({ 0, d / b, 1, 0, e / b }, 1))
            + 2. * G({ 0, 0, 1, a / b, d / b, e / b }, 1)
            + 2. * G({ 0, 0, 1, d / b, a / b, e / b }, 1)
            + 2. * G({ 0, 0, 1, d / b, e / b, a / b }, 1)
            + 2. * G({ 0, 0, a / b, 1, d / b, e / b }, 1)
            + 2. * G({ 0, 0, a / b, d / b, 1, e / b }, 1)
            + 2. * G({ 0, 0, a / b, d / b, e / b, 1 }, 1)
            + 2. * G({ 0, 0, d / b, 1, a / b, e / b }, 1)
            + 2. * G({ 0, 0, d / b, 1, e / b, a / b }, 1)
            + 2. * G({ 0, 0, d / b, a / b, 1, e / b }, 1)
            + 2. * G({ 0, 0, d / b, a / b, e / b, 1 }, 1)
            + 2. * G({ 0, 0, d / b, e / b, 1, a / b }, 1)
            + 2. * G({ 0, 0, d / b, e / b, a / b, 1 }, 1)
            + G({ 0, 1, 0, a / b, d / b, e / b }, 1)
            + G({ 0, 1, 0, d / b, a / b, e / b }, 1)
            + G({ 0, 1, 0, d / b, e / b, a / b }, 1)
            + G({ 0, 1, a / b, 0, d / b, e / b }, 1)
            + G({ 0, 1, a / b, d / b, 0, e / b }, 1)
            + G({ 0, 1, a / b, d / b, e / b, x / b }, 1)
            + G({ 0, 1, d / b, 0, a / b, e / b }, 1)
            + G({ 0, 1, d / b, 0, e / b, a / b }, 1)
            + G({ 0, 1, d / b, a / b, 0, e / b }, 1)
            + G({ 0, 1, d / b, a / b, e / b, x / b }, 1)
            + G({ 0, 1, d / b, e / b, 0, a / b }, 1)
            + G({ 0, 1, d / b, e / b, a / b, x / b }, 1)
            + 2. * G({ 0, a / b, 0, 1, d / b, e / b }, 1)
            + 2. * G({ 0, a / b, 0, d / b, 1, e / b }, 1)
            + 2. * G({ 0, a / b, 0, d / b, e / b, 1 }, 1)
            + G({ 0, a / b, 1, 0, d / b, e / b }, 1)
            + G({ 0, a / b, 1, d / b, 0, e / b }, 1)
            + G({ 0, a / b, 1, d / b, e / b, x / b }, 1)
            + G({ 0, a / b, d / b, 0, 1, e / b }, 1)
            + G({ 0, a / b, d / b, 0, e / b, 1 }, 1)
            + G({ 0, a / b, d / b, 1, 0, e / b }, 1)
            + G({ 0, a / b, d / b, 1, e / b, x / b }, 1)
            + G({ 0, a / b, d / b, e / b, 1, x / b }, 1)
            + G({ 0, a / b, d / b, e / b, x / b, 1 }, 1)
            + G({ 0, d / b, 0, 1, a / b, e / b }, 1)
            + G({ 0, d / b, 0, 1, e / b, a / b }, 1)
            + G({ 0, d / b, 0, a / b, 1, e / b }, 1)
            + G({ 0, d / b, 0, a / b, e / b, 1 }, 1)
            + G({ 0, d / b, 0, e / b, 1, a / b }, 1)
            + G({ 0, d / b, 0, e / b, a / b, 1 }, 1)
            + G({ 0, d / b, 1, 0, a / b, e / b }, 1)
            + G({ 0, d / b, 1, 0, e / b, a / b }, 1)
            + G({ 0, d / b, 1, a / b, 0, e / b }, 1)
            + G({ 0, d / b, 1, a / b, e / b, x / b }, 1)
            + G({ 0, d / b, 1, e / b, 0, a / b }, 1)
            + G({ 0, d / b, 1, e / b, a / b, x / b }, 1)
            + G({ 0, d / b, a / b, 0, 1, e / b }, 1)
            + G({ 0, d / b, a / b, 0, e / b, 1 }, 1)
            + G({ 0, d / b, a / b, 1, 0, e / b }, 1)
            + G({ 0, d / b, a / b, 1, e / b, x / b }, 1)
            + G({ 0, d / b, a / b, e / b, 1, x / b }, 1)
            + G({ 0, d / b, a / b, e / b, x / b, 1 }, 1)
            + G({ 0, d / b, e / b, 0, 1, a / b }, 1)
            + G({ 0, d / b, e / b, 0, a / b, 1 }, 1)
            + G({ 0, d / b, e / b, 1, 0, a / b }, 1)
            + G({ 0, d / b, e / b, 1, a / b, x / b }, 1)
            + G({ 0, d / b, e / b, a / b, 1, x / b }, 1)
            + G({ 0, d / b, e / b, a / b, x / b, 1 }, 1)
            + 2. * G({ a / b, 0, 0, 1, d / b, e / b }, 1)
            + 2. * G({ a / b, 0, 0, d / b, 1, e / b }, 1)
            + 2. * G({ a / b, 0, 0, d / b, e / b, 1 }, 1)
            + G({ a / b, 0, 1, 0, d / b, e / b }, 1)
            + G({ a / b, 0, 1, d / b, 0, e / b }, 1)
            + G({ a / b, 0, 1, d / b, e / b, x / b }, 1)
            + G({ a / b, 0, d / b, 0, 1, e / b }, 1)
            + G({ a / b, 0, d / b, 0, e / b, 1 }, 1)
            + G({ a / b, 0, d / b, 1, 0, e / b }, 1)
            + G({ a / b, 0, d / b, 1, e / b, x / b }, 1)
            + G({ a / b, 0, d / b, e / b, 1, x / b }, 1)
            + G({ a / b, 0, d / b, e / b, x / b, 1 }, 1)
            + (sy[11] + sy[12] + sy[13] + sy[14] + sy[15] + sy[16] + sy[17] + sy[6])
                * G({ 0, e }, { 1, se }, x)
            + G({ a / b, 1 }, 1) * (sy[18] - G({ a, 0, d, e }, { sa, 1, sd, se }, x))
            + G({ a / b }, 1)
                * (-G({ 0, b, 0, d, e }, { 1, sb, 1, sd, se }, x)
                    + G({ a, b, 0, d, e }, { sa, sb, 1, sd, se }, x))
            + G({ a, 0, b, 0, d, e }, { sa, 1, sb, 1, sd, se }, x)
            + (-sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27]
                  - sy[28] - sy[29] - sy[30] - sy[31] - sy[32] - sy[33] - sy[34]
                  + sy[19] * (sy[9] + sy[10] + sy[8]))
                * Log(b, sb)
            + (sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27]
                  + sy[28] + sy[29] + sy[30] + sy[31] + sy[32] + sy[33] + sy[34]
                  + sy[19] * (-sy[9] - sy[10] - sy[8]))
                * Log(-x, sb) };
        if (d != x) {
            res += (-sy[13] - sy[14] - sy[15] - sy[16] - sy[17]
                       + G({ 0, 1, a / b, x / b }, 1) + G({ 0, a / b, 1, x / b }, 1)
                       + G({ 0, a / b, x / b, 1 }, 1) + G({ a / b, 0, 1, x / b }, 1)
                       + G({ a / b, 0, x / b, 1 }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (-sy[24] - sy[25] - sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31]
                       - sy[32] - sy[33] - sy[34] + G({ 0, 1, a / b, d / b, x / b }, 1)
                       + G({ 0, 1, d / b, a / b, x / b }, 1)
                       + G({ 0, a / b, 1, d / b, x / b }, 1)
                       + G({ 0, a / b, d / b, 1, x / b }, 1)
                       + G({ 0, a / b, d / b, x / b, 1 }, 1)
                       + G({ 0, d / b, 1, a / b, x / b }, 1)
                       + G({ 0, d / b, a / b, 1, x / b }, 1)
                       + G({ 0, d / b, a / b, x / b, 1 }, 1)
                       + G({ a / b, 0, 1, d / b, x / b }, 1)
                       + G({ a / b, 0, d / b, 1, x / b }, 1)
                       + G({ a / b, 0, d / b, x / b, 1 }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    } else { // abcde
        const vector<complex<double>> sy = { G({ a, d, e }, { sa, sd, se }, x),
            G({ c / b, 0, a / b }, 1), G({ c / b, 0, d / b }, 1), G({ a / b, c / b }, 1),
            G({ 0, 0, d, e }, { 1, 1, sd, se }, x), G({ c / b, a / b }, 1),
            G({ a, e }, { sa, se }, x), G({ c / b, 0, d / b, a / b }, 1),
            G({ c / b, 0, d / b, e / b }, 1), G({ a / b, c / b, 0, d / b }, 1),
            G({ c / b, 0, a / b, d / b }, 1), G({ c / b, a / b, 0, d / b }, 1),
            G({ a, c, 0, d, e }, { sa, sc, 1, sd, se }, x), G({ a }, { sa }, x),
            G({ c / b, 0, d / b, e / b, a / b }, 1),
            G({ a / b, c / b, 0, d / b, e / b }, 1),
            G({ c / b, 0, a / b, d / b, e / b }, 1),
            G({ c / b, 0, d / b, a / b, e / b }, 1),
            G({ c / b, a / b, 0, d / b, e / b }, 1) };
        complex<double> res { sy[0] * sy[1] - sy[0] * sy[2] - sy[3] * sy[4]
            - sy[4] * sy[5] + sy[6] * sy[7] - sy[6] * sy[8]
            + sy[13]
                * (sy[14] + G({ 0, c / b, 0, d / b, e / b }, 1)
                    + 2. * G({ c / b, 0, 0, d / b, e / b }, 1)
                    + G({ c / b, 0, d / b, 0, e / b }, 1))
            - G({ 0, a / b, c / b, 0, d / b, e / b }, 1)
            - G({ 0, c / b, 0, a / b, d / b, e / b }, 1)
            - G({ 0, c / b, 0, d / b, a / b, e / b }, 1)
            - G({ 0, c / b, 0, d / b, e / b, a / b }, 1)
            - G({ 0, c / b, a / b, 0, d / b, e / b }, 1)
            - G({ a / b, 0, c / b, 0, d / b, e / b }, 1)
            - 2. * G({ a / b, c / b, 0, 0, d / b, e / b }, 1)
            - G({ a / b, c / b, 0, d / b, 0, e / b }, 1)
            - G({ a / b, c / b, 0, d / b, e / b, x / b }, 1)
            - 2. * G({ c / b, 0, 0, a / b, d / b, e / b }, 1)
            - 2. * G({ c / b, 0, 0, d / b, a / b, e / b }, 1)
            - 2. * G({ c / b, 0, 0, d / b, e / b, a / b }, 1)
            - 2. * G({ c / b, 0, a / b, 0, d / b, e / b }, 1)
            - G({ c / b, 0, a / b, d / b, 0, e / b }, 1)
            - G({ c / b, 0, a / b, d / b, e / b, x / b }, 1)
            - G({ c / b, 0, d / b, 0, a / b, e / b }, 1)
            - G({ c / b, 0, d / b, 0, e / b, a / b }, 1)
            - G({ c / b, 0, d / b, a / b, 0, e / b }, 1)
            - G({ c / b, 0, d / b, a / b, e / b, x / b }, 1)
            - G({ c / b, 0, d / b, e / b, 0, a / b }, 1)
            - G({ c / b, 0, d / b, e / b, a / b, x / b }, 1)
            - 2. * G({ c / b, a / b, 0, 0, d / b, e / b }, 1)
            - G({ c / b, a / b, 0, d / b, 0, e / b }, 1)
            - G({ c / b, a / b, 0, d / b, e / b, x / b }, 1)
            - sy[8] * G({ 0, a }, { 1, sa }, x)
            + (-sy[9] - sy[10] - sy[11] - sy[7]) * G({ 0, e }, { 1, se }, x)
            + (-sy[1] - G({ a / b, c / b, x / b }, 1) - G({ c / b, a / b, x / b }, 1))
                * G({ 0, d, e }, { 1, sd, se }, x)
            + sy[2] * G({ a, 0, e }, { sa, 1, se }, x)
            + sy[5] * G({ a, 0, d, e }, { sa, 1, sd, se }, x)
            + G({ a / b }, 1) * (sy[12] - G({ 0, c, 0, d, e }, { 1, sc, 1, sd, se }, x))
            + G({ c / b }, 1) * (-sy[12] + G({ a, 0, 0, d, e }, { sa, 1, 1, sd, se }, x))
            + G({ a, 0, c, 0, d, e }, { sa, 1, sc, 1, sd, se }, x)
            + (sy[14] + sy[15] + sy[16] + sy[17] + sy[18] - sy[13] * sy[8]) * Log(b, sb)
            + (-sy[14] - sy[15] - sy[16] - sy[17] - sy[18] + sy[13] * sy[8])
                * Log(-x, sb) };
        if (c != x) {
            res += (sy[3] - G({ a / b, x / b }, 1))
                * G({ c, 0, d, e }, { sc, 1, sd, se }, x);
        }
        if (d != x) {
            res += (sy[9] + sy[10] + sy[11] - G({ a / b, c / b, 0, x / b }, 1)
                       - G({ c / b, 0, a / b, x / b }, 1)
                       - G({ c / b, a / b, 0, x / b }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (sy[15] + sy[16] + sy[17] + sy[18]
                       - G({ a / b, c / b, 0, d / b, x / b }, 1)
                       - G({ c / b, 0, a / b, d / b, x / b }, 1)
                       - G({ c / b, 0, d / b, a / b, x / b }, 1)
                       - G({ c / b, a / b, 0, d / b, x / b }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_abc0de_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (a == b) {
        if (a == c) {
            // aaade
            const vector<complex<double>> sy = { G({ 0, 1, 1, d / a }, 1),
                G({ 0, 1, d / a, 1 }, 1), G({ 0, d / a, 1, 1 }, 1),
                G({ 0, 1, 1, d / a, e / a }, 1), G({ 0, 1, d / a, 1, e / a }, 1),
                G({ 0, 1, d / a, e / a, 1 }, 1), G({ 0, d / a, 1, 1, e / a }, 1),
                G({ 0, d / a, 1, e / a, 1 }, 1), G({ 0, d / a, e / a, 1, 1 }, 1) };
            complex<double> res { 2. * G({ 0, 0, 1, 1, d / a, e / a }, 1)
                + 2. * G({ 0, 0, 1, d / a, 1, e / a }, 1)
                + 2. * G({ 0, 0, 1, d / a, e / a, 1 }, 1)
                + 2. * G({ 0, 0, d / a, 1, 1, e / a }, 1)
                + 2. * G({ 0, 0, d / a, 1, e / a, 1 }, 1)
                + 2. * G({ 0, 0, d / a, e / a, 1, 1 }, 1)
                + G({ 0, 1, 0, 1, d / a, e / a }, 1) + G({ 0, 1, 0, d / a, 1, e / a }, 1)
                + G({ 0, 1, 0, d / a, e / a, 1 }, 1) + G({ 0, 1, 1, 0, d / a, e / a }, 1)
                + G({ 0, 1, 1, d / a, 0, e / a }, 1)
                + G({ 0, 1, 1, d / a, e / a, x / a }, 1)
                + G({ 0, 1, d / a, 0, 1, e / a }, 1) + G({ 0, 1, d / a, 0, e / a, 1 }, 1)
                + G({ 0, 1, d / a, 1, 0, e / a }, 1)
                + G({ 0, 1, d / a, 1, e / a, x / a }, 1)
                + G({ 0, 1, d / a, e / a, 1, x / a }, 1)
                + G({ 0, 1, d / a, e / a, x / a, 1 }, 1)
                + G({ 0, d / a, 0, 1, 1, e / a }, 1) + G({ 0, d / a, 0, 1, e / a, 1 }, 1)
                + G({ 0, d / a, 0, e / a, 1, 1 }, 1) + G({ 0, d / a, 1, 0, 1, e / a }, 1)
                + G({ 0, d / a, 1, 0, e / a, 1 }, 1) + G({ 0, d / a, 1, 1, 0, e / a }, 1)
                + G({ 0, d / a, 1, 1, e / a, x / a }, 1)
                + G({ 0, d / a, 1, e / a, 1, x / a }, 1)
                + G({ 0, d / a, 1, e / a, x / a, 1 }, 1)
                + G({ 0, d / a, e / a, 1, 1, x / a }, 1)
                + G({ 0, d / a, e / a, 1, x / a, 1 }, 1)
                + G({ 0, d / a, e / a, x / a, 1, 1 }, 1)
                + (sy[0] + sy[1] + sy[2]) * G({ 0, e }, { 1, se }, x)
                + G({ x / a, 1, 1 }, 1) * G({ 0, d, e }, { 1, sd, se }, x)
                - G({ x / a, 1 }, 1) * G({ a, 0, d, e }, { sa, 1, sd, se }, x)
                + G({ x / a }, 1) * G({ a, a, 0, d, e }, { sa, sa, 1, sd, se }, x)
                + G({ 0, a, a, 0, d, e }, { 1, sa, sa, 1, sd, se }, x)
                + (-sy[3] - sy[4] - sy[5] - sy[6] - sy[7] - sy[8]) * Log(a, sa)
                + (sy[3] + sy[4] + sy[5] + sy[6] + sy[7] + sy[8]) * Log(-x, sa) };
            if (d != x) {
                res += (-sy[0] - sy[1] - sy[2] + G({ 0, 1, 1, x / a }, 1)
                           + G({ 0, 1, x / a, 1 }, 1) + G({ 0, x / a, 1, 1 }, 1))
                    * G({ d, e }, { sd, se }, x);
            }
            if (e != x) {
                res += (-sy[3] - sy[4] - sy[5] - sy[6] - sy[7] - sy[8]
                           + G({ 0, 1, 1, d / a, x / a }, 1)
                           + G({ 0, 1, d / a, 1, x / a }, 1)
                           + G({ 0, 1, d / a, x / a, 1 }, 1)
                           + G({ 0, d / a, 1, 1, x / a }, 1)
                           + G({ 0, d / a, 1, x / a, 1 }, 1)
                           + G({ 0, d / a, x / a, 1, 1 }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        } else { // aacde
            const vector<complex<double>> sy = { G({ c / a, 1 }, 1),
                G({ c / a, 0, 1, d / a }, 1), G({ c / a, 0, d / a, 1 }, 1),
                G({ c / a, 1, 0, d / a }, 1), G({ c / a, 0, 1, d / a, e / a }, 1),
                G({ c / a, 0, d / a, 1, e / a }, 1), G({ c / a, 0, d / a, e / a, 1 }, 1),
                G({ c / a, 1, 0, d / a, e / a }, 1) };
            complex<double> res { -G({ 0, c / a, 0, 1, d / a, e / a }, 1)
                - G({ 0, c / a, 0, d / a, 1, e / a }, 1)
                - G({ 0, c / a, 0, d / a, e / a, 1 }, 1)
                - G({ 0, c / a, 1, 0, d / a, e / a }, 1)
                - 2. * G({ c / a, 0, 0, 1, d / a, e / a }, 1)
                - 2. * G({ c / a, 0, 0, d / a, 1, e / a }, 1)
                - 2. * G({ c / a, 0, 0, d / a, e / a, 1 }, 1)
                - 2. * G({ c / a, 0, 1, 0, d / a, e / a }, 1)
                - G({ c / a, 0, 1, d / a, 0, e / a }, 1)
                - G({ c / a, 0, 1, d / a, e / a, x / a }, 1)
                - G({ c / a, 0, d / a, 0, 1, e / a }, 1)
                - G({ c / a, 0, d / a, 0, e / a, 1 }, 1)
                - G({ c / a, 0, d / a, 1, 0, e / a }, 1)
                - G({ c / a, 0, d / a, 1, e / a, x / a }, 1)
                - G({ c / a, 0, d / a, e / a, 1, x / a }, 1)
                - G({ c / a, 0, d / a, e / a, x / a, 1 }, 1)
                - 2. * G({ c / a, 1, 0, 0, d / a, e / a }, 1)
                - G({ c / a, 1, 0, d / a, 0, e / a }, 1)
                - G({ c / a, 1, 0, d / a, e / a, x / a }, 1)
                + (-sy[1] - sy[2] - sy[3]) * G({ 0, e }, { 1, se }, x)
                + (-G({ c / a, 1, x / a }, 1) - G({ c / a, x / a, 1 }, 1))
                    * G({ 0, d, e }, { 1, sd, se }, x)
                - sy[0] * G({ 0, 0, d, e }, { 1, 1, sd, se }, x)
                + G({ x / a }, 1) * G({ a, c, 0, d, e }, { sa, sc, 1, sd, se }, x)
                + G({ 0, a, c, 0, d, e }, { 1, sa, sc, 1, sd, se }, x)
                + (sy[4] + sy[5] + sy[6] + sy[7]) * Log(a, sa)
                + (-sy[4] - sy[5] - sy[6] - sy[7]) * Log(-x, sa) };
            if (c != x) {
                res += (sy[0] - G({ x / a, 1 }, 1))
                    * G({ c, 0, d, e }, { sc, 1, sd, se }, x);
            }
            if (d != x) {
                res += (sy[1] + sy[2] + sy[3] - G({ c / a, 0, 1, x / a }, 1)
                           - G({ c / a, 0, x / a, 1 }, 1) - G({ c / a, 1, 0, x / a }, 1))
                    * G({ d, e }, { sd, se }, x);
            }
            if (e != x) {
                res += (sy[4] + sy[5] + sy[6] + sy[7]
                           - G({ c / a, 0, 1, d / a, x / a }, 1)
                           - G({ c / a, 0, d / a, 1, x / a }, 1)
                           - G({ c / a, 0, d / a, x / a, 1 }, 1)
                           - G({ c / a, 1, 0, d / a, x / a }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else { // abcde
        const vector<complex<double>> sy
            = { G({ b / a, c / a }, 1), G({ b / a, c / a, 0, d / a }, 1), G({ b / a }, 1),
                  G({ b / a, c / a, 0, d / a, e / a }, 1) };
        complex<double> res { G({ 0, b / a, c / a, 0, d / a, e / a }, 1)
            + G({ b / a, 0, c / a, 0, d / a, e / a }, 1)
            + 2. * G({ b / a, c / a, 0, 0, d / a, e / a }, 1)
            + G({ b / a, c / a, 0, d / a, 0, e / a }, 1)
            + G({ b / a, c / a, 0, d / a, e / a, x / a }, 1)
            + sy[1] * G({ 0, e }, { 1, se }, x)
            + G({ b / a, c / a, x / a }, 1) * G({ 0, d, e }, { 1, sd, se }, x)
            + sy[0] * G({ 0, 0, d, e }, { 1, 1, sd, se }, x)
            + sy[2] * G({ 0, c, 0, d, e }, { 1, sc, 1, sd, se }, x)
            + G({ 0, b, c, 0, d, e }, { 1, sb, sc, 1, sd, se }, x) - sy[3] * Log(a, sa)
            + sy[3] * Log(-x, sa) };
        if (b != x) {
            res += (-sy[2] + G({ x / a }, 1))
                * G({ b, c, 0, d, e }, { sb, sc, 1, sd, se }, x);
        }
        if (c != x) {
            res += (-sy[0] + G({ b / a, x / a }, 1))
                * G({ c, 0, d, e }, { sc, 1, sd, se }, x);
        }
        if (d != x) {
            res += (-sy[1] + G({ b / a, c / a, 0, x / a }, 1))
                * G({ d, e }, { sd, se }, x);
        }
        if (e != x) {
            res += (-sy[3] + G({ b / a, c / a, 0, d / a, x / a }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}

complex<double> G6_abcd0e_e(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    const vector<complex<double>> sy
        = { G({ a, b, c }, { sa, sb, sc }, x), G({ 0, d / e, c / e }, 1),
              G({ a, b }, { sa, sb }, x), G({ 0, d / e, c / e, b / e }, 1),
              G({ 0, d / e }, 1), G({ a, b, c, d }, { sa, sb, sc, sd }, x),
              G({ x / e }, 1), G({ 0, a, b, c, d }, { 1, sa, sb, sc, sd }, x),
              G({ a }, { sa }, x), G({ 0, d / e, c / e, b / e, a / e }, 1),
              G({ a, 0, b, c, d }, { sa, 1, sb, sc, sd }, x),
              G({ a, b, 0, c, d }, { sa, sb, 1, sc, sd }, x),
              G({ a, b, c, 0, d }, { sa, sb, sc, 1, sd }, x) };
    complex<double> res { (sy[10] + sy[11] + sy[12]) * sy[6] + sy[6] * sy[7]
        + 2. * sy[0] * G({ 0, 0, d / e }, 1)
        + sy[2]
            * (-sy[3] - 2. * G({ 0, 0, d / e, c / e }, 1) - G({ 0, d / e, 0, c / e }, 1))
        + sy[8]
            * (sy[9] + 2. * G({ 0, 0, d / e, c / e, b / e }, 1)
                + G({ 0, d / e, 0, c / e, b / e }, 1)
                + G({ 0, d / e, c / e, 0, b / e }, 1))
        - 2. * G({ 0, 0, d / e, c / e, b / e, a / e }, 1)
        - G({ 0, d / e, 0, c / e, b / e, a / e }, 1)
        - G({ 0, d / e, c / e, 0, b / e, a / e }, 1)
        - G({ 0, d / e, c / e, b / e, 0, a / e }, 1)
        - G({ 0, d / e, c / e, b / e, a / e, x / e }, 1)
        + (-sy[10] - sy[11] - sy[12] - sy[7]) * G({ e }, { se }, x)
        - sy[3] * G({ 0, a }, { 1, sa }, x)
        + sy[5] * (G({ 0, x / e }, 1) + G({ 0, e }, { 1, se }, x))
        + sy[1]
            * (sy[0] + G({ 0, a, b }, { 1, sa, sb }, x)
                + G({ a, 0, b }, { sa, 1, sb }, x))
        + sy[4]
            * (-sy[5] - G({ 0, a, b, c }, { 1, sa, sb, sc }, x)
                - G({ a, 0, b, c }, { sa, 1, sb, sc }, x)
                - G({ a, b, 0, c }, { sa, sb, 1, sc }, x))
        + G({ 0, 0, a, b, c, d }, { 1, 1, sa, sb, sc, sd }, x)
        + G({ 0, a, 0, b, c, d }, { 1, sa, 1, sb, sc, sd }, x)
        + G({ 0, a, b, 0, c, d }, { 1, sa, sb, 1, sc, sd }, x)
        + G({ 0, a, b, c, 0, d }, { 1, sa, sb, sc, 1, sd }, x)
        + G({ a, 0, 0, b, c, d }, { sa, 1, 1, sb, sc, sd }, x)
        + G({ a, 0, b, 0, c, d }, { sa, 1, sb, 1, sc, sd }, x)
        + G({ a, 0, b, c, 0, d }, { sa, 1, sb, sc, 1, sd }, x)
        + G({ a, b, 0, 0, c, d }, { sa, sb, 1, 1, sc, sd }, x)
        + G({ a, b, 0, c, 0, d }, { sa, sb, 1, sc, 1, sd }, x)
        + G({ a, b, c, 0, 0, d }, { sa, sb, sc, 1, 1, sd }, x)
        + (sy[9] + sy[1] * sy[2] - sy[0] * sy[4] - sy[3] * sy[8]) * Log(e, se)
        + (-sy[9] - sy[1] * sy[2] + sy[0] * sy[4] + sy[3] * sy[8]) * Log(-x, se) };
    return res;
}
complex<double> G6_abcd0e_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;

    const vector<complex<double>> sy = { G({ a, b, c }, { sa, sb, sc }, x),
        G({ 0, e / d, c / d }, 1), G({ a, 0, b }, { sa, 1, sb }, x),
        G({ a, 0, e }, { sa, 1, se }, x), G({ a, b, e }, { sa, sb, se }, x),
        G({ 0, c / d, e / d }, 1), G({ 0, a, b }, { 1, sa, sb }, x),
        G({ c / d, 0, e / d }, 1), G({ c / d, 0, b / d }, 1),
        G({ c / d, b / d, a / d }, 1), G({ a, e }, { sa, se }, x),
        G({ 0, c / d, b / d, a / d }, 1), G({ 0, a }, { 1, sa }, x),
        G({ 0, c / d, b / d, e / d }, 1), G({ a, b, 0, e }, { sa, sb, 1, se }, x),
        G({ 0, e / d }, 1), G({ a, b, c, e }, { sa, sb, sc, se }, x),
        G({ a, b }, { sa, sb }, x), G({ 0, c / d, e / d, b / d }, 1),
        G({ 0, e / d, c / d, b / d }, 1), G({ c / d, 0, e / d, b / d }, 1),
        G({ c / d, 0, b / d, e / d }, 1), G({ c / d, b / d, 0, e / d }, 1),
        G({ c / d, 0, b / d, a / d }, 1), G({ c / d, b / d, 0, a / d }, 1),
        G({ a }, { sa }, x), G({ 0, c / d, b / d, e / d, a / d }, 1),
        G({ 0, c / d, e / d, b / d, a / d }, 1), G({ 0, e / d, c / d, b / d, a / d }, 1),
        G({ c / d, 0, b / d, e / d, a / d }, 1), G({ c / d, 0, e / d, b / d, a / d }, 1),
        G({ c / d, b / d, 0, e / d, a / d }, 1), G({ 0, c / d, b / d, a / d, e / d }, 1),
        G({ c / d, 0, b / d, a / d, e / d }, 1), G({ c / d, b / d, 0, a / d, e / d }, 1),
        G({ c / d, b / d, a / d, 0, e / d }, 1) };
    complex<double> res { sy[10] * sy[11] - sy[12] * sy[13]
        + sy[12] * (-sy[18] - sy[19] - sy[20] - sy[21] - sy[22])
        + sy[10] * (-sy[13] - sy[21] - sy[22] + sy[23] + sy[24]) + sy[1] * (sy[0] + sy[2])
        + sy[5] * (sy[2] + sy[4] + sy[6]) + sy[2] * sy[7] + sy[6] * (sy[1] + sy[7])
        + sy[4] * (sy[7] - sy[8]) + sy[3] * (sy[9] + sy[8])
        + (-sy[14] + sy[16]) * G({ 0, c / d }, 1) + 2. * sy[0] * G({ 0, 0, e / d }, 1)
        + (sy[3] - sy[4]) * G({ 0, c / d, b / d }, 1)
        + sy[17]
            * (-sy[18] - sy[19] - sy[20] - 2. * G({ 0, 0, c / d, e / d }, 1)
                - 2. * G({ 0, 0, e / d, c / d }, 1) - 2. * G({ 0, c / d, 0, e / d }, 1)
                - G({ 0, e / d, 0, c / d }, 1) - 2. * G({ c / d, 0, 0, e / d }, 1))
        + sy[25]
            * (sy[26] + sy[27] + sy[28] + sy[29] + sy[30] + sy[31]
                + 2. * G({ 0, 0, c / d, b / d, e / d }, 1)
                + 2. * G({ 0, 0, c / d, e / d, b / d }, 1)
                + 2. * G({ 0, 0, e / d, c / d, b / d }, 1)
                + 2. * G({ 0, c / d, 0, b / d, e / d }, 1)
                + 2. * G({ 0, c / d, 0, e / d, b / d }, 1)
                + 2. * G({ 0, c / d, b / d, 0, e / d }, 1)
                + G({ 0, c / d, e / d, 0, b / d }, 1)
                + G({ 0, e / d, 0, c / d, b / d }, 1)
                + G({ 0, e / d, c / d, 0, b / d }, 1)
                + 2. * G({ c / d, 0, 0, b / d, e / d }, 1)
                + 2. * G({ c / d, 0, 0, e / d, b / d }, 1)
                + 2. * G({ c / d, 0, b / d, 0, e / d }, 1)
                + G({ c / d, 0, e / d, 0, b / d }, 1)
                + 2. * G({ c / d, b / d, 0, 0, e / d }, 1))
        - 2. * G({ 0, 0, c / d, b / d, a / d, e / d }, 1)
        - 2. * G({ 0, 0, c / d, b / d, e / d, a / d }, 1)
        - 2. * G({ 0, 0, c / d, e / d, b / d, a / d }, 1)
        - 2. * G({ 0, 0, e / d, c / d, b / d, a / d }, 1)
        - 2. * G({ 0, c / d, 0, b / d, a / d, e / d }, 1)
        - 2. * G({ 0, c / d, 0, b / d, e / d, a / d }, 1)
        - 2. * G({ 0, c / d, 0, e / d, b / d, a / d }, 1)
        - 2. * G({ 0, c / d, b / d, 0, a / d, e / d }, 1)
        - 2. * G({ 0, c / d, b / d, 0, e / d, a / d }, 1)
        - 2. * G({ 0, c / d, b / d, a / d, 0, e / d }, 1)
        - G({ 0, c / d, b / d, a / d, e / d, x / d }, 1)
        - G({ 0, c / d, b / d, e / d, 0, a / d }, 1)
        - G({ 0, c / d, b / d, e / d, a / d, x / d }, 1)
        - G({ 0, c / d, e / d, 0, b / d, a / d }, 1)
        - G({ 0, c / d, e / d, b / d, 0, a / d }, 1)
        - G({ 0, c / d, e / d, b / d, a / d, x / d }, 1)
        - G({ 0, e / d, 0, c / d, b / d, a / d }, 1)
        - G({ 0, e / d, c / d, 0, b / d, a / d }, 1)
        - G({ 0, e / d, c / d, b / d, 0, a / d }, 1)
        - G({ 0, e / d, c / d, b / d, a / d, x / d }, 1)
        - 2. * G({ c / d, 0, 0, b / d, a / d, e / d }, 1)
        - 2. * G({ c / d, 0, 0, b / d, e / d, a / d }, 1)
        - 2. * G({ c / d, 0, 0, e / d, b / d, a / d }, 1)
        - 2. * G({ c / d, 0, b / d, 0, a / d, e / d }, 1)
        - 2. * G({ c / d, 0, b / d, 0, e / d, a / d }, 1)
        - 2. * G({ c / d, 0, b / d, a / d, 0, e / d }, 1)
        - G({ c / d, 0, b / d, a / d, e / d, x / d }, 1)
        - G({ c / d, 0, b / d, e / d, 0, a / d }, 1)
        - G({ c / d, 0, b / d, e / d, a / d, x / d }, 1)
        - G({ c / d, 0, e / d, 0, b / d, a / d }, 1)
        - G({ c / d, 0, e / d, b / d, 0, a / d }, 1)
        - G({ c / d, 0, e / d, b / d, a / d, x / d }, 1)
        - 2. * G({ c / d, b / d, 0, 0, a / d, e / d }, 1)
        - 2. * G({ c / d, b / d, 0, 0, e / d, a / d }, 1)
        - 2. * G({ c / d, b / d, 0, a / d, 0, e / d }, 1)
        - G({ c / d, b / d, 0, a / d, e / d, x / d }, 1)
        - G({ c / d, b / d, 0, e / d, 0, a / d }, 1)
        - G({ c / d, b / d, 0, e / d, a / d, x / d }, 1)
        - 2. * G({ c / d, b / d, a / d, 0, 0, e / d }, 1)
        - G({ c / d, b / d, a / d, 0, e / d, x / d }, 1)
        + (-sy[11] - sy[23] - sy[24] - G({ c / d, b / d, a / d, x / d }, 1))
            * G({ 0, e }, { 1, se }, x)
        - sy[9] * G({ 0, 0, e }, { 1, 1, se }, x)
        + G({ c / d, b / d }, 1) * (-sy[14] + G({ a, 0, 0, e }, { sa, 1, 1, se }, x))
        + sy[15]
            * (-sy[16] - G({ 0, a, b, c }, { 1, sa, sb, sc }, x)
                - G({ a, 0, b, c }, { sa, 1, sb, sc }, x)
                - G({ a, b, 0, c }, { sa, sb, 1, sc }, x))
        + G({ c / d }, 1)
            * (-G({ a, b, 0, 0, e }, { sa, sb, 1, 1, se }, x)
                + G({ a, b, c, 0, e }, { sa, sb, sc, 1, se }, x))
        + G({ a, b, c, 0, 0, e }, { sa, sb, sc, 1, 1, se }, x)
        + (-(sy[0] * sy[15])
              + (-sy[13] - sy[18] - sy[19] - sy[20] - sy[21] - sy[22]) * sy[25] + sy[26]
              + sy[27] + sy[28] + sy[29] + sy[30] + sy[31] + sy[32] + sy[33] + sy[34]
              + sy[35] + sy[17] * (sy[1] + sy[5] + sy[7]))
            * Log(d, sd)
        + (sy[0] * sy[15] + (sy[13] + sy[18] + sy[19] + sy[20] + sy[21] + sy[22]) * sy[25]
              - sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31] - sy[32] - sy[33]
              - sy[34] - sy[35] + sy[17] * (-sy[1] - sy[5] - sy[7]))
            * Log(-x, sd) };
    if (e != x) {
        res += (sy[32] + sy[33] + sy[34] + sy[35]
                   - G({ 0, c / d, b / d, a / d, x / d }, 1)
                   - G({ c / d, 0, b / d, a / d, x / d }, 1)
                   - G({ c / d, b / d, 0, a / d, x / d }, 1)
                   - G({ c / d, b / d, a / d, 0, x / d }, 1))
            * G({ e }, { se }, x);
    }
    return res;
}
complex<double> G6_abcd0e_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (c == d) {
        const vector<complex<double>> sy = { G({ 0, a, b }, { 1, sa, sb }, x),
            G({ 0, e / c, 1 }, 1), G({ a, 0, e }, { sa, 1, se }, x),
            G({ a, b, e }, { sa, sb, se }, x), G({ 0, b / c, 1 }, 1),
            G({ a, 0, b }, { sa, 1, sb }, x), G({ 0, 1, e / c }, 1),
            G({ b / c, a / c, 1 }, 1), G({ c, 0, e }, { sc, 1, se }, x),
            G({ a, e }, { sa, se }, x), G({ 0, 1, b / c, a / c }, 1),
            G({ 0, a }, { 1, sa }, x), G({ 0, 1, b / c, e / c }, 1),
            G({ a, b }, { sa, sb }, x), G({ 0, 1, e / c, b / c }, 1),
            G({ 0, e / c, 1, b / c }, 1), G({ 0, e / c, b / c, 1 }, 1),
            G({ a, c, 0, e }, { sa, sc, 1, se }, x), G({ 0, b / c, 1, e / c }, 1),
            G({ 0, b / c, e / c, 1 }, 1), G({ b / c, 0, 1, e / c }, 1),
            G({ b / c, 0, e / c, 1 }, 1), G({ 0, b / c, 1, a / c }, 1),
            G({ 0, b / c, a / c, 1 }, 1), G({ b / c, 0, 1, a / c }, 1),
            G({ b / c, 0, a / c, 1 }, 1), G({ a }, { sa }, x),
            G({ 0, 1, b / c, e / c, a / c }, 1), G({ 0, 1, e / c, b / c, a / c }, 1),
            G({ 0, b / c, 1, e / c, a / c }, 1), G({ 0, b / c, e / c, 1, a / c }, 1),
            G({ 0, b / c, e / c, a / c, 1 }, 1), G({ 0, e / c, 1, b / c, a / c }, 1),
            G({ 0, e / c, b / c, 1, a / c }, 1), G({ 0, e / c, b / c, a / c, 1 }, 1),
            G({ b / c, 0, 1, e / c, a / c }, 1), G({ b / c, 0, e / c, 1, a / c }, 1),
            G({ b / c, 0, e / c, a / c, 1 }, 1), G({ 0, 1, b / c, a / c, e / c }, 1),
            G({ 0, b / c, 1, a / c, e / c }, 1), G({ 0, b / c, a / c, 1, e / c }, 1),
            G({ 0, b / c, a / c, e / c, 1 }, 1), G({ b / c, 0, 1, a / c, e / c }, 1),
            G({ b / c, 0, a / c, 1, e / c }, 1), G({ b / c, 0, a / c, e / c, 1 }, 1),
            G({ b / c, a / c, 0, 1, e / c }, 1), G({ b / c, a / c, 0, e / c, 1 }, 1) };
        complex<double> res { sy[9] * sy[10] - sy[11] * sy[12] + sy[0] * sy[1]
            + sy[11] * (-sy[14] - sy[15] - sy[16] - sy[18] - sy[19] - sy[20] - sy[21])
            + sy[9]
                * (-sy[12] - sy[18] - sy[19] - sy[20] - sy[21] + sy[22] + sy[23] + sy[24]
                    + sy[25])
            - sy[3] * sy[4] + sy[1] * (sy[3] + sy[5]) + (sy[0] + sy[3] + sy[5]) * sy[6]
            + sy[2] * (sy[4] + sy[7]) - sy[7] * sy[8]
            + (sy[2] - sy[3]) * G({ 0, 1, b / c }, 1)
            + sy[8] * G({ b / c, a / c, x / c }, 1)
            + sy[13]
                * (-sy[14] - sy[15] - sy[16] - 2. * G({ 0, 0, 1, e / c }, 1)
                    - 2. * G({ 0, 0, e / c, 1 }, 1) - G({ 0, 1, 0, e / c }, 1))
            + sy[26]
                * (sy[27] + sy[28] + sy[29] + sy[30] + sy[31] + sy[32] + sy[33] + sy[34]
                    + sy[35] + sy[36] + sy[37] + 2. * G({ 0, 0, 1, b / c, e / c }, 1)
                    + 2. * G({ 0, 0, 1, e / c, b / c }, 1)
                    + 2. * G({ 0, 0, b / c, 1, e / c }, 1)
                    + 2. * G({ 0, 0, b / c, e / c, 1 }, 1)
                    + 2. * G({ 0, 0, e / c, 1, b / c }, 1)
                    + 2. * G({ 0, 0, e / c, b / c, 1 }, 1)
                    + G({ 0, 1, 0, b / c, e / c }, 1) + G({ 0, 1, 0, e / c, b / c }, 1)
                    + G({ 0, 1, b / c, 0, e / c }, 1) + G({ 0, 1, e / c, 0, b / c }, 1)
                    + 2. * G({ 0, b / c, 0, 1, e / c }, 1)
                    + 2. * G({ 0, b / c, 0, e / c, 1 }, 1)
                    + G({ 0, b / c, 1, 0, e / c }, 1) + G({ 0, e / c, 0, 1, b / c }, 1)
                    + G({ 0, e / c, 0, b / c, 1 }, 1) + G({ 0, e / c, 1, 0, b / c }, 1)
                    + 2. * G({ b / c, 0, 0, 1, e / c }, 1)
                    + 2. * G({ b / c, 0, 0, e / c, 1 }, 1)
                    + G({ b / c, 0, 1, 0, e / c }, 1))
            - 2. * G({ 0, 0, 1, b / c, a / c, e / c }, 1)
            - 2. * G({ 0, 0, 1, b / c, e / c, a / c }, 1)
            - 2. * G({ 0, 0, 1, e / c, b / c, a / c }, 1)
            - 2. * G({ 0, 0, b / c, 1, a / c, e / c }, 1)
            - 2. * G({ 0, 0, b / c, 1, e / c, a / c }, 1)
            - 2. * G({ 0, 0, b / c, a / c, 1, e / c }, 1)
            - 2. * G({ 0, 0, b / c, a / c, e / c, 1 }, 1)
            - 2. * G({ 0, 0, b / c, e / c, 1, a / c }, 1)
            - 2. * G({ 0, 0, b / c, e / c, a / c, 1 }, 1)
            - 2. * G({ 0, 0, e / c, 1, b / c, a / c }, 1)
            - 2. * G({ 0, 0, e / c, b / c, 1, a / c }, 1)
            - 2. * G({ 0, 0, e / c, b / c, a / c, 1 }, 1)
            - G({ 0, 1, 0, b / c, a / c, e / c }, 1)
            - G({ 0, 1, 0, b / c, e / c, a / c }, 1)
            - G({ 0, 1, 0, e / c, b / c, a / c }, 1)
            - G({ 0, 1, b / c, 0, a / c, e / c }, 1)
            - G({ 0, 1, b / c, 0, e / c, a / c }, 1)
            - G({ 0, 1, b / c, a / c, 0, e / c }, 1)
            - G({ 0, 1, b / c, a / c, e / c, x / c }, 1)
            - G({ 0, 1, b / c, e / c, 0, a / c }, 1)
            - G({ 0, 1, b / c, e / c, a / c, x / c }, 1)
            - G({ 0, 1, e / c, 0, b / c, a / c }, 1)
            - G({ 0, 1, e / c, b / c, 0, a / c }, 1)
            - G({ 0, 1, e / c, b / c, a / c, x / c }, 1)
            - 2. * G({ 0, b / c, 0, 1, a / c, e / c }, 1)
            - 2. * G({ 0, b / c, 0, 1, e / c, a / c }, 1)
            - 2. * G({ 0, b / c, 0, a / c, 1, e / c }, 1)
            - 2. * G({ 0, b / c, 0, a / c, e / c, 1 }, 1)
            - 2. * G({ 0, b / c, 0, e / c, 1, a / c }, 1)
            - 2. * G({ 0, b / c, 0, e / c, a / c, 1 }, 1)
            - G({ 0, b / c, 1, 0, a / c, e / c }, 1)
            - G({ 0, b / c, 1, 0, e / c, a / c }, 1)
            - G({ 0, b / c, 1, a / c, 0, e / c }, 1)
            - G({ 0, b / c, 1, a / c, e / c, x / c }, 1)
            - G({ 0, b / c, 1, e / c, 0, a / c }, 1)
            - G({ 0, b / c, 1, e / c, a / c, x / c }, 1)
            - 2. * G({ 0, b / c, a / c, 0, 1, e / c }, 1)
            - 2. * G({ 0, b / c, a / c, 0, e / c, 1 }, 1)
            - G({ 0, b / c, a / c, 1, 0, e / c }, 1)
            - G({ 0, b / c, a / c, 1, e / c, x / c }, 1)
            - G({ 0, b / c, a / c, e / c, 1, x / c }, 1)
            - G({ 0, b / c, a / c, e / c, x / c, 1 }, 1)
            - G({ 0, b / c, e / c, 0, 1, a / c }, 1)
            - G({ 0, b / c, e / c, 0, a / c, 1 }, 1)
            - G({ 0, b / c, e / c, 1, 0, a / c }, 1)
            - G({ 0, b / c, e / c, 1, a / c, x / c }, 1)
            - G({ 0, b / c, e / c, a / c, 1, x / c }, 1)
            - G({ 0, b / c, e / c, a / c, x / c, 1 }, 1)
            - G({ 0, e / c, 0, 1, b / c, a / c }, 1)
            - G({ 0, e / c, 0, b / c, 1, a / c }, 1)
            - G({ 0, e / c, 0, b / c, a / c, 1 }, 1)
            - G({ 0, e / c, 1, 0, b / c, a / c }, 1)
            - G({ 0, e / c, 1, b / c, 0, a / c }, 1)
            - G({ 0, e / c, 1, b / c, a / c, x / c }, 1)
            - G({ 0, e / c, b / c, 0, 1, a / c }, 1)
            - G({ 0, e / c, b / c, 0, a / c, 1 }, 1)
            - G({ 0, e / c, b / c, 1, 0, a / c }, 1)
            - G({ 0, e / c, b / c, 1, a / c, x / c }, 1)
            - G({ 0, e / c, b / c, a / c, 1, x / c }, 1)
            - G({ 0, e / c, b / c, a / c, x / c, 1 }, 1)
            - 2. * G({ b / c, 0, 0, 1, a / c, e / c }, 1)
            - 2. * G({ b / c, 0, 0, 1, e / c, a / c }, 1)
            - 2. * G({ b / c, 0, 0, a / c, 1, e / c }, 1)
            - 2. * G({ b / c, 0, 0, a / c, e / c, 1 }, 1)
            - 2. * G({ b / c, 0, 0, e / c, 1, a / c }, 1)
            - 2. * G({ b / c, 0, 0, e / c, a / c, 1 }, 1)
            - G({ b / c, 0, 1, 0, a / c, e / c }, 1)
            - G({ b / c, 0, 1, 0, e / c, a / c }, 1)
            - G({ b / c, 0, 1, a / c, 0, e / c }, 1)
            - G({ b / c, 0, 1, a / c, e / c, x / c }, 1)
            - G({ b / c, 0, 1, e / c, 0, a / c }, 1)
            - G({ b / c, 0, 1, e / c, a / c, x / c }, 1)
            - 2. * G({ b / c, 0, a / c, 0, 1, e / c }, 1)
            - 2. * G({ b / c, 0, a / c, 0, e / c, 1 }, 1)
            - G({ b / c, 0, a / c, 1, 0, e / c }, 1)
            - G({ b / c, 0, a / c, 1, e / c, x / c }, 1)
            - G({ b / c, 0, a / c, e / c, 1, x / c }, 1)
            - G({ b / c, 0, a / c, e / c, x / c, 1 }, 1)
            - G({ b / c, 0, e / c, 0, 1, a / c }, 1)
            - G({ b / c, 0, e / c, 0, a / c, 1 }, 1)
            - G({ b / c, 0, e / c, 1, 0, a / c }, 1)
            - G({ b / c, 0, e / c, 1, a / c, x / c }, 1)
            - G({ b / c, 0, e / c, a / c, 1, x / c }, 1)
            - G({ b / c, 0, e / c, a / c, x / c, 1 }, 1)
            - 2. * G({ b / c, a / c, 0, 0, 1, e / c }, 1)
            - 2. * G({ b / c, a / c, 0, 0, e / c, 1 }, 1)
            - G({ b / c, a / c, 0, 1, 0, e / c }, 1)
            - G({ b / c, a / c, 0, 1, e / c, x / c }, 1)
            - G({ b / c, a / c, 0, e / c, 1, x / c }, 1)
            - G({ b / c, a / c, 0, e / c, x / c, 1 }, 1)
            + (-sy[10] - sy[22] - sy[23] - sy[24] - sy[25]
                  - G({ b / c, a / c, x / c, 1 }, 1))
                * G({ 0, e }, { 1, se }, x)
            + G({ b / c, a / c }, 1) * (-sy[17] + G({ 0, c, 0, e }, { 1, sc, 1, se }, x))
            + G({ b / c, 1 }, 1) * (sy[17] - G({ a, b, 0, e }, { sa, sb, 1, se }, x))
            + G({ b / c }, 1)
                * (-G({ a, 0, c, 0, e }, { sa, 1, sc, 1, se }, x)
                    + G({ a, b, c, 0, e }, { sa, sb, sc, 1, se }, x))
            + G({ a, b, 0, c, 0, e }, { sa, sb, 1, sc, 1, se }, x)
            + ((-sy[12] - sy[14] - sy[15] - sy[16] - sy[18] - sy[19] - sy[20] - sy[21])
                      * sy[26]
                  + sy[27] + sy[28] + sy[29] + sy[30] + sy[31] + sy[32] + sy[33] + sy[34]
                  + sy[35] + sy[36] + sy[37] + sy[38] + sy[39] + sy[40] + sy[41] + sy[42]
                  + sy[43] + sy[44] + sy[45] + sy[46] + sy[13] * (sy[1] + sy[6]))
                * Log(c, sc)
            + ((sy[12] + sy[14] + sy[15] + sy[16] + sy[18] + sy[19] + sy[20] + sy[21])
                      * sy[26]
                  - sy[27] - sy[28] - sy[29] - sy[30] - sy[31] - sy[32] - sy[33] - sy[34]
                  - sy[35] - sy[36] - sy[37] - sy[38] - sy[39] - sy[40] - sy[41] - sy[42]
                  - sy[43] - sy[44] - sy[45] - sy[46] + sy[13] * (-sy[1] - sy[6]))
                * Log(-x, sc) };
        if (e != x) {
            res += (sy[38] + sy[39] + sy[40] + sy[41] + sy[42] + sy[43] + sy[44] + sy[45]
                       + sy[46] - G({ 0, 1, b / c, a / c, x / c }, 1)
                       - G({ 0, b / c, 1, a / c, x / c }, 1)
                       - G({ 0, b / c, a / c, 1, x / c }, 1)
                       - G({ 0, b / c, a / c, x / c, 1 }, 1)
                       - G({ b / c, 0, 1, a / c, x / c }, 1)
                       - G({ b / c, 0, a / c, 1, x / c }, 1)
                       - G({ b / c, 0, a / c, x / c, 1 }, 1)
                       - G({ b / c, a / c, 0, 1, x / c }, 1)
                       - G({ b / c, a / c, 0, x / c, 1 }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    } else {
        const vector<complex<double>> sy = { G({ a, 0, e }, { sa, 1, se }, x),
            G({ b / c, d / c, a / c }, 1), G({ a, b, e }, { sa, sb, se }, x),
            G({ d / c, 0, b / c }, 1), G({ d / c, 0, e / c }, 1),
            G({ d / c, b / c, a / c }, 1), G({ b / c, a / c, d / c }, 1),
            G({ d / c, b / c }, 1), G({ a, 0, 0, e }, { sa, 1, 1, se }, x),
            G({ a, d, 0, e }, { sa, sd, 1, se }, x), G({ a, e }, { sa, se }, x),
            G({ b / c, d / c, 0, a / c }, 1), G({ 0, a }, { 1, sa }, x),
            G({ b / c, d / c, 0, e / c }, 1), G({ a, b }, { sa, sb }, x),
            G({ d / c, 0, e / c, b / c }, 1), G({ d / c, 0, b / c, e / c }, 1),
            G({ d / c, b / c, 0, e / c }, 1), G({ d / c, 0, b / c, a / c }, 1),
            G({ d / c, b / c, 0, a / c }, 1),
            G({ a, b, d, 0, e }, { sa, sb, sd, 1, se }, x), G({ a }, { sa }, x),
            G({ b / c, d / c, 0, e / c, a / c }, 1),
            G({ d / c, 0, b / c, e / c, a / c }, 1),
            G({ d / c, 0, e / c, b / c, a / c }, 1),
            G({ d / c, b / c, 0, e / c, a / c }, 1),
            G({ b / c, a / c, d / c, 0, e / c }, 1),
            G({ b / c, d / c, 0, a / c, e / c }, 1),
            G({ b / c, d / c, a / c, 0, e / c }, 1),
            G({ d / c, 0, b / c, a / c, e / c }, 1),
            G({ d / c, b / c, 0, a / c, e / c }, 1),
            G({ d / c, b / c, a / c, 0, e / c }, 1) };
        complex<double> res { -(sy[10] * sy[11]) + sy[12] * sy[13]
            + sy[12] * (sy[15] + sy[16] + sy[17]) - sy[0] * sy[1]
            + sy[10] * (sy[13] + sy[16] + sy[17] - sy[18] - sy[19]) + sy[2] * sy[3]
            - sy[2] * sy[4] + sy[0] * (-sy[3] - sy[5]) - sy[7] * sy[8]
            + (sy[9] - sy[8]) * G({ b / c, d / c }, 1)
            + sy[14]
                * (sy[15] + G({ 0, d / c, 0, e / c }, 1)
                    + 2. * G({ d / c, 0, 0, e / c }, 1))
            + sy[21]
                * (-sy[22] - sy[23] - sy[24] - sy[25]
                    - G({ 0, b / c, d / c, 0, e / c }, 1)
                    - G({ 0, d / c, 0, b / c, e / c }, 1)
                    - G({ 0, d / c, 0, e / c, b / c }, 1)
                    - G({ 0, d / c, b / c, 0, e / c }, 1)
                    - G({ b / c, 0, d / c, 0, e / c }, 1)
                    - 2. * G({ b / c, d / c, 0, 0, e / c }, 1)
                    - 2. * G({ d / c, 0, 0, b / c, e / c }, 1)
                    - 2. * G({ d / c, 0, 0, e / c, b / c }, 1)
                    - 2. * G({ d / c, 0, b / c, 0, e / c }, 1)
                    - G({ d / c, 0, e / c, 0, b / c }, 1)
                    - 2. * G({ d / c, b / c, 0, 0, e / c }, 1))
            + G({ 0, b / c, a / c, d / c, 0, e / c }, 1)
            + G({ 0, b / c, d / c, 0, a / c, e / c }, 1)
            + G({ 0, b / c, d / c, 0, e / c, a / c }, 1)
            + G({ 0, b / c, d / c, a / c, 0, e / c }, 1)
            + G({ 0, d / c, 0, b / c, a / c, e / c }, 1)
            + G({ 0, d / c, 0, b / c, e / c, a / c }, 1)
            + G({ 0, d / c, 0, e / c, b / c, a / c }, 1)
            + G({ 0, d / c, b / c, 0, a / c, e / c }, 1)
            + G({ 0, d / c, b / c, 0, e / c, a / c }, 1)
            + G({ 0, d / c, b / c, a / c, 0, e / c }, 1)
            + G({ b / c, 0, a / c, d / c, 0, e / c }, 1)
            + G({ b / c, 0, d / c, 0, a / c, e / c }, 1)
            + G({ b / c, 0, d / c, 0, e / c, a / c }, 1)
            + G({ b / c, 0, d / c, a / c, 0, e / c }, 1)
            + G({ b / c, a / c, 0, d / c, 0, e / c }, 1)
            + 2. * G({ b / c, a / c, d / c, 0, 0, e / c }, 1)
            + G({ b / c, a / c, d / c, 0, e / c, x / c }, 1)
            + 2. * G({ b / c, d / c, 0, 0, a / c, e / c }, 1)
            + 2. * G({ b / c, d / c, 0, 0, e / c, a / c }, 1)
            + 2. * G({ b / c, d / c, 0, a / c, 0, e / c }, 1)
            + G({ b / c, d / c, 0, a / c, e / c, x / c }, 1)
            + G({ b / c, d / c, 0, e / c, 0, a / c }, 1)
            + G({ b / c, d / c, 0, e / c, a / c, x / c }, 1)
            + 2. * G({ b / c, d / c, a / c, 0, 0, e / c }, 1)
            + G({ b / c, d / c, a / c, 0, e / c, x / c }, 1)
            + 2. * G({ d / c, 0, 0, b / c, a / c, e / c }, 1)
            + 2. * G({ d / c, 0, 0, b / c, e / c, a / c }, 1)
            + 2. * G({ d / c, 0, 0, e / c, b / c, a / c }, 1)
            + 2. * G({ d / c, 0, b / c, 0, a / c, e / c }, 1)
            + 2. * G({ d / c, 0, b / c, 0, e / c, a / c }, 1)
            + 2. * G({ d / c, 0, b / c, a / c, 0, e / c }, 1)
            + G({ d / c, 0, b / c, a / c, e / c, x / c }, 1)
            + G({ d / c, 0, b / c, e / c, 0, a / c }, 1)
            + G({ d / c, 0, b / c, e / c, a / c, x / c }, 1)
            + G({ d / c, 0, e / c, 0, b / c, a / c }, 1)
            + G({ d / c, 0, e / c, b / c, 0, a / c }, 1)
            + G({ d / c, 0, e / c, b / c, a / c, x / c }, 1)
            + 2. * G({ d / c, b / c, 0, 0, a / c, e / c }, 1)
            + 2. * G({ d / c, b / c, 0, 0, e / c, a / c }, 1)
            + 2. * G({ d / c, b / c, 0, a / c, 0, e / c }, 1)
            + G({ d / c, b / c, 0, a / c, e / c, x / c }, 1)
            + G({ d / c, b / c, 0, e / c, 0, a / c }, 1)
            + G({ d / c, b / c, 0, e / c, a / c, x / c }, 1)
            + 2. * G({ d / c, b / c, a / c, 0, 0, e / c }, 1)
            + G({ d / c, b / c, a / c, 0, e / c, x / c }, 1)
            + (sy[11] + sy[18] + sy[19] + G({ b / c, a / c, d / c, x / c }, 1)
                  + G({ b / c, d / c, a / c, x / c }, 1)
                  + G({ d / c, b / c, a / c, x / c }, 1))
                * G({ 0, e }, { 1, se }, x)
            + (sy[1] + sy[5] + sy[6]) * G({ 0, 0, e }, { 1, 1, se }, x)
            + sy[4]
                * (-G({ 0, a, b }, { 1, sa, sb }, x) - G({ a, 0, b }, { sa, 1, sb }, x))
            + G({ b / c, a / c }, 1) * (-sy[9] + G({ 0, d, 0, e }, { 1, sd, 1, se }, x))
            + sy[7] * G({ a, b, 0, e }, { sa, sb, 1, se }, x)
            + G({ b / c }, 1) * (sy[20] - G({ a, 0, d, 0, e }, { sa, 1, sd, 1, se }, x))
            + G({ d / c }, 1) * (-sy[20] + G({ a, b, 0, 0, e }, { sa, sb, 1, 1, se }, x))
            + G({ a, b, 0, d, 0, e }, { sa, sb, 1, sd, 1, se }, x)
            + ((sy[13] + sy[15] + sy[16] + sy[17]) * sy[21] - sy[22] - sy[23] - sy[24]
                  - sy[25] - sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31]
                  - sy[14] * sy[4])
                * Log(c, sc)
            + ((-sy[13] - sy[15] - sy[16] - sy[17]) * sy[21] + sy[22] + sy[23] + sy[24]
                  + sy[25] + sy[26] + sy[27] + sy[28] + sy[29] + sy[30] + sy[31]
                  + sy[14] * sy[4])
                * Log(-x, sc) };
        if (d != x) {
            res += (-sy[6] + G({ b / c, a / c, x / c }, 1))
                * G({ d, 0, e }, { sd, 1, se }, x);
        }
        if (e != x) {
            res += (-sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31]
                       + G({ b / c, a / c, d / c, 0, x / c }, 1)
                       + G({ b / c, d / c, 0, a / c, x / c }, 1)
                       + G({ b / c, d / c, a / c, 0, x / c }, 1)
                       + G({ d / c, 0, b / c, a / c, x / c }, 1)
                       + G({ d / c, b / c, 0, a / c, x / c }, 1)
                       + G({ d / c, b / c, a / c, 0, x / c }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_abcd0e_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (b == c) {
        if (b == d) {
            // abbbe
            const vector<complex<double>> sy = { G({ a / b, 1, 1 }, 1),
                G({ b, 0, e }, { sb, 1, se }, x), G({ a, e }, { sa, se }, x),
                G({ 0, 1, 1, a / b }, 1), G({ 0, a }, { 1, sa }, x),
                G({ 0, 1, 1, e / b }, 1), G({ 0, 1, e / b, 1 }, 1),
                G({ 0, e / b, 1, 1 }, 1), G({ 0, 1, a / b, 1 }, 1),
                G({ 0, a / b, 1, 1 }, 1), G({ b, b, 0, e }, { sb, sb, 1, se }, x),
                G({ a }, { sa }, x), G({ 0, 1, 1, e / b, a / b }, 1),
                G({ 0, 1, e / b, 1, a / b }, 1), G({ 0, 1, e / b, a / b, 1 }, 1),
                G({ 0, e / b, 1, 1, a / b }, 1), G({ 0, e / b, 1, a / b, 1 }, 1),
                G({ 0, e / b, a / b, 1, 1 }, 1), G({ 0, 1, 1, a / b, e / b }, 1),
                G({ 0, 1, a / b, 1, e / b }, 1), G({ 0, 1, a / b, e / b, 1 }, 1),
                G({ 0, a / b, 1, 1, e / b }, 1), G({ 0, a / b, 1, e / b, 1 }, 1),
                G({ 0, a / b, e / b, 1, 1 }, 1), G({ a / b, 0, 1, 1, e / b }, 1),
                G({ a / b, 0, 1, e / b, 1 }, 1), G({ a / b, 0, e / b, 1, 1 }, 1) };
            complex<double> res { -(sy[0] * sy[1]) + sy[2] * sy[3] - sy[4] * sy[5]
                + sy[4] * (-sy[6] - sy[7])
                + sy[2] * (sy[9] - sy[5] - sy[6] - sy[7] + sy[8])
                - sy[10] * G({ a / b, x / b }, 1) + sy[1] * G({ a / b, x / b, 1 }, 1)
                + sy[11]
                    * (sy[12] + sy[13] + sy[14] + sy[15] + sy[16] + sy[17]
                        + 2. * G({ 0, 0, 1, 1, e / b }, 1)
                        + 2. * G({ 0, 0, 1, e / b, 1 }, 1)
                        + 2. * G({ 0, 0, e / b, 1, 1 }, 1) + G({ 0, 1, 0, 1, e / b }, 1)
                        + G({ 0, 1, 0, e / b, 1 }, 1) + G({ 0, 1, 1, 0, e / b }, 1))
                - 2. * G({ 0, 0, 1, 1, a / b, e / b }, 1)
                - 2. * G({ 0, 0, 1, 1, e / b, a / b }, 1)
                - 2. * G({ 0, 0, 1, a / b, 1, e / b }, 1)
                - 2. * G({ 0, 0, 1, a / b, e / b, 1 }, 1)
                - 2. * G({ 0, 0, 1, e / b, 1, a / b }, 1)
                - 2. * G({ 0, 0, 1, e / b, a / b, 1 }, 1)
                - 2. * G({ 0, 0, a / b, 1, 1, e / b }, 1)
                - 2. * G({ 0, 0, a / b, 1, e / b, 1 }, 1)
                - 2. * G({ 0, 0, a / b, e / b, 1, 1 }, 1)
                - 2. * G({ 0, 0, e / b, 1, 1, a / b }, 1)
                - 2. * G({ 0, 0, e / b, 1, a / b, 1 }, 1)
                - 2. * G({ 0, 0, e / b, a / b, 1, 1 }, 1)
                - G({ 0, 1, 0, 1, a / b, e / b }, 1) - G({ 0, 1, 0, 1, e / b, a / b }, 1)
                - G({ 0, 1, 0, a / b, 1, e / b }, 1) - G({ 0, 1, 0, a / b, e / b, 1 }, 1)
                - G({ 0, 1, 0, e / b, 1, a / b }, 1) - G({ 0, 1, 0, e / b, a / b, 1 }, 1)
                - G({ 0, 1, 1, 0, a / b, e / b }, 1) - G({ 0, 1, 1, 0, e / b, a / b }, 1)
                - G({ 0, 1, 1, a / b, 0, e / b }, 1)
                - G({ 0, 1, 1, a / b, e / b, x / b }, 1)
                - G({ 0, 1, 1, e / b, 0, a / b }, 1)
                - G({ 0, 1, 1, e / b, a / b, x / b }, 1)
                - G({ 0, 1, a / b, 0, 1, e / b }, 1) - G({ 0, 1, a / b, 0, e / b, 1 }, 1)
                - G({ 0, 1, a / b, 1, 0, e / b }, 1)
                - G({ 0, 1, a / b, 1, e / b, x / b }, 1)
                - G({ 0, 1, a / b, e / b, 1, x / b }, 1)
                - G({ 0, 1, a / b, e / b, x / b, 1 }, 1)
                - G({ 0, 1, e / b, 0, 1, a / b }, 1) - G({ 0, 1, e / b, 0, a / b, 1 }, 1)
                - G({ 0, 1, e / b, 1, 0, a / b }, 1)
                - G({ 0, 1, e / b, 1, a / b, x / b }, 1)
                - G({ 0, 1, e / b, a / b, 1, x / b }, 1)
                - G({ 0, 1, e / b, a / b, x / b, 1 }, 1)
                - 2. * G({ 0, a / b, 0, 1, 1, e / b }, 1)
                - 2. * G({ 0, a / b, 0, 1, e / b, 1 }, 1)
                - 2. * G({ 0, a / b, 0, e / b, 1, 1 }, 1)
                - G({ 0, a / b, 1, 0, 1, e / b }, 1) - G({ 0, a / b, 1, 0, e / b, 1 }, 1)
                - G({ 0, a / b, 1, 1, 0, e / b }, 1)
                - G({ 0, a / b, 1, 1, e / b, x / b }, 1)
                - G({ 0, a / b, 1, e / b, 1, x / b }, 1)
                - G({ 0, a / b, 1, e / b, x / b, 1 }, 1)
                - G({ 0, a / b, e / b, 1, 1, x / b }, 1)
                - G({ 0, a / b, e / b, 1, x / b, 1 }, 1)
                - G({ 0, a / b, e / b, x / b, 1, 1 }, 1)
                - G({ 0, e / b, 0, 1, 1, a / b }, 1) - G({ 0, e / b, 0, 1, a / b, 1 }, 1)
                - G({ 0, e / b, 0, a / b, 1, 1 }, 1) - G({ 0, e / b, 1, 0, 1, a / b }, 1)
                - G({ 0, e / b, 1, 0, a / b, 1 }, 1) - G({ 0, e / b, 1, 1, 0, a / b }, 1)
                - G({ 0, e / b, 1, 1, a / b, x / b }, 1)
                - G({ 0, e / b, 1, a / b, 1, x / b }, 1)
                - G({ 0, e / b, 1, a / b, x / b, 1 }, 1)
                - G({ 0, e / b, a / b, 1, 1, x / b }, 1)
                - G({ 0, e / b, a / b, 1, x / b, 1 }, 1)
                - G({ 0, e / b, a / b, x / b, 1, 1 }, 1)
                - 2. * G({ a / b, 0, 0, 1, 1, e / b }, 1)
                - 2. * G({ a / b, 0, 0, 1, e / b, 1 }, 1)
                - 2. * G({ a / b, 0, 0, e / b, 1, 1 }, 1)
                - G({ a / b, 0, 1, 0, 1, e / b }, 1) - G({ a / b, 0, 1, 0, e / b, 1 }, 1)
                - G({ a / b, 0, 1, 1, 0, e / b }, 1)
                - G({ a / b, 0, 1, 1, e / b, x / b }, 1)
                - G({ a / b, 0, 1, e / b, 1, x / b }, 1)
                - G({ a / b, 0, 1, e / b, x / b, 1 }, 1)
                - G({ a / b, 0, e / b, 1, 1, x / b }, 1)
                - G({ a / b, 0, e / b, 1, x / b, 1 }, 1)
                - G({ a / b, 0, e / b, x / b, 1, 1 }, 1)
                + (-sy[9] - sy[3] - sy[8] - G({ a / b, x / b, 1, 1 }, 1))
                    * G({ 0, e }, { 1, se }, x)
                + sy[0] * G({ a, 0, e }, { sa, 1, se }, x)
                + G({ a / b, 1 }, 1) * (sy[10] - G({ a, b, 0, e }, { sa, sb, 1, se }, x))
                + G({ a / b }, 1)
                    * (-G({ 0, b, b, 0, e }, { 1, sb, sb, 1, se }, x)
                        + G({ a, b, b, 0, e }, { sa, sb, sb, 1, se }, x))
                + G({ a, 0, b, b, 0, e }, { sa, 1, sb, sb, 1, se }, x)
                + (sy[12] + sy[13] + sy[14] + sy[15] + sy[16] + sy[17] + sy[18] + sy[19]
                      + sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25] + sy[26]
                      + sy[11] * (-sy[5] - sy[6] - sy[7]))
                    * Log(b, sb)
                + (-sy[12] - sy[13] - sy[14] - sy[15] - sy[16] - sy[17] - sy[18] - sy[19]
                      - sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25] - sy[26]
                      + sy[11] * (sy[5] + sy[6] + sy[7]))
                    * Log(-x, sb) };

            if (e != x) {
                res += (sy[18] + sy[19] + sy[20] + sy[21] + sy[22] + sy[23] + sy[24]
                           + sy[25] + sy[26] - G({ 0, 1, 1, a / b, x / b }, 1)
                           - G({ 0, 1, a / b, 1, x / b }, 1)
                           - G({ 0, 1, a / b, x / b, 1 }, 1)
                           - G({ 0, a / b, 1, 1, x / b }, 1)
                           - G({ 0, a / b, 1, x / b, 1 }, 1)
                           - G({ 0, a / b, x / b, 1, 1 }, 1)
                           - G({ a / b, 0, 1, 1, x / b }, 1)
                           - G({ a / b, 0, 1, x / b, 1 }, 1)
                           - G({ a / b, 0, x / b, 1, 1 }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        } else { // abbde
            const vector<complex<double>> sy = { G({ a, 0, e }, { sa, 1, se }, x),
                G({ d / b, 1, a / b }, 1), G({ d / b, a / b, 1 }, 1),
                G({ a / b, d / b, 1 }, 1), G({ a, d, 0, e }, { sa, sd, 1, se }, x),
                G({ b, d, 0, e }, { sb, sd, 1, se }, x), G({ a, e }, { sa, se }, x),
                G({ d / b, 0, 1, a / b }, 1), G({ 0, a }, { 1, sa }, x),
                G({ d / b, 0, 1, e / b }, 1), G({ d / b, 0, e / b, 1 }, 1),
                G({ d / b, 1, 0, e / b }, 1), G({ d / b, 0, a / b, 1 }, 1),
                G({ d / b, 1, 0, a / b }, 1), G({ a }, { sa }, x),
                G({ d / b, 0, 1, e / b, a / b }, 1), G({ d / b, 0, e / b, 1, a / b }, 1),
                G({ d / b, 0, e / b, a / b, 1 }, 1), G({ d / b, 1, 0, e / b, a / b }, 1),
                G({ a / b, d / b, 0, 1, e / b }, 1), G({ a / b, d / b, 0, e / b, 1 }, 1),
                G({ a / b, d / b, 1, 0, e / b }, 1), G({ d / b, 0, 1, a / b, e / b }, 1),
                G({ d / b, 0, a / b, 1, e / b }, 1), G({ d / b, 0, a / b, e / b, 1 }, 1),
                G({ d / b, 1, 0, a / b, e / b }, 1), G({ d / b, 1, a / b, 0, e / b }, 1),
                G({ d / b, a / b, 0, 1, e / b }, 1), G({ d / b, a / b, 0, e / b, 1 }, 1),
                G({ d / b, a / b, 1, 0, e / b }, 1) };
            complex<double> res { -(sy[0] * sy[1]) - sy[0] * sy[2]
                + (sy[9] + sy[10] + sy[11] - sy[12] - sy[13]) * sy[6] - sy[6] * sy[7]
                + sy[9] * sy[8] + (sy[10] + sy[11]) * sy[8]
                + (-sy[4] + sy[5]) * G({ a / b, 1 }, 1) - sy[5] * G({ a / b, x / b }, 1)
                + sy[14]
                    * (-sy[15] - sy[16] - sy[17] - sy[18]
                        - G({ 0, d / b, 0, 1, e / b }, 1)
                        - G({ 0, d / b, 0, e / b, 1 }, 1)
                        - G({ 0, d / b, 1, 0, e / b }, 1)
                        - 2. * G({ d / b, 0, 0, 1, e / b }, 1)
                        - 2. * G({ d / b, 0, 0, e / b, 1 }, 1)
                        - 2. * G({ d / b, 0, 1, 0, e / b }, 1)
                        - 2. * G({ d / b, 1, 0, 0, e / b }, 1))
                + G({ 0, a / b, d / b, 0, 1, e / b }, 1)
                + G({ 0, a / b, d / b, 0, e / b, 1 }, 1)
                + G({ 0, a / b, d / b, 1, 0, e / b }, 1)
                + G({ 0, d / b, 0, 1, a / b, e / b }, 1)
                + G({ 0, d / b, 0, 1, e / b, a / b }, 1)
                + G({ 0, d / b, 0, a / b, 1, e / b }, 1)
                + G({ 0, d / b, 0, a / b, e / b, 1 }, 1)
                + G({ 0, d / b, 0, e / b, 1, a / b }, 1)
                + G({ 0, d / b, 0, e / b, a / b, 1 }, 1)
                + G({ 0, d / b, 1, 0, a / b, e / b }, 1)
                + G({ 0, d / b, 1, 0, e / b, a / b }, 1)
                + G({ 0, d / b, 1, a / b, 0, e / b }, 1)
                + G({ 0, d / b, a / b, 0, 1, e / b }, 1)
                + G({ 0, d / b, a / b, 0, e / b, 1 }, 1)
                + G({ 0, d / b, a / b, 1, 0, e / b }, 1)
                + G({ a / b, 0, d / b, 0, 1, e / b }, 1)
                + G({ a / b, 0, d / b, 0, e / b, 1 }, 1)
                + G({ a / b, 0, d / b, 1, 0, e / b }, 1)
                + 2. * G({ a / b, d / b, 0, 0, 1, e / b }, 1)
                + 2. * G({ a / b, d / b, 0, 0, e / b, 1 }, 1)
                + 2. * G({ a / b, d / b, 0, 1, 0, e / b }, 1)
                + G({ a / b, d / b, 0, 1, e / b, x / b }, 1)
                + G({ a / b, d / b, 0, e / b, 1, x / b }, 1)
                + G({ a / b, d / b, 0, e / b, x / b, 1 }, 1)
                + 2. * G({ a / b, d / b, 1, 0, 0, e / b }, 1)
                + G({ a / b, d / b, 1, 0, e / b, x / b }, 1)
                + 2. * G({ d / b, 0, 0, 1, a / b, e / b }, 1)
                + 2. * G({ d / b, 0, 0, 1, e / b, a / b }, 1)
                + 2. * G({ d / b, 0, 0, a / b, 1, e / b }, 1)
                + 2. * G({ d / b, 0, 0, a / b, e / b, 1 }, 1)
                + 2. * G({ d / b, 0, 0, e / b, 1, a / b }, 1)
                + 2. * G({ d / b, 0, 0, e / b, a / b, 1 }, 1)
                + 2. * G({ d / b, 0, 1, 0, a / b, e / b }, 1)
                + 2. * G({ d / b, 0, 1, 0, e / b, a / b }, 1)
                + 2. * G({ d / b, 0, 1, a / b, 0, e / b }, 1)
                + G({ d / b, 0, 1, a / b, e / b, x / b }, 1)
                + G({ d / b, 0, 1, e / b, 0, a / b }, 1)
                + G({ d / b, 0, 1, e / b, a / b, x / b }, 1)
                + 2. * G({ d / b, 0, a / b, 0, 1, e / b }, 1)
                + 2. * G({ d / b, 0, a / b, 0, e / b, 1 }, 1)
                + 2. * G({ d / b, 0, a / b, 1, 0, e / b }, 1)
                + G({ d / b, 0, a / b, 1, e / b, x / b }, 1)
                + G({ d / b, 0, a / b, e / b, 1, x / b }, 1)
                + G({ d / b, 0, a / b, e / b, x / b, 1 }, 1)
                + G({ d / b, 0, e / b, 0, 1, a / b }, 1)
                + G({ d / b, 0, e / b, 0, a / b, 1 }, 1)
                + G({ d / b, 0, e / b, 1, 0, a / b }, 1)
                + G({ d / b, 0, e / b, 1, a / b, x / b }, 1)
                + G({ d / b, 0, e / b, a / b, 1, x / b }, 1)
                + G({ d / b, 0, e / b, a / b, x / b, 1 }, 1)
                + 2. * G({ d / b, 1, 0, 0, a / b, e / b }, 1)
                + 2. * G({ d / b, 1, 0, 0, e / b, a / b }, 1)
                + 2. * G({ d / b, 1, 0, a / b, 0, e / b }, 1)
                + G({ d / b, 1, 0, a / b, e / b, x / b }, 1)
                + G({ d / b, 1, 0, e / b, 0, a / b }, 1)
                + G({ d / b, 1, 0, e / b, a / b, x / b }, 1)
                + 2. * G({ d / b, 1, a / b, 0, 0, e / b }, 1)
                + G({ d / b, 1, a / b, 0, e / b, x / b }, 1)
                + 2. * G({ d / b, a / b, 0, 0, 1, e / b }, 1)
                + 2. * G({ d / b, a / b, 0, 0, e / b, 1 }, 1)
                + 2. * G({ d / b, a / b, 0, 1, 0, e / b }, 1)
                + G({ d / b, a / b, 0, 1, e / b, x / b }, 1)
                + G({ d / b, a / b, 0, e / b, 1, x / b }, 1)
                + G({ d / b, a / b, 0, e / b, x / b, 1 }, 1)
                + 2. * G({ d / b, a / b, 1, 0, 0, e / b }, 1)
                + G({ d / b, a / b, 1, 0, e / b, x / b }, 1)
                + (sy[12] + sy[13] + sy[7] + G({ a / b, d / b, 1, x / b }, 1)
                      + G({ a / b, d / b, x / b, 1 }, 1)
                      + G({ d / b, 1, a / b, x / b }, 1)
                      + G({ d / b, a / b, 1, x / b }, 1)
                      + G({ d / b, a / b, x / b, 1 }, 1))
                    * G({ 0, e }, { 1, se }, x)
                + (sy[1] + sy[2] + sy[3]) * G({ 0, 0, e }, { 1, 1, se }, x)
                + G({ d / b, 1 }, 1) * (sy[4] - G({ a, 0, 0, e }, { sa, 1, 1, se }, x))
                + G({ a / b }, 1)
                    * (-G({ 0, b, d, 0, e }, { 1, sb, sd, 1, se }, x)
                        + G({ a, b, d, 0, e }, { sa, sb, sd, 1, se }, x))
                + G({ a, 0, b, d, 0, e }, { sa, 1, sb, sd, 1, se }, x)
                + ((sy[9] + sy[10] + sy[11]) * sy[14] - sy[15] - sy[16] - sy[17] - sy[18]
                      - sy[19] - sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25]
                      - sy[26] - sy[27] - sy[28] - sy[29])
                    * Log(b, sb)
                + ((-sy[9] - sy[10] - sy[11]) * sy[14] + sy[15] + sy[16] + sy[17] + sy[18]
                      + sy[19] + sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25]
                      + sy[26] + sy[27] + sy[28] + sy[29])
                    * Log(-x, sb) };
            if (d != x) {
                res += (-sy[3] + G({ a / b, x / b, 1 }, 1))
                    * G({ d, 0, e }, { sd, 1, se }, x);
            }
            if (e != x) {
                res += (-sy[19] - sy[20] - sy[21] - sy[22] - sy[23] - sy[24] - sy[25]
                           - sy[26] - sy[27] - sy[28] - sy[29]
                           + G({ a / b, d / b, 0, 1, x / b }, 1)
                           + G({ a / b, d / b, 0, x / b, 1 }, 1)
                           + G({ a / b, d / b, 1, 0, x / b }, 1)
                           + G({ d / b, 0, 1, a / b, x / b }, 1)
                           + G({ d / b, 0, a / b, 1, x / b }, 1)
                           + G({ d / b, 0, a / b, x / b, 1 }, 1)
                           + G({ d / b, 1, 0, a / b, x / b }, 1)
                           + G({ d / b, 1, a / b, 0, x / b }, 1)
                           + G({ d / b, a / b, 0, 1, x / b }, 1)
                           + G({ d / b, a / b, 0, x / b, 1 }, 1)
                           + G({ d / b, a / b, 1, 0, x / b }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else { // abcde
        const vector<complex<double>> sy
            = { G({ a / b, c / b, d / b }, 1), G({ c / b, a / b, d / b }, 1),
                  G({ c / b, d / b, a / b }, 1), G({ a / b, c / b }, 1),
                  G({ 0, d, 0, e }, { 1, sd, 1, se }, x), G({ c / b, a / b }, 1),
                  G({ a, d, 0, e }, { sa, sd, 1, se }, x), G({ a, e }, { sa, se }, x),
                  G({ c / b, d / b, 0, a / b }, 1), G({ c / b, d / b, 0, e / b }, 1),
                  G({ a, c, d, 0, e }, { sa, sc, sd, 1, se }, x), G({ a }, { sa }, x),
                  G({ c / b, d / b, 0, e / b, a / b }, 1),
                  G({ a / b, c / b, d / b, 0, e / b }, 1),
                  G({ c / b, a / b, d / b, 0, e / b }, 1),
                  G({ c / b, d / b, 0, a / b, e / b }, 1),
                  G({ c / b, d / b, a / b, 0, e / b }, 1) };
        complex<double> res { -(sy[3] * sy[4]) - sy[4] * sy[5] + sy[5] * sy[6]
            - sy[9] * sy[7] + sy[7] * sy[8]
            + sy[11]
                * (sy[12] + G({ 0, c / b, d / b, 0, e / b }, 1)
                    + G({ c / b, 0, d / b, 0, e / b }, 1)
                    + 2. * G({ c / b, d / b, 0, 0, e / b }, 1))
            - G({ 0, a / b, c / b, d / b, 0, e / b }, 1)
            - G({ 0, c / b, a / b, d / b, 0, e / b }, 1)
            - G({ 0, c / b, d / b, 0, a / b, e / b }, 1)
            - G({ 0, c / b, d / b, 0, e / b, a / b }, 1)
            - G({ 0, c / b, d / b, a / b, 0, e / b }, 1)
            - G({ a / b, 0, c / b, d / b, 0, e / b }, 1)
            - G({ a / b, c / b, 0, d / b, 0, e / b }, 1)
            - 2. * G({ a / b, c / b, d / b, 0, 0, e / b }, 1)
            - G({ a / b, c / b, d / b, 0, e / b, x / b }, 1)
            - G({ c / b, 0, a / b, d / b, 0, e / b }, 1)
            - G({ c / b, 0, d / b, 0, a / b, e / b }, 1)
            - G({ c / b, 0, d / b, 0, e / b, a / b }, 1)
            - G({ c / b, 0, d / b, a / b, 0, e / b }, 1)
            - G({ c / b, a / b, 0, d / b, 0, e / b }, 1)
            - 2. * G({ c / b, a / b, d / b, 0, 0, e / b }, 1)
            - G({ c / b, a / b, d / b, 0, e / b, x / b }, 1)
            - 2. * G({ c / b, d / b, 0, 0, a / b, e / b }, 1)
            - 2. * G({ c / b, d / b, 0, 0, e / b, a / b }, 1)
            - 2. * G({ c / b, d / b, 0, a / b, 0, e / b }, 1)
            - G({ c / b, d / b, 0, a / b, e / b, x / b }, 1)
            - G({ c / b, d / b, 0, e / b, 0, a / b }, 1)
            - G({ c / b, d / b, 0, e / b, a / b, x / b }, 1)
            - 2. * G({ c / b, d / b, a / b, 0, 0, e / b }, 1)
            - G({ c / b, d / b, a / b, 0, e / b, x / b }, 1)
            - sy[9] * G({ 0, a }, { 1, sa }, x)
            + (-sy[8] - G({ a / b, c / b, d / b, x / b }, 1)
                  - G({ c / b, a / b, d / b, x / b }, 1)
                  - G({ c / b, d / b, a / b, x / b }, 1))
                * G({ 0, e }, { 1, se }, x)
            + (-sy[0] - sy[1] - sy[2]) * G({ 0, 0, e }, { 1, 1, se }, x)
            + sy[2] * G({ a, 0, e }, { sa, 1, se }, x)
            + G({ c / b, d / b }, 1) * (-sy[6] + G({ a, 0, 0, e }, { sa, 1, 1, se }, x))
            + G({ a / b }, 1) * (sy[10] - G({ 0, c, d, 0, e }, { 1, sc, sd, 1, se }, x))
            + G({ c / b }, 1) * (-sy[10] + G({ a, 0, d, 0, e }, { sa, 1, sd, 1, se }, x))
            + G({ a, 0, c, d, 0, e }, { sa, 1, sc, sd, 1, se }, x)
            + (-(sy[9] * sy[11]) + sy[12] + sy[13] + sy[14] + sy[15] + sy[16])
                * Log(b, sb)
            + (sy[9] * sy[11] - sy[12] - sy[13] - sy[14] - sy[15] - sy[16])
                * Log(-x, sb) };
        if (c != x) {
            res += (sy[3] - G({ a / b, x / b }, 1))
                * G({ c, d, 0, e }, { sc, sd, 1, se }, x);
        }
        if (d != x) {
            res += (sy[0] + sy[1] - G({ a / b, c / b, x / b }, 1)
                       - G({ c / b, a / b, x / b }, 1))
                * G({ d, 0, e }, { sd, 1, se }, x);
        }
        if (e != x) {
            res += (sy[13] + sy[14] + sy[15] + sy[16]
                       - G({ a / b, c / b, d / b, 0, x / b }, 1)
                       - G({ c / b, a / b, d / b, 0, x / b }, 1)
                       - G({ c / b, d / b, 0, a / b, x / b }, 1)
                       - G({ c / b, d / b, a / b, 0, x / b }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}
complex<double> G6_abcd0e_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, int sa, int sb, int sc, int sd, int se,
    double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    double x = 1.;
    if (a == b) {
        if (a == c) {
            if (a == d) { // aaaae
                const vector<complex<double>> sy
                    = { G({ 0, 1, 1, 1, e / a }, 1), G({ 0, 1, 1, e / a, 1 }, 1),
                          G({ 0, 1, e / a, 1, 1 }, 1), G({ 0, e / a, 1, 1, 1 }, 1) };
                complex<double> res { -2. * G({ 0, 0, 1, 1, 1, e / a }, 1)
                    - 2. * G({ 0, 0, 1, 1, e / a, 1 }, 1)
                    - 2. * G({ 0, 0, 1, e / a, 1, 1 }, 1)
                    - 2. * G({ 0, 0, e / a, 1, 1, 1 }, 1) - G({ 0, 1, 0, 1, 1, e / a }, 1)
                    - G({ 0, 1, 0, 1, e / a, 1 }, 1) - G({ 0, 1, 0, e / a, 1, 1 }, 1)
                    - G({ 0, 1, 1, 0, 1, e / a }, 1) - G({ 0, 1, 1, 0, e / a, 1 }, 1)
                    - G({ 0, 1, 1, 1, 0, e / a }, 1) - G({ 0, 1, 1, 1, e / a, x / a }, 1)
                    - G({ 0, 1, 1, e / a, 1, x / a }, 1)
                    - G({ 0, 1, 1, e / a, x / a, 1 }, 1)
                    - G({ 0, 1, e / a, 1, 1, x / a }, 1)
                    - G({ 0, 1, e / a, 1, x / a, 1 }, 1)
                    - G({ 0, 1, e / a, x / a, 1, 1 }, 1)
                    - G({ 0, e / a, 1, 1, 1, x / a }, 1)
                    - G({ 0, e / a, 1, 1, x / a, 1 }, 1)
                    - G({ 0, e / a, 1, x / a, 1, 1 }, 1)
                    - G({ 0, e / a, x / a, 1, 1, 1 }, 1)
                    - G({ x / a, 1, 1, 1 }, 1) * G({ 0, e }, { 1, se }, x)
                    + G({ x / a, 1, 1 }, 1) * G({ a, 0, e }, { sa, 1, se }, x)
                    - G({ x / a, 1 }, 1) * G({ a, a, 0, e }, { sa, sa, 1, se }, x)
                    + G({ x / a }, 1) * G({ a, a, a, 0, e }, { sa, sa, sa, 1, se }, x)
                    + G({ 0, a, a, a, 0, e }, { 1, sa, sa, sa, 1, se }, x)
                    + (sy[0] + sy[1] + sy[2] + sy[3]) * Log(a, sa)
                    + (-sy[0] - sy[1] - sy[2] - sy[3]) * Log(-x, sa) };
                if (e != x) {
                    res += (sy[0] + sy[1] + sy[2] + sy[3] - G({ 0, 1, 1, 1, x / a }, 1)
                               - G({ 0, 1, 1, x / a, 1 }, 1) - G({ 0, 1, x / a, 1, 1 }, 1)
                               - G({ 0, x / a, 1, 1, 1 }, 1))
                        * G({ e }, { se }, x);
                }
                return res;
            } else { // aaade
                const vector<complex<double>> sy = { G({ d / a, 1, 1 }, 1),
                    G({ d / a, 0, 1, 1, e / a }, 1), G({ d / a, 0, 1, e / a, 1 }, 1),
                    G({ d / a, 0, e / a, 1, 1 }, 1), G({ d / a, 1, 0, 1, e / a }, 1),
                    G({ d / a, 1, 0, e / a, 1 }, 1), G({ d / a, 1, 1, 0, e / a }, 1) };
                complex<double> res { G({ 0, d / a, 0, 1, 1, e / a }, 1)
                    + G({ 0, d / a, 0, 1, e / a, 1 }, 1)
                    + G({ 0, d / a, 0, e / a, 1, 1 }, 1)
                    + G({ 0, d / a, 1, 0, 1, e / a }, 1)
                    + G({ 0, d / a, 1, 0, e / a, 1 }, 1)
                    + G({ 0, d / a, 1, 1, 0, e / a }, 1)
                    + 2. * G({ d / a, 0, 0, 1, 1, e / a }, 1)
                    + 2. * G({ d / a, 0, 0, 1, e / a, 1 }, 1)
                    + 2. * G({ d / a, 0, 0, e / a, 1, 1 }, 1)
                    + 2. * G({ d / a, 0, 1, 0, 1, e / a }, 1)
                    + 2. * G({ d / a, 0, 1, 0, e / a, 1 }, 1)
                    + 2. * G({ d / a, 0, 1, 1, 0, e / a }, 1)
                    + G({ d / a, 0, 1, 1, e / a, x / a }, 1)
                    + G({ d / a, 0, 1, e / a, 1, x / a }, 1)
                    + G({ d / a, 0, 1, e / a, x / a, 1 }, 1)
                    + G({ d / a, 0, e / a, 1, 1, x / a }, 1)
                    + G({ d / a, 0, e / a, 1, x / a, 1 }, 1)
                    + G({ d / a, 0, e / a, x / a, 1, 1 }, 1)
                    + 2. * G({ d / a, 1, 0, 0, 1, e / a }, 1)
                    + 2. * G({ d / a, 1, 0, 0, e / a, 1 }, 1)
                    + 2. * G({ d / a, 1, 0, 1, 0, e / a }, 1)
                    + G({ d / a, 1, 0, 1, e / a, x / a }, 1)
                    + G({ d / a, 1, 0, e / a, 1, x / a }, 1)
                    + G({ d / a, 1, 0, e / a, x / a, 1 }, 1)
                    + 2. * G({ d / a, 1, 1, 0, 0, e / a }, 1)
                    + G({ d / a, 1, 1, 0, e / a, x / a }, 1)
                    + (G({ d / a, 1, 1, x / a }, 1) + G({ d / a, 1, x / a, 1 }, 1)
                          + G({ d / a, x / a, 1, 1 }, 1))
                        * G({ 0, e }, { 1, se }, x)
                    + sy[0] * G({ 0, 0, e }, { 1, 1, se }, x)
                    - G({ x / a, 1 }, 1) * G({ a, d, 0, e }, { sa, sd, 1, se }, x)
                    + G({ x / a }, 1) * G({ a, a, d, 0, e }, { sa, sa, sd, 1, se }, x)
                    + G({ 0, a, a, d, 0, e }, { 1, sa, sa, sd, 1, se }, x)
                    + (-sy[1] - sy[2] - sy[3] - sy[4] - sy[5] - sy[6]) * Log(a, sa)
                    + (sy[1] + sy[2] + sy[3] + sy[4] + sy[5] + sy[6]) * Log(-x, sa) };
                if (d != x) {
                    res += (-sy[0] + G({ x / a, 1, 1 }, 1))
                        * G({ d, 0, e }, { sd, 1, se }, x);
                }
                if (e != x) {
                    res += (-sy[1] - sy[2] - sy[3] - sy[4] - sy[5] - sy[6]
                               + G({ d / a, 0, 1, 1, x / a }, 1)
                               + G({ d / a, 0, 1, x / a, 1 }, 1)
                               + G({ d / a, 0, x / a, 1, 1 }, 1)
                               + G({ d / a, 1, 0, 1, x / a }, 1)
                               + G({ d / a, 1, 0, x / a, 1 }, 1)
                               + G({ d / a, 1, 1, 0, x / a }, 1))
                        * G({ e }, { se }, x);
                }
                return res;
            }
        } else { // aacde
            const vector<complex<double>> sy = { G({ c / a, 1, d / a }, 1),
                G({ c / a, d / a, 1 }, 1), G({ c / a, 1 }, 1),
                G({ c / a, 1, d / a, 0, e / a }, 1), G({ c / a, d / a, 0, 1, e / a }, 1),
                G({ c / a, d / a, 0, e / a, 1 }, 1),
                G({ c / a, d / a, 1, 0, e / a }, 1) };
            complex<double> res { -G({ 0, c / a, 1, d / a, 0, e / a }, 1)
                - G({ 0, c / a, d / a, 0, 1, e / a }, 1)
                - G({ 0, c / a, d / a, 0, e / a, 1 }, 1)
                - G({ 0, c / a, d / a, 1, 0, e / a }, 1)
                - G({ c / a, 0, 1, d / a, 0, e / a }, 1)
                - G({ c / a, 0, d / a, 0, 1, e / a }, 1)
                - G({ c / a, 0, d / a, 0, e / a, 1 }, 1)
                - G({ c / a, 0, d / a, 1, 0, e / a }, 1)
                - G({ c / a, 1, 0, d / a, 0, e / a }, 1)
                - 2. * G({ c / a, 1, d / a, 0, 0, e / a }, 1)
                - G({ c / a, 1, d / a, 0, e / a, x / a }, 1)
                - 2. * G({ c / a, d / a, 0, 0, 1, e / a }, 1)
                - 2. * G({ c / a, d / a, 0, 0, e / a, 1 }, 1)
                - 2. * G({ c / a, d / a, 0, 1, 0, e / a }, 1)
                - G({ c / a, d / a, 0, 1, e / a, x / a }, 1)
                - G({ c / a, d / a, 0, e / a, 1, x / a }, 1)
                - G({ c / a, d / a, 0, e / a, x / a, 1 }, 1)
                - 2. * G({ c / a, d / a, 1, 0, 0, e / a }, 1)
                - G({ c / a, d / a, 1, 0, e / a, x / a }, 1)
                + (-G({ c / a, 1, d / a, x / a }, 1) - G({ c / a, d / a, 1, x / a }, 1)
                      - G({ c / a, d / a, x / a, 1 }, 1))
                    * G({ 0, e }, { 1, se }, x)
                + (-sy[0] - sy[1]) * G({ 0, 0, e }, { 1, 1, se }, x)
                - sy[2] * G({ 0, d, 0, e }, { 1, sd, 1, se }, x)
                + G({ x / a }, 1) * G({ a, c, d, 0, e }, { sa, sc, sd, 1, se }, x)
                + G({ 0, a, c, d, 0, e }, { 1, sa, sc, sd, 1, se }, x)
                + (sy[3] + sy[4] + sy[5] + sy[6]) * Log(a, sa)
                + (-sy[3] - sy[4] - sy[5] - sy[6]) * Log(-x, sa) };
            if (c != x) {
                res += (sy[2] - G({ x / a, 1 }, 1))
                    * G({ c, d, 0, e }, { sc, sd, 1, se }, x);
            }
            if (d != x) {
                res += (sy[0] + sy[1] - G({ c / a, 1, x / a }, 1)
                           - G({ c / a, x / a, 1 }, 1))
                    * G({ d, 0, e }, { sd, 1, se }, x);
            }
            if (e != x) {
                res += (sy[3] + sy[4] + sy[5] + sy[6]
                           - G({ c / a, 1, d / a, 0, x / a }, 1)
                           - G({ c / a, d / a, 0, 1, x / a }, 1)
                           - G({ c / a, d / a, 0, x / a, 1 }, 1)
                           - G({ c / a, d / a, 1, 0, x / a }, 1))
                    * G({ e }, { se }, x);
            }
            return res;
        }
    } else { // abcde
        const vector<complex<double>> sy
            = { G({ b / a, c / a, d / a }, 1), G({ b / a, c / a }, 1), G({ b / a }, 1),
                  G({ b / a, c / a, d / a, 0, e / a }, 1) };
        complex<double> res { G({ 0, b / a, c / a, d / a, 0, e / a }, 1)
            + G({ b / a, 0, c / a, d / a, 0, e / a }, 1)
            + G({ b / a, c / a, 0, d / a, 0, e / a }, 1)
            + 2. * G({ b / a, c / a, d / a, 0, 0, e / a }, 1)
            + G({ b / a, c / a, d / a, 0, e / a, x / a }, 1)
            + G({ b / a, c / a, d / a, x / a }, 1) * G({ 0, e }, { 1, se }, x)
            + sy[0] * G({ 0, 0, e }, { 1, 1, se }, x)
            + sy[1] * G({ 0, d, 0, e }, { 1, sd, 1, se }, x)
            + sy[2] * G({ 0, c, d, 0, e }, { 1, sc, sd, 1, se }, x)
            + G({ 0, b, c, d, 0, e }, { 1, sb, sc, sd, 1, se }, x) - sy[3] * Log(a, sa)
            + sy[3] * Log(-x, sa) };
        if (b != x) {
            res += (-sy[2] + G({ x / a }, 1))
                * G({ b, c, d, 0, e }, { sb, sc, sd, 1, se }, x);
        }
        if (c != x) {
            res += (-sy[1] + G({ b / a, x / a }, 1))
                * G({ c, d, 0, e }, { sc, sd, 1, se }, x);
        }
        if (d != x) {
            res += (-sy[0] + G({ b / a, c / a, x / a }, 1))
                * G({ d, 0, e }, { sd, 1, se }, x);
        }
        if (e != x) {
            res += (-sy[3] + G({ b / a, c / a, d / a, 0, x / a }, 1))
                * G({ e }, { se }, x);
        }
        return res;
    }
}

complex<double> G6_0abcde(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{
    // e is smallest
    if (abs(e) < abs(d) && abs(e) < abs(c) && abs(e) < abs(b) && abs(e) < abs(a))
        return G6_0abcde_e(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_0abcde_d(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_0abcde_c(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_0abcde_a(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // b is smallest
    return G6_0abcde_b(a, b, c, d, e, sa, sb, sc, sd, se, x);
}
complex<double> G6_a0bcde(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{
    // e is smallest
    if (abs(e) < abs(d) && abs(e) < abs(c) && abs(e) < abs(b) && abs(e) < abs(a))
        return G6_a0bcde_e(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_a0bcde_d(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_a0bcde_c(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_a0bcde_a(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // b is smallest
    return G6_a0bcde_b(a, b, c, d, e, sa, sb, sc, sd, se, x);
}
complex<double> G6_ab0cde(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{
    // e is smallest
    if (abs(e) < abs(d) && abs(e) < abs(c) && abs(e) < abs(b) && abs(e) < abs(a))
        return G6_ab0cde_e(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_ab0cde_d(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_ab0cde_c(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_ab0cde_a(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // b is smallest
    return G6_ab0cde_b(a, b, c, d, e, sa, sb, sc, sd, se, x);
}
complex<double> G6_abc0de(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{
    // e is smallest
    if (abs(e) < abs(d) && abs(e) < abs(c) && abs(e) < abs(b) && abs(e) < abs(a))
        return G6_abc0de_e(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_abc0de_d(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_abc0de_c(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_abc0de_a(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // b is smallest
    return G6_abc0de_b(a, b, c, d, e, sa, sb, sc, sd, se, x);
}
complex<double> G6_abcd0e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x)
{
    // e is smallest
    if (abs(e) < abs(d) && abs(e) < abs(c) && abs(e) < abs(b) && abs(e) < abs(a))
        return G6_abcd0e_e(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_abcd0e_d(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_abcd0e_c(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_abcd0e_a(a, b, c, d, e, sa, sb, sc, sd, se, x);

    // b is smallest
    return G6_abcd0e_b(a, b, c, d, e, sa, sb, sc, sd, se, x);
}
