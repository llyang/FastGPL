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

complex<double> G6_abcdef_f(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, complex<double> f1, int sa, int sb, int sc,
    int sd, int se, int sf, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    complex<double> f = f1 / x1;
    double x = 1.;
    const vector<complex<double>> sy
        = { G({ e / f, d / f, c / f }, 1), G({ a, b, c }, { sa, sb, sc }, x),
              G({ e / f, d / f }, 1), G({ a, b, c, d }, { sa, sb, sc, sd }, x),
              G({ a, b }, { sa, sb }, x), G({ e / f, d / f, c / f, b / f }, 1),
              G({ a, b, c, d, e }, { sa, sb, sc, sd, se }, x), G({ e / f }, 1),
              G({ a }, { sa }, x), G({ e / f, d / f, c / f, b / f, a / f }, 1) };
    complex<double> res { -(sy[3] * G({ 0, e / f }, 1))
        + sy[1] * (sy[0] + G({ 0, e / f, d / f }, 1) + G({ e / f, 0, d / f }, 1))
        + sy[4]
            * (-sy[5] - G({ 0, e / f, d / f, c / f }, 1)
                - G({ e / f, 0, d / f, c / f }, 1) - G({ e / f, d / f, 0, c / f }, 1))
        + sy[8]
            * (sy[9] + G({ 0, e / f, d / f, c / f, b / f }, 1)
                + G({ e / f, 0, d / f, c / f, b / f }, 1)
                + G({ e / f, d / f, 0, c / f, b / f }, 1)
                + G({ e / f, d / f, c / f, 0, b / f }, 1))
        - G({ 0, e / f, d / f, c / f, b / f, a / f }, 1)
        - G({ e / f, 0, d / f, c / f, b / f, a / f }, 1)
        - G({ e / f, d / f, 0, c / f, b / f, a / f }, 1)
        - G({ e / f, d / f, c / f, 0, b / f, a / f }, 1)
        - G({ e / f, d / f, c / f, b / f, 0, a / f }, 1)
        - G({ e / f, d / f, c / f, b / f, a / f, x / f }, 1)
        + sy[6] * (-G({ x / f }, 1) + G({ f }, { sf }, x))
        - sy[5] * G({ 0, a }, { 1, sa }, x)
        + sy[0] * (G({ 0, a, b }, { 1, sa, sb }, x) + G({ a, 0, b }, { sa, 1, sb }, x))
        + sy[2]
            * (-sy[3] - G({ 0, a, b, c }, { 1, sa, sb, sc }, x)
                - G({ a, 0, b, c }, { sa, 1, sb, sc }, x)
                - G({ a, b, 0, c }, { sa, sb, 1, sc }, x))
        + sy[7]
            * (sy[6] + G({ 0, a, b, c, d }, { 1, sa, sb, sc, sd }, x)
                + G({ a, 0, b, c, d }, { sa, 1, sb, sc, sd }, x)
                + G({ a, b, 0, c, d }, { sa, sb, 1, sc, sd }, x)
                + G({ a, b, c, 0, d }, { sa, sb, sc, 1, sd }, x))
        - G({ 0, a, b, c, d, e }, { 1, sa, sb, sc, sd, se }, x)
        - G({ a, 0, b, c, d, e }, { sa, 1, sb, sc, sd, se }, x)
        - G({ a, b, 0, c, d, e }, { sa, sb, 1, sc, sd, se }, x)
        - G({ a, b, c, 0, d, e }, { sa, sb, sc, 1, sd, se }, x)
        - G({ a, b, c, d, 0, e }, { sa, sb, sc, sd, 1, se }, x)
        + (sy[9] - sy[1] * sy[2] + sy[0] * sy[4] + sy[3] * sy[7] - sy[5] * sy[8])
            * Log(f, sf)
        + (-sy[9] + sy[1] * sy[2] - sy[0] * sy[4] - sy[3] * sy[7] + sy[5] * sy[8])
            * Log(-x, sf) };
    return res;
}
complex<double> G6_abcdef_e(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, complex<double> f1, int sa, int sb, int sc,
    int sd, int se, int sf, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    complex<double> f = f1 / x1;
    double x = 1.;
    if (e == f) { // abcdee
        const vector<complex<double>> sy = { G({ a, b, e }, { sa, sb, se }, x),
            G({ d / e, c / e, 1 }, 1), G({ d / e, c / e, b / e }, 1),
            G({ a, b, c, e }, { sa, sb, sc, se }, x),
            G({ a, b, c, d }, { sa, sb, sc, sd }, x), G({ d / e, c / e, b / e, 1 }, 1),
            G({ a, e }, { sa, se }, x), G({ d / e, c / e, b / e, a / e }, 1),
            G({ e }, { se }, x), G({ d / e, c / e, b / e, a / e, 1 }, 1) };
        complex<double> res { -(sy[0] * sy[1]) + sy[0] * sy[2] + sy[5] * sy[6]
            - sy[6] * sy[7] - sy[9] * sy[8] + (sy[3] - sy[4]) * G({ d / e, 1 }, 1)
            + sy[8] * G({ d / e, c / e, b / e, a / e, x / e }, 1)
            + G({ d / e, c / e, b / e, a / e, 0, 1 }, 1)
            - G({ d / e, c / e, b / e, a / e, x / e, 1 }, 1)
            + (sy[9] - G({ d / e, c / e, b / e, 0, 1 }, 1)) * G({ a }, { sa }, x)
            + sy[7] * G({ 0, e }, { 1, se }, x)
            + (-sy[5] + G({ d / e, c / e, 0, 1 }, 1)) * G({ a, b }, { sa, sb }, x)
            - sy[2] * G({ a, 0, e }, { sa, 1, se }, x)
            + (sy[1] - G({ d / e, 0, 1 }, 1)) * G({ a, b, c }, { sa, sb, sc }, x)
            + G({ d / e, c / e }, 1) * (-sy[3] + G({ a, b, 0, e }, { sa, sb, 1, se }, x))
            + G({ d / e }, 1)
                * (-G({ a, b, c, 0, e }, { sa, sb, sc, 1, se }, x)
                    + G({ a, b, c, d, e }, { sa, sb, sc, sd, se }, x))
            + G({ a, b, c, d, 0, e }, { sa, sb, sc, sd, 1, se }, x) - sy[4] * Zeta(2) };
        return res;
    } else { // abcdef
        const vector<complex<double>> sy = { G({ d / e, c / e, b / e }, 1),
            G({ a, b, f }, { sa, sb, sf }, x), G({ 0, a, b }, { 1, sa, sb }, x),
            G({ a, 0, b }, { sa, 1, sb }, x), G({ d / e, c / e, f / e }, 1),
            G({ d / e, f / e, c / e }, 1), G({ f / e, d / e, c / e }, 1),
            G({ a, b, c }, { sa, sb, sc }, x), G({ a, b, c, d }, { sa, sb, sc, sd }, x),
            G({ f / e, d / e }, 1), G({ 0, a, b, c }, { 1, sa, sb, sc }, x),
            G({ a, 0, b, c }, { sa, 1, sb, sc }, x),
            G({ a, b, 0, c }, { sa, sb, 1, sc }, x),
            G({ a, b, c, f }, { sa, sb, sc, sf }, x), G({ d / e, f / e }, 1),
            G({ d / e, c / e, b / e, a / e }, 1), G({ a, f }, { sa, sf }, x),
            G({ 0, a }, { 1, sa }, x), G({ d / e, c / e, b / e, f / e }, 1),
            G({ d / e, c / e, f / e, b / e }, 1), G({ d / e, f / e, c / e, b / e }, 1),
            G({ f / e, d / e, c / e, b / e }, 1), G({ a, b }, { sa, sb }, x),
            G({ f / e }, 1), G({ a, b, c, d, f }, { sa, sb, sc, sd, sf }, x),
            G({ a }, { sa }, x), G({ d / e, c / e, b / e, a / e, f / e }, 1),
            G({ d / e, c / e, b / e, f / e, a / e }, 1),
            G({ d / e, c / e, f / e, b / e, a / e }, 1),
            G({ d / e, f / e, c / e, b / e, a / e }, 1),
            G({ f / e, d / e, c / e, b / e, a / e }, 1) };
        complex<double> res { (sy[10] + sy[11] + sy[12] + sy[13]) * sy[14]
            - sy[15] * sy[16] + sy[16] * sy[18] + sy[17] * sy[18] + sy[0] * sy[1]
            + sy[17] * (sy[19] + sy[20] + sy[21]) - sy[1] * sy[4]
            + (-sy[2] - sy[3]) * sy[4] + sy[2] * (-sy[5] - sy[6])
            + sy[3] * (-sy[5] - sy[6]) + sy[9] * (sy[10] + sy[11] + sy[12] + sy[8])
            + sy[8] * G({ 0, f / e }, 1)
            + sy[7]
                * (-sy[5] - sy[6] - G({ 0, d / e, f / e }, 1) - G({ 0, f / e, d / e }, 1)
                    - G({ d / e, 0, f / e }, 1) - G({ f / e, 0, d / e }, 1))
            + sy[22]
                * (sy[19] + sy[20] + sy[21] + G({ 0, d / e, c / e, f / e }, 1)
                    + G({ 0, d / e, f / e, c / e }, 1) + G({ 0, f / e, d / e, c / e }, 1)
                    + G({ d / e, 0, c / e, f / e }, 1) + G({ d / e, 0, f / e, c / e }, 1)
                    + G({ d / e, c / e, 0, f / e }, 1) + G({ d / e, f / e, 0, c / e }, 1)
                    + G({ f / e, 0, d / e, c / e }, 1) + G({ f / e, d / e, 0, c / e }, 1))
            + sy[25]
                * (-sy[27] - sy[28] - sy[29] - sy[30]
                    - G({ 0, d / e, c / e, b / e, f / e }, 1)
                    - G({ 0, d / e, c / e, f / e, b / e }, 1)
                    - G({ 0, d / e, f / e, c / e, b / e }, 1)
                    - G({ 0, f / e, d / e, c / e, b / e }, 1)
                    - G({ d / e, 0, c / e, b / e, f / e }, 1)
                    - G({ d / e, 0, c / e, f / e, b / e }, 1)
                    - G({ d / e, 0, f / e, c / e, b / e }, 1)
                    - G({ d / e, c / e, 0, b / e, f / e }, 1)
                    - G({ d / e, c / e, 0, f / e, b / e }, 1)
                    - G({ d / e, c / e, b / e, 0, f / e }, 1)
                    - G({ d / e, c / e, f / e, 0, b / e }, 1)
                    - G({ d / e, f / e, 0, c / e, b / e }, 1)
                    - G({ d / e, f / e, c / e, 0, b / e }, 1)
                    - G({ f / e, 0, d / e, c / e, b / e }, 1)
                    - G({ f / e, d / e, 0, c / e, b / e }, 1)
                    - G({ f / e, d / e, c / e, 0, b / e }, 1))
            + G({ 0, d / e, c / e, b / e, a / e, f / e }, 1)
            + G({ 0, d / e, c / e, b / e, f / e, a / e }, 1)
            + G({ 0, d / e, c / e, f / e, b / e, a / e }, 1)
            + G({ 0, d / e, f / e, c / e, b / e, a / e }, 1)
            + G({ 0, f / e, d / e, c / e, b / e, a / e }, 1)
            + G({ d / e, 0, c / e, b / e, a / e, f / e }, 1)
            + G({ d / e, 0, c / e, b / e, f / e, a / e }, 1)
            + G({ d / e, 0, c / e, f / e, b / e, a / e }, 1)
            + G({ d / e, 0, f / e, c / e, b / e, a / e }, 1)
            + G({ d / e, c / e, 0, b / e, a / e, f / e }, 1)
            + G({ d / e, c / e, 0, b / e, f / e, a / e }, 1)
            + G({ d / e, c / e, 0, f / e, b / e, a / e }, 1)
            + G({ d / e, c / e, b / e, 0, a / e, f / e }, 1)
            + G({ d / e, c / e, b / e, 0, f / e, a / e }, 1)
            + G({ d / e, c / e, b / e, a / e, 0, f / e }, 1)
            + G({ d / e, c / e, b / e, a / e, f / e, x / e }, 1)
            + G({ d / e, c / e, b / e, f / e, 0, a / e }, 1)
            + G({ d / e, c / e, b / e, f / e, a / e, x / e }, 1)
            + G({ d / e, c / e, f / e, 0, b / e, a / e }, 1)
            + G({ d / e, c / e, f / e, b / e, 0, a / e }, 1)
            + G({ d / e, c / e, f / e, b / e, a / e, x / e }, 1)
            + G({ d / e, f / e, 0, c / e, b / e, a / e }, 1)
            + G({ d / e, f / e, c / e, 0, b / e, a / e }, 1)
            + G({ d / e, f / e, c / e, b / e, 0, a / e }, 1)
            + G({ d / e, f / e, c / e, b / e, a / e, x / e }, 1)
            + G({ f / e, 0, d / e, c / e, b / e, a / e }, 1)
            + G({ f / e, d / e, 0, c / e, b / e, a / e }, 1)
            + G({ f / e, d / e, c / e, 0, b / e, a / e }, 1)
            + G({ f / e, d / e, c / e, b / e, 0, a / e }, 1)
            + G({ f / e, d / e, c / e, b / e, a / e, x / e }, 1)
            + sy[15] * G({ 0, f }, { 1, sf }, x)
            - sy[0] * G({ a, 0, f }, { sa, 1, sf }, x)
            + G({ d / e, c / e }, 1) * (-sy[13] + G({ a, b, 0, f }, { sa, sb, 1, sf }, x))
            + sy[23]
                * (-sy[24] - G({ 0, a, b, c, d }, { 1, sa, sb, sc, sd }, x)
                    - G({ a, 0, b, c, d }, { sa, 1, sb, sc, sd }, x)
                    - G({ a, b, 0, c, d }, { sa, sb, 1, sc, sd }, x)
                    - G({ a, b, c, 0, d }, { sa, sb, sc, 1, sd }, x))
            + G({ d / e }, 1) * (sy[24] - G({ a, b, c, 0, f }, { sa, sb, sc, 1, sf }, x))
            + G({ a, b, c, d, 0, f }, { sa, sb, sc, sd, 1, sf }, x)
            + ((sy[18] + sy[19] + sy[20] + sy[21]) * sy[25] - sy[26] - sy[27] - sy[28]
                  - sy[29] - sy[30] + sy[22] * (-sy[4] - sy[5] - sy[6]) + sy[9] * sy[7]
                  + sy[14] * sy[7] - sy[23] * sy[8])
                * Log(e, se)
            + ((-sy[18] - sy[19] - sy[20] - sy[21]) * sy[25] + sy[26] + sy[27] + sy[28]
                  + sy[29] + sy[30] + sy[22] * (sy[4] + sy[5] + sy[6]) - sy[9] * sy[7]
                  - sy[14] * sy[7] + sy[23] * sy[8])
                * Log(-x, se) };
        if (f != x) {
            res += (-sy[26] + G({ d / e, c / e, b / e, a / e, x / e }, 1))
                * G({ f }, { sf }, x);
        }
        return res;
    }
}
complex<double> G6_abcdef_d(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, complex<double> f1, int sa, int sb, int sc,
    int sd, int se, int sf, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    complex<double> f = f1 / x1;
    double x = 1.;
    if (d == e) {
        if (d == f) { // abcddd
            const vector<complex<double>> sy = { G({ a, b, d }, { sa, sb, sd }, x),
                G({ c / d, 1, 1 }, 1), G({ c / d, b / d, 1 }, 1),
                G({ a, d, d }, { sa, sd, sd }, x), G({ c / d, b / d, a / d }, 1),
                G({ a, b, d, d }, { sa, sb, sd, sd }, x), G({ c / d, b / d, 1, 1 }, 1),
                G({ a, d }, { sa, sd }, x), G({ c / d, b / d, a / d, 1 }, 1),
                G({ d, d }, { sd, sd }, x), G({ d }, { sd }, x),
                G({ c / d, b / d, a / d, 1, 1 }, 1) };
            complex<double> res { -(sy[10] * sy[11]) - sy[0] * sy[1] + sy[0] * sy[2]
                - sy[2] * sy[3] + sy[3] * sy[4] + sy[6] * sy[7] + sy[9] * sy[8]
                - sy[7] * sy[8] - sy[9] * G({ c / d, b / d, a / d, x / d }, 1)
                + sy[10] * G({ c / d, b / d, a / d, x / d, 1 }, 1)
                + G({ c / d, b / d, a / d, 0, 1, 1 }, 1)
                - G({ c / d, b / d, a / d, x / d, 1, 1 }, 1)
                + (sy[11] - G({ c / d, b / d, 0, 1, 1 }, 1)) * G({ a }, { sa }, x)
                + (-sy[6] + G({ c / d, 0, 1, 1 }, 1)) * G({ a, b }, { sa, sb }, x)
                - sy[4] * G({ 0, d, d }, { 1, sd, sd }, x)
                + G({ c / d, b / d }, 1)
                    * (-sy[5] + G({ a, 0, d, d }, { sa, 1, sd, sd }, x))
                + G({ c / d, 1 }, 1) * (sy[5] - G({ a, b, c, d }, { sa, sb, sc, sd }, x))
                + G({ c / d }, 1)
                    * (-G({ a, b, 0, d, d }, { sa, sb, 1, sd, sd }, x)
                        + G({ a, b, c, d, d }, { sa, sb, sc, sd, sd }, x))
                + G({ a, b, c, 0, d, d }, { sa, sb, sc, 1, sd, sd }, x)
                + G({ a, b, c }, { sa, sb, sc }, x) * (sy[1] - Zeta(3)) };
            return res;
        } else { // abcddf
            const vector<complex<double>> sy = { G({ a, d, f }, { sa, sd, sf }, x),
                G({ c / d, b / d, 1 }, 1), G({ c / d, b / d, a / d }, 1),
                G({ c / d, f / d, 1 }, 1), G({ 0, a, b }, { 1, sa, sb }, x),
                G({ a, 0, b }, { sa, 1, sb }, x), G({ f / d, 1, c / d }, 1),
                G({ f / d, c / d, 1 }, 1), G({ a, b, c }, { sa, sb, sc }, x),
                G({ f / d, 1 }, 1), G({ a, b, c, f }, { sa, sb, sc, sf }, x),
                G({ a, b, d, f }, { sa, sb, sd, sf }, x), G({ d, f }, { sd, sf }, x),
                G({ c / d, b / d, a / d, 1 }, 1), G({ 0, a }, { 1, sa }, x),
                G({ c / d, b / d, f / d, 1 }, 1), G({ c / d, f / d, 1, b / d }, 1),
                G({ c / d, f / d, b / d, 1 }, 1), G({ f / d, 1, c / d, b / d }, 1),
                G({ f / d, c / d, 1, b / d }, 1), G({ f / d, c / d, b / d, 1 }, 1),
                G({ a, b }, { sa, sb }, x), G({ a }, { sa }, x),
                G({ c / d, b / d, a / d, f / d, 1 }, 1),
                G({ c / d, b / d, f / d, 1, a / d }, 1),
                G({ c / d, b / d, f / d, a / d, 1 }, 1),
                G({ c / d, f / d, 1, b / d, a / d }, 1),
                G({ c / d, f / d, b / d, 1, a / d }, 1),
                G({ c / d, f / d, b / d, a / d, 1 }, 1),
                G({ f / d, 1, c / d, b / d, a / d }, 1),
                G({ f / d, c / d, 1, b / d, a / d }, 1),
                G({ f / d, c / d, b / d, 1, a / d }, 1),
                G({ f / d, c / d, b / d, a / d, 1 }, 1) };
            complex<double> res { sy[12] * sy[13] + sy[14] * sy[15] - sy[0] * sy[1]
                + sy[14] * (sy[16] + sy[17] + sy[18] + sy[19] + sy[20]) + sy[0] * sy[2]
                + sy[3] * (-sy[4] - sy[5]) + sy[4] * (-sy[6] - sy[7])
                + sy[5] * (-sy[6] - sy[7]) + (-sy[10] + sy[11]) * G({ c / d, 1 }, 1)
                + sy[8] * (-sy[6] - sy[7] - G({ 0, f / d, 1 }, 1))
                - sy[12] * G({ c / d, b / d, a / d, x / d }, 1)
                + sy[21]
                    * (sy[16] + sy[17] + sy[18] + sy[19] + sy[20]
                        + G({ 0, c / d, f / d, 1 }, 1) + G({ 0, f / d, 1, c / d }, 1)
                        + G({ 0, f / d, c / d, 1 }, 1) + G({ c / d, 0, f / d, 1 }, 1)
                        + G({ f / d, 0, 1, c / d }, 1) + G({ f / d, 0, c / d, 1 }, 1)
                        + G({ f / d, 1, 0, c / d }, 1))
                + sy[22]
                    * (-sy[24] - sy[25] - sy[26] - sy[27] - sy[28] - sy[29] - sy[30]
                        - sy[31] - sy[32] - G({ 0, c / d, b / d, f / d, 1 }, 1)
                        - G({ 0, c / d, f / d, 1, b / d }, 1)
                        - G({ 0, c / d, f / d, b / d, 1 }, 1)
                        - G({ 0, f / d, 1, c / d, b / d }, 1)
                        - G({ 0, f / d, c / d, 1, b / d }, 1)
                        - G({ 0, f / d, c / d, b / d, 1 }, 1)
                        - G({ c / d, 0, b / d, f / d, 1 }, 1)
                        - G({ c / d, 0, f / d, 1, b / d }, 1)
                        - G({ c / d, 0, f / d, b / d, 1 }, 1)
                        - G({ c / d, b / d, 0, f / d, 1 }, 1)
                        - G({ c / d, f / d, 0, 1, b / d }, 1)
                        - G({ c / d, f / d, 0, b / d, 1 }, 1)
                        - G({ c / d, f / d, 1, 0, b / d }, 1)
                        - G({ f / d, 0, 1, c / d, b / d }, 1)
                        - G({ f / d, 0, c / d, 1, b / d }, 1)
                        - G({ f / d, 0, c / d, b / d, 1 }, 1)
                        - G({ f / d, 1, 0, c / d, b / d }, 1)
                        - G({ f / d, 1, c / d, 0, b / d }, 1)
                        - G({ f / d, c / d, 0, 1, b / d }, 1)
                        - G({ f / d, c / d, 0, b / d, 1 }, 1)
                        - G({ f / d, c / d, 1, 0, b / d }, 1))
                + G({ 0, c / d, b / d, a / d, f / d, 1 }, 1)
                + G({ 0, c / d, b / d, f / d, 1, a / d }, 1)
                + G({ 0, c / d, b / d, f / d, a / d, 1 }, 1)
                + G({ 0, c / d, f / d, 1, b / d, a / d }, 1)
                + G({ 0, c / d, f / d, b / d, 1, a / d }, 1)
                + G({ 0, c / d, f / d, b / d, a / d, 1 }, 1)
                + G({ 0, f / d, 1, c / d, b / d, a / d }, 1)
                + G({ 0, f / d, c / d, 1, b / d, a / d }, 1)
                + G({ 0, f / d, c / d, b / d, 1, a / d }, 1)
                + G({ 0, f / d, c / d, b / d, a / d, 1 }, 1)
                + G({ c / d, 0, b / d, a / d, f / d, 1 }, 1)
                + G({ c / d, 0, b / d, f / d, 1, a / d }, 1)
                + G({ c / d, 0, b / d, f / d, a / d, 1 }, 1)
                + G({ c / d, 0, f / d, 1, b / d, a / d }, 1)
                + G({ c / d, 0, f / d, b / d, 1, a / d }, 1)
                + G({ c / d, 0, f / d, b / d, a / d, 1 }, 1)
                + G({ c / d, b / d, 0, a / d, f / d, 1 }, 1)
                + G({ c / d, b / d, 0, f / d, 1, a / d }, 1)
                + G({ c / d, b / d, 0, f / d, a / d, 1 }, 1)
                + G({ c / d, b / d, a / d, 0, f / d, 1 }, 1)
                + G({ c / d, b / d, a / d, f / d, 1, x / d }, 1)
                + G({ c / d, b / d, a / d, f / d, x / d, 1 }, 1)
                + G({ c / d, b / d, f / d, 0, 1, a / d }, 1)
                + G({ c / d, b / d, f / d, 0, a / d, 1 }, 1)
                + G({ c / d, b / d, f / d, 1, 0, a / d }, 1)
                + G({ c / d, b / d, f / d, 1, a / d, x / d }, 1)
                + G({ c / d, b / d, f / d, a / d, 1, x / d }, 1)
                + G({ c / d, b / d, f / d, a / d, x / d, 1 }, 1)
                + G({ c / d, f / d, 0, 1, b / d, a / d }, 1)
                + G({ c / d, f / d, 0, b / d, 1, a / d }, 1)
                + G({ c / d, f / d, 0, b / d, a / d, 1 }, 1)
                + G({ c / d, f / d, 1, 0, b / d, a / d }, 1)
                + G({ c / d, f / d, 1, b / d, 0, a / d }, 1)
                + G({ c / d, f / d, 1, b / d, a / d, x / d }, 1)
                + G({ c / d, f / d, b / d, 0, 1, a / d }, 1)
                + G({ c / d, f / d, b / d, 0, a / d, 1 }, 1)
                + G({ c / d, f / d, b / d, 1, 0, a / d }, 1)
                + G({ c / d, f / d, b / d, 1, a / d, x / d }, 1)
                + G({ c / d, f / d, b / d, a / d, 1, x / d }, 1)
                + G({ c / d, f / d, b / d, a / d, x / d, 1 }, 1)
                + G({ f / d, 0, 1, c / d, b / d, a / d }, 1)
                + G({ f / d, 0, c / d, 1, b / d, a / d }, 1)
                + G({ f / d, 0, c / d, b / d, 1, a / d }, 1)
                + G({ f / d, 0, c / d, b / d, a / d, 1 }, 1)
                + G({ f / d, 1, 0, c / d, b / d, a / d }, 1)
                + G({ f / d, 1, c / d, 0, b / d, a / d }, 1)
                + G({ f / d, 1, c / d, b / d, 0, a / d }, 1)
                + G({ f / d, 1, c / d, b / d, a / d, x / d }, 1)
                + G({ f / d, c / d, 0, 1, b / d, a / d }, 1)
                + G({ f / d, c / d, 0, b / d, 1, a / d }, 1)
                + G({ f / d, c / d, 0, b / d, a / d, 1 }, 1)
                + G({ f / d, c / d, 1, 0, b / d, a / d }, 1)
                + G({ f / d, c / d, 1, b / d, 0, a / d }, 1)
                + G({ f / d, c / d, 1, b / d, a / d, x / d }, 1)
                + G({ f / d, c / d, b / d, 0, 1, a / d }, 1)
                + G({ f / d, c / d, b / d, 0, a / d, 1 }, 1)
                + G({ f / d, c / d, b / d, 1, 0, a / d }, 1)
                + G({ f / d, c / d, b / d, 1, a / d, x / d }, 1)
                + G({ f / d, c / d, b / d, a / d, 1, x / d }, 1)
                + G({ f / d, c / d, b / d, a / d, x / d, 1 }, 1)
                + (-sy[13] + sy[15]) * G({ a, f }, { sa, sf }, x)
                - sy[2] * G({ 0, d, f }, { 1, sd, sf }, x)
                + (sy[1] - sy[3]) * G({ a, b, f }, { sa, sb, sf }, x)
                + G({ c / d, b / d }, 1)
                    * (-sy[11] + G({ a, 0, d, f }, { sa, 1, sd, sf }, x))
                + sy[9]
                    * (sy[10] + G({ 0, a, b, c }, { 1, sa, sb, sc }, x)
                        + G({ a, 0, b, c }, { sa, 1, sb, sc }, x)
                        + G({ a, b, 0, c }, { sa, sb, 1, sc }, x))
                + G({ c / d }, 1)
                    * (-G({ a, b, 0, d, f }, { sa, sb, 1, sd, sf }, x)
                        + G({ a, b, c, d, f }, { sa, sb, sc, sd, sf }, x))
                + G({ a, b, c, 0, d, f }, { sa, sb, sc, 1, sd, sf }, x)
                + ((sy[15] + sy[16] + sy[17] + sy[18] + sy[19] + sy[20]) * sy[22] - sy[23]
                      - sy[24] - sy[25] - sy[26] - sy[27] - sy[28] - sy[29] - sy[30]
                      - sy[31] - sy[32] + sy[21] * (-sy[3] - sy[6] - sy[7])
                      + sy[9] * sy[8])
                    * Log(d, sd)
                + ((-sy[15] - sy[16] - sy[17] - sy[18] - sy[19] - sy[20]) * sy[22]
                      + sy[23] + sy[24] + sy[25] + sy[26] + sy[27] + sy[28] + sy[29]
                      + sy[30] + sy[31] + sy[32] + sy[21] * (sy[3] + sy[6] + sy[7])
                      - sy[9] * sy[8])
                    * Log(-x, sd) };
            if (f != x) {
                res += (-sy[23] + G({ c / d, b / d, a / d, x / d, 1 }, 1))
                    * G({ f }, { sf }, x);
            }
            return res;
        }
    } else { // abcdef
        const vector<complex<double>> sy = { G({ c / d, b / d, a / d }, 1),
            G({ a, e, f }, { sa, se, sf }, x), G({ c / d, b / d, e / d }, 1),
            G({ a, b, f }, { sa, sb, sf }, x), G({ c / d, e / d, b / d }, 1),
            G({ 0, a, b }, { 1, sa, sb }, x), G({ a, 0, b }, { sa, 1, sb }, x),
            G({ c / d, e / d, f / d }, 1), G({ e / d, c / d, b / d }, 1),
            G({ e / d, c / d, f / d }, 1), G({ a, b, c }, { sa, sb, sc }, x),
            G({ e / d, f / d, c / d }, 1), G({ e / d, c / d }, 1),
            G({ a, b, 0, f }, { sa, sb, 1, sf }, x), G({ e / d, f / d }, 1),
            G({ a, b, c, f }, { sa, sb, sc, sf }, x),
            G({ a, b, e, f }, { sa, sb, se, sf }, x), G({ a, f }, { sa, sf }, x),
            G({ c / d, b / d, e / d, a / d }, 1), G({ 0, a }, { 1, sa }, x),
            G({ c / d, b / d, e / d, f / d }, 1), G({ c / d, b / d, a / d, e / d }, 1),
            G({ c / d, e / d, b / d, a / d }, 1), G({ e / d, c / d, b / d, a / d }, 1),
            G({ c / d, e / d, b / d, f / d }, 1), G({ e / d, c / d, b / d, f / d }, 1),
            G({ c / d, e / d, f / d, b / d }, 1), G({ e / d, c / d, f / d, b / d }, 1),
            G({ e / d, f / d, c / d, b / d }, 1), G({ a, b }, { sa, sb }, x),
            G({ a, b, c, e, f }, { sa, sb, sc, se, sf }, x), G({ a }, { sa }, x),
            G({ c / d, b / d, a / d, e / d, f / d }, 1),
            G({ c / d, b / d, e / d, a / d, f / d }, 1),
            G({ c / d, b / d, e / d, f / d, a / d }, 1),
            G({ c / d, e / d, b / d, a / d, f / d }, 1),
            G({ c / d, e / d, b / d, f / d, a / d }, 1),
            G({ c / d, e / d, f / d, b / d, a / d }, 1),
            G({ e / d, c / d, b / d, a / d, f / d }, 1),
            G({ e / d, c / d, b / d, f / d, a / d }, 1),
            G({ e / d, c / d, f / d, b / d, a / d }, 1),
            G({ e / d, f / d, c / d, b / d, a / d }, 1) };
        complex<double> res { -(sy[12] * sy[13]) + sy[12] * sy[15] + sy[17] * sy[18]
            + sy[0] * sy[1] - sy[19] * sy[20]
            + sy[17] * (-sy[20] + sy[22] + sy[23] - sy[24] - sy[25])
            + sy[19] * (-sy[24] - sy[25] - sy[26] - sy[27] - sy[28]) - sy[1] * sy[2]
            - sy[3] * sy[4] + (sy[9] + sy[11]) * sy[5] + (sy[9] + sy[11]) * sy[6]
            + (sy[5] + sy[6]) * sy[7] + sy[3] * (sy[9] + sy[7] - sy[8])
            + (-sy[13] + sy[16]) * G({ c / d, e / d }, 1)
            + sy[10] * (sy[11] + G({ 0, e / d, f / d }, 1) + G({ e / d, 0, f / d }, 1))
            + sy[29]
                * (-sy[26] - sy[27] - sy[28] - G({ 0, c / d, e / d, f / d }, 1)
                    - G({ 0, e / d, c / d, f / d }, 1) - G({ 0, e / d, f / d, c / d }, 1)
                    - G({ c / d, 0, e / d, f / d }, 1) - G({ c / d, e / d, 0, f / d }, 1)
                    - G({ e / d, 0, c / d, f / d }, 1) - G({ e / d, 0, f / d, c / d }, 1)
                    - G({ e / d, c / d, 0, f / d }, 1) - G({ e / d, f / d, 0, c / d }, 1))
            + sy[31]
                * (sy[34] + sy[36] + sy[37] + sy[39] + sy[40] + sy[41]
                    + G({ 0, c / d, b / d, e / d, f / d }, 1)
                    + G({ 0, c / d, e / d, b / d, f / d }, 1)
                    + G({ 0, c / d, e / d, f / d, b / d }, 1)
                    + G({ 0, e / d, c / d, b / d, f / d }, 1)
                    + G({ 0, e / d, c / d, f / d, b / d }, 1)
                    + G({ 0, e / d, f / d, c / d, b / d }, 1)
                    + G({ c / d, 0, b / d, e / d, f / d }, 1)
                    + G({ c / d, 0, e / d, b / d, f / d }, 1)
                    + G({ c / d, 0, e / d, f / d, b / d }, 1)
                    + G({ c / d, b / d, 0, e / d, f / d }, 1)
                    + G({ c / d, b / d, e / d, 0, f / d }, 1)
                    + G({ c / d, e / d, 0, b / d, f / d }, 1)
                    + G({ c / d, e / d, 0, f / d, b / d }, 1)
                    + G({ c / d, e / d, b / d, 0, f / d }, 1)
                    + G({ c / d, e / d, f / d, 0, b / d }, 1)
                    + G({ e / d, 0, c / d, b / d, f / d }, 1)
                    + G({ e / d, 0, c / d, f / d, b / d }, 1)
                    + G({ e / d, 0, f / d, c / d, b / d }, 1)
                    + G({ e / d, c / d, 0, b / d, f / d }, 1)
                    + G({ e / d, c / d, 0, f / d, b / d }, 1)
                    + G({ e / d, c / d, b / d, 0, f / d }, 1)
                    + G({ e / d, c / d, f / d, 0, b / d }, 1)
                    + G({ e / d, f / d, 0, c / d, b / d }, 1)
                    + G({ e / d, f / d, c / d, 0, b / d }, 1))
            - G({ 0, c / d, b / d, a / d, e / d, f / d }, 1)
            - G({ 0, c / d, b / d, e / d, a / d, f / d }, 1)
            - G({ 0, c / d, b / d, e / d, f / d, a / d }, 1)
            - G({ 0, c / d, e / d, b / d, a / d, f / d }, 1)
            - G({ 0, c / d, e / d, b / d, f / d, a / d }, 1)
            - G({ 0, c / d, e / d, f / d, b / d, a / d }, 1)
            - G({ 0, e / d, c / d, b / d, a / d, f / d }, 1)
            - G({ 0, e / d, c / d, b / d, f / d, a / d }, 1)
            - G({ 0, e / d, c / d, f / d, b / d, a / d }, 1)
            - G({ 0, e / d, f / d, c / d, b / d, a / d }, 1)
            - G({ c / d, 0, b / d, a / d, e / d, f / d }, 1)
            - G({ c / d, 0, b / d, e / d, a / d, f / d }, 1)
            - G({ c / d, 0, b / d, e / d, f / d, a / d }, 1)
            - G({ c / d, 0, e / d, b / d, a / d, f / d }, 1)
            - G({ c / d, 0, e / d, b / d, f / d, a / d }, 1)
            - G({ c / d, 0, e / d, f / d, b / d, a / d }, 1)
            - G({ c / d, b / d, 0, a / d, e / d, f / d }, 1)
            - G({ c / d, b / d, 0, e / d, a / d, f / d }, 1)
            - G({ c / d, b / d, 0, e / d, f / d, a / d }, 1)
            - G({ c / d, b / d, a / d, 0, e / d, f / d }, 1)
            - G({ c / d, b / d, a / d, e / d, 0, f / d }, 1)
            - G({ c / d, b / d, a / d, e / d, f / d, x / d }, 1)
            - G({ c / d, b / d, e / d, 0, a / d, f / d }, 1)
            - G({ c / d, b / d, e / d, 0, f / d, a / d }, 1)
            - G({ c / d, b / d, e / d, a / d, 0, f / d }, 1)
            - G({ c / d, b / d, e / d, a / d, f / d, x / d }, 1)
            - G({ c / d, b / d, e / d, f / d, 0, a / d }, 1)
            - G({ c / d, b / d, e / d, f / d, a / d, x / d }, 1)
            - G({ c / d, e / d, 0, b / d, a / d, f / d }, 1)
            - G({ c / d, e / d, 0, b / d, f / d, a / d }, 1)
            - G({ c / d, e / d, 0, f / d, b / d, a / d }, 1)
            - G({ c / d, e / d, b / d, 0, a / d, f / d }, 1)
            - G({ c / d, e / d, b / d, 0, f / d, a / d }, 1)
            - G({ c / d, e / d, b / d, a / d, 0, f / d }, 1)
            - G({ c / d, e / d, b / d, a / d, f / d, x / d }, 1)
            - G({ c / d, e / d, b / d, f / d, 0, a / d }, 1)
            - G({ c / d, e / d, b / d, f / d, a / d, x / d }, 1)
            - G({ c / d, e / d, f / d, 0, b / d, a / d }, 1)
            - G({ c / d, e / d, f / d, b / d, 0, a / d }, 1)
            - G({ c / d, e / d, f / d, b / d, a / d, x / d }, 1)
            - G({ e / d, 0, c / d, b / d, a / d, f / d }, 1)
            - G({ e / d, 0, c / d, b / d, f / d, a / d }, 1)
            - G({ e / d, 0, c / d, f / d, b / d, a / d }, 1)
            - G({ e / d, 0, f / d, c / d, b / d, a / d }, 1)
            - G({ e / d, c / d, 0, b / d, a / d, f / d }, 1)
            - G({ e / d, c / d, 0, b / d, f / d, a / d }, 1)
            - G({ e / d, c / d, 0, f / d, b / d, a / d }, 1)
            - G({ e / d, c / d, b / d, 0, a / d, f / d }, 1)
            - G({ e / d, c / d, b / d, 0, f / d, a / d }, 1)
            - G({ e / d, c / d, b / d, a / d, 0, f / d }, 1)
            - G({ e / d, c / d, b / d, a / d, f / d, x / d }, 1)
            - G({ e / d, c / d, b / d, f / d, 0, a / d }, 1)
            - G({ e / d, c / d, b / d, f / d, a / d, x / d }, 1)
            - G({ e / d, c / d, f / d, 0, b / d, a / d }, 1)
            - G({ e / d, c / d, f / d, b / d, 0, a / d }, 1)
            - G({ e / d, c / d, f / d, b / d, a / d, x / d }, 1)
            - G({ e / d, f / d, 0, c / d, b / d, a / d }, 1)
            - G({ e / d, f / d, c / d, 0, b / d, a / d }, 1)
            - G({ e / d, f / d, c / d, b / d, 0, a / d }, 1)
            - G({ e / d, f / d, c / d, b / d, a / d, x / d }, 1)
            + (-sy[18] - sy[21] - sy[22] - sy[23]) * G({ 0, f }, { 1, sf }, x)
            - sy[0] * G({ 0, e, f }, { 1, se, sf }, x)
            + (sy[2] + sy[4] + sy[8]) * G({ a, 0, f }, { sa, 1, sf }, x)
            + G({ c / d, b / d }, 1) * (-sy[16] + G({ a, 0, e, f }, { sa, 1, se, sf }, x))
            + sy[14]
                * (-sy[15] - G({ 0, a, b, c }, { 1, sa, sb, sc }, x)
                    - G({ a, 0, b, c }, { sa, 1, sb, sc }, x)
                    - G({ a, b, 0, c }, { sa, sb, 1, sc }, x))
            + G({ c / d }, 1) * (sy[30] - G({ a, b, 0, e, f }, { sa, sb, 1, se, sf }, x))
            + G({ e / d }, 1) * (-sy[30] + G({ a, b, c, 0, f }, { sa, sb, sc, 1, sf }, x))
            + G({ a, b, c, 0, e, f }, { sa, sb, sc, 1, se, sf }, x)
            + (-(sy[10] * sy[14])
                  + (-sy[20] - sy[24] - sy[25] - sy[26] - sy[27] - sy[28]) * sy[31]
                  + sy[32] + sy[33] + sy[34] + sy[35] + sy[36] + sy[37] + sy[38] + sy[39]
                  + sy[40] + sy[41] + sy[29] * (sy[9] + sy[11] + sy[7]))
                * Log(d, sd)
            + (sy[10] * sy[14]
                  + (sy[20] + sy[24] + sy[25] + sy[26] + sy[27] + sy[28]) * sy[31]
                  - sy[32] - sy[33] - sy[34] - sy[35] - sy[36] - sy[37] - sy[38] - sy[39]
                  - sy[40] - sy[41] + sy[29] * (-sy[9] - sy[11] - sy[7]))
                * Log(-x, sd) };
        if (e != x) {
            res += (sy[21] - G({ c / d, b / d, a / d, x / d }, 1))
                * G({ e, f }, { se, sf }, x);
        }
        if (f != x) {
            res += (sy[32] + sy[33] + sy[35] + sy[38]
                       - G({ c / d, b / d, a / d, e / d, x / d }, 1)
                       - G({ c / d, b / d, e / d, a / d, x / d }, 1)
                       - G({ c / d, e / d, b / d, a / d, x / d }, 1)
                       - G({ e / d, c / d, b / d, a / d, x / d }, 1))
                * G({ f }, { sf }, x);
        }
        return res;
    }
}
complex<double> G6_abcdef_c(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, complex<double> f1, int sa, int sb, int sc,
    int sd, int se, int sf, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    complex<double> f = f1 / x1;
    double x = 1.;
    if (c == d) {
        if (c == e) {
            if (c == f) { // abcccc
                const vector<complex<double>> sy
                    = { G({ b / c, 1, 1 }, 1), G({ a, c, c }, { sa, sc, sc }, x),
                          G({ b / c, a / c, 1 }, 1), G({ c, c, c }, { sc, sc, sc }, x),
                          G({ a, c, c, c }, { sa, sc, sc, sc }, x),
                          G({ a, c }, { sa, sc }, x), G({ b / c, 1, 1, 1 }, 1),
                          G({ b / c, a / c, 1, 1 }, 1), G({ c, c }, { sc, sc }, x),
                          G({ c }, { sc }, x), G({ b / c, a / c, 1, 1, 1 }, 1) };
                complex<double> res { -(sy[9] * sy[10]) - sy[0] * sy[1] + sy[1] * sy[2]
                    - sy[2] * sy[3] + sy[5] * sy[6] - sy[5] * sy[7] + sy[7] * sy[8]
                    + sy[3] * G({ b / c, a / c, x / c }, 1)
                    - sy[8] * G({ b / c, a / c, x / c, 1 }, 1)
                    + sy[9] * G({ b / c, a / c, x / c, 1, 1 }, 1)
                    + G({ b / c, a / c, 0, 1, 1, 1 }, 1)
                    - G({ b / c, a / c, x / c, 1, 1, 1 }, 1)
                    + (sy[10] - G({ b / c, 0, 1, 1, 1 }, 1)) * G({ a }, { sa }, x)
                    + sy[0] * G({ a, b, c }, { sa, sb, sc }, x)
                    + G({ b / c, a / c }, 1)
                        * (-sy[4] + G({ 0, c, c, c }, { 1, sc, sc, sc }, x))
                    + G({ b / c, 1 }, 1)
                        * (sy[4] - G({ a, b, c, c }, { sa, sb, sc, sc }, x))
                    + G({ b / c }, 1)
                        * (-G({ a, 0, c, c, c }, { sa, 1, sc, sc, sc }, x)
                            + G({ a, b, c, c, c }, { sa, sb, sc, sc, sc }, x))
                    + G({ a, b, 0, c, c, c }, { sa, sb, 1, sc, sc, sc }, x)
                    + G({ a, b }, { sa, sb }, x) * (-sy[6] - Zeta(4)) };
                return res;
            } else { // abcccf
                const vector<complex<double>> sy
                    = { G({ a, c, f }, { sa, sc, sf }, x), G({ b / c, 1, 1 }, 1),
                          G({ b / c, a / c, 1 }, 1), G({ c, c, f }, { sc, sc, sf }, x),
                          G({ f / c, 1, 1 }, 1), G({ a, c, c, f }, { sa, sc, sc, sf }, x),
                          G({ c, f }, { sc, sf }, x), G({ b / c, a / c, 1, 1 }, 1),
                          G({ 0, a }, { 1, sa }, x), G({ b / c, f / c, 1, 1 }, 1),
                          G({ f / c, 1, 1, b / c }, 1), G({ f / c, 1, b / c, 1 }, 1),
                          G({ f / c, b / c, 1, 1 }, 1), G({ a, b }, { sa, sb }, x),
                          G({ a }, { sa }, x), G({ b / c, a / c, f / c, 1, 1 }, 1),
                          G({ b / c, f / c, 1, 1, a / c }, 1),
                          G({ b / c, f / c, 1, a / c, 1 }, 1),
                          G({ b / c, f / c, a / c, 1, 1 }, 1),
                          G({ f / c, 1, 1, b / c, a / c }, 1),
                          G({ f / c, 1, b / c, 1, a / c }, 1),
                          G({ f / c, 1, b / c, a / c, 1 }, 1),
                          G({ f / c, b / c, 1, 1, a / c }, 1),
                          G({ f / c, b / c, 1, a / c, 1 }, 1),
                          G({ f / c, b / c, a / c, 1, 1 }, 1) };
                complex<double> res { -(sy[0] * sy[1]) + sy[0] * sy[2] - sy[2] * sy[3]
                    + sy[6] * sy[7] + sy[9] * sy[8] + (sy[10] + sy[11] + sy[12]) * sy[8]
                    + sy[3] * G({ b / c, a / c, x / c }, 1)
                    + sy[13] * (sy[10] + sy[11] + sy[12] + G({ 0, f / c, 1, 1 }, 1))
                    - sy[6] * G({ b / c, a / c, x / c, 1 }, 1)
                    + sy[14]
                        * (-sy[16] - sy[17] - sy[18] - sy[19] - sy[20] - sy[21] - sy[22]
                            - sy[23] - sy[24] - G({ 0, b / c, f / c, 1, 1 }, 1)
                            - G({ 0, f / c, 1, 1, b / c }, 1)
                            - G({ 0, f / c, 1, b / c, 1 }, 1)
                            - G({ 0, f / c, b / c, 1, 1 }, 1)
                            - G({ b / c, 0, f / c, 1, 1 }, 1)
                            - G({ f / c, 0, 1, 1, b / c }, 1)
                            - G({ f / c, 0, 1, b / c, 1 }, 1)
                            - G({ f / c, 0, b / c, 1, 1 }, 1)
                            - G({ f / c, 1, 0, 1, b / c }, 1)
                            - G({ f / c, 1, 0, b / c, 1 }, 1)
                            - G({ f / c, 1, 1, 0, b / c }, 1))
                    + G({ 0, b / c, a / c, f / c, 1, 1 }, 1)
                    + G({ 0, b / c, f / c, 1, 1, a / c }, 1)
                    + G({ 0, b / c, f / c, 1, a / c, 1 }, 1)
                    + G({ 0, b / c, f / c, a / c, 1, 1 }, 1)
                    + G({ 0, f / c, 1, 1, b / c, a / c }, 1)
                    + G({ 0, f / c, 1, b / c, 1, a / c }, 1)
                    + G({ 0, f / c, 1, b / c, a / c, 1 }, 1)
                    + G({ 0, f / c, b / c, 1, 1, a / c }, 1)
                    + G({ 0, f / c, b / c, 1, a / c, 1 }, 1)
                    + G({ 0, f / c, b / c, a / c, 1, 1 }, 1)
                    + G({ b / c, 0, a / c, f / c, 1, 1 }, 1)
                    + G({ b / c, 0, f / c, 1, 1, a / c }, 1)
                    + G({ b / c, 0, f / c, 1, a / c, 1 }, 1)
                    + G({ b / c, 0, f / c, a / c, 1, 1 }, 1)
                    + G({ b / c, a / c, 0, f / c, 1, 1 }, 1)
                    + G({ b / c, a / c, f / c, 1, 1, x / c }, 1)
                    + G({ b / c, a / c, f / c, 1, x / c, 1 }, 1)
                    + G({ b / c, a / c, f / c, x / c, 1, 1 }, 1)
                    + G({ b / c, f / c, 0, 1, 1, a / c }, 1)
                    + G({ b / c, f / c, 0, 1, a / c, 1 }, 1)
                    + G({ b / c, f / c, 0, a / c, 1, 1 }, 1)
                    + G({ b / c, f / c, 1, 0, 1, a / c }, 1)
                    + G({ b / c, f / c, 1, 0, a / c, 1 }, 1)
                    + G({ b / c, f / c, 1, 1, 0, a / c }, 1)
                    + G({ b / c, f / c, 1, 1, a / c, x / c }, 1)
                    + G({ b / c, f / c, 1, a / c, 1, x / c }, 1)
                    + G({ b / c, f / c, 1, a / c, x / c, 1 }, 1)
                    + G({ b / c, f / c, a / c, 1, 1, x / c }, 1)
                    + G({ b / c, f / c, a / c, 1, x / c, 1 }, 1)
                    + G({ b / c, f / c, a / c, x / c, 1, 1 }, 1)
                    + G({ f / c, 0, 1, 1, b / c, a / c }, 1)
                    + G({ f / c, 0, 1, b / c, 1, a / c }, 1)
                    + G({ f / c, 0, 1, b / c, a / c, 1 }, 1)
                    + G({ f / c, 0, b / c, 1, 1, a / c }, 1)
                    + G({ f / c, 0, b / c, 1, a / c, 1 }, 1)
                    + G({ f / c, 0, b / c, a / c, 1, 1 }, 1)
                    + G({ f / c, 1, 0, 1, b / c, a / c }, 1)
                    + G({ f / c, 1, 0, b / c, 1, a / c }, 1)
                    + G({ f / c, 1, 0, b / c, a / c, 1 }, 1)
                    + G({ f / c, 1, 1, 0, b / c, a / c }, 1)
                    + G({ f / c, 1, 1, b / c, 0, a / c }, 1)
                    + G({ f / c, 1, 1, b / c, a / c, x / c }, 1)
                    + G({ f / c, 1, b / c, 0, 1, a / c }, 1)
                    + G({ f / c, 1, b / c, 0, a / c, 1 }, 1)
                    + G({ f / c, 1, b / c, 1, 0, a / c }, 1)
                    + G({ f / c, 1, b / c, 1, a / c, x / c }, 1)
                    + G({ f / c, 1, b / c, a / c, 1, x / c }, 1)
                    + G({ f / c, 1, b / c, a / c, x / c, 1 }, 1)
                    + G({ f / c, b / c, 0, 1, 1, a / c }, 1)
                    + G({ f / c, b / c, 0, 1, a / c, 1 }, 1)
                    + G({ f / c, b / c, 0, a / c, 1, 1 }, 1)
                    + G({ f / c, b / c, 1, 0, 1, a / c }, 1)
                    + G({ f / c, b / c, 1, 0, a / c, 1 }, 1)
                    + G({ f / c, b / c, 1, 1, 0, a / c }, 1)
                    + G({ f / c, b / c, 1, 1, a / c, x / c }, 1)
                    + G({ f / c, b / c, 1, a / c, 1, x / c }, 1)
                    + G({ f / c, b / c, 1, a / c, x / c, 1 }, 1)
                    + G({ f / c, b / c, a / c, 1, 1, x / c }, 1)
                    + G({ f / c, b / c, a / c, 1, x / c, 1 }, 1)
                    + G({ f / c, b / c, a / c, x / c, 1, 1 }, 1)
                    + (sy[9] - sy[7]) * G({ a, f }, { sa, sf }, x)
                    + sy[4]
                        * (-G({ 0, a, b }, { 1, sa, sb }, x)
                            - G({ a, 0, b }, { sa, 1, sb }, x))
                    + (sy[1] - sy[4]) * G({ a, b, f }, { sa, sb, sf }, x)
                    + G({ b / c, a / c }, 1)
                        * (-sy[5] + G({ 0, c, c, f }, { 1, sc, sc, sf }, x))
                    + G({ b / c, 1 }, 1)
                        * (sy[5] - G({ a, b, c, f }, { sa, sb, sc, sf }, x))
                    + G({ b / c }, 1)
                        * (-G({ a, 0, c, c, f }, { sa, 1, sc, sc, sf }, x)
                            + G({ a, b, c, c, f }, { sa, sb, sc, sc, sf }, x))
                    + G({ a, b, 0, c, c, f }, { sa, sb, 1, sc, sc, sf }, x)
                    + ((sy[9] + sy[10] + sy[11] + sy[12]) * sy[14] - sy[15] - sy[16]
                          - sy[17] - sy[18] - sy[19] - sy[20] - sy[21] - sy[22] - sy[23]
                          - sy[24] - sy[13] * sy[4])
                        * Log(c, sc)
                    + ((-sy[9] - sy[10] - sy[11] - sy[12]) * sy[14] + sy[15] + sy[16]
                          + sy[17] + sy[18] + sy[19] + sy[20] + sy[21] + sy[22] + sy[23]
                          + sy[24] + sy[13] * sy[4])
                        * Log(-x, sc) };

                if (f != x) {
                    res += (-sy[15] + G({ b / c, a / c, x / c, 1, 1 }, 1))
                        * G({ f }, { sf }, x);
                }
                return res;
            }
        } else { // abccef
            const vector<complex<double>> sy = { G({ b / c, a / c, 1 }, 1),
                G({ b / c, e / c, 1 }, 1), G({ c, e, f }, { sc, se, sf }, x),
                G({ a, b, f }, { sa, sb, sf }, x), G({ e / c, 1, b / c }, 1),
                G({ 0, a, b }, { 1, sa, sb }, x), G({ a, 0, b }, { sa, 1, sb }, x),
                G({ e / c, 1, f / c }, 1), G({ e / c, b / c, 1 }, 1),
                G({ e / c, f / c, 1 }, 1), G({ a, b, e, f }, { sa, sb, se, sf }, x),
                G({ a, c, e, f }, { sa, sc, se, sf }, x), G({ a, f }, { sa, sf }, x),
                G({ b / c, e / c, 1, a / c }, 1), G({ 0, a }, { 1, sa }, x),
                G({ b / c, e / c, 1, f / c }, 1), G({ b / c, a / c, e / c, 1 }, 1),
                G({ b / c, e / c, a / c, 1 }, 1), G({ e / c, 1, b / c, a / c }, 1),
                G({ e / c, b / c, 1, a / c }, 1), G({ e / c, b / c, a / c, 1 }, 1),
                G({ b / c, e / c, f / c, 1 }, 1), G({ e / c, 1, b / c, f / c }, 1),
                G({ e / c, b / c, 1, f / c }, 1), G({ e / c, b / c, f / c, 1 }, 1),
                G({ a, b }, { sa, sb }, x), G({ e / c, 1, f / c, b / c }, 1),
                G({ e / c, f / c, 1, b / c }, 1), G({ e / c, f / c, b / c, 1 }, 1),
                G({ a }, { sa }, x), G({ b / c, a / c, e / c, 1, f / c }, 1),
                G({ b / c, a / c, e / c, f / c, 1 }, 1),
                G({ b / c, e / c, 1, a / c, f / c }, 1),
                G({ b / c, e / c, 1, f / c, a / c }, 1),
                G({ b / c, e / c, a / c, 1, f / c }, 1),
                G({ b / c, e / c, a / c, f / c, 1 }, 1),
                G({ b / c, e / c, f / c, 1, a / c }, 1),
                G({ b / c, e / c, f / c, a / c, 1 }, 1),
                G({ e / c, 1, b / c, a / c, f / c }, 1),
                G({ e / c, 1, b / c, f / c, a / c }, 1),
                G({ e / c, 1, f / c, b / c, a / c }, 1),
                G({ e / c, b / c, 1, a / c, f / c }, 1),
                G({ e / c, b / c, 1, f / c, a / c }, 1),
                G({ e / c, b / c, a / c, 1, f / c }, 1),
                G({ e / c, b / c, a / c, f / c, 1 }, 1),
                G({ e / c, b / c, f / c, 1, a / c }, 1),
                G({ e / c, b / c, f / c, a / c, 1 }, 1),
                G({ e / c, f / c, 1, b / c, a / c }, 1),
                G({ e / c, f / c, b / c, 1, a / c }, 1),
                G({ e / c, f / c, b / c, a / c, 1 }, 1) };
            complex<double> res { sy[12] * sy[13] - sy[14] * sy[15]
                + sy[12]
                    * (-sy[15] + sy[17] + sy[18] + sy[19] + sy[20] - sy[21] - sy[22]
                        - sy[23] - sy[24])
                + sy[14] * (-sy[21] - sy[22] - sy[23] - sy[24] - sy[26] - sy[27] - sy[28])
                - sy[0] * sy[2] - sy[3] * sy[4] + sy[9] * sy[5] + sy[9] * sy[6]
                + (sy[5] + sy[6]) * sy[7] + sy[3] * (sy[9] + sy[7] - sy[8])
                + (-sy[10] + sy[11]) * G({ b / c, 1 }, 1)
                + sy[2] * G({ b / c, a / c, x / c }, 1)
                + sy[25]
                    * (-sy[26] - sy[27] - sy[28] - G({ 0, e / c, 1, f / c }, 1)
                        - G({ 0, e / c, f / c, 1 }, 1) - G({ e / c, 0, 1, f / c }, 1)
                        - G({ e / c, 0, f / c, 1 }, 1) - G({ e / c, 1, 0, f / c }, 1))
                + sy[29]
                    * (sy[33] + sy[36] + sy[37] + sy[39] + sy[40] + sy[42] + sy[45]
                        + sy[46] + sy[47] + sy[48] + sy[49]
                        + G({ 0, b / c, e / c, 1, f / c }, 1)
                        + G({ 0, b / c, e / c, f / c, 1 }, 1)
                        + G({ 0, e / c, 1, b / c, f / c }, 1)
                        + G({ 0, e / c, 1, f / c, b / c }, 1)
                        + G({ 0, e / c, b / c, 1, f / c }, 1)
                        + G({ 0, e / c, b / c, f / c, 1 }, 1)
                        + G({ 0, e / c, f / c, 1, b / c }, 1)
                        + G({ 0, e / c, f / c, b / c, 1 }, 1)
                        + G({ b / c, 0, e / c, 1, f / c }, 1)
                        + G({ b / c, 0, e / c, f / c, 1 }, 1)
                        + G({ b / c, e / c, 0, 1, f / c }, 1)
                        + G({ b / c, e / c, 0, f / c, 1 }, 1)
                        + G({ b / c, e / c, 1, 0, f / c }, 1)
                        + G({ e / c, 0, 1, b / c, f / c }, 1)
                        + G({ e / c, 0, 1, f / c, b / c }, 1)
                        + G({ e / c, 0, b / c, 1, f / c }, 1)
                        + G({ e / c, 0, b / c, f / c, 1 }, 1)
                        + G({ e / c, 0, f / c, 1, b / c }, 1)
                        + G({ e / c, 0, f / c, b / c, 1 }, 1)
                        + G({ e / c, 1, 0, b / c, f / c }, 1)
                        + G({ e / c, 1, 0, f / c, b / c }, 1)
                        + G({ e / c, 1, b / c, 0, f / c }, 1)
                        + G({ e / c, 1, f / c, 0, b / c }, 1)
                        + G({ e / c, b / c, 0, 1, f / c }, 1)
                        + G({ e / c, b / c, 0, f / c, 1 }, 1)
                        + G({ e / c, b / c, 1, 0, f / c }, 1)
                        + G({ e / c, f / c, 0, 1, b / c }, 1)
                        + G({ e / c, f / c, 0, b / c, 1 }, 1)
                        + G({ e / c, f / c, 1, 0, b / c }, 1))
                - G({ 0, b / c, a / c, e / c, 1, f / c }, 1)
                - G({ 0, b / c, a / c, e / c, f / c, 1 }, 1)
                - G({ 0, b / c, e / c, 1, a / c, f / c }, 1)
                - G({ 0, b / c, e / c, 1, f / c, a / c }, 1)
                - G({ 0, b / c, e / c, a / c, 1, f / c }, 1)
                - G({ 0, b / c, e / c, a / c, f / c, 1 }, 1)
                - G({ 0, b / c, e / c, f / c, 1, a / c }, 1)
                - G({ 0, b / c, e / c, f / c, a / c, 1 }, 1)
                - G({ 0, e / c, 1, b / c, a / c, f / c }, 1)
                - G({ 0, e / c, 1, b / c, f / c, a / c }, 1)
                - G({ 0, e / c, 1, f / c, b / c, a / c }, 1)
                - G({ 0, e / c, b / c, 1, a / c, f / c }, 1)
                - G({ 0, e / c, b / c, 1, f / c, a / c }, 1)
                - G({ 0, e / c, b / c, a / c, 1, f / c }, 1)
                - G({ 0, e / c, b / c, a / c, f / c, 1 }, 1)
                - G({ 0, e / c, b / c, f / c, 1, a / c }, 1)
                - G({ 0, e / c, b / c, f / c, a / c, 1 }, 1)
                - G({ 0, e / c, f / c, 1, b / c, a / c }, 1)
                - G({ 0, e / c, f / c, b / c, 1, a / c }, 1)
                - G({ 0, e / c, f / c, b / c, a / c, 1 }, 1)
                - G({ b / c, 0, a / c, e / c, 1, f / c }, 1)
                - G({ b / c, 0, a / c, e / c, f / c, 1 }, 1)
                - G({ b / c, 0, e / c, 1, a / c, f / c }, 1)
                - G({ b / c, 0, e / c, 1, f / c, a / c }, 1)
                - G({ b / c, 0, e / c, a / c, 1, f / c }, 1)
                - G({ b / c, 0, e / c, a / c, f / c, 1 }, 1)
                - G({ b / c, 0, e / c, f / c, 1, a / c }, 1)
                - G({ b / c, 0, e / c, f / c, a / c, 1 }, 1)
                - G({ b / c, a / c, 0, e / c, 1, f / c }, 1)
                - G({ b / c, a / c, 0, e / c, f / c, 1 }, 1)
                - G({ b / c, a / c, e / c, 0, 1, f / c }, 1)
                - G({ b / c, a / c, e / c, 0, f / c, 1 }, 1)
                - G({ b / c, a / c, e / c, 1, 0, f / c }, 1)
                - G({ b / c, a / c, e / c, 1, f / c, x / c }, 1)
                - G({ b / c, a / c, e / c, f / c, 1, x / c }, 1)
                - G({ b / c, a / c, e / c, f / c, x / c, 1 }, 1)
                - G({ b / c, e / c, 0, 1, a / c, f / c }, 1)
                - G({ b / c, e / c, 0, 1, f / c, a / c }, 1)
                - G({ b / c, e / c, 0, a / c, 1, f / c }, 1)
                - G({ b / c, e / c, 0, a / c, f / c, 1 }, 1)
                - G({ b / c, e / c, 0, f / c, 1, a / c }, 1)
                - G({ b / c, e / c, 0, f / c, a / c, 1 }, 1)
                - G({ b / c, e / c, 1, 0, a / c, f / c }, 1)
                - G({ b / c, e / c, 1, 0, f / c, a / c }, 1)
                - G({ b / c, e / c, 1, a / c, 0, f / c }, 1)
                - G({ b / c, e / c, 1, a / c, f / c, x / c }, 1)
                - G({ b / c, e / c, 1, f / c, 0, a / c }, 1)
                - G({ b / c, e / c, 1, f / c, a / c, x / c }, 1)
                - G({ b / c, e / c, a / c, 0, 1, f / c }, 1)
                - G({ b / c, e / c, a / c, 0, f / c, 1 }, 1)
                - G({ b / c, e / c, a / c, 1, 0, f / c }, 1)
                - G({ b / c, e / c, a / c, 1, f / c, x / c }, 1)
                - G({ b / c, e / c, a / c, f / c, 1, x / c }, 1)
                - G({ b / c, e / c, a / c, f / c, x / c, 1 }, 1)
                - G({ b / c, e / c, f / c, 0, 1, a / c }, 1)
                - G({ b / c, e / c, f / c, 0, a / c, 1 }, 1)
                - G({ b / c, e / c, f / c, 1, 0, a / c }, 1)
                - G({ b / c, e / c, f / c, 1, a / c, x / c }, 1)
                - G({ b / c, e / c, f / c, a / c, 1, x / c }, 1)
                - G({ b / c, e / c, f / c, a / c, x / c, 1 }, 1)
                - G({ e / c, 0, 1, b / c, a / c, f / c }, 1)
                - G({ e / c, 0, 1, b / c, f / c, a / c }, 1)
                - G({ e / c, 0, 1, f / c, b / c, a / c }, 1)
                - G({ e / c, 0, b / c, 1, a / c, f / c }, 1)
                - G({ e / c, 0, b / c, 1, f / c, a / c }, 1)
                - G({ e / c, 0, b / c, a / c, 1, f / c }, 1)
                - G({ e / c, 0, b / c, a / c, f / c, 1 }, 1)
                - G({ e / c, 0, b / c, f / c, 1, a / c }, 1)
                - G({ e / c, 0, b / c, f / c, a / c, 1 }, 1)
                - G({ e / c, 0, f / c, 1, b / c, a / c }, 1)
                - G({ e / c, 0, f / c, b / c, 1, a / c }, 1)
                - G({ e / c, 0, f / c, b / c, a / c, 1 }, 1)
                - G({ e / c, 1, 0, b / c, a / c, f / c }, 1)
                - G({ e / c, 1, 0, b / c, f / c, a / c }, 1)
                - G({ e / c, 1, 0, f / c, b / c, a / c }, 1)
                - G({ e / c, 1, b / c, 0, a / c, f / c }, 1)
                - G({ e / c, 1, b / c, 0, f / c, a / c }, 1)
                - G({ e / c, 1, b / c, a / c, 0, f / c }, 1)
                - G({ e / c, 1, b / c, a / c, f / c, x / c }, 1)
                - G({ e / c, 1, b / c, f / c, 0, a / c }, 1)
                - G({ e / c, 1, b / c, f / c, a / c, x / c }, 1)
                - G({ e / c, 1, f / c, 0, b / c, a / c }, 1)
                - G({ e / c, 1, f / c, b / c, 0, a / c }, 1)
                - G({ e / c, 1, f / c, b / c, a / c, x / c }, 1)
                - G({ e / c, b / c, 0, 1, a / c, f / c }, 1)
                - G({ e / c, b / c, 0, 1, f / c, a / c }, 1)
                - G({ e / c, b / c, 0, a / c, 1, f / c }, 1)
                - G({ e / c, b / c, 0, a / c, f / c, 1 }, 1)
                - G({ e / c, b / c, 0, f / c, 1, a / c }, 1)
                - G({ e / c, b / c, 0, f / c, a / c, 1 }, 1)
                - G({ e / c, b / c, 1, 0, a / c, f / c }, 1)
                - G({ e / c, b / c, 1, 0, f / c, a / c }, 1)
                - G({ e / c, b / c, 1, a / c, 0, f / c }, 1)
                - G({ e / c, b / c, 1, a / c, f / c, x / c }, 1)
                - G({ e / c, b / c, 1, f / c, 0, a / c }, 1)
                - G({ e / c, b / c, 1, f / c, a / c, x / c }, 1)
                - G({ e / c, b / c, a / c, 0, 1, f / c }, 1)
                - G({ e / c, b / c, a / c, 0, f / c, 1 }, 1)
                - G({ e / c, b / c, a / c, 1, 0, f / c }, 1)
                - G({ e / c, b / c, a / c, 1, f / c, x / c }, 1)
                - G({ e / c, b / c, a / c, f / c, 1, x / c }, 1)
                - G({ e / c, b / c, a / c, f / c, x / c, 1 }, 1)
                - G({ e / c, b / c, f / c, 0, 1, a / c }, 1)
                - G({ e / c, b / c, f / c, 0, a / c, 1 }, 1)
                - G({ e / c, b / c, f / c, 1, 0, a / c }, 1)
                - G({ e / c, b / c, f / c, 1, a / c, x / c }, 1)
                - G({ e / c, b / c, f / c, a / c, 1, x / c }, 1)
                - G({ e / c, b / c, f / c, a / c, x / c, 1 }, 1)
                - G({ e / c, f / c, 0, 1, b / c, a / c }, 1)
                - G({ e / c, f / c, 0, b / c, 1, a / c }, 1)
                - G({ e / c, f / c, 0, b / c, a / c, 1 }, 1)
                - G({ e / c, f / c, 1, 0, b / c, a / c }, 1)
                - G({ e / c, f / c, 1, b / c, 0, a / c }, 1)
                - G({ e / c, f / c, 1, b / c, a / c, x / c }, 1)
                - G({ e / c, f / c, b / c, 0, 1, a / c }, 1)
                - G({ e / c, f / c, b / c, 0, a / c, 1 }, 1)
                - G({ e / c, f / c, b / c, 1, 0, a / c }, 1)
                - G({ e / c, f / c, b / c, 1, a / c, x / c }, 1)
                - G({ e / c, f / c, b / c, a / c, 1, x / c }, 1)
                - G({ e / c, f / c, b / c, a / c, x / c, 1 }, 1)
                + (-sy[13] - sy[16] - sy[17] - sy[18] - sy[19] - sy[20])
                    * G({ 0, f }, { 1, sf }, x)
                + (sy[1] + sy[4] + sy[8]) * G({ a, 0, f }, { sa, 1, sf }, x)
                + (sy[0] - sy[1]) * G({ a, e, f }, { sa, se, sf }, x)
                + G({ b / c, a / c }, 1)
                    * (-sy[11] + G({ 0, c, e, f }, { 1, sc, se, sf }, x))
                + G({ e / c, 1 }, 1) * (sy[10] - G({ a, b, 0, f }, { sa, sb, 1, sf }, x))
                + G({ b / c }, 1)
                    * (-G({ a, 0, c, e, f }, { sa, 1, sc, se, sf }, x)
                        + G({ a, b, c, e, f }, { sa, sb, sc, se, sf }, x))
                + G({ a, b, 0, c, e, f }, { sa, sb, 1, sc, se, sf }, x)
                + ((-sy[15] - sy[21] - sy[22] - sy[23] - sy[24] - sy[26] - sy[27]
                       - sy[28])
                          * sy[29]
                      + sy[30] + sy[31] + sy[32] + sy[33] + sy[34] + sy[35] + sy[36]
                      + sy[37] + sy[38] + sy[39] + sy[40] + sy[41] + sy[42] + sy[43]
                      + sy[44] + sy[45] + sy[46] + sy[47] + sy[48] + sy[49]
                      + sy[25] * (sy[9] + sy[7]))
                    * Log(c, sc)
                + ((sy[15] + sy[21] + sy[22] + sy[23] + sy[24] + sy[26] + sy[27] + sy[28])
                          * sy[29]
                      - sy[30] - sy[31] - sy[32] - sy[33] - sy[34] - sy[35] - sy[36]
                      - sy[37] - sy[38] - sy[39] - sy[40] - sy[41] - sy[42] - sy[43]
                      - sy[44] - sy[45] - sy[46] - sy[47] - sy[48] - sy[49]
                      + sy[25] * (-sy[9] - sy[7]))
                    * Log(-x, sc) };
            if (e != x) {
                res += (sy[16] - G({ b / c, a / c, x / c, 1 }, 1))
                    * G({ e, f }, { se, sf }, x);
            }
            if (f != x) {
                res += (sy[30] + sy[31] + sy[32] + sy[34] + sy[35] + sy[38] + sy[41]
                           + sy[43] + sy[44] - G({ b / c, a / c, e / c, 1, x / c }, 1)
                           - G({ b / c, a / c, e / c, x / c, 1 }, 1)
                           - G({ b / c, e / c, 1, a / c, x / c }, 1)
                           - G({ b / c, e / c, a / c, 1, x / c }, 1)
                           - G({ b / c, e / c, a / c, x / c, 1 }, 1)
                           - G({ e / c, 1, b / c, a / c, x / c }, 1)
                           - G({ e / c, b / c, 1, a / c, x / c }, 1)
                           - G({ e / c, b / c, a / c, 1, x / c }, 1)
                           - G({ e / c, b / c, a / c, x / c, 1 }, 1))
                    * G({ f }, { sf }, x);
            }
            return res;
        }
    } else { // abcdef
        const vector<complex<double>> sy = { G({ a, e, f }, { sa, se, sf }, x),
            G({ b / c, d / c, a / c }, 1), G({ b / c, a / c, d / c }, 1),
            G({ d / c, b / c, a / c }, 1), G({ b / c, d / c, e / c }, 1),
            G({ d / c, b / c, e / c }, 1), G({ d / c, e / c, b / c }, 1),
            G({ a, b, f }, { sa, sb, sf }, x), G({ d / c, e / c, f / c }, 1),
            G({ d / c, b / c }, 1), G({ a, 0, e, f }, { sa, 1, se, sf }, x),
            G({ a, b, e, f }, { sa, sb, se, sf }, x),
            G({ a, d, e, f }, { sa, sd, se, sf }, x), G({ a, f }, { sa, sf }, x),
            G({ b / c, d / c, e / c, a / c }, 1), G({ 0, a }, { 1, sa }, x),
            G({ b / c, d / c, e / c, f / c }, 1), G({ b / c, a / c, d / c, e / c }, 1),
            G({ b / c, d / c, a / c, e / c }, 1), G({ d / c, b / c, a / c, e / c }, 1),
            G({ d / c, b / c, e / c, a / c }, 1), G({ d / c, e / c, b / c, a / c }, 1),
            G({ d / c, b / c, e / c, f / c }, 1), G({ d / c, e / c, b / c, f / c }, 1),
            G({ a, b }, { sa, sb }, x), G({ d / c, e / c, f / c, b / c }, 1),
            G({ a, b, d, e, f }, { sa, sb, sd, se, sf }, x), G({ a }, { sa }, x),
            G({ b / c, a / c, d / c, e / c, f / c }, 1),
            G({ b / c, d / c, a / c, e / c, f / c }, 1),
            G({ b / c, d / c, e / c, a / c, f / c }, 1),
            G({ b / c, d / c, e / c, f / c, a / c }, 1),
            G({ d / c, b / c, a / c, e / c, f / c }, 1),
            G({ d / c, b / c, e / c, a / c, f / c }, 1),
            G({ d / c, b / c, e / c, f / c, a / c }, 1),
            G({ d / c, e / c, b / c, a / c, f / c }, 1),
            G({ d / c, e / c, b / c, f / c, a / c }, 1),
            G({ d / c, e / c, f / c, b / c, a / c }, 1) };
        complex<double> res { -(sy[9] * sy[10]) + sy[9] * sy[11] - sy[13] * sy[14]
            + sy[15] * sy[16] - sy[0] * sy[1]
            + sy[13] * (sy[16] - sy[20] - sy[21] + sy[22] + sy[23])
            + sy[15] * (sy[22] + sy[23] + sy[25]) + sy[0] * (-sy[3] + sy[4] + sy[5])
            + sy[6] * sy[7] - sy[7] * sy[8] + (-sy[10] + sy[12]) * G({ b / c, d / c }, 1)
            + sy[24]
                * (sy[25] + G({ 0, d / c, e / c, f / c }, 1)
                    + G({ d / c, 0, e / c, f / c }, 1) + G({ d / c, e / c, 0, f / c }, 1))
            + sy[27]
                * (-sy[31] - sy[34] - sy[36] - sy[37]
                    - G({ 0, b / c, d / c, e / c, f / c }, 1)
                    - G({ 0, d / c, b / c, e / c, f / c }, 1)
                    - G({ 0, d / c, e / c, b / c, f / c }, 1)
                    - G({ 0, d / c, e / c, f / c, b / c }, 1)
                    - G({ b / c, 0, d / c, e / c, f / c }, 1)
                    - G({ b / c, d / c, 0, e / c, f / c }, 1)
                    - G({ b / c, d / c, e / c, 0, f / c }, 1)
                    - G({ d / c, 0, b / c, e / c, f / c }, 1)
                    - G({ d / c, 0, e / c, b / c, f / c }, 1)
                    - G({ d / c, 0, e / c, f / c, b / c }, 1)
                    - G({ d / c, b / c, 0, e / c, f / c }, 1)
                    - G({ d / c, b / c, e / c, 0, f / c }, 1)
                    - G({ d / c, e / c, 0, b / c, f / c }, 1)
                    - G({ d / c, e / c, 0, f / c, b / c }, 1)
                    - G({ d / c, e / c, b / c, 0, f / c }, 1)
                    - G({ d / c, e / c, f / c, 0, b / c }, 1))
            + G({ 0, b / c, a / c, d / c, e / c, f / c }, 1)
            + G({ 0, b / c, d / c, a / c, e / c, f / c }, 1)
            + G({ 0, b / c, d / c, e / c, a / c, f / c }, 1)
            + G({ 0, b / c, d / c, e / c, f / c, a / c }, 1)
            + G({ 0, d / c, b / c, a / c, e / c, f / c }, 1)
            + G({ 0, d / c, b / c, e / c, a / c, f / c }, 1)
            + G({ 0, d / c, b / c, e / c, f / c, a / c }, 1)
            + G({ 0, d / c, e / c, b / c, a / c, f / c }, 1)
            + G({ 0, d / c, e / c, b / c, f / c, a / c }, 1)
            + G({ 0, d / c, e / c, f / c, b / c, a / c }, 1)
            + G({ b / c, 0, a / c, d / c, e / c, f / c }, 1)
            + G({ b / c, 0, d / c, a / c, e / c, f / c }, 1)
            + G({ b / c, 0, d / c, e / c, a / c, f / c }, 1)
            + G({ b / c, 0, d / c, e / c, f / c, a / c }, 1)
            + G({ b / c, a / c, 0, d / c, e / c, f / c }, 1)
            + G({ b / c, a / c, d / c, 0, e / c, f / c }, 1)
            + G({ b / c, a / c, d / c, e / c, 0, f / c }, 1)
            + G({ b / c, a / c, d / c, e / c, f / c, x / c }, 1)
            + G({ b / c, d / c, 0, a / c, e / c, f / c }, 1)
            + G({ b / c, d / c, 0, e / c, a / c, f / c }, 1)
            + G({ b / c, d / c, 0, e / c, f / c, a / c }, 1)
            + G({ b / c, d / c, a / c, 0, e / c, f / c }, 1)
            + G({ b / c, d / c, a / c, e / c, 0, f / c }, 1)
            + G({ b / c, d / c, a / c, e / c, f / c, x / c }, 1)
            + G({ b / c, d / c, e / c, 0, a / c, f / c }, 1)
            + G({ b / c, d / c, e / c, 0, f / c, a / c }, 1)
            + G({ b / c, d / c, e / c, a / c, 0, f / c }, 1)
            + G({ b / c, d / c, e / c, a / c, f / c, x / c }, 1)
            + G({ b / c, d / c, e / c, f / c, 0, a / c }, 1)
            + G({ b / c, d / c, e / c, f / c, a / c, x / c }, 1)
            + G({ d / c, 0, b / c, a / c, e / c, f / c }, 1)
            + G({ d / c, 0, b / c, e / c, a / c, f / c }, 1)
            + G({ d / c, 0, b / c, e / c, f / c, a / c }, 1)
            + G({ d / c, 0, e / c, b / c, a / c, f / c }, 1)
            + G({ d / c, 0, e / c, b / c, f / c, a / c }, 1)
            + G({ d / c, 0, e / c, f / c, b / c, a / c }, 1)
            + G({ d / c, b / c, 0, a / c, e / c, f / c }, 1)
            + G({ d / c, b / c, 0, e / c, a / c, f / c }, 1)
            + G({ d / c, b / c, 0, e / c, f / c, a / c }, 1)
            + G({ d / c, b / c, a / c, 0, e / c, f / c }, 1)
            + G({ d / c, b / c, a / c, e / c, 0, f / c }, 1)
            + G({ d / c, b / c, a / c, e / c, f / c, x / c }, 1)
            + G({ d / c, b / c, e / c, 0, a / c, f / c }, 1)
            + G({ d / c, b / c, e / c, 0, f / c, a / c }, 1)
            + G({ d / c, b / c, e / c, a / c, 0, f / c }, 1)
            + G({ d / c, b / c, e / c, a / c, f / c, x / c }, 1)
            + G({ d / c, b / c, e / c, f / c, 0, a / c }, 1)
            + G({ d / c, b / c, e / c, f / c, a / c, x / c }, 1)
            + G({ d / c, e / c, 0, b / c, a / c, f / c }, 1)
            + G({ d / c, e / c, 0, b / c, f / c, a / c }, 1)
            + G({ d / c, e / c, 0, f / c, b / c, a / c }, 1)
            + G({ d / c, e / c, b / c, 0, a / c, f / c }, 1)
            + G({ d / c, e / c, b / c, 0, f / c, a / c }, 1)
            + G({ d / c, e / c, b / c, a / c, 0, f / c }, 1)
            + G({ d / c, e / c, b / c, a / c, f / c, x / c }, 1)
            + G({ d / c, e / c, b / c, f / c, 0, a / c }, 1)
            + G({ d / c, e / c, b / c, f / c, a / c, x / c }, 1)
            + G({ d / c, e / c, f / c, 0, b / c, a / c }, 1)
            + G({ d / c, e / c, f / c, b / c, 0, a / c }, 1)
            + G({ d / c, e / c, f / c, b / c, a / c, x / c }, 1)
            + (sy[14] + sy[17] + sy[18] + sy[19] + sy[20] + sy[21])
                * G({ 0, f }, { 1, sf }, x)
            + (sy[1] + sy[2] + sy[3]) * G({ 0, e, f }, { 1, se, sf }, x)
            + sy[8]
                * (-G({ 0, a, b }, { 1, sa, sb }, x) - G({ a, 0, b }, { sa, 1, sb }, x))
            + (-sy[4] - sy[5] - sy[6]) * G({ a, 0, f }, { sa, 1, sf }, x)
            + G({ b / c, a / c }, 1) * (-sy[12] + G({ 0, d, e, f }, { 1, sd, se, sf }, x))
            + G({ d / c, e / c }, 1) * (-sy[11] + G({ a, b, 0, f }, { sa, sb, 1, sf }, x))
            + G({ b / c }, 1) * (sy[26] - G({ a, 0, d, e, f }, { sa, 1, sd, se, sf }, x))
            + G({ d / c }, 1) * (-sy[26] + G({ a, b, 0, e, f }, { sa, sb, 1, se, sf }, x))
            + G({ a, b, 0, d, e, f }, { sa, sb, 1, sd, se, sf }, x)
            + ((sy[16] + sy[22] + sy[23] + sy[25]) * sy[27] - sy[28] - sy[29] - sy[30]
                  - sy[31] - sy[32] - sy[33] - sy[34] - sy[35] - sy[36] - sy[37]
                  - sy[24] * sy[8])
                * Log(c, sc)
            + ((-sy[16] - sy[22] - sy[23] - sy[25]) * sy[27] + sy[28] + sy[29] + sy[30]
                  + sy[31] + sy[32] + sy[33] + sy[34] + sy[35] + sy[36] + sy[37]
                  + sy[24] * sy[8])
                * Log(-x, sc) };
        if (d != x) {
            res += (-sy[2] + G({ b / c, a / c, x / c }, 1))
                * G({ d, e, f }, { sd, se, sf }, x);
        }
        if (e != x) {
            res += (-sy[17] - sy[18] - sy[19] + G({ b / c, a / c, d / c, x / c }, 1)
                       + G({ b / c, d / c, a / c, x / c }, 1)
                       + G({ d / c, b / c, a / c, x / c }, 1))
                * G({ e, f }, { se, sf }, x);
        }
        if (f != x) {
            res += (-sy[28] - sy[29] - sy[30] - sy[32] - sy[33] - sy[35]
                       + G({ b / c, a / c, d / c, e / c, x / c }, 1)
                       + G({ b / c, d / c, a / c, e / c, x / c }, 1)
                       + G({ b / c, d / c, e / c, a / c, x / c }, 1)
                       + G({ d / c, b / c, a / c, e / c, x / c }, 1)
                       + G({ d / c, b / c, e / c, a / c, x / c }, 1)
                       + G({ d / c, e / c, b / c, a / c, x / c }, 1))
                * G({ f }, { sf }, x);
        }
        return res;
    }
}
complex<double> G6_abcdef_b(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, complex<double> f1, int sa, int sb, int sc,
    int sd, int se, int sf, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    complex<double> f = f1 / x1;
    double x = 1.;
    if (b == c) {
        if (b == d) {
            if (b == e) {
                if (b == f) { // abbbbb
                    const vector<complex<double>> sy = { G({ a / b, 1, 1 }, 1),
                        G({ b, b, b }, { sb, sb, sb }, x), G({ a / b, 1 }, 1),
                        G({ a / b, 1, 1, 1 }, 1), G({ b, b }, { sb, sb }, x),
                        G({ b, b, b, b }, { sb, sb, sb, sb }, x), G({ a / b }, 1),
                        G({ a }, { sa }, x), G({ a / b, 1, 1, 1, 1 }, 1),
                        G({ b }, { sb }, x) };
                    complex<double> res { -(sy[0] * sy[1]) + sy[3] * sy[4] + sy[2] * sy[5]
                        - sy[9] * sy[8] + sy[7] * sy[8] - sy[5] * G({ a / b, x / b }, 1)
                        + sy[1] * G({ a / b, x / b, 1 }, 1)
                        - sy[4] * G({ a / b, x / b, 1, 1 }, 1)
                        + sy[9] * G({ a / b, x / b, 1, 1, 1 }, 1)
                        + G({ a / b, 0, 1, 1, 1, 1 }, 1)
                        - G({ a / b, x / b, 1, 1, 1, 1 }, 1)
                        - sy[3] * G({ a, b }, { sa, sb }, x)
                        + sy[0] * G({ a, b, b }, { sa, sb, sb }, x)
                        - sy[2] * G({ a, b, b, b }, { sa, sb, sb, sb }, x)
                        - sy[6] * G({ 0, b, b, b, b }, { 1, sb, sb, sb, sb }, x)
                        + sy[6] * G({ a, b, b, b, b }, { sa, sb, sb, sb, sb }, x)
                        + G({ a, 0, b, b, b, b }, { sa, 1, sb, sb, sb, sb }, x)
                        - sy[7] * Zeta(5) };

                    return res;
                } else { // abbbbf
                    const vector<complex<double>> sy = { G({ a / b, 1, 1 }, 1),
                        G({ b, b, f }, { sb, sb, sf }, x), G({ b, f }, { sb, sf }, x),
                        G({ a / b, 1, 1, 1 }, 1),
                        G({ b, b, b, f }, { sb, sb, sb, sf }, x),
                        G({ f / b, 1, 1, 1 }, 1), G({ a }, { sa }, x),
                        G({ f / b, 1, 1, 1, a / b }, 1), G({ f / b, 1, 1, a / b, 1 }, 1),
                        G({ f / b, 1, a / b, 1, 1 }, 1), G({ f / b, a / b, 1, 1, 1 }, 1),
                        G({ a / b, f / b, 1, 1, 1 }, 1) };
                    complex<double> res { -(sy[0] * sy[1]) + sy[2] * sy[3]
                        - sy[4] * G({ a / b, x / b }, 1)
                        + sy[1] * G({ a / b, x / b, 1 }, 1)
                        - sy[2] * G({ a / b, x / b, 1, 1 }, 1)
                        + sy[6]
                            * (-sy[9] - sy[10] - sy[7] - sy[8]
                                - G({ 0, f / b, 1, 1, 1 }, 1))
                        + G({ 0, a / b, f / b, 1, 1, 1 }, 1)
                        + G({ 0, f / b, 1, 1, 1, a / b }, 1)
                        + G({ 0, f / b, 1, 1, a / b, 1 }, 1)
                        + G({ 0, f / b, 1, a / b, 1, 1 }, 1)
                        + G({ 0, f / b, a / b, 1, 1, 1 }, 1)
                        + G({ a / b, 0, f / b, 1, 1, 1 }, 1)
                        + G({ a / b, f / b, 1, 1, 1, x / b }, 1)
                        + G({ a / b, f / b, 1, 1, x / b, 1 }, 1)
                        + G({ a / b, f / b, 1, x / b, 1, 1 }, 1)
                        + G({ a / b, f / b, x / b, 1, 1, 1 }, 1)
                        + G({ f / b, 0, 1, 1, 1, a / b }, 1)
                        + G({ f / b, 0, 1, 1, a / b, 1 }, 1)
                        + G({ f / b, 0, 1, a / b, 1, 1 }, 1)
                        + G({ f / b, 0, a / b, 1, 1, 1 }, 1)
                        + G({ f / b, 1, 0, 1, 1, a / b }, 1)
                        + G({ f / b, 1, 0, 1, a / b, 1 }, 1)
                        + G({ f / b, 1, 0, a / b, 1, 1 }, 1)
                        + G({ f / b, 1, 1, 0, 1, a / b }, 1)
                        + G({ f / b, 1, 1, 0, a / b, 1 }, 1)
                        + G({ f / b, 1, 1, 1, 0, a / b }, 1)
                        + G({ f / b, 1, 1, 1, a / b, x / b }, 1)
                        + G({ f / b, 1, 1, a / b, 1, x / b }, 1)
                        + G({ f / b, 1, 1, a / b, x / b, 1 }, 1)
                        + G({ f / b, 1, a / b, 1, 1, x / b }, 1)
                        + G({ f / b, 1, a / b, 1, x / b, 1 }, 1)
                        + G({ f / b, 1, a / b, x / b, 1, 1 }, 1)
                        + G({ f / b, a / b, 1, 1, 1, x / b }, 1)
                        + G({ f / b, a / b, 1, 1, x / b, 1 }, 1)
                        + G({ f / b, a / b, 1, x / b, 1, 1 }, 1)
                        + G({ f / b, a / b, x / b, 1, 1, 1 }, 1)
                        + sy[5] * G({ 0, a }, { 1, sa }, x)
                        + (-sy[3] + sy[5]) * G({ a, f }, { sa, sf }, x)
                        + sy[0] * G({ a, b, f }, { sa, sb, sf }, x)
                        + G({ a / b, 1 }, 1)
                            * (sy[4] - G({ a, b, b, f }, { sa, sb, sb, sf }, x))
                        + G({ a / b }, 1)
                            * (-G({ 0, b, b, b, f }, { 1, sb, sb, sb, sf }, x)
                                + G({ a, b, b, b, f }, { sa, sb, sb, sb, sf }, x))
                        + G({ a, 0, b, b, b, f }, { sa, 1, sb, sb, sb, sf }, x)
                        + (-sy[9] - sy[10] - sy[11] + sy[5] * sy[6] - sy[7] - sy[8])
                            * Log(b, sb)
                        + (sy[9] + sy[10] + sy[11] - sy[5] * sy[6] + sy[7] + sy[8])
                            * Log(-x, sb) };
                    if (f != x) {
                        res += (-sy[11] + G({ a / b, x / b, 1, 1, 1 }, 1))
                            * G({ f }, { sf }, x);
                    }
                    return res;
                }
            } else { // abbbef
                const vector<complex<double>> sy = { G({ a / b, 1, 1 }, 1),
                    G({ b, e, f }, { sb, se, sf }, x), G({ e / b, 1, 1 }, 1),
                    G({ b, b, e, f }, { sb, sb, se, sf }, x), G({ a, f }, { sa, sf }, x),
                    G({ e / b, 1, 1, a / b }, 1), G({ 0, a }, { 1, sa }, x),
                    G({ e / b, 1, 1, f / b }, 1), G({ a / b, e / b, 1, 1 }, 1),
                    G({ e / b, 1, a / b, 1 }, 1), G({ e / b, a / b, 1, 1 }, 1),
                    G({ e / b, 1, f / b, 1 }, 1), G({ e / b, f / b, 1, 1 }, 1),
                    G({ a }, { sa }, x), G({ e / b, 1, 1, f / b, a / b }, 1),
                    G({ e / b, 1, f / b, 1, a / b }, 1),
                    G({ e / b, 1, f / b, a / b, 1 }, 1),
                    G({ e / b, f / b, 1, 1, a / b }, 1),
                    G({ e / b, f / b, 1, a / b, 1 }, 1),
                    G({ e / b, f / b, a / b, 1, 1 }, 1),
                    G({ a / b, e / b, 1, 1, f / b }, 1),
                    G({ a / b, e / b, 1, f / b, 1 }, 1),
                    G({ a / b, e / b, f / b, 1, 1 }, 1),
                    G({ e / b, 1, 1, a / b, f / b }, 1),
                    G({ e / b, 1, a / b, 1, f / b }, 1),
                    G({ e / b, 1, a / b, f / b, 1 }, 1),
                    G({ e / b, a / b, 1, 1, f / b }, 1),
                    G({ e / b, a / b, 1, f / b, 1 }, 1),
                    G({ e / b, a / b, f / b, 1, 1 }, 1) };
                complex<double> res { -(sy[0] * sy[1]) + sy[4] * sy[5]
                    + (-sy[11] - sy[12]) * sy[6]
                    + sy[4] * (sy[9] + sy[10] - sy[11] - sy[12] - sy[7]) - sy[6] * sy[7]
                    - sy[3] * G({ a / b, x / b }, 1) + sy[1] * G({ a / b, x / b, 1 }, 1)
                    + sy[13]
                        * (sy[14] + sy[15] + sy[16] + sy[17] + sy[18] + sy[19]
                            + G({ 0, e / b, 1, 1, f / b }, 1)
                            + G({ 0, e / b, 1, f / b, 1 }, 1)
                            + G({ 0, e / b, f / b, 1, 1 }, 1)
                            + G({ e / b, 0, 1, 1, f / b }, 1)
                            + G({ e / b, 0, 1, f / b, 1 }, 1)
                            + G({ e / b, 0, f / b, 1, 1 }, 1)
                            + G({ e / b, 1, 0, 1, f / b }, 1)
                            + G({ e / b, 1, 0, f / b, 1 }, 1)
                            + G({ e / b, 1, 1, 0, f / b }, 1))
                    - G({ 0, a / b, e / b, 1, 1, f / b }, 1)
                    - G({ 0, a / b, e / b, 1, f / b, 1 }, 1)
                    - G({ 0, a / b, e / b, f / b, 1, 1 }, 1)
                    - G({ 0, e / b, 1, 1, a / b, f / b }, 1)
                    - G({ 0, e / b, 1, 1, f / b, a / b }, 1)
                    - G({ 0, e / b, 1, a / b, 1, f / b }, 1)
                    - G({ 0, e / b, 1, a / b, f / b, 1 }, 1)
                    - G({ 0, e / b, 1, f / b, 1, a / b }, 1)
                    - G({ 0, e / b, 1, f / b, a / b, 1 }, 1)
                    - G({ 0, e / b, a / b, 1, 1, f / b }, 1)
                    - G({ 0, e / b, a / b, 1, f / b, 1 }, 1)
                    - G({ 0, e / b, a / b, f / b, 1, 1 }, 1)
                    - G({ 0, e / b, f / b, 1, 1, a / b }, 1)
                    - G({ 0, e / b, f / b, 1, a / b, 1 }, 1)
                    - G({ 0, e / b, f / b, a / b, 1, 1 }, 1)
                    - G({ a / b, 0, e / b, 1, 1, f / b }, 1)
                    - G({ a / b, 0, e / b, 1, f / b, 1 }, 1)
                    - G({ a / b, 0, e / b, f / b, 1, 1 }, 1)
                    - G({ a / b, e / b, 0, 1, 1, f / b }, 1)
                    - G({ a / b, e / b, 0, 1, f / b, 1 }, 1)
                    - G({ a / b, e / b, 0, f / b, 1, 1 }, 1)
                    - G({ a / b, e / b, 1, 0, 1, f / b }, 1)
                    - G({ a / b, e / b, 1, 0, f / b, 1 }, 1)
                    - G({ a / b, e / b, 1, 1, 0, f / b }, 1)
                    - G({ a / b, e / b, 1, 1, f / b, x / b }, 1)
                    - G({ a / b, e / b, 1, f / b, 1, x / b }, 1)
                    - G({ a / b, e / b, 1, f / b, x / b, 1 }, 1)
                    - G({ a / b, e / b, f / b, 1, 1, x / b }, 1)
                    - G({ a / b, e / b, f / b, 1, x / b, 1 }, 1)
                    - G({ a / b, e / b, f / b, x / b, 1, 1 }, 1)
                    - G({ e / b, 0, 1, 1, a / b, f / b }, 1)
                    - G({ e / b, 0, 1, 1, f / b, a / b }, 1)
                    - G({ e / b, 0, 1, a / b, 1, f / b }, 1)
                    - G({ e / b, 0, 1, a / b, f / b, 1 }, 1)
                    - G({ e / b, 0, 1, f / b, 1, a / b }, 1)
                    - G({ e / b, 0, 1, f / b, a / b, 1 }, 1)
                    - G({ e / b, 0, a / b, 1, 1, f / b }, 1)
                    - G({ e / b, 0, a / b, 1, f / b, 1 }, 1)
                    - G({ e / b, 0, a / b, f / b, 1, 1 }, 1)
                    - G({ e / b, 0, f / b, 1, 1, a / b }, 1)
                    - G({ e / b, 0, f / b, 1, a / b, 1 }, 1)
                    - G({ e / b, 0, f / b, a / b, 1, 1 }, 1)
                    - G({ e / b, 1, 0, 1, a / b, f / b }, 1)
                    - G({ e / b, 1, 0, 1, f / b, a / b }, 1)
                    - G({ e / b, 1, 0, a / b, 1, f / b }, 1)
                    - G({ e / b, 1, 0, a / b, f / b, 1 }, 1)
                    - G({ e / b, 1, 0, f / b, 1, a / b }, 1)
                    - G({ e / b, 1, 0, f / b, a / b, 1 }, 1)
                    - G({ e / b, 1, 1, 0, a / b, f / b }, 1)
                    - G({ e / b, 1, 1, 0, f / b, a / b }, 1)
                    - G({ e / b, 1, 1, a / b, 0, f / b }, 1)
                    - G({ e / b, 1, 1, a / b, f / b, x / b }, 1)
                    - G({ e / b, 1, 1, f / b, 0, a / b }, 1)
                    - G({ e / b, 1, 1, f / b, a / b, x / b }, 1)
                    - G({ e / b, 1, a / b, 0, 1, f / b }, 1)
                    - G({ e / b, 1, a / b, 0, f / b, 1 }, 1)
                    - G({ e / b, 1, a / b, 1, 0, f / b }, 1)
                    - G({ e / b, 1, a / b, 1, f / b, x / b }, 1)
                    - G({ e / b, 1, a / b, f / b, 1, x / b }, 1)
                    - G({ e / b, 1, a / b, f / b, x / b, 1 }, 1)
                    - G({ e / b, 1, f / b, 0, 1, a / b }, 1)
                    - G({ e / b, 1, f / b, 0, a / b, 1 }, 1)
                    - G({ e / b, 1, f / b, 1, 0, a / b }, 1)
                    - G({ e / b, 1, f / b, 1, a / b, x / b }, 1)
                    - G({ e / b, 1, f / b, a / b, 1, x / b }, 1)
                    - G({ e / b, 1, f / b, a / b, x / b, 1 }, 1)
                    - G({ e / b, a / b, 0, 1, 1, f / b }, 1)
                    - G({ e / b, a / b, 0, 1, f / b, 1 }, 1)
                    - G({ e / b, a / b, 0, f / b, 1, 1 }, 1)
                    - G({ e / b, a / b, 1, 0, 1, f / b }, 1)
                    - G({ e / b, a / b, 1, 0, f / b, 1 }, 1)
                    - G({ e / b, a / b, 1, 1, 0, f / b }, 1)
                    - G({ e / b, a / b, 1, 1, f / b, x / b }, 1)
                    - G({ e / b, a / b, 1, f / b, 1, x / b }, 1)
                    - G({ e / b, a / b, 1, f / b, x / b, 1 }, 1)
                    - G({ e / b, a / b, f / b, 1, 1, x / b }, 1)
                    - G({ e / b, a / b, f / b, 1, x / b, 1 }, 1)
                    - G({ e / b, a / b, f / b, x / b, 1, 1 }, 1)
                    - G({ e / b, f / b, 0, 1, 1, a / b }, 1)
                    - G({ e / b, f / b, 0, 1, a / b, 1 }, 1)
                    - G({ e / b, f / b, 0, a / b, 1, 1 }, 1)
                    - G({ e / b, f / b, 1, 0, 1, a / b }, 1)
                    - G({ e / b, f / b, 1, 0, a / b, 1 }, 1)
                    - G({ e / b, f / b, 1, 1, 0, a / b }, 1)
                    - G({ e / b, f / b, 1, 1, a / b, x / b }, 1)
                    - G({ e / b, f / b, 1, a / b, 1, x / b }, 1)
                    - G({ e / b, f / b, 1, a / b, x / b, 1 }, 1)
                    - G({ e / b, f / b, a / b, 1, 1, x / b }, 1)
                    - G({ e / b, f / b, a / b, 1, x / b, 1 }, 1)
                    - G({ e / b, f / b, a / b, x / b, 1, 1 }, 1)
                    + (-sy[9] - sy[10] - sy[5] - sy[8]) * G({ 0, f }, { 1, sf }, x)
                    + sy[2] * G({ a, 0, f }, { sa, 1, sf }, x)
                    + (sy[0] - sy[2]) * G({ a, e, f }, { sa, se, sf }, x)
                    + G({ a / b, 1 }, 1)
                        * (sy[3] - G({ a, b, e, f }, { sa, sb, se, sf }, x))
                    + G({ a / b }, 1)
                        * (-G({ 0, b, b, e, f }, { 1, sb, sb, se, sf }, x)
                            + G({ a, b, b, e, f }, { sa, sb, sb, se, sf }, x))
                    + G({ a, 0, b, b, e, f }, { sa, 1, sb, sb, se, sf }, x)
                    + (sy[14] + sy[15] + sy[16] + sy[17] + sy[18] + sy[19] + sy[20]
                          + sy[21] + sy[22] + sy[23] + sy[24] + sy[25] + sy[26] + sy[27]
                          + sy[28] + sy[13] * (-sy[11] - sy[12] - sy[7]))
                        * Log(b, sb)
                    + (-sy[14] - sy[15] - sy[16] - sy[17] - sy[18] - sy[19] - sy[20]
                          - sy[21] - sy[22] - sy[23] - sy[24] - sy[25] - sy[26] - sy[27]
                          - sy[28] + sy[13] * (sy[11] + sy[12] + sy[7]))
                        * Log(-x, sb) };
                if (e != x) {
                    res += (sy[8] - G({ a / b, x / b, 1, 1 }, 1))
                        * G({ e, f }, { se, sf }, x);
                }
                if (f != x) {
                    res += (sy[20] + sy[21] + sy[22] + sy[23] + sy[24] + sy[25] + sy[26]
                               + sy[27] + sy[28] - G({ a / b, e / b, 1, 1, x / b }, 1)
                               - G({ a / b, e / b, 1, x / b, 1 }, 1)
                               - G({ a / b, e / b, x / b, 1, 1 }, 1)
                               - G({ e / b, 1, 1, a / b, x / b }, 1)
                               - G({ e / b, 1, a / b, 1, x / b }, 1)
                               - G({ e / b, 1, a / b, x / b, 1 }, 1)
                               - G({ e / b, a / b, 1, 1, x / b }, 1)
                               - G({ e / b, a / b, 1, x / b, 1 }, 1)
                               - G({ e / b, a / b, x / b, 1, 1 }, 1))
                        * G({ f }, { sf }, x);
                }
                return res;
            }
        } else { // abbdef
            const vector<complex<double>> sy = { G({ a, e, f }, { sa, se, sf }, x),
                G({ d / b, 1, a / b }, 1), G({ a / b, d / b, 1 }, 1),
                G({ d / b, a / b, 1 }, 1), G({ d / b, 1, e / b }, 1),
                G({ d / b, e / b, 1 }, 1), G({ a, d, e, f }, { sa, sd, se, sf }, x),
                G({ b, d, e, f }, { sb, sd, se, sf }, x), G({ a, f }, { sa, sf }, x),
                G({ d / b, 1, e / b, a / b }, 1), G({ 0, a }, { 1, sa }, x),
                G({ d / b, 1, e / b, f / b }, 1), G({ a / b, d / b, 1, e / b }, 1),
                G({ a / b, d / b, e / b, 1 }, 1), G({ d / b, 1, a / b, e / b }, 1),
                G({ d / b, a / b, 1, e / b }, 1), G({ d / b, a / b, e / b, 1 }, 1),
                G({ d / b, e / b, 1, a / b }, 1), G({ d / b, e / b, a / b, 1 }, 1),
                G({ d / b, e / b, 1, f / b }, 1), G({ d / b, e / b, f / b, 1 }, 1),
                G({ a }, { sa }, x), G({ d / b, 1, e / b, f / b, a / b }, 1),
                G({ d / b, e / b, 1, f / b, a / b }, 1),
                G({ d / b, e / b, f / b, 1, a / b }, 1),
                G({ d / b, e / b, f / b, a / b, 1 }, 1),
                G({ a / b, d / b, 1, e / b, f / b }, 1),
                G({ a / b, d / b, e / b, 1, f / b }, 1),
                G({ a / b, d / b, e / b, f / b, 1 }, 1),
                G({ d / b, 1, a / b, e / b, f / b }, 1),
                G({ d / b, 1, e / b, a / b, f / b }, 1),
                G({ d / b, a / b, 1, e / b, f / b }, 1),
                G({ d / b, a / b, e / b, 1, f / b }, 1),
                G({ d / b, a / b, e / b, f / b, 1 }, 1),
                G({ d / b, e / b, 1, a / b, f / b }, 1),
                G({ d / b, e / b, a / b, 1, f / b }, 1),
                G({ d / b, e / b, a / b, f / b, 1 }, 1) };
            complex<double> res { sy[10] * sy[11] - sy[0] * sy[1]
                + sy[10] * (sy[19] + sy[20]) + sy[0] * (-sy[3] + sy[4] + sy[5])
                - sy[9] * sy[8] + (sy[11] - sy[17] - sy[18] + sy[19] + sy[20]) * sy[8]
                + (-sy[6] + sy[7]) * G({ a / b, 1 }, 1) - sy[7] * G({ a / b, x / b }, 1)
                + sy[21]
                    * (-sy[22] - sy[23] - sy[24] - sy[25]
                        - G({ 0, d / b, 1, e / b, f / b }, 1)
                        - G({ 0, d / b, e / b, 1, f / b }, 1)
                        - G({ 0, d / b, e / b, f / b, 1 }, 1)
                        - G({ d / b, 0, 1, e / b, f / b }, 1)
                        - G({ d / b, 0, e / b, 1, f / b }, 1)
                        - G({ d / b, 0, e / b, f / b, 1 }, 1)
                        - G({ d / b, 1, 0, e / b, f / b }, 1)
                        - G({ d / b, 1, e / b, 0, f / b }, 1)
                        - G({ d / b, e / b, 0, 1, f / b }, 1)
                        - G({ d / b, e / b, 0, f / b, 1 }, 1)
                        - G({ d / b, e / b, 1, 0, f / b }, 1))
                + G({ 0, a / b, d / b, 1, e / b, f / b }, 1)
                + G({ 0, a / b, d / b, e / b, 1, f / b }, 1)
                + G({ 0, a / b, d / b, e / b, f / b, 1 }, 1)
                + G({ 0, d / b, 1, a / b, e / b, f / b }, 1)
                + G({ 0, d / b, 1, e / b, a / b, f / b }, 1)
                + G({ 0, d / b, 1, e / b, f / b, a / b }, 1)
                + G({ 0, d / b, a / b, 1, e / b, f / b }, 1)
                + G({ 0, d / b, a / b, e / b, 1, f / b }, 1)
                + G({ 0, d / b, a / b, e / b, f / b, 1 }, 1)
                + G({ 0, d / b, e / b, 1, a / b, f / b }, 1)
                + G({ 0, d / b, e / b, 1, f / b, a / b }, 1)
                + G({ 0, d / b, e / b, a / b, 1, f / b }, 1)
                + G({ 0, d / b, e / b, a / b, f / b, 1 }, 1)
                + G({ 0, d / b, e / b, f / b, 1, a / b }, 1)
                + G({ 0, d / b, e / b, f / b, a / b, 1 }, 1)
                + G({ a / b, 0, d / b, 1, e / b, f / b }, 1)
                + G({ a / b, 0, d / b, e / b, 1, f / b }, 1)
                + G({ a / b, 0, d / b, e / b, f / b, 1 }, 1)
                + G({ a / b, d / b, 0, 1, e / b, f / b }, 1)
                + G({ a / b, d / b, 0, e / b, 1, f / b }, 1)
                + G({ a / b, d / b, 0, e / b, f / b, 1 }, 1)
                + G({ a / b, d / b, 1, 0, e / b, f / b }, 1)
                + G({ a / b, d / b, 1, e / b, 0, f / b }, 1)
                + G({ a / b, d / b, 1, e / b, f / b, x / b }, 1)
                + G({ a / b, d / b, e / b, 0, 1, f / b }, 1)
                + G({ a / b, d / b, e / b, 0, f / b, 1 }, 1)
                + G({ a / b, d / b, e / b, 1, 0, f / b }, 1)
                + G({ a / b, d / b, e / b, 1, f / b, x / b }, 1)
                + G({ a / b, d / b, e / b, f / b, 1, x / b }, 1)
                + G({ a / b, d / b, e / b, f / b, x / b, 1 }, 1)
                + G({ d / b, 0, 1, a / b, e / b, f / b }, 1)
                + G({ d / b, 0, 1, e / b, a / b, f / b }, 1)
                + G({ d / b, 0, 1, e / b, f / b, a / b }, 1)
                + G({ d / b, 0, a / b, 1, e / b, f / b }, 1)
                + G({ d / b, 0, a / b, e / b, 1, f / b }, 1)
                + G({ d / b, 0, a / b, e / b, f / b, 1 }, 1)
                + G({ d / b, 0, e / b, 1, a / b, f / b }, 1)
                + G({ d / b, 0, e / b, 1, f / b, a / b }, 1)
                + G({ d / b, 0, e / b, a / b, 1, f / b }, 1)
                + G({ d / b, 0, e / b, a / b, f / b, 1 }, 1)
                + G({ d / b, 0, e / b, f / b, 1, a / b }, 1)
                + G({ d / b, 0, e / b, f / b, a / b, 1 }, 1)
                + G({ d / b, 1, 0, a / b, e / b, f / b }, 1)
                + G({ d / b, 1, 0, e / b, a / b, f / b }, 1)
                + G({ d / b, 1, 0, e / b, f / b, a / b }, 1)
                + G({ d / b, 1, a / b, 0, e / b, f / b }, 1)
                + G({ d / b, 1, a / b, e / b, 0, f / b }, 1)
                + G({ d / b, 1, a / b, e / b, f / b, x / b }, 1)
                + G({ d / b, 1, e / b, 0, a / b, f / b }, 1)
                + G({ d / b, 1, e / b, 0, f / b, a / b }, 1)
                + G({ d / b, 1, e / b, a / b, 0, f / b }, 1)
                + G({ d / b, 1, e / b, a / b, f / b, x / b }, 1)
                + G({ d / b, 1, e / b, f / b, 0, a / b }, 1)
                + G({ d / b, 1, e / b, f / b, a / b, x / b }, 1)
                + G({ d / b, a / b, 0, 1, e / b, f / b }, 1)
                + G({ d / b, a / b, 0, e / b, 1, f / b }, 1)
                + G({ d / b, a / b, 0, e / b, f / b, 1 }, 1)
                + G({ d / b, a / b, 1, 0, e / b, f / b }, 1)
                + G({ d / b, a / b, 1, e / b, 0, f / b }, 1)
                + G({ d / b, a / b, 1, e / b, f / b, x / b }, 1)
                + G({ d / b, a / b, e / b, 0, 1, f / b }, 1)
                + G({ d / b, a / b, e / b, 0, f / b, 1 }, 1)
                + G({ d / b, a / b, e / b, 1, 0, f / b }, 1)
                + G({ d / b, a / b, e / b, 1, f / b, x / b }, 1)
                + G({ d / b, a / b, e / b, f / b, 1, x / b }, 1)
                + G({ d / b, a / b, e / b, f / b, x / b, 1 }, 1)
                + G({ d / b, e / b, 0, 1, a / b, f / b }, 1)
                + G({ d / b, e / b, 0, 1, f / b, a / b }, 1)
                + G({ d / b, e / b, 0, a / b, 1, f / b }, 1)
                + G({ d / b, e / b, 0, a / b, f / b, 1 }, 1)
                + G({ d / b, e / b, 0, f / b, 1, a / b }, 1)
                + G({ d / b, e / b, 0, f / b, a / b, 1 }, 1)
                + G({ d / b, e / b, 1, 0, a / b, f / b }, 1)
                + G({ d / b, e / b, 1, 0, f / b, a / b }, 1)
                + G({ d / b, e / b, 1, a / b, 0, f / b }, 1)
                + G({ d / b, e / b, 1, a / b, f / b, x / b }, 1)
                + G({ d / b, e / b, 1, f / b, 0, a / b }, 1)
                + G({ d / b, e / b, 1, f / b, a / b, x / b }, 1)
                + G({ d / b, e / b, a / b, 0, 1, f / b }, 1)
                + G({ d / b, e / b, a / b, 0, f / b, 1 }, 1)
                + G({ d / b, e / b, a / b, 1, 0, f / b }, 1)
                + G({ d / b, e / b, a / b, 1, f / b, x / b }, 1)
                + G({ d / b, e / b, a / b, f / b, 1, x / b }, 1)
                + G({ d / b, e / b, a / b, f / b, x / b, 1 }, 1)
                + G({ d / b, e / b, f / b, 0, 1, a / b }, 1)
                + G({ d / b, e / b, f / b, 0, a / b, 1 }, 1)
                + G({ d / b, e / b, f / b, 1, 0, a / b }, 1)
                + G({ d / b, e / b, f / b, 1, a / b, x / b }, 1)
                + G({ d / b, e / b, f / b, a / b, 1, x / b }, 1)
                + G({ d / b, e / b, f / b, a / b, x / b, 1 }, 1)
                + (sy[9] + sy[12] + sy[13] + sy[14] + sy[15] + sy[16] + sy[17] + sy[18])
                    * G({ 0, f }, { 1, sf }, x)
                + (sy[1] + sy[2] + sy[3]) * G({ 0, e, f }, { 1, se, sf }, x)
                + (-sy[4] - sy[5]) * G({ a, 0, f }, { sa, 1, sf }, x)
                + G({ d / b, 1 }, 1) * (sy[6] - G({ a, 0, e, f }, { sa, 1, se, sf }, x))
                + G({ a / b }, 1)
                    * (-G({ 0, b, d, e, f }, { 1, sb, sd, se, sf }, x)
                        + G({ a, b, d, e, f }, { sa, sb, sd, se, sf }, x))
                + G({ a, 0, b, d, e, f }, { sa, 1, sb, sd, se, sf }, x)
                + ((sy[11] + sy[19] + sy[20]) * sy[21] - sy[22] - sy[23] - sy[24] - sy[25]
                      - sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31] - sy[32]
                      - sy[33] - sy[34] - sy[35] - sy[36])
                    * Log(b, sb)
                + ((-sy[11] - sy[19] - sy[20]) * sy[21] + sy[22] + sy[23] + sy[24]
                      + sy[25] + sy[26] + sy[27] + sy[28] + sy[29] + sy[30] + sy[31]
                      + sy[32] + sy[33] + sy[34] + sy[35] + sy[36])
                    * Log(-x, sb) };
            if (d != x) {
                res += (-sy[2] + G({ a / b, x / b, 1 }, 1))
                    * G({ d, e, f }, { sd, se, sf }, x);
            }
            if (e != x) {
                res += (-sy[12] - sy[13] - sy[14] - sy[15] - sy[16]
                           + G({ a / b, d / b, 1, x / b }, 1)
                           + G({ a / b, d / b, x / b, 1 }, 1)
                           + G({ d / b, 1, a / b, x / b }, 1)
                           + G({ d / b, a / b, 1, x / b }, 1)
                           + G({ d / b, a / b, x / b, 1 }, 1))
                    * G({ e, f }, { se, sf }, x);
            }
            if (f != x) {
                res += (-sy[26] - sy[27] - sy[28] - sy[29] - sy[30] - sy[31] - sy[32]
                           - sy[33] - sy[34] - sy[35] - sy[36]
                           + G({ a / b, d / b, 1, e / b, x / b }, 1)
                           + G({ a / b, d / b, e / b, 1, x / b }, 1)
                           + G({ a / b, d / b, e / b, x / b, 1 }, 1)
                           + G({ d / b, 1, a / b, e / b, x / b }, 1)
                           + G({ d / b, 1, e / b, a / b, x / b }, 1)
                           + G({ d / b, a / b, 1, e / b, x / b }, 1)
                           + G({ d / b, a / b, e / b, 1, x / b }, 1)
                           + G({ d / b, a / b, e / b, x / b, 1 }, 1)
                           + G({ d / b, e / b, 1, a / b, x / b }, 1)
                           + G({ d / b, e / b, a / b, 1, x / b }, 1)
                           + G({ d / b, e / b, a / b, x / b, 1 }, 1))
                    * G({ f }, { sf }, x);
            }
            return res;
        }
    } else { // abcdef
        const vector<complex<double>> sy = { G({ a / b, c / b, d / b }, 1),
            G({ c / b, a / b, d / b }, 1), G({ c / b, d / b, a / b }, 1),
            G({ a, e, f }, { sa, se, sf }, x), G({ c / b, d / b, e / b }, 1),
            G({ a / b, c / b }, 1), G({ 0, d, e, f }, { 1, sd, se, sf }, x),
            G({ c / b, a / b }, 1), G({ a, d, e, f }, { sa, sd, se, sf }, x),
            G({ a / b, c / b, d / b, e / b }, 1), G({ c / b, a / b, d / b, e / b }, 1),
            G({ c / b, d / b, a / b, e / b }, 1), G({ c / b, d / b, e / b, a / b }, 1),
            G({ a, f }, { sa, sf }, x), G({ c / b, d / b, e / b, f / b }, 1),
            G({ a, c, d, e, f }, { sa, sc, sd, se, sf }, x), G({ a }, { sa }, x),
            G({ c / b, d / b, e / b, f / b, a / b }, 1),
            G({ a / b, c / b, d / b, e / b, f / b }, 1),
            G({ c / b, a / b, d / b, e / b, f / b }, 1),
            G({ c / b, d / b, a / b, e / b, f / b }, 1),
            G({ c / b, d / b, e / b, a / b, f / b }, 1) };
        complex<double> res { sy[12] * sy[13] - sy[13] * sy[14] + sy[2] * sy[3]
            - sy[3] * sy[4] - sy[5] * sy[6] - sy[6] * sy[7] + sy[7] * sy[8]
            + sy[16]
                * (sy[17] + G({ 0, c / b, d / b, e / b, f / b }, 1)
                    + G({ c / b, 0, d / b, e / b, f / b }, 1)
                    + G({ c / b, d / b, 0, e / b, f / b }, 1)
                    + G({ c / b, d / b, e / b, 0, f / b }, 1))
            - G({ 0, a / b, c / b, d / b, e / b, f / b }, 1)
            - G({ 0, c / b, a / b, d / b, e / b, f / b }, 1)
            - G({ 0, c / b, d / b, a / b, e / b, f / b }, 1)
            - G({ 0, c / b, d / b, e / b, a / b, f / b }, 1)
            - G({ 0, c / b, d / b, e / b, f / b, a / b }, 1)
            - G({ a / b, 0, c / b, d / b, e / b, f / b }, 1)
            - G({ a / b, c / b, 0, d / b, e / b, f / b }, 1)
            - G({ a / b, c / b, d / b, 0, e / b, f / b }, 1)
            - G({ a / b, c / b, d / b, e / b, 0, f / b }, 1)
            - G({ a / b, c / b, d / b, e / b, f / b, x / b }, 1)
            - G({ c / b, 0, a / b, d / b, e / b, f / b }, 1)
            - G({ c / b, 0, d / b, a / b, e / b, f / b }, 1)
            - G({ c / b, 0, d / b, e / b, a / b, f / b }, 1)
            - G({ c / b, 0, d / b, e / b, f / b, a / b }, 1)
            - G({ c / b, a / b, 0, d / b, e / b, f / b }, 1)
            - G({ c / b, a / b, d / b, 0, e / b, f / b }, 1)
            - G({ c / b, a / b, d / b, e / b, 0, f / b }, 1)
            - G({ c / b, a / b, d / b, e / b, f / b, x / b }, 1)
            - G({ c / b, d / b, 0, a / b, e / b, f / b }, 1)
            - G({ c / b, d / b, 0, e / b, a / b, f / b }, 1)
            - G({ c / b, d / b, 0, e / b, f / b, a / b }, 1)
            - G({ c / b, d / b, a / b, 0, e / b, f / b }, 1)
            - G({ c / b, d / b, a / b, e / b, 0, f / b }, 1)
            - G({ c / b, d / b, a / b, e / b, f / b, x / b }, 1)
            - G({ c / b, d / b, e / b, 0, a / b, f / b }, 1)
            - G({ c / b, d / b, e / b, 0, f / b, a / b }, 1)
            - G({ c / b, d / b, e / b, a / b, 0, f / b }, 1)
            - G({ c / b, d / b, e / b, a / b, f / b, x / b }, 1)
            - G({ c / b, d / b, e / b, f / b, 0, a / b }, 1)
            - G({ c / b, d / b, e / b, f / b, a / b, x / b }, 1)
            - sy[14] * G({ 0, a }, { 1, sa }, x)
            + (-sy[9] - sy[10] - sy[11] - sy[12]) * G({ 0, f }, { 1, sf }, x)
            + (-sy[0] - sy[1] - sy[2]) * G({ 0, e, f }, { 1, se, sf }, x)
            + sy[4] * G({ a, 0, f }, { sa, 1, sf }, x)
            + G({ c / b, d / b }, 1) * (-sy[8] + G({ a, 0, e, f }, { sa, 1, se, sf }, x))
            + G({ a / b }, 1) * (sy[15] - G({ 0, c, d, e, f }, { 1, sc, sd, se, sf }, x))
            + G({ c / b }, 1) * (-sy[15] + G({ a, 0, d, e, f }, { sa, 1, sd, se, sf }, x))
            + G({ a, 0, c, d, e, f }, { sa, 1, sc, sd, se, sf }, x)
            + (-(sy[14] * sy[16]) + sy[17] + sy[18] + sy[19] + sy[20] + sy[21])
                * Log(b, sb)
            + (sy[14] * sy[16] - sy[17] - sy[18] - sy[19] - sy[20] - sy[21])
                * Log(-x, sb) };
        if (c != x) {
            res += (sy[5] - G({ a / b, x / b }, 1))
                * G({ c, d, e, f }, { sc, sd, se, sf }, x);
        }
        if (d != x) {
            res += (sy[0] + sy[1] - G({ a / b, c / b, x / b }, 1)
                       - G({ c / b, a / b, x / b }, 1))
                * G({ d, e, f }, { sd, se, sf }, x);
        }
        if (e != x) {
            res += (sy[9] + sy[10] + sy[11] - G({ a / b, c / b, d / b, x / b }, 1)
                       - G({ c / b, a / b, d / b, x / b }, 1)
                       - G({ c / b, d / b, a / b, x / b }, 1))
                * G({ e, f }, { se, sf }, x);
        }
        if (f != x) {
            res += (sy[18] + sy[19] + sy[20] + sy[21]
                       - G({ a / b, c / b, d / b, e / b, x / b }, 1)
                       - G({ c / b, a / b, d / b, e / b, x / b }, 1)
                       - G({ c / b, d / b, a / b, e / b, x / b }, 1)
                       - G({ c / b, d / b, e / b, a / b, x / b }, 1))
                * G({ f }, { sf }, x);
        }
        return res;
    }
}
complex<double> G6_abcdef_a(complex<double> a1, complex<double> b1, complex<double> c1,
    complex<double> d1, complex<double> e1, complex<double> f1, int sa, int sb, int sc,
    int sd, int se, int sf, double x1)
{
    complex<double> a = a1 / x1;
    complex<double> b = b1 / x1;
    complex<double> c = c1 / x1;
    complex<double> d = d1 / x1;
    complex<double> e = e1 / x1;
    complex<double> f = f1 / x1;
    double x = 1.;
    if (a == b) {
        if (a == c) {
            if (a == d) {
                if (a == e) { // aaaaaf
                    const vector<complex<double>> sy = { G({ f / a, 1, 1, 1, 1 }, 1) };
                    complex<double> res { G({ 0, f / a, 1, 1, 1, 1 }, 1)
                        + G({ f / a, 1, 1, 1, 1, x / a }, 1)
                        + G({ f / a, 1, 1, 1, x / a, 1 }, 1)
                        + G({ f / a, 1, 1, x / a, 1, 1 }, 1)
                        + G({ f / a, 1, x / a, 1, 1, 1 }, 1)
                        + G({ f / a, x / a, 1, 1, 1, 1 }, 1)
                        - G({ x / a, 1, 1, 1 }, 1) * G({ a, f }, { sa, sf }, x)
                        + G({ x / a, 1, 1 }, 1) * G({ a, a, f }, { sa, sa, sf }, x)
                        - G({ x / a, 1 }, 1) * G({ a, a, a, f }, { sa, sa, sa, sf }, x)
                        + G({ x / a }, 1)
                            * G({ a, a, a, a, f }, { sa, sa, sa, sa, sf }, x)
                        + G({ 0, a, a, a, a, f }, { 1, sa, sa, sa, sa, sf }, x)
                        - sy[0] * Log(a, sa) + sy[0] * Log(-x, sa) };
                    if (f != x) {
                        res += (-sy[0] + G({ x / a, 1, 1, 1, 1 }, 1))
                            * G({ f }, { sf }, x);
                    }
                    return res;
                } else { // aaaaef
                    const vector<complex<double>> sy = { G({ e / a, 1, 1, 1 }, 1),
                        G({ e / a, 1, 1, 1, f / a }, 1), G({ e / a, 1, 1, f / a, 1 }, 1),
                        G({ e / a, 1, f / a, 1, 1 }, 1),
                        G({ e / a, f / a, 1, 1, 1 }, 1) };
                    complex<double> res { -G({ 0, e / a, 1, 1, 1, f / a }, 1)
                        - G({ 0, e / a, 1, 1, f / a, 1 }, 1)
                        - G({ 0, e / a, 1, f / a, 1, 1 }, 1)
                        - G({ 0, e / a, f / a, 1, 1, 1 }, 1)
                        - G({ e / a, 0, 1, 1, 1, f / a }, 1)
                        - G({ e / a, 0, 1, 1, f / a, 1 }, 1)
                        - G({ e / a, 0, 1, f / a, 1, 1 }, 1)
                        - G({ e / a, 0, f / a, 1, 1, 1 }, 1)
                        - G({ e / a, 1, 0, 1, 1, f / a }, 1)
                        - G({ e / a, 1, 0, 1, f / a, 1 }, 1)
                        - G({ e / a, 1, 0, f / a, 1, 1 }, 1)
                        - G({ e / a, 1, 1, 0, 1, f / a }, 1)
                        - G({ e / a, 1, 1, 0, f / a, 1 }, 1)
                        - G({ e / a, 1, 1, 1, 0, f / a }, 1)
                        - G({ e / a, 1, 1, 1, f / a, x / a }, 1)
                        - G({ e / a, 1, 1, f / a, 1, x / a }, 1)
                        - G({ e / a, 1, 1, f / a, x / a, 1 }, 1)
                        - G({ e / a, 1, f / a, 1, 1, x / a }, 1)
                        - G({ e / a, 1, f / a, 1, x / a, 1 }, 1)
                        - G({ e / a, 1, f / a, x / a, 1, 1 }, 1)
                        - G({ e / a, f / a, 1, 1, 1, x / a }, 1)
                        - G({ e / a, f / a, 1, 1, x / a, 1 }, 1)
                        - G({ e / a, f / a, 1, x / a, 1, 1 }, 1)
                        - G({ e / a, f / a, x / a, 1, 1, 1 }, 1)
                        - sy[0] * G({ 0, f }, { 1, sf }, x)
                        + G({ x / a, 1, 1 }, 1) * G({ a, e, f }, { sa, se, sf }, x)
                        - G({ x / a, 1 }, 1) * G({ a, a, e, f }, { sa, sa, se, sf }, x)
                        + G({ x / a }, 1)
                            * G({ a, a, a, e, f }, { sa, sa, sa, se, sf }, x)
                        + G({ 0, a, a, a, e, f }, { 1, sa, sa, sa, se, sf }, x)
                        + (sy[1] + sy[2] + sy[3] + sy[4]) * Log(a, sa)
                        + (-sy[1] - sy[2] - sy[3] - sy[4]) * Log(-x, sa) };
                    if (e != x) {
                        res += (sy[0] - G({ x / a, 1, 1, 1 }, 1))
                            * G({ e, f }, { se, sf }, x);
                    }
                    if (f != x) {
                        res += (sy[1] + sy[2] + sy[3] + sy[4]
                                   - G({ e / a, 1, 1, 1, x / a }, 1)
                                   - G({ e / a, 1, 1, x / a, 1 }, 1)
                                   - G({ e / a, 1, x / a, 1, 1 }, 1)
                                   - G({ e / a, x / a, 1, 1, 1 }, 1))
                            * G({ f }, { sf }, x);
                    }
                    return res;
                }
            } else { // aaadef
                const vector<complex<double>> sy = { G({ d / a, 1, 1 }, 1),
                    G({ d / a, 1, 1, e / a }, 1), G({ d / a, 1, e / a, 1 }, 1),
                    G({ d / a, e / a, 1, 1 }, 1), G({ d / a, 1, 1, e / a, f / a }, 1),
                    G({ d / a, 1, e / a, 1, f / a }, 1),
                    G({ d / a, 1, e / a, f / a, 1 }, 1),
                    G({ d / a, e / a, 1, 1, f / a }, 1),
                    G({ d / a, e / a, 1, f / a, 1 }, 1),
                    G({ d / a, e / a, f / a, 1, 1 }, 1) };
                complex<double> res { G({ 0, d / a, 1, 1, e / a, f / a }, 1)
                    + G({ 0, d / a, 1, e / a, 1, f / a }, 1)
                    + G({ 0, d / a, 1, e / a, f / a, 1 }, 1)
                    + G({ 0, d / a, e / a, 1, 1, f / a }, 1)
                    + G({ 0, d / a, e / a, 1, f / a, 1 }, 1)
                    + G({ 0, d / a, e / a, f / a, 1, 1 }, 1)
                    + G({ d / a, 0, 1, 1, e / a, f / a }, 1)
                    + G({ d / a, 0, 1, e / a, 1, f / a }, 1)
                    + G({ d / a, 0, 1, e / a, f / a, 1 }, 1)
                    + G({ d / a, 0, e / a, 1, 1, f / a }, 1)
                    + G({ d / a, 0, e / a, 1, f / a, 1 }, 1)
                    + G({ d / a, 0, e / a, f / a, 1, 1 }, 1)
                    + G({ d / a, 1, 0, 1, e / a, f / a }, 1)
                    + G({ d / a, 1, 0, e / a, 1, f / a }, 1)
                    + G({ d / a, 1, 0, e / a, f / a, 1 }, 1)
                    + G({ d / a, 1, 1, 0, e / a, f / a }, 1)
                    + G({ d / a, 1, 1, e / a, 0, f / a }, 1)
                    + G({ d / a, 1, 1, e / a, f / a, x / a }, 1)
                    + G({ d / a, 1, e / a, 0, 1, f / a }, 1)
                    + G({ d / a, 1, e / a, 0, f / a, 1 }, 1)
                    + G({ d / a, 1, e / a, 1, 0, f / a }, 1)
                    + G({ d / a, 1, e / a, 1, f / a, x / a }, 1)
                    + G({ d / a, 1, e / a, f / a, 1, x / a }, 1)
                    + G({ d / a, 1, e / a, f / a, x / a, 1 }, 1)
                    + G({ d / a, e / a, 0, 1, 1, f / a }, 1)
                    + G({ d / a, e / a, 0, 1, f / a, 1 }, 1)
                    + G({ d / a, e / a, 0, f / a, 1, 1 }, 1)
                    + G({ d / a, e / a, 1, 0, 1, f / a }, 1)
                    + G({ d / a, e / a, 1, 0, f / a, 1 }, 1)
                    + G({ d / a, e / a, 1, 1, 0, f / a }, 1)
                    + G({ d / a, e / a, 1, 1, f / a, x / a }, 1)
                    + G({ d / a, e / a, 1, f / a, 1, x / a }, 1)
                    + G({ d / a, e / a, 1, f / a, x / a, 1 }, 1)
                    + G({ d / a, e / a, f / a, 1, 1, x / a }, 1)
                    + G({ d / a, e / a, f / a, 1, x / a, 1 }, 1)
                    + G({ d / a, e / a, f / a, x / a, 1, 1 }, 1)
                    + (sy[1] + sy[2] + sy[3]) * G({ 0, f }, { 1, sf }, x)
                    + sy[0] * G({ 0, e, f }, { 1, se, sf }, x)
                    - G({ x / a, 1 }, 1) * G({ a, d, e, f }, { sa, sd, se, sf }, x)
                    + G({ x / a }, 1) * G({ a, a, d, e, f }, { sa, sa, sd, se, sf }, x)
                    + G({ 0, a, a, d, e, f }, { 1, sa, sa, sd, se, sf }, x)
                    + (-sy[9] - sy[4] - sy[5] - sy[6] - sy[7] - sy[8]) * Log(a, sa)
                    + (sy[9] + sy[4] + sy[5] + sy[6] + sy[7] + sy[8]) * Log(-x, sa) };
                if (d != x) {
                    res += (-sy[0] + G({ x / a, 1, 1 }, 1))
                        * G({ d, e, f }, { sd, se, sf }, x);
                }
                if (e != x) {
                    res += (-sy[1] - sy[2] - sy[3] + G({ d / a, 1, 1, x / a }, 1)
                               + G({ d / a, 1, x / a, 1 }, 1)
                               + G({ d / a, x / a, 1, 1 }, 1))
                        * G({ e, f }, { se, sf }, x);
                }
                if (f != x) {
                    res += (-sy[9] - sy[4] - sy[5] - sy[6] - sy[7] - sy[8]
                               + G({ d / a, 1, 1, e / a, x / a }, 1)
                               + G({ d / a, 1, e / a, 1, x / a }, 1)
                               + G({ d / a, 1, e / a, x / a, 1 }, 1)
                               + G({ d / a, e / a, 1, 1, x / a }, 1)
                               + G({ d / a, e / a, 1, x / a, 1 }, 1)
                               + G({ d / a, e / a, x / a, 1, 1 }, 1))
                        * G({ f }, { sf }, x);
                }
                return res;
            }
        } else { // aacdef
            const vector<complex<double>> sy = { G({ c / a, 1, d / a }, 1),
                G({ c / a, d / a, 1 }, 1), G({ c / a, 1 }, 1),
                G({ c / a, 1, d / a, e / a }, 1), G({ c / a, d / a, 1, e / a }, 1),
                G({ c / a, d / a, e / a, 1 }, 1), G({ c / a, 1, d / a, e / a, f / a }, 1),
                G({ c / a, d / a, 1, e / a, f / a }, 1),
                G({ c / a, d / a, e / a, 1, f / a }, 1),
                G({ c / a, d / a, e / a, f / a, 1 }, 1) };
            complex<double> res { -G({ 0, c / a, 1, d / a, e / a, f / a }, 1)
                - G({ 0, c / a, d / a, 1, e / a, f / a }, 1)
                - G({ 0, c / a, d / a, e / a, 1, f / a }, 1)
                - G({ 0, c / a, d / a, e / a, f / a, 1 }, 1)
                - G({ c / a, 0, 1, d / a, e / a, f / a }, 1)
                - G({ c / a, 0, d / a, 1, e / a, f / a }, 1)
                - G({ c / a, 0, d / a, e / a, 1, f / a }, 1)
                - G({ c / a, 0, d / a, e / a, f / a, 1 }, 1)
                - G({ c / a, 1, 0, d / a, e / a, f / a }, 1)
                - G({ c / a, 1, d / a, 0, e / a, f / a }, 1)
                - G({ c / a, 1, d / a, e / a, 0, f / a }, 1)
                - G({ c / a, 1, d / a, e / a, f / a, x / a }, 1)
                - G({ c / a, d / a, 0, 1, e / a, f / a }, 1)
                - G({ c / a, d / a, 0, e / a, 1, f / a }, 1)
                - G({ c / a, d / a, 0, e / a, f / a, 1 }, 1)
                - G({ c / a, d / a, 1, 0, e / a, f / a }, 1)
                - G({ c / a, d / a, 1, e / a, 0, f / a }, 1)
                - G({ c / a, d / a, 1, e / a, f / a, x / a }, 1)
                - G({ c / a, d / a, e / a, 0, 1, f / a }, 1)
                - G({ c / a, d / a, e / a, 0, f / a, 1 }, 1)
                - G({ c / a, d / a, e / a, 1, 0, f / a }, 1)
                - G({ c / a, d / a, e / a, 1, f / a, x / a }, 1)
                - G({ c / a, d / a, e / a, f / a, 1, x / a }, 1)
                - G({ c / a, d / a, e / a, f / a, x / a, 1 }, 1)
                + (-sy[3] - sy[4] - sy[5]) * G({ 0, f }, { 1, sf }, x)
                + (-sy[0] - sy[1]) * G({ 0, e, f }, { 1, se, sf }, x)
                - sy[2] * G({ 0, d, e, f }, { 1, sd, se, sf }, x)
                + G({ x / a }, 1) * G({ a, c, d, e, f }, { sa, sc, sd, se, sf }, x)
                + G({ 0, a, c, d, e, f }, { 1, sa, sc, sd, se, sf }, x)
                + (sy[9] + sy[6] + sy[7] + sy[8]) * Log(a, sa)
                + (-sy[9] - sy[6] - sy[7] - sy[8]) * Log(-x, sa) };
            if (c != x) {
                res += (sy[2] - G({ x / a, 1 }, 1))
                    * G({ c, d, e, f }, { sc, sd, se, sf }, x);
            }
            if (d != x) {
                res += (sy[0] + sy[1] - G({ c / a, 1, x / a }, 1)
                           - G({ c / a, x / a, 1 }, 1))
                    * G({ d, e, f }, { sd, se, sf }, x);
            }
            if (e != x) {
                res += (sy[3] + sy[4] + sy[5] - G({ c / a, 1, d / a, x / a }, 1)
                           - G({ c / a, d / a, 1, x / a }, 1)
                           - G({ c / a, d / a, x / a, 1 }, 1))
                    * G({ e, f }, { se, sf }, x);
            }
            if (f != x) {
                res += (sy[9] + sy[6] + sy[7] + sy[8]
                           - G({ c / a, 1, d / a, e / a, x / a }, 1)
                           - G({ c / a, d / a, 1, e / a, x / a }, 1)
                           - G({ c / a, d / a, e / a, 1, x / a }, 1)
                           - G({ c / a, d / a, e / a, x / a, 1 }, 1))
                    * G({ f }, { sf }, x);
            }
            return res;
        }
    } else { // abcdef
        const vector<complex<double>> sy = { G({ b / a, c / a, d / a }, 1),
            G({ b / a, c / a }, 1), G({ b / a, c / a, d / a, e / a }, 1), G({ b / a }, 1),
            G({ b / a, c / a, d / a, e / a, f / a }, 1) };
        complex<double> res { G({ 0, b / a, c / a, d / a, e / a, f / a }, 1)
            + G({ b / a, 0, c / a, d / a, e / a, f / a }, 1)
            + G({ b / a, c / a, 0, d / a, e / a, f / a }, 1)
            + G({ b / a, c / a, d / a, 0, e / a, f / a }, 1)
            + G({ b / a, c / a, d / a, e / a, 0, f / a }, 1)
            + G({ b / a, c / a, d / a, e / a, f / a, x / a }, 1)
            + sy[2] * G({ 0, f }, { 1, sf }, x) + sy[0] * G({ 0, e, f }, { 1, se, sf }, x)
            + sy[1] * G({ 0, d, e, f }, { 1, sd, se, sf }, x)
            + sy[3] * G({ 0, c, d, e, f }, { 1, sc, sd, se, sf }, x)
            + G({ 0, b, c, d, e, f }, { 1, sb, sc, sd, se, sf }, x) - sy[4] * Log(a, sa)
            + sy[4] * Log(-x, sa) };
        if (b != x) {
            res += (-sy[3] + G({ x / a }, 1))
                * G({ b, c, d, e, f }, { sb, sc, sd, se, sf }, x);
        }
        if (c != x) {
            res += (-sy[1] + G({ b / a, x / a }, 1))
                * G({ c, d, e, f }, { sc, sd, se, sf }, x);
        }
        if (d != x) {
            res += (-sy[0] + G({ b / a, c / a, x / a }, 1))
                * G({ d, e, f }, { sd, se, sf }, x);
        }
        if (e != x) {
            res += (-sy[2] + G({ b / a, c / a, d / a, x / a }, 1))
                * G({ e, f }, { se, sf }, x);
        }
        if (f != x) {
            res += (-sy[4] + G({ b / a, c / a, d / a, e / a, x / a }, 1))
                * G({ f }, { sf }, x);
        }
        return res;
    }
}

complex<double> G6_abcdef(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, complex<double> f, int sa, int sb, int sc,
    int sd, int se, int sf, double x)
{
    // e is smallest
    if (abs(f) < abs(e) && abs(f) < abs(d) && abs(f) < abs(c) && abs(f) < abs(b)
        && abs(f) < abs(a))
        return G6_abcdef_f(a, b, c, d, e, f, sa, sb, sc, sd, se, sf, x);

    // e is smallest
    if (abs(e) < abs(d) && abs(e) < abs(c) && abs(e) < abs(b) && abs(e) < abs(a))
        return G6_abcdef_e(a, b, c, d, e, f, sa, sb, sc, sd, se, sf, x);

    // d is smallest
    if (abs(d) < abs(c) && abs(d) < abs(b) && abs(d) < abs(a))
        return G6_abcdef_d(a, b, c, d, e, f, sa, sb, sc, sd, se, sf, x);

    // c is smallest
    if (abs(c) < abs(b) && abs(c) < abs(a))
        return G6_abcdef_c(a, b, c, d, e, f, sa, sb, sc, sd, se, sf, x);

    // a is smallest
    if (abs(a) <= abs(b))
        return G6_abcdef_a(a, b, c, d, e, f, sa, sb, sc, sd, se, sf, x);

    // b is smallest
    return G6_abcdef_b(a, b, c, d, e, f, sa, sb, sc, sd, se, sf, x);
}
