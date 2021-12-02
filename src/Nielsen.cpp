#include <array>
#include <cmath>
#include <complex>

#include "FastGPL.h"

#include "FastGPL_internal.h"

using std::array;
using std::complex;
using std::string;

using FastGPL::FastGPL_error;
using FastGPL::Log;
using FastGPL::PolyLog;

using FastGPL_internal::Zeta;

complex<double> FastGPL::S22(complex<double> x)
{
    int s { -1 };
    if (x.imag() > 0)
        s = 1;
    return S22(x, s);
}

complex<double> FastGPL::S22(complex<double> x, int s)
{
    if (x == 0.0)
        return 0;

    if (x == 1.0)
        return Zeta(4) / 4.0;

    if (x.real() <= 0.5 && abs(x) <= 1.0) {
        return FastGPL_internal::S22_series2(x);
    } else if (x.real() > 0.5 && abs(x - 1.0) <= 1.0) {
        return FastGPL_internal::S22_series3(x, s);
    } else {
        return FastGPL_internal::S22_inversion(x, s);
    }
}

// Series expansion of S22[1-Exp[-alpha]]
// Applicable when Re[x]<=1/2 && Abs[x]<=1
complex<double> FastGPL_internal::S22_series2(complex<double> x)
{
    constexpr array<double, 40> C { { 0.125, -0.0694444444444444444,
        0.0182291666666666667, -0.00243055555555555556, 0.0000675154320987654321,
        0.0000248015873015873016, -1.72233245149911817e-6, -4.21014599255339996e-7,
        4.29434891240446796e-8, 8.1419352252685586e-9, -1.05543604771999834e-9,
        -1.67281706633558485e-10, 2.5744345773117864e-11, 3.56473144276537239e-12,
        -6.25827037844391082e-13, -7.79007947558238045e-14, 1.51973086355922995e-14,
        1.7344855317440117e-15, -3.69122553744208134e-16, -3.91824260504742907e-17,
        8.97337227588407234e-18, 8.95449892357568294e-19, -2.18406852756017018e-19,
        -2.06587606063508607e-20, 5.32311626792523144e-21, 4.80388741558123429e-22,
        -1.29919307887332759e-22, -1.12454760086733879e-23, 3.1752931223297456e-24,
        2.64754547573690901e-25, -7.77100909127361322e-26, -6.26407321905136578e-27,
        1.90427148194954478e-27, 1.48849618171058467e-28, -4.67207275929987443e-29,
        -3.55053344275398719e-30, 1.14759930183546983e-30, 8.49785158385502234e-32,
        -2.8219025185892045e-32 } };

    complex<double> alpha { -log(1. - x) };

    size_t j { 0 };
    complex<double> power { alpha * alpha };
    complex<double> term { power * C[0] };
    complex<double> res { term };
    while (abs(term) > abs(res) * 1e-16) {
        ++j;
        power *= alpha;
        if (j == C.size())
            throw FastGPL_error { "S22_series2 converges too slowly" };
        if (C[j] == 0)
            continue;
        term = power * C[j];
        res += term;
    }
    return res;
}

// Series expansion of S22[Exp[-alpha]]
// Applicable when Re[x]>1/2 && Abs[x-1]<=1
complex<double> FastGPL_internal::S22_series3(complex<double> x, int s)
{
    complex<double> alpha { -log(x) };

    complex<double> Lalpha { Log(alpha, -s) };

    complex<double> res { Zeta(4) / 4. - alpha * Zeta(3)
        + alpha * alpha * Lalpha * Lalpha / 4. };

    constexpr array<double, 40> C { { -0.75, -0.0833333333333333333,
        0.00347222222222222222, 0, -0.0000115740740740740741, 0, 9.8418997228521038e-8, 0,
        -1.14822163433274544e-9, 0, 1.58157249908091659e-11, 0, -2.41950097925251519e-13,
        0, 3.98289777698948775e-15, 0, -6.92336661830592906e-17, 0,
        1.25527223044997728e-18, 0, -2.35375400276846523e-20, 0, 4.53639890345868702e-22,
        0, -8.94516967039264317e-24, 0, 1.79828400469549627e-25, 0,
        -3.67549976479373844e-27, 0, 7.62080797156479523e-29, 0, -1.60004196436948598e-30,
        0, 3.39676114756037559e-32, 0, -7.2822722867577647e-34, 0,
        1.57502264795800349e-35, 0 } };

    constexpr array<double, 40> D { { 0.875, 0.0694444444444444444, 0.0083912037037037037,
        -0.00104166666666666667, 0.0000331790123456790123, 4.1335978835978836e-6,
        -2.84712099125364431e-7, -3.82740544777581815e-8, 3.46380193023711542e-9,
        4.69727032227032227e-10, -4.9873450965335722e-11, -6.69126826534233942e-12,
        7.94446060563646951e-13, 1.04845042434275658e-13, -1.35529160467003403e-14,
        -1.75716078396595048e-15, 2.43097318868525724e-16, 3.09729559240002089e-17,
        -4.53165832274831275e-18, -5.67861247108323053e-19, 8.71027688189545228e-20,
        1.07453987082908195e-20, -1.71652329657411147e-21, -2.08674349559099603e-22,
        3.45371807219537732e-23, 4.14128225481140887e-24, -7.07211001730474564e-25,
        -8.37132209082386195e-26, 1.47009173000715507e-26, 1.71918537385513572e-27,
        -3.09600683514091908e-28, -3.57947041088649473e-29, 6.59499725647344126e-30,
        7.5430549748847196e-31, -1.41906088090354445e-31, -1.60657621844071818e-32,
        3.08089430710768301e-33, 3.45441121294919607e-34, -6.74270543587236185e-35,
        -7.49096137443440683e-36 } };

    complex<double> power { alpha * alpha };
    size_t i { 0 };
    complex<double> term { power * (C[i] * Lalpha + D[i]) };
    res += term;
    while (abs(term) > abs(res) * 1e-16) {
        i += 1;
        if (i >= C.size())
            throw FastGPL_error { "S22_series3 converges too slowly" };
        power *= alpha;
        term = power * (C[i] * Lalpha + D[i]);
        res += term;
    }
    return res;
}

// Relate S22(x) to S22(1/x)
// Should be used when |x| > 1 && |x-1| > 1
complex<double> FastGPL_internal::S22_inversion(complex<double> x, int s)
{
    complex<double> S22_inv { S22_series2(1. / x) };
    complex<double> l { Log(-x, -s) };
    complex<double> res { S22_inv - 2. * PolyLog(4, 1. / x) - l * PolyLog(3, 1. / x)
        - Zeta(4) * 7. / 4. + l * Zeta(3) + positive_int_power(l, 4) / factorial(4) };
    return res;
}
