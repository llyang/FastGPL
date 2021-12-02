#include <complex>
#include <vector>

#include "FastGPL.h"

using std::complex;

complex<double> FastGPL::EllipticK(complex<double> z)
{
    if (z.real() > 1. && z.imag() == 0)
        throw std::runtime_error("z > 1 in EllipticK");

    constexpr double eps { 1e-16 };

    complex<double> s { sqrt(1.0 - z) };
    complex<double> k { (1.0 - s) / (1.0 + s) }; // k1

    complex<double> r { 1.0 + k };

    do {
        s = sqrt((1.0 + k) * (1.0 - k));
        k = (1.0 - s) / (1.0 + s);
        r *= 1.0 + k;
    } while (abs(k) > eps);

    r *= M_PI / 2.0;

    return r;
}

complex<double> FastGPL::EllipticE(complex<double> z)
{
    if (z.real() > 1. && z.imag() == 0)
        throw std::runtime_error("z > 1 in EllipticE");

    constexpr double eps { 1e-16 };

    complex<double> s { sqrt(1.0 - z) };
    complex<double> b { s };
    complex<double> a { 1.0 + s };
    complex<double> k { (1.0 - s) / (1.0 + s) };

    do {
        s = sqrt((1.0 + k) * (1.0 - k)); // s[n+1]=sqrt(1-k[n]^2)
        b = b * (1.0 + k) + a * s; // b[n+1]=b[n]*(1+k[n])+a[n]*s[n+1]
        a *= 1.0 + s; // a[n+1]=a[n]*(1+s[n+1])
        k = (1.0 - s) / (1.0 + s); // k[n+1]=(1-s[n+1])/(1+s[n+1])
    } while (abs(k) > eps);

    return (a - (1.0 + k) * b) * M_PI / 2.0;
}
