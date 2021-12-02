#include <complex>
#include <vector>

#include "FastGPL.h"

#include "FastGPL_internal.h"

using std::complex;

double FastGPL_internal::factorial(unsigned int n)
{
    double res { 1 };
    for (unsigned int i = 2; i <= n; ++i)
        res *= i;
    return res;
}

template <class T> T FastGPL_internal::positive_int_power(T x, int n)
{
    if (n < 0)
        throw FastGPL::FastGPL_error("negative n for positive_int_power");

    if (n == 0)
        return 1.0;

    T res { 1.0 };
    while (n > 0) {
        if (n % 2 != 0)
            res *= x;
        n /= 2;
        x *= x;
    }
    return res;
}

template double FastGPL_internal::positive_int_power<double>(double, int);
template complex<double> FastGPL_internal::positive_int_power<complex<double>>(
    complex<double>, int);
