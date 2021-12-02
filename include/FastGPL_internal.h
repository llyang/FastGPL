#include <complex>
#include <vector>

#include "FastGPL.h"

namespace FastGPL_internal {

inline std::complex<double> G(const std::vector<std::complex<double>>& a, double x)
{
    return FastGPL::G(a, x, true);
}

inline std::complex<double> G(
    const std::vector<std::complex<double>>& a, std::vector<int> s, double x)
{
    return FastGPL::G(a, s, x, true);
}

// Riemann zeta values for 2<=n<=8 and n=0
constexpr double Zeta(int n)
{
    if (n == 0)
        return -0.5;
    constexpr double z[7] { 1.6449340668482264365, 1.2020569031595942854,
        1.0823232337111381915, 1.0369277551433699263, 1.0173430619844491397,
        1.0083492773819228268, 1.0040773561979443394 };
    return z[n - 2];
}

// Harmonic numbers for n >= 1
constexpr double HarmonicNumber(int n)
{
    double res { 1 };
    for (int i = 2; i <= n; ++i) {
        res += 1. / i;
    }
    return res;
}

double factorial(unsigned int n);

template <class T> T positive_int_power(T x, int n);

std::complex<double> PolyLog_series2(int n, std::complex<double> x);
std::complex<double> PolyLog_series3(int n, std::complex<double> x, int s);

std::complex<double> PolyLog_inversion(int n, std::complex<double> x, int s);

std::complex<double> S22_series2(std::complex<double> x);
std::complex<double> S22_series3(std::complex<double> x, int s);

std::complex<double> S22_inversion(std::complex<double> x, int s);

bool is_convergent(const std::vector<std::complex<double>>& a, double x);
bool is_well_convergent(const std::vector<std::complex<double>>& a, double x);

std::complex<double> G_as_mLi(const std::vector<std::complex<double>>& a, double x);
template <class T>
std::complex<T> mLi_series(
    const std::vector<int>& m, const std::vector<std::complex<double>>& x);

std::complex<double> G_Hoelder(
    const std::vector<std::complex<double>>& a, const std::vector<int>& s, double x);

std::complex<double> G_remove_trailing_zeros(
    std::vector<std::complex<double>> a, std::vector<int> s, double x);

std::complex<double> G1(std::complex<double> a, int s, double x);

std::complex<double> G2_explicit(
    const std::vector<std::complex<double>>& a, const std::vector<int>& s, double x);

std::complex<double> G3_dispatch(
    const std::vector<std::complex<double>>& a, const std::vector<int>& s, double x);

std::complex<double> G4_dispatch(
    const std::vector<std::complex<double>>& a, const std::vector<int>& s, double x);
};
