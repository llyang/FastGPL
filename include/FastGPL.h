/*
 Library for fast evaluation of generalized polylogarithms

 Usage: G({a1,...,an}, x) or G({a1,...,an}, {s1,...,sn}, x)
 * Input: complex<double> ai
 * Input: int si (1 or -1, the sign of the imaginary part of ai; only useful when
 Im(ai)==0, otherwise will be overwritten according to Im(ai))
 * Input: double x
 * Output: complex<double>

 For best performance and precision, it is recommended to analytically transform the
 expressions to
 * Remove trailing zeros;
 * Normalize the argument to x=1;
 * Set Im(ai) to zero and set si appropriately if Im(ai) comes from i*eps prescription

 Usage: PolyLog(n, z) or PolyLog(n, z, s)
*/

#ifndef _FASTGPL_H
#define _FASTGPL_H

#include <complex>
#include <vector>

namespace FastGPL {

std::complex<double> Log(std::complex<double> x);
std::complex<double> Log(std::complex<double> x, int s);

std::complex<double> PolyLog(int n, std::complex<double> x);
std::complex<double> PolyLog(int n, std::complex<double> x, int s);

std::complex<double> S22(std::complex<double> x);
std::complex<double> S22(std::complex<double> x, int s);

constexpr double Hoelder_threshold { 1.1 };

std::complex<double> G(
    const std::vector<std::complex<double>>& a, double x, bool is_clean = false);
std::complex<double> G(const std::vector<std::complex<double>>& a, std::vector<int> s,
    double x, bool is_clean = false);

std::complex<double> EllipticK(std::complex<double> z);
std::complex<double> EllipticE(std::complex<double> z);

class FastGPL_error : public std::runtime_error {
public:
    explicit FastGPL_error(const std::string& s)
        : std::runtime_error(s)
    {
    }
};

};

#endif //_FASTGPL_H
