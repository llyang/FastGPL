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
//#include<iostream>

complex<double> G6_explicit_0000aa(complex<double> a, int sa, double x);
complex<double> G6_explicit_a0000a(complex<double> a, int sa, double x);
complex<double> G6_explicit_000a0a(complex<double> a, int sa, double x);
complex<double> G6_explicit_0a000a(complex<double> a, int sa, double x);
complex<double> G6_explicit_00a00a(complex<double> a, int sa, double x);

complex<double> G6_explicit_000aaa(complex<double> a, int sa, double x);
complex<double> G6_explicit_a000aa(complex<double> a, int sa, double x);
complex<double> G6_explicit_aa000a(complex<double> a, int sa, double x);
complex<double> G6_explicit_0a00aa(complex<double> a, int sa, double x);
complex<double> G6_explicit_00a0aa(complex<double> a, int sa, double x);
complex<double> G6_explicit_0aa00a(complex<double> a, int sa, double x);
complex<double> G6_explicit_00aa0a(complex<double> a, int sa, double x);
complex<double> G6_explicit_a0a00a(complex<double> a, int sa, double x);
complex<double> G6_explicit_a00a0a(complex<double> a, int sa, double x);
complex<double> G6_explicit_0a0a0a(complex<double> a, int sa, double x);

complex<double> G6_explicit_00aaaa(complex<double> a, int sa, double x);
complex<double> G6_explicit_a00aaa(complex<double> a, int sa, double x);
complex<double> G6_explicit_aa00aa(complex<double> a, int sa, double x);
complex<double> G6_explicit_aaa00a(complex<double> a, int sa, double x);
complex<double> G6_explicit_0a0aaa(complex<double> a, int sa, double x);
complex<double> G6_explicit_0aa0aa(complex<double> a, int sa, double x);
complex<double> G6_explicit_0aaa0a(complex<double> a, int sa, double x);
complex<double> G6_explicit_a0a0aa(complex<double> a, int sa, double x);
complex<double> G6_explicit_a0aa0a(complex<double> a, int sa, double x);
complex<double> G6_explicit_aa0a0a(complex<double> a, int sa, double x);

complex<double> G6_explicit_0aaaaa(complex<double> a, int sa, double x);
complex<double> G6_explicit_a0aaaa(complex<double> a, int sa, double x);
complex<double> G6_explicit_aa0aaa(complex<double> a, int sa, double x);
complex<double> G6_explicit_aaa0aa(complex<double> a, int sa, double x);
complex<double> G6_explicit_aaaa0a(complex<double> a, int sa, double x);

complex<double> G6_0000ab(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G6_a0000b(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G6_000a0b(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G6_0a000b(complex<double> a, complex<double> b, int sa, int sb, double x);
complex<double> G6_00a00b(complex<double> a, complex<double> b, int sa, int sb, double x);

complex<double> G6_000abc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G6_a000bc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G6_ab000c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G6_00a0bc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G6_0a00bc(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G6_00ab0c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G6_0ab00c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G6_a00b0c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G6_a0b00c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);
complex<double> G6_0a0b0c(complex<double> a, complex<double> b, complex<double> c, int sa,
    int sb, int sc, double x);

complex<double> G6_00abcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_00abcd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_00abcd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_00abcd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_00abcd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_a00bcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a00bcd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a00bcd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a00bcd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a00bcd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_ab00cd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_ab00cd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_ab00cd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_ab00cd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_ab00cd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_abc00d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_abc00d_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_abc00d_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_abc00d_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_abc00d_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_0a0bcd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0a0bcd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0a0bcd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0a0bcd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0a0bcd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_0ab0cd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0ab0cd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0ab0cd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0ab0cd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0ab0cd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_0abc0d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0abc0d_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0abc0d_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0abc0d_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_0abc0d_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_a0b0cd(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a0b0cd_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a0b0cd_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a0b0cd_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a0b0cd_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_a0bc0d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a0bc0d_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a0bc0d_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a0bc0d_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_a0bc0d_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_ab0c0d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_ab0c0d_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_ab0c0d_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_ab0c0d_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);
complex<double> G6_ab0c0d_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, int sa, int sb, int sc, int sd, double x);

complex<double> G6_0abcde(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_0abcde_e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_0abcde_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_0abcde_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_0abcde_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_0abcde_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);

complex<double> G6_a0bcde(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_a0bcde_e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_a0bcde_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_a0bcde_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_a0bcde_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_a0bcde_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);

complex<double> G6_ab0cde(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_ab0cde_e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_ab0cde_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_ab0cde_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_ab0cde_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_ab0cde_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);

complex<double> G6_abc0de(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abc0de_e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abc0de_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abc0de_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abc0de_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abc0de_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);

complex<double> G6_abcd0e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abcd0e_e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abcd0e_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abcd0e_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abcd0e_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);
complex<double> G6_abcd0e_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, int sa, int sb, int sc, int sd, int se,
    double x);

complex<double> G6_abcdef(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, complex<double> f, int sa, int sb, int sc,
    int sd, int se, int sf, double x);
complex<double> G6_abcdef_f(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, complex<double> f, int sa, int sb, int sc,
    int sd, int se, int sf, double x);
complex<double> G6_abcdef_e(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, complex<double> f, int sa, int sb, int sc,
    int sd, int se, int sf, double x);
complex<double> G6_abcdef_d(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, complex<double> f, int sa, int sb, int sc,
    int sd, int se, int sf, double x);
complex<double> G6_abcdef_c(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, complex<double> f, int sa, int sb, int sc,
    int sd, int se, int sf, double x);
complex<double> G6_abcdef_b(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, complex<double> f, int sa, int sb, int sc,
    int sd, int se, int sf, double x);
complex<double> G6_abcdef_a(complex<double> a, complex<double> b, complex<double> c,
    complex<double> d, complex<double> e, complex<double> f, int sa, int sb, int sc,
    int sd, int se, int sf, double x);

complex<double> FastGPL_internal::G6_dispatch(
    const vector<complex<double>>& a, const vector<int>& s, double x)
{
    if (a[0] == 0.0) {
        if (a[1] == 0.0) {
            if (a[2] == 0.0) {
                if (a[3] == 0.0) { // 0000ab
                    if (a[4] == a[5])
                        if (s[4] == s[5])
                            return G6_explicit_0000aa(a[5], s[5], x);
                        else
                            throw FastGPL_error { "GPL: a[4]==a[5] but s[4]!=s[5]" };
                    else
                        return G6_0000ab(a[4], a[5], s[4], s[5], x);
                } else if (a[4] == 0.0) { // 000a0b
                    if (a[3] == a[5])
                        if (s[3] == s[5])
                            return G6_explicit_000a0a(a[5], s[5], x);
                        else
                            throw FastGPL_error { "GPL: a[3]==a[5] but s[3]!=s[5]" };
                    else
                        return G6_000a0b(a[3], a[5], s[3], s[5], x);
                } else { // 000abc//345
                    if (a[3] == a[4] && a[3] == a[5])
                        if (s[3] == s[4] && s[3] == s[5])
                            return G6_explicit_000aaa(a[5], s[5], x);
                        else
                            throw FastGPL_error {
                                "GPL: a[3]==a[4]==a[5] but s[3]!=s[4]!=s[5]"
                            };
                    else
                        return G6_000abc(a[3], a[4], a[5], s[3], s[4], s[5], x);
                }
            } else if (a[3] == 0.0) {
                if (a[4] == 0.0) { // 00a00b
                    if (a[2] == a[5])
                        if (s[2] == s[5])
                            return G6_explicit_00a00a(a[5], s[5], x);
                        else
                            throw FastGPL_error { "GPL: a[2]==a[5] but s[2]!=s[5]" };
                    else
                        return G6_00a00b(a[2], a[5], s[2], s[5], x);
                } else { // 00a0bc//245
                    if (a[2] == a[4] && a[2] == a[5])
                        if (s[2] == s[4] && s[2] == s[5])
                            return G6_explicit_00a0aa(a[5], s[5], x);
                        else
                            throw FastGPL_error {
                                "GPL: a[2]==a[4]==a[5] but s[2]!=s[4]!=s[5]"
                            };
                    else
                        return G6_00a0bc(a[2], a[4], a[5], s[2], s[4], s[5], x);
                }
            } else if (a[4] == 0.0) { // 00ab0c//235
                if (a[2] == a[3] && a[2] == a[5])
                    if (s[2] == s[3] && s[2] == s[5])
                        return G6_explicit_00aa0a(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[2]==a[3]==a[5] but s[2]!=s[3]!=s[5]"
                        };
                else
                    return G6_00ab0c(a[2], a[3], a[5], s[2], s[3], s[5], x);

            } else { // 00abcd//2345
                if (a[2] == a[3] && a[2] == a[4] && a[2] == a[5])
                    if (s[2] == s[3] && s[2] == s[4] && s[2] == s[5])
                        return G6_explicit_00aaaa(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[2]=a[3]==a[4]==a[5] but s[2]!=s[3]!=s[4]!=s[5]"
                        };
                else
                    return G6_00abcd(a[2], a[3], a[4], a[5], s[2], s[3], s[4], s[5], x);
            }
        } else if (a[2] == 0.0) {
            if (a[3] == 0.0) { // 0a00()()
                if (a[4] == 0.0) { // 0a000b//15
                    if (a[1] == a[5])
                        if (s[1] == s[5])
                            return G6_explicit_0a000a(a[5], s[5], x);
                        else
                            throw FastGPL_error { "GPL: a[1]==a[5] but s[1]!=s[5]" };
                    else
                        return G6_0a000b(a[1], a[5], s[1], s[5], x);
                } else { // 0a00bc//145
                    if (a[1] == a[4] && a[1] == a[5])
                        if (s[1] == s[4] && s[1] == s[5])
                            return G6_explicit_0a00aa(a[5], s[5], x);
                        else
                            throw FastGPL_error {
                                "GPL: a[1]==a[4]==a[5] but s[1]!=s[4]!=s[5]"
                            };
                    else
                        return G6_0a00bc(a[1], a[4], a[5], s[1], s[4], s[5], x);
                }
            } else if (a[4] == 0.0) { // 0a0b0c//135
                if (a[1] == a[3] && a[1] == a[5])
                    if (s[1] == s[3] && s[1] == s[5])
                        return G6_explicit_0a0a0a(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[1]==a[3]==a[5] but s[1]!=s[3]!=s[5]"
                        };
                else
                    return G6_0a0b0c(a[1], a[3], a[5], s[1], s[3], s[5], x);
            } else { // 0a0bcd//1345
                if (a[1] == a[3] && a[1] == a[4] && a[1] == a[5])
                    if (s[1] == s[3] && s[1] == s[4] && s[1] == s[5])
                        return G6_explicit_0a0aaa(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[1]==a[3]==a[4]==a[5] but s[1]!=s[3]!=s[4]!=s[5]"
                        };
                else
                    return G6_0a0bcd(a[1], a[3], a[4], a[5], s[1], s[3], s[4], s[5], x);
            }
        } else if (a[3] == 0.0) { // 0ab0()()
            if (a[4] == 0.0) { // 0ab00c//125
                if (a[1] == a[2] && a[1] == a[5])
                    if (s[1] == s[2] && s[1] == s[5])
                        return G6_explicit_0aa00a(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[1]==a[2]==a[5] but s[1]!=s[2]!=s[5]"
                        };
                else
                    return G6_0ab00c(a[1], a[2], a[5], s[1], s[2], s[5], x);
            } else { // 0ab0cd//1245
                if (a[1] == a[2] && a[1] == a[4] && a[1] == a[5])
                    if (s[1] == s[2] && s[1] == s[4] && s[1] == s[5])
                        return G6_explicit_0aa0aa(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[1]==a[2]==a[4]==a[5] but s[1]!=s[2]!=s[4]!=s[5]"
                        };
                else
                    return G6_0ab0cd(a[1], a[2], a[4], a[5], s[1], s[2], s[4], s[5], x);
            }

        } else if (a[4] == 0.0) { // 0abc0d//1235
            if (a[1] == a[2] && a[1] == a[3] && a[1] == a[5])
                if (s[1] == s[2] && s[1] == s[3] && s[1] == s[5])
                    return G6_explicit_0aaa0a(a[5], s[5], x);
                else
                    throw FastGPL_error {
                        "GPL: a[1]==a[2]==a[3]==a[5] but s[1]!=s[2]!=s[3]!=s[5]"
                    };
            else
                return G6_0abc0d(a[1], a[2], a[3], a[5], s[1], s[2], s[3], s[5], x);
        } else { // 0abcde//12345
            if (a[1] == a[2] && a[1] == a[3] && a[1] == a[4] && a[1] == a[5])
                if (s[1] == s[2] && s[1] == s[3] && s[1] == s[4] && s[1] == s[5])
                    return G6_explicit_0aaaaa(a[5], s[5], x);
                else
                    throw FastGPL_error { "GPL: a[1]==a[2]==a[3]==a[4]==a[5] but "
                                          "s[1]!=s[2]!=s[3]!=s[4]!=s[5]" };
            else
                return G6_0abcde(
                    a[1], a[2], a[3], a[4], a[5], s[1], s[2], s[3], s[4], s[5], x);
        }
    } else if (a[1] == 0.0) { // a0()()()()
        if (a[2] == 0.0) { // a00()()()
            if (a[3] == 0.0) { // a000()()
                if (a[4] == 0.0) { // a0000b//05
                    if (a[0] == a[5])
                        if (s[0] == s[5])
                            return G6_explicit_a0000a(a[5], s[5], x);
                        else
                            throw FastGPL_error { "GPL: a[0]==a[5] but s[0]!=s[5]" };
                    else
                        return G6_a0000b(a[0], a[5], s[0], s[5], x);
                } else { // a000bc//045
                    if (a[0] == a[4] && a[0] == a[5])
                        if (s[0] == s[4] && s[0] == s[5])
                            return G6_explicit_a000aa(a[5], s[5], x);
                        else
                            throw FastGPL_error {
                                "GPL: a[0]==a[4]==a[5] but s[0]!=s[4]!=s[5]"
                            };
                    else
                        return G6_a000bc(a[0], a[4], a[5], s[0], s[4], s[5], x);
                }
            } else if (a[4] == 0.0) { // a00b0c//035
                if (a[0] == a[3] && a[0] == a[5])
                    if (s[0] == s[3] && s[0] == s[5])
                        return G6_explicit_a00a0a(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[0]==a[3]==a[5] but s[0]!=s[3]!=s[5]"
                        };
                else
                    return G6_a00b0c(a[0], a[3], a[5], s[0], s[3], s[5], x);
            } else { // a00bcd//0345
                if (a[0] == a[3] && a[0] == a[4] && a[0] == a[5])
                    if (s[0] == s[3] && s[0] == s[4] && s[0] == s[5])
                        return G6_explicit_a00aaa(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[0]==a[3]==a[4]==a[5] but s[0]!=s[3]!=s[4]!=s[5]"
                        };
                else
                    return G6_a00bcd(a[0], a[3], a[4], a[5], s[0], s[3], s[4], s[5], x);
            }

        } else if (a[3] == 0.0) { // a0b0()()
            if (a[4] == 0.0) { // a0b00c//025
                if (a[0] == a[2] && a[0] == a[5])
                    if (s[0] == s[2] && s[0] == s[5])
                        return G6_explicit_a0a00a(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[0]==a[2]==a[5] but s[0]!=s[2]!=s[5]"
                        };
                else
                    return G6_a0b00c(a[0], a[2], a[5], s[0], s[2], s[5], x);
            } else { // a0b0cd//0245
                if (a[0] == a[2] && a[0] == a[4] && a[0] == a[5])
                    if (s[0] == s[2] && s[0] == s[4] && s[0] == s[5])
                        return G6_explicit_a0a0aa(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[0]==a[2]==a[4]==a[5] but s[0]!=s[2]!=s[4]!=s[5]"
                        };
                else
                    return G6_a0b0cd(a[0], a[2], a[4], a[5], s[0], s[2], s[4], s[5], x);
            }
        } else if (a[4] == 0.0) { // a0bc0d//0235
            if (a[0] == a[2] && a[0] == a[3] && a[0] == a[5])
                if (s[0] == s[2] && s[0] == s[3] && s[0] == s[5])
                    return G6_explicit_a0aa0a(a[5], s[5], x);
                else
                    throw FastGPL_error {
                        "GPL: a[0]==a[2]==a[3]==a[5] but s[0]!=s[2]!=s[3]!=s[5]"
                    };
            else
                return G6_a0bc0d(a[0], a[2], a[3], a[5], s[0], s[2], s[3], s[5], x);
        } else { // a0bcde//02345
            if (a[0] == a[2] && a[0] == a[3] && a[0] == a[4] && a[0] == a[5])
                if (s[0] == s[2] && s[0] == s[3] && s[0] == s[4] && s[0] == s[5])
                    return G6_explicit_a0aaaa(a[5], s[5], x);
                else
                    throw FastGPL_error { "GPL: a[0]==a[2]==a[3]==a[4]==a[5] but "
                                          "s[0]!=s[2]!=s[3]!=s[4]!=s[5]" };
            else
                return G6_a0bcde(
                    a[0], a[2], a[3], a[4], a[5], s[0], s[2], s[3], s[4], s[5], x);
        }

    } else if (a[2] == 0.0) { // ab0()()()
        if (a[3] == 0.0) { // ab00()()
            if (a[4] == 0.0) { // ab000c//015
                if (a[0] == a[1] && a[0] == a[5])
                    if (s[0] == s[1] && s[0] == s[5])
                        return G6_explicit_aa000a(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[0]==a[1]==a[5] but s[0]!=s[1]!=s[5]"
                        };
                else
                    return G6_ab000c(a[0], a[1], a[5], s[0], s[1], s[5], x);
            } else { // ab00cd//0145
                if (a[0] == a[1] && a[0] == a[4] && a[0] == a[5])
                    if (s[0] == s[1] && s[0] == s[4] && s[0] == s[5])
                        return G6_explicit_aa00aa(a[5], s[5], x);
                    else
                        throw FastGPL_error {
                            "GPL: a[0]==a[1]==a[4]==a[5] but s[0]!=s[1]!=s[4]!=s[5]"
                        };
                else
                    return G6_ab00cd(a[0], a[1], a[4], a[5], s[0], s[1], s[4], s[5], x);
            }

        } else if (a[4] == 0.0) { // ab0c0d //0135
            if (a[0] == a[1] && a[0] == a[3] && a[0] == a[5])
                if (s[0] == s[1] && s[0] == s[3] && s[0] == s[5])
                    return G6_explicit_aa0a0a(a[5], s[5], x);
                else
                    throw FastGPL_error {
                        "GPL: a[0]==a[1]==a[3]==a[5] but s[0]!=s[1]!=s[3]!=s[5]"
                    };
            else
                return G6_ab0c0d(a[0], a[1], a[3], a[5], s[0], s[1], s[3], s[5], x);
        } else { // ab0cde//01345
            if (a[0] == a[1] && a[0] == a[3] && a[0] == a[4] && a[0] == a[5])
                if (s[0] == s[1] && s[0] == s[3] && s[0] == s[4] && s[0] == s[5])
                    return G6_explicit_aa0aaa(a[5], s[5], x);
                else
                    throw FastGPL_error { "GPL: a[0]==a[1]==a[3]==a[4]==a[5] but "
                                          "s[0]!=s[1]!=s[3]!=s[4]!=s[5]" };
            else
                return G6_ab0cde(
                    a[0], a[1], a[3], a[4], a[5], s[0], s[1], s[3], s[4], s[5], x);
        }

    } else if (a[3] == 0.0) { // abc0()()
        if (a[4] == 0.0) { // abc00d//0125
            if (a[0] == a[1] && a[0] == a[2] && a[0] == a[5])
                if (s[0] == s[1] && s[0] == s[2] && s[0] == s[5])
                    return G6_explicit_aaa00a(a[5], s[5], x);
                else
                    throw FastGPL_error {
                        "GPL: a[0]==a[1]==a[2]==a[5] but s[0]!=s[1]!=s[2]!=s[5]"
                    };
            else
                return G6_abc00d(a[0], a[1], a[2], a[5], s[0], s[1], s[2], s[5], x);
        } else { // abc0de//01245
            if (a[0] == a[1] && a[0] == a[2] && a[0] == a[4] && a[0] == a[5])
                if (s[0] == s[1] && s[0] == s[2] && s[0] == s[4] && s[0] == s[5])
                    return G6_explicit_aaa0aa(a[5], s[5], x);
                else
                    throw FastGPL_error { "GPL: a[0]==a[1]==a[2]==a[4]==a[5] but "
                                          "s[0]!=s[1]!=s[2]!=s[4]!=s[5]" };
            else
                return G6_abc0de(
                    a[0], a[1], a[2], a[4], a[5], s[0], s[1], s[2], s[4], s[5], x);
        }

    } else if (a[4] == 0.0) { // abcd0e//01235
        if (a[0] == a[1] && a[0] == a[2] && a[0] == a[3] && a[0] == a[5])
            if (s[0] == s[1] && s[0] == s[2] && s[0] == s[3] && s[0] == s[5])
                return G6_explicit_aaaa0a(a[5], s[5], x);
            else
                throw FastGPL_error {
                    "GPL: a[0]==a[1]==a[2]==a[3]==a[5] but s[0]!=s[1]!=s[2]!=s[3]!=s[5]"
                };
        else
            return G6_abcd0e(
                a[0], a[1], a[2], a[3], a[5], s[0], s[1], s[2], s[3], s[5], x);
    } else { // abcdef
        return G6_abcdef(
            a[0], a[1], a[2], a[3], a[4], a[5], s[0], s[1], s[2], s[3], s[4], s[5], x);
    }
}
