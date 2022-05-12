#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

#include <chrono>
#include <random>

#include "FastGPL.h"

using std::complex;
using std::cout;
using std::vector;

void speed_test()
{
    constexpr double para_max { 1000 };
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_int_distribution<> dice(0, 6);
    std::uniform_real_distribution<> para(-para_max, para_max);

    auto start { std::chrono::steady_clock::now() };
    for (int i = 0; i < 10000; ++i) {
        vector<complex<double>> a;

        double x { fabs(para(gen)) };

        for (int j = 0; j < 4; ++j) {
            if (dice(gen) < 2)
                a.push_back(0);
            else
                a.push_back({ para(gen), para(gen) });
        }
        FastGPL::G(a, x);
    }
    auto end { std::chrono::steady_clock::now() };
    std::chrono::duration<double> diff { end - start };
    cout << "Time used for 10,000 weight-4 GPLs: " << diff.count() << " s" << std::endl;
}

void speed6_test()
{
    constexpr double para_max { 1000 };
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_int_distribution<> dice(0, 6);
    std::uniform_real_distribution<> para(-para_max, para_max);

    auto start { std::chrono::steady_clock::now() };
    for (int i = 0; i < 10000; ++i) {
        vector<complex<double>> a;

        double x { fabs(para(gen)) };

        for (int j = 0; j < 6; ++j) {
            if (dice(gen) < 2)
                a.push_back(0);
            else
                a.push_back({ para(gen), para(gen) });
        }
        FastGPL::G(a, x);
    }
    auto end { std::chrono::steady_clock::now() };
    std::chrono::duration<double> diff { end - start };
    cout << "Time used for 10,000 weight-6 GPLs: " << diff.count() << " s" << std::endl;
}

void simple_test()
{
    vector<complex<double>> a { { 0.1, 0.2 }, 0.4, 0.5, { -0.2, -0.3 } };
    double x { 0.8 };

    auto start { std::chrono::steady_clock::now() };

    cout << "G(0.1+0.2i, 0.4+0i, 0.5+0i, -0.2-0.3i; 0.8) = " << FastGPL::G(a, x) << '\n';

    vector<int> s { 1, -1, 1, -1 };
    cout << "G(0.1+0.2i, 0.4-0i, 0.5+0i, -0.2-0.3i; 0.8) = " << FastGPL::G(a, s, x) << '\n';

    auto end { std::chrono::steady_clock::now() };

    std::chrono::duration<double> diff { end - start };

    cout << "Time used: " << diff.count() << "s\n\n";
}

void weight6_test()
{
    vector<complex<double>> a { { 0.1, 0.2 }, 0.4, 0.5, { -0.2, -0.3 }, 0.6, 0.7 };
    double x { 0.8 };

    auto start { std::chrono::steady_clock::now() };

    cout << "G(0.1+0.2i, 0.4+0i, 0.5+0i, -0.2-0.3i, 0.6+0i, 0.7+0i; 0.8) = " << FastGPL::G(a, x) << '\n';

    vector<int> s { 1, -1, 1, -1, -1, -1 };
    cout << "G(0.1+0.2i, 0.4-0i, 0.5+0i, -0.2-0.3i, 0.6-0i, 0.7-0i; 0.8) = " << FastGPL::G(a, s, x) << '\n';

    auto end { std::chrono::steady_clock::now() };

    std::chrono::duration<double> diff { end - start };

    cout << "Time used: " << diff.count() << "s\n\n";
}

int main()
{

    cout << std::setprecision(16);

    simple_test();

    weight6_test();

    speed_test();

    speed6_test();
}
