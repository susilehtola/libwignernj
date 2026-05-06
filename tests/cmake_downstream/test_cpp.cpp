// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Susi Lehtola
//
// Smoke test for the C++ wrapper of an installed libwignernj.

#include <cmath>
#include <iostream>
#include "wignernj.hpp"

int main()
{
    // 2*j integer form
    double v = wignernj::symbol3j<double>(2, 2, 0,  0, 0, 0);
    // Real-valued half-integer form
    double w = wignernj::symbol3j(1.0, 1.0, 0.0,  0.0, 0.0, 0.0);
    double ref = -1.0 / std::sqrt(3.0);
    if (std::fabs(v - ref) > 1e-14 || std::fabs(w - ref) > 1e-14) {
        std::cerr << "C++: symbol3j returned " << v << " / " << w
                  << ", expected " << ref << "\n";
        return 1;
    }
    std::cout << "C++: symbol3j(1,1,0;0,0,0) = " << v
              << " / " << w << "  (expected " << ref << ")\n";
    return 0;
}
