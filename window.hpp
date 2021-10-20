//
// Created by emb on 10/20/2021.
//

#ifndef DOWSER_WINDOW_HPP
#define DOWSER_WINDOW_HPP

#include <cmath>
#include <iostream>
#include <functional>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace dowser {

    namespace window {

        enum class shape_t: int {
            RECT = 0,
            BARTLETT,   // triangular
            WELCH,      // parabolic
            HANN,       // squared cosine
            HAMMING,       // raised squared cosine
            MLT,        // sine (modular lapped transform)
        };
        static constexpr int numShapes = 6;

        template<shape_t shape, int n>
        static std::array<float, n> data() {
            std::array<float, n> y;
            for(int i=0; i<n; ++i) {
                y[i] = calc<shape>(i, n);
            }
            return std::move(y);
        };

        template<shape_t shape>
        static float calc(int i, int n);
    };


    template<>
    float window::calc<window::shape_t::RECT>(int i, int n) {
        (void)i;
        (void)n;
        return 1.f;
    }

    template<>
    float window::calc<window::shape_t::BARTLETT>(int i, int n) {
        i = i - (n - 1) / 2;
        return static_cast<float>(1.f - (std::abs(i) / ((double)(n - 1) * 0.5f)));
    }

    template<>
    float window::calc<window::shape_t::WELCH>(int i, int n) {
        double y = ((i - (0.5 * (n - 1))) / (0.5 * (n - 1)));
        return static_cast<float>(1.0 - y * y);
    }

    template<>
    float window::calc<window::shape_t::HANN>(int i, int n) {
        double t = 2.0 * M_PI * (double) i / (double) (n - 1);
        double y = 0.5 - 0.5 * std::cos(t);
        return static_cast<float>(y);
    }

    template<>
    float window::calc<window::shape_t::HAMMING>(int i, int n) {
        double t = 2.0 * M_PI * (double) i / (double) (n - 1);
        double y = 0.54 - 0.46 * std::cos(t);
        // endpoints should be halved for COLA
        if (i == 0 || i == (n - 1)) { y *= 0.5; }
        return static_cast<float>(y);
    }

    template<>
    float window::calc<window::shape_t::MLT>(int i, int n) {
        double t = 2.0 * M_PI * (double) i / (double) (n - 1);
        double y = 0.5 - 0.5 * std::cos(t);
        return static_cast<float>(std::sqrt(y));
    }
}

#endif //DOWSER_WINDOW_HPP
