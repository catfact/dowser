//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_ANALYSIS_HPP
#define DOWSER_ANALYSIS_HPP

#include <memory>

#include "process.hpp"

namespace dowser {

    template<int fftOrder, int overlap>
    class analysis {
    public:
        struct results {

            struct frame {
                std::vector<std::pair<float,float>> peaks;
                float papr;
                float centroid;
                float flatness;
                float totalPower;
            };

            std::vector<struct frame> frames;
        };

        static std::unique_ptr<struct results> perform(std::unique_ptr<typename process<fftOrder, overlap>::data> data) {
            auto res = std::make_unique<struct results>();
            // TODO
            return std::move(res);
        }
    };
}

#endif //DOWSER_ANALYSIS_HPP