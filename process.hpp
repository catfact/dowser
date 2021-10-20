//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_PROCESS_HPP
#define DOWSER_PROCESS_HPP

#include <complex>
#include <memory>

#include <juce_core/juce_core.h>


namespace dowser {
    namespace process {

        struct data {
            typedef std::vector<std::complex<float>> spectral_frame;
            std::vector<spectral_frame> frames;
        };

        std::unique_ptr<struct data> perform(juce::File soundfile) {
            auto data = std::make_unique<struct data>();
            // TODO
            return std::move(data);
        }

    };
}
#endif //DOWSER_PROCESS_HPP
