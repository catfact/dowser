//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_OUTPUT_HPP
#define DOWSER_OUTPUT_HPP


#include "analysis.hpp"

namespace dowser {
    class output {

    public:
        enum class format_t:int {
            csv, python, supercollider, lua
        };

        struct config {
            juce::File file;
            format_t fmt;
        };

        static void perform(std::unique_ptr<analysis::results> results) {
            // TODO
        }
    };
}

#endif //DOWSER_OUTPUT_HPP
