//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_OUTPUT_HPP
#define DOWSER_OUTPUT_HPP

#include <juce_core/juce_core.h>

#include "analysis.hpp"

namespace dowser {
    class output {

    public:
        enum class format_t : int {
            csv, python, supercollider, lua
        };

        struct config {
            juce::File file;
            format_t fmt;
        };

        static void perform(std::unique_ptr<analysis::results> results, const juce::File& outfile) {
            outfile.deleteFile();
            juce::FileOutputStream fos(outfile);

            if (fos.failedToOpen()) {
                std::cerr << "failed to open file for output: " << outfile.getFullPathName() << std::endl;
                return;
            }

            fos << "[\n";
            for (auto frame: results->frames) {
                fos << "  ( \n";
                fos << "    papr:" << frame.papr << ", centroid:" << frame.centroid << ", flatness:" << frame.flatness
                << ", meanMag:" << frame.meanMag << ", maxMag:" << frame.maxMag << ",\n";
                fos << "    peaks: [\n";
                for (auto peak: frame.magPeaks) {
                    fos << "      (hz:" << peak.first << ", mag:" << sqrt(peak.second) << "),\n";
                }
                fos << "    ]\n";
                fos << "  ),\n";
            }
            fos << "]\n\n";
        }
    };
}

#endif //DOWSER_OUTPUT_HPP
