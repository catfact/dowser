//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_OUTPUT_HPP
#define DOWSER_OUTPUT_HPP

#include <juce_core/juce_core.h>

#include "analysis.hpp"

namespace dowser {
    class output {

    private:
        /// helper
        static void log_field(juce::FileOutputStream &fos, const char *name, float val) {
            if (!std::isnan(val)) {
                fos <<  name << ":" << val << ", ";
            } 
        }
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
                fos << "  ( \n    ";
                log_field(fos, "papr", frame.papr);
                log_field(fos, "centroid", frame.centroid);
                log_field(fos, "flatness", frame.flatness);
                log_field(fos, "meanMag", frame.meanMag);
                log_field(fos, "maxMag", frame.maxMag);
                fos << "\n    ";
                fos << "peaks: [\n";
                for (auto peak: frame.magPeaks) {
                    fos << "      (hz:" << peak.first << ", mag:" << sqrt(peak.second) << "),\n";
                }

#if INCLUDE_AUTOCORR
                fos << "    acPeaks: [\n";
                for (auto peak: frame.acPeaks) {
                    fos << "      (hz:" << peak.first << ", ac:" << sqrt(peak.second) << "),\n";
                }
#endif
                fos << "    ]\n";
                fos << "  ),\n";
            }
            fos << "]\n\n";
        }
    };
}

#endif //DOWSER_OUTPUT_HPP
