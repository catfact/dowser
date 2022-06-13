//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_OUTPUT_HPP
#define DOWSER_OUTPUT_HPP

#include <juce_core/juce_core.h>

#include "analysis.hpp"

namespace dowser {

    enum class output_format_t : int {
        csv, python, supercollider, lua
    };

    template <output_format_t fmt>
    class output {

    private:
        juce::FileOutputStream *fos;

    private:
        void log_field(const char *name, float val) {
            if constexpr(fmt == output_format_t::supercollider) {
                if (!std::isnan(val)) {
                    *fos << name << ":" << val << ", ";
                } else {
                    *fos << name << ": nan" << ", ";
                }
            }
            else if constexpr(fmt == output_format_t::lua) {
                if (!std::isnan(val)) {
                    *fos << name << " = " << val << ", ";
                } else {
                    *fos << name << " = nan" << ", ";
                }
            }
            else if constexpr(fmt == output_format_t::python) {
                if (!std::isnan(val)) {
//                    *fos << name << " = " << val << ", ";
                } else {
//                    *fos << name << " = nan" << ", ";
                }
            }
        }

        void log_all(const analysis::results *results) {
            if constexpr(fmt == output_format_t::supercollider) {
                *fos << "[\n";
                for (auto frame: results->frames) {
                    *fos << "  ( \n    ";
                    log_field("papr", frame.papr);
                    log_field("centroid", frame.centroid);
                    log_field("flatness", frame.flatness);
                    log_field("meanMag", frame.meanMag);
                    log_field("maxMag", frame.maxMag);
                    log_field("fluxPositive", frame.fluxPositive);
                    log_field("fluxNegative", frame.fluxNegative);
                    *fos << "\n    ";
                    *fos << "peaks: [\n";
                    for (auto peak: frame.magPeaks) {
                        *fos << "      (hz:" << peak.hz << ", mag:" << sqrt(peak.pow) << ", persistence: "
                             << peak.persistence << "),\n";
                    }
                    *fos << "    ]\n";
                    *fos << "  ),\n";
                }
                *fos << "]\n\n";
            }

            else if constexpr(fmt == output_format_t::lua) {
                *fos << "{\n";
                for (auto frame: results->frames) {
                    *fos << "  { \n    ";
                    log_field("papr", frame.papr);
                    log_field("centroid", frame.centroid);
                    log_field("flatness", frame.flatness);
                    log_field("meanMag", frame.meanMag);
                    log_field("maxMag", frame.maxMag);
                    log_field("fluxPositive", frame.fluxPositive);
                    log_field("fluxNegative", frame.fluxNegative);
                    *fos << "\n    ";
                    *fos << "peaks: {\n";
                    for (auto peak: frame.magPeaks) {
                        *fos << "      (hz:" << peak.hz << ", mag:" << sqrt(peak.pow) << ", persistence: "
                             << peak.persistence << "),\n";
                    }
                    *fos << "    }\n";
                    *fos << "  },\n";
                }
                *fos << "}\n\n";
            }

            else if constexpr(fmt == output_format_t::python) {
//                *fos << "[\n";
                for (auto frame: results->frames) {
//...
                }
                //           *fos << "]\n\n";
            }


        }

    public:
        void perform(const analysis::results *results, const juce::File& outfile) {
            outfile.deleteFile();
            if (fos != nullptr)
                delete fos;

            fos = new juce::FileOutputStream(outfile);

            if (fos->failedToOpen()) {
                std::cerr << "failed to open file for output: " << outfile.getFullPathName() << std::endl;
                delete fos;
                return;
            }
            log_all(results);
        }
    };


}

#endif //DOWSER_OUTPUT_HPP
