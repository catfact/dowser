#pragma once

#include <cmath>

#include <juce_core/juce_core.h>

#include "analysis.hpp"

#include "output_format.hpp"

//class results;

namespace dowser
{
    class output
    {

    private:
        template<output_format_t fmt>
        class logger
        {
        public:
            static void log_field(juce::FileOutputStream *fos, const char *name, float val)
            {
                if constexpr (fmt == output_format_t::supercollider)
                {
                    if (!std::isnan(val))
                    {
                        *fos << name << ":" << val << ", ";
                    } else
                    {
                        *fos << name << ": nan"
                             << ", ";
                    }
                } else if constexpr (fmt == output_format_t::lua)
                {
                    if (!std::isnan(val))
                    {
                        *fos << name << " = " << val << ", ";
                    } else
                    {
                        *fos << name << " = math.nan"
                             << ", ";
                    }
                } else if constexpr (fmt == output_format_t::python)
                {
                    if (!std::isnan(val))
                    {
                        *fos << "\"" << name << "\": " << val << ", ";
                    } else
                    {
                        *fos << name << ": NaN"
                             << ", ";
                    }
                }
            }

            // FIXME: could be DRYd
            static void log_all(juce::FileOutputStream *fos, const analysis::results *results)
            {
                if constexpr (fmt == output_format_t::supercollider)
                {
                    *fos << "[\n";
                    for (auto frame: results->frames)
                    {
                        *fos << "  ( \n    ";
                        log_field(fos, "papr", frame.papr);
                        log_field(fos, "centroid", frame.centroid);
                        log_field(fos, "flatness", frame.flatness);
                        log_field(fos, "meanMag", frame.meanMag);
                        log_field(fos, "maxMag", frame.maxMag);
                        log_field(fos, "fluxPositive", frame.fluxPositive);
                        log_field(fos, "fluxNegative", frame.fluxNegative);
                        *fos << "\n    ";
                        *fos << "peaks: [\n";
                        for (auto peak: frame.magPeaks)
                        {
                            *fos << "( ";
                            log_field(fos, "hz", peak.hz);
                            log_field(fos, "mag", sqrt(peak.pow));
                            log_field(fos, "persist", peak.persistence);
                            *fos << " ),\n";
                        }
                        *fos << "    ]\n";
                        *fos << "  ),\n";
                    }
                    *fos << "]\n\n";
                } else if constexpr (fmt == output_format_t::lua)
                {
                    *fos << "{\n";
                    for (auto frame: results->frames)
                    {
                        *fos << "  { \n    ";
                        log_field(fos, "papr", frame.papr);
                        log_field(fos, "centroid", frame.centroid);
                        log_field(fos, "flatness", frame.flatness);
                        log_field(fos, "meanMag", frame.meanMag);
                        log_field(fos, "maxMag", frame.maxMag);
                        log_field(fos, "fluxPositive", frame.fluxPositive);
                        log_field(fos, "fluxNegative", frame.fluxNegative);
                        *fos << "\n    ";
                        *fos << "peaks = {\n";
                        for (auto peak: frame.magPeaks)
                        {
                            *fos << "{ ";
                            log_field(fos, "hz", peak.hz);
                            log_field(fos, "mag", sqrt(peak.pow));
                            log_field(fos, "persist", peak.persistence);
                            *fos << " }\n";
                        }
                        *fos << "    }\n";
                        *fos << "  },\n";
                    }
                    *fos << "}\n\n";
                } else if constexpr (fmt == output_format_t::python)
                {
                    *fos << "[\n";
                    for (auto frame: results->frames)
                    {
                        *fos << "  { \n    ";
                        log_field(fos, "papr", frame.papr);
                        log_field(fos, "centroid", frame.centroid);
                        log_field(fos, "flatness", frame.flatness);
                        log_field(fos, "meanMag", frame.meanMag);
                        log_field(fos, "maxMag", frame.maxMag);
                        log_field(fos, "fluxPositive", frame.fluxPositive);
                        log_field(fos, "fluxNegative", frame.fluxNegative);
                        *fos << "\n    ";
                        *fos << "peaks: [\n";
                        for (auto peak: frame.magPeaks)
                        {
                            *fos << "{ ";
                            log_field(fos, "hz", peak.hz);
                            log_field(fos, "mag", sqrt(peak.pow));
                            log_field(fos, "persist", peak.persistence);
                            *fos << " },\n";
                        }
                        *fos << "    ]\n";
                        *fos << "  },\n";
                    }
                    *fos << "]\n\n";
                }
            }
        };

    private:
        juce::FileOutputStream *fos{nullptr};

    public:
        void perform(output_format_t fmt, const analysis::results *results, const juce::File &outfile)
        {
            outfile.deleteFile();
            if (fos != nullptr)
                delete fos;

            fos = new juce::FileOutputStream(outfile);

            if (fos->failedToOpen())
            {
                std::cerr << "failed to open file for output: " << outfile.getFullPathName() << std::endl;
                delete fos;
                return;
            }
            switch (fmt)
            {
                case output_format_t::supercollider:
                    logger<output_format_t::supercollider>::log_all(fos, results);
                    break;
                case output_format_t::python:
                    logger<output_format_t::python>::log_all(fos, results);
                    break;
                case output_format_t::lua:
                    logger<output_format_t::lua>::log_all(fos, results);
                    break;
            }
            delete fos;
        }
    };
}
