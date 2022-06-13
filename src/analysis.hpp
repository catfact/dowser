//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_ANALYSIS_HPP
#define DOWSER_ANALYSIS_HPP

#include <memory>

#include "process.hpp"
//#include "findpeaks/findpeaks.hpp"
//#include "findpeaks/persistence.hpp"
#include "peaks_watershed.hpp"

#include "process_config.hpp"

namespace dowser {

    class analysis {

    public:
        //typedef std::pair<float, float> peak_t;
        struct MagPeak {
            float hz;
            float pow;
            float persistence;
        };
        struct results {

            struct frame {
                std::vector<MagPeak> magPeaks;

#if INCLUDE_AUTOCORR
                std::vector<peak_t> acPeaks;
#endif
                float papr{};
                float centroid{};
                float flatness{};
                float meanMag{};
                float maxMag{};
                float fluxPositive;
                float fluxNegative;
            };

            std::vector<struct frame> frames;
        };

        template<int fftOrder>
        static std::unique_ptr<struct results>
        perform(std::unique_ptr<typename process<fftOrder>::analysisData> analysisData,
                const dowser::ProcessConfig &config) {
            static constexpr int fftSize = process<fftOrder>::fftSizePadded;
            static constexpr int numRealBins = fftSize / 2 + 1;

            auto res = std::make_unique<struct results>();

            auto sr = analysisData->sampleRate;
            int minBin = std::max(1, static_cast<int>(config.minHz * fftSize / sr));
            int maxBin = std::min(numRealBins - 2, static_cast<int>(config.maxHz * fftSize / sr));

            for (size_t i = 0; i < analysisData->specPowFrames.size(); ++i) {
                const auto &powFrame = analysisData->specPowFrames[i];
#if INCLUDE_AUTOCORR
                const auto &acFrame = analysisData->autoCorrFrames[i];
#endif
                typename results::frame resFrame;

                resFrame.magPeaks = getPowPeaks<fftSize>(powFrame.data(), sr,
                                                         config.maxPeaksPerFrame,
                                                         powf(10, config.minPowDb / 20),
                                                         config.minPersistence,
                                                         static_cast<unsigned int>(minBin),
                                                         static_cast<unsigned int>(maxBin));

#if INCLUDE_AUTOCORR
                resFrame.acPeaks = getMagPeaksFromPowFrames<fftSize>(acFrame.data(), sr, maxPeaks, 0,
                                                                     static_cast<unsigned int>(minBin),
                                                                     static_cast<unsigned int>(maxBin));
#endif
                computeStats<fftSize>(resFrame, powFrame.data(), sr, minBin, maxBin);

                if (i > 0) {
                    const auto &powFrame0 = analysisData->specPowFrames[i - 1];
                    double fluxUp = 0;
                    double fluxDown = 0;
                    for (size_t bin = 0; bin < powFrame.size(); ++bin) {
                        double diff = sqrt(powFrame[bin]) - sqrt(powFrame0[bin]);
                        if (diff > 0) {
                            fluxUp += diff;
                        } else {
                            fluxDown += diff;
                        }
                    }

                    resFrame.fluxPositive = fluxUp;
                    resFrame.fluxNegative = fluxDown;
                } else {
                    resFrame.fluxPositive = 0.f;
                    resFrame.fluxNegative = 0.f;
                }

                res->frames.push_back(resFrame);
            }

            return res;
        }

    private:
        template<int fftSize>
        static std::vector<MagPeak> getPowPeaks(const double *powBuf, double sr,
                                                unsigned int maxPeaks,
                                                double minPow, double minPersist,
                                                unsigned int minBin, unsigned int maxBin) {
            std::vector<MagPeak> y; // results vector

            std::vector<double> workBuffer;
            unsigned int bin = minBin;
            while (bin <= maxBin) {
                workBuffer.push_back(log10(powBuf[bin]) * 20);
                bin++;
            }
            auto watershedPeaks = dowser::peaks::Watershed<double>::findPeaks(
                    workBuffer.data(), static_cast<unsigned int>(maxBin - minBin + 1));

            // normalize persistence
            std::vector<unsigned int> idx;
            for (auto &p: watershedPeaks) {
                unsigned int i = static_cast<unsigned int>(p.index) + minBin;
                if (powBuf[i] > minPow && p.persistence > minPersist) {
                    idx.push_back(i);
                }
            }

            // interpolate true locations / magnitudes
            y.reserve(idx.size());
            for (auto &pos: idx) {
                auto peak = refinePeak<fftSize>(powBuf, static_cast<int>(pos), sr);
                peak.persistence = watershedPeaks[pos - minBin].persistence;
                if (peak.pow > 0) {
                    y.push_back(peak);
                }
            }

            std::sort(y.begin(), y.end(), [](MagPeak a, MagPeak b) {
                return a.persistence > b.persistence;
            });

            // normalize persistence
            float maxPersistence = y[1].persistence;
            for (auto &p: y) {
                if (p.persistence <= maxPersistence) {
                    p.persistence /= maxPersistence;
                }
            }

            if (y.size() > maxPeaks) {
                y = {y.begin(), y.begin() + maxPeaks};
            }
            return y;
        }

        // approximate true peak location by quadratic fit
        // assumption: pos-1, pos+1 are in range
        template<int fftSize>
        static MagPeak refinePeak(const double *powBuf, int pos, double sr) {
            MagPeak y;
            auto a = powBuf[pos - 1];
            auto b = powBuf[pos];
            auto c = powBuf[pos + 1];
            auto p = (a - c) / (a - 2 * b + c) * 0.5;
            auto h = b - 0.5 * (a - c) * p;
            if (h < 0) {
                // std::cerr << "peak interpolation returned negative power (probably not a true peak)" << std::endl;
                h = 0;
            }
            y.hz = static_cast<float>((pos + p) / fftSize * sr);
            y.pow = static_cast<float>(h);
            return y;
        }

        // compute spectral statistics of current frame, over given hz range
        template<int fftSize>
        static void computeStats(typename results::frame &dst, const double *powBuf, double sr,
                                 int minBin, int maxBin) {
            // static const double normScale = pow(1.0 / static_cast<double>(fftSize), 2);
            double meanLogMag = 0;
            double meanMag = 0;
            double meanPow = 0;
            double maxMag = 0;
            double meanHzW = 0; // weighted hz for centroid
            int n = 0;
            for (int bin = minBin; bin < maxBin; ++bin) {
                double pow = powBuf[bin];
                // double pow = powBuf[bin] * normScale;
                double mag = sqrt(pow);
                double hz = static_cast<double>(bin) / static_cast<double>(fftSize) * sr;
                meanLogMag += std::log(mag);
                meanMag += mag;
                meanHzW += mag * hz;
                meanPow += pow;
                maxMag = mag > maxMag ? mag : maxMag;
                n++;
            }
            double nscale = 1.0 / static_cast<double>(n);
            dst.papr = static_cast<float>((maxMag * maxMag) / meanPow);
            dst.centroid = static_cast<float>(meanHzW / meanMag);
            meanLogMag *= nscale;
            meanMag *= nscale;
            dst.flatness = static_cast<float>(exp(meanLogMag) / meanMag);
            dst.meanMag = static_cast<float>(meanMag);
            dst.maxMag = static_cast<float>(maxMag);
        }
    };
}

#endif // DOWSER_ANALYSIS_HPP
