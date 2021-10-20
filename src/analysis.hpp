//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_ANALYSIS_HPP
#define DOWSER_ANALYSIS_HPP

#include <memory>

#include "process.hpp"

namespace dowser {

    class analysis {

    public:

        typedef std::pair<float, float> peak_t;
        struct results {

            struct frame {
                std::vector<peak_t> peaks;
                float papr{};
                float centroid{};
                float flatness{};
                float meanPower{};
            };

            std::vector<struct frame> frames;
        };


        template<int fftOrder>
        static std::unique_ptr<struct results>
        perform(std::unique_ptr<typename process<fftOrder>::data> data,
                float minHz, float maxHz) {
            static constexpr int fftSize = 1 << fftOrder;
            static constexpr int numRealBins = fftSize / 2 + 1;

            auto res = std::make_unique<struct results>();

            auto sr = data->sampleRate;
            int minBin = std::max(1, static_cast<int>(minHz * fftSize / sr));
            int maxBin = std::min(numRealBins - 2, static_cast<int>(maxHz * fftSize / sr));
            unsigned int maxPeaks = 16;
            double minMag = 0.0001;
            for (auto &frame: data->specPowFrames) {
                typename results::frame resFrame;
                resFrame.peaks = findPeaks<fftSize>(frame.data(), sr, maxPeaks, minMag, minBin, maxBin);
                computeStats<fftSize>(resFrame, frame.data(), sr, minBin, maxBin);
                res->frames.push_back(resFrame);
            }

            return res;
        }


    private:
        // return vector of peaks for the current magnitudes

        template<int fftSize>
        static std::vector<peak_t> findPeaks(const double *powBuf, double sr,
                                             unsigned int maxPeaks, double minPow,
                                             int minBin, int maxBin) {
            std::vector<peak_t> y; // results vector


            // find peak indices
            std::vector<int> idx;
            for (int bin = minBin; bin <= maxBin; ++bin) {
                if (powBuf[bin] > minPow) {
                    // assumption: bin-1 and bin+1 are in range
                    if (powBuf[bin] > powBuf[bin - 1] && powBuf[bin] > powBuf[bin + 1]) {
                        idx.push_back(bin);
                    }
                }
            }
            // interpolate true locations / magnitudes
            y.reserve(idx.size());
            for (auto &pos: idx) {
                auto peak = refinePeak<fftSize>(powBuf, pos, sr);
                y.push_back(peak);
                // .. could compute peak "width" here but not sure it's useful
            }
            // sort by peak power
            std::sort(y.begin(), y.end(), [](peak_t a, peak_t b) { return a.second > b.second; });
            // retain only highest N peaks
            if (y.size() > maxPeaks) {
                y = {y.begin(), y.begin() + maxPeaks};
            }
            return y;
        }

        // approximate true peak location by quadratic fit
        // assumption: pos-1, pos+1 are in range
        template<int fftSize>
        static peak_t refinePeak(const double *powBuf, int pos, double sr) {
            peak_t y;
            auto a = powBuf[pos - 1];
            auto b = powBuf[pos];
            auto c = powBuf[pos + 1];
            auto p = (a - c) / (a - 2 * b + c) * 0.5;
            auto h = b - 0.5 * (a - c) * p;
            y.first = static_cast<float>((pos + p) / fftSize * sr);
            y.second = static_cast<float>(h);
            return y;
        }


        // compute spectral flatness of current frame, over given hz range
        template<int fftSize>
        static void computeStats(typename results::frame &dst, const double *powBuf, double sr,
                                 int minBin, int maxBin) {
            static const double normScale = pow(1.0 / static_cast<double>(fftSize), 2);
            double meanLogMag = 0;
            double meanMag = 0;
            double meanPow = 0;
            double maxPow = 0;
            double meanHzW = 0; // weighted hz for centroid
            int n = 0;
            for (int bin = minBin; bin < maxBin; ++bin) {
                double pow = powBuf[bin] * normScale;
                double mag = sqrt(pow);
                double hz = static_cast<double>(bin) / static_cast<double>(fftSize) * sr;
                meanLogMag += std::log(mag);
                meanMag += mag;
                meanHzW += mag * hz;
                meanPow += pow;
                maxPow = pow > maxPow ? pow : maxPow;
                n++;
            }
            double nscale = 1.0 / static_cast<double>(n);
            meanPow *= nscale;
            dst.papr = static_cast<float>(maxPow / meanPow);
            dst.centroid = static_cast<float>(meanHzW / meanMag);
            meanLogMag *= nscale;
            meanMag *= nscale;
            dst.flatness = static_cast<float>(exp(meanLogMag) / meanMag);
            dst.meanPower = static_cast<float>(meanPow);
        }
    };
}

#endif //DOWSER_ANALYSIS_HPP
