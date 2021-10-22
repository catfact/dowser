//
// Created by ezra on 10/21/21.
//
#ifndef DOWSER_PEAKS_HPP
#define DOWSER_PEAKS_HPP

#include <array>
#include <limits>

namespace dowser {
    namespace peaks {

        // find peaks in an array of values, returning a vector of indices
        template<typename T, int N>
        static std::vector<unsigned int> find_peaks(const T *src, T minValue,
                                                    unsigned int maxPeaks, unsigned int minBin, unsigned int maxBin,
                                                    T floorValue=0) {
            const size_t n = maxBin+1;
            std::vector<T> tmp(n);
            std::vector<bool> isPeak(n);

            for (unsigned int i = 0; i <n; ++i) {
                if ((src[i] < minValue) || (i < minBin)) {
                    tmp[i] = floorValue;
                } else {
                    tmp[i] = src[i];
                }
                isPeak[i] = false;
            }

            unsigned int count = 0;
            while (true) {
                // 1. find the position with the maximum value
                unsigned int maxIdx = minBin;
                double maxVal = floorValue;
                for (unsigned int i=minBin; i<maxBin; ++i) {
                    if (tmp[i] > maxVal) {
                        maxVal = tmp[i];
                        maxIdx = i;
                    }
                }
                if (maxVal <= minValue) {
                    // all the values are at/below threshold (maybe b/c we zeroed them), so we're done
                    // std::cout << "no more peak candidates; stopping search" << std::endl;
                    break;
                }

                // otherwise, yes this is a peak
                isPeak[maxIdx] = true;
                count++;

                // 2. zero out the region around the peak, "rolling down the slopes"
                double peakVal = tmp[maxIdx];
                tmp[maxIdx] = floorValue;
                // roll one way
                double val = peakVal;
                unsigned int idx = maxIdx;
                while (true) {
                    if (++idx >= maxBin) { break; }
                    double delta = tmp[idx] - val;
                    if (delta > 0) { break; }
                    val = tmp[idx];
                    tmp[idx] = floorValue;
                }
                // roll the other way
                val = peakVal;
                idx = maxIdx;
                while (true) {
                    if (--idx <= minBin) { break; }
                    double delta = tmp[idx] - val;
                    if (delta > 0) { break; }
                    val = tmp[idx];
                    tmp[idx] = floorValue;
                }

                // 3. repeat until enough peaks are found
                if (count >= maxPeaks) {
                    // std::cout << "found max peaks; stopping search" << std::endl;
                    break;
                }
            } // not-done loop
            std::vector<unsigned int> peakIdx;
            for (unsigned int i = 0; i < n; ++i) {
                if (isPeak[i]) {
                    peakIdx.push_back(i);
                }
            }

//            std::cout << "found peak indices:" << std::endl;
//            for (auto &idx: peakIdx) {
//                std::cout << idx << std::endl;
//            }
//            std::cout << std::endl;

            return peakIdx;
        } // find_peaks
    } // peaks
} // dowser
#endif //DOWSER_PEAKS_HPP
