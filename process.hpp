//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_PROCESS_HPP
#define DOWSER_PROCESS_HPP

#include <complex>
#include <memory>

#include <juce_core/juce_core.h>
#include <juce_dsp/juce_dsp.h>

#include "window.hpp"

namespace dowser {
    template<int fftOrder, int overlap>
    class process {

    public:
        static constexpr int fftSize = 1 << fftOrder;
        static constexpr int hopSize = fftSize / overlap;
        static constexpr int numRealBins = fftSize / 2 + 1;

        struct data {
            typedef std::vector<float> spec_frame;
            std::vector<spec_frame> frames;
        };

        static std::unique_ptr<struct data> perform(const juce::File& soundfile, float minHz, float maxHz ) {
            juce::AudioFormatManager formatManager;
            formatManager.registerBasicFormats();
            auto reader = formatManager.createReaderFor(soundfile.createInputStream());
            double sampleRate = reader->sampleRate;

            // analysis window
            std::array<float, fftSize> win = window::data<window::shape_t::HANN, fftSize>();

            // working buffer for FFT
            // (format is required by juce API)
            std::array<float, fftSize*2> buf;
            for (auto &x: buf) { x = 0.f; }

            // working buffer for squared bin magnitudes
            std::array<double, numRealBins> binMag2;

            int minBin = std::max(1, static_cast<int>(minHz * fftSize / sampleRate));
            int maxBin = std::min(numRealBins - 2, static_cast<int>(maxHz * fftSize / sampleRate));

            auto data = std::make_unique<struct data>();
            juce::dsp::FFT fft(fftOrder);

            juce::int64 numSamples = reader->lengthInSamples;
            int offset = 0;
            int maxOffset = findMaxFftFrameOffset(numSamples);


            auto *src = new float[numSamples];
            //-------------------
            /// TODO:
            /// read the entire soundfile into `src`
            //-------------------

            while (offset < maxOffset)
            {
                const float *x = &(src[offset]);
                for (int i = 0; i < fftSize; ++i)
                {
                    buf[i] = x[i] * win[i];
                }

                fft.performRealOnlyForwardTransform(buf.data(), true);
                // result is placed in the IO buffer with real/imag interleaved
                x = &(buf[0]);
                for (int i=0; i<numRealBins; ++i) {
                    float re = *x++;
                    float im = *x++;
                    binMag2[i] = re*re + (im*im);
                }
                offset += hopSize;
            }

            delete[] src;
            return std::move(data);
        }

    private:

        static juce::int64 findMaxFftFrameOffset(juce::int64 numSampleFrames)
        {
            juce::int64 used = 0;
            juce::int64 idx = 0;
            juce::int64 idx0 = 0;
            while (used < numSampleFrames)
            {
                idx0 = idx;
                idx += hopSize;
                used = idx + fftSize;
            }
            return idx0;
        }


    };
}
#endif //DOWSER_PROCESS_HPP
