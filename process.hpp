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

    template<int fftOrder>
    class process {

    public:
        static constexpr int fftSize = 1 << fftOrder;
        static constexpr int numRealBins = fftSize / 2 + 1;

        struct data {
            double sampleRate;
            std::vector<std::array<double, numRealBins>> specMagFrames;
        };

        template<int overlap>
        static std::unique_ptr<struct data> perform(const juce::File& soundfile) {
            static constexpr int hopSize = fftSize / overlap;
            using juce::int64;
            juce::AudioFormatManager formatManager;
            formatManager.registerBasicFormats();
            auto reader = formatManager.createReaderFor(soundfile.createInputStream());

            // analysis window
            std::array<float, fftSize> win = window::data<window::shape_t::HANN, fftSize>();

            // working buffer for FFT
            // (format is required by juce API)
            std::array<float, fftSize*2> buf;
            for (auto &x: buf) { x = 0.f; }

            // working buffer for squared bin magnitudes
            std::array<double, numRealBins> binMag2;

            int64 numSampleFrames = reader->lengthInSamples;
            int64 offset = 0;
            int64 maxOffset = findMaxFftFrameOffset<hopSize>(numSampleFrames);

            //-- perform STFT
            auto data = std::make_unique<struct data>();
            data->sampleRate = reader->sampleRate;

            juce::dsp::FFT fft(fftOrder);

            while (offset < maxOffset)
            {
                // read an FFT frame's worth of input from the soundfile
                /// FIXME: we are doing some unnecessary disk access by reading each frame in full.
                /// would be better to use a circular buffer and only pull one hop's worth
                float* src = buf.data();
                reader->read(&src, 1, offset, fftSize);

                for (int i = 0; i < fftSize; ++i)
                {
                    buf[i] *= win[i];
                }

                fft.performRealOnlyForwardTransform(buf.data(), true);
                // result is placed in the IO buffer with real/imag interleaved
                const float *x = &(buf[0]);
                for (int i=0; i<numRealBins; ++i) {
                    float re = *x++;
                    float im = *x++;
                    binMag2[i] = re*re + (im*im);
                }
                data->specMagFrames.push_back(binMag2);
                offset += hopSize;
            }
            delete reader;
            return std::move(data);
        }

    private:
        template<int hopSize>
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
