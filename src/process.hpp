//
// Created by emb on 10/19/2021.
//

#ifndef DOWSER_PROCESS_HPP
#define DOWSER_PROCESS_HPP

#include <memory>
#include <numeric>

#include <juce_core/juce_core.h>
#include <juce_dsp/juce_dsp.h>

#include "window.hpp"

namespace dowser {

    template<int fftOrder>
    class process {

    public:
        static constexpr int fftSize = 1 << fftOrder;
        static constexpr int fftPadRatio = 2;
        static constexpr int fftSizePadded = fftSize * fftPadRatio;
        static constexpr int numRealBins = fftSizePadded / 2 + 1;

        struct analysisData {
            double sampleRate;
            std::vector<std::array<double, numRealBins>> specPowFrames;
            std::vector<std::array<double, numRealBins>> autoCorrFrames;
        };

        template<int overlap>
        static std::unique_ptr<struct analysisData> perform(const juce::File& soundfile) {
            static constexpr int hopSize = fftSize / overlap;
            using juce::int64;

            juce::AudioFormatManager formatManager;
            formatManager.registerBasicFormats();
            auto reader = formatManager.createReaderFor(soundfile.createInputStream());

            // analysis window
            std::array<float, fftSize> win = window::data<window::shape_t::HANN, fftSize>();
            const double fftNormScale = fftPadRatio * 1.0 / std::accumulate(win.begin(), win.end(), 0.0);

            // working buffer for FFT
            // (interleaved format is required by juce API)
            std::array<float, fftSizePadded*2> buf;
            for (auto &x: buf) { x = 0.f; }

            // output buffer for squared bin magnitudes
            std::array<double, numRealBins> binPow;
            // output buffer for autocorrelation
            std::array<double, numRealBins> binAc;
            
            int64 numSampleFrames = reader->lengthInSamples;
            int64 offset = 0;
            int64 maxOffset = findMaxFftFrameOffset<hopSize>(numSampleFrames);

            //------------------------------------
            //-- perform the STFT
            juce::dsp::FFT fft(fftOrder);
            

            auto results = std::make_unique<struct analysisData>();
            results->sampleRate = reader->sampleRate * fftPadRatio;

            while (offset < maxOffset)
            {
                // read an FFT frame's worth of input from the soundfile
                /// FIXME: we are doing some unnecessary disk access by reading each frame in full.
                /// would be better to use a circular buffer and only pull one hop's worth
                float* src = buf.data();
                reader->read(&src, 1, offset, fftSize);

                unsigned int i=0;
                while (i<fftSize) {
                    buf[i] *= win[i];
                    i++;
                }
                while (i<buf.size()) {
                    buf[i++] = 0.f;
                }

                fft.performRealOnlyForwardTransform(buf.data(), true);
                // result is placed in the IO buffer with real/imag interleaved
                float *pRe = &(buf[0]);
                float *pIm = &(buf[1]);

                for (i=0; i<numRealBins; ++i) {
                    double re = *pRe * fftNormScale;
                    double im = *pIm * fftNormScale;
                    double pow =(re * re) + (im * im);
                    binPow[i] = pow;
#if INCLUDE_AUTOCORR
                    // FIXME: echh.. don't want windowing on the signal for AC.
                    // also need different circularity correction, etc;
                    /// anyway it needs work, so disabling it for now.

                    // store calculated power for autocorrelation
                    *pRe = static_cast<float>(pow);
                    *pIm = 0.f;
#endif
                    pRe += 2;
                    pIm += 2;
                }

#if INCLUDE_AUTOCORR
                // perform autocorrelation
                fft.performRealOnlyInverseTransform(buf.data());
                // (i believe the AC data is now non-interleaved in the I/O buf?)
                for (i=0; i<numRealBins; ++i) {
                    binAc[(size_t)i] = buf.data()[i];
                }
                results->autoCorrFrames.push_back(binAc);
#endif
                results->specPowFrames.push_back(binPow);
                offset += hopSize;
            }
            delete reader;
            return std::move(results);
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