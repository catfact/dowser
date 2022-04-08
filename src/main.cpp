#include <juce_core/juce_core.h>

#include "process.hpp"
#include "analysis.hpp"
#include "output.hpp"

// FIXME: make these into args...
static constexpr int fftOrder = 14;
static constexpr int overlap = 2;

static void process(const juce::File infile, const juce::File outfile) {
    // FIXME: definitely make these into args
    static const float minHz = 10.f;
    static const float maxHz = 10000.f;

    std::cout << "performing import / STFT... (" <<infile.getFullPathName() << ")" << std::endl;
    auto data = dowser::process<fftOrder>::perform<overlap>(infile);

    std::cout << "performing spectral analysis..." << std::endl;
    auto results = dowser::analysis::perform<fftOrder>(std::move(data), minHz, maxHz);

    std::cout << "performing output... (" << outfile.getFullPathName() << ")" << std::endl;
    dowser::output::perform(std::move(results), outfile);

    std::cout << "done." << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "infile soundfile required" << std::endl;
        return 1;
    }
    juce::File infile = juce::File::getCurrentWorkingDirectory().getChildFile(argv[1]);
    juce::File outfile = argc > 2 ? juce::File::getCurrentWorkingDirectory().getChildFile(argv[2])
            : juce::File::getCurrentWorkingDirectory().getChildFile("dowser-output.scd");

    process(infile, outfile);
}
