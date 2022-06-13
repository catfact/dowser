#include <juce_core/juce_core.h>

#include "argh.h"

#include "process.hpp"
#include "analysis.hpp"
#include "output.hpp"
#include "process_config.hpp"

// FIXME: make these into args...
static constexpr int fftOrder = 14;
static constexpr int overlap = 2;

static void process(const juce::File infile, const juce::File outfile, const dowser::ProcessConfig &config) {

    std::cout << "performing import / STFT... (" <<infile.getFullPathName() << ")" << std::endl;
    auto data = dowser::process<fftOrder>::perform<overlap>(infile);

    std::cout << "performing spectral analysis..." << std::endl;
    auto results = dowser::analysis::perform<fftOrder>(std::move(data), config);

    std::cout << "performing output... (" << outfile.getFullPathName() << ")" << std::endl;
    dowser::output<dowser::output_format_t::supercollider> out;
    out.perform(std::move(results).get(), outfile);

    std::cout << "done." << std::endl;
}

int main(int argc, char* argv[]) {

    argh::parser cmdl(argv);
    dowser::ProcessConfig config;
    cmdl("--min_pow_db", -200.f) >> config.minPowDb;
    cmdl("--min_hz", 30.f) >> config.minHz;
    cmdl("--max_hz", 6000.f) >> config.maxHz;
    cmdl("--max_peaks", 48) >> config.maxPeaksPerFrame;
    cmdl("--min_persistence", 0.5) >> config.minPersistence;

    if (argc < 2) {
        std::cerr << "infile soundfile required" << std::endl;
        return 1;
    }
    juce::File infile = juce::File::getCurrentWorkingDirectory().getChildFile(argv[1]);
    juce::File outfile = argc > 2 ? juce::File::getCurrentWorkingDirectory().getChildFile(argv[2])
            : juce::File::getCurrentWorkingDirectory().getChildFile("dowser-output.scd");

    process(infile, outfile, config);
}
