#include <juce_core/juce_core.h>

#include "process.hpp"
#include "analysis.hpp"
#include "output.hpp"

// FIXME: make these into args
static constexpr int fftOrder = 11;
static constexpr int overlap = 2;

static void process(const juce::ArgumentList &args) {
    juce::File file = args.getExistingFileForOption("input");
    auto data = dowser::process<fftOrder, overlap>::perform(file, 1, 10000);
    auto results = dowser::analysis<fftOrder, overlap>::perform(std::move(data));
    dowser::output<fftOrder, overlap>::perform(std::move(results));
}

int main(int argc, char* argv[]) {

    juce::ConsoleApplication app;

    app.addHelpCommand ("--help|-h", "Usage:", true);
    app.addVersionCommand ("--version|-v", "MyApp version 1.2.3");


    app.addCommand ({ "--process",
                      "--process filename",
                      "process a soundfile",
                      "performs spectral analysis on given soundfile",
                      [] (const auto& args) { process (args); }});

    return app.findAndRunCommand (argc, argv);
}
