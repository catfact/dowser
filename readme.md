## dowser

`dowser` is (presently) a quick and simple low-level utility for performing spectral analysis of sound files.

### usage

`dowser [in] [out] <options>`

options include:

- `--minhz`, lower frequency bound of analysis, in hz
- `--maxhz`, upper frequency bound of analysis, in hz
- `--minmag`, peak magnitude threshold
- `--maxpeaks`, max number of peaks per frame

the output is a supercollider script, defining a list of dictionaries. each dictionary contains data for a single spectral analysis frame.

per-frame measurements include:

- `peaks`: list of events with `(hz, mag)` keys, listing spectral peaks for the given frame
- `papr`: peak-to-average-power ratio, a measure of "tonal-ness"
- `flatness`: AKA weiner entropy, geometric mean / arithmetic mean. another measure of "tonalness"
- `centroid`: spectral centroid, a measure of "brightness"
- `fluxPositive`: spectral flux, in a positive direction
- `fluxNegative`: spectral flux, in a negative direction

### requirements

`cmake` and a suitable c++ compilation toolset (requires c++17)

### build

1. fetch dependencies (JUCE, which is lame but so it goes):

`git submodule update --init --recursive`

2. create build folder:

`mkdir cmake-build-release`

3. configure and build:

`cd cmake-build-release`

`cmake -DCMAKE_BUILD_TYPE=Release ..` (note the `..`)

`cmake --build .`

the freshly compiled executable should then be located at:
`cmake_build_release/dowser_artefacts/Release/dowser` 
(or `.../dowser.exe` on windows)

### TODO

- other output formats (csv, python, etc)
- expose more parameters (hz range, max peaks per frame, FFT config)
- more analyses? (crest factor? autocorrelation peaks?)
- drag/drop GUI option
- maybe some higher-level analysis (e.g. clustering)
