#include <fstream>
#include <iterator>
#include <filesystem>
#include <iostream>
#include "peaks_watershed.hpp"

const char path[] = "/Users/emb/code/dowser/test_signals/spectrum_frame.bin";
constexpr size_t valueSize = sizeof(double);
constexpr int N = 256;

int main()
{
    union {
        char buf[8];
        double value;
    } u;
    char buf[8];

    std::vector<double> X;

    auto fileSize = std::filesystem::file_size(path);
    assert(N <= fileSize/valueSize);
    X.reserve(N);

    std::ifstream input(path, std::ios::binary );

    for (int i=0; i<N; ++i) {
        input.read(buf, valueSize);
        u.buf[0] = buf[7];
        u.buf[1] = buf[6];
        u.buf[2] = buf[5];
        u.buf[3] = buf[4];
        u.buf[4] = buf[3];
        u.buf[5] = buf[2];
        u.buf[6] = buf[1];
        u.buf[7] = buf[0];

        X.push_back(u.value);
    }

//    std::cout << "data = [" << std::endl;
//    for (auto x: X) {
//        std::cout << x << ", ";
//    }
//    std::cout << std::endl << "];" << std::endl;

    auto peaks = dowser::peaks::Watershed<double>::findPeaks(X.data(), N);

    int n = peaks.size();

    std::cout << "x = [" << std::endl;
    for (int i=0; i<n; ++i) {
        std::cout << peaks[i].index << ", ";
    }
    std::cout << std::endl << "]; " << std::endl << std::endl;

    std::cout << "y = [" << std::endl;
    for (int i=0; i<n; ++i) {
        std::cout << peaks[i].value << ", ";
    }
    std::cout << std::endl << "]; " << std::endl << std::endl;

    std::cout << "z = [" << std::endl;
    for (int i=0; i<n; ++i) {
        std::cout << peaks[i].persistence << ", ";
    }
    std::cout << std::endl << "]; " << std::endl;
//
//    for (auto &p: peaks) {
//        std::cout << "idx = "<<p.index << "; y = " << p.value << "; p = " << p.persistence << std::endl;
//    }
}