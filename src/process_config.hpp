#pragma once

#include "output.hpp"

namespace dowser
{
    struct ProcessConfig
    {
        dowser::output_format_t outputFormat;
        float minHz, maxHz, minPowDb, minPersistence;
        int maxPeaksPerFrame;
    };
}