#pragma once

#include "output_format.hpp"

namespace dowser
{
    struct ProcessConfig
    {
        output_format_t outputFormat;
        float minHz, maxHz, minPowDb, minPersistence;
        int maxPeaksPerFrame;
    };
}