#pragma once

namespace dowser
{
    struct ProcessConfig
    {
        float minPowDb, minHz, maxHz, minPersistence;
        int maxPeaksPerFrame;
    };
}