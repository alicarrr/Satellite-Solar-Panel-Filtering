#pragma once

#include <string>
#include "ephem.hpp"

typedef struct 
{
    double time_step;
    double t_final;
    double et_0;
    double tYear_0;
    std::string initialUtc;
    SolarSystemEphemeris* ephemeris;
    EarthFrameCache earthFrameCache;

} SimParams_t;
