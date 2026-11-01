#include "ephem.hpp"

#include <stdio.h>
#include <iostream>

#include <stdexcept>
#include <cstdio>
#include <cmath>
#include <vector>

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/MagneticModel.hpp>

#include "spice/cspice/include/SpiceUsr.h"

#include "time_conv.hpp"
#include "transforms.hpp"
#include "constants.hpp"
#include <string>


double decimalYearFromEtUtc(const std::string& date_0, double t)
{
    Time time = timeFromDateUTCdt(date_0, t);

    const double secondsOfDay = 3600.0 * time.hour + 60.0 * time.minute + time.second;
    const double daysInYear   = isLeapYear(time.year) ? 366.0 : 365.0;

    return static_cast<double>(time.year)
         + (static_cast<double>(dayOfYear(time) - 1) + secondsOfDay / 86400.0) / daysInYear;
}

namespace {

void requireFieldModel(const void* ptr, const char* name)
{
    if (!ptr) {
        throw std::runtime_error(std::string(name) +
                                 " is not loaded. Call loadEarthFieldModels(...) first.");
    }
}

} // namespace


// ---------------------------------------------------------------------------
// Construction / destruction
// ---------------------------------------------------------------------------

SolarSystemEphemeris::SolarSystemEphemeris(const std::string& lskPath,
                                           const std::string& spkPath,
                                           const std::string& initialUtc)
    : et0_(0.0), initialUtc_(initialUtc)
{
    furnsh_c(lskPath.c_str());
    checkSpiceError();

    furnsh_c(spkPath.c_str());
    checkSpiceError();

    str2et_c(initialUtc_.c_str(), &et0_);
    checkSpiceError();
}

SolarSystemEphemeris::~SolarSystemEphemeris()
{
    kclear_c();
}


// ---------------------------------------------------------------------------
// Time helpers
// ---------------------------------------------------------------------------

void SolarSystemEphemeris::setInitialUtc(const std::string& utc)
{
    initialUtc_ = utc;
    str2et_c(initialUtc_.c_str(), &et0_);
    checkSpiceError();

    // Invalidate the frame cache – it was built for the old epoch.
    frameCache_.simSeconds = std::numeric_limits<double>::quiet_NaN();
}

double SolarSystemEphemeris::initialEt() const
{
    return et0_;
}

double SolarSystemEphemeris::etAtSeconds(double simSeconds)
{
    return et0_ + simSeconds;
}


// ---------------------------------------------------------------------------
// Internal frame cache
//
// Rebuilds the J2000 <-> ITRF93 rotation matrices and decimal year whenever
// simSeconds differs from the cached value.  An exact floating-point
// comparison is intentional: the simulation typically calls these functions
// with the same time step value repeatedly, so bit-identical reuse is both
// correct and desirable.
// ---------------------------------------------------------------------------

EarthFrameCache& SolarSystemEphemeris::getOrRefreshFrameCache(double simSeconds)
{
    if (frameCache_.simSeconds == simSeconds) {
        return frameCache_;   // cache hit – nothing to recompute
    }

    const double et = etAtSeconds(simSeconds);

    SpiceDouble R_spice[3][3];
    pxform_c("J2000", "ITRF93", et, R_spice);
    checkSpiceError();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            frameCache_.R_j2e(i, j) = R_spice[i][j];

    frameCache_.R_e2j     = frameCache_.R_j2e.transpose();
    frameCache_.tYear     = decimalYearFromEtUtc(initialUtc_, simSeconds);
    frameCache_.simSeconds = simSeconds;

    return frameCache_;
}


// ---------------------------------------------------------------------------
// ECI (J2000) <-> ECEF (ITRF93) conversions
// ---------------------------------------------------------------------------

Eigen::Vector3d SolarSystemEphemeris::eciToEcef(const Eigen::Vector3d& vecEci,
                                                 double simSeconds)
{
    const EarthFrameCache& c = getOrRefreshFrameCache(simSeconds);
    return c.R_j2e * vecEci;
}

Eigen::Vector3d SolarSystemEphemeris::ecefToEci(const Eigen::Vector3d& vecEcef,
                                                 double simSeconds)
{
    const EarthFrameCache& c = getOrRefreshFrameCache(simSeconds);
    return c.R_e2j * vecEcef;
}


// ---------------------------------------------------------------------------
// SPICE ephemeris queries
// ---------------------------------------------------------------------------

Eigen::VectorXd SolarSystemEphemeris::positionM(const std::string& target,
                                                 const std::string& observer,
                                                 double simSeconds,
                                                 const std::string& frame,
                                                 const std::string& abcorr)
{
    const double et = etAtSeconds(simSeconds);

    SpiceDouble pos_km[3];
    SpiceDouble lt = 0.0;

    spkpos_c(target.c_str(), et, frame.c_str(), abcorr.c_str(),
             observer.c_str(), pos_km, &lt);
    checkSpiceError();

    Eigen::VectorXd pos_m(3);
    pos_m(0) = pos_km[0] * 1000.0;
    pos_m(1) = pos_km[1] * 1000.0;
    pos_m(2) = pos_km[2] * 1000.0;

    return pos_m;
}

Eigen::VectorXd SolarSystemEphemeris::stateMAndMps(const std::string& target,
                                                    const std::string& observer,
                                                    double simSeconds,
                                                    const std::string& frame,
                                                    const std::string& abcorr)
{
    const double et = etAtSeconds(simSeconds);

    SpiceDouble state_km_kmps[6];
    SpiceDouble lt = 0.0;

    spkezr_c(target.c_str(), et, frame.c_str(), abcorr.c_str(),
             observer.c_str(), state_km_kmps, &lt);
    checkSpiceError();

    Eigen::VectorXd state_m_mps(6);
    for (int i = 0; i < 6; ++i)
        state_m_mps(i) = state_km_kmps[i] * 1000.0;

    return state_m_mps;
}

Eigen::VectorXd SolarSystemEphemeris::sunRelativeToEarthM(double simSeconds,
                                                           const std::string& frame,
                                                           const std::string& abcorr)
{
    return positionM("SUN", "EARTH", simSeconds, frame, abcorr);
}

std::string SolarSystemEphemeris::utcStringAtSeconds(double simSeconds)
{
    const double et = etAtSeconds(simSeconds);

    SpiceChar out[64];
    timout_c(et, "YYYY-MM-DD HR:MN:SC.###### ::UTC ::RND", 64, out);
    checkSpiceError();

    return std::string(out);
}


// ---------------------------------------------------------------------------
// Field-model loading
// ---------------------------------------------------------------------------

void SolarSystemEphemeris::loadEarthFieldModels(const std::string& earthPckPath,
                                                 const std::string& gravityModelName,
                                                 const std::string& gravityModelPath,
                                                 int gravityNmax,
                                                 int gravityMmax,
                                                 const std::string& magneticModelName,
                                                 const std::string& magneticModelPath,
                                                 int magneticNmax,
                                                 int magneticMmax)
{
    if (!earthPckPath.empty()) {
        furnsh_c(earthPckPath.c_str());
        checkSpiceError();
    }

    gravityModel_ = std::make_unique<GeographicLib::GravityModel>(
        gravityModelName, gravityModelPath, gravityNmax, gravityMmax);

    magneticModel_ = std::make_unique<GeographicLib::MagneticModel>(
        magneticModelName,
        magneticModelPath,
        GeographicLib::Geocentric::WGS84(),
        magneticNmax,
        magneticMmax);
}


// ---------------------------------------------------------------------------
// Earth gravity field (J2000)
//
// Uses the internal frame cache.  A SPICE rotation-matrix call is avoided
// whenever simSeconds matches the previously cached value.
// ---------------------------------------------------------------------------

Eigen::Vector3d SolarSystemEphemeris::earthGravityFieldJ2000Mps2(
    const Eigen::Vector3d& rEarthToPointJ2000M,
    const std::string& /*date_utc_0*/,
    double simSeconds)
{
    requireFieldModel(gravityModel_.get(), "Gravity model");

    const EarthFrameCache& c = getOrRefreshFrameCache(simSeconds);

    const Eigen::Vector3d rItrf = c.R_j2e * rEarthToPointJ2000M;

    double GX = 0.0, GY = 0.0, GZ = 0.0;
    gravityModel_->V(rItrf.x(), rItrf.y(), rItrf.z(), GX, GY, GZ);

    const Eigen::Vector3d gItrf(GX, GY, GZ);
    return c.R_e2j * gItrf;
}


// ---------------------------------------------------------------------------
// Earth magnetic field (J2000)
//
// Uses the internal frame cache in the same way as the gravity method.
// ---------------------------------------------------------------------------

Eigen::Vector3d SolarSystemEphemeris::earthMagneticFieldJ2000Tesla(
    const Eigen::Vector3d& rEarthToPointJ2000M,
    double simSeconds)
{
    requireFieldModel(magneticModel_.get(), "Magnetic model");

    const EarthFrameCache& c = getOrRefreshFrameCache(simSeconds);

    const Eigen::Vector3d rItrf = c.R_j2e * rEarthToPointJ2000M;

    double BX_nT = 0.0, BY_nT = 0.0, BZ_nT = 0.0;
    double BXt   = 0.0, BYt   = 0.0, BZt   = 0.0;

    magneticModel_->FieldGeocentric(
        c.tYear,
        rItrf.x(), rItrf.y(), rItrf.z(),
        BX_nT, BY_nT, BZ_nT,
        BXt,   BYt,   BZt);

    const Eigen::Vector3d Bitrf_T = 1e-9 * Eigen::Vector3d(BX_nT, BY_nT, BZ_nT);
    return c.R_e2j * Bitrf_T;
}


// ---------------------------------------------------------------------------
// SPICE error checking
// ---------------------------------------------------------------------------

void SolarSystemEphemeris::checkSpiceError()
{
    if (failed_c()) {
        SpiceChar msg[1841];
        getmsg_c("LONG", 1840, msg);
        reset_c();
        throw std::runtime_error(std::string("SPICE error: ") + msg);
    }
}


// ---------------------------------------------------------------------------
// Free-function helper – retained for callers that manage their own cache.
//
// tYear_0 is the decimal year corresponding to simTimeSeconds = 0 for the
// ephemeris; pass the same value you would have used before.  The returned
// cache now includes simSeconds so it can be compared directly against the
// class-internal cache.
// ---------------------------------------------------------------------------

EarthFrameCache makeEarthFrameCache(SolarSystemEphemeris& eph,
                                    double simTimeSeconds,
                                    double tYear_0)
{
    EarthFrameCache c;

    const double et = eph.etAtSeconds(simTimeSeconds);

    SpiceDouble R_spice[3][3];
    pxform_c("J2000", "ITRF93", et, R_spice);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            c.R_j2e(i, j) = R_spice[i][j];

    c.R_e2j     = c.R_j2e.transpose();
    c.tYear     = tYear_0 + simTimeSeconds / year;
    c.simSeconds = simTimeSeconds;

    return c;
}


Eigen::Vector3d get_sun_location_eci(double t, SolarSystemEphemeris* ephem)
{
    Eigen::Vector3d s_sol = ephem->sunRelativeToEarthM(t);

    return s_sol;
}