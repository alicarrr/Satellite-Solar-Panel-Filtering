#pragma once

#include <limits>
#include <memory>
#include <string>
#include <eigen3/Eigen/Dense>

#include "time_conv.hpp"

namespace GeographicLib {
class GravityModel;
class MagneticModel;
}

// ---------------------------------------------------------------------------
// EarthFrameCache
//
// Stores the J2000 <-> ITRF93 rotation matrices and the decimal year at a
// particular simulation time.  The simSeconds field records when the cache
// was built so callers (and the class itself) can decide whether a refresh
// is needed.
// ---------------------------------------------------------------------------
struct EarthFrameCache
{
    Eigen::Matrix3d R_j2e;                                         // J2000  -> ECEF (ITRF93)
    Eigen::Matrix3d R_e2j;                                         // ECEF   -> J2000
    double tYear      = std::numeric_limits<double>::quiet_NaN();  // decimal year of cache
    double simSeconds = std::numeric_limits<double>::quiet_NaN();  // sim time of cache
};

double decimalYearFromEtUtc(const std::string& date_0, double t);

class SolarSystemEphemeris {
public:
    SolarSystemEphemeris(const std::string& lskPath,
                         const std::string& spkPath,
                         const std::string& initialUtc);
    ~SolarSystemEphemeris();

    void setInitialUtc(const std::string& utc);
    double initialEt() const;
    double etAtSeconds(double simSeconds);

    Eigen::VectorXd positionM(const std::string& target,
                              const std::string& observer,
                              double simSeconds,
                              const std::string& frame  = "J2000",
                              const std::string& abcorr = "NONE");

    Eigen::VectorXd stateMAndMps(const std::string& target,
                                 const std::string& observer,
                                 double simSeconds,
                                 const std::string& frame  = "J2000",
                                 const std::string& abcorr = "NONE");

    Eigen::VectorXd sunRelativeToEarthM(double simSeconds,
                                        const std::string& frame  = "J2000",
                                        const std::string& abcorr = "NONE");

    std::string utcStringAtSeconds(double simSeconds);

    void loadEarthFieldModels(const std::string& earthPckPath,
                              const std::string& gravityModelName  = "egm2008",
                              const std::string& gravityModelPath  = "",
                              int gravityNmax                      = 100,
                              int gravityMmax                      = 100,
                              const std::string& magneticModelName = "wmm2025",
                              const std::string& magneticModelPath = "",
                              int magneticNmax                     = 100,
                              int magneticMmax                     = 100);

    Eigen::Vector3d eciToEcef(const Eigen::Vector3d& vecEci,  double simSeconds);
    Eigen::Vector3d ecefToEci(const Eigen::Vector3d& vecEcef, double simSeconds);

    Eigen::Vector3d earthGravityFieldJ2000Mps2(const Eigen::Vector3d& rEarthToPointJ2000M,
                                               const std::string& date_utc_0,
                                               double simSeconds);

    Eigen::Vector3d earthMagneticFieldJ2000Tesla(const Eigen::Vector3d& rEarthToPointJ2000M,
                                                 double simSeconds);

private:
    double      et0_;
    std::string initialUtc_;

    std::unique_ptr<GeographicLib::GravityModel>  gravityModel_;
    std::unique_ptr<GeographicLib::MagneticModel> magneticModel_;

    EarthFrameCache frameCache_;

    EarthFrameCache& getOrRefreshFrameCache(double simSeconds);

    static void checkSpiceError();
};

// ---------------------------------------------------------------------------
// Free-function helper retained for callers that manage their own cache.
// The returned cache has simSeconds set so it can be compared directly.
// ---------------------------------------------------------------------------
EarthFrameCache makeEarthFrameCache(SolarSystemEphemeris& eph,
                                    double simTimeSeconds,
                                    double tYear_0);



Eigen::Vector3d get_sun_location_eci(double t, SolarSystemEphemeris* ephem);
