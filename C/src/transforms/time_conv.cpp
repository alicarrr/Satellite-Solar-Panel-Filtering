#include "time_conv.hpp"
#include "vectorMaths.hpp"
#include "constants.hpp"

#include <stdexcept>
#include <iostream>
#include <string>
#include <stdio.h>
#include <iomanip>

bool isLeapYear(int year)
{
    return (year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0));
}




std::string dateFromDateUTdt(const std::string& date_utc, double t)
{
    SpiceDouble et_date;
    str2et_c(date_utc.c_str(), &et_date);
    

    SpiceDouble et = et_date + t;

    SpiceChar out[64];
    timout_c(et, "YYYY-MM-DD HR:MN:SC.###### ::UTC ::RND", 64, out);

    std::string date = std::string(out);


    return date;
}


Time timeFromDateUTC(const std::string& date_utc)
{
    Time time;
    
    int year = 0;
    int month = 0;
    int day = 0;
    int hour = 0;
    int minute = 0;
    double second = 0.0;

    if (std::sscanf(date_utc.c_str(), "%d-%d-%d %d:%d:%lf", &year, &day, &month, &hour, &minute, &second) != 6) {
        throw std::runtime_error("Failed to parse UTC time from SPICE for decimal-year conversion");
    }
    

    time.year = year;
    time.day = day;
    time.month = month;
    time.hour = hour;
    time.minute = minute;
    time.second = second;


    return time;

}


int dayOfYear(Time& t) {
    static const int days_before_month[12] = {
        0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334
    };
    bool is_leap = (t.year % 4 == 0 && t.year % 100 != 0) || (t.year % 400 == 0);
    int doy = days_before_month[t.month - 1] + t.day;
    if (is_leap && t.month > 2) doy += 1;
    return doy;
}

std::string dateUTCfromTime(Time time)
{
    char date[64];
    std::sprintf(date, "%d-%d-%d %d:%d:%lf", time.year, time.month, time.day, time.hour, time.minute, time.second);
    return std::string(date);
}



Time timeFromDateUTCdt(const std::string& date_utc, double t)
{
    
    std::string date = dateFromDateUTdt(date_utc, t);

    return timeFromDateUTC(date);

    
}


double convertUTtoJ0(const std::string& date_utc)
{
    Time time = timeFromDateUTC(date_utc);

    double J0 = 367.0 * time.year - floorf64(7.0 * (time.year + floorf64((time.month + 9.0) / 12.0)) / 4.0) + floorf64(275.0 * time.month / 9.0) + time.day + 1721013.5;

    return J0;
}


double convertUTtoJDfrac(const std::string& date_utc)
{
    Time time = timeFromDateUTC(date_utc);

    double JD_UT_frac = (time.hour * hour + time.minute * minute + time.second * second);

    return JD_UT_frac;
}


double epoch_to_GMST(const std::string& date_utc_0, double t)
{
    std::string date_utc = dateFromDateUTdt(date_utc_0, t);

    double J0 = convertUTtoJ0(date_utc);
    double JD_frac = convertUTtoJDfrac(date_utc);

    double T0 = (J0 - J2000_JD) / 36525;

    double angular_GMST_0 = 100.4606184 + 36000.77004 * T0 + 0.000387933 * T0*T0 - 2.583e-8 * T0*T0*T0;

    double angular_GMST = angular_GMST_0 + 360.98564724 * JD_frac;

    double radian_GMST = DEG2RAD * angular_GMST;
    double GMST = radian_GMST / angular_vel_earth;

    return GMST;
}


double epoch_to_GMST(const std::string& date_utc)
{

    double J0 = convertUTtoJ0(date_utc);
    double JD_frac = convertUTtoJDfrac(date_utc);

    double T0 = (J0 - J2000_JD) / 36525;

    double angular_GMST_0 = 100.4606184 + 36000.77004 * T0 + 0.000387933 * T0*T0 - 2.583e-8 * T0*T0*T0;

    double angular_GMST = angular_GMST_0 + 360.98564724 * JD_frac;

    double radian_GMST = DEG2RAD * angular_GMST;
    double GMST = radian_GMST / angular_vel_earth;

    return GMST;
}