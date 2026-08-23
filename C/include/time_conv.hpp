#pragma once

#include <math.h>
#include "lib/spice/cspice/include/SpiceUsr.h"
#include <string>


typedef struct
{
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
} Time;





double wrap_number(double N, double lower_bound, double upper_bound);

bool isLeapYear(int year);

std::string dateFromDateUTdt(const std::string& date_utc, double t);

int dayOfYear(Time& t);

Time timeFromDateUTC(const std::string& date_utc);

Time timeFromDateUTCdt(const std::string& date_utc, double t);

std::string dateUTCfromTime(Time time);

double convertUTtoJ0(const std::string& date_utc);

double convertUTtoJDfrac(const std::string& date_utc);

double epoch_to_GMST(const std::string& date_utc_0, double t);

double epoch_to_GMST(const std::string& date_utc_0);

