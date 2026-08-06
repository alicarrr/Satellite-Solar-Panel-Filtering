// In our CPP file we have to include the corresponding header file. To include our own files we use "" which tells the compiler to look for the file we are after here first before searching the system files.
#include "transforms/coordTransforms.hpp"
#include <math.h> 

// Now we can actually define what the functions in our header file does, think of this as defining what the code in those functions is, so that we can compile them if we need them.
Eigen::Vector2d polarToCartesian(Eigen::Vector2d r_polar)
{
    // Initialise as zeros.
    Eigen::Vector2d r_cartesian = Eigen::Vector2d::Zero();

    double r     = r_polar(0);
    double theta = r_polar(1);
    
    r_cartesian(0) = r * cos(theta);
    r_cartesian(1) = r * sin(theta);

    return r_cartesian;
}