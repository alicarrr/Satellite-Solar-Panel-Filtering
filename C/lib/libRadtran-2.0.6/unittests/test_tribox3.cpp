#include "catch.hpp"
#include "libsrc_c/common_math.h"
#include "libsrc_c/tribox3.h"

TEST_CASE ("Testing bounding box intersection", "[AABB|Triangle intersection]") {
  double boxcenter[3]   = {1, 1, 1};
  double boxhalfsize[3] = {1, 1, 1};

  SECTION ("Large flat triangle") {
    double A[3] = {-10, -10, 1};
    double B[3] = {10, -10, 1};
    double C[3] = {5, 10, 1};

    REQUIRE (triBoxOverlap (boxcenter, boxhalfsize, A, B, C) != 0);
  }

  SECTION ("Small flat triangle") {
    double A[3] = {1, 1, 1};
    double B[3] = {1.1, 1, 1};
    double C[3] = {1.05, 1.1, 1};

    REQUIRE (triBoxOverlap (boxcenter, boxhalfsize, A, B, C) != 0);
  }

  SECTION ("Outside on the right side") {
    double A[3] = {2.00001, 1, 1};
    double B[3] = {2.1, 1, 1};
    double C[3] = {2.00001, 1, 1.1};

    REQUIRE (triBoxOverlap (boxcenter, boxhalfsize, A, B, C) == 0);
  }

  SECTION ("Outside on the front side") {
    double A[3] = {1, 2.00001, 1};
    double B[3] = {1, 2.1, 1};
    double C[3] = {1, 2.00001, 1.1};

    REQUIRE (triBoxOverlap (boxcenter, boxhalfsize, A, B, C) == 0);
  }

  SECTION ("Outside on the top side") {
    double A[3] = {1, 1, 2.00001};
    double B[3] = {1, 1, 2.1};
    double C[3] = {1, 1.1, 2.00001};

    REQUIRE (triBoxOverlap (boxcenter, boxhalfsize, A, B, C) == 0);
  }
}
