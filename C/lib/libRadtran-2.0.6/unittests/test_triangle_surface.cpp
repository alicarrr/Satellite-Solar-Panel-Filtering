#include "catch.hpp"
#include "libsrc_c/common_math.h"
#include "libsrc_c/triangle_surface.h"

TEST_CASE ("Testing simple intersection", "[woop intersection]") {
  double A[3] = {0, 0, 0};
  double B[3] = {1, 0, 0};
  double C[3] = {0, 1, 0};

  double tnear = 0;
  double tfar  = std::numeric_limits<double>::infinity();

  double distance;
  double hit[3];
  int    lhit;

  SECTION ("Hitting triangle from above, straight down") {
    double origin[3]    = {0.5, 0.5, 1};
    double direction[3] = {0, 0, -1};
    lhit                = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);

    REQUIRE (lhit != 0);
    REQUIRE (distance == 1);
    double target[3] = {0.5, 0.5, 0};
    for (size_t i = 0; i < 3; ++i)
      REQUIRE (hit[i] == Approx (target[i]).epsilon (1e-6));
  }

  SECTION ("Hitting triangle from above, straight down, not normalized direction") {
    double origin[3]    = {0.5, 0.5, 1};
    double direction[3] = {0, 0, -100};
    lhit                = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);

    REQUIRE (lhit != 0);
    REQUIRE (distance == 1. / vec3_norm2 (direction));
    double target[3] = {0.5, 0.5, 0};
    for (size_t i = 0; i < 3; ++i)
      REQUIRE (hit[i] == Approx (target[i]).epsilon (1e-6));
  }

  SECTION ("Hitting triangle from below, straight up") {
    double origin[3]    = {0.5, 0.5, -1};
    double direction[3] = {0, 0, +1};
    lhit                = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);

    REQUIRE (lhit != 0);
    REQUIRE (distance == 1);
    double target[3] = {0.5, 0.5, 0};
    for (size_t i = 0; i < 3; ++i)
      REQUIRE (hit[i] == Approx (target[i]).epsilon (1e-6));
  }

  SECTION ("Miss triangle flying away from it towards top") {
    double origin[3]    = {0.5, 0.5, 1e-30};
    double direction[3] = {0, 0, +1};
    lhit                = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);
    REQUIRE (lhit == 0);
  }

  SECTION ("Miss triangle flying away from it towards bot") {
    double origin[3]    = {0.5, 0.5, -1e-30};
    double direction[3] = {0, 0, -1};
    lhit                = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);
    REQUIRE (lhit == 0);
  }

  SECTION ("Hit triangle tiny distance") {
    double origin[3]    = {0.5, 0.5, -1e-10};
    double direction[3] = {0, 0, 1};
    lhit                = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);
    REQUIRE (lhit != 0);
    REQUIRE (distance == 1e-10);
  }

  SECTION ("Hit triangle tilted angle") {
    double origin[3]    = {0., 0.5, .1};
    double direction[3] = {1, 0, -1};
    vec_scale (3, 1. / vec3_norm2 (direction), direction);  // normalize direction
    lhit = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);
    REQUIRE (lhit != 0);
    REQUIRE (distance == Approx (origin[2] * sqrt (2)));
  }
}

TEST_CASE ("Testing large triangle intersections", "[woop intersection, large_triangle]") {
  double dx   = 1e100;
  double A[3] = {-dx, -dx, 0};
  double B[3] = {dx, -dx, 0};
  double C[3] = {0, dx, 0};

  double tnear = 0;
  double tfar  = std::numeric_limits<double>::infinity();

  double distance;
  double hit[3];
  int    lhit;

  SECTION ("Hitting triangle from above, straight down") {
    double origin[3]    = {0.5, 0, 1};
    double direction[3] = {0, 0, -1};
    lhit                = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);

    REQUIRE (lhit != 0);
    REQUIRE (distance == 1);
    REQUIRE (hit[0] == origin[0]);
    REQUIRE (hit[1] == origin[1]);
    REQUIRE (hit[2] == 0.0);
  }
}

TEST_CASE ("Testing simple vertical intersection", "[woop intersection]") {
  double A[3] = {0, 0, 0};
  double B[3] = {1, 0, 0};
  double C[3] = {0, 0, 1};

  double tnear = 0;
  double tfar  = std::numeric_limits<double>::infinity();

  double distance;
  double hit[3];
  int    lhit;

  SECTION ("Hitting triangle from front") {
    double origin[3]    = {0.5, -1, 0.5};
    double direction[3] = {0, 1, 0};
    lhit                = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);

    REQUIRE (lhit != 0);
    REQUIRE (distance == 1);
    double target[3] = {0.5, 0, 0.5};
    for (size_t i = 0; i < 3; ++i)
      REQUIRE (hit[i] == Approx (target[i]).epsilon (1e-6));
  }

  SECTION ("Hitting triangle from back") {
    double origin[3]    = {0.5, 1, 0.5};
    double direction[3] = {0, -1, 0};
    lhit                = triangle_intersection_woop (origin, direction, A, B, C, tnear, tfar, &distance, hit);

    REQUIRE (lhit != 0);
    REQUIRE (distance == 1);
    double target[3] = {0.5, 0, 0.5};
    for (size_t i = 0; i < 3; ++i)
      REQUIRE (hit[i] == Approx (target[i]).epsilon (1e-6));
  }
}
