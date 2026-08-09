#include "catch.hpp"
#include "libsrc_c/common_math.h"

TEST_CASE ("Testing max_loc", "[common_math,max_loc]") {
  int ierr;

  SECTION ("Searching max_loc in increasing array") {
    double arr[4] = {0, 1, 2, 3};
    REQUIRE (max_loc (arr, 0, 3, &ierr) == 3);
    REQUIRE (ierr == 0);
    REQUIRE (max_loc (arr, 0, 2, &ierr) == 2);
    REQUIRE (ierr == 0);
  }

  SECTION ("Searching max_loc in decreasing array") {
    double arr[4] = {-0, -1, -2, -3};
    REQUIRE (max_loc (arr, 0, 3, &ierr) == 0);
    REQUIRE (ierr == 0);
    REQUIRE (max_loc (arr, 1, 3, &ierr) == 1);
    REQUIRE (ierr == 0);
  }

  SECTION ("Searching max_loc in ambiguous array") {
    double arr[4] = {0, 0, 1, 1};
    REQUIRE (max_loc (arr, 0, 3, &ierr) == 3);
    REQUIRE (ierr == -1);
  }

  SECTION ("Searching max_loc in single entry array") {
    double arr[1] = {1};
    REQUIRE (max_loc (arr, 0, 0, &ierr) == 0);
    REQUIRE (ierr == 0);
  }
}

TEST_CASE ("Testing sign function", "[common_math,sign]") {
  REQUIRE (sign (-1) == -1);
  REQUIRE (sign (1) == 1);
}

TEST_CASE ("Testing vec_AXPBY", "[common_math, vec_AXPBY]") {
  const size_t N    = 2;
  const double x[N] = {1, 3};
  const double y[N] = {5, 7};
  double       r[N];
  SECTION ("alpha=2, beta=0") {
    for (size_t i = 0; i < N; ++i)
      r[i] = y[i];
    const double alpha = 2;
    const double beta  = 0;
    REQUIRE (vec_AXPBY (N, alpha, x, beta, r) == 0);
    for (size_t i = 0; i < N; ++i)
      REQUIRE (r[i] == alpha * x[i] + beta * y[i]);
  }
  SECTION ("alpha=0, beta=1") {
    for (size_t i = 0; i < N; ++i)
      r[i] = y[i];
    const double alpha = 0;
    const double beta  = 1;
    REQUIRE (vec_AXPBY (N, alpha, x, beta, r) == 0);
    for (size_t i = 0; i < N; ++i)
      REQUIRE (r[i] == alpha * x[i] + beta * y[i]);
  }
  SECTION ("alpha=10, beta=10") {
    for (size_t i = 0; i < N; ++i)
      r[i] = y[i];
    const double alpha = 10;
    const double beta  = 10;
    REQUIRE (vec_AXPBY (N, alpha, x, beta, r) == 0);
    for (size_t i = 0; i < N; ++i)
      REQUIRE (r[i] == alpha * x[i] + beta * y[i]);
  }
}

TEST_CASE ("Testing vec_AXPBYPCZ", "[common_math, vec_AXPBYPCZ]") {
  const size_t N    = 2;
  const double x[N] = {1, 3};
  const double y[N] = {5, 7};
  const double z[N] = {13, 17};
  double       r[N];
  SECTION ("alpha=2, beta=1, gamma=0") {
    for (size_t i = 0; i < N; ++i)
      r[i] = z[i];
    const double alpha = 2;
    const double beta  = 1;
    const double gamma = 0;
    REQUIRE (vec_AXPBYPCZ (N, alpha, x, beta, y, gamma, r) == 0);
    for (size_t i = 0; i < N; ++i)
      REQUIRE (r[i] == alpha * x[i] + beta * y[i] + gamma * z[i]);
  }

  SECTION ("alpha=1, beta=1, gamma=0") {
    for (size_t i = 0; i < N; ++i)
      r[i] = z[i];
    const double alpha = 1;
    const double beta  = 1;
    const double gamma = 0;
    REQUIRE (vec_AXPBYPCZ (N, alpha, x, beta, y, gamma, r) == 0);
    for (size_t i = 0; i < N; ++i)
      REQUIRE (r[i] == alpha * x[i] + beta * y[i] + gamma * z[i]);
  }

  SECTION ("alpha=1, beta=1, gamma=1") {
    for (size_t i = 0; i < N; ++i)
      r[i] = z[i];
    const double alpha = 1;
    const double beta  = 1;
    const double gamma = 1;
    REQUIRE (vec_AXPBYPCZ (N, alpha, x, beta, y, gamma, r) == 0);
    for (size_t i = 0; i < N; ++i)
      REQUIRE (r[i] == alpha * x[i] + beta * y[i] + gamma * z[i]);
  }

  SECTION ("alpha=10, beta=10, gamma=10") {
    for (size_t i = 0; i < N; ++i)
      r[i] = z[i];
    const double alpha = 1;
    const double beta  = 1;
    const double gamma = 1;
    REQUIRE (vec_AXPBYPCZ (N, alpha, x, beta, y, gamma, r) == 0);
    for (size_t i = 0; i < N; ++i)
      REQUIRE (r[i] == alpha * x[i] + beta * y[i] + gamma * z[i]);
  }
}

TEST_CASE ("Testing reverse functions", "[common_math, reverse]") {
  SECTION ("float_reverse functions") {
    float x[2] = {1, 2};
    float y[3] = {1, 2, 3};

    SECTION ("test reverse of an even numbered array") {
      reverse_1d_f4 (2, x);
      REQUIRE (x[0] == 2);
      REQUIRE (x[1] == 1);
    }
    SECTION ("test reverse of an odd numbered array") {
      reverse_1d_f4 (3, y);
      REQUIRE (y[0] == 3);
      REQUIRE (y[1] == 2);
      REQUIRE (y[2] == 1);
    }
  }
  SECTION ("double reverse functions") {
    double x[2] = {1, 2};
    double y[3] = {1, 2, 3};

    SECTION ("test reverse of an even numbered array") {
      reverse_1d_f8 (2, x);
      REQUIRE (x[0] == 2);
      REQUIRE (x[1] == 1);
    }
    SECTION ("test reverse of an odd numbered array") {
      reverse_1d_f8 (3, y);
      REQUIRE (y[0] == 3);
      REQUIRE (y[1] == 2);
      REQUIRE (y[2] == 1);
    }
  }
}

TEST_CASE ("Testing MIN/MAX macros", "[common_math, MIN/MAX macros]") {
  SECTION ("Integer checks") {
    int x = 0;
    int y = 2;
    REQUIRE (MIN (x, y) == x);
    REQUIRE (MIN (y, x) == x);

    REQUIRE (MAX (x, y) == y);
    REQUIRE (MAX (y, x) == y);

    REQUIRE (MAX (x, x) == MIN (x, x));

    SECTION ("Check for side effects evaluation in MIN") {
      REQUIRE (MIN (++x, ++y) == 1);
      REQUIRE (x == 1);
      REQUIRE (y == 3);
    }
    SECTION ("Check for side effects evaluation in MAX") {
      REQUIRE (MAX (++x, ++y) == 3);
      REQUIRE (x == 1);
      REQUIRE (y == 3);
    }
  }
  SECTION ("Double checks") {
    double x = 0.5;
    double y = 2.5;
    REQUIRE (MIN (x, y) == x);
    REQUIRE (MIN (y, x) == x);

    REQUIRE (MAX (x, y) == y);
    REQUIRE (MAX (y, x) == y);

    REQUIRE (MAX (x, x) == MIN (x, x));

    SECTION ("Check for side effects evaluation in MIN") {
      REQUIRE (MIN (++x, ++y) == Approx (1.5));
      REQUIRE (x == Approx (1.5));
      REQUIRE (y == Approx (3.5));
    }
    SECTION ("Check for side effects evaluation in MAX") {
      REQUIRE (MAX (++x, ++y) == Approx (3.5));
      REQUIRE (x == Approx (1.5));
      REQUIRE (y == Approx (3.5));
    }
  }
  SECTION ("Mixed precision checks") {
    int   x = 1;
    float y = 1.1;
    REQUIRE (MIN (x, y) == x);
    REQUIRE (MIN (y, x) == x);

    REQUIRE (MAX (x, y) == y);
    REQUIRE (MAX (y, x) == y);

    REQUIRE (MAX (x, x) == MIN (x, x));

    SECTION ("Check for side effects evaluation in MIN") {
      REQUIRE (MIN (++x, ++y) == 2);
      REQUIRE (x == 2);
      REQUIRE (y == Approx (2.1));
    }
    SECTION ("Check for side effects evaluation in MAX") {
      REQUIRE (MAX (++x, ++y) == Approx (2.1));
      REQUIRE (x == 2);
      REQUIRE (y == Approx (2.1));
    }
  }
}

TEST_CASE ("Testing vec_count", "[common_math, vec_count]") {
  SECTION ("count_true") {
    const size_t N    = 3;
    const int    a[N] = {0, 1, 0};
    REQUIRE (vec_count_true (N, a) == 1);
    REQUIRE (vec_count_false (N, a) == 2);
  }
}

TEST_CASE ("Testing triangle_area", "[common_math, triangle_area]") {
  const double A[3] = {0, 0, 0};
  const double B[3] = {1, 0, 0};
  const double C[3] = {1, 1, 0};

  REQUIRE (vec_distance_between_two_points (2, A, B) == Approx (1));
  REQUIRE (vec_distance_between_two_points (2, B, C) == Approx (1));
  REQUIRE (vec_distance_between_two_points (2, A, C) == Approx (sqrt (2)));

  REQUIRE (triangle_area_by_vertices (A, B, C) == Approx (0.5));
  REQUIRE (triangle_area_by_vertices (A, C, B) == Approx (0.5));
  REQUIRE (triangle_area_by_vertices (B, A, C) == Approx (0.5));
}

TEST_CASE ("Testing is_inrange/is_between", "[common_math,is_inrange,is_between]") {
  REQUIRE (is_inrange (2, 1, 3) != 0);
  REQUIRE (is_inrange (2, 3, 1) != 0);
  REQUIRE (is_inrange (4, 1, 3) == 0);
  REQUIRE (is_inrange (-1, 1, 3) == 0);

  REQUIRE (is_between (2, 1, 3) != 0);
  REQUIRE (is_between (2, 3, 1) != 0);
  REQUIRE (is_between (4, 1, 3) == 0);
  REQUIRE (is_between (-1, 1, 3) == 0);

  // tests on the interval edges
  REQUIRE (is_inrange (1, 1, 3) != 0);
  REQUIRE (is_inrange (1, 3, 1) != 0);
  REQUIRE (is_inrange (3, 1, 3) != 0);
  REQUIRE (is_inrange (3, 3, 1) != 0);
  REQUIRE (is_between (1, 1, 3) == 0);
  REQUIRE (is_between (1, 3, 1) == 0);
  REQUIRE (is_between (3, 1, 3) == 0);
  REQUIRE (is_between (3, 3, 1) == 0);
}

TEST_CASE ("Testing cosd_sind", "[common_math,cosd_sind]") {
  for (int i = -360 * 200; i < 361 * 200; ++i) {
    const double phi = i / 100.;
    double       cosphi, sinphi;
    int          status = cosd_sind (phi, &cosphi, &sinphi);
    INFO ("phi is " << phi);
    INFO ("cos(phi) is " << cosphi << " should be " << cosd (phi));
    INFO ("sin(phi) is " << sinphi << " should be " << sind (phi));
    REQUIRE (status == 0);
    REQUIRE (cosphi == Approx (cosd (phi)).margin (std::numeric_limits<float>::epsilon()));
    REQUIRE (sinphi == Approx (sind (phi)).margin (std::numeric_limits<double>::epsilon()));
  }
}
