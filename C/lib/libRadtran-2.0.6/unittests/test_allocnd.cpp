#include "catch.hpp"
#include "../libsrc_c/allocnd.h"

TEST_CASE ("apply_array_offsets_1D", "[allocnd, apply_array_offsets]") {

  SECTION ("Test 1D array offsets") {
    const size_t dimlen[1] = {4};
    int          data[4]   = {0, 1, 2, 3};
    int*         arr       = &data[0];

    {  // move left
      const int offsets[1] = {-2};

      int ierr = apply_array_offsets_int_1D (dimlen, offsets, &arr);
      REQUIRE (ierr == 0);

      REQUIRE (arr[2] == 0);
      REQUIRE (arr[3] == 1);
      REQUIRE (arr[4] == 2);
      REQUIRE (arr[5] == 3);
    }

    {  // move right
      const int offsets[1] = {4};

      int ierr = apply_array_offsets_int_1D (dimlen, offsets, &arr);
      REQUIRE (ierr == 0);

      REQUIRE (arr[-2] == 0);
      REQUIRE (arr[-1] == 1);
      REQUIRE (arr[0] == 2);
      REQUIRE (arr[1] == 3);
    }
  }
}

TEST_CASE ("apply_array_offsets_2D", "[allocnd, apply_array_offsets]") {
  SECTION ("Test 2D array offsets") {
    const size_t Nx = 2;
    const size_t Ny = 3;

    const size_t dimlen[2] = {Nx, Ny};
    int**        arr       = calloc_int_2D (Nx, Ny, "test_data");

    size_t cnt = 0;
    for (size_t i = 0; i < Nx; ++i) {
      for (size_t j = 0; j < Ny; ++j, ++cnt) {
        arr[i][j] = cnt;
      }
    }

    REQUIRE (arr[0][0] == 0);
    REQUIRE (arr[0][1] == 1);
    REQUIRE (arr[0][2] == 2);
    REQUIRE (arr[1][0] == 3);
    REQUIRE (arr[1][1] == 4);
    REQUIRE (arr[1][2] == 5);

    {  // move left
      const int offsets[2] = {-1, 0};

      int ierr = apply_array_offsets_int_2D (dimlen, offsets, &arr);
      REQUIRE (ierr == 0);

      REQUIRE (arr[0 - offsets[0]][0 - offsets[1]] == 0);
      REQUIRE (arr[0 - offsets[0]][1 - offsets[1]] == 1);
      REQUIRE (arr[0 - offsets[0]][2 - offsets[1]] == 2);
      REQUIRE (arr[1 - offsets[0]][0 - offsets[1]] == 3);
      REQUIRE (arr[1 - offsets[0]][1 - offsets[1]] == 4);
      REQUIRE (arr[1 - offsets[0]][2 - offsets[1]] == 5);

      // move back
      arr -= offsets[0];
      REQUIRE (arr[0][0] == 0);
      REQUIRE (arr[0][1] == 1);
      REQUIRE (arr[0][2] == 2);
      REQUIRE (arr[1][0] == 3);
      REQUIRE (arr[1][1] == 4);
      REQUIRE (arr[1][2] == 5);
    }

    free_int_2D (arr);
  }
}

TEST_CASE ("calloc_nd", "[allocnd, calloc_nd]") {
  SECTION ("Test 1D array allocation") {
    const size_t Nx = 2;

    int* arr = calloc_int_1D (Nx, "test_data");

    size_t cnt = 0;
    for (size_t i = 0; i < Nx; ++i, cnt++) {
      arr[i] = cnt;
    }

    REQUIRE (arr[0] == 0);
    REQUIRE (arr[1] == 1);

    free_int_1D (arr);
  }

  SECTION ("Test 2D array allocation") {
    const size_t Nx = 2;
    const size_t Ny = 3;

    int** arr = calloc_int_2D (Nx, Ny, "test_data");

    size_t cnt = 0;
    for (size_t i = 0; i < Nx; ++i) {
      for (size_t j = 0; j < Ny; ++j, ++cnt) {
        arr[i][j] = cnt;
      }
    }

    REQUIRE (arr[0][0] == 0);
    REQUIRE (arr[0][1] == 1);
    REQUIRE (arr[0][2] == 2);
    REQUIRE (arr[1][0] == 3);
    REQUIRE (arr[1][1] == 4);
    REQUIRE (arr[1][2] == 5);

    free_int_2D (arr);
  }

  SECTION ("Test 3D array allocation") {
    const size_t Nx = 3;
    const size_t Ny = 2;
    const size_t Nz = 5;

    double*** arr = calloc_double_3D (Nx, Ny, Nz, "test_data_3D");

    {  // set vals
      size_t cnt = 0;
      for (size_t i = 0; i < Nx; ++i) {
        for (size_t j = 0; j < Ny; ++j) {
          for (size_t k = 0; k < Nz; ++k, ++cnt) {
            arr[i][j][k] = cnt;
          }
        }
      }
    }

    {  //check vals again
      size_t cnt = 0;
      for (size_t i = 0; i < Nx; ++i)
        for (size_t j = 0; j < Ny; ++j)
          for (size_t k = 0; k < Nz; ++k, ++cnt)
            REQUIRE (arr[i][j][k] == cnt);
    }

    free_double_3D (arr);
  }

  SECTION ("Test 4D array allocation") {
    const size_t N1 = 3;
    const size_t N2 = 5;
    const size_t N3 = 2;
    const size_t N4 = 1;

    int**** arr = calloc_int_4D (N1, N2, N3, N4, "test_data_4D");

    {  // set vals
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = 0; l < N4; ++l, cnt++) {
              arr[i][j][k][l] = cnt;
            }
    }

    {  //check vals again
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = 0; l < N4; ++l, cnt++) {
              int trgt = static_cast<int> (cnt);
              REQUIRE (arr[i][j][k][l] == trgt);
            }
    }

    free_int_4D (arr);
  }

  SECTION ("Test 5D array allocation") {
    const size_t N1 = 3;
    const size_t N2 = 5;
    const size_t N3 = 7;
    const size_t N4 = 11;
    const size_t N5 = 13;

    size_t***** arr = calloc_size_t_5D (N1, N2, N3, N4, N5, "test_data_5D");

    {  // set vals
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = 0; l < N4; ++l)
              for (size_t m = 0; m < N5; ++m)
                arr[i][j][k][l][m] = ++cnt;
    }

    {  //check vals again
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = 0; l < N4; ++l)
              for (size_t m = 0; m < N5; ++m) {
                REQUIRE (arr[i][j][k][l][m] == ++cnt);
              }
    }

    free_size_t_5D (arr);
  }

  SECTION ("Test 6D array allocation") {
    const size_t N1 = 7;
    const size_t N2 = 3;
    const size_t N3 = 1;
    const size_t N4 = 7;
    const size_t N5 = 5;
    const size_t N6 = 3;

    float****** arr = calloc_float_6D (N1, N2, N3, N4, N5, N6, "test_data_6D");

    {  // set vals
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = 0; l < N4; ++l)
              for (size_t m = 0; m < N5; ++m)
                for (size_t n = 0; n < N6; ++n)
                  arr[i][j][k][l][m][n] = ++cnt;
    }

    {  //check vals again
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = 0; l < N4; ++l)
              for (size_t m = 0; m < N5; ++m)
                for (size_t n = 0; n < N6; ++n) {
                  float trgt = static_cast<float> (++cnt);
                  REQUIRE (arr[i][j][k][l][m][n] == trgt);
                }
    }

    free_float_6D (arr);
  }

  SECTION ("Test 7D array allocation") {
    const size_t N1 = 2;
    const size_t N2 = 3;
    const size_t N3 = 2;
    const size_t N4 = 3;
    const size_t N5 = 2;
    const size_t N6 = 3;
    const size_t N7 = 2;

    float******* arr = calloc_float_7D (N1, N2, N3, N4, N5, N6, N7, "test_data_7D");

    {  // set vals
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = 0; l < N4; ++l)
              for (size_t m = 0; m < N5; ++m)
                for (size_t n = 0; n < N6; ++n)
                  for (size_t o = 0; o < N7; ++o)
                    arr[i][j][k][l][m][n][o] = ++cnt;
    }

    {  //check vals again
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = 0; l < N4; ++l)
              for (size_t m = 0; m < N5; ++m)
                for (size_t n = 0; n < N6; ++n)
                  for (size_t o = 0; o < N7; ++o)
                    REQUIRE (arr[i][j][k][l][m][n][o] == ++cnt);
    }

    free_float_7D (arr);
  }

  SECTION ("Test 5D restricted array allocation") {
    const size_t N1   = 3;
    const size_t N2   = 5;
    const size_t N3   = 7;
    const size_t N4_s = 11;
    const size_t N4_e = 13;
    int          N5[N4_e + 1];

    for (size_t l = 0; l <= N4_e; ++l)
      N5[l] = -2147483648;  // default value to inf to produce hopefully an error if we ever use it

    for (size_t l = N4_s; l <= N4_e; ++l)
      N5[l] = l;

    int***** arr = calloc_int_5D_restricted (N1, N2, N3, N4_s, N4_e, N5, "test_data_5D_restricted");

    {  // set vals
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = N4_s; l <= N4_e; ++l)
              for (int m = 0; m < N5[l]; ++m) {
                REQUIRE (arr[i][j][k][l][m] == 0);
                arr[i][j][k][l][m] = ++cnt;
              }
    }

    {  //check vals again
      size_t cnt = 0;
      for (size_t i = 0; i < N1; ++i)
        for (size_t j = 0; j < N2; ++j)
          for (size_t k = 0; k < N3; ++k)
            for (size_t l = N4_s; l <= N4_e; ++l)
              for (int m = 0; m < N5[l]; ++m) {
                int trgt = static_cast<int> (++cnt);
                REQUIRE (arr[i][j][k][l][m] == trgt);
              }
    }

    free_int_5D_restricted (arr);
  }
}
