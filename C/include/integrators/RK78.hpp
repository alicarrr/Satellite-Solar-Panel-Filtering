#pragma once

/* Derivative function type */
typedef void (*DerivFunc)(double t, double* x, double* x_dot, void* params);

int rkf78_step(
    double t,
    double h,
    double* x,
    int n,
    void* params,
    DerivFunc f,
    double* x_next,
    double* err_est
);
