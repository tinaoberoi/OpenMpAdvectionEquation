#pragma once
#include <stdio.h>


double calc_c_n1_first_order_upwind(double **C_n1, int i, int j, double dt, double dx, double u, double v, int N);

void first_order_upwind(double **C_n1, double **C_n2, double dt, double dx, double u, double v, int N, int NT, int nti, int ntj);

double calc_c_n1_second_order_upwind(double **C_n1, int i, int j, double dt, double dx, double u, double v, int N);

void second_order_upwind(double **C_n1, double **C_n2, double dt, double dx, double u, double v, int N, int NT, int nti, int ntj);

void lax_method(double **C_n1, double **C_n2, double dt, double dx, double u, double v, int N, int NT, int nti, int ntj);

void populate_guassian(double **C_n1, double L, double dt, double dx, double u, double v, int N, char* algo, int nti, int ntj);

void write_to_file(FILE *fptr, double **C_n1, char *filename, int N);

int fix_boundary_index(int idx, int N);
