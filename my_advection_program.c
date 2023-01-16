#include "my_advection_program.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    int N, NT, nti, ntj;
    double L, T, u, v;
    char* algo;

    N = atoi(argv[1]);
    printf("N : %d\n", N);
    NT = atoi(argv[2]);
    printf("NT : %d\n", NT);
    L = atof(argv[3]);
    printf("L : %f\n", L);
    T = atof(argv[4]);
    printf("T : %f\n", T);
    u = atof(argv[5]);
    printf("u : %.10f\n", u);
    v = atof(argv[6]);
    printf("v : %.10f\n", v);

    algo = argv[7];
    printf("Algorithm : %s\n", algo);

    nti = atoi(argv[8]);
    printf("Number of threads for outer loops: %d\n", nti);

    ntj = 1;

    double **C_n1 = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++)
    {
        C_n1[i] = (double *)malloc(N * sizeof(double));
    }
    double **C_n2 = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++)
    {
        C_n2[i] = (double *)malloc(N * sizeof(double));
    }
    double dx = L / N;
    double dt = T / NT;
    double stability_param = dx / sqrt(2 * (u * u + v * v));
    assert(dt <= stability_param);

    populate_guassian(C_n1, L, dt, dx, u, v, N, algo, nti, ntj);

    /*Choose method to call*/
    if (strcmp(algo, "lax") == 0) 
    {
        printf("Executing Lax Method.\n");
        lax_method(C_n1, C_n2, dt, dx, u, v, N, NT, nti, ntj);
    }
    else if (strcmp(algo, "first_order") == 0) 
    {
        printf("Executing First Order Method.\n");
        first_order_upwind(C_n1, C_n2, dt, dx, u, v, N, NT, nti, ntj);
    }
    else {
        printf("Executing Second Order Method.\n");
        second_order_upwind(C_n1, C_n2, dt, dx, u, v, N, NT, nti, ntj);
    }

    return 0;
}

void write_to_file(FILE *fptr, double **C_n1, char *filename, int N)
{
    fptr = fopen(filename, "w");
    if (fptr == NULL)
    {
        printf("Error!");
        exit(1);
    }
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
            fprintf(fptr, "%e, ", C_n1[x][y]);
        }
        fprintf(fptr, "\n");
    }
}

int fix_boundary_index(int idx, int N)
{
    idx = (((idx % N) + N) % N);
    return idx;
}

void populate_guassian(double **C_n1, double L, double dt, double dx, double u, double v, int N, char* algo, int nti, int ntj)
{
    int i, j;
    double start = omp_get_wtime();
    #ifdef isParallel
    #pragma omp parallel for default(none) shared(C_n1, L, dt, dx, u, v, N) private(i, j) num_threads(8) schedule(static)
    #endif
    for (i = 0; i < N; i++)
    {
        #ifdef isParallel
        #pragma omp parallel for default(none) shared(C_n1, L, dt, dx, u, v, N, i) private(j) num_threads(8) schedule(static)
        #endif
        for (j = 0; j < N; j++)
        {
            double valx = pow(i * dx - L / 2, 2);
            double valy = pow(j * dx - L / 2, 2);
            double ran = (valx + valy) / (2 * pow(L / 4, 2));
            double temp = exp(-ran);
            C_n1[i][j] = temp;
        }
    }
    double end = omp_get_wtime(); 
    printf("Time taken for computing guassian for algo {%s} : %lf \n", algo, end-start);
    #ifdef isParallel
    FILE *fptr;
    char* guassian_file = (char*)malloc(1 + strlen(algo) + strlen("_guassian.txt"));
    strcpy(guassian_file, algo);
    strcat(guassian_file, "/serial_guassian.txt");
    printf("Printing to : %s\n", guassian_file);
    write_to_file(fptr, C_n1, guassian_file, N);
    #endif
}

void lax_method(double **C_n1, double **C_n2, double dt, double dx, double u, double v, int N, int NT, int nti, int ntj)
{
    int i, j, n;
    
    double start = omp_get_wtime(); 
    for (n = 0; n < NT; n++)
    {
        #ifdef isParallel
        #pragma omp parallel for default(none) shared(C_n1, C_n2, dt, dx, u, v, N) private(i, j) num_threads(8) schedule(static)
        #endif
        for (i = 0; i < N; i++)
        {
            #ifdef isParallel
            #pragma omp parallel for default(none) shared(C_n1, C_n2, dt, dx, u, v, N, i) private(j) num_threads(8) schedule(static)
            #endif
            for (j = 0; j < N; j++)
            {
                double up = C_n1[(((i - 1) % N) + N) % N][j];
                double down = C_n1[(((i + 1) % N) + N) % N][j];
                double left = C_n1[i][(((j - 1) % N) + N) % N];
                double right = C_n1[i][((j + 1 % N) + N) % N];
                C_n2[i][j] = (up + down + left + right) / 4;
                C_n2[i][j] -= (dt / (2 * dx)) * (u * (down - up) + v * (right - left));
            }
        }

        double** helper = C_n2;
        C_n2 = C_n1;
        C_n1 = helper;

        #ifdef isParallel
        FILE *fptr;
        if (n == 400)
        {
            char str_n[20];
            snprintf(str_n, sizeof(str_n), "%d", n);
            char* filename = (char*)malloc(1 + strlen("/lax/") + strlen(str_n) + strlen(".txt"));
            strcpy(filename, "lax/");
            strcat(filename, str_n);
            strcat(filename, "_serial.txt");
            printf("Printing to : %s\n", filename);
            write_to_file(fptr, C_n1, filename, N);
        }
        #endif
    }
    double end = omp_get_wtime(); 
    printf("Time Taken: %lf \n", end-start);
}

double calc_c_n1_first_order_upwind(double **C_n1, int i, int j, double dt, double dx, double u, double v, int N)
{
    double output;
    if (u > 0 && v > 0)
    {
        output = C_n1[i][j] + dt *
                                  (-v * ((C_n1[i][j] - C_n1[i][fix_boundary_index(j - 1, N)]) / dx) - u * ((C_n1[i][j] - C_n1[fix_boundary_index(i - 1, N)][j]) / dx));
    }
    else if (u < 0 && v < 0)
    {
        output = C_n1[i][j] + dt *
                                  (-v * ((C_n1[i][fix_boundary_index(j + 1, N)] - C_n1[i][j]) / dx) - u * ((C_n1[fix_boundary_index(i + 1, N)][j] - C_n1[i][j]) / dx));
    }
    else if (u < 0)
    {
        output = C_n1[i][j] + dt *
                                  (-v * ((C_n1[i][j] - C_n1[i][fix_boundary_index(j - 1, N)]) / dx) - u * ((C_n1[fix_boundary_index(i + 1, N)][j] - C_n1[i][j]) / dx));
    }
    else
    {
        output = C_n1[i][j] + dt *
                                  (-v * ((C_n1[i][fix_boundary_index(j + 1, N)] - C_n1[i][j]) / dx) - u * ((C_n1[i][j] - C_n1[fix_boundary_index(i - 1, N)][j]) / dx));
    }
    return output;
}

void first_order_upwind(double **C_n1, double **C_n2, double dt, double dx, double u, double v, int N, int NT, int nti, int ntj)
{
    int i, j, n;

    double start = omp_get_wtime(); 
    for (n = 0; n < NT; n++)
    {
        #ifdef isParallel
        #pragma omp parallel for default(none) shared(C_n1, C_n2, dt, dx, u, v, N) private(i, j) num_threads(8) schedule(static)
        #endif
        for (i = 0; i < N; i++)
        {
            #ifdef isParallel
            #pragma omp parallel for default(none) shared(C_n1, C_n2, dt, dx, u, v, N, i) private(j) num_threads(8) schedule(static)
            #endif
            for (j = 0; j < N; j++)
            {
                C_n2[i][j] = calc_c_n1_first_order_upwind(C_n1, i, j, dt, dx, u, v, N);
            }
        }

        double** helper = C_n2;
        C_n2 = C_n1;
        C_n1 = helper;

        #ifdef isParallel
        FILE *fptr;
        if (n == 400)
        {
            char str_n[20];
            snprintf(str_n, sizeof(str_n), "%d", n);
            char* filename = (char*)malloc(1 + strlen("/first_order/") + strlen(str_n) + strlen(".txt"));
            strcpy(filename, "first_order/");
            strcat(filename, str_n);
            strcat(filename, "_serial.txt");
            printf("Printing to : %s\n", filename);
            write_to_file(fptr, C_n1, filename, N);
        }
        #endif
    }
    double end = omp_get_wtime(); 
    printf("Time Taken: %lf \n", end-start);
}

double calc_c_n1_second_order_upwind(double **C_n1, int i, int j, double dt, double dx, double u, double v, int N)
{
    double output;
    if (u > 0 && v > 0)
    {
        output = C_n1[i][j] + dt *
                                  (-v * ((3 * C_n1[i][j] - 4 * C_n1[i][(((j - 1) % N) + N) % N] + C_n1[i][(((j - 2) % N) + N) % N]) / (2 * dx)) - u * ((3 * C_n1[i][j] - 4 * C_n1[(((i - 1) % N) + N) % N][j] + C_n1[(((i - 2) % N) + N) % N][j]) / (2 * dx)));
    }
    else if (u < 0 && v < 0)
    {
        output = C_n1[i][j] + dt *
                                  (-v * ((-1 * C_n1[i][(((j + 2) % N) + N) % N] + 4 * C_n1[i][(((j + 1) % N) + N) % N] - 3 * C_n1[i][j]) / (2 * dx)) - u * ((-1 * C_n1[(((i + 2) % N) + N) % N][j] + 4 * C_n1[(((i + 1) % N) + N) % N][j] - 3 * C_n1[i][j]) / (2 * dx)));
    }
    else if (u < 0)
    {
        output = C_n1[i][j] + dt *
                                  (-v * ((3 * C_n1[i][j] - 4 * C_n1[i][(((j - 1) % N) + N) % N] + C_n1[i][(((j - 2) % N) + N) % N]) / (2 * dx)) - u * ((-1 * C_n1[(((i + 2) % N) + N) % N][j] + 4 * C_n1[(((i + 1) % N) + N) % N][j] - 3 * C_n1[i][j]) / (2 * dx)));
    }
    else
    {
        output = C_n1[i][j] + dt *
                                  (-v * ((-1 * C_n1[i][(((j + 2) % N) + N) % N] + 4 * C_n1[i][(((j + 1) % N) + N) % N] - 3 * C_n1[i][j]) / (2 * dx)) - u * ((3 * C_n1[i][j] - 4 * C_n1[(((i - 1) % N) + N) % N][j] + C_n1[(((i - 2) % N) + N) % N][j]) / (2 * dx)));
    }
    return output;
}

void second_order_upwind(double **C_n1, double **C_n2, double dt, double dx, double u, double v, int N, int NT, int nti, int ntj)
{
    int i, j, n;

    double start = omp_get_wtime(); 
    for (n = 0; n < NT; n++)
    {
        #ifdef isParallel
        #pragma omp parallel for default(none) shared(C_n1, C_n2, dt, dx, u, v, N, NT) private(i, j) num_threads(8) schedule(static)
        #endif
        for (i = 0; i < N; i++)
        {
            #ifdef isParallel
            #pragma omp parallel for default(none) shared(C_n1, C_n2, dt, dx, u, v, N, i) private(j) num_threads(8) schedule(static)
            #endif
            for (j = 0; j < N; j++)
            {
                C_n2[i][j] = calc_c_n1_second_order_upwind(C_n1, i, j, dt, dx, u, v, N);
            }
        }

        double** helper = C_n2;
        C_n2 = C_n1;
        C_n1 = helper;

        #ifdef isParallel
        FILE *fptr;
        if (n == 400)
        {
            char str_n[20];
            snprintf(str_n, sizeof(str_n), "%d", n);
            char* filename = (char*)malloc(1 + strlen("/first_order/") + strlen(str_n) + strlen(".txt"));
            strcpy(filename, "second_order/");
            strcat(filename, str_n);
            strcat(filename, "_serial.txt");
            printf("Printing to : %s\n", filename);
            write_to_file(fptr, C_n1, filename, N);
        }
        #endif
    }
    double end = omp_get_wtime(); 
    printf("Time Taken: %lf \n", end-start);
}