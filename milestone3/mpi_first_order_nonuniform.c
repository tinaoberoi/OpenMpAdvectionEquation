#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <time.h>
#define MASTER 0

int modulo_index(int idx, int N);
double set_guassian_matrix(int i, int j, double L, double dx, double u, double v);
double calc_first_order_value(int CHUNK_SIZE, int i, int j, double m1[CHUNK_SIZE + 2][CHUNK_SIZE + 2], double dt, double dx);
void first_order(int CHUNK_SIZE, int x, int y, int top, int bottom, int left, int right, 
                    double m1[CHUNK_SIZE + 2][CHUNK_SIZE + 2], double m2[CHUNK_SIZE + 2][CHUNK_SIZE + 2], 
                    double dt, double dx, double u, double v, int NT);
int numtasks, taskid, len, src;
char hostname[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;

int main(int argc, char *argv[])
{
    int N, NT;
    double L, T, u, v;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(hostname, &len);
    if(taskid == 0){
        N = atoi(argv[1]);
        NT = atoi(argv[2]);
        L = atof(argv[3]);
        T = atof(argv[2]);
        u = atof(argv[5]);
        v = atof(argv[6]);
        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&NT, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&L, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&u, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&v, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    } else {
        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&NT, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&L, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&u, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&v, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    double dx = L / N;
    double dt = T / NT;
    double stability_param = dx / sqrt(2 * (u * u + v * v));
    assert(dt <= stability_param);

    int r, c, top, bottom, left, right;
    
    int CHUNKSIZE = sqrt(N * N / numtasks);
    int N_CHUNK = N / CHUNKSIZE;
    r = taskid / N_CHUNK;
    c = taskid % N_CHUNK;

    top = (modulo_index(r - 1, N_CHUNK)) * N_CHUNK + c;
    bottom = (modulo_index(r + 1, N_CHUNK)) * N_CHUNK + c;
    left = r * N_CHUNK + (modulo_index(c - 1, N_CHUNK));
    right = r * N_CHUNK + (modulo_index(c + 1, N_CHUNK));
    double m1[CHUNKSIZE + 2][CHUNKSIZE + 2];
    double m2[CHUNKSIZE + 2][CHUNKSIZE + 2];
    int x = taskid / N_CHUNK;
    int y = taskid % N_CHUNK;

    double start = omp_get_wtime();
    #pragma omp parallel num_threads(4)
    for (int i = 0; i < CHUNKSIZE; i++)
    {
        for (int j = 0; j < CHUNKSIZE; j++)
        {
            m2[i+1][j+1] = set_guassian_matrix(i + x * CHUNKSIZE, j + y * CHUNKSIZE, L, dx, u, v);
        }
    }

    

    first_order(CHUNKSIZE, x, y, top, bottom, left, right, m1, m2, dt, dx, u, v, NT);

    if (taskid == 0)
    {
        double final_arr[N][N];

        //Update the array by master
        for (int m = 1; m < CHUNKSIZE+1; m++)
        {
            for (int n = 1; n < CHUNKSIZE+1; n++)
            {
                final_arr[m-1][n-1] = m2[m][n];
            }
        }

        // Update rest of the ranks
        for (int i = 0; i < N_CHUNK; i++)
        {
            for (int j = 0; j < N_CHUNK; j++)
            {
                if (i + j == 0)
                    continue;
                // Receive array from all ranks
                double temp[CHUNKSIZE+2][CHUNKSIZE+2];
                
                MPI_Recv(temp, (CHUNKSIZE+2) * (CHUNKSIZE+2), MPI_DOUBLE, i * N_CHUNK + j, 0, MPI_COMM_WORLD, &status);

                // Update data from ranks in final array
                for (int m = 1; m < CHUNKSIZE+1; m++)
                {
                    for (int n = 1; n < CHUNKSIZE+1; n++)
                    {
                        final_arr[(m-1) + i * CHUNKSIZE][(n-1) + j * CHUNKSIZE] = temp[m][n];
                    }
                }
            }
        }

        // for (int i = 0; i < N; i++)
        // {
        //     for (int j = 0; j < N; j++)
        //     {
        //         printf("%f      ", final_arr[i][j]);
        //     }
        //     printf("\n");
        // }
        double end = omp_get_wtime(); 
        printf("Time taken for first order non uniform: %lf \n", end-start);
    }
    else
    {
        MPI_Send(m2, (CHUNKSIZE+2) * (CHUNKSIZE+2), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

int modulo_index(int idx, int n)
{
    return (((idx % n) + n) % n);
}

double set_guassian_matrix(int i, int j, double L, double dx, double u, double v)
{
    double valx = pow(i * dx - L / 2, 2);
    double valy = pow(j * dx - L / 2, 2);
    double ran = (valx + valy) / (2 * pow(L / 4, 2));
    double temp = exp(-ran);
    return temp;
}

double calc_first_order_value(int CHUNK_SIZE, int i, int j, double m1[CHUNK_SIZE + 2][CHUNK_SIZE + 2], double dt, double dx)
{
    double output;
    double u = i*10e-9;
    double v = j*10e-9;
    if (u > 0 && v > 0)
    {
        output = m1[i][j] + dt *
                                  (-v * ((m1[i][j] - m1[i][j - 1]) / dx) - u * ((m1[i][j] - m1[i - 1][j]) / dx));
    }
    else if (u < 0 && v < 0)
    {
        output = m1[i][j] + dt *
                                  (-v * ((m1[i][j + 1] - m1[i][j]) / dx) - u * ((m1[i + 1][j] - m1[i][j]) / dx));
    }
    else if (u < 0)
    {
        output = m1[i][j] + dt *
                                  (-v * ((m1[i][j] - m1[i][j - 1]) / dx) - u * ((m1[i + 1][j] - m1[i][j]) / dx));
    }
    else
    {
        output = m1[i][j] + dt *
                                  (-v * ((m1[i][j + 1] - m1[i][j]) / dx) - u * ((m1[i][j] - m1[i - 1][j]) / dx));
    }
    return output;
}

void update(int CHUNK_SIZE, double m1[CHUNK_SIZE + 2][CHUNK_SIZE + 2], double m2[CHUNK_SIZE + 2][CHUNK_SIZE + 2],
            double top_recv_arr[CHUNK_SIZE], double bottom_recv_arr[CHUNK_SIZE], double left_recv_arr[CHUNK_SIZE],
            double right_recv_arr[CHUNK_SIZE], double dt, double dx, double u, double v, int taskid)
{
    // Pad Top;
    for(int j = 0; j < CHUNK_SIZE; j++){
        m2[0][j+1] = top_recv_arr[j];
        m1[0][j+1] = top_recv_arr[j];
    }

    // // Pad Bottom;
    for(int j = 0; j < CHUNK_SIZE; j++){
        m2[CHUNK_SIZE + 1][j+1] = bottom_recv_arr[j];
        m1[CHUNK_SIZE + 1][j+1] = bottom_recv_arr[j];
    }

    // Pad Left;
    for(int i = 0; i < CHUNK_SIZE; i++){
            // Need to flip columns before sending
            m2[i+1][0] = left_recv_arr[i];
            m1[i+1][0] = left_recv_arr[i];
    }

    // Pad Right;
    for(int i = 0; i < CHUNK_SIZE; i++){
        m2[i+1][CHUNK_SIZE + 1] = right_recv_arr[i];
        m1[i+1][CHUNK_SIZE + 1] = right_recv_arr[i];
    }

    // Update the inner matrix
    for (int i = 1; i < CHUNK_SIZE+1 ; i++)
    {
        for (int j = 1; j < CHUNK_SIZE+1 ; j++)
        {
            m2[i][j] = calc_first_order_value(CHUNK_SIZE, i, j, m1, dt, dx);
        }
    }
}

void first_order(int CHUNK_SIZE, int x, int y, int top, int bottom, int left, int right, 
                    double m1[CHUNK_SIZE + 2][CHUNK_SIZE + 2], double m2[CHUNK_SIZE + 2][CHUNK_SIZE + 2], 
                    double dt, double dx, double u, double v, int NT)
{
    double top_arr[CHUNK_SIZE], bottom_arr[CHUNK_SIZE], left_arr[CHUNK_SIZE], right_arr[CHUNK_SIZE];
    for (int t = 0; t <NT; t++)
    {
        #pragma omp parallel num_threads(4)
        for(int i = 0; i<CHUNK_SIZE + 2; i++){
            for(int j = 0; j<CHUNK_SIZE + 2; j++){
                m1[i][j] = m2[i][j];
            }
        }
        if ((x + y) % 2 == 0)
        {
            for (int l = 0; l < CHUNK_SIZE; l++)
            {
                left_arr[l] = m1[l+1][1];
                right_arr[l] = m1[l+1][CHUNK_SIZE];
                
            }

            for (int k = 0; k < CHUNK_SIZE; k++)
            {
                top_arr[k] = m1[1][k+1];
                bottom_arr[k] = m1[CHUNK_SIZE][k+1];
            }

            //Sending
            MPI_Send(top_arr, CHUNK_SIZE, MPI_DOUBLE, top, 0, MPI_COMM_WORLD);
            MPI_Send(bottom_arr, CHUNK_SIZE, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD);
            MPI_Send(left_arr, CHUNK_SIZE, MPI_DOUBLE, left, 0, MPI_COMM_WORLD);
            MPI_Send(right_arr, CHUNK_SIZE, MPI_DOUBLE, right, 0, MPI_COMM_WORLD);

            // Receiving
            MPI_Recv(bottom_arr, CHUNK_SIZE, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(top_arr, CHUNK_SIZE, MPI_DOUBLE, top, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(right_arr, CHUNK_SIZE, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(left_arr, CHUNK_SIZE, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);

            ////////////////////////////////////////////////////////////////////////////////////////////
            update(CHUNK_SIZE, m1, m2, top_arr, bottom_arr, left_arr, right_arr, dt, dx, u, v, taskid);
            ////////////////////////////////////////////////////////////////////////////////////////////
        }
        else if ((x + y) % 2 != 0)
        {
            // Receiving
            MPI_Recv(bottom_arr, CHUNK_SIZE, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(top_arr, CHUNK_SIZE, MPI_DOUBLE, top, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(right_arr, CHUNK_SIZE, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(left_arr, CHUNK_SIZE, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);

            ////////////////////////////////////////////////////////////////////////////////////////////
            update(CHUNK_SIZE, m1, m2, top_arr, bottom_arr, left_arr, right_arr, dt, dx, u, v, taskid);
            ////////////////////////////////////////////////////////////////////////////////////////////

            for (int l = 0; l < CHUNK_SIZE; l++)
            {   
                left_arr[l] = m1[l+1][1];
                right_arr[l] = m1[l+1][CHUNK_SIZE];   
            }

            for (int k = 0; k < CHUNK_SIZE; k++)
            {
                top_arr[k] = m1[1][k+1];
                bottom_arr[k] = m1[CHUNK_SIZE][k+1];
            }
                

            // Sending
            MPI_Send(top_arr, CHUNK_SIZE, MPI_DOUBLE, top, 0, MPI_COMM_WORLD);
            MPI_Send(bottom_arr, CHUNK_SIZE, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD);
            MPI_Send(left_arr, CHUNK_SIZE, MPI_DOUBLE, left, 0, MPI_COMM_WORLD);
            MPI_Send(right_arr, CHUNK_SIZE, MPI_DOUBLE, right, 0, MPI_COMM_WORLD);
        }
    }
}