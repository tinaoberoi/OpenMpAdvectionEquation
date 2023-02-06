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
double calc_lax_value(double up, double down, double left, double right, double u, double v, double dt, double dx);
void update(int CHUNK_SIZE, double m1[CHUNK_SIZE][CHUNK_SIZE], double m2[CHUNK_SIZE][CHUNK_SIZE],
            double top_recv_arr[CHUNK_SIZE], double bottom_recv_arr[CHUNK_SIZE], double left_recv_arr[CHUNK_SIZE],
            double right_recv_arr[CHUNK_SIZE], double dt, double dx, double u, double v);

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
        T = atof(argv[4]);
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
    double top_arr[CHUNKSIZE], bottom_arr[CHUNKSIZE], left_arr[CHUNKSIZE], right_arr[CHUNKSIZE];
    double top_recv_arr[CHUNKSIZE], bottom_recv_arr[CHUNKSIZE], left_recv_arr[CHUNKSIZE], right_recv_arr[CHUNKSIZE];
    r = taskid / N_CHUNK;
    c = taskid % N_CHUNK;

    top = (modulo_index(r - 1, N_CHUNK)) * N_CHUNK + c;
    bottom = (modulo_index(r + 1, N_CHUNK)) * N_CHUNK + c;
    left = r * N_CHUNK + (modulo_index(c - 1, N_CHUNK));
    right = r * N_CHUNK + (modulo_index(c + 1, N_CHUNK));
    double m1[CHUNKSIZE][CHUNKSIZE];
    double m2[CHUNKSIZE][CHUNKSIZE];
    int x = taskid / N_CHUNK;
    int y = taskid % N_CHUNK;

    for (int i = 0; i < CHUNKSIZE; i++)
    {
        for (int j = 0; j < CHUNKSIZE; j++)
        {
            m1[i][j] = set_guassian_matrix(i + x * CHUNKSIZE, j + y * CHUNKSIZE, L, dx, u, v);
        }
    }

    for (int t = 0; t < NT; t++)
    {
        if ((x + y) % 2 == 0)
        {
            for (int k = 0; k < CHUNKSIZE; k++)
            {
                top_arr[k] = m1[0][k];
                bottom_arr[k] = m1[CHUNKSIZE - 1][k];
                left_arr[k] = m1[k][0];
                right_arr[k] = m1[k][CHUNKSIZE - 1];
            }

            // Sending
            MPI_Send(top_arr, CHUNKSIZE, MPI_DOUBLE, top, 0, MPI_COMM_WORLD);
            MPI_Send(bottom_arr, CHUNKSIZE, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD);
            MPI_Send(left_arr, CHUNKSIZE, MPI_DOUBLE, left, 0, MPI_COMM_WORLD);
            MPI_Send(right_arr, CHUNKSIZE, MPI_DOUBLE, right, 0, MPI_COMM_WORLD);

            ////////////////////////////////////////////////////////////////////////////////////////////////////

            // Receiving
            MPI_Recv(bottom_arr, CHUNKSIZE, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(top_arr, CHUNKSIZE, MPI_DOUBLE, top, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(right_arr, CHUNKSIZE, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(left_arr, CHUNKSIZE, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);

            ////////////////////////////////////////////////////////////////////////////////////
            update(CHUNKSIZE, m1, m2, top_arr, bottom_arr, left_arr, right_arr, dt, dx, u, v);
            ////////////////////////////////////////////////////////////////////////////////////////////////////
        }
        else if ((x + y) % 2 != 0)
        {
            // Receiving
            MPI_Recv(bottom_arr, CHUNKSIZE, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(top_arr, CHUNKSIZE, MPI_DOUBLE, top, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(right_arr, CHUNKSIZE, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(left_arr, CHUNKSIZE, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);

            ////////////////////////////////////////////////////////////////////////////////////
            update(CHUNKSIZE, m1, m2, top_arr, bottom_arr, left_arr, right_arr, dt, dx, u, v);
            ////////////////////////////////////////////////////////////////////////////////////

            for (int k = 0; k < CHUNKSIZE; k++)
            {
                top_arr[k] = m1[0][k];
                bottom_arr[k] = m1[CHUNKSIZE - 1][k];
                left_arr[k] = m1[k][0];
                right_arr[k] = m1[k][CHUNKSIZE - 1];
            }

            // Sending
            MPI_Send(top_arr, CHUNKSIZE, MPI_DOUBLE, top, 0, MPI_COMM_WORLD);
            MPI_Send(bottom_arr, CHUNKSIZE, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD);
            MPI_Send(left_arr, CHUNKSIZE, MPI_DOUBLE, left, 0, MPI_COMM_WORLD);
            MPI_Send(right_arr, CHUNKSIZE, MPI_DOUBLE, right, 0, MPI_COMM_WORLD);
        }

        if (taskid == 0)
        {
            double final_arr[N][N];

            // Update the array by master
            for (int i = 0; i < N_CHUNK; i++)
            {
                for (int j = 0; j < N_CHUNK; j++)
                {
                    final_arr[i][j] = m2[i][j];
                }
            }

            for (int i = 0; i < N_CHUNK; i++)
            {
                for (int j = 0; j < N_CHUNK; j++)
                {
                    if (i + j == 0)
                        continue;
                    // Receive array from all ranks
                    double temp[CHUNKSIZE][CHUNKSIZE];
                    MPI_Recv(temp, CHUNKSIZE * CHUNKSIZE, MPI_DOUBLE, i * N_CHUNK + j, 0, MPI_COMM_WORLD, &status);

                    // Update data from ranks in final array
                    for (int m = 0; m < CHUNKSIZE; m++)
                    {
                        for (int n = 0; n < CHUNKSIZE; n++)
                        {
                            final_arr[m + i * CHUNKSIZE][n + j * CHUNKSIZE] = temp[m][n];
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
        }
        else
        {
            MPI_Send(m2, CHUNKSIZE * CHUNKSIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
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

double calc_lax_value(double up, double down, double left, double right, double u, double v, double dt, double dx)
{
    double lax_val = (up + down + left + right) / 4 - ((dt / (2 * dx)) * (u * (down - up) + v * (right - left)));
    return lax_val;
}

void update(int CHUNK_SIZE, double m1[CHUNK_SIZE][CHUNK_SIZE], double m2[CHUNK_SIZE][CHUNK_SIZE],
            double top_recv_arr[CHUNK_SIZE], double bottom_recv_arr[CHUNK_SIZE], double left_recv_arr[CHUNK_SIZE],
            double right_recv_arr[CHUNK_SIZE], double dt, double dx, double u, double v)
{

    double up, down, left, right;
    // Update the inner matrix
    for (int i = 1; i < CHUNK_SIZE - 1; i++)
    {
        for (int j = 1; j < CHUNK_SIZE - 1; j++)
        {
            up = top_recv_arr[j];
            down = m1[i - 1][j];
            left = m1[i][j - 1];
            right = m1[i][j + 1];
            m2[i][j] = calc_lax_value(up, down, left, right, u, v, dt, dx);
        }
    }

    // Updating boundaries without corners
    for (int j = 1; j < CHUNK_SIZE - 1; j++)
    {
        up = top_recv_arr[j];
        down = m1[1][j];
        left = m1[0][j - 1];
        right = m1[0][j + 1];
        m2[0][j] = calc_lax_value(up, down, left, right, u, v, dt, dx);

        up = m1[CHUNK_SIZE - 2][j];
        down = bottom_recv_arr[j];
        left = m1[CHUNK_SIZE - 1][j - 1];
        right = m1[CHUNK_SIZE - 1][j + 1];
        m2[CHUNK_SIZE - 1][j] = calc_lax_value(up, down, left, right, u, v, dt, dx);
    }

    for (int i = 1; i < CHUNK_SIZE - 1; i++)
    {
        up = m1[i - 1][0];
        down = m1[i + 1][0];
        left = left_recv_arr[i];
        right = m1[i][1];
        m2[i][0] = calc_lax_value(up, down, left, right, u, v, dt, dx);

        up = m1[i - 1][CHUNK_SIZE - 1];
        down = m1[i + 1][CHUNK_SIZE - 1];
        left = m1[i][CHUNK_SIZE - 2];
        right = right_recv_arr[i];
        m2[i][CHUNK_SIZE - 1] = calc_lax_value(up, down, left, right, u, v, dt, dx);
    }

    // Update corners
    up = top_recv_arr[0];
    down = CHUNK_SIZE == 1 ? bottom_recv_arr[0] : m1[1][0];
    left = left_recv_arr[0];
    right = CHUNK_SIZE == 1 ? right_recv_arr[0] : m1[0][1];
    m2[0][0] = calc_lax_value(up, down, left, right, u, v, dt, dx);
    // printf("TASK: %d :: top = %f bottom = %f left = %f right = %f \n", taskid, top_recv_arr[0], m1[1][0], left_recv_arr[0], m1[0][1]);

    up = CHUNK_SIZE == 1 ? top_recv_arr[0] : m1[CHUNK_SIZE - 2][0];
    down = bottom_recv_arr[0];
    left = left_recv_arr[CHUNK_SIZE - 1];
    right = CHUNK_SIZE == 1 ? right_recv_arr[0] : m1[CHUNK_SIZE - 1][1];
    m2[CHUNK_SIZE - 1][0] = calc_lax_value(up, down, left, right, u, v, dt, dx);

    // printf("TASK: %d :: t = %f b = %f l = %f r = %f \n", taskid, m1[CHUNKSIZE - 2][0], bottom_recv_arr[0], left_recv_arr[CHUNKSIZE - 1], m1[CHUNKSIZE - 1][1]);
    up = top_recv_arr[CHUNK_SIZE - 1];
    down = CHUNK_SIZE == 1 ? bottom_recv_arr[0] : m1[1][CHUNK_SIZE - 1];
    left = CHUNK_SIZE == 1 ? left_recv_arr[0] : m1[0][CHUNK_SIZE - 2];
    right = right_recv_arr[0];
    m2[0][CHUNK_SIZE - 1] = calc_lax_value(up, down, left, right, u, v, dt, dx);
    // printf("TASK: %d :: t = %f b = %f l = %f r = %f \n", taskid, top_recv_arr[CHUNKSIZE - 1], m1[1][CHUNKSIZE - 1], m1[0][CHUNKSIZE - 2], right_recv_arr[0]);

    up = CHUNK_SIZE == 1 ? top_recv_arr[CHUNK_SIZE - 1] : m1[CHUNK_SIZE - 2][CHUNK_SIZE - 1];
    down = bottom_recv_arr[CHUNK_SIZE - 1];
    left = CHUNK_SIZE == 1 ? left_recv_arr[CHUNK_SIZE - 1] : m1[CHUNK_SIZE - 1][CHUNK_SIZE - 2];
    right = right_recv_arr[CHUNK_SIZE - 1];
    m2[CHUNK_SIZE - 1][CHUNK_SIZE - 1] = calc_lax_value(up, down, left, right, u, v, dt, dx);
}